# Known values are at 2 decimals, items_averaged = 1, rounding = "up_or_down".

agg <- function(data, ...) {
  grim_aggregate(data, digits_x = 2, items_averaged = 1, ...)
}

pass_rate <- function(n, x = 1, ...) {
  agg(data.frame(x = x, n = n), ...)$data$p_pass
}

get <- function(res, quantity, col = "value") {
  res$summary[[col]][res$summary$quantity == quantity]
}

test_that("pass rates match the known values", {
  expect_equal(pass_rate(30), 0.30)
  expect_equal(pass_rate(50), 0.50)
  # N = 40 and 60: S/N lands on .x25 / .x75 ties, which up_or_down accepts
  # in both directions (confirmed against scrutiny 0.6.2)
  expect_equal(pass_rate(40), 0.60)
  expect_equal(pass_rate(60), 0.60)
  expect_equal(pass_rate(c(100, 137, 250)), c(1, 1, 1))
})

test_that("pass rate uses items_averaged, as an argument or a column", {
  two <- grim_aggregate(
    data.frame(x = 1, n = 15),
    digits_x = 2,
    items_averaged = 2
  )
  expect_equal(two$data$p_pass, pass_rate(30))
  dat <- data.frame(
    x = c(1, 1),
    n = c(15, 30),
    items_averaged = c(2, 1),
    digits_x = 2
  )
  res <- grim_aggregate(dat)
  expect_equal(res$data$p_pass, c(0.30, 0.30))
})

test_that("digits_x and items_averaged: one way, not both, not neither", {
  dat <- data.frame(x = 3.4, n = 30, digits_x = 2)
  expect_error(
    grim_aggregate(dat, digits_x = 2, items_averaged = 1),
    "both a column"
  )
  expect_error(
    grim_aggregate(data.frame(x = 3.4, n = 30), items_averaged = 1),
    "`digits_x` is required"
  )
  expect_error(grim_aggregate(dat), "`items_averaged` is required")
  expect_error(
    grim_aggregate(transform(dat, items = 1)),
    "rename it to `items_averaged`"
  )
  expect_error(grim_aggregate(dat, items = 1), "renamed to `items_averaged`")
})

test_that("F = 0 gives zero lower bounds throughout", {
  dat <- data.frame(x = c(3.40, 2.50, 1.20), n = c(30, 20, 25))
  res <- agg(dat)
  expect_equal(get(res, "F"), 0)
  expect_equal(get(res, "theta", "lower"), 0)
  expect_equal(get(res, "theta", "upper"), 1)
  expect_equal(get(res, "O_E"), 0)
  expect_equal(get(res, "lambda_mle"), 0)
  expect_equal(get(res, "lambda_support_8", "lower"), 0)
  expect_equal(get(res, "lambda_support_32", "lower"), 0)
  expect_equal(get(res, "p_more_than_random"), 1)
})

test_that("all informative values failing gives lambda_mle = 1", {
  # 5.19 fails at N = 28; 3.14 fails at N = 30; 2.47 fails at N = 25
  dat <- data.frame(x = c(5.19, 3.14, 2.47), n = c(28, 30, 25))
  res <- agg(dat)
  expect_equal(get(res, "F"), 3)
  expect_equal(get(res, "k"), 3)
  expect_equal(get(res, "lambda_mle"), 1)
  expect_equal(get(res, "theta", "lower"), 1)
  expect_equal(get(res, "lambda_support_8", "upper"), 1)
})

test_that("equal q gives lambda_mle = min(1, O/E)", {
  # N = 30 throughout, so q = 0.7 for every value; one failure in four
  dat <- data.frame(x = c(5.10, 5.13, 5.17, 5.14), n = 30)
  res <- agg(dat)
  expect_equal(get(res, "F"), 1)
  expect_equal(get(res, "lambda_mle"), min(1, get(res, "O_E")),
               tolerance = 1e-8)
})

test_that("a failing value with q = 0 is flagged, warned about and excluded", {
  # At N = 100 every 2-decimal ending is attainable (q = 0): uninformative
  dat <- data.frame(x = c(3.14, 2.47, 5.19), n = c(100, 25, 28))
  res0 <- agg(dat)
  expect_false(res0$data$flag_impossible[1])
  expect_false(res0$data$informative[1])
  expect_equal(get(res0, "n_uninformative"), 1)

  # scrutiny itself never fails a value on a full grid, so such a failure
  # can only reach grim_aggregate() from outside (a pasted verdict, a wrong
  # N). Stub the observed verdicts to exercise the flag.
  local_mocked_bindings(grim_observed = function(...) c(FALSE, FALSE, FALSE))
  expect_warning(res <- agg(dat), "extraction or reporting error")
  expect_true(res$data$flag_impossible[1])
  expect_false(res$data$informative[1])
  expect_equal(get(res, "n_flagged_impossible"), 1)
  expect_equal(get(res, "k"), 2)
  expect_equal(get(res, "F"), 2)
})

test_that("q and mean_q are exactly 0 when every ending passes", {
  res <- agg(data.frame(x = c(3.14, 2.47), n = c(100, 250)))
  expect_identical(res$data$q, c(0, 0))
  expect_identical(get(res, "mean_q"), 0)
})

test_that("Poisson-binomial pmf sums to 1 and matches dbinom for equal q", {
  pmf <- strait:::poisbinom_pmf(c(0.1, 0.5, 0.73, 0.2, 0.99))
  expect_equal(sum(pmf), 1)
  expect_length(pmf, 6)
  for (q in c(0.05, 0.3, 0.7)) {
    expect_equal(strait:::poisbinom_pmf(rep(q, 12)), dbinom(0:12, 12, q))
  }
})

test_that("O/E is calibrated: mean approx lambda under the mixture", {
  skip_on_cran() # a 600-run simulation; CI and local runs still do it
  set.seed(20260927)
  ns <- c(12, 17, 23, 28, 33, 39, 41, 46, 52, 57,
          63, 68, 71, 77, 84, 88, 93, 97, 30, 45)
  sim <- function(lambda, reps = 200) {
    replicate(reps, {
      fake <- runif(length(ns)) < lambda
      # non-genuine: uniform on the reporting grid
      x_fake <- sample(100:599, length(ns), replace = TRUE) / 100
      # genuine: an integer sum divided by N, rounded
      x_real <- vapply(ns, function(nn) sample(nn:(5 * nn), 1) / nn, 0)
      x <- round(ifelse(fake, x_fake, x_real), 2)
      get(agg(data.frame(x = x, n = ns)), "O_E")
    })
  }
  expect_equal(mean(sim(1)), 1, tolerance = 0.08)
  expect_equal(mean(sim(0.5)), 0.5, tolerance = 0.08)
  expect_equal(mean(sim(0)), 0)
})

test_that("malformed input errors instead of passing GRIM", {
  # "3,14" used to become NA, which then counted as consistent
  expect_error(agg(data.frame(x = c("3.14", "3,14"), n = 30)), "decimal_places")
  expect_error(agg(data.frame(x = c(3.14, NA), n = 30)), "missing")
  expect_error(agg(data.frame(x = 3.14, n = c(30, NA))), "`n`")
  expect_error(agg(data.frame(x = 3.14, n = Inf)), "`n`")
  expect_error(agg(data.frame(x = 3.14, n = 30.5)), "`n`")
  expect_error(
    grim_aggregate(
      data.frame(x = 3.14, n = 30),
      digits_x = NA,
      items_averaged = 1
    ),
    "`digits_x`"
  )
  expect_error(
    grim_aggregate(
      data.frame(x = 3.14, n = 30),
      digits_x = 2,
      items_averaged = NA
    ),
    "`items_averaged`"
  )
  expect_error(agg(data.frame(x = numeric(0), n = numeric(0))), "no rows")
})

test_that("conf_level and support_c are validated; CI label keeps decimals", {
  dat <- data.frame(x = c(5.19, 3.14, 2.40, 3.33), n = c(28, 30, 25, 27))
  expect_error(agg(dat, ci = TRUE, conf_level = 90), "conf_level")
  expect_error(agg(dat, ci = TRUE, conf_level = 0), "conf_level")
  expect_error(agg(dat, support_c = 1), "support_c")
  expect_error(agg(dat, support_c = -8), "support_c")
  res <- agg(dat, ci = TRUE, conf_level = 0.975)
  expect_true("lambda_ci_97.5" %in% res$summary$quantity)
})

test_that("the CI is NA with a warning when even lambda = 1 is rejected", {
  # three failures where each fails with probability 0.05 at most
  expect_warning(
    ci <- strait:::lambda_ci(c(1, 1, 1), rep(0.05, 3), 0.9),
    "rejected"
  )
  expect_equal(ci, c(NA_real_, NA_real_))
})

test_that("support intervals nest and contain the MLE; CI is ordered", {
  dat <- data.frame(
    x = c(5.19, 3.14, 2.40, 3.33, 1.25, 4.44),
                    n = c(28, 30, 25, 27, 36, 18))
  res <- agg(dat, ci = TRUE)
  mle <- get(res, "lambda_mle")
  s8 <- c(get(res, "lambda_support_8", "lower"),
          get(res, "lambda_support_8", "upper"))
  s32 <- c(get(res, "lambda_support_32", "lower"),
           get(res, "lambda_support_32", "upper"))
  expect_true(s32[1] <= s8[1] && s8[1] <= mle && mle <= s8[2] && s8[2] <= s32[2])
  ci <- c(get(res, "lambda_ci_90", "lower"), get(res, "lambda_ci_90", "upper"))
  expect_true(0 <= ci[1] && ci[1] <= ci[2] && ci[2] <= 1)
})

test_that("custom g: final-digit function and grid vector agree", {
  dat <- data.frame(x = c(5.19, 3.14, 2.40), n = c(28, 30, 25))
  gf <- function(d) ifelse(d %in% c(0, 5), 3, 1)
  gv <- gf((0:99) %% 10)
  a <- agg(dat, g = gf)
  b <- agg(dat, g = gv)
  expect_equal(a$data$p_pass, b$data$p_pass)
  expect_error(agg(dat, g = 1:10), "length 10\\^d")
})

test_that("baseline: SIR = (O/E) / tau when eps = tau * q", {
  dat <- data.frame(x = c(5.19, 3.14, 2.40, 3.33), n = c(28, 30, 25, 27))
  tau <- 0.05
  calls <- 0
  res <- agg(dat, eps = function(q, n, d) {
    calls <<- calls + 1
    tau * q
  })
  expect_equal(calls, 1)
  expect_equal(get(res, "SIR"), get(res, "O_E") / tau)
  expect_true(get(res, "p_vs_baseline") >= 0 && get(res, "p_vs_baseline") <= 1)
})

test_that("the eps column is there when no value is informative", {
  res <- agg(data.frame(x = 3.14, n = 100), eps = 0.01)
  expect_equal(get(res, "k"), 0)
  expect_equal(res$data$eps, 0.01)
})

test_that("grim_aggregate_by and grim_sensitivity return one row per unit", {
  dat <- data.frame(
    paper = rep(c("A", "B"), each = 3),
    x = c(5.19, 3.14, 2.40, 3.40, 2.50, 1.20),
    n = c(28, 30, 25, 30, 20, 25)
  )
  by <- grim_aggregate_by(dat, by = "paper", digits_x = 2, items_averaged = 1)
  expect_s3_class(by, "tbl_df")
  expect_equal(nrow(by), 2)
  expect_true(all(c("paper", "k", "F", "O_E", "theta_lower",
                    "lambda_support_8_upper") %in% names(by)))
  sens <- grim_sensitivity(
    dat,
    list(heap = function(d) ifelse(d == 0, 5, 1)),
    digits_x = 2,
    items_averaged = 1
  )
  expect_s3_class(sens, "tbl_df")
  expect_equal(sens$counterfactual, c("uniform", "heap"))
  expect_equal(sens$O_E[1], get(agg(dat), "O_E"))
})

test_that("grim_aggregate_by keeps NA groups and rejects clashing names", {
  dat <- data.frame(
    paper = c("A", "A", NA, NA, "B"),
    x = c(5.19, 3.14, 2.47, 3.40, 2.50),
    n = c(28, 30, 25, 30, 20)
  )
  by <- grim_aggregate_by(dat, by = "paper", digits_x = 2, items_averaged = 1)
  expect_equal(by$paper, c("A", "B", NA))
  expect_equal(sum(by$n_values), 5)
  names(dat)[1] <- "F"
  expect_error(
    grim_aggregate_by(dat, by = "F", digits_x = 2, items_averaged = 1),
    "clash"
  )
})
