test_that("chisq2 reproduces the reference table (BDI-II, M = 25, SD = 8)", {
  r <- sd_plausibility(
    0,
    63,
    n = c(20, 50, 100, 500),
    mean = 25,
    sd = 8,
    p = c(0.01, 0.001),
    condition_on_mean = FALSE
  )
  expect_s3_class(r, "tbl_df")
  expect_equal(round(r$rho[1], 3), 0.052)
  got <- matrix(round(c(r$lower, r$upper), 2), ncol = 2)
  # rows alternate p = 0.01, 0.001 within each n
  ref <- rbind(
    c(5.24, 10.86),
    c(4.47, 11.90),
    c(6.26, 9.78),
    c(5.74, 10.40),
               c(6.77, 9.25), c(6.39, 9.68), c(7.45, 8.55), c(7.27, 8.74))
  expect_equal(got, ref, tolerance = 0.006)
  expect_equal(round(bb_cumulants(63, 25 / 63, r$rho[1])$g2, 3), -0.242)
})

test_that("the pmf matches extraDistr and the moments match closed forms", {
  skip_if_not_installed("extraDistr")
  m <- 20; mu <- 0.3; rho <- 0.1; phi <- 1 / rho - 1
  expect_equal(strait:::bb_pmf(m, mu, rho),
               extraDistr::dbbinom(0:m, m, mu * phi, (1 - mu) * phi))
  cu <- bb_cumulants(m, mu, rho)
  expect_equal(cu$mean, m * mu)
  expect_equal(cu$k2, m * mu * (1 - mu) * (1 + (m - 1) * rho))
})

test_that("the exact DP matches brute-force enumeration", {
  m <- 4
  mu <- 0.4
  rho <- 0.1
  n <- 5
  pr <- strait:::bb_pmf(m, mu, rho)
  g <- as.matrix(expand.grid(rep(list(0:m), n)))
  w <- apply(g, 1, function(r) prod(pr[r + 1]))
  s2 <- apply(g, 1, stats::var)
  bf <- tapply(w, round(s2 * n * (n - 1)), sum)
  ex <- s2_exact(m, mu, rho, n = n)[["5"]]
  expect_s3_class(ex, "tbl_df")
  expect_equal(unname(as.numeric(bf)), ex$prob)
  expect_equal(sort(unique(round(s2, 10))), round(ex$s2, 10))

  # conditioning on the sum
  keep <- rowSums(g) == 8
  bfc <- tapply(w[keep], round(s2[keep] * n * (n - 1)), sum)
  exc <- s2_exact(m, mu, rho, n = n, sums = 8)[["5"]]
  expect_equal(unname(as.numeric(bfc)) / sum(bfc), exc$prob)
  expect_error(s2_exact(m, mu, rho, n = 1), "`n`")
})

test_that("cumulants of S^2 match the exact distribution, including kappa_3", {
  for (n in c(2, 3, 7)) {
    ex <- s2_exact(5, 0.35, 0.2, n = n)[[1]]
    sc <- strait:::s2_cumulants(
      strait:::cumulants_from_pmf(strait:::bb_pmf(5, 0.35, 0.2)),
      n
    )
    e1 <- sum(ex$s2 * ex$prob)
    expect_equal(e1, sc$k1)
    expect_equal(sum((ex$s2 - e1)^2 * ex$prob), sc$k2)
    expect_equal(sum((ex$s2 - e1)^3 * ex$prob), sc$k3)
  }
})

test_that("normal-theory checks: nu = n - 1 and kappa_3 = 8 sigma^6/(n-1)^2", {
  # a normal population has g2 = 0 and kappa_3 = kappa_4 = kappa_6 = 0
  cu <- list(k2 = 4, k3 = 0, k4 = 0, k6 = 0, g2 = 0)
  sc <- strait:::s2_cumulants(cu, 11)
  expect_equal(2 * sc$k1^2 / sc$k2, 10)
  expect_equal(sc$k3, 8 * 4^3 / 10^2)
})

test_that("approximate methods agree with the exact method at moderate n", {
  e <- sd_plausibility(
    0,
    10,
    n = 12,
    mean = 4,
    rho = 0.1,
    method = "exact",
    condition_on_mean = FALSE
  )
  for (meth in c("chisq2", "chisq3", "cf")) {
    a <- sd_plausibility(
      0,
      10,
      n = 12,
      mean = 4,
      rho = 0.1,
      method = meth,
      condition_on_mean = FALSE
    )
    expect_equal(a$lower, e$lower, tolerance = 0.08)
    expect_equal(a$upper, e$upper, tolerance = 0.08)
  }
})

test_that("conditioning on the mean is the default and narrows the bounds", {
  u <- sd_plausibility(
    0,
    20,
    n = 15,
    mean = 2,
    rho = 0.1,
    condition_on_mean = FALSE
  )
  c <- sd_plausibility(0, 20, n = 15, mean = 2, rho = 0.1)
  expect_true(c$conditional)
  expect_lt(c$upper - c$lower, u$upper - u$lower)
  ec <- sd_plausibility(0, 20, n = 15, mean = 2, rho = 0.1, method = "exact")
  expect_equal(c$lower, ec$lower, tolerance = 0.1)
  expect_equal(c$upper, ec$upper, tolerance = 0.1)
  expect_error(
    sd_plausibility(0, 20, n = 15, mean = 2, rho = 0.1, method = "cf"),
    "condition_on_mean = FALSE"
  )
})

test_that("the mixture over rho widens the bounds and tau = 0 is the point", {
  a <- sd_plausibility(0, 63, n = 50, mean = 25, rho = 0.05)
  b <- sd_plausibility(0, 63, n = 50, mean = 25, rho = 0.05, tau = 1e-9)
  w <- sd_plausibility(0, 63, n = 50, mean = 25, rho = 0.05, tau = 0.03)
  expect_equal(a$lower, b$lower, tolerance = 1e-6)
  expect_lt(w$lower, a$lower)
  expect_gt(w$upper, a$upper)
})

test_that("truncated-normal nodes stay in [0, 1) and have the right centre", {
  nodes <- strait:::rho_nodes(0.01, 0.05, 41)
  expect_true(all(nodes >= 0 & nodes < 1))
  expect_length(nodes, 41)
  expect_equal(strait:::rho_nodes(0.3, 0, 41), 0.3)
})

test_that("hard bounds are the integer floor and the sharp ceiling", {
  r <- sd_plausibility(1, 7, n = 30, mean = 3.5, rho = 0.1)
  expect_equal(r$hard_min, sd_min_integer(2.5, 30))
  expect_equal(r$hard_max, sd_max_structure_s(3.5, 30, 1, 7))
  # with digits_mean, the floor is the minimum over every sum that rounds
  r2 <- sd_plausibility(1, 7, n = 30, mean = 3.5, rho = 0.1, digits_mean = 1)
  expect_lte(r2$hard_min, r$hard_min)
})

test_that("digits_mean is recycled with mean, one value per row", {
  both <- sd_plausibility(
    1,
    7,
    n = 30,
    mean = c(3.5, 3.52),
    rho = 0.1,
    digits_mean = c(1, 2)
  )
  one <- sd_plausibility(1, 7, n = 30, mean = 3.52, rho = 0.1, digits_mean = 2)
  expect_equal(both[2, ], one)
  expect_error(
    sd_plausibility(
      1,
      7,
      n = 30,
      mean = c(3.5, 3.52, 3.6),
      rho = 0.1,
      digits_mean = c(1, 2)
    ),
    "common length"
  )
})

test_that("lower quantiles below the integer floor are clipped and flagged", {
  # a mean half-way between integers forces SD >= ~0.5, while a nearly
  # binomial DGP on a 2-step scale at tiny n can go lower
  r <- sd_plausibility(
    0,
    2,
    n = 4,
    mean = 0.5,
    rho = 0,
    p = 0.01,
    condition_on_mean = FALSE
  )
  expect_gte(r$lower, r$hard_min)
  expect_match(r$clipped, "lower")
  expect_match(r$note, "clipped")
})

test_that("edge cases: mean at or near the limits", {
  z <- sd_plausibility(0, 63, n = 30, mean = 0, rho = 0.05)
  expect_equal(c(z$lower, z$upper, z$hard_max), c(0, 0, 0))
  expect_match(z$note, "scale limit")
  # the sd path used to divide by mu (1 - mu) = 0 and report missing input
  zs <- sd_plausibility(0, 63, n = 30, mean = 63, sd = 0)
  expect_equal(c(zs$lower, zs$upper, zs$hard_max), c(0, 0, 0))
  expect_match(zs$note, "scale limit")
  near <- sd_plausibility(0, 63, n = 30, mean = 0.2, rho = 0.05)
  expect_true(is.finite(near$lower) && is.finite(near$upper))
  expect_lte(near$upper, near$hard_max)
  out <- sd_plausibility(0, 63, n = 30, mean = 70, rho = 0.05)
  expect_true(is.na(out$lower))
  expect_match(out$note, "outside")
})

test_that("edge cases: rho <= 0 and rho >= 1", {
  neg <- sd_plausibility(0, 10, n = 20, mean = 5, sd = 1) # under-dispersed
  expect_lt(neg$rho, 1e-12)
  expect_match(neg$note, "under-dispersed")
  zero <- sd_plausibility(0, 10, n = 20, mean = 5, rho = 0)
  expect_equal(neg$lower, zero$lower)
  one <- sd_plausibility(0, 10, n = 20, mean = 5, rho = 1)
  expect_match(one$note, "endpoints")
  expect_equal(one$sigma, 5)
})

test_that("a two-point scale is Bernoulli: rho is not identified", {
  a <- sd_plausibility(0, 1, n = 30, mean = 0.4, rho = 0.3)
  b <- sd_plausibility(0, 1, n = 30, mean = 0.4, sd = 0.5)
  expect_true(is.na(a$rho))
  expect_match(a$note, "binary")
  expect_false(grepl("endpoints", a$note))
  expect_equal(a$sigma, sqrt(0.4 * 0.6))
  expect_equal(a[c("lower", "upper")], b[c("lower", "upper")])
  expect_error(
    rho_from_summary(0, 1, n = 30, mean = 0.4, sd = 0.5),
    "not identified"
  )
})

test_that("edge cases: n = 2, missing inputs", {
  two <- sd_plausibility(0, 10, n = 2, mean = 5, rho = 0.1, method = "exact")
  expect_true(is.finite(two$lower) && is.finite(two$upper))
  two_a <- sd_plausibility(0, 10, n = 2, mean = 5, rho = 0.1)
  expect_true(is.finite(two_a$upper))
  miss <- sd_plausibility(0, 10, n = 20, mean = c(5, NA), rho = 0.1)
  expect_equal(nrow(miss), 2)
  expect_true(is.na(miss$lower[2]))
  expect_match(miss$note[2], "missing")
  # a non-integer l is fine when u - l is whole
  expect_true(is.finite(
    sd_plausibility(0.5, 3.5, n = 20, mean = 2, rho = 0.1)$upper
  ))
})

test_that("invalid input errors clearly", {
  expect_error(
    sd_plausibility(1, 5.5, n = 20, mean = 3, rho = 0.1),
               "whole number")
  expect_error(sd_plausibility(0, 10, n = 20, mean = 5), "exactly one")
  expect_error(
    sd_plausibility(0, 10, n = 20, mean = 5, sd = 2, rho = 0.1),
    "exactly one"
  )
  expect_error(
    sd_plausibility(c(0, 1), 10, n = 20, mean = 5, rho = 0.1),
    "single"
  )
  expect_error(sd_plausibility(0, 10, n = 1, mean = 5, rho = 0.1), "`n`")
  expect_error(sd_plausibility(0, 10, n = Inf, mean = 5, rho = 0.1), "`n`")
  expect_error(sd_plausibility(0, 10, n = 20.5, mean = 5, rho = 0.1), "`n`")
  expect_error(
    sd_plausibility(0, 10, n = 20, mean = 5, rho = 0.1, p = 0.6),
    "`p`"
  )
  expect_error(
    sd_plausibility(0, 10, N = 20, M = 5, rho = 0.1),
    "`N` was renamed to `n`"
  )
  expect_error(
    rho_from_summary(c(0, 1), 10, n = 20, mean = 5, sd = 2),
    "single"
  )
})

test_that("the exact DP runs once per rho node, not once per p", {
  real <- strait:::sq_joint
  calls <- 0
  local_mocked_bindings(sq_joint = function(...) {
    calls <<- calls + 1
    real(...)
  })
  r <- sd_plausibility(
    0,
    6,
    n = 8,
    mean = 2.5,
    rho = 0.1,
    p = c(0.005, 0.01, 0.05),
    method = "exact"
  )
  expect_equal(nrow(r), 3)
  expect_equal(calls, 1)
  calls <- 0
  sd_plausibility(
    0,
    6,
    n = 8,
    mean = 2.5,
    rho = 0.1,
    tau = 0.02,
    K = 5,
    p = c(0.005, 0.01, 0.05),
    method = "exact"
  )
  expect_equal(calls, 5)
})

test_that("exact method refuses infeasible sizes", {
  expect_error(
    sd_plausibility(0, 63, n = 200, mean = 25, rho = 0.05, method = "exact"),
    "approximate method"
  )
})

test_that("rho_from_summary inverts the moment formula and keeps negatives", {
  m <- 63; mu <- 0.4; rho <- 0.07
  sd <- sqrt(m * mu * (1 - mu) * (1 + (m - 1) * rho))
  r <- rho_from_summary(0, 63, n = 50, mean = m * mu, sd = sd)
  expect_equal(r$rho_hat, rho)
  expect_equal(r$mu_hat, mu)
  expect_lt(rho_from_summary(0, 10, n = 50, mean = 5, sd = 1)$rho_hat, 0)
  # var_rho is Var(S^2) / [m (m - 1) mu (1 - mu)]^2
  sc <- strait:::s2_cumulants(strait:::cumulants_from_pmf(
    strait:::bb_pmf(m, mu, rho)), 50)
  expect_equal(r$var_rho, sc$k2 / (m * (m - 1) * mu * (1 - mu))^2)
})

test_that("rho_corpus_fit recovers a flat rho and predicts tau", {
  skip_if_not_installed("metafor")
  skip_if_not_installed("extraDistr")
  set.seed(3)
  k <- 80; m <- 30
  mu <- runif(k, 0.15, 0.85)
  n <- sample(40:300, k, replace = TRUE)
  sd <- vapply(seq_len(k), function(i) {
      x <- extraDistr::rbbinom(
        n[i],
        m,
        mu[i] * (1 / 0.08 - 1),
                             (1 - mu[i]) * (1 / 0.08 - 1))
    stats::sd(x)
    },
    0
  )
  # population means, as a corpus would report sample means
  corpus <- data.frame(l = 0, u = m, n = n, mean = m * mu, sd = sd)
  fit <- rho_corpus_fit(corpus)
  expect_s3_class(fit, "rho_corpus")
  pr <- predict(fit, c(0.3, 0.5, 0.7))
  expect_equal(pr$rho, rep(0.08, 3), tolerance = 0.3)
  expect_true(all(pr$tau_total >= pr$tau))
  old <- data.frame(L = 0, U = m, M = m * mu, SD = sd, N = n)
  expect_error(rho_corpus_fit(old), "were renamed")
})
