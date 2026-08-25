# The analytic certification engine (R/certify-analytic.R). Exactness rests
# on the completeness of the demotion decomposition, so the tests are
# agreement tests: against brute force where brute force is affordable, and
# against the two older, independent engines (constructive search / corridor
# DP, and the full lattice) everywhere else. brimmest(engine = "both") errors
# on any disagreement, which makes whole sweeps one-line assertions.

brute_states <- function(n, W) {
  grid <- as.matrix(expand.grid(rep(list(0:W), n)))
  data.frame(S = rowSums(grid), Q = rowSums(grid^2))
}

test_that(".hit_deficit reproduces brute-force attainability on small designs", {
  for (n in c(3L, 5L)) for (W in c(4L, 6L)) {
    b <- brute_states(n, W)
    for (s in unique(b$S)) {
      qs <- unique(b$Q[b$S == s])
      Qx <- .q_max_cnt(s, W, n)
      for (q in seq(min(qs), max(qs))) {
        expect_identical(
          .hit_deficit(s, W, n, Qx - q, Qx - q),
          q %in% qs
        )
      }
    }
  }
})

test_that(".analytic_window agrees with .witness_search run to exhaustion", {
  set.seed(42)
  W <- 9L; n <- 14L
  for (rep in 1:300) {
    S <- sample.int(n * W - 1L, 1L)
    Qn <- .q_min_int(S, n); Qx <- .q_max_cnt(S, W, n)
    k_lo <- sample(seq(Qn, Qx), 1L)
    k_hi <- min(Qx, k_lo + sample(0:10, 1L))
    got_s <- .witness_search(W, n, S, k_lo, k_hi, budget = 1e7)$possible
    expect_identical(.analytic_window(W, n, S, k_lo, k_hi), got_s)
  }
})

test_that("the analytic engine settles the canonical blind-spot tuples", {
  # passes bounds, GRIM and GRIMMER, yet no integer sample produces it
  expect_false(brimmest(1, 5, 9, mean = 1.3, sd = 0.9, digits = 1)$possible)
  expect_true(brimmest(1, 5, 9, mean = 3.0, sd = 1.0, digits = 1)$possible)
  # top-centre phantoms at two decimals (impossible despite passing GRIMMER)
  r <- brimmest(1, 7, 25, mean = c(4.00, 4.00, 4.00), digits = 2,
                sd = c(2.93, 2.94, 3.00))
  expect_identical(r$possible, c(TRUE, FALSE, TRUE))
})

test_that("engine = 'both' finds no disagreement across rounding rules (lattice route)", {
  set.seed(7)
  rules <- c("half_up", "native", "trunc")     # incl. bankers and truncation
  means <- round(runif(120, 1, 7), 1)
  sds <- round(runif(120, 0, 2.3), 1)
  expect_silent(
    r <- brimmest(1, 7, 9, mean = means, sd = sds, digits = 1,
                  rounding = rules, engine = "both")
  )
  # and the verdicts equal the default engine's
  a <- brimmest(1, 7, 9, mean = means, sd = sds, digits = 1,
                rounding = rules, engine = "analytic")
  expect_identical(r$possible, a$possible)
  expect_identical(r$rules, a$rules)
})

test_that("engine = 'both' finds no disagreement on a wide design (targeted route)", {
  set.seed(8)
  W <- 63; n <- 50
  sd_max_of <- function(S) {
    a <- S %/% W; b <- S - a * W
    sqrt((a * W^2 + b^2 - S^2 / n) / (n - 1))
  }
  S <- sample(seq(2, n * W - 2), 40)
  means <- round(S / n, 2)
  sds <- round(pmax(0, vapply(S, sd_max_of, 0) - runif(40, 0, 2) * 0.01), 2)
  expect_silent(
    r <- brimmest(0, 63, 50, mean = means, sd = sds, digits = 2,
                  rounding = c("half_up", "native", "trunc"), engine = "both")
  )
  s <- brimmest(0, 63, 50, mean = means, sd = sds, digits = 2,
                rounding = c("half_up", "native", "trunc"), engine = "search")
  expect_identical(r$possible, s$possible)
  expect_identical(r$rules, s$rules)
})

test_that("brimmest_multiple passes engine through as a call constant", {
  reports <- data.frame(mean = c(3.0, 1.3), sd = c(1.0, 0.9))
  r <- brimmest_multiple(reports, l = 1, u = 5, n = 9, digits = 1,
                         engine = "both")
  expect_identical(r$possible, c(TRUE, FALSE))
  expect_error(
    brimmest_multiple(cbind(reports, engine = "analytic"),
                      l = 1, u = 5, n = 9, digits = 1),
    "constant"
  )
})

test_that("lattice boundary bug is fixed: trunc-rounded reports on exact rationals", {
  # 14 sevens, 1 five, 10 ones: mean = 113/25 = 4.52 exactly, sd = 2.9597...
  # trunc-rounds to (4.52, 2.95). The lattice route used to file 4.52 under
  # 4.51 (its scaled double is 451.99999999999994) and declare this real
  # sample impossible; engine = "both" caught the disagreement.
  x <- c(rep(7, 14), 5, rep(1, 10))
  expect_identical(round(mean(x), 10), 4.52)
  r <- brimmest(1, 7, 25, mean = 4.52, sd = 2.95, digits = 2,
                rounding = "trunc", engine = "both")
  expect_true(r$possible)
})

test_that(".round_scaled matches scrutiny's forward rounding away from boundaries", {
  set.seed(11)
  x <- runif(400, -8, 8)                       # generic values: no boundaries
  for (d in 0:2) {
    expect_identical(.round_scaled(x, d, "half_up"),
                     round(scrutiny::round_up_from(x, d, 5, symmetric = FALSE) * 10^d))
    expect_identical(.round_scaled(x, d, "native"), round(round(x, d) * 10^d))
    expect_identical(.round_scaled(x, d, "floor"),
                     round(scrutiny::round_floor(x, d) * 10^d))
    expect_identical(.round_scaled(x, d, "trunc"),
                     round(scrutiny::round_trunc(x, d) * 10^d))
  }
})

test_that(".round_scaled is exact on the boundary values scrutiny misfiles", {
  # 113/25 = 4.52 exactly; its scaled double 451.99999999999994 makes
  # scrutiny::round_trunc(113/25, 2) return 4.51. The snapped scaled rounding
  # must file every rule's boundary case where the exact rational belongs.
  expect_identical(.round_scaled(113/25, 2, "trunc"), 452)
  expect_identical(.round_scaled(113/25, 2, "floor"), 452)
  expect_identical(.round_scaled(-113/25, 2, "trunc"), -452)
  expect_identical(.round_scaled(-113/25, 2, "ceiling"), -452)
  # half boundaries: toward +Inf (half_up), toward -Inf (half_down), to even
  expect_identical(.round_scaled(4.525, 2, "half_up"), 453)
  expect_identical(.round_scaled(-1.25, 1, "half_up"), -12)
  expect_identical(.round_scaled(-1.25, 1, "half_down"), -13)
  expect_identical(.round_scaled(4.525, 2, "native"), 452)
})
