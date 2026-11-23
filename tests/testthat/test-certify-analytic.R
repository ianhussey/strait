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

test_that("hit_deficit reproduces brute-force attainability on small designs", {
  for (n in c(3L, 5L)) for (W in c(4L, 6L)) {
    b <- brute_states(n, W)
    for (s in unique(b$S)) {
      qs <- unique(b$Q[b$S == s])
      Qx <- q_max_cnt(s, W, n)
      for (q in seq(min(qs), max(qs))) {
        expect_identical(
          hit_deficit(s, W, n, Qx - q, Qx - q),
          q %in% qs
        )
      }
    }
  }
})

test_that("analytic_window agrees with witness_search run to exhaustion", {
  set.seed(42)
  W <- 9L; n <- 14L
  for (rep in 1:300) {
    S <- sample.int(n * W - 1L, 1L)
    Qn <- q_min_int(S, n); Qx <- q_max_cnt(S, W, n)
    k_lo <- sample(seq(Qn, Qx), 1L)
    k_hi <- min(Qx, k_lo + sample(0:10, 1L))
    got_s <- witness_search(W, n, S, k_lo, k_hi, budget = 1e7)$possible
    expect_identical(analytic_window(W, n, S, k_lo, k_hi), got_s)
  }
})

test_that("the analytic engine settles the canonical blind-spot tuples", {
  # passes bounds, GRIM and GRIMMER, yet no integer sample produces it
  expect_false(brimmest(1, 5, 9, mean = 1.3, sd = 0.9, digits_mean = 1, digits_sd = 1)$possible)
  expect_true(brimmest(1, 5, 9, mean = 3.0, sd = 1.0, digits_mean = 1, digits_sd = 1)$possible)
  # top-centre phantoms at two decimals (impossible despite passing GRIMMER)
  r <- brimmest(1, 7, 25, mean = c(4.00, 4.00, 4.00), digits_mean = 2, digits_sd = 2,
                sd = c(2.93, 2.94, 3.00))
  expect_identical(r$possible, c(TRUE, FALSE, TRUE))
})

test_that("engine = 'both' finds no disagreement across rounding rules (lattice route)", {
  set.seed(7)
  rules <- c("up", "even", "trunc")     # incl. bankers and truncation
  means <- round(runif(120, 1, 7), 1)
  sds <- round(runif(120, 0, 2.3), 1)
  expect_silent(
    r <- brimmest(1, 7, 9, mean = means, sd = sds, digits_mean = 1, digits_sd = 1,
                  rounding = rules, engine = "both")
  )
  # and the verdicts equal the default engine's
  a <- brimmest(1, 7, 9, mean = means, sd = sds, digits_mean = 1, digits_sd = 1,
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
    r <- brimmest(0, 63, 50, mean = means, sd = sds, digits_mean = 2, digits_sd = 2,
                  rounding = c("up", "even", "trunc"), engine = "both")
  )
  s <- brimmest(0, 63, 50, mean = means, sd = sds, digits_mean = 2, digits_sd = 2,
                rounding = c("up", "even", "trunc"), engine = "search")
  expect_identical(r$possible, s$possible)
  expect_identical(r$rules, s$rules)
})

test_that("brimmest_map passes engine through as a call constant", {
  reports <- data.frame(mean = c(3.0, 1.3), sd = c(1.0, 0.9))
  r <- brimmest_map(reports, l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
                         engine = "both")
  expect_identical(r$possible, c(TRUE, FALSE))
  expect_error(
    brimmest_map(cbind(reports, engine = "analytic"),
                      l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1),
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
  r <- brimmest(1, 7, 25, mean = 4.52, sd = 2.95, digits_mean = 2, digits_sd = 2,
                rounding = "trunc", engine = "both")
  expect_true(r$possible)
})

test_that("round_exact matches roundwork's forward rounding away from boundaries", {
  set.seed(11)
  num <- round(runif(400, -8000, 8000)) * 7 + 3  # whole numbers, odd thirteenths
  den <- 13
  x <- num / den
  for (d in 0:2) {
    p <- num * 10^d
    expect_identical(round_exact(p, den, "up")$lo,
                     round(roundwork::round_up(x, d) * 10^d))
    expect_identical(round_exact(p, den, "floor")$lo,
                     round(roundwork::round_floor(x, d) * 10^d))
    expect_identical(round_exact(p, den, "trunc")$lo,
                     round(roundwork::round_trunc(x, d) * 10^d))
    expect_identical(round_exact(p, den, "anti_trunc")$lo,
                     round(roundwork::round_anti_trunc(x, d) * 10^d))
  }
})

test_that("round_exact files boundary values where the exact value belongs", {
  # 113/25 = 4.52 exactly; scaled as a double it is 451.99999999999994
  expect_identical(round_exact(113 * 100, 25, "trunc")$lo, 452)
  expect_identical(round_exact(113 * 100, 25, "floor")$lo, 452)
  expect_identical(round_exact(-113 * 100, 25, "trunc")$lo, -452)
  expect_identical(round_exact(-113 * 100, 25, "ceiling")$lo, -452)
  # ties: toward +Inf (up), toward -Inf (down); "even" lets them go either way
  expect_identical(round_exact(4525, 10, "up")$lo, 453)
  expect_identical(round_exact(-125, 10, "up")$lo, -12)
  expect_identical(round_exact(-125, 10, "down")$lo, -13)
  expect_identical(unlist(round_exact(4525, 10, "even")), c(lo = 452, hi = 453))
  # an SD as sqrt(num / den): sd = 0.475 exactly is a tie at 2 decimals
  expect_identical(round_exact(4 * 47.5^2, 4, "up", root = TRUE)$lo, 48)
  expect_identical(round_exact(4 * 47.5^2, 4, "down", root = TRUE)$lo, 47)
  # a hair above the tie is no tie: (47.5^2 * 4 + 1) / 4 rounds up either way
  expect_identical(round_exact(4 * 47.5^2 + 1, 4, "down", root = TRUE)$lo, 48)
})

test_that("an SD a hair off a rounding tie is not snapped onto it", {
  # binary scale, n = 33017, S = 11353: sd = 0.475000000005, above the tie at
  # 2 decimals, so "down" must still give 0.48. A 1e-9 snap used to file it
  # under the tie and round it down to 0.47.
  n <- 33017; S <- 11353
  nss <- n * S - S^2
  expect_identical(round_exact(nss * 10^4, n * (n - 1), "down",
                                root = TRUE)$lo, 48)
})

test_that("both routes agree under 'even' ties and at a truncated 0", {
  # a mean of exactly 0.25 rounds to 0.3 under one tie direction
  expect_true(brimmest(0, 1, 4, mean = 0.3, sd = 0.5, digits_mean = 1,
                       digits_sd = 1, rounding = "even", engine = "both")$possible)
  expect_identical(
    brimmest(0, 1, 4, mean = 0.3, sd = 0.5, digits_mean = 1, digits_sd = 1,
             rounding = "even")$possible,
    brimmest(0, 1, 4, mean = 0.3, sd = 0.5, digits_mean = 1, digits_sd = 1,
             rounding = "even", max_cells = 1)$possible
  )
  # c(-3, 1, 1): mean -0.33, sd 2.31, which truncate to 0 and 2
  for (mc in c(2e7, 1)) {
    expect_true(brimmest(-3, 1, 3, mean = 0, sd = 2, digits_mean = 0,
                         digits_sd = 0, rounding = "trunc",
                         max_cells = mc)$possible)
    # under anti_trunc only an exact 0 reports as 0, and no pair on -2..2
    # has mean 0 with an SD that rounds away from zero to 1
    expect_false(brimmest(-2, 2, 2, mean = 0, sd = 1, digits_mean = 0,
                          digits_sd = 0, rounding = "anti_trunc",
                          max_cells = mc)$possible)
  }
})

test_that("brimmest() validates its inputs on every route", {
  expect_error(brimmest(1, 5.6, 3, mean = 6, sd = 0, digits_mean = 1,
                        digits_sd = 1), "response grid")
  expect_error(brimmest(1, 5, 4.5, mean = 3, sd = 0, digits_mean = 1,
                        digits_sd = 1), "whole number")
  expect_error(brimmest(1, 5, 2, mean = 2, sd = -0.1, digits_mean = 1,
                        digits_sd = 1, max_cells = 1), "negative")
  expect_error(brimmest(1, 5, 9, mean = 3.04, sd = 1, digits_mean = 1,
                        digits_sd = 1), "more decimal places")
  expect_error(brimmest(1, 5, 9, mean = 3, sd = 1, digits_mean = 1,
                        digits_sd = 1, rounding = "half_up"), '"up"')
  r <- brimmest(1, 5, 9, mean = c(3, NA), sd = c(1, 1), digits_mean = 1,
                digits_sd = 1)
  expect_identical(r$possible, c(TRUE, NA))
})
