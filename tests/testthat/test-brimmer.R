test_that("brimmer() returns a verdict with POMP columns", {
  r <- brimmer(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2,
               sd = 2.83, sd_digits = 2, Z = "integer")
  expect_true(all(c("consistent", "failed_tests", "min_sd", "max_sd",
                    "feasible", "in_scale_range", "grim", "grimmer",
                    "sd_in_bounds", "pomp_mean", "pomp_sd_parity",
                    "pomp_sd_sharp", "note") %in% names(r)))
  expect_true(r$consistent)
  expect_equal(r$failed_tests, "")
})

test_that("an SD below the floor fails the bounds test", {
  r <- brimmer(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2,
               sd = 0.10, sd_digits = 2, Z = "quasiinteger")
  expect_false(r$consistent)
  expect_true(grepl("bounds", r$failed_tests))
  expect_lt(r$pomp_sd_sharp, 0)          # below the band -> negative
})

test_that("a GRIM-impossible mean fails grim, not feasibility", {
  r <- brimmer(l = 1, u = 7, n = 30, mean = 3.51, mean_digits = 2,
               sd = 1.00, sd_digits = 2, Z = "integer")
  expect_false(r$consistent)
  expect_true(grepl("grim", r$failed_tests))
  # 3.51 is inside [1, 7]; only its granularity is impossible, so the
  # generic feasibility residual must not also fire
  expect_false(grepl("feasibility", r$failed_tests))
  expect_true(r$in_scale_range)
})

test_that("pomp_sd_parity is linear and mean-independent", {
  p1 <- brimmer(l = 1, u = 7, n = 30, mean = 4, sd = 1.0, sd_digits = 6,
                rounding = NULL, Z = "continuous")
  p2 <- brimmer(l = 1, u = 7, n = 30, mean = 4, sd = 2.0, sd_digits = 6,
                rounding = NULL, Z = "continuous")
  p3 <- brimmer(l = 1, u = 7, n = 30, mean = 6, sd = 1.0, sd_digits = 6,
                rounding = NULL, Z = "continuous")
  expect_equal(p2$pomp_sd_parity, 2 * p1$pomp_sd_parity)
  expect_equal(p3$pomp_sd_parity, p1$pomp_sd_parity)
})

test_that("pomp_sd_sharp is 1 at the ceiling and 0 at the floor", {
  mx <- sd_bounds(l = 1, u = 5, n = 9, mean = 3)$max_sd
  mn <- sd_bounds(l = 1, u = 5, n = 9, mean = 2.44, Z = "quasiinteger")$min_sd
  sh1 <- brimmer(l = 1, u = 5, n = 9, mean = 3, sd = mx, sd_digits = 9,
                 rounding = NULL, Z = "continuous")$pomp_sd_sharp
  sh0 <- brimmer(l = 1, u = 5, n = 9, mean = 2.44, sd = mn, sd_digits = 9,
                 rounding = NULL, Z = "quasiinteger")$pomp_sd_sharp
  expect_equal(sh1, 1, tolerance = 1e-6)
  expect_equal(sh0, 0, tolerance = 1e-6)
})

test_that("brimmer_multiple() equals a row-wise loop, incl. duplicates", {
  df <- expand.grid(mean = c(2.90, 2.97, 3.50), sd = c(0.10, 1.50, 2.83))
  df <- rbind(df, df[1:3, ])
  mult <- brimmer_multiple(df, l = 1, u = 7, n = 30, mean_digits = 2,
                           sd_digits = 2, Z = "integer")
  roww <- do.call(rbind, lapply(seq_len(nrow(df)), function(i)
    brimmer(l = 1, u = 7, n = 30, mean = df$mean[i], mean_digits = 2,
            sd = df$sd[i], sd_digits = 2, Z = "integer")))
  expect_equal(mult$consistent, roww$consistent)
  expect_equal(mult$min_sd, roww$min_sd)
  expect_equal(mult$pomp_sd_sharp, roww$pomp_sd_sharp)
})

test_that("brimmer_multiple() honours include_inputs and column/constant collisions", {
  df <- data.frame(mean = c(3, 4), sd = c(1.5, 1.6))
  keep <- brimmer_multiple(df, l = 1, u = 7, n = 30, mean_digits = 1,
                           sd_digits = 1, Z = "integer")
  drop <- brimmer_multiple(df, l = 1, u = 7, n = 30, mean_digits = 1,
                           sd_digits = 1, Z = "integer",
                           include_inputs = FALSE)
  expect_true(all(c("mean", "sd") %in% names(keep)))
  expect_false("mean" %in% setdiff(names(drop), names(keep)))
  expect_equal(nrow(drop), nrow(df))
  df$n <- 30
  expect_error(brimmer_multiple(df, n = 30, l = 1, u = 7,
                                mean_digits = 1, sd_digits = 1, sd = 1))
})

test_that("brim() runs the mean-side tests only", {
  ok <- brim(l = 1, u = 7, n = 30, mean = 2.97, mean_digits = 2)
  expect_true(ok$consistent)
  expect_identical(ok$failed_tests, "")
  expect_true(ok$in_scale_range)
  # SD-side columns are omitted, not returned as NA
  expect_false(any(c("min_sd", "max_sd", "sd_in_bounds", "grimmer",
                     "pomp_sd_parity", "pomp_sd_sharp") %in% names(ok)))
  expect_true(all(c("in_scale_range", "band_lo", "band_hi",
                    "pomp_mean") %in% names(ok)))
})

test_that("brim() rejects a mean outside the feasible band", {
  bad <- brim(l = 1, u = 7, n = 30, mean = 7.5, mean_digits = 1)
  expect_false(bad$consistent)
  expect_identical(bad$failed_tests, "in_scale_range")
  expect_false(bad$in_scale_range)
  expect_equal(c(bad$band_lo, bad$band_hi), c(1, 7))
})

test_that("brim() reports the band narrowed by attained extremes", {
  # both extremes attained pins one observation at each, so the band is
  # [a + W/n, b - W/n] = [1.2, 6.8] for a = 1, b = 7, n = 30
  r <- brim(a = 1, b = 7, n = 30, mean = 1.10, mean_digits = 2)
  expect_equal(c(r$band_lo, r$band_hi), c(1.2, 6.8))
  expect_false(r$consistent)
  # the same mean passes against bare scale limits, where the band is [1, 7]
  expect_true(brim(l = 1, u = 7, n = 30, mean = 1.10, mean_digits = 2)$consistent)
})

test_that("brim() agrees with brimmer()'s mean-side verdict", {
  for (mu in c(2.97, 3.51, 7.5)) {
    b  <- brim(l = 1, u = 7, n = 30, mean = mu, mean_digits = 2, Z = "integer")
    bm <- brimmer(l = 1, u = 7, n = 30, mean = mu, mean_digits = 2,
                  sd = NULL, Z = "integer")
    expect_identical(b$in_scale_range, bm$in_scale_range)
    expect_identical(b$grim, bm$grim)
  }
})

test_that("in_scale_range and grim are independent, and neither double-counts", {
  # the four combinations on a 1-7 scale at n = 30, under Z = "integer":
  #   3.50  in range,     GRIM ok      -> passes
  #   3.51  in range,     30*x = 105.3 -> grim only
  #   7.50  out of range, 30*x = 225   -> in_scale_range only
  #   7.51  out of range, 30*x = 225.3 -> both
  f <- function(mu) brim(l = 1, u = 7, n = 30, mean = mu, mean_digits = 2,
                         Z = "integer")
  expect_identical(f(3.50)$failed_tests, "")
  expect_identical(f(3.51)$failed_tests, "grim")
  expect_identical(f(7.50)$failed_tests, "in_scale_range")
  expect_identical(f(7.51)$failed_tests, "in_scale_range,grim")
  # a GRIM failure must not be reported as an out-of-range mean
  expect_true(f(3.51)$in_scale_range)
  # and an out-of-range mean must not be reported as a GRIM failure
  expect_true(f(7.50)$grim)
})

test_that("brimmer() nests the mean-side tests under the SD-side ones", {
  f <- function(mu) brimmer(l = 1, u = 7, n = 30, mean = mu, mean_digits = 2,
                            sd = 1.00, sd_digits = 2, Z = "integer")
  # every brim() failure is still reported by brimmer(), plus the SD-side ones
  expect_identical(f(3.51)$failed_tests, "grim,grimmer")
  expect_identical(f(7.50)$failed_tests, "in_scale_range,grimmer")
  expect_identical(f(7.51)$failed_tests, "in_scale_range,grim,grimmer")
  for (mu in c(3.50, 3.51, 7.50, 7.51)) {
    expect_identical(f(mu)$in_scale_range,
                     brim(l = 1, u = 7, n = 30, mean = mu, mean_digits = 2,
                          Z = "integer")$in_scale_range)
  }
})

test_that("in_scale_range is NA when there is no mean, and Z-independent", {
  sd_only <- brimmer(l = 1, u = 7, n = 30, sd = 2.00, sd_digits = 2)
  expect_true(is.na(sd_only$in_scale_range))
  # the predicate is bounds-only, so granularity must not change it
  for (z in c("continuous", "integer", "quasiinteger")) {
    expect_false(brimmer(l = 1, u = 7, n = 30, mean = 7.5, mean_digits = 1,
                         sd = 1.0, sd_digits = 1, Z = z)$in_scale_range)
    expect_true(brimmer(l = 1, u = 7, n = 30, mean = 3.5, mean_digits = 1,
                        sd = 1.0, sd_digits = 1, Z = z)$in_scale_range)
  }
})

test_that("brim() requires a mean and brimmer() requires mean or sd", {
  expect_error(brim(l = 1, u = 7, n = 30), "mean is required")
  expect_error(brimmer(l = 1, u = 7, n = 30), "reported sd, a reported mean")
})

test_that("feasibility is a residual, fired only when no named test explains it", {
  # a 5-item composite on [5, 10] at n = 6 with alpha = 0.90: the
  # alpha-amplified floor exceeds the alpha-conditional ceiling, so no sample
  # exists for a reason that is neither an out-of-range mean nor a GRIM failure.
  # SD 1.1 rather than 1.0: integer sums with mean 6.0 and SD 1.0 at n = 6 fail
  # GRIMMER's parity test (sum 36, sum of squares 221).
  r <- brimmer(l = 5, u = 10, n = 6, mean = 6.0, mean_digits = 1,
               sd = 1.1, sd_digits = 1, Z = "integer",
               scoring = "sumscored", n_items = 5, alpha = 0.90)
  expect_false(r$consistent)
  expect_identical(r$failed_tests, "feasibility")
  expect_false(r$feasible)
  # the mean itself is fine; feasibility is carrying the alpha contradiction
  expect_true(r$in_scale_range)
})
