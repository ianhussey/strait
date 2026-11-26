# The quasi-integer floor against brute force: every multiset of n - 1
# integers in a window around the mean, with the free observation forced by
# the mean. The window is wide enough that wider placements only increase
# the sum of squares.

brute_quasi_min_sd <- function(n, mean, window = 5) {
  v <- floor(mean)
  # one row per multiset, also when n - 1 = 1 (multisets() then returns a
  # single row)
  fixed <- matrix(multisets((v - window):(v + window + 1), n - 1), ncol = n - 1)
  free <- n * mean - rowSums(fixed)
  ss <- rowSums((fixed - mean)^2) + (free - mean)^2

  sqrt(min(ss) / (n - 1))
}

test_that("sd_min_quasi_integer() is the brute-force minimum at generic means", {
  skip_on_cran() # 2,406 enumerations; CI and local runs still do it
  for (n in 2:7) {
    means <- 0.013 + (0:400) * 0.0137
    expect_equal(
      sd_min_quasi_integer(means, n),
      vapply(means, \(mean) brute_quasi_min_sd(n, mean), numeric(1)),
      tolerance = 1e-8
    )
  }
})

test_that("sd_bounds_sample() builds a sample that attains the quasi-integer floor", {
  for (n in 2:7) {
    for (mean in 0.013 + (0:40) * 0.137) {
      x <- sd_bounds_sample(l = 0, u = 7, n = n, mean = mean, which = "min")
      expect_equal(mean(x), mean, tolerance = 1e-9)
      expect_equal(sd(x), sd_min_quasi_integer(mean, n), tolerance = 1e-9)
      expect_lte(sum(abs(x - round(x)) > 1e-9), 1) # at most one non-integer
    }
  }
})

test_that("at GRIM-consistent means the quasi-integer floor is the integer floor", {
  for (n in 2:7) {
    means <- (0:42) / n
    expect_equal(sd_min_quasi_integer(means, n), sd_min_integer(means, n),
                 tolerance = 1e-9)
    expect_equal(
      sd_min_quasi_integer(means, n),
      vapply(means, \(mean) brute_quasi_min_sd(n, mean), numeric(1)),
      tolerance = 1e-8
    )
  }
})

test_that("the integer floor interpolated off the grid overstates it by g(1 - g)", {
  # in sums of squares: n d (1 - d) - [n d (1 - d) - g (1 - g)] = g (1 - g),
  # with g the fractional part of n * mean
  for (n in 2:29) {
    means <- 0.017 + (0:299) * 0.0173
    g <- n * means - floor(n * means)
    g[g > 1 - 1e-9] <- 0
    naive_ss <- sd_min_integer(means, n)^2 * (n - 1)
    quasi_ss <- sd_min_quasi_integer(means, n)^2 * (n - 1)
    expect_true(all(quasi_ss <= naive_ss + 1e-12))
    expect_equal(naive_ss - quasi_ss, g * (1 - g), tolerance = 1e-8)
  }
})
