consistent_within <- function(l, u, ...) {
  brimmer(l = l, u = u, n = 30, digits_mean = 2, digits_sd = 2, ...)$consistent
}

test_that("an SD near the ceiling forces both scale limits into the data", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 2.97,
    digits_mean = 2,
    sd = 2.83,
    digits_sd = 2,
    Z = "integer"
  )
  expect_equal(r$min_at_most, 1)
  expect_equal(r$max_at_least, 7)
})

test_that("a small SD forces only the neighbouring values", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 4.00,
    digits_mean = 2,
    sd = 0.53,
    digits_sd = 2,
    Z = "integer"
  )
  expect_equal(r$min_at_most, 3)
  expect_equal(r$max_at_least, 5)
})

test_that("integer results are exactly where brimmer() changes its verdict", {
  reports <- expand.grid(
    mean = c(2.10, 3.50, 4.80, 6.20),
    sd = c(0.80, 1.40, 2.10)
  )
  for (i in seq_len(nrow(reports))) {
    m <- reports$mean[i]
    s <- reports$sd[i]
    r <- implied_range(
      l = 1,
      u = 7,
      n = 30,
      mean = m,
      digits_mean = 2,
      sd = s,
      digits_sd = 2,
      Z = "integer"
    )
    if (is.na(r$min_at_most)) {
      next
    }
    lo <- r$min_at_most
    hi <- r$max_at_least
    expect_true(consistent_within(lo, 7, mean = m, sd = s, Z = "integer"))
    expect_true(consistent_within(1, hi, mean = m, sd = s, Z = "integer"))
    if (lo + 1 < m) {
      expect_false(consistent_within(
        lo + 1,
        7,
        mean = m,
        sd = s,
        Z = "integer"
      ))
    }
    if (hi - 1 > m) {
      expect_false(consistent_within(
        1,
        hi - 1,
        mean = m,
        sd = s,
        Z = "integer"
      ))
    }
  }
})

test_that("continuous results sit at the threshold of the bounds test", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 4.00,
    digits_mean = 2,
    sd = 1.90,
    digits_sd = 2
  )
  expect_false(consistent_within(r$min_at_most, 7, mean = 4, sd = 1.9))
  expect_true(consistent_within(r$min_at_most - 1e-4, 7, mean = 4, sd = 1.9))
  # symmetric about the scale midpoint
  expect_equal(r$min_at_most + r$max_at_least, 8, tolerance = 1e-5)
})

test_that("a mean-scored composite is bounded on its 1/n_items grid", {
  r <- implied_range(
    l = 1,
    u = 5,
    n = 20,
    mean = 3.00,
    digits_mean = 2,
    sd = 1.20,
    digits_sd = 2,
    Z = "integer",
    scoring = "meanscored",
    n_items = 3
  )
  # on [7/3, 5] the ceiling at mean 3 is 1.185, below 1.195: some response
  # is at most 2, the grid value below 7/3
  expect_equal(r$min_at_most, 2)
  expect_equal(r$max_at_least, 4)
})

test_that("nothing forced: the bounds fall back to the mean", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 4.00,
    digits_mean = 2,
    sd = 0.00,
    digits_sd = 2,
    Z = "integer"
  )
  expect_equal(r$min_at_most, 4)
  expect_equal(r$max_at_least, 4)
})

test_that("attained extremes are returned as given", {
  r <- implied_range(
    a = 2,
    u = 7,
    n = 30,
    mean = 4.00,
    digits_mean = 2,
    sd = 1.50,
    digits_sd = 2
  )
  expect_equal(r$min_at_most, 2)
})

test_that("an inconsistent report has no implied range", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 3.51,
    digits_mean = 2,
    sd = 1.00,
    digits_sd = 2,
    Z = "integer"
  )
  expect_true(is.na(r$min_at_most) && is.na(r$max_at_least))
  expect_match(r$note, "inconsistent")
})

test_that("missing inputs are errors", {
  expect_error(implied_range(l = 1, u = 7, n = 30, mean = 4), "sd")
  expect_error(implied_range(l = 1, n = 30, mean = 4, sd = 1), "limits")
})

test_that("mean-scored data move the limit along the 1/n_items grid", {
  # SD 1.95 at mean 3 on 1-5 needs a 1: on [4/3, 5] the ceiling is about 1.87
  r <- implied_range(
    l = 1,
    u = 5,
    n = 20,
    mean = 3.00,
    digits_mean = 2,
    sd = 1.95,
    digits_sd = 2,
    Z = "integer",
    scoring = "meanscored",
    n_items = 3
  )
  expect_equal(r$min_at_most, 1)
  expect_equal(r$max_at_least, 5)
})

test_that("a missing input gives NA results", {
  r <- implied_range(
    l = 1,
    u = 7,
    n = 30,
    mean = 3,
    digits_mean = 1,
    sd = NA,
    digits_sd = 1
  )
  expect_true(is.na(r$min_at_most) && is.na(r$max_at_least))
  expect_match(r$note, "missing input")
})
