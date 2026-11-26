# The unimodal ceiling against brute force: every multiset of n values on the
# scale, kept if its frequency distribution is unimodal (unobserved values
# count as frequency 0).

# Rises, then falls, with ties allowed: no decrease followed by an increase.
is_unimodal_counts <- function(counts) {
  steps <- sign(diff(counts))
  steps <- steps[steps != 0]
  !any(head(steps, -1) < 0 & tail(steps, -1) > 0)
}

# The largest and smallest sample SD among the unimodal samples, and the
# smallest among all samples, at each attainable mean.
unimodal_brute_force <- function(values, n) {
  samples <- multisets(values, n)
  unimodal <- apply(samples, 1, \(x) {
    is_unimodal_counts(tabulate(match(x, values), length(values)))
  })
  means <- round(rowMeans(samples), 10)
  sds <- apply(samples, 1, stats::sd)
  data.frame(
    mean = sort(unique(means)),
    max_unimodal = tapply(sds[unimodal], means[unimodal], max),
    min_unimodal = tapply(sds[unimodal], means[unimodal], min),
    min_all = tapply(sds, means, min)
  )
}

test_that("sd_max_unimodal() is the largest SD of a unimodal sample at every mean", {
  for (design in list(c(2, 2), c(2, 9), c(4, 2), c(4, 7), c(4, 12), c(6, 8))) {
    top <- design[1]
    n <- design[2]
    bf <- unimodal_brute_force(0:top, n)
    expect_equal(
      sd_max_unimodal(bf$mean, n, lower = 0, upper = top),
      unname(bf$max_unimodal),
      tolerance = 1e-9
    )
  }
})

test_that("without the mean, sd_max_unimodal() is the largest over all means", {
  for (top in 1:5) {
    for (n in 2:12) {
      bf <- unimodal_brute_force(0:top, n)
      expect_equal(
        sd_max_unimodal(n = n, lower = 0, upper = top),
        max(bf$max_unimodal),
        tolerance = 1e-9
      )
    }
  }
})

test_that("sd_max_unimodal() is NA at a mean no sample of whole numbers has", {
  expect_true(is.na(sd_max_unimodal(mean = 2.55, n = 10, lower = 1, upper = 5)))
  expect_false(is.na(sd_max_unimodal(mean = 2.5, n = 10, lower = 1, upper = 5)))
})

test_that("shape = 'unimodal' lowers the ceiling and keeps the floor", {
  bf <- unimodal_brute_force(1:5, 11)
  for (i in seq_len(nrow(bf))) {
    any_shape <- sd_bounds(l = 1, u = 5, n = 11, mean = bf$mean[i],
                           granularity = "integer")
    unimodal <- sd_bounds(l = 1, u = 5, n = 11, mean = bf$mean[i],
                          granularity = "integer", shape = "unimodal")
    expect_equal(unimodal$max_sd, unname(bf$max_unimodal[i]), tolerance = 1e-9)
    expect_equal(unimodal$min_sd, any_shape$min_sd)
    expect_equal(unimodal$min_sd, unname(bf$min_unimodal[i]), tolerance = 1e-9)
    expect_lte(unimodal$max_sd, any_shape$max_sd + 1e-9)
  }
})

test_that("shape = 'unimodal' works on the grid of mean-scored items", {
  # means of two items on a 1-4 scale lie on a grid of spacing 1/2
  bf <- unimodal_brute_force(seq(1, 4, by = 0.5), 8)
  ceiling <- vapply(bf$mean, \(mean) {
    sd_bounds(l = 1, u = 4, n = 8, mean = mean, granularity = "integer",
              scoring = "meanscored", n_items = 2, shape = "unimodal")$max_sd
  }, numeric(1))
  expect_equal(ceiling, unname(bf$max_unimodal), tolerance = 1e-9)
})

test_that("a rounded mean gives the largest ceiling over its rounding interval", {
  # n = 30 on a 1-5 scale, mean 4.5 reported with one decimal, so the sum
  # of the observations minus n is 104, 105, or 106
  b <- sd_bounds(l = 1, u = 5, n = 30, mean = 4.5, digits_mean = 1,
                 rounding = "up_or_down", granularity = "integer",
                 shape = "unimodal")
  each <- sd_max_unimodal(mean = 1 + (104:106) / 30, n = 30, lower = 1,
                          upper = 5)
  expect_equal(b$max_sd, max(each), tolerance = 1e-9)
  expect_equal(round(each, 2), c(1.11, 1.01, 1.01))
})

test_that("sd_max_unimodal() reproduces known values", {
  expect_equal(round(sd_max_unimodal(n = 10, lower = 1, upper = 5), 2), 1.49)
  expect_equal(round(sd_max_unimodal(4.5, n = 30, lower = 1, upper = 5), 2), 1.01)
  # the counts 1, 2, 2, 2, 3, 2, 2 on a 1-7 scale attain the ceiling at
  # their mean
  counts <- c(1, 2, 2, 2, 3, 2, 2)
  x <- rep(1:7, counts)
  expect_equal(sd_max_unimodal(mean(x), n = 14, lower = 1, upper = 7), sd(x))
})

test_that("brimmer() flags an SD that only a non-unimodal sample can have", {
  # c(rep(1, 20), rep(7, 9), 6): mean 2.97, sd 2.83, bimodal
  args <- list(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
               sd = 2.83, digits_sd = 2, granularity = "integer")
  expect_true(do.call(brimmer, args)$consistent)
  unimodal <- do.call(brimmer, c(args, shape = "unimodal"))
  expect_false(unimodal$consistent)
  expect_identical(unimodal$failed_tests, "bounds")
})

test_that("shape = 'unimodal' refuses constraints it does not support", {
  expect_error(sd_bounds(l = 1, u = 5, n = 10, mean = 3, shape = "unimodal"),
               "granularity")
  expect_error(sd_bounds(l = 1, n = 10, mean = 3, granularity = "integer",
                         shape = "unimodal"), "requires")
  expect_error(sd_bounds(l = 1, u = 5, a = 2, n = 10, mean = 3,
                         granularity = "integer", shape = "unimodal"),
               "not supported")
  expect_error(sd_bounds(l = 3, u = 15, n = 20, mean = 9,
                         granularity = "integer", scoring = "sumscored",
                         n_items = 3, alpha = 0.8, shape = "unimodal"),
               "not supported")
})

test_that("sd_max_unimodal() validates its inputs and snaps a near-whole range", {
  expect_error(sd_max_unimodal(3, n = 10, lower = 5, upper = 1), "at least")
  expect_error(sd_max_unimodal(3, n = 1, lower = 1, upper = 5), "at least")
  expect_error(sd_max_unimodal(3, n = 10.5, lower = 1, upper = 5), "whole")
  expect_identical(
    sd_max_unimodal(c(3, NA), n = 10, lower = 1, upper = 5)[2],
    NA_real_
  )
  expect_equal(
    sd_bounds(
      l = 1 + 5e-10,
      u = 5,
      n = 10,
      mean = 3,
      granularity = "integer",
      shape = "unimodal"
    )$max_sd,
    sd_max_unimodal(3, n = 10, lower = 1, upper = 5)
  )
})

test_that("shape = 'unimodal' works without a mean and in brimmer_map()", {
  expect_equal(
    sd_bounds(
      l = 1,
      u = 5,
      n = 10,
      granularity = "integer",
      shape = "unimodal"
    )$max_sd,
    sd_max_unimodal(n = 10, lower = 1, upper = 5)
  )

  reports <- data.frame(mean = 2.97, sd = 2.83)
  out <- brimmer_map(
    reports,
    l = 1,
    u = 7,
    n = 30,
    digits_mean = 2,
    digits_sd = 2,
    granularity = "integer",
    shape = "unimodal"
  )
  expect_false(out$consistent)
})

test_that("umbrella_data() tests cells against the unimodal ceiling on the full grid", {
  any_shape <- umbrella_data(l = 1, u = 5, n = 10, digits = 1)
  unimodal <- umbrella_data(l = 1, u = 5, n = 10, digits = 1, shape = "unimodal")

  # same cells, up to the ordinary ceiling
  expect_identical(unimodal[c("mean", "sd")], any_shape[c("mean", "sd")])
  # fewer pass, and some fail only for want of a unimodal sample
  expect_true(all(!unimodal$consistent | any_shape$consistent))
  expect_true(any(any_shape$consistent & !unimodal$in_bounds))
  # the bounds reported per cell are the unimodal ones
  first <- unimodal[1, ]
  expect_equal(
    first$max_sd,
    sd_bounds(l = 1, u = 5, n = 10, mean = first$mean, digits_mean = 1,
              rounding = "up_or_down", granularity = "integer",
              shape = "unimodal")$max_sd
  )
})

test_that("sd_region_data() takes shape = 'unimodal' for the integer rule only", {
  any_shape <- sd_region_data(l = 1, u = 5, n = 10, rule = "integer", digits = 1)
  unimodal <- sd_region_data(l = 1, u = 5, n = 10, rule = "integer", digits = 1,
                             shape = "unimodal")
  key <- \(d) paste(d$mean, d$sd)
  expect_true(all(key(unimodal) %in% key(any_shape)))
  expect_lt(nrow(unimodal), nrow(any_shape))
  expect_error(
    sd_region_data(l = 1, u = 5, n = 10, rule = "quasi", shape = "unimodal"),
    "only for"
  )
})
