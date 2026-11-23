test_that("sd_bounds_curve() is hole-free under quasi-integer", {
  cv <- sd_bounds_curve(l = 1, u = 7, n = 30, Z = "quasiinteger")
  expect_true(all(c("mean", "min_sd", "max_sd", "feasible", "pomp_mean",
                    "parity_max", "ceil_parity", "floor_parity") %in% names(cv)))
  expect_false(any(is.na(cv$min_sd[cv$feasible])))
  expect_false(any(is.na(cv$max_sd[cv$feasible])))
})

test_that("umbrella_data() passing set matches an independent GRIMMER + bounds screen", {
  um <- umbrella_data(n = 12, l = 1, u = 3, digits = 1)
  expect_true(all(c("mean", "sd", "min_sd", "max_sd", "in_bounds",
                    "grimmer", "consistent") %in% names(um)))
  # all three verdict categories must be present, or this proves little
  expect_true(any(um$consistent))
  expect_true(any(um$in_bounds & !is.na(um$grimmer) & !um$grimmer))
  expect_true(any(!um$in_bounds))
  indep <- with(um, {
    gm <- as.logical(scrutiny::grimmer(x = mean, sd = sd, n = 12,
                                       digits_x = 1, digits_sd = 1,
                                       items_averaged = 1,
                                       data_type = "integer",
                                       rounding = "up_or_down"))
    ib <- (sd + 0.05) >= min_sd - 1e-9 & (sd - 0.05) <= max_sd + 1e-9
    gm & ib
  })
  expect_equal(um$consistent, indep)
})

test_that("the plot functions return ggplot objects", {
  cv  <- sd_bounds_curve(l = 1, u = 7, n = 30, Z = "quasiinteger")
  pts <- brimmer_multiple(
    data.frame(mean = c(3, 4, 5), sd = c(2.5, 0.2, 1.8)),
    l = 1, u = 7, n = 30, digits_mean = 1, digits_sd = 1, Z = "quasiinteger")
  expect_s3_class(sd_bounds_plot(cv, reports = pts), "ggplot")
  expect_s3_class(sd_bounds_plot(sd_region_data(1, 7, 30), reports = pts), "ggplot")
  expect_s3_class(sd_bounds_pomp_plot(cv, reports = pts, reference = "parity"), "ggplot")
  expect_s3_class(sd_bounds_pomp_plot(cv, reports = pts, reference = "sharp"), "ggplot")
  um <- umbrella_data(n = 12, l = 1, u = 7, digits = 2)
  expect_s3_class(umbrella_plot(um, bounds = sd_bounds_curve(l = 1, u = 7, n = 12,
                                                             Z = "quasiinteger")), "ggplot")
})

test_that("umbrella_plot() fills a band", {
  cv <- sd_bounds_curve(l = 1, u = 5, n = 12)
  p <- umbrella_plot(cv, fill = "red")
  expect_s3_class(p, "ggplot")
  band <- ggplot2::layer_data(p, 2)
  expect_equal(band$ymin, ifelse(cv$feasible, cv$min_sd, NA))
  expect_equal(band$ymax, ifelse(cv$feasible, cv$max_sd, NA))
  expect_equal(unique(band$fill), "red")
  # the region is the outline of the points: every consistent tuple lies in it
  # once rounding is undone. A reported mean stands for the interval of half a
  # step around it (1.08 at n = 12 is 13 / 12, where the SD reaches 0.2887 but
  # is at most 0.277 at 1.08 itself), and so does a reported SD.
  um <- umbrella_data(n = 12, l = 1, u = 5, digits = 2)
  ok <- um[um$consistent, ]
  cv <- cv[cv$feasible, ]
  h <- 0.005
  around <- function(m, col, f) vapply(m, function(x) {
    inside <- cv$mean >= x - h & cv$mean <= x + h
    f(c(stats::approx(cv$mean, cv[[col]], c(x - h, x + h), rule = 2)$y,
        cv[[col]][inside]))
  }, numeric(1))
  lo <- around(ok$mean, "min_sd", min)
  hi <- around(ok$mean, "max_sd", max)
  expect_true(all(ok$sd >= lo - h - 1e-9 & ok$sd <= hi + h + 1e-9))
  # a band from sd_region_data() works the same way
  sm <- sd_region_data(l = 1, u = 5, n = 12, rule = "mean_naive_floor")
  band <- ggplot2::layer_data(umbrella_plot(sm), 2)
  expect_equal(band$ymin, sm$min_sd)
  expect_equal(band$ymax, sm$max_sd)
  # bounds are dashed lines on top, for a band as for a lattice
  dashed <- ggplot2::layer_data(umbrella_plot(sm, bounds = cv), 3)
  expect_true(all(dashed$linetype == "dashed"))
  expect_error(umbrella_plot(cv, style = "tiles"), "needs a lattice")
  expect_error(umbrella_plot(data.frame(mean = 1, max_sd = 1)), "needs `mean` and `sd`")
  lat <- sd_region_data(1, 5, 12, rule = "integer", digits = 1)
  expect_error(umbrella_plot(lat, style = "tiles"), "full umbrella_data")
})

# Regression: the plot functions used to drop the rows of a band where no SD
# is possible, and ggplot2 joins the rows on either side of a missing one, so a
# ribbon filled in the gap and a line crossed it. ggplot2 breaks a ribbon or
# line at an NA row instead, so the rows must still be there, as NA.
test_that("bands are not drawn across a gap", {
  # a reported alpha leaves 87 means next to each limit with no feasible SD
  d <- sd_region_data(0, 3, 7, rule = "alpha", n_items = 2, alpha = 0.70,
                      by = 0.001)
  gap <- d$mean[is.na(d$max_sd)]
  breaks_at_gap <- function(p, i, col) {
    ld <- ggplot2::layer_data(p, i)
    in_gap <- ld$x %in% gap
    sum(in_gap) == length(gap) && all(is.na(ld[[col]][in_gap]))
  }
  expect_true(breaks_at_gap(umbrella_plot(d), 2, "ymax"))
  expect_true(breaks_at_gap(umbrella_plot(d, bounds = d), 3, "y"))
  expect_true(breaks_at_gap(sd_bounds_plot(d, shade = "inside"), 1, "ymax"))
  expect_true(breaks_at_gap(sd_bounds_plot(d, shade = "inside"), 2, "y"))
  quasi <- sd_region_data(0, 3, 7, rule = "quasi", n_items = 2)
  expect_true(breaks_at_gap(sd_region_plot(quasi, reference = d), 3, "y"))
})

test_that("sd_region_plot() draws only the consistent tuples of a lattice", {
  um <- umbrella_data(n = 12, l = 1, u = 3, digits = 1)
  pts <- ggplot2::layer_data(sd_region_plot(um), 1)
  expect_equal(nrow(pts), sum(um$consistent))
})

# Regression: sd_bounds_plot() used to infer the mean-grid spacing from
# median(diff(curve$mean)). sd_bounds_curve()'s grid is deliberately not
# uniform - it adds every kink of the 1/(n * n_items) lattice plus a pair of
# neighbours 1e-9 away - so once the kinks outnumber the uniform grid (about
# n * (u - l) > 333 at the default `by`) the median fell BELOW the plain grid
# spacing. band_polygon() then read every ordinary interval as a gap and the
# feasible region was drawn as thousands of one-column slivers instead of one
# band. sd_bounds_curve() now records the spacing it used and sd_bounds_plot()
# reads it.

# the white knock-out rings sd_bounds_plot() builds under shade = "outside"
n_rings <- function(p) {
  for (ly in p$layers) {
    if (
      inherits(ly$geom, "GeomPolygon") &&
        is.data.frame(ly$data) &&
        "ring" %in% names(ly$data)
    ) {
      return(length(unique(ly$data$ring)))
    }
  }
  NA_integer_
}

test_that("sd_bounds_curve() records the grid spacing it used", {
  expect_equal(attr(sd_bounds_curve(l = 1, u = 7, n = 30), "step"), 6 / 1000)
  expect_equal(attr(sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.01), "step"), 0.01)
})

test_that("a hole-free band is one ring, however dense the kink lattice", {
  # n = 30 was already fine; 60 and 200 are the cases the median guess broke,
  # at 921 and 2001 rings respectively
  for (n in c(30, 60, 200)) {
    cv <- sd_bounds_curve(l = 1, u = 7, n = n, Z = "quasiinteger")
    expect_equal(n_rings(sd_bounds_plot(cv)), 1L,
                 info = sprintf("l = 1, u = 7, n = %d", n))
  }
  # and on a wider scale, where the threshold is crossed sooner
  cv <- sd_bounds_curve(l = 0, u = 10, n = 80, Z = "quasiinteger")
  expect_equal(n_rings(sd_bounds_plot(cv)), 1L)
})

test_that("genuine gaps in the band still split it into separate rings", {
  # strictly integer data is feasible only at GRIM-consistent means, so the
  # band is a comb: merging those into one ring would shade the impossible
  # means between the teeth as feasible, the error band_polygon() exists to
  # prevent. The fix must not over-merge.
  cv <- sd_bounds_curve(l = 1, u = 7, n = 30, Z = "integer")
  expect_gt(n_rings(sd_bounds_plot(cv)), 50L)
})

test_that("a hand-built uniform curve still works without the attribute", {
  cv <- sd_bounds_curve(l = 1, u = 7, n = 60, Z = "quasiinteger")
  flat <- data.frame(mean = seq(1, 7, by = 0.006))
  bd <- do.call(rbind, lapply(flat$mean, function(m)
    sd_bounds(l = 1, u = 7, n = 60, mean = m, Z = "quasiinteger")))
  flat$min_sd <- bd$min_sd
  flat$max_sd <- bd$max_sd
  flat$feasible <- bd$feasible
  expect_null(attr(flat, "step"))
  expect_equal(n_rings(sd_bounds_plot(flat)), 1L)
})
