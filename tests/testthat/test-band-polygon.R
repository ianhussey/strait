test_that("band_polygon() splits a gapped band into separate rings", {
  # the alpha rule leaves stretches near each limit where no composite with the
  # reported alpha exists; sd_region_data() returns NA there
  d <- sd_region_data(0, 3, 7, rule = "alpha", n_items = 2, alpha = 0.70,
                      by = 0.001)
  expect_equal(nrow(d), 3001L)
  expect_equal(sum(is.na(d$min_sd)), 174L)

  rings <- band_polygon(d, by = 0.001)
  expect_equal(length(unique(rings$ring)), 3L)
  expect_named(rings, c("mean", "y", "ring"))
  # each ring is closed: first and last vertex coincide
  for (k in unique(rings$ring)) {
    g <- rings[rings$ring == k, ]
    expect_equal(g$mean[1], g$mean[nrow(g)])
    expect_equal(g$y[1], g$y[nrow(g)])
  }
})

test_that("band_polygon() leaves an ungapped band as one ring", {
  d <- sd_region_data(1, 5, 7, rule = "quasi", by = 0.01)
  rings <- band_polygon(d, by = 0.01)
  expect_equal(length(unique(rings$ring)), 1L)
  # a ring traces the ceiling out and the floor back
  expect_equal(nrow(rings), 2 * nrow(d) + 1)
})

test_that("band_polygon() handles degenerate input", {
  expect_null(band_polygon(data.frame(mean = 1, min_sd = NA, max_sd = NA), by = 0.1))
  expect_null(band_polygon(data.frame(mean = numeric(0), min_sd = numeric(0),
                                      max_sd = numeric(0)), by = 0.1))
})

test_that("sd_delta() is exactly the gap between the two ceilings", {
  # sd_max_structure_s^2 == muilwijk^2 * delta, two independent implementations
  for (n in c(5, 7, 10, 30, 101)) {
    for (lu in list(c(1, 5), c(1, 7), c(0, 6))) {
      m <- seq(lu[1], lu[2], length.out = 401)
      expect_equal(sd_max_structure_s(m, n, lu[1], lu[2])^2,
                   sd_max_muilwijk(m, n, lu[1], lu[2])^2 *
                     sd_delta(m, n, lu[1], lu[2]),
                   tolerance = 1e-9)
    }
  }
})

test_that("sd_delta() is 1 at whole counts and falls toward the limits", {
  # midpoint of a 1-5 scale at n = 7: n_l = n_u = 3.5, delta = 6/7
  expect_equal(sd_delta(3, 7, 1, 5), 6 / 7, tolerance = 1e-9)
  # whole counts give exactly 1 (no epsilon shortfall)
  expect_equal(sd_delta(3, 6, 1, 5), 1)
  # barely attainable near a limit, and never outside [0, 1]
  expect_lt(sd_delta(1.07, 7, 1, 5), 0.2)
  d <- sd_delta(seq(1, 5, by = 0.01), 7, 1, 5)
  expect_true(all(d >= 0 & d <= 1))
})

test_that("the muilwijk alias is the mean rule", {
  expect_identical(sd_region_data(1, 7, 7, rule = "muilwijk"),
                   sd_region_data(1, 7, 7, rule = "mean"))
  # and it is genuinely the uncorrected form: mestdagh never exceeds it
  a <- sd_region_data(1, 7, 7, rule = "mean")
  b <- sd_region_data(1, 7, 7, rule = "mestdagh")
  expect_true(all(b$max_sd <= a$max_sd + 1e-9))
  expect_equal(round(max(a$max_sd - b$max_sd), 2), 0.64)
})

test_that("plot_umbrella() supports both styles", {
  um <- suppressWarnings(umbrella_data(n = 12, l = 1, u = 3, digits = 1))
  cur <- sd_bounds_curve(l = 1, u = 3, n = 12, by = 0.1)
  expect_s3_class(plot_umbrella(um), "ggplot")
  expect_s3_class(plot_umbrella(um, style = "tiles"), "ggplot")
  expect_s3_class(plot_umbrella(um, bounds = cur), "ggplot")
  expect_error(plot_umbrella(um, style = "blobs"))
  # points style accepts an already-filtered lattice with no `consistent` column
  lat <- sd_region_data(1, 3, 12, rule = "integer", digits = 1)
  expect_false("consistent" %in% names(lat))
  expect_s3_class(plot_umbrella(lat), "ggplot")
})

test_that("plot_sd_region() supports both shading conventions", {
  quasi <- sd_region_data(1, 5, 7, rule = "quasi")
  expect_s3_class(plot_sd_region(quasi), "ggplot")
  expect_s3_class(plot_sd_region(quasi, shade = "inside"), "ggplot")
  # the gapped rule is the one shade = "outside" must not misdraw
  expect_s3_class(plot_sd_region(sd_region_data(0, 3, 7, rule = "alpha",
                                                n_items = 2, alpha = 0.70)),
                  "ggplot")
})
