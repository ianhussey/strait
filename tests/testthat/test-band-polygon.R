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
  empty <- tibble::tibble(mean = numeric(0), y = numeric(0), ring = numeric(0))
  expect_identical(band_polygon(data.frame(mean = 1, min_sd = NA, max_sd = NA), by = 0.1), empty)
  expect_identical(band_polygon(data.frame(mean = numeric(0), min_sd = numeric(0),
                                           max_sd = numeric(0)), by = 0.1), empty)
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

test_that("umbrella_plot() supports both styles", {
  um <- suppressWarnings(umbrella_data(n = 12, l = 1, u = 3, digits = 1))
  cur <- sd_bounds_curve(l = 1, u = 3, n = 12, by = 0.1)
  expect_s3_class(umbrella_plot(um), "ggplot")
  expect_s3_class(umbrella_plot(um, style = "tiles"), "ggplot")
  expect_s3_class(umbrella_plot(um, bounds = cur), "ggplot")
  expect_error(umbrella_plot(um, style = "blobs"))
  # points style accepts an already-filtered lattice with no `consistent` column
  lat <- sd_region_data(1, 3, 12, rule = "integer", digits = 1)
  expect_false("consistent" %in% names(lat))
  expect_s3_class(umbrella_plot(lat), "ggplot")
})

test_that("sd_region_plot() supports both shading conventions", {
  quasi <- sd_region_data(1, 5, 7, rule = "quasi")
  expect_s3_class(sd_region_plot(quasi), "ggplot")
  expect_s3_class(sd_region_plot(quasi, shade = "inside"), "ggplot")
  # the gapped rule is the one shade = "outside" must not misdraw
  expect_s3_class(sd_region_plot(sd_region_data(0, 3, 7, rule = "alpha",
                                                n_items = 2, alpha = 0.70)),
                  "ggplot")
})

test_that("umbrella_contour() is the envelope of the consistent tuples", {
  um <- suppressWarnings(umbrella_data(l = 1, u = 3, n = 12, digits = 1))
  pts <- um[um$consistent, c("mean", "sd")]
  rings <- umbrella_contour(um)
  # one region, spanning the means that have a consistent tuple, and reaching
  # the extreme SDs those tuples take
  length(unique(rings$ring)) |> expect_equal(1L)
  range(rings$mean) |> expect_equal(range(pts$mean))
  range(rings$y) |> expect_equal(range(pts$sd))
  # it filters by `consistent` itself, so an already-filtered lattice gives
  # the same rings
  umbrella_contour(pts) |> expect_equal(rings)
  nrow(umbrella_contour(um[0, ])) |> expect_equal(0L)
})

test_that("the contour joins GRIM's stripes but splits a genuine gap", {
  # at two decimals only the means an integer sum rounds to survive, so the
  # reporting grid's own step would make every stripe a zero-width ring
  striped <- suppressWarnings(umbrella_data(l = 1, u = 5, n = 25, digits = 2))
  pts <- striped[striped$consistent, c("mean", "sd")]
  expect_gt(min(diff(sort(unique(pts$mean)))), 0.01)
  expect_equal(length(unique(umbrella_contour(pts)$ring)), 1L)
  # forcing the reporting step back in splits it into one ring per stripe
  by_grid <- umbrella_contour(pts, by = 0.01)
  expect_equal(length(unique(by_grid$ring)), length(unique(pts$mean)))
  # a reported alpha leaves stretches of infeasible means, which stay split
  gapped <- suppressWarnings(umbrella_data(
    l = 0, u = 3, n = 7, digits = 2, scoring = "meanscored", n_items = 2,
    alpha = 0.70
  ))
  gapped <- gapped[gapped$consistent, c("mean", "sd")]
  expect_gt(length(unique(umbrella_contour(gapped)$ring)), 1L)
})

test_that("umbrella_plot() draws the contour, shaded or not", {
  um <- suppressWarnings(umbrella_data(l = 1, u = 3, n = 12, digits = 1))
  umbrella_plot(um, style = "contour") |> expect_s3_class("ggplot")
  # the shading is a layer, so "none" has one layer fewer
  n_layers <- function(p) length(p$layers)
  expect_equal(
    n_layers(umbrella_plot(um, shade = "none")),
    n_layers(umbrella_plot(um)) - 1L
  )
  # and the ring is left unfilled, so what is underneath shows through
  contour_none <- umbrella_plot(um, style = "contour", shade = "none")
  expect_true(is.na(contour_none$layers[[1]]$aes_params$fill))
  expect_identical(
    umbrella_plot(um, style = "contour")$layers[[2]]$aes_params$fill,
    "white"
  )
  expect_error(umbrella_plot(sd_bounds_curve(1, 3, 12), style = "contour"))
})
