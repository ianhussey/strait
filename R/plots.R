# ggplot2 visualisations of the SD bounds, POMP-normalised bounds, and umbrella.

.bounds_point_layer <- function(pts, xvar, yvar) {
  if (is.null(pts[["consistent"]])) {
    pts$consistent <- NA
  }
  pts$.consistent <- as.character(pts$consistent)
  list(
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(
        x = .data[[xvar]],
        y = .data[[yvar]],
        fill = .data$.consistent
      ),
      shape = 21,
      colour = "black",
      size = 2.4,
      na.rm = TRUE
    ),
    ggplot2::scale_fill_manual(
      values = c("TRUE" = COLOUR_CONSISTENT, "FALSE" = COLOUR_INCONSISTENT),
      na.value = "grey50",
      name = "Consistent",
      labels = c("TRUE" = "consistent", "FALSE" = "inconsistent")
    )
  )
}

#' Plot SD bounds on the native scale
#'
#' The feasible SD band (floor to ceiling) against the mean, optionally with
#' reported points coloured by consistency (green/red outlined dots).
#'
#' `shade = "outside"` (default) shades the infeasible region and leaves the
#' feasible one clear, matching [sd_bounds_pomp_plot()] with
#' `reference = "sharp"`, so the two scales read the same way: shaded means
#' unreachable. `shade = "inside"` fills the feasible band instead, which was
#' the behaviour before this argument existed.
#'
#' @param data Output of [sd_bounds_curve()], or a band from
#'   [sd_region_data()]. Its `"step"` attribute, when present, is the mean-grid
#'   spacing used to tell a sampling gap from a genuine one under `shade =
#'   "outside"`; a band built by hand, without that attribute, should be on a
#'   uniform grid, from which the spacing is inferred.
#' @param reports Optional data frame of reported values, with `mean`, `sd`,
#'   and (optionally) `consistent`; e.g. the output of [brimmer_map()].
#' @param title Optional plot title.
#' @param fill,line_colour Band fill and outline colours. `fill` is used only
#'   by `shade = "inside"`; the infeasible shading has its own fixed grey.
#' @param shade `"outside"` (default) shades the infeasible region;
#'   `"inside"` fills the feasible band.
#' @param expand Padding around the plotted region, as a proportion of the
#'   scale width `u - l` rather than a fixed number of SD units, so that the
#'   margin looks the same on a 1-5 scale and a 0-100 one. Applied to both
#'   axes. The limits always stretch to include `reports`, so an out-of-bounds
#'   report is never clipped out of view.
#' @return A ggplot object.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#' sd_bounds_plot(curve, title = "Feasible SDs, 1-7 scale, n = 30")
#'
#' # overlay reported values, coloured by consistency (the second is
#' # above the ceiling, so it plots as inconsistent)
#' reports <- data.frame(mean = c(2.97, 3.51), sd = c(2.83, 3.50))
#' checked <- brimmer_map(reports, l = 1, u = 7, n = 30,
#'                             digits_mean = 2, digits_sd = 2)
#' sd_bounds_plot(curve, reports = checked)
#'
#' # the previous look, with the feasible band filled
#' sd_bounds_plot(curve, reports = checked, shade = "inside")
#' @export
sd_bounds_plot <- function(
  data,
  reports = NULL,
  title = NULL,
  fill = "grey85",
  line_colour = "grey30",
  shade = c("outside", "inside"),
  expand = 0.03
) {
  shade <- match.arg(shade)
  cur <- .band_mask(data)

  # Padding is a proportion of the scale width, not a fixed SD amount, so the
  # margin is visually constant across scales of very different widths.
  drawn <- cur$mean[!is.na(cur$max_sd)]
  lo_m <- min(drawn)
  hi_m <- max(drawn)
  pad <- expand * (hi_m - lo_m)
  # Shading the outside needs finite limits, and finite limits will silently
  # clip an out-of-bounds point - exactly the case the plot exists to show -
  # so the ceiling of the view must account for the reported points as well.
  y_hi <- max(c(cur$max_sd, reports$sd), na.rm = TRUE)

  p <- ggplot2::ggplot(cur, ggplot2::aes(x = .data$mean))
  if (shade == "inside") {
    p <- p +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$min_sd, ymax = .data$max_sd),
        fill = fill
      )
  } else {
    # Shade the whole panel, then knock the feasible rings out of it. Doing it
    # this way rather than assembling the grey from side rectangles and ribbons
    # means a mean at which the band is undefined stays shaded by construction:
    # there is nothing there to knock out. Assembling the grey instead leaves
    # such a gap unshaded, which would assert that any SD at all is possible
    # there. See band_polygon().
    # band_polygon() needs the grid spacing to tell a sampling gap from a
    # genuine one; see .band_step().
    rings <- band_polygon(cur, by = .band_step(data))
    p <- p +
      ggplot2::annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = "grey10",
        alpha = 0.12
      )
    if (!is.null(rings)) {
      p <- p +
        ggplot2::geom_polygon(
          data = rings,
          ggplot2::aes(x = .data$mean, y = .data$y, group = .data$ring),
          inherit.aes = FALSE,
          fill = "white"
        )
    }
  }
  p <- p +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$max_sd),
      colour = line_colour,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$min_sd),
      colour = line_colour,
      na.rm = TRUE
    ) +
    ggplot2::coord_cartesian(
      xlim = c(lo_m - pad, hi_m + pad),
      ylim = c(-pad, y_hi + pad),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Mean", y = "SD", title = title) +
    ggplot2::theme_minimal()
  if (!is.null(reports)) {
    p <- p + .bounds_point_layer(reports, "mean", "sd")
  }
  p
}

#' Plot SD bounds on a percent-of-maximum-possible (POMP) scale
#'
#' `reference = "parity"` normalises every SD by the mean-agnostic parity
#' (Popoviciu) ceiling: a linear rescaling, so the Structure S ceiling appears
#' as a dome under 1 and the umbrella geometry is undistorted. `reference =
#' "sharp"` normalises each SD by its own sharp mean-conditional band, so the
#' feasible region is exactly the unit square and a point's height is its
#' position within the band (`pomp_sd_sharp`); regions outside `[0, 1]^2` are
#' shaded infeasible.
#'
#' @param data Output of [sd_bounds_curve()] (used for the parity band).
#' @param reports Optional data frame of reported values, with `pomp_mean` and
#'   `pomp_sd_parity` / `pomp_sd_sharp` and `consistent` (e.g. from
#'   [brimmer_map()]).
#' @param reference `"sharp"` (default) or `"parity"`.
#' @param title Optional plot title.
#' @return A ggplot object.
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#'
#' # "sharp" makes the feasible region exactly the unit square
#' sd_bounds_pomp_plot(curve)
#'
#' # "parity" is a linear rescaling, so the umbrella keeps its shape
#' sd_bounds_pomp_plot(curve, reference = "parity")
#' @export
sd_bounds_pomp_plot <- function(
  data,
  reports = NULL,
  reference = c("sharp", "parity"),
  title = NULL
) {
  reference <- match.arg(reference)
  if (reference == "parity") {
    cur <- .band_mask(data)
    p <- ggplot2::ggplot(cur, ggplot2::aes(x = .data$pomp_mean)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$floor_parity, ymax = .data$ceil_parity),
        fill = "grey85"
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$ceil_parity),
        colour = "grey30",
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        ggplot2::aes(y = .data$floor_parity),
        colour = "grey30",
        na.rm = TRUE
      ) +
      ggplot2::labs(
        x = "Relative location (POMP mean)",
        y = "Relative dispersion (parity-normalised SD)",
        title = title
      )
    if (!is.null(reports)) {
      p <- p + .bounds_point_layer(reports, "pomp_mean", "pomp_sd_parity")
    }
  } else {
    shade <- data.frame(
      xmin = c(-Inf, 1, 0, 0),
      xmax = c(0, Inf, 1, 1),
      ymin = c(-Inf, -Inf, 1, -Inf),
      ymax = c(Inf, Inf, Inf, 0)
    )
    p <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = shade,
        ggplot2::aes(
          xmin = .data$xmin,
          xmax = .data$xmax,
          ymin = .data$ymin,
          ymax = .data$ymax
        ),
        fill = "grey10",
        alpha = 0.12
      ) +
      ggplot2::geom_rect(
        data = data.frame(x = 0),
        ggplot2::aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1),
        fill = NA,
        colour = "black",
        linewidth = 0.3
      ) +
      ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1)) +
      ggplot2::labs(
        x = "Relative location (POMP mean)",
        y = "Position in sharp SD band",
        title = title
      )
    if (!is.null(reports)) {
      p <- p + .bounds_point_layer(reports, "pomp_mean", "pomp_sd_sharp")
    }
  }
  p + ggplot2::theme_minimal()
}

#' Plot the umbrella
#'
#' Draws the umbrella of a design: its reportable `(mean, sd)` tuples, or the
#' continuous region they are drawn from. What is drawn follows from `data`.
#'
#' A lattice of tuples, such as [umbrella_data()] output, is drawn in one of
#' three styles. `style = "points"` (default) greys the whole panel and draws only the
#' consistent tuples. Nothing else is drawn, because nothing else exists: every
#' other cell of the grid is a value that cannot be reported. This matches the
#' convention of [sd_bounds_plot()] and [sd_region_plot()], where shading marks
#' what is ruled out, and it makes the real structure legible — the vertical
#' striping at means an integer sum can round to.
#'
#' `style = "tiles"` draws every cell of the reporting grid, coloured as
#' consistent, GRIMMER-inconsistent, or out of bounds, so it needs the full
#' [umbrella_data()] output. That is useful for methods exposition, since it
#' separates what the bounds rule out from what GRIMMER additionally rules out,
#' but it inverts the emphasis: at `n = 14` on a 1-7 scale at two decimal places
#' it spends two thirds of its ink, in the most saturated colour, on impossible
#' tuples, and GRIMMER's verdict alternating between adjacent SDs produces
#' interference banding that obscures the striping.
#'
#' `style = "contour"` draws the outline of the consistent tuples instead of the
#' tuples themselves: at each mean, the lowest and highest consistent SD, from
#' [umbrella_contour()]. Where the striping of the lattice would only add noise,
#' this shows the shape of the umbrella.
#'
#' A band, with `min_sd` and `max_sd` instead of `sd`, is the umbrella itself,
#' filled in `fill`. [sd_bounds_curve()] gives the sharp one: under its default
#' `granularity = "quasiinteger"`, every `(mean, SD)` of `n` values in `[l, u]` that are
#' whole numbers, except possibly one, which lets the mean take any value. The
#' tuples of a lattice are this region rounded to the reporting grid and
#' thinned by GRIM and GRIMMER, so the region is their outline, free of the
#' grid. Its edges are smooth between kinks at multiples of `1 / n`, which
#' `sd_bounds_curve()` evaluates exactly.
#'
#' For an umbrella without the kinks, use the band of
#' `sd_region_data(l, u, n, rule = "mean_naive_floor")`. Both of its edges are
#' variances of two-point distributions, `(mean - a) * (b - mean)` times
#' `n / (n - 1)`: the ceiling takes `a` and `b` to be the scale limits (the
#' Bhatia-Davis inequality), and each arch of the floor takes them to be the
#' whole numbers on either side of the mean. They meet the sharp bounds
#' wherever these can be attained: the floor at every mean that integer data
#' can have, the ceiling at every multiple of `(u - l) / n` above `l`. For data
#' of only 0s and 1s, the two coincide in DEBIT's curve.
#'
#' Note what the tuples represent. The consistent set is the GRIM- and
#' GRIMMER-consistent one, which is what [brimmer()] applies. It is strictly
#' larger than the set of attainable tuples, so an umbrella plot shows what the
#' test admits rather than what exists; use [brimmest()] to certify a tuple.
#' The gap is not always small — 9\% at `l = 1, u = 5, n = 10`, and nearly a
#' factor of three at `l = 0, u = 6, n = 23` with two items.
#'
#' @param data What to draw: a lattice with `mean` and `sd`, such as
#'   [umbrella_data()] output or `sd_region_data(rule = "integer")`, filtered by
#'   `consistent` when that column is present; or a band with `mean`, `min_sd`,
#'   and `max_sd`, such as [sd_bounds_curve()] output.
#' @param bounds Optional band, such as [sd_bounds_curve()] output, drawn as
#'   dashed lines on top of `data`.
#' @param title Optional plot title.
#' @param style For a lattice, `"points"` (default), `"tiles"`, or
#'   `"contour"`; see Details.
#' @param fill Fill colour of a band.
#' @param point_colour,point_size Point appearance. `point_size` defaults to a
#'   value chosen from the number of points, since a size that reads well for a
#'   few hundred is a solid mass at twenty thousand.
#' @param line_colour Colour of the outline under `style = "contour"`.
#' @param contour_by Mean spacing that separates a true gap in the umbrella
#'   from its striping under `style = "contour"`; see [umbrella_contour()].
#' @param reference_colour Colour of the `bounds` lines.
#' @param shade `"outside"` (default) shades the panel where the umbrella is
#'   not, as [sd_bounds_plot()] does; `"none"` leaves it clear, so the plot can
#'   be layered over another.
#' @param expand Padding as a proportion of the scale width, as in
#'   [sd_bounds_plot()].
#' @return A ggplot object.
#' @examples
#' grid <- umbrella_data(l = 1, u = 3, n = 12, digits = 1)
#' umbrella_plot(grid, title = "n = 12, 1-3 scale")
#'
#' # the previous look, separating the two ways a tuple can be ruled out
#' umbrella_plot(grid, style = "tiles")
#'
#' # the outline of the consistent tuples
#' umbrella_plot(grid, style = "contour")
#'
#' # the continuous umbrella that the tuples are drawn from
#' umbrella_plot(sd_bounds_curve(l = 1, u = 3, n = 12))
#'
#' # the same without kinks: its edges are smooth arcs
#' umbrella_plot(sd_region_data(l = 1, u = 3, n = 12, rule = "mean_naive_floor"))
#' @export
umbrella_plot <- function(
  data,
  bounds = NULL,
  title = NULL,
  style = c("points", "tiles", "contour"),
  fill = "grey30",
  point_colour = "black",
  point_size = NULL,
  line_colour = "grey30",
  contour_by = NULL,
  reference_colour = "grey35",
  shade = c("outside", "none"),
  expand = 0.03
) {
  style <- match.arg(style)
  shade <- match.arg(shade)
  is_band <- !"sd" %in% names(data)
  if (is_band && !all(c("min_sd", "max_sd") %in% names(data))) {
    stop(
      "`data` needs `mean` and `sd` (a lattice) or `mean`, `min_sd`, and `max_sd` (a band)"
    )
  }
  if (is_band && style != "points") {
    stop(
      "style = '",
      style,
      "' needs a lattice from umbrella_data(), not a band"
    )
  }

  bnd <- if (!is.null(bounds)) .band_mask(bounds)
  bounds_layers <- if (!is.null(bnd)) {
    lapply(c("max_sd", "min_sd"), function(col) {
      ggplot2::geom_line(
        data = bnd,
        ggplot2::aes(.data$mean, .data[[col]]),
        inherit.aes = FALSE,
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.35,
        na.rm = TRUE
      )
    })
  }

  if (style == "tiles") {
    if (!all(c("consistent", "in_bounds", "grimmer") %in% names(data))) {
      stop(
        "style = 'tiles' needs the full umbrella_data() output, with `consistent`, `in_bounds`, and `grimmer`"
      )
    }
    data$category <- ifelse(
      data$consistent,
      "consistent",
      ifelse(
        data$in_bounds & !is.na(data$grimmer) & !data$grimmer,
        "GRIMMER-inconsistent",
        "out of bounds"
      )
    )
    return(
      ggplot2::ggplot(
        data,
        ggplot2::aes(x = .data$mean, y = .data$sd, fill = .data$category)
      ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_manual(
          values = c(
            "consistent" = COLOUR_CONSISTENT,
            "GRIMMER-inconsistent" = COLOUR_INCONSISTENT,
            "out of bounds" = "grey80"
          ),
          name = NULL
        ) +
        bounds_layers +
        ggplot2::labs(x = "Mean", y = "SD", title = title) +
        ggplot2::theme_minimal()
    )
  }

  if (is_band) {
    cur <- .band_mask(data)
    region <- ggplot2::geom_ribbon(
      data = cur,
      ggplot2::aes(x = .data$mean, ymin = .data$min_sd, ymax = .data$max_sd),
      fill = fill
    )
    means <- range(cur$mean[!is.na(cur$max_sd)])
    y_hi <- max(c(cur$max_sd, bnd$max_sd), na.rm = TRUE)
  } else {
    pts <- if ("consistent" %in% names(data)) {
      data[
        !is.na(data$consistent) & data$consistent,
        c("mean", "sd"),
        drop = FALSE
      ]
    } else {
      data[, c("mean", "sd"), drop = FALSE]
    }
    if (style == "contour") {
      # rings knocked out of the shaded panel rather than a ribbon, for the
      # reason band_polygon() exists: where the umbrella genuinely stops, the
      # shading must stay
      region <- ggplot2::geom_polygon(
        data = umbrella_contour(pts, by = contour_by),
        ggplot2::aes(x = .data$mean, y = .data$y, group = .data$ring),
        fill = if (shade == "outside") "white" else NA,
        colour = line_colour,
        linewidth = 0.4
      )
    } else {
      if (is.null(point_size)) {
        point_size <- if (nrow(pts) > 12000) {
          0.045
        } else if (nrow(pts) > 2000) {
          0.12
        } else {
          0.35
        }
      }
      region <- ggplot2::geom_point(
        data = pts,
        ggplot2::aes(.data$mean, .data$sd),
        colour = point_colour,
        size = point_size,
        shape = 16,
        na.rm = TRUE
      )
    }
    means <- range(data$mean)
    y_hi <- max(c(pts$sd, bnd$max_sd), na.rm = TRUE)
  }
  pad <- expand * diff(means)

  # the panel is infeasible everywhere except where the umbrella is drawn
  shading <- if (shade == "outside") {
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey10",
      alpha = 0.12
    )
  }

  ggplot2::ggplot() +
    shading +
    region +
    bounds_layers +
    ggplot2::coord_cartesian(
      xlim = means + c(-pad, pad),
      ylim = c(-pad, y_hi + pad),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Mean", y = "Sample standard deviation", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
