# One plotting entry point for the feasible-region panels: each `rule` is a
# constraint set from the nested framework, drawn either as a continuous band or
# (for strictly integer data) as the lattice of attainable reported tuples.

# Internal: the alpha-free quasi-integer sharp band. `mg` is the granularity
# multiplier: the affine map w = mg * (x - l) sends the reported grid to the
# integers, so mg = n_items for mean-scored composites and 1 otherwise.
quasi_band <- function(mean, n, l, u, mg) {
  tibble::tibble(
    mean = mean,
    min_sd = sd_min_quasi_integer(mg * (mean - l), n) / mg,
    max_sd = sd_max_structure_s(mg * (mean - l), n, 0, mg * (u - l)) / mg
  )
}

# Internal: resolve a `scoring` choice into the three quantities the region
# formulas need. `mg` is the granularity multiplier (above); `to_sum` maps a
# reported mean to the sum-score mean that the alpha bounds are stated in;
# `sd_div` converts a sum-score SD back to the reported SD's units.
scoring_geometry <- function(scoring, k, l, u) {
  switch(
    scoring,
    singleitem = {
      if (k != 1L) {
        cli::cli_abort(
          "{.code scoring = \"singleitem\"} requires {.code n_items = 1}"
        )
      }
      list(mg = 1, item_l = l, item_u = u, to_sum = function(m) m, sd_div = 1)
    },
    meanscored = list(
      mg = k,
      item_l = l,
      item_u = u,
      to_sum = function(m) k * m,
      sd_div = k
    ),
    sumscored = {
      if (
        abs(l / k - round(l / k)) > TOLERANCE ||
          abs(u / k - round(u / k)) > TOLERANCE
      ) {
        cli::cli_abort(
          "{.code scoring = \"sumscored\"} needs {.arg l} and {.arg u} divisible by {.arg n_items} (they are the composite's limits, k times the per-item limits)"
        )
      }
      list(
        mg = 1,
        item_l = l / k,
        item_u = u / k,
        to_sum = function(m) m,
        sd_div = 1
      )
    }
  )
}

gcd <- function(a, b) {
  while (b) {
    remainder <- a %% b
    a <- b
    b <- remainder
  }
  a
}

# Internal: round a lattice of (mean, sd) pairs to reporting precision and
# collapse the duplicates that creates. Distinct attainable pairs routinely
# round into one reported cell, which is exactly why a rounded lattice can look
# solid where the exact one is full of holes.
round_lattice <- function(d, digits, rounding) {
  if (is.null(digits)) {
    return(d)
  }
  d$mean <- roundwork::reround(d$mean, digits, rounding)
  d$sd <- roundwork::reround(d$sd, digits, rounding)
  # some rounding rules return negative zero at 0, which compares equal but
  # prints and string-matches differently; normalise it away
  d$mean[d$mean == 0] <- 0
  d$sd[d$sd == 0] <- 0
  d <- unique(d)
  d <- d[order(d$mean, d$sd), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# Internal: the EXACT attainable (mean, SD) pairs for strictly integer data,
# with no reporting grid. A dynamic program over the n observations accumulates
# the reachable (sum, sum-of-squares) pairs on the shifted grid
# y = mg * (x - l) in 0..W, from which both statistics are recovered. Unlike the
# rounded `"integer"` lattice this shows the true interior holes, which rounding
# smears shut.
#
# Five things keep the state space small enough for this to stay fast in R:
#
#  (1) The second axis is NOT the sum of squares Q but R = W*S - Q = sum
#      y(W - y). Every term is non-negative and at most floor(W^2/4), so R
#      spans n*floor(W^2/4) rather than n*W^2: about a quarter of the cells.
#  (2) R is then divided by g, the gcd of the achievable y(W - y). For odd W
#      every y(W - y) is even, so g = 2 halves the axis again.
#  (3) Only the frontier reachable after t items is live at step t, so the
#      active window grows instead of the full grid being swept n times.
#  (4) The reachable set is symmetric under y -> W - y, which maps S -> tW - S
#      and leaves R fixed. Only the lower half of the S axis is computed; the
#      upper half is a reversed copy.
#  (5) Cells are `raw` (one byte) rather than `logical` (four).
#
# Shifted whole-block ORs are used rather than scattering into which(...)
# indices: the block form is what lets (3) and (4) restrict the work.
attainable_lattice <- function(
  l,
  u,
  n,
  mg = 1,
  max_cells = MAX_LATTICE_CELLS
) {
  W <- mg * (u - l)
  if (abs(W - round(W)) > TOLERANCE) {
    cli::cli_abort(
      "{.arg l}, {.arg u} and {.arg n_items} must put the scale limits on the integer grid"
    )
  }
  W <- as.integer(round(W))
  if (W < 1L) {
    cli::cli_abort("the scale limits must span at least one grid step")
  }
  ys <- 0:W
  pr <- ys * (W - ys) # each item's contribution to R
  g <- Reduce(gcd, pr[pr > 0])
  if (!length(g) || is.na(g) || g < 1) {
    g <- 1
  }
  pm <- max(pr) / g # per-item span of the R axis
  Smax <- n * W
  Rmax <- n * pm
  if ((Smax + 1) * (Rmax + 1) > max_cells) {
    cli::cli_abort(c(
      "the exact lattice is too large to enumerate here ({format(Smax + 1)} x {format(Rmax + 1)} cells)",
      "i" = "Use {.code rule = \"integer\"} with a reporting precision instead."
    ))
  }

  z <- as.raw(0)
  A <- matrix(z, Smax + 1L, Rmax + 1L)
  B <- matrix(z, Smax + 1L, Rmax + 1L)
  A[1L, 1L] <- as.raw(1)
  drs <- pr / g
  for (t in seq_len(n)) {
    Sh <- (t - 1L) * W
    Rh <- (t - 1L) * pm # frontier before this item
    So <- t * W
    Ro <- t * pm # frontier after it
    hl <- floor(So / 2) + 1L # rows computed; rest mirrored
    B[1L:(So + 1L), 1L:(Ro + 1L)] <- z
    for (i in seq_along(ys)) {
      v <- ys[i]
      dr <- drs[i]
      r_hi <- min(Sh + 1L + v, hl)
      if (r_hi < 1L + v) {
        next
      }
      B[(1L + v):r_hi, (1L + dr):(Rh + 1L + dr)] <-
        B[(1L + v):r_hi, (1L + dr):(Rh + 1L + dr)] |
        A[1L:(r_hi - v), 1L:(Rh + 1L), drop = FALSE]
    }
    if (hl < So + 1L) {
      # y -> W - y symmetry
      B[(hl + 1L):(So + 1L), 1L:(Ro + 1L)] <-
        B[(So + 1L - hl):1L, 1L:(Ro + 1L), drop = FALSE]
    }
    swap <- A
    A <- B
    B <- swap
  }

  live <- which(A != z) - 1L
  n_rows <- Smax + 1L
  S <- live %% n_rows
  Q <- W * S - (live %/% n_rows) * g
  # n times the sum of squares, a whole number, so that brimmest() can round
  # the SD exactly; ss = nss / n is translation-free
  nss <- n * Q - S^2
  means <- l + S / (n * mg)
  sds <- sqrt(pmax(0, nss) / (n * (n - 1))) / mg
  # new_tibble() rather than tibble(): the lattice can run to millions of rows
  ord <- order(means, sds)
  tibble::new_tibble(
    list(mean = means[ord], sd = sds[ord], S = S[ord], nss = nss[ord]),
    nrow = length(ord)
  )
}

# Internal: give a region the same columns whatever its rule, so that
# sd_region_data() is type-stable: a band leaves `sd` NA, a lattice
# `min_sd` and `max_sd`.
region_columns <- function(region, type) {
  template <- c("mean", "min_sd", "max_sd", "sd")
  for (column in setdiff(template, names(region))) {
    region[[column]] <- rep(NA_real_, nrow(region))
  }
  out <- tibble::as_tibble(region[template])
  attr(out, "type") <- type
  out
}

# Internal: is `data` a band (`min_sd`, `max_sd`) rather than a lattice of
# tuples (`sd`)? sd_region_data() has both columns and says which in its
# `type` attribute; sd_bounds_curve() output has no `sd` at all.
is_band <- function(data) {
  sd <- data[["sd"]]
  is.null(sd) ||
    identical(attr(data, "type"), "band") ||
    (length(sd) > 0 && all(is.na(sd)))
}

#' Feasible-region data for one constraint set
#'
#' The `(mean, min_sd, max_sd)` band, or the lattice of attainable reported
#' `(mean, sd)` tuples, for a chosen constraint set. [sd_region_plot()] draws
#' it, and [umbrella_plot()] draws the bands of the umbrella-shaped rules.
#'
#' Rules, in the order they enter the framework (`n_items > 1` puts every rule in
#' mean-score units for a composite of that many integer items):
#'
#' \describe{
#'   \item{`"range"`}{`(u - l)/sqrt(2)`, no floor; ignores `n` and the mean
#'     (Popoviciu 1935, range only).}
#'   \item{`"range_n"`}{the parity ceiling, no floor (Popoviciu, as restored by
#'     Petocz 2005).}
#'   \item{`"mean"`, `"muilwijk"`}{Muilwijk's mean-conditional ceiling as
#'     originally stated, over a floor of zero (Muilwijk 1966; Bhatia-Davis
#'     2000). A valid ceiling at every mean, but attained only where the counts
#'     at each limit come out whole, so it is conservative elsewhere: see
#'     `"mestdagh"` for the same constraint set with the count-parity
#'     correction applied, and [sd_delta()] for the correction itself. At
#'     `n = 7` on a 1-7 scale the two differ by up to 0.64 SD units. The two
#'     names are aliases.}
#'   \item{`"mean_naive_floor"`}{that arch plus the naive Bernoulli floor
#'     (Fuenderich et al. 2025).}
#'   \item{`"pesant_regin"`}{integer minimum with a loose (parity) ceiling
#'     (Pesant and Regin 2005).}
#'   \item{`"mestdagh"`}{the sharp integer maximum, no floor (Mestdagh et al.
#'     2018): `"mean"` with the count-parity correction, hence never above it.}
#'   \item{`"quasi"`}{both bounds sharp and GRIM-free (this package's default).
#'     This is the umbrella, which [umbrella_plot()] fills when given this
#'     band.}
#'   \item{`"alpha"`}{additionally conditioning on a reported Cronbach's `alpha`;
#'     needs `n_items >= 2`.}
#'   \item{`"integer"`}{strictly integer data: the lattice of reported
#'     `(mean, sd)` tuples passing GRIM, GRIMMER and the bounds — that is,
#'     exactly what [brimmer()] admits. Built via [umbrella_data()], so the
#'     pairs are rounded to `digits`. A set of points, because no band is
#'     defined at means and SDs integer data cannot produce.}
#'   \item{`"integer_alpha"`}{that lattice, additionally inside the
#'     alpha-conditional bounds.}
#'   \item{`"attainable"`}{the EXACT attainable `(mean, sd)` tuples of strictly
#'     integer data, with no reporting grid and so no `digits`. Unlike
#'     `"integer"` this shows the true interior holes, which rounding to a
#'     reporting grid smears shut; it is enumerated by dynamic programming and
#'     errors if the lattice is too large.}
#'   \item{`"attainable_alpha"`}{those exact tuples, additionally inside the
#'     alpha-conditional bounds.}
#' }
#'
#' The two lattice rules are easy to confuse and the difference is not small.
#' GRIMMER is necessary but not sufficient for an integer sample to exist, so
#' the attainable set is a strict subset of the GRIMMER-consistent one. At
#' `l = 1, u = 5, n = 10, digits = 2` there are 491 GRIMMER-consistent reported
#' pairs against 447 attainable ones rounded to the same grid: 9 percent pass
#' `"integer"` without being attainable, and none go the other way. At
#' `l = 0, u = 6, n = 23` with two items it is 49,030 against 47,871. Use
#' [brimmest()] to certify a single tuple rather than rounding a lattice to
#' compare against.
#'
#' Note also that `"integer"` returns pairs already rounded to `digits`, so
#' testing membership means rounding the candidate first: an exact sample SD of
#' 0.5270463 matches nothing until it becomes 0.53.
#'
#' @param l,u Numeric scalars, the scale limits (mean-score units when
#'   `n_items > 1` and `scoring` is left at its default).
#' @param n Integer scalar, sample size.
#' @param rule Which constraint set; see Details.
#' @param scoring One of `"singleitem"`, `"sumscored"`, `"meanscored"`, as in
#'   [sd_bounds()]. Defaults to `"meanscored"` when `n_items > 1` and
#'   `"singleitem"` otherwise. Under `"sumscored"`, `l` and `u` are the
#'   composite's limits and the reported values live on the unit grid, so
#'   `n_items` enters only as Cronbach's alpha's `k`.
#' @param n_items Integer, number of items in the composite (default 1).
#' @param alpha Reported Cronbach's alpha; required by the alpha rules.
#' @param digits Reported decimal places, used by the rounded lattice rules
#'   (default 2); ignored by the exact `"attainable"` rules.
#' @param round_digits Optionally round the returned `(mean, sd)` lattice to
#'   this many decimal places, collapsing pairs that round together. `NULL`
#'   (default) returns them unrounded, which for the `"attainable"` rules means
#'   exactly. Band rules are unaffected: rounding a bound can turn it
#'   anti-conservative, so bounds are always returned at full precision.
#' @param rounding How to round when `round_digits` is given: one of
#'   scrutiny's single rounding methods, such as `"up"` (default), `"down"`,
#'   `"even"`, `"ceiling"`, `"floor"`, `"trunc"`, or `"anti_trunc"`, applied
#'   by [roundwork::reround()].
#' @param by Optional mean-grid spacing (default `(u - l) / 1000`); ignored by
#'   the lattice rules, which step by `10^-digits`.
#' @return A tibble with `mean`, `min_sd`, `max_sd` and `sd`, and a `type`
#'   attribute of `"band"` or `"points"`. Band rules fill `min_sd` and
#'   `max_sd`, which are `NA` where no SD is possible, leave `sd` `NA`, and
#'   record the grid spacing `by` in a `"step"` attribute. The lattice rules
#'   (`"integer"`, `"integer_alpha"`, `"attainable"`, `"attainable_alpha"`)
#'   fill `sd` with the tuples passing every applicable test and leave
#'   `min_sd` and `max_sd` `NA`.
#' @examples
#' # a band rule returns (mean, min_sd, max_sd) across the scale
#' band <- sd_region_data(l = 1, u = 5, n = 7, rule = "quasi")
#' head(band)
#' attr(band, "type")
#'
#' # a lattice rule returns the attainable (mean, sd) tuples themselves
#' pts <- sd_region_data(l = 1, u = 5, n = 7, rule = "attainable")
#' head(pts)
#' attr(pts, "type")
#' @export
sd_region_data <- function(
  l,
  u,
  n,
  rule = c(
    "quasi",
    "range",
    "range_n",
    "mean",
    "muilwijk",
    "mean_naive_floor",
    "mestdagh",
    "pesant_regin",
    "alpha",
    "integer",
    "integer_alpha",
    "attainable",
    "attainable_alpha"
  ),
  scoring = NULL,
  n_items = 1,
  alpha = NULL,
  digits = 2,
  by = NULL,
  round_digits = NULL,
  rounding = "up"
) {
  rule <- rlang::arg_match(rule)
  # alias: named for its author
  if (rule == "muilwijk") {
    rule <- "mean"
  }
  check_rounding(
    rounding,
    setdiff(ROUNDING_METHODS, c("up_or_down", "ceiling_or_floor"))
  )
  if (
    rule %in% c("alpha", "integer_alpha", "attainable_alpha") && is.null(alpha)
  ) {
    cli::cli_abort("{.arg rule} = {.val {rule}} requires {.arg alpha}")
  }
  k <- n_items
  # back-compatible default: a composite is in mean-score units unless told
  # otherwise, and one item is a single item.
  if (is.null(scoring)) {
    scoring <- if (k > 1) "meanscored" else "singleitem"
  }
  scoring <- rlang::arg_match(
    scoring,
    c("singleitem", "sumscored", "meanscored")
  )
  geometry <- scoring_geometry(scoring, k, l, u)

  # lattice rules: only GRIM/GRIMMER-attainable reported tuples exist, so the
  # region is a set of points, not a band (no curve is defined at the means and
  # SDs that strictly integer data cannot produce).
  if (rule %in% c("integer", "integer_alpha")) {
    umbrella <- umbrella_data(
      l = l,
      u = u,
      n = n,
      digits = digits,
      granularity = "integer",
      scoring = scoring,
      n_items = k,
      alpha = if (rule == "integer_alpha") alpha else NULL
    )
    out <- umbrella[which(umbrella$consistent), c("mean", "sd")] |>
      tibble::as_tibble() |>
      round_lattice(round_digits, rounding)
    return(region_columns(out, "points"))
  }

  # exact lattice rules: the true attainable tuples, with no reporting grid.
  if (rule %in% c("attainable", "attainable_alpha")) {
    out <- attainable_lattice(l, u, n, geometry$mg)[c("mean", "sd")]
    if (rule == "attainable_alpha") {
      bounds <- unique(out$mean) |>
        purrr::map(function(mean_current) {
          bound <- sd_bounds(
            l = l,
            u = u,
            n = n,
            mean = mean_current,
            granularity = "integer",
            scoring = scoring,
            n_items = k,
            alpha = alpha
          )
          tibble::tibble(
            mean = mean_current,
            min_sd = bound$min_sd,
            max_sd = bound$max_sd,
            ok = isTRUE(bound$feasible) && !is.na(bound$min_sd)
          )
        }) |>
        dplyr::bind_rows()
      # the lattice is already sorted by mean and then sd, and the join keeps
      # its row order
      out <- out |>
        dplyr::inner_join(bounds, by = "mean")
      out <- out[
        out$ok &
          out$sd >= out$min_sd - TOLERANCE &
          out$sd <= out$max_sd + TOLERANCE,
        c("mean", "sd")
      ]
    }
    out <- round_lattice(out, round_digits, rounding)
    return(region_columns(out, "points"))
  }

  if (is.null(by)) {
    by <- (u - l) / MEAN_GRID_DIVISIONS
  }
  grid <- seq(l, u, by = by)
  mg <- geometry$mg
  quasi <- quasi_band(grid, n, l, u, mg)
  naive <- sd_min_integer(mg * grid, n) / mg # strict Bernoulli floor, off-grid too

  out <- switch(
    rule,
    range = tibble::tibble(
      mean = grid,
      min_sd = 0,
      max_sd = sd_max_span(l, u)
    ),
    range_n = tibble::tibble(
      mean = grid,
      min_sd = 0,
      max_sd = sd_max_span_n(l, u, n)
    ),
    mean = tibble::tibble(
      mean = grid,
      min_sd = 0,
      max_sd = sd_max_muilwijk(grid, n, l, u)
    ),
    mean_naive_floor = tibble::tibble(
      mean = grid,
      min_sd = naive,
      max_sd = sd_max_muilwijk(grid, n, l, u)
    ),
    mestdagh = tibble::tibble(mean = grid, min_sd = 0, max_sd = quasi$max_sd),
    pesant_regin = tibble::tibble(
      mean = grid,
      min_sd = quasi$min_sd,
      max_sd = sd_max_span_n(l, u, n)
    ),
    quasi = tibble::tibble(
      mean = grid,
      min_sd = quasi$min_sd,
      max_sd = quasi$max_sd
    ),
    alpha = {
      D <- 1 - (k - 1) / k * alpha
      if (k < 2) {
        cli::cli_abort(
          "{.code rule = \"alpha\"} needs {.code n_items >= 2} (alpha is inert at one item)"
        )
      }
      if (D <= TOLERANCE_TIGHT) {
        cli::cli_abort("{.arg alpha} too high for this {.arg n_items}")
      }
      v_max <- v_max_alpha(
        geometry$to_sum(grid),
        k,
        n,
        geometry$item_l,
        geometry$item_u
      )
      ceil_alpha <- sqrt((n / (n - 1)) * v_max / D) / geometry$sd_div
      tibble::tibble(
        mean = grid,
        min_sd = quasi$min_sd / sqrt(D), # alpha-amplified quasi-integer floor
        max_sd = pmin(ceil_alpha, quasi$max_sd) # alpha can only tighten
      )
    }
  )
  # an empty band (floor above ceiling) is infeasible, not a negative region
  bad <- out$min_sd > out$max_sd + TOLERANCE_TIGHT
  out$min_sd[bad] <- NA_real_
  out$max_sd[bad] <- NA_real_
  out <- region_columns(out, "band")
  attr(out, "step") <- by
  out
}

#' Plot the feasible SD region for a constraint set
#'
#' Draws what one constraint set from [sd_region_data()] asserts to be
#' possible: a band of sample SDs as a function of the mean, or a lattice of
#' attainable `(mean, sd)` tuples. A second band can be drawn as a dashed
#' reference, such as the sharp `"quasi"` band that the other rules are judged
#' against. Together these reproduce the panels of the nested-constraints
#' figure in the STRAIT article.
#'
#' @param data Output of [sd_region_data()]: a band, with `mean`, `min_sd`, and
#'   `max_sd`, or a lattice, with `mean` and `sd`, filtered by `consistent`
#'   when that column is present, as in [umbrella_data()] output.
#' @param reference Optional band to draw as a dashed reference, such as
#'   `sd_region_data(l, u, n, rule = "quasi")` for the same design.
#'   [sd_bounds_curve()] output works too.
#' @param title Optional plot title.
#' @param shade For a band, `"outside"` (default) shades the infeasible region
#'   and leaves the feasible one clear, matching [sd_bounds_plot()]; `"inside"`
#'   fills the band. The outside form is drawn with [band_polygon()], so a band
#'   that is undefined at some means — the alpha rules near each limit — leaves
#'   those means shaded rather than blank.
#' @param expand Padding around the plotted region, as a proportion of the
#'   range of the means (the scale width) rather than a fixed number of SD
#'   units, so the margin looks the same on a 1-5 scale and a 0-100 one.
#'   Applied to both axes.
#' @param fill,line_colour,reference_colour,point_colour,point_size Appearance.
#' @return A ggplot object.
#' @examples
#' quasi <- sd_region_data(l = 1, u = 5, n = 7, rule = "quasi")
#' sd_region_plot(quasi)                                   # sharp, GRIM-free
#' sd_region_plot(sd_region_data(l = 1, u = 5, n = 7, rule = "mean"),
#'                reference = quasi)                       # Bhatia-Davis
#' \donttest{
#' # digits = 1 keeps the reporting grid small: the lattice rules screen every
#' # cell with GRIMMER
#' sd_region_plot(sd_region_data(l = 1, u = 5, n = 7, rule = "integer",
#'                               digits = 1),
#'                reference = quasi)                       # attainable tuples
#' sd_region_plot(sd_region_data(l = 1, u = 5, n = 7, rule = "alpha",
#'                               n_items = 2, alpha = 0.7),
#'                reference = sd_region_data(l = 1, u = 5, n = 7,
#'                                           rule = "quasi", n_items = 2))
#' }
#' @export
sd_region_plot <- function(
  data,
  reference = NULL,
  title = NULL,
  shade = c("outside", "inside"),
  expand = 0.03,
  fill = "grey85",
  line_colour = "black",
  reference_colour = "grey45",
  point_colour = "#1d4ed8",
  point_size = 0.5
) {
  shade <- rlang::arg_match(shade)
  is_band <- is_band(data)
  if (is_band) {
    data <- band_mask(data)
  } else if ("consistent" %in% names(data)) {
    # full umbrella_data() output: only the consistent tuples exist
    data <- data[!is.na(data$consistent) & data$consistent, , drop = FALSE]
  }
  reference_masked <- if (!is.null(reference)) band_mask(reference)

  p <- ggplot2::ggplot()
  # the region itself first, so the dashed reference stays visible on top of it
  if (is_band) {
    if (shade == "inside") {
      p <- p +
        ggplot2::geom_ribbon(
          data = data,
          ggplot2::aes(
            x = .data$mean,
            ymin = .data$min_sd,
            ymax = .data$max_sd
          ),
          fill = fill,
          na.rm = TRUE
        )
    } else {
      # Knock the feasible rings out of a shaded panel rather than assembling
      # the shading around them. The alpha rules leave stretches near each
      # limit where no composite exists and `min_sd`/`max_sd` are NA; a ribbon
      # draws nothing there, so assembled shading would leave those means
      # unshaded and imply every SD is possible at them. See band_polygon().
      rings <- band_polygon(data, by = band_step(data))
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
      if (nrow(rings)) {
        p <- p +
          ggplot2::geom_polygon(
            data = rings,
            ggplot2::aes(x = .data$mean, y = .data$y, group = .data$ring),
            inherit.aes = FALSE,
            fill = "white"
          )
      }
    }
  }

  if (!is.null(reference_masked)) {
    p <- p +
      ggplot2::geom_line(
        data = reference_masked,
        ggplot2::aes(.data$mean, .data$max_sd),
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.25,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = reference_masked,
        ggplot2::aes(.data$mean, .data$min_sd),
        colour = reference_colour,
        linetype = "dashed",
        linewidth = 0.25,
        na.rm = TRUE
      )
  }

  if (is_band) {
    p <- p +
      ggplot2::geom_line(
        data = data,
        ggplot2::aes(.data$mean, .data$max_sd),
        colour = line_colour,
        linewidth = 0.4,
        na.rm = TRUE
      ) +
      ggplot2::geom_line(
        data = data,
        ggplot2::aes(.data$mean, .data$min_sd),
        colour = line_colour,
        linewidth = 0.4,
        na.rm = TRUE
      )
  } else {
    p <- p +
      ggplot2::geom_point(
        data = data,
        ggplot2::aes(.data$mean, .data$sd),
        colour = point_colour,
        size = point_size,
        na.rm = TRUE
      )
  }

  # limits are set from the data rather than left to ggplot2, because the
  # outside shading needs finite ones; padding is a proportion of the scale
  # width so the margin is constant across scales of very different widths.
  means <- range(c(data$mean, reference$mean), na.rm = TRUE)
  pad <- expand * diff(means)
  y_hi <- max(
    c(if (is_band) data$max_sd else data$sd, reference_masked$max_sd),
    na.rm = TRUE
  )

  p +
    ggplot2::coord_cartesian(
      xlim = means + c(-pad, pad),
      ylim = c(-pad, y_hi + pad),
      expand = FALSE
    ) +
    ggplot2::labs(x = "Mean", y = "Sample standard deviation", title = title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
