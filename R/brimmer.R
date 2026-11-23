# Report-level consistency checking (single row and batch) and the data
# builders (bound curves, umbrella grid) that drive the plots.

# ---- Layer 5: report-level consistency check ---------------------------------

#' BRIMMER: check a reported (mean, SD, n) against the SD bounds
#'
#' Bounds-Related Inconsistency of Means and Errors Reported. The SD-side
#' bounds test, and the package's main entry point: it asks whether a reported
#' standard deviation is arithmetically attainable given the scale limits, the
#' sample size and the reported mean. [brim()] is the mean-side test, which
#' needs no SD.
#'
#' The report-checking wrapper above [sd_bounds()]: give it the reported
#' summary statistics as they appear in a paper (numeric values plus their
#' reported decimal places), and it unrounds them, computes the bounds
#' envelope, runs the applicable consistency tests, and returns a single
#' consistent/inconsistent verdict with the reasons.
#'
#' BRIMMER is nested on top of BRIM, exactly as GRIMMER is nested on top of
#' GRIM: it applies the mean-side test and then adds the SD-side one. A report
#' is consistent only if every applicable test passes, and each names a
#' distinct defect.
#'
#' Tests applied (each only when its inputs are present):
#' * `in_scale_range` — the BRIM predicate: can the reported mean's rounding
#'   interval meet the feasible mean band at all? Bounds only, independent of
#'   granularity; with attained extremes the band narrows (see
#'   [feasible_mean_band()]). Run by [brim()] on its own.
#' * `bounds` — does the reported SD's own rounding interval overlap
#'   `[min_sd, max_sd]`? (The package's validated mathematics.)
#' * `feasibility` — a residual: the constraint set admits no sample for a
#'   reason none of the other tests already accounts for, such as an alpha
#'   floor exceeding the alpha ceiling. An out-of-range mean reports as
#'   `in_scale_range` and a granular-impossible one as `grim`, so neither
#'   double-counts here.
#' * `grim`, `grimmer` — scrutiny's verdicts, deferred verbatim (only under
#'   `granularity = "integer"` with `rounding`). NOTE scrutiny's GRIMMER (and, under
#'   `rounding = "up"`, GRIM) have documented floating-point boundary bugs,
#'   including false flags under the default rounding. When GRIMMER is the ONLY
#'   failing test, the note says so explicitly, since that pattern matches
#'   the documented false-flag family.
#'
#' `consistent` is TRUE iff every applicable test passes. `failed_tests`
#' names the failures (empty string if none). As in scrutiny, a missing
#' value among the inputs gives a missing verdict: `consistent` and
#' `failed_tests` are `NA`, and `note` names the missing input.
#'
#' Unlike [sd_bounds()], `rounding` here defaults to `"up_or_down"` — a
#' report-checking context virtually always deals with rounded values. Pass
#' `rounding = NULL` to treat the inputs as exact instead (digits then
#' unused).
#'
#' @param l,u Numeric scalars or NULL, logical limits.
#' @param a,b Numeric scalars or NULL, observed (attained) extremes.
#' @param n Integer scalar, sample size.
#' @param mean Numeric scalar or NULL, the reported mean.
#' @param digits_mean Integer scalar or NULL, its reported decimal places
#'   (required when `mean` and `rounding` are both given).
#' @param sd Numeric scalar or NULL, the reported SD — the thing being
#'   checked. NULL runs the mean-side tests only, which is what [brim()]
#'   does; the SD columns are then `NA`.
#' @param digits_sd Integer scalar or NULL, its reported decimal places
#'   (required when `rounding` is given).
#' @param rounding Rounding method for unrounding the reported values, named
#'   as in scrutiny (default `"up_or_down"`); NULL treats inputs as exact. See
#'   [sd_bounds()] and [roundwork::unround()].
#' @param granularity "continuous" (default), "integer", or "quasiinteger".
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1).
#' @param alpha Numeric scalar or NULL, reported Cronbach's alpha.
#'
#' @return A one-row tibble: `consistent` (logical), `failed_tests`
#'   (comma-separated names, "" if none), then the [sd_bounds()] columns
#'   `min_sd`, `max_sd`, `feasible`, the BRIM predicate `in_scale_range`, then
#'   `grim`, `grimmer`, `sd_in_bounds`, the
#'   percent-of-maximum-possible transforms `pomp_mean`, `pomp_sd_parity`
#'   (linear, against the mean-agnostic parity ceiling) and `pomp_sd_sharp`
#'   (position within the sharp mean-conditional band), and `note`.
#'
#' @seealso [brim()] for the mean-side test, [brimmest()] for the exact
#'   certificate that settles what this test can only screen, and
#'   [brimmer_map()] to apply this across a data frame.
#' @examples
#' # a perfectly ordinary report on a 1-7 scale
#' # (witnessed by c(rep(1, 20), rep(7, 9), 6): mean 2.9667, sd 2.8343)
#' brimmer(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'         sd = 2.83, digits_sd = 2, granularity = "integer")
#'
#' # SD below the quasi-integer floor: inconsistent via the bounds test
#' brimmer(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'         sd = 0.10, digits_sd = 2, granularity = "quasiinteger")
#'
#' # GRIM-impossible mean: inconsistent via feasibility (and grim)
#' brimmer(l = 1, u = 7, n = 30, mean = 3.51, digits_mean = 2,
#'         sd = 1.00, digits_sd = 2, granularity = "integer")
#' @export
brimmer <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  digits_mean = NULL,
  sd = NULL,
  digits_sd = NULL,
  rounding = "up_or_down",
  granularity = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  alpha = NULL
) {
  if (is.null(sd) && is.null(mean)) {
    stop("brimmer() needs a reported sd, a reported mean, or both")
  }
  if (is.null(n)) {
    stop("brimmer() requires n")
  }

  # scrutiny's NA rule: a missing input gives a missing verdict, rather than a
  # pass (an NA SD used to pass every SD-side test) or an opaque error
  supplied <- list(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = sd,
    digits_sd = digits_sd,
    n_items = n_items,
    alpha = alpha
  )
  na_in <- names(supplied)[vapply(
    supplied,
    function(v) {
      length(v) == 1L && is.na(v)
    },
    logical(1)
  )]
  if (length(na_in)) {
    return(.brimmer_na(paste("missing input:", paste(na_in, collapse = ", "))))
  }

  # a mean-free check is legitimate (e.g. sd vs the parity ceiling), but
  # rounding then applies to the sd only
  r <- sd_bounds(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = sd,
    digits_sd = digits_sd,
    rounding = if (is.null(mean)) NULL else rounding,
    granularity = granularity,
    scoring = scoring,
    n_items = n_items,
    alpha = alpha
  )

  # mean-free path: sd_bounds() skipped unrounding; redo the sd overlap with
  # the sd's own interval if rounding is in force
  if (is.null(mean) && !is.null(rounding)) {
    if (is.null(digits_sd)) {
      stop("rounding requires digits_sd for the reported sd")
    }
    iv <- .unround(sd, digits_sd, rounding, "sd", "digits_sd")
    r$sd_in_bounds <- if (is.na(r$min_sd)) {
      NA
    } else {
      (iv$hi >= r$min_sd - TOLERANCE && iv$lo <= r$max_sd + TOLERANCE)
    }
  }

  # POMP transforms on the reported scale (effective limits: a supersedes l,
  # b supersedes u). Descriptive, so they use the reported point mean/sd, not
  # the unrounding interval.
  eff_lower <- if (!is.null(a)) a else l
  eff_upper <- if (!is.null(b)) b else u

  # The BRIM predicate, computed independently of sd_bounds(). sd_bounds()
  # legitimately folds GRIM into `feasible` under granularity = "integer" (with no
  # GRIM-consistent mean there is no sample, so no bounds exist), but an
  # out-of-range mean and a granular-impossible mean are different defects and
  # must not be reported as one. This asks only: can the reported mean's
  # rounding interval meet the feasible mean band at all?
  in_scale_range <- NA
  if (
    !is.null(mean) &&
      (!is.null(l) || !is.null(u) || !is.null(a) || !is.null(b))
  ) {
    band <- feasible_mean_band(
      lower = eff_lower,
      upper = eff_upper,
      lower_attained = !is.null(a),
      upper_attained = !is.null(b),
      n = n
    )
    iv <- if (is.null(rounding)) {
      list(lo = mean, hi = mean)
    } else {
      .unround(mean, digits_mean, rounding, "mean", "digits_mean")
    }
    in_scale_range <- (iv$hi >= band[1] - TOLERANCE &&
      iv$lo <= band[2] + TOLERANCE)
  }

  failed <- character(0)
  if (isFALSE(in_scale_range)) {
    failed <- c(failed, "in_scale_range")
  }
  if (isFALSE(r$sd_in_bounds)) {
    failed <- c(failed, "bounds")
  }
  if (isFALSE(r$grim)) {
    failed <- c(failed, "grim")
  }
  if (isFALSE(r$grimmer)) {
    failed <- c(failed, "grimmer")
  }
  # `feasibility` is now the residual: the constraint set admits no sample for
  # a reason none of the named tests above already accounts for (e.g. an alpha
  # floor exceeding the alpha ceiling, or a GRIM divergence from scrutiny).
  if (!isTRUE(r$feasible) && !length(failed)) {
    failed <- c(failed, "feasibility")
  }

  note <- r$note
  if (identical(failed, "grimmer")) {
    caveat <- "only GRIMMER failed while feasibility and bounds pass; scrutiny's GRIMMER has documented false-flag cases - verify before flagging"
    note <- if (is.na(note)) caveat else paste(note, caveat, sep = "; ")
  }

  pomp <- .pomp_cols(
    mean,
    sd,
    r$min_sd,
    r$max_sd,
    eff_lower,
    eff_upper,
    n,
    has_mean = !is.null(mean)
  )

  # new_tibble(), as in sd_bounds(), since brimmer_map() calls this per row
  tibble::new_tibble(
    list(
      consistent = length(failed) == 0,
      failed_tests = paste(failed, collapse = ","),
      min_sd = r$min_sd,
      max_sd = r$max_sd,
      feasible = r$feasible,
      in_scale_range = in_scale_range,
      grim = r$grim,
      grimmer = r$grimmer,
      sd_in_bounds = r$sd_in_bounds,
      pomp_mean = pomp$pomp_mean,
      pomp_sd_parity = pomp$pomp_sd_parity,
      pomp_sd_sharp = pomp$pomp_sd_sharp,
      note = note
    ),
    nrow = 1L
  )
}

# Internal: brimmer()'s row for a report that cannot be judged.
.brimmer_na <- function(note) {
  tibble::new_tibble(
    list(
      consistent = NA,
      failed_tests = NA_character_,
      min_sd = NA_real_,
      max_sd = NA_real_,
      feasible = NA,
      in_scale_range = NA,
      grim = NA,
      grimmer = NA,
      sd_in_bounds = NA,
      pomp_mean = NA_real_,
      pomp_sd_parity = NA_real_,
      pomp_sd_sharp = NA_real_,
      note = note
    ),
    nrow = 1L
  )
}

#' BRIM: check a reported mean against the scale bounds
#'
#' Bounds-Related Inconsistency of Means. The mean-side test: is the reported
#' mean attainable at all, given the scale limits, the sample size and, when
#' supplied, the attained extremes? Where GRIM asks whether a mean is
#' attainable by strictly integer data at a given `n`, BRIM asks whether it is
#' attainable within the reporting range at all — a weaker condition that
#' applies to continuous data too, and one that tightens sharply once an
#' observed minimum or maximum is reported.
#'
#' A thin wrapper on [brimmer()] with no reported SD, so only the mean-side
#' tests apply: `in_scale_range` (the reported mean's rounding interval must
#' meet the feasible mean band, see [feasible_mean_band()]) and, under
#' `granularity = "integer"` with `rounding`, scrutiny's `grim` verdict. GRIMMER is never
#' invoked, and the SD-side columns are omitted rather than returned as `NA`.
#'
#' The two tests are independent and name different defects: a mean can be out
#' of range but granular-attainable (7.50 on a 1-7 scale at `n = 30`), inside
#' the range but GRIM-impossible (3.51 at `n = 30`), or both (7.51).
#'
#' @inheritParams brimmer
#' @return A one-row tibble: `consistent` (logical), `failed_tests`
#'   (comma-separated, "" if none), `in_scale_range`, `grim`, the feasible mean
#'   band `band_lo` and `band_hi`, `pomp_mean`, and `note`.
#' @seealso [brimmer()] for the SD-side test, [brimmest()] for the exact
#'   certificate.
#' @examples
#' # an ordinary mean on a 1-7 scale
#' brim(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2)
#'
#' # a mean above the scale maximum cannot be attained
#' brim(l = 1, u = 7, n = 30, mean = 7.5, digits_mean = 1)
#'
#' # attained extremes narrow the band, and can exclude a mean the bare
#' # scale limits would allow
#' brim(a = 1, b = 7, n = 30, mean = 1.10, digits_mean = 2)
#' @export
brim <- function(
  l = NULL,
  u = NULL,
  a = NULL,
  b = NULL,
  n = NULL,
  mean = NULL,
  digits_mean = NULL,
  rounding = "up_or_down",
  granularity = c("continuous", "integer", "quasiinteger"),
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1
) {
  if (is.null(mean)) {
    stop("brim() checks a reported mean: mean is required")
  }
  # sd = NULL, not NA: NULL is this package's "not supplied", whereas an NA
  # would propagate through the bounds arithmetic instead of switching the
  # SD-side tests off
  r <- brimmer(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = NULL,
    digits_sd = NULL,
    rounding = rounding,
    granularity = granularity,
    scoring = scoring,
    n_items = n_items
  )
  # same resolution of sides and attainment that sd_bounds() applies, so the
  # reported band is the one feasibility was actually tested against
  band <- feasible_mean_band(
    lower = if (!is.null(a)) a else l,
    upper = if (!is.null(b)) b else u,
    lower_attained = !is.null(a),
    upper_attained = !is.null(b),
    n = n
  )
  tibble::new_tibble(
    c(
      as.list(r[c("consistent", "failed_tests", "in_scale_range", "grim")]),
      list(band_lo = band[1], band_hi = band[2]),
      as.list(r[c("pomp_mean", "note")])
    ),
    nrow = 1L
  )
}

# ---- Layer 6: batch report checking ------------------------------------------

#' Check many reported (mean, SD, n) rows against the SD bounds
#'
#' Applies [brimmer()] to each row of a data frame, de-duplicating
#' identical constraint tuples so the (possibly expensive) rounding-envelope is
#' computed once per distinct input and reused. Columns of `data` whose names
#' match [brimmer()] arguments are taken per row; any argument given in
#' `...` is a constant broadcast to every row. Supplying one name both ways is
#' an error.
#'
#' Recognised names (column or constant): `l`, `u`, `a`, `b`, `n`, `mean`,
#' `digits_mean`, `sd`, `digits_sd`, `rounding`, `granularity`, `scoring`, `n_items`,
#' `alpha`.
#'
#' @param data A data frame, one reported statistic set per row.
#' @param ... Constant arguments applied to all rows.
#' @param include_inputs If TRUE (default), returns `data` column-bound to the
#'   results; if FALSE, only the result columns (same row order), for drop-in
#'   use inside a `dplyr::mutate()`/`purrr` pipeline.
#' @return A tibble of the [brimmer()] columns, one row per input row
#'   (optionally with the inputs prepended). A row with a missing input gets
#'   a missing verdict, as in [brimmer()].
#' @examples
#' reports <- data.frame(mean = c(2.97, 3.51, 4.20),
#'                       sd   = c(2.83, 3.50, 0.90),
#'                       n    = c(30, 30, 30))
#'
#' # the scale limits and reported precision are constant across rows;
#' # row 2 reports an SD above the ceiling for a 1-7 scale at n = 30
#' out <- brimmer_map(reports, l = 1, u = 7,
#'                                 digits_mean = 2, digits_sd = 2)
#' out[, c("mean", "sd", "n", "consistent", "failed_tests")]
#' @export
brimmer_map <- function(data, ..., include_inputs = TRUE) {
  if (!is.data.frame(data)) {
    stop("data must be a data frame")
  }
  arg_names <- c(
    "l",
    "u",
    "a",
    "b",
    "n",
    "mean",
    "digits_mean",
    "sd",
    "digits_sd",
    "rounding",
    "granularity",
    "scoring",
    "n_items",
    "alpha"
  )
  consts <- list(...)
  unknown <- setdiff(names(consts), arg_names)
  if (length(unknown)) {
    stop("unknown constant argument(s): ", paste(unknown, collapse = ", "))
  }
  N <- nrow(data)
  resolve <- function(nm) {
    incol <- nm %in% names(data)
    incon <- nm %in% names(consts) && !is.null(consts[[nm]])
    if (incol && incon) {
      stop(sprintf("'%s' supplied as both a column and a constant", nm))
    }
    v <- if (incol) {
      data[[nm]]
    } else if (incon) {
      rep(consts[[nm]], length.out = N)
    }
    if (is.factor(v)) as.character(v) else v
  }
  cols <- lapply(arg_names, resolve)
  names(cols) <- arg_names
  present <- arg_names[!vapply(cols, is.null, logical(1))]
  # a constant given as NULL on purpose, such as rounding = NULL for exact
  # inputs, is passed on rather than replaced by brimmer()'s default
  null_consts <- names(consts)[vapply(consts, is.null, logical(1))]
  if (!("sd" %in% present)) {
    stop("a reported sd is required (as a column of data or a constant)")
  }

  # de-duplicate identical input tuples; compute once per unique tuple
  key <- do.call(
    paste,
    c(
      lapply(present, function(nm) {
        format(cols[[nm]], nsmall = 6, trim = TRUE)
      }),
      sep = "\r"
    )
  )
  uk_idx <- which(!duplicated(key))
  back <- match(key, key[uk_idx])
  res_uni <- do.call(
    rbind,
    lapply(uk_idx, function(i) {
      args <- lapply(present, function(nm) cols[[nm]][i])
      names(args) <- present
      args[null_consts] <- list(NULL)
      do.call(brimmer, args)
    })
  )
  res <- res_uni[back, , drop = FALSE]
  if (include_inputs) .bind_inputs(data, res) else res
}

# Internal: the inputs, then the results, as one tibble. Names are kept as
# they are, so an input column that shares a result's name stays visible.
.bind_inputs <- function(data, res) {
  tibble::as_tibble(
    cbind(as.data.frame(data), as.data.frame(res)),
    .name_repair = "minimal"
  )
}

# ---- Layer 7: bound curves, umbrella grid, and plots -------------------------

#' SD-bound curves across the mean (hole-free under granularity)
#'
#' Traces the minimum and maximum sample SD as the mean sweeps `[l, u]`, on a
#' grid dense enough to resolve the piecewise floors and the Structure S ceiling
#' (their kinks fall on the `1/(n * m)` lattice, `m = n_items` for mean-scored
#' data else 1). Under a granularity constraint the quasi-integer floor is
#' defined at every mean, so the curves have no gaps — unlike a strict-integer
#' floor, which exists only at GRIM means. Returns raw bounds and their
#' parity-normalised counterparts for the POMP plot.
#'
#' @param l,u Numeric scalars, limits.
#' @param n Integer scalar, sample size.
#' @param granularity `"quasiinteger"` (default) gives a hole-free floor.
#' @param scoring,n_items,alpha As in [sd_bounds()].
#' @param by Numeric or NULL; mean-grid spacing (default ~ `(u - l) / 1000`).
#' @return A tibble: `mean`, `min_sd`, `max_sd`, `feasible`, `pomp_mean`,
#'   `parity_max`, `ceil_parity` (= `max_sd / parity_max`), `floor_parity`.
#'   The `"step"` attribute records the uniform grid spacing `by` actually
#'   used. The mean grid is deliberately NOT uniform - the kinks and their
#'   neighbourhoods are sampled far more densely than `by` - so the spacing
#'   cannot be recovered from the returned means, and a consumer that needs to
#'   tell a sampling gap from a genuine gap in the band (as [band_polygon()]
#'   and hence [sd_bounds_plot()] must) has to be told. See [sd_bounds_plot()].
#' @examples
#' curve <- sd_bounds_curve(l = 1, u = 7, n = 30, by = 0.1)
#' head(curve[, c("mean", "min_sd", "max_sd")])
#'
#' # the ceiling peaks near the scale midpoint
#' curve[which.max(curve$max_sd), c("mean", "min_sd", "max_sd")]
#' @export
sd_bounds_curve <- function(
  l,
  u,
  n,
  granularity = "quasiinteger",
  scoring = "singleitem",
  n_items = 1,
  alpha = NULL,
  by = NULL
) {
  m <- if (scoring == "meanscored") n_items else 1L
  if (is.null(by)) {
    by <- (u - l) / MEAN_GRID_DIVISIONS
  }
  nm <- n * m
  kinks <- unique(c(
    seq(ceiling(l * nm), floor(u * nm)) / nm,
    seq(ceiling(l * m), floor(u * m)) / m
  ))
  means <- sort(unique(pmin(
    u,
    pmax(l, c(seq(l, u, by = by), kinks, kinks + 1e-9, kinks - 1e-9))
  )))
  # one row per mean, collected column-wise: binding thousands of one-row
  # tables took longer than computing the bounds
  rows <- means |>
    lapply(\(mu) {
      sd_bounds(
        l = l,
        u = u,
        n = n,
        mean = mu,
        granularity = granularity,
        scoring = scoring,
        n_items = n_items,
        alpha = alpha
      )
    })
  pick <- \(name, type) vapply(rows, \(d) d[[name]], type)
  min_sd <- pick("min_sd", numeric(1))
  max_sd <- pick("max_sd", numeric(1))
  parity_max <- sd_max_span_n(l, u, n)
  out <- tibble::tibble(
    mean = means,
    min_sd = min_sd,
    max_sd = max_sd,
    feasible = pick("feasible", logical(1)),
    pomp_mean = (means - l) / (u - l),
    parity_max = parity_max,
    ceil_parity = max_sd / parity_max,
    floor_parity = min_sd / parity_max
  )
  attr(out, "step") <- by
  out
}

#' Build the GRIM x GRIMMER x bounds umbrella grid (full, with verdicts)
#'
#' For each reported mean on the `10^-digits` grid over `[l, u]`, computes the
#' SD-bounds envelope once, then evaluates every reported SD on the same grid
#' from 0 up to that mean's ceiling, tagging each `(mean, sd)` cell with the
#' bounds-overlap and (under `granularity = "integer"`) GRIMMER verdicts. This returns the
#' FULL grid — feasible means x candidate SDs — with verdict columns, rather
#' than only the passing tuples, so the plotting layer can render failures as
#' well as the passing "umbrella".
#'
#' @param l,u Numeric scalars, limits.
#' @param n Integer scalar, sample size.
#' @param digits Integer, reported decimal places (grid step `10^-digits`).
#' @param granularity Response granularity (default `"integer"`).
#' @param scoring,n_items As in [sd_bounds()].
#' @param alpha Optional reported Cronbach's alpha; when supplied the bounds each
#'   cell is tested against are the alpha-conditional ones, so the grid shows the
#'   jointly GRIM-, GRIMMER- and alpha-consistent tuples.
#' @param rounding Rounding method for mean and SD, named as in scrutiny
#'   (default `"up_or_down"`).
#' @return A tibble: `mean`, `sd`, `min_sd`, `max_sd`, `in_bounds`,
#'   `grimmer`, `consistent`. `grimmer` is evaluated only for SDs inside the
#'   sharp bounds (elsewhere the tuple is already inconsistent and `grimmer`
#'   is `NA`); GRIM-inconsistent means are pruned before any SD is tested,
#'   which is why they are absent from the grid.
#' @examples
#' # the full grid of reportable (mean, sd) pairs for a small design
#' # (a 3-point scale at 1 decimal place, to keep the example quick)
#' grid <- umbrella_data(l = 1, u = 3, n = 12, digits = 1)
#' head(grid)
#'
#' # how many pairs survive every test
#' table(grid$consistent)
#' @export
umbrella_data <- function(
  l,
  u,
  n,
  digits = 2,
  granularity = "integer",
  scoring = "singleitem",
  n_items = 1,
  alpha = NULL,
  rounding = "up_or_down"
) {
  .check_rounding(rounding)
  step <- 10^(-digits)
  # round() the grids back onto the nearest double to each decimal. seq() by a
  # decimal step accumulates error (seq(0, 6, by = 0.1)[4] is 0.3 + 5.6e-17),
  # and the granularity tests read the value as a decimal string, so an
  # off-by-one-ulp grid point can flip a GRIMMER verdict.
  means <- round(seq(l, u, by = step), digits)
  # items averaged per observation, as in sd_bounds()
  m <- if (scoring == "meanscored") n_items else 1L
  use_grimmer <- granularity == "integer"
  rows <- list()
  for (mu in means) {
    # sd_bounds() embeds the rounding-aware GRIM prefilter: a mean whose
    # rounding interval admits no integer sum is infeasible under strict granularity,
    # so none of its SDs need testing.
    d <- sd_bounds(
      l = l,
      u = u,
      n = n,
      mean = mu,
      digits_mean = digits,
      rounding = rounding,
      granularity = granularity,
      scoring = scoring,
      n_items = n_items,
      alpha = alpha
    )
    if (!isTRUE(d$feasible) || is.na(d$max_sd)) {
      next
    }
    sds <- round(seq(0, ceiling(d$max_sd / step) * step, by = step), digits)
    # each SD's own rounding interval, as sd_bounds() unrounds a reported SD
    iv <- .unround(sds, digits, rounding, "sd", "digits")
    in_bounds <- iv$hi >= d$min_sd - TOLERANCE & iv$lo <= d$max_sd + TOLERANCE
    # GRIMMER (the expensive per-tuple test) runs only where the SD is inside
    # the sharp bounds; outside, the tuple is already inconsistent and the
    # verdict is left NA.
    grimmer <- rep(NA, length(sds))
    if (use_grimmer && any(in_bounds)) {
      grimmer[in_bounds] <- as.logical(scrutiny::grimmer(
        x = mu,
        sd = sds[in_bounds],
        n = n,
        digits_x = digits,
        digits_sd = digits,
        items_averaged = m,
        data_type = "integer",
        rounding = rounding
      ))
    }
    consistent <- in_bounds &
      (if (use_grimmer) !is.na(grimmer) & grimmer else TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      mean = mu,
      sd = sds,
      min_sd = d$min_sd,
      max_sd = d$max_sd,
      in_bounds = in_bounds,
      grimmer = grimmer,
      consistent = consistent
    )
  }
  tibble::as_tibble(do.call(rbind, rows))
}

# Internal: shared point layer (green/red outlined dots)
