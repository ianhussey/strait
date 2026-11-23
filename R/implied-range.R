#' The range a reported (mean, SD) forces the data to cover
#'
#' How polarised must the data behind a report be? A sample whose values all
#' lie in `[x, y]` has an SD no larger than the ceiling for that narrower
#' range: the arch between `x` and `y` in the umbrella plot, which lies inside
#' the arch between the scale limits (see [sd_max_muilwijk()]). So a reported
#' SD can force the smallest observation down toward `l`, and the largest up
#' toward `u`. At the extreme, the data must contain the scale limits
#' themselves: Heathers' "horns", read off the summary statistics without
#' simulating a single sample.
#'
#' Each result is found by moving one scale limit inward until [brimmer()] no
#' longer calls the report consistent, so it can never disagree with the bounds
#' test. Moving a limit can only narrow the bounds, so the search is a
#' bisection. Under `granularity = "integer"` and `"quasiinteger"` the limit moves along
#' the response grid (the integers, or the `1 / n_items` lattice of a
#' mean-scored composite); under `"continuous"`, to within `1e-6` of the scale
#' range.
#'
#' Every result is a bound the data must meet, but not always the tightest
#' one: [brimmer()] is a necessary test, so it can accept a narrowed range that
#' no real sample fits. The range is exactly the one [brimmer()] can prove, and
#' real data may be forced to cover more.
#'
#' Each result is a constraint on its own: the smallest observation is at most
#' `min_at_most` (with the largest anywhere in range), and the largest is at
#' least `max_at_least`. `min_at_most == l` means the data must include `l`.
#' Where the report forces nothing, `min_at_most` is the top of the reported
#' mean's interval and `max_at_least` its bottom, since the smallest
#' observation never exceeds the mean, nor the largest falls short of it. An
#' attained extreme `a` or `b` already is the smallest or largest observation,
#' and is returned as given.
#'
#' This describes a report; it never flags one. A report that is inconsistent
#' to begin with has no sample to describe, and both results are `NA`. A
#' reported alpha is not supported: its bounds come from the item limits, which
#' narrowing the composite's range would misstate.
#'
#' @inheritParams brimmer
#' @param sd Numeric scalar, the reported SD (required).
#' @return A one-row tibble: `min_at_most`, `max_at_least`, and `note`. A
#'   missing input gives `NA` results, as in [brimmer()].
#' @seealso [brimmer()], whose question this extends from "is the report
#'   possible?" to "what must the data look like?".
#' @examples
#' # an SD near the ceiling of a 1-7 scale: the data must include both a 1
#' # and a 7
#' implied_range(l = 1, u = 7, n = 30, mean = 2.97, digits_mean = 2,
#'               sd = 2.83, digits_sd = 2, granularity = "integer")
#'
#' # a small SD forces only neighbouring values: at least one response of 3 or
#' # lower, and one of 5 or higher
#' implied_range(l = 1, u = 7, n = 30, mean = 4.00, digits_mean = 2,
#'               sd = 0.53, digits_sd = 2, granularity = "integer")
#' @export
implied_range <- function(
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
  n_items = 1
) {
  granularity <- match.arg(granularity)
  scoring <- match.arg(scoring)
  if (is.null(sd)) {
    stop("implied_range() needs a reported sd")
  }
  lower <- if (!is.null(a)) a else l
  upper <- if (!is.null(b)) b else u
  if (is.null(lower) || is.null(upper)) {
    stop("implied_range() needs both scale limits (l or a, and u or b)")
  }

  consistent <- function(l, u, a, b) {
    brimmer(
      l = l,
      u = u,
      a = a,
      b = b,
      n = n,
      mean = mean,
      digits_mean = digits_mean,
      sd = sd,
      digits_sd = digits_sd,
      rounding = rounding,
      granularity = granularity,
      scoring = scoring,
      n_items = n_items
    )$consistent
  }
  result <- function(min_at_most, max_at_least, note = NA_character_) {
    tibble::tibble(
      min_at_most = min_at_most,
      max_at_least = max_at_least,
      note = note
    )
  }
  first <- brimmer(
    l = l,
    u = u,
    a = a,
    b = b,
    n = n,
    mean = mean,
    digits_mean = digits_mean,
    sd = sd,
    digits_sd = digits_sd,
    rounding = rounding,
    granularity = granularity,
    scoring = scoring,
    n_items = n_items
  )
  if (is.na(first$consistent)) {
    return(result(NA_real_, NA_real_, first$note))
  }
  if (!first$consistent) {
    return(result(
      NA_real_,
      NA_real_,
      "the report is inconsistent, so no sample can produce it"
    ))
  }

  iv <- if (is.null(mean)) {
    list(lo = lower, hi = upper)
  } else if (is.null(rounding)) {
    list(lo = mean, hi = mean)
  } else {
    .unround(mean, digits_mean, rounding, "mean", "digits_mean")
  }
  granular <- granularity != "continuous"
  span <- upper - lower
  # the response grid a granular limit moves along
  unit <- if (scoring == "meanscored") 1 / n_items else 1

  # The first distance d in (0, d_max] by which a limit can move inward before
  # the report turns inconsistent, or NA if it stays consistent throughout;
  # ok(0) is known to hold. In grid steps for granular data.
  first_failure <- function(ok, d_max) {
    if (granular) {
      k_max <- floor(d_max / unit + 1e-9)
      if (k_max < 1 || ok(k_max * unit)) {
        return(NA)
      }
      lo <- 0
      hi <- k_max
      while (hi - lo > 1) {
        mid <- (lo + hi) %/% 2
        if (ok(mid * unit)) lo <- mid else hi <- mid
      }
      return(hi * unit)
    }
    if (ok(d_max)) {
      return(NA)
    }
    lo <- 0
    hi <- d_max
    while (hi - lo > 1e-6 * span) {
      mid <- (lo + hi) / 2
      if (ok(mid)) lo <- mid else hi <- mid
    }
    hi
  }
  # A limit stops short of the far end of the mean's interval, and of the
  # opposite limit. Once the report fails with the limit moved d inward, no
  # observation lies that far in: on the grid of integer data, the bound is
  # the grid value just outside d.
  d_room <- function(to_mean) {
    min(to_mean, span - if (granular) unit else 1e-9 * span)
  }
  back <- if (granularity == "integer") unit else 0
  # where nothing is forced, the mean's interval, snapped inward to that grid
  snap <- function(d) if (back > 0) floor(d / back + 1e-9) * back else d

  min_at_most <- if (!is.null(a)) {
    a
  } else {
    d <- first_failure(
      function(d) consistent(lower + d, u, NULL, b),
      d_room(iv$hi - lower)
    )
    if (is.na(d)) lower + snap(iv$hi - lower) else lower + d - back
  }
  max_at_least <- if (!is.null(b)) {
    b
  } else {
    d <- first_failure(
      function(d) consistent(l, upper - d, a, NULL),
      d_room(upper - iv$lo)
    )
    if (is.na(d)) upper - snap(upper - iv$lo) else upper - d + back
  }
  result(min_at_most, max_at_least)
}
