# bounds-primitives.R
#
# The bounds of the sample standard deviation under nested constraints, as a
# family of small single-purpose functions governed by the sd_bounds()
# dispatcher (see sd_bounds.R). The mathematics is the closed-form (not
# constructive) derivation developed in the STRAIT article; each primitive here
# implements one formula.
#
# ARCHITECTURE
#   Layer 0  utilities: bessel_factor, frac, is_grim_consistent, sd_from_ss
#   Layer 1  bound primitives (one formula each):
#            ceilings:  sd_max_span, sd_max_span_n, sd_max_structure_s
#            floors:    sd_min_integer, sd_min_quasi_integer   (range-free)
#                       sd_min_two_pin, sd_min_one_pin         (attained)
#            feasibility: feasible_mean_band
#   Layer 2  alpha composite layer: v_max_alpha, sd_bounds_alpha
#   Layer 3  rounding layer: .check_rounding, .unround (via roundwork)
#
# CONVENTIONS
#   - All SDs are sample SDs (Bessel, n - 1 divisor).
#   - l, u are LOGICAL limits (walls: membership only); a, b are OBSERVED
#     extremes (attained: at least one observation equals each). a supersedes l
#     and b supersedes u, because a wall beyond an attained extreme is vacuous.
#   - Z is the granularity constraint: "continuous" (none), "integer" (all
#     observations on the grid), "quasiinteger" (all but one; GRIM-free).
#   - Sharp mean-conditional ceilings are Structure S in every case: it equals
#     the delta-corrected Muilwijk bound (verified to 1e-13), is valid for
#     continuous data (endpoint counts must be whole regardless), and is
#     unchanged by attainment inside the feasible mean band.
#   - A reported mean may be treated as exact (default) or as a rounded or
#     truncated value via `digits_mean` + `rounding` (scrutiny's options); the
#     bounds are then the ENVELOPE over all exact means consistent with the
#     report: min over the interval of the floor, max over it of the ceiling.

# ---- Layer 0: utilities ------------------------------------------------------

#' Bessel factor
#'
#' The factor `sqrt(n / (n - 1))` converting a population SD of an attaining
#' configuration into the sample SD.
#'
#' @param n Integer scalar or vector, sample size(s), `n >= 2`.
#' @return Numeric, `sqrt(n / (n - 1))`.
#' @keywords internal
bessel_factor <- function(n) sqrt(n / (n - 1))

#' Fractional part, robust to floating-point dust near an integer
#'
#' `x - floor(x)`, but a value within `tol` of an integer (on either side) is
#' snapped to `0`. This matters because scaling a rounded reported mean onto the
#' integer grid (`w = m * x`, mean-scored data) can land a hair away from an
#' integer — e.g. `round(4/3, 12) * 3 = 3.999999999999` (fractional part `~1`)
#' or `round(5/3, 12) * 3 = 5.000000000100` (fractional part `~1e-10`) — where a
#' naive fractional part produces a spurious granularity floor
#' `sqrt(d(1 - d))`. The snap is safe: the smallest legitimate on-grid
#' fractional part is `1/n`, orders of magnitude larger than `tol`, so a
#' genuinely off-grid mean such as `2.999` (fractional part `0.999`) is left
#' untouched.
#'
#' @param x Numeric vector.
#' @param tol Numeric tolerance for snapping to the nearest integer.
#' @return `x - floor(x)`, in `[0, 1)`, dust near either integer snapped to 0.
#' @keywords internal
frac <- function(x, tol = 1e-9) {
  r <- x - floor(x)
  ifelse(r < tol | r > 1 - tol, 0, r)
}

#' Is a mean GRIM-consistent?
#'
#' A mean is GRIM-consistent at sample size `n` when `n * mean` is an integer,
#' i.e. some sample of `n` integers has exactly that mean.
#'
#' @param mean Numeric scalar, the (exact) mean.
#' @param n Integer scalar, sample size.
#' @param tol Numeric tolerance on `n * mean`'s distance from an integer.
#' @return Logical scalar.
#' @keywords internal
is_grim_consistent <- function(mean, n, tol = 1e-9) {
  abs(n * mean - round(n * mean)) < tol
}

#' Sample SD from a sum of squared deviations
#'
#' @param ss Numeric vector, sum of squared deviations about the mean.
#' @param n Integer scalar, sample size.
#' @return `sqrt(pmax(0, ss) / (n - 1))`.
#' @keywords internal
sd_from_ss <- function(ss, n) sqrt(pmax(0, ss) / (n - 1))

# ---- Layer 1: ceilings -------------------------------------------------------

#' Maximum SD given two limits only (no n)
#'
#' The n-free maximum over all sample sizes: attained at `n = 2` with one
#' observation at each limit (Popoviciu's configuration).
#' Valid whether the limits are walls or attained extremes.
#'
#' @param lower,upper Numeric scalars, the effective lower and upper limits.
#' @return `(upper - lower) / sqrt(2)`.
#' @examples
#' # the n-free maximum for a 1-7 scale, attained at n = 2
#' sd_max_span(lower = 1, upper = 7)
#' @export
sd_max_span <- function(lower, upper) (upper - lower) / sqrt(2)

#' Maximum SD given two limits and n (parity ceiling)
#'
#' The mean-agnostic maximum for a given sample size: half the observations at
#' each limit, with the odd-n parity correction (Popoviciu 1935 sec 4,
#' Petocz 2005). Unchanged by attainment (the attaining
#' configuration occupies both limits).
#'
#' @param lower,upper Numeric scalars, the effective limits.
#' @param n Integer vector, sample size(s), `n >= 2`.
#' @return Numeric, the maximum sample SD.
#' @examples
#' # an even n splits exactly; an odd n pays the parity correction
#' sd_max_span_n(lower = 1, upper = 7, n = 30)
#' sd_max_span_n(lower = 1, upper = 7, n = 31)
#'
#' # both sit below the n-free maximum
#' sd_max_span(lower = 1, upper = 7)
#' @export
sd_max_span_n <- function(lower, upper, n) {
  (upper - lower) /
    2 *
    ifelse(n %% 2 == 0, sqrt(n / (n - 1)), sqrt((n + 1) / n))
}

#' Smooth mean-conditional maximum SD (Muilwijk / Bhatia-Davis)
#'
#' The mean-conditional ceiling `sqrt(n/(n-1)) * sqrt((upper - mean)(mean - lower))`,
#' due to Muilwijk (1966) and rediscovered by Bhatia and Davis (2000). It is the
#' smooth arch obtained by allowing fractional counts at the two limits, so it is
#' *not* sharp: the sharp ceiling applies the count-parity correction and is
#' [sd_max_structure_s()], which never exceeds this. Retained because it is the
#' form most often cited in the literature and is useful as a reference curve.
#'
#' For any two values `lower < upper`, not just the scale limits, it is the SD
#' of a sample on those two values alone, as a function of its mean: an arch
#' from `lower` to `upper` (a semicircle for the population SD, with mean and SD
#' in the same units). The umbrella plot is made of these arches. Its ceiling is
#' the arch between the scale limits; its integer floor, [sd_min_integer()], is
#' the arch between the two whole numbers either side of the mean; and the arch
#' between any other two values is the ceiling for data confined to them, which
#' [implied_range()] uses to find the range a report forces the data to cover.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Numeric scalars, the limits.
#' @return Numeric vector, the Muilwijk / Bhatia-Davis maximum sample SD.
#' @examples
#' # the smooth arch, highest at the scale midpoint
#' sd_max_muilwijk(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#'
#' # it is not sharp: the parity-corrected ceiling never exceeds it
#' sd_max_structure_s(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#' @export
sd_max_muilwijk <- function(mean, n, lower, upper) {
  bessel_factor(n) * sqrt(pmax(0, (upper - mean) * (mean - lower)))
}

#' Sharp mean-conditional maximum SD (Structure S)
#'
#' The sharp maximum given limits, n, and the mean: `n_u` observations at the
#' upper limit, `n_l` at the lower, at most one interior remainder (Mestdagh et
#' al. 2018 "Structure S"). Equal to the count-parity (delta) corrected Muilwijk
#' bound, so it is the sharp ceiling for continuous data too (endpoint counts
#' must be whole numbers regardless of data continuity), and it is unchanged by
#' attainment of the limits inside the feasible mean band. Under the
#' quasi-integer constraint the remainder observation is the one free real.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Numeric scalars, the effective limits.
#' @return Numeric vector, the maximum sample SD at each mean.
#' @examples
#' # the sharp mean-conditional ceiling across a 1-7 scale
#' sd_max_structure_s(mean = c(2, 4, 6), n = 30, lower = 1, upper = 7)
#'
#' # a mean at a limit forces every observation there, so the SD is 0
#' sd_max_structure_s(mean = 1, n = 30, lower = 1, upper = 7)
#' @export
sd_max_structure_s <- function(mean, n, lower, upper) {
  R <- upper - lower
  S_excess <- n * (mean - lower)
  n_u <- pmin(floor(S_excess / R + 1e-9), n - 1)
  n_l <- n - n_u - 1
  # deviations from the mean, each found from e = mean - lower: the sum of
  # squares is tiny next to a wall, and as a difference of squared raw values,
  # n_l * lower^2 + x_r^2 + n_u * upper^2 - n * mean^2, it cancels to 0
  e <- mean - lower
  ss <- n_l * e^2 + (S_excess - n_u * R - e)^2 + n_u * (R - e)^2
  sd_from_ss(ss, n)
}

# ---- Layer 1: range-free floors (granularity only) ---------------------------

#' Minimum SD for strictly integer data (Bernoulli floor)
#'
#' All observations on the two integers adjacent to the mean (Pesant & Regin
#' 2005; Bernoulli form). Range-free. Defined only at
#' GRIM-consistent means; the caller is responsible for that check (see
#' [is_grim_consistent()]), otherwise use [sd_min_quasi_integer()].
#'
#' It is the same arch as the ceiling, [sd_max_muilwijk()], taken between the
#' two integers either side of the mean instead of between the scale limits:
#' the floor and the ceiling of the umbrella plot are both SDs of two-value
#' samples, at the narrowest and the widest spacing.
#'
#' @param mean Numeric vector, GRIM-consistent mean(s).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @return Numeric vector, the minimum sample SD at each mean.
#' @examples
#' # zero at a whole-number mean, largest half way between two integers
#' sd_min_integer(mean = 3, n = 30)
#' sd_min_integer(mean = 3.5, n = 30)
#' @export
sd_min_integer <- function(mean, n) {
  d <- frac(mean)
  bessel_factor(n) * sqrt(d * (1 - d))
}

#' Minimum SD for quasi-integer data (GRIM-free floor)
#'
#' All but one observation an integer (the quasi-integer constraint, i.e. all
#' but one observation on the grid). Defined at every real mean; coincides with
#' [sd_min_integer()] at GRIM-consistent means and dips between them, so it is
#' a valid floor under either the strict or the relaxed hypothesis. Range-free.
#'
#' @param mean Numeric vector, mean(s).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @return Numeric vector, the minimum sample SD at each mean.
#' @examples
#' # at a GRIM-consistent mean it agrees with the strict-integer floor
#' sd_min_quasi_integer(mean = 3.5, n = 30)
#' sd_min_integer(mean = 3.5, n = 30)
#'
#' # between GRIM means it dips below, so it is defined at every mean
#' sd_min_quasi_integer(mean = 3.51, n = 30)
#' @export
sd_min_quasi_integer <- function(mean, n) {
  d <- frac(mean)
  g <- frac(n * mean)
  sd_from_ss(n * d * (1 - d) - g * (1 - g), n)
}

# ---- Layer 1: attained-extremes floors ---------------------------------------

#' Minimum SD with both observed extremes attained (two-pin floor)
#'
#' One observation pinned at each observed extreme. Without a mean this is the
#' Nagy (1918) / Thomson (1955) floor `W / sqrt(2 (n - 1))`. With a mean it is
#' the pinned-plus-interior decomposition of the attained-extremes analysis:
#' `SS = p^2 + q^2 + (p - q)^2 / (n - 2)` plus, for integer data, the
#' Pesant-Regin clustering of the `n - 2` interior observations, and for
#' quasi-integer data the `g (1 - g)` relaxation. The mean must lie in the
#' feasible band `[a + W/n, b - W/n]` (see [feasible_mean_band()]); this
#' function does not check it.
#'
#' @param a,b Numeric scalars, observed minimum and maximum (attained).
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param mean Numeric scalar or NULL. If NULL, the mean-agnostic floor.
#' @param Z One of "continuous", "integer", "quasiinteger".
#' @return Numeric scalar, the minimum sample SD.
#' @examples
#' # mean-agnostic floor: one observation pinned at each observed extreme
#' sd_min_two_pin(a = 1, b = 7, n = 30)
#'
#' # knowing the mean sharpens it, and granularity sharpens it again
#' sd_min_two_pin(a = 1, b = 7, n = 30, mean = 3)
#' sd_min_two_pin(a = 1, b = 7, n = 30, mean = 3, Z = "integer")
#' @export
sd_min_two_pin <- function(a, b, n, mean = NULL, Z = "continuous") {
  W <- b - a
  if (is.null(mean)) {
    return(W / sqrt(2 * (n - 1)))
  }
  if (n == 2) {
    return(W / sqrt(2))
  } # band is the single midpoint mean
  k <- n - 2
  p <- mean - a
  q <- b - mean
  ss_cont <- p^2 + q^2 + (p - q)^2 / k
  ss <- ss_cont
  if (Z %in% c("integer", "quasiinteger")) {
    T_int <- n * mean - a - b # interior sum
    d <- frac(T_int / k) # interior mean's fractional part
    ss <- ss_cont + k * d * (1 - d)
    if (Z == "quasiinteger") {
      ss <- ss - frac(T_int) * (1 - frac(T_int))
    }
    ss <- max(ss, ss_cont)
  }
  sd_from_ss(ss, n)
}

#' Minimum SD with one observed extreme attained (one-pin floor)
#'
#' One observation pinned at a single attained extreme, the other side at most
#' walled (the attained-extremes analysis, mixed constraint sets). NOT the
#' two-pin formula with a pin deleted: the interior count is `n - 1`, not
#' `n - 2`. Continuous floor `SS = n q^2 / (n - 1)` (strictly sharper than the
#' Laguerre-Samuelson corollary), plus the same granularity terms on the
#' `n - 1` free observations. Requires a mean (without one the floor is zero:
#' all observations may sit at the pin). The mirror case is handled by
#' reflection.
#'
#' @param pin Numeric scalar, the attained extreme.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param mean Numeric scalar.
#' @param Z One of "continuous", "integer", "quasiinteger".
#' @param side "max" if the pin is the observed maximum, "min" if the minimum.
#' @return Numeric scalar, the minimum sample SD.
#' @examples
#' # an observed maximum of 7 alongside a mean of 3 forces some spread
#' sd_min_one_pin(pin = 7, n = 30, mean = 3)
#'
#' # the mirror case: an observed minimum of 1 alongside a mean of 5
#' sd_min_one_pin(pin = 1, n = 30, mean = 5, side = "min")
#' @export
sd_min_one_pin <- function(
  pin,
  n,
  mean,
  Z = "continuous",
  side = c("max", "min")
) {
  side <- match.arg(side)
  # reflect the observed-minimum case onto the observed-maximum formulas
  if (side == "min") {
    return(sd_min_one_pin(-pin, n, -mean, Z, side = "max"))
  }
  k <- n - 1
  q <- pin - mean # distance from mean up to the pin
  ss_cont <- n * q^2 / (n - 1)
  ss <- ss_cont
  if (Z %in% c("integer", "quasiinteger")) {
    T_int <- n * mean - pin # sum of the n - 1 free observations
    d <- frac(T_int / k)
    ss <- ss_cont + k * d * (1 - d)
    if (Z == "quasiinteger") {
      ss <- ss - frac(T_int) * (1 - frac(T_int))
    }
    ss <- max(ss, ss_cont)
  }
  sd_from_ss(ss, n)
}

# ---- Layer 1: feasibility ----------------------------------------------------

#' Feasible mean interval given the side constraints
#'
#' Walls admit any mean in `[lower, upper]`. Each attained extreme pins one
#' observation and narrows the band: both attained gives
#' `[a + W/n, b - W/n]`; an attained maximum with (at most) a walled minimum
#' gives `[lower + (b - lower)/n, b]` (and the mirror for an attained
#' minimum). `-Inf`/`Inf` when a side is absent.
#'
#' @param lower,upper Numeric scalars or NULL, the effective limits.
#' @param lower_attained,upper_attained Logical, is each limit attained?
#' @param n Integer scalar or NULL.
#' @return Numeric length-2 vector, the closed feasible mean interval.
#' @examples
#' # walls alone admit any mean on the scale
#' feasible_mean_band(lower = 1, upper = 7, n = 30)
#'
#' # an attained maximum pins one observation at 7, lifting the lowest mean
#' feasible_mean_band(lower = 1, upper = 7, upper_attained = TRUE, n = 30)
#'
#' # both extremes attained narrows the band at each end
#' feasible_mean_band(lower = 1, upper = 7, lower_attained = TRUE,
#'                    upper_attained = TRUE, n = 30)
#' @export
feasible_mean_band <- function(
  lower = NULL,
  upper = NULL,
  lower_attained = FALSE,
  upper_attained = FALSE,
  n = NULL
) {
  lo <- if (is.null(lower)) -Inf else lower
  hi <- if (is.null(upper)) Inf else upper
  if (!is.null(n)) {
    if (lower_attained && upper_attained) {
      W <- upper - lower
      return(c(lower + W / n, upper - W / n))
    }
    if (upper_attained && is.finite(lo)) {
      lo <- max(lo, lower + (upper - lower) / n)
    }
    if (lower_attained && is.finite(hi)) {
      hi <- min(hi, upper - (upper - lower) / n)
    }
  }
  c(lo, hi)
}

# ---- Layer 2: alpha composite layer ------------------------------------------

#' Allocation maximum of the item-variance sum (Theorem H3)
#'
#' The quasi-integer maximum of the summed item variances of a k-item integer
#' battery at a given sum-score mean (the STRAIT article, Appendix H).
#'
#' @param mean_sum Numeric vector, sum-score mean(s).
#' @param k Integer scalar, number of items.
#' @param n Integer scalar, sample size.
#' @param item_l,item_u Numeric scalars, the per-item limits.
#' @return Numeric vector, the maximum of the population item-variance sum.
#' @examples
#' # a 3-item 1-5 battery (sum score 3-15) at a sum-score mean of 9
#' v_max_alpha(mean_sum = 9, k = 3, n = 50, item_l = 1, item_u = 5)
#' @keywords internal
#' @export
v_max_alpha <- function(mean_sum, k, n, item_l, item_u) {
  h <- (item_u - item_l) / n
  T_ <- (mean_sum - k * item_l) / h
  x <- T_ / k
  j <- pmin(floor(x + 1e-12), n - 1)
  th <- x - j
  phi <- T_ - floor(T_ + 1e-12)
  phi <- ifelse(phi > 1 - 1e-9, 0, phi)
  A <- function(jj) h^2 * jj * (n - jj)
  pmax(
    0,
    k * ((1 - th) * A(j) + th * A(j + 1)) - h^2 * (n - 1) * phi * (1 - phi)
  )
}

# Internal: all non-negative integer count vectors (c_0, ..., c_{slots-1}) that
# sum to `total`. Rows enumerate compositions; used for the exact Gini envelope.
.count_vectors <- function(slots, total) {
  if (slots == 1L) {
    return(matrix(total, ncol = 1))
  }
  do.call(
    rbind,
    lapply(0:total, function(k) {
      cbind(k, .count_vectors(slots - 1L, total - k))
    })
  )
}

#' Sharp alpha-conditional composite floor via the Gini mean difference
#'
#' The exact alpha-conditional minimum SD for a *strictly integer* composite
#' (STRAIT article, Theorem H5 and Corollaries H5a-H5b): the least sample SD
#' among integer sum-score profiles \eqn{S} with the given mean, integer values
#' in \eqn{[l, u]}, whose design-factor cap
#' \eqn{m_{\max}(S) = 2 n\, SS_S / \sum_{s,t}|S_s - S_t|} is at least the reported
#' design factor \eqn{m = 1/(1 - c\alpha)}. By inequality (\eqn{\star}) of the
#' article, every integer item set achieving the reported alpha has such a
#' profile, so this is a valid lower bound on the composite SD, and it is the
#' *sharp* one when the full logical composite window (spread \eqn{u - l}) is
#' enumerated. Unlike the amplified floor it is strictly positive whenever
#' \eqn{m > 1}, including at whole-number sum-score means (Corollary H5b).
#'
#' Enumeration spans the full composite range, so the returned floor is exact.
#' Returns `NULL` (caller falls back to the proven amplified floor) when no
#' integer composite exists (\eqn{n \cdot mean} not an integer) or the
#' enumeration `choose(n + W, W)`, \eqn{W = u - l}, would exceed `max_profiles`.
#'
#' @param l,u Integer scalars, sum-score limits.
#' @param n Integer scalar, sample size.
#' @param mean Numeric scalar, sum-score mean (with \eqn{n \cdot mean} an integer).
#' @param m Numeric scalar, the reported design factor \eqn{1/(1 - c\alpha)}, \eqn{\ge 1}.
#' @param max_profiles Numeric, enumeration budget (default 5e5).
#' @return Numeric scalar, the exact floor in sum-score SD units, or `NULL`.
#' @keywords internal
sd_min_alpha_gini <- function(l, u, n, mean, m, max_profiles = 5e5) {
  N <- round(n * mean)
  if (abs(n * mean - N) > 1e-9) {
    return(NULL)
  } # no integer composite
  target <- N - n * l # required sum of (value - l)
  env <- .alpha_gini_envelope(l, u, n, m, max_profiles)
  if (is.null(env)) {
    return(NULL)
  } # over budget -> fall back
  i <- match(target, env$target)
  if (is.na(i)) {
    return(NULL)
  }
  sqrt(env$ss[i] / (n - 1))
}

# Internal: the whole Gini envelope in ONE enumeration pass, memoized.
# The per-mean entry point above used to filter a full enumeration of the
# profile space down to a single sum, so sweeping every mean re-enumerated the
# same 10^5-10^6 profiles once per mean. Grouping by sum instead yields every
# mean's floor from a single pass, which is what makes mean-sweeps (figures,
# umbrella grids) affordable. Cached on (l, u, n, m, budget) because `m`
# decides which profiles can support the reported alpha.
.gini_envelope_cache <- new.env(parent = emptyenv())

.alpha_gini_envelope <- function(l, u, n, m, max_profiles = 5e5) {
  key <- paste(l, u, n, signif(m, 12), max_profiles, sep = "|")
  hit <- .gini_envelope_cache[[key]]
  if (!is.null(hit)) {
    return(if (identical(hit, NA)) NULL else hit)
  }
  W <- as.integer(round(u - l)) # composite spread
  if (
    W < 1L ||
      !is.finite(suppressWarnings(choose(n + W, W))) ||
      choose(n + W, W) > max_profiles
  ) {
    assign(key, NA, envir = .gini_envelope_cache)
    return(NULL)
  }
  vals <- 0:W # deviations above l
  cv <- .count_vectors(W + 1L, n)
  tot <- as.vector(cv %*% vals)
  SS <- as.vector(cv %*% (vals^2)) - tot^2 / n # SS_S (shift-invariant)
  # Gini double-sum over unordered value pairs: sum_{a<b} c_a c_b (v_b - v_a)
  G2 <- numeric(nrow(cv))
  for (a in seq_len(W)) {
    for (b in (a + 1L):(W + 1L)) {
      G2 <- G2 + cv[, a] * cv[, b] * (vals[b] - vals[a])
    }
  }
  mmax <- ifelse(G2 > 0, n * SS / G2, 0) # m_max(S) = n SS / (n V_min)
  ok <- mmax >= m - 1e-9 & SS > 1e-12 # non-constant, supports alpha
  res <- if (!any(ok)) {
    list(target = integer(0), ss = numeric(0))
  } else {
    agg <- tapply(SS[ok], tot[ok], min)
    list(target = as.integer(names(agg)), ss = as.numeric(agg))
  }
  assign(key, res, envir = .gini_envelope_cache)
  res
}

#' SD bounds for a k-item composite with reported Cronbach's alpha
#'
#' Bounds of the sum score's sample SD given the composite limits, n, sum-score
#' mean, granularity, and a reported alpha (the STRAIT article, Appendix H).
#' Ceiling: the mean-conditional ceiling divided by
#' `sqrt(k - (k - 1) alpha)`; under granularity, additionally the sharper
#' quasi-integer `V_max` form, intersected with the alpha-free ceiling (alpha
#' can only tighten). Floor under granularity: the alpha-amplified quasi-integer
#' floor `s_min / sqrt(1 - c alpha)`, sharpened for strictly integer composites
#' to the exact Gini envelope of [sd_min_alpha_gini()] (positive even at
#' whole-number means; Theorem H5) when that enumeration is affordable, else 0.
#' The floor can exceed the ceiling for small k (a genuinely infeasible
#' constraint set); this is reported rather than clipped.
#'
#' @param l,u Numeric scalars, limits of the SUM score (composite units).
#' @param n Integer scalar, sample size.
#' @param mean Numeric scalar, sum-score mean.
#' @param Z One of "continuous", "integer", "quasiinteger" (item-level grid).
#' @param alpha Numeric scalar, reported Cronbach's alpha, `alpha < 1`.
#' @param k_items Integer scalar, number of items, `k_items >= 1`.
#' @return List with `min_sd`, `max_sd`, `feasible`, `min_rule`, `note`.
#' @examples
#' # a 3-item 1-5 battery, so the sum score runs 3-15
#' sd_bounds_alpha(l = 3, u = 15, n = 50, mean = 9.4, Z = "continuous",
#'                 alpha = 0.8, k_items = 3)
#'
#' # integer items lift the floor off zero and tighten the ceiling
#' sd_bounds_alpha(l = 3, u = 15, n = 50, mean = 9.4, Z = "integer",
#'                 alpha = 0.8, k_items = 3)
#' @export
sd_bounds_alpha <- function(l, u, n, mean, Z, alpha, k_items) {
  k <- k_items
  cc <- (k - 1) / k
  D <- 1 - cc * alpha # k - (k-1) alpha = k * D
  if (D <= 1e-12) {
    stop("alpha too high for this k: 1 - alpha*(k-1)/k must be positive")
  }

  ceil_free <- sd_max_structure_s(mean, n, l, u) # alpha-free sharp ceiling
  ceil_smooth <- bessel_factor(n) *
    sqrt(pmax(0, (u - mean) * (mean - l)) / (k * D))
  ceiling <- min(ceil_free, ceil_smooth)
  floor_ <- 0
  floor_rule <- "s >= 0"
  note <- NA_character_
  if (Z %in% c("integer", "quasiinteger")) {
    item_l <- l / k
    item_u <- u / k
    ceil_vmax <- sqrt(
      (n / (n - 1)) * v_max_alpha(mean, k, n, item_l, item_u) / D
    )
    ceiling <- min(ceiling, ceil_vmax)
    floor_ <- sd_min_quasi_integer(mean, n) / sqrt(D) # proven amplified floor
    floor_rule <- "alpha-amplified quasi-integer floor"
    # sharpen with the exact Gini envelope for strictly integer composites; it
    # is >= the amplified floor and strictly positive at whole-number means.
    if (Z == "integer" && alpha > 0) {
      env <- sd_min_alpha_gini(l, u, n, mean, m = 1 / D)
      if (is.null(env)) {
        note <- "exact alpha floor (Theorem H5) not evaluated (n*mean non-integer or composite window over budget); proven amplified floor used, conservative near whole-number means"
      } else if (env > floor_ + 1e-12) {
        floor_ <- env
        floor_rule <- "alpha-conditional Gini envelope floor (Theorem H5)"
      }
    }
  }
  list(
    min_sd = floor_,
    max_sd = ceiling,
    feasible = floor_ <= ceiling + 1e-9,
    min_rule = floor_rule,
    note = note
  )
}

# ---- Layer 3: rounding / truncation of reported inputs -----------------------

# The rounding methods strait accepts, named as in scrutiny and roundwork.
# The threshold methods ("up_from", ...) are left out: strait has no
# `threshold` argument, and at its default they equal "up", "down" and
# "up_or_down".
# fmt: skip
ROUNDING_METHODS <- c(
  "up_or_down", "ceiling_or_floor", "up", "down", "even",
  "ceiling", "floor", "trunc", "anti_trunc",
  "ties_up", "ties_down", "ties_away", "ties_zero", "ties_even"
)

# Internal: check a `rounding` value against `allowed`, pointing the names
# strait used before 0.5.0 to their scrutiny equivalents.
.check_rounding <- function(
  rounding,
  allowed = ROUNDING_METHODS,
  several = FALSE
) {
  old <- c(half_up = "up", half_down = "down", native = "even")
  if (
    !is.character(rounding) ||
      !length(rounding) ||
      (!several && length(rounding) != 1L)
  ) {
    stop(
      "`rounding` must be ",
      if (several) {
        "a character vector"
      } else {
        "a single string"
      },
      " naming a rounding method",
      call. = FALSE
    )
  }
  bad <- setdiff(rounding, allowed)
  if (length(bad)) {
    hint <- intersect(bad, names(old))
    stop(
      "unknown rounding method(s): ",
      paste(bad, collapse = ", "),
      ". ",
      if (length(hint)) {
        paste0(
          "strait now uses scrutiny's names: ",
          paste0('"', hint, '" is "', old[hint], '"', collapse = ", "),
          ". "
        )
      },
      "Choose from: ",
      paste(allowed, collapse = ", "),
      call. = FALSE
    )
  }
  rounding
}

# Internal: the interval of exact values behind reported values `x` with
# `digits` decimal places, from roundwork::unround(), so that strait unrounds
# exactly as scrutiny does. Vectorised over `x`. `arg` and `arg_digits` name
# the arguments in errors.
.unround <- function(x, digits, rounding, arg = "x", arg_digits = "digits") {
  if (any(abs(x - round(x, digits)) > 1e-9 * pmax(1, abs(x)))) {
    stop(
      sprintf(
        "`%s` has more decimal places than `%s` states (%s)",
        arg,
        arg_digits,
        paste(digits, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  r <- roundwork::unround(
    round(x, digits),
    rounding = rounding,
    digits = digits
  )
  list(
    lo = r$lower,
    hi = r$upper,
    lo_incl = r$incl_lower,
    hi_incl = r$incl_upper
  )
}
