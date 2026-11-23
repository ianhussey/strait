# Plausibility bounds for the sample SD of bounded integer scores.
#
# Hard bounds (sd_min_quasi_integer(), sd_max_structure_s()) say which SDs are
# arithmetically possible. Plausibility bounds say which SDs a realistic
# data-generating process would produce: here, scores X in {0, ..., m},
# m = u - l, drawn iid from BetaBinomial(m, mu, rho), with mu = alpha /
# (alpha + beta) and rho = 1 / (alpha + beta + 1). The bounds are quantiles of
# the sampling distribution of the sample SD S, obtained by one of:
#
#   chisq2  two-moment (Box / Satterthwaite) scaled chi-square for S^2
#   chisq3  three-moment shifted chi-square (Pearson type III) for S^2
#   cf      Cornish-Fisher expansion for S^2, using its third cumulant
#   exact   the exact distribution of S^2, by dynamic programming over the
#           joint (sum X, sum X^2) state space (small N only)
#
# All population cumulants are computed from the pmf, not closed forms. With
# between-study heterogeneity tau on the rho scale, the bound is the quantile
# of the mixture over rho ~ N(rho, tau^2) truncated to [0, 1): the per-rho
# CDFs are averaged, then inverted. See paper/sd-plausibility-bounds/.

# ---- model -------------------------------------------------------------------

# Internal: beta-binomial pmf on 0:m. rho <= 0 is the binomial (the floor of
# the family); rho >= 1 is the limit with all mass at the two endpoints.
bb_pmf <- function(m, mu, rho) {
  x <- 0:m
  if (rho <= 1e-10) {
    return(stats::dbinom(x, m, mu))
  }
  if (rho >= 1) {
    return(c(1 - mu, rep(0, m - 1), mu))
  }
  phi <- 1 / rho - 1
  a <- mu * phi
  b <- (1 - mu) * phi
  p <- exp(lchoose(m, x) + lbeta(x + a, m - x + b) - lbeta(a, b))
  p / sum(p)
}

# Internal: population cumulants of X from the pmf.
cumulants_from_pmf <- function(pmf) {
  x <- seq_along(pmf) - 1
  mean <- sum(pmf * x)
  central_moments <- purrr::map_dbl(2:6, function(order) {
    sum(pmf * (x - mean)^order)
  })
  m2 <- central_moments[1]
  m3 <- central_moments[2]
  m4 <- central_moments[3]
  m5 <- central_moments[4]
  m6 <- central_moments[5]
  k4 <- m4 - 3 * m2^2
  list(
    mean = mean,
    k2 = m2,
    k3 = m3,
    k4 = k4,
    k5 = m5 - 10 * m3 * m2,
    k6 = m6 - 15 * m4 * m2 - 10 * m3^2 + 30 * m2^3,
    g1 = m3 / m2^1.5,
    g2 = k4 / m2^2
  )
}

#' Beta-binomial cumulants
#'
#' Population mean, cumulants \eqn{\kappa_2, \dots, \kappa_6}, skewness
#' \eqn{\gamma_1} and excess kurtosis \eqn{\gamma_2} of
#' BetaBinomial\eqn{(m, \mu, \rho)}, computed exactly from the pmf
#' (\eqn{\log p(x) = \log\binom{m}{x} + \log B(x + \alpha, m - x + \beta) -
#' \log B(\alpha, \beta)}) rather than from closed forms.
#'
#' @param m Whole number of scale steps, `u - l`.
#' @param mu Mean proportion in \eqn{(0, 1)}.
#' @param rho Intra-class correlation in \eqn{[0, 1]}. \eqn{\rho \le 0} is
#'   the binomial; \eqn{\rho = 1} puts all mass on the two endpoints.
#' @return A one-row tibble with `mean`, `k2` to `k6`, `g1` and `g2`.
#' @examples
#' bb_cumulants(63, 25 / 63, 0.052)
#' @export
bb_cumulants <- function(m, mu, rho) {
  tibble::as_tibble(cumulants_from_pmf(bb_pmf(m, mu, rho)))
}

# Internal: cumulants of the unbiased sample variance S^2 of N iid draws.
# conditional = TRUE gives the leading-order moments of S^2 given that the
# sample mean equals the population mean, removing the part of Var(S^2) that
# Cov(xbar, S^2) = kappa_3 / N explains.
s2_cumulants <- function(cumulants, N, conditional = FALSE) {
  k2 <- cumulants$k2
  k3 <- cumulants$k3
  k4 <- cumulants$k4
  k6 <- cumulants$k6
  variance <- k2^2 * (2 / (N - 1) + cumulants$g2 / N)
  if (conditional) {
    variance <- variance - k3^2 / (N * k2)
  }
  third <- k6 /
    N^2 +
    12 * k2 * k4 / (N * (N - 1)) +
    4 * (N - 2) * k3^2 / (N * (N - 1)^2) +
    8 * k2^3 / (N - 1)^2
  list(k1 = k2, k2 = variance, k3 = third)
}

# ---- approximate CDFs of S ---------------------------------------------------

# Internal: CDF of S at s for one rho, under an approximation of S^2.
cdf_s_approx <- function(s, cumulants, method) {
  s2 <- s^2
  switch(
    method,
    chisq2 = {
      nu <- 2 * cumulants$k1^2 / cumulants$k2
      stats::pchisq(s2 * nu / cumulants$k1, nu)
    },
    chisq3 = {
      if (cumulants$k3 <= 0) {
        return(cdf_s_approx(s, cumulants, "chisq2"))
      }
      b <- cumulants$k3 / (4 * cumulants$k2)
      nu <- 8 * cumulants$k2^3 / cumulants$k3^2
      a <- cumulants$k1 - b * nu
      stats::pchisq(pmax(0, (s2 - a) / b), nu)
    },
    cf = {
      # invert the Cornish-Fisher quantile k1 + sqrt(k2) (z + g/6 (z^2 - 1))
      # for z; the + root is the monotone branch
      standardized <- (s2 - cumulants$k1) / sqrt(cumulants$k2)
      skew <- cumulants$k3 / cumulants$k2^1.5 / 6
      if (abs(skew) < 1e-8) {
        return(stats::pnorm(standardized))
      }
      discriminant <- 1 + 4 * skew * (skew + standardized)
      # `s` is a scalar from uniroot(), so a scalar `if` (an NA discriminant
      # still gives NA) avoids vectorised-branch overhead in this hot loop
      z <- if (isTRUE(discriminant < 0)) {
        if (skew > 0) -Inf else Inf
      } else {
        (-1 + sqrt(pmax(discriminant, 0))) / (2 * skew)
      }
      stats::pnorm(z)
    },
    cli::cli_abort("Unknown method: {method}.", .internal = TRUE)
  )
}

# ---- exact distribution by dynamic programming -------------------------------

# Internal: minimum and maximum of sum x^2 over n integers in [0, m] summing
# to s.
q_min <- function(n, s) {
  if (n == 0) {
    return(0)
  }
  a <- s %/% n
  remainder <- s %% n
  (n - remainder) * a^2 + remainder * (a + 1)^2
}
q_max <- function(n, s, m) {
  k <- s %/% m
  remainder <- s %% m
  k * m^2 + remainder^2
}

# Internal: exact joint distribution of (sum X, sum X^2) for N iid draws from
# the pmf `pr` on 0:m, by dynamic programming over observations. State: for
# each sum s, a vector over q in [q_min(n, s), q_max(n, s)]. Values with
# p(x) < prune and sums whose total mass is below prune are dropped; the
# discarded mass is returned as `lost`. `keep` lists the sample sizes whose
# joint distribution is returned (all <= N), so one run serves several N.
sq_joint <- function(pmf, N, keep = N, prune = 1e-17, max_cells = 4e7) {
  m <- length(pmf) - 1
  # vapply() rather than purrr: this is the exact method's hot path
  cells <- sum(vapply(
    0:(N * m),
    function(s) q_max(N, s, m) - q_min(N, s) + 1,
    numeric(1)
  ))
  if (cells > max_cells) {
    cells_shown <- format(cells, big.mark = ",")
    max_cells_shown <- format(max_cells, big.mark = ",")
    cli::cli_abort(
      "The exact method needs about {cells_shown} cells for n = {N} on a
      {m}-step scale, more than max_cells = {max_cells_shown}. Use an
      approximate method for this n.",
      call = NULL
    )
  }
  support <- which(pmf >= prune) - 1
  state <- list(1) # n = 0: s = 0, q = 0
  lost <- 0
  out <- list()
  for (n in 0:(N - 1)) {
    sum_max_next <- (n + 1) * m
    state_next <- vector("list", sum_max_next + 1)
    for (sum_next in 0:sum_max_next) {
      base <- q_min(n + 1, sum_next)
      probs_next <- NULL
      for (x in support) {
        s <- sum_next - x
        if (s < 0 || s > n * m) {
          next
        }
        probs <- state[[s + 1]]
        if (is.null(probs)) {
          next
        }
        if (is.null(probs_next)) {
          probs_next <- numeric(q_max(n + 1, sum_next, m) - base + 1)
        }
        index <- q_min(n, s) + x^2 - base + seq_along(probs)
        probs_next[index] <- probs_next[index] + pmf[x + 1] * probs
      }
      if (!is.null(probs_next)) {
        total <- sum(probs_next)
        if (total < prune) {
          lost <- lost + total
        } else {
          state_next[[sum_next + 1]] <- probs_next
        }
      }
    }
    state <- state_next
    if ((n + 1) %in% keep) {
      out[[as.character(n + 1)]] <- list(state = state, N = n + 1, m = m)
    }
  }
  structure(out, lost = lost + (1 - sum(pmf[support + 1])) * N)
}

# Internal: distribution of S^2 from a joint state, optionally restricted to
# the sums in `sums` (conditioning on the sample mean). Returns a tibble
# (s2, prob) sorted by s2, with prob normalised within the condition.
s2_from_joint <- function(joint, sums = NULL) {
  N <- joint$N
  state <- joint$state
  # vapply() rather than purrr: part of the exact method's hot path
  sums_attained <- which(!vapply(state, is.null, logical(1))) - 1
  if (!is.null(sums)) {
    sums_attained <- intersect(sums_attained, sums)
  }
  if (length(sums_attained) == 0) {
    cli::cli_abort(
      "No attainable sample has a mean in the conditioning set.",
      call = NULL
    )
  }
  scaled_s2 <- list()
  probs <- list()
  for (s in sums_attained) {
    probs_current <- state[[s + 1]]
    nonzero <- which(probs_current > 0)
    q <- q_min(N, s) + nonzero - 1
    scaled_s2[[length(scaled_s2) + 1]] <- N * q - s^2 # N (N - 1) S^2, an integer
    probs[[length(probs) + 1]] <- probs_current[nonzero]
  }
  totals <- rowsum(unlist(probs), unlist(scaled_s2), reorder = TRUE)
  tibble::tibble(
    s2 = as.numeric(rownames(totals)) / (N * (N - 1)),
    prob = unname(totals[, 1] / sum(totals[, 1]))
  )
}

#' Exact sampling distribution of the sample variance of beta-binomial scores
#'
#' The exact distribution of the unbiased sample variance \eqn{S^2} of \eqn{n}
#' iid BetaBinomial\eqn{(m, \mu, \rho)} scores, by dynamic programming over
#' the joint state \eqn{(\sum X, \sum X^2)}. Feasible because \eqn{X} is a
#' bounded integer; cost grows roughly as \eqn{n^2 m^3}, so this is for small
#' \eqn{n} (about \eqn{n \le 25} on a 0-63 scale).
#'
#' @inheritParams bb_cumulants
#' @param n Sample size(s), each a whole number \eqn{\ge 2}. A vector is
#'   computed in one pass up to `max(n)`.
#' @param sums Optional integer vector of attainable sums \eqn{\sum X}
#'   (on the \eqn{0..m} scale) to condition on, i.e. conditioning on the
#'   sample mean. Applied to every `n`.
#' @return A named list (by `n`) of tibbles with columns `s2` and `prob`,
#'   sorted by `s2`. The attribute `lost` is the probability mass discarded
#'   by pruning (values and sums below \eqn{10^{-17}}).
#' @examples
#' d <- s2_exact(m = 4, mu = 0.4, rho = 0.1, n = 5)[["5"]]
#' sum(d$prob)
#' sum(d$s2 * d$prob)   # equals the population variance
#' @export
s2_exact <- function(m, mu, rho, n, sums = NULL) {
  plausibility_check_n(n, allow_na = FALSE)
  joint <- sq_joint(bb_pmf(m, mu, rho), max(n), keep = n)
  out <- purrr::map(joint, function(joint_current) {
    s2_from_joint(joint_current, sums)
  })
  attr(out, "lost") <- attr(joint, "lost")
  out
}

# Internal: stop unless n holds whole numbers >= 2 (NA allowed if allow_na).
plausibility_check_n <- function(n, allow_na = TRUE) {
  values <- if (allow_na) n[!is.na(n)] else n
  if (
    !is.numeric(n) ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < 2) ||
      any(values != round(values))
  ) {
    cli::cli_abort(
      "{.arg n} must contain finite whole numbers >= 2.",
      call = NULL
    )
  }
}
# ---- hard bounds for a reported mean -----------------------------------------

# Internal: sums t (on the 0..m scale) whose mean l + t / n rounds to `mean`
# at `digits` decimals (half-width tolerance), or the exact sum when digits
# is NULL and the mean is GRIM-consistent. `digits` is a single value.
consistent_sums <- function(l, mean, n, m, digits = NULL) {
  if (is.null(digits)) {
    t <- round(n * (mean - l))
    if (abs(t - n * (mean - l)) > 1e-7) {
      return(integer(0))
    }
    return(as.integer(t))
  }
  half_width <- 0.5 * 10^-digits + TOLERANCE
  lowest <- ceiling(n * (mean - l - half_width))
  highest <- floor(n * (mean - l + half_width))
  if (highest < lowest) {
    return(integer(0))
  }
  t <- lowest:highest
  t[t >= 0 & t <= n * m]
}

# Internal: hard min / max SD for integer data on [l, u] with reported mean.
hard_bounds <- function(l, u, mean, n, m, digits) {
  t <- consistent_sums(l, mean, n, m, digits)
  if (length(t) == 0) {
    # no integer sample has this mean: fall back to the quasi-integer floor,
    # which is a valid floor at every real mean
    return(c(
      sd_min_quasi_integer(mean - l, n),
      sd_max_structure_s(mean, n, l, u)
    ))
  }
  means <- l + t / n
  c(min(sd_min_integer(means - l, n)), max(sd_max_structure_s(means, n, l, u)))
}

# ---- the main function -------------------------------------------------------

# Internal: rho nodes of N(rho, tau^2) truncated to [0, 1), equal weights.
rho_nodes <- function(rho, tau, K) {
  if (tau <= 0) {
    return(rho)
  }
  p_lower <- stats::pnorm(0, rho, tau)
  p_upper <- stats::pnorm(1, rho, tau)
  u <- (seq_len(K) - 0.5) / K
  stats::qnorm(p_lower + (p_upper - p_lower) * u, rho, tau) |>
    pmin(1 - TOLERANCE)
}

# Internal: lower and upper quantiles of S for one row, at every tail
# probability in `p`. The distribution (exact) or the cumulants (approximate)
# are computed once per rho node and shared by all `p`.
sd_plausibility_row <- function(
  l,
  u,
  mean,
  n,
  p,
  rho,
  tau,
  method,
  conditional,
  digits,
  K
) {
  m <- u - l
  mu <- (mean - l) / m
  rhos <- rho_nodes(rho, tau, K)
  if (method == "exact") {
    sums <- if (conditional) {
      t <- consistent_sums(l, mean, n, m, digits)
      if (length(t) == 0) {
        cli::cli_abort(
          "Conditioning on the mean with the exact method needs a
          GRIM-consistent mean, or {.arg digits_mean} to define the sums that
          round to it.",
          call = NULL
        )
      }
      t
    }
    distribution <- rhos |>
      purrr::map(function(rho_current) {
        joint <- sq_joint(bb_pmf(m, mu, rho_current), n)[[1]]
        s2_from_joint(joint, sums)
      }) |>
      dplyr::bind_rows()
    totals <- rowsum(
      distribution$prob / length(rhos),
      round(distribution$s2 * n * (n - 1)),
      reorder = TRUE
    )
    s <- sqrt(as.numeric(rownames(totals)) / (n * (n - 1)))
    cdf <- cumsum(totals[, 1])
    q <- function(pp) s[which(cdf >= pp - TOLERANCE_TIGHT)[1]]
  } else {
    s2_cumulants_all <- purrr::map(rhos, function(rho_current) {
      bb_pmf(m, mu, rho_current) |>
        cumulants_from_pmf() |>
        s2_cumulants(n, conditional)
    })
    # vapply() rather than purrr: uniroot() calls this many times per bound
    cdf <- function(s) {
      cdf_nodes <- vapply(
        s2_cumulants_all,
        function(cumulants) cdf_s_approx(s, cumulants, method),
        numeric(1)
      )
      mean(cdf_nodes)
    }
    upper_start <- 2 * sqrt(max(purrr::map_dbl(s2_cumulants_all, "k1"))) + 1
    q <- function(pp) {
      if (cdf(0) >= pp) {
        return(0)
      }
      stats::uniroot(
        function(s) cdf(s) - pp,
        c(0, upper_start),
        extendInt = "upX",
        tol = 1e-10
      )$root
    }
  }
  list(lower = purrr::map_dbl(p, q), upper = purrr::map_dbl(1 - p, q))
}

#' Plausibility bounds for a reported SD on a bounded integer scale
#'
#' Quantiles of the sampling distribution of the sample SD under a
#' beta-binomial data-generating process: the range of SDs that samples of
#' size `n` would produce at mean `mean` if scores were
#' BetaBinomial\eqn{(u - l, \mu, \rho)} with \eqn{\mu = (mean - l)/(u - l)}.
#' Returned next to the hard (arithmetic) bounds for integer data.
#'
#' @section Model:
#' \eqn{X = \text{score} - l \in \{0, \dots, m\}}, \eqn{m = u - l}, with
#' \eqn{E[X] = m\mu} and \eqn{\sigma^2 = m\mu(1 - \mu)[1 + (m - 1)\rho]}.
#' \eqn{\rho = 0} is the binomial, the least dispersed member of the family.
#' All moments are computed from the pmf. On a two-point scale (\eqn{m = 1})
#' scores are Bernoulli\eqn{(\mu)} whatever \eqn{\rho}, so \eqn{\rho} is not
#' identified and is returned as `NA`.
#'
#' @section Methods:
#' For iid data with finite fourth moment, \eqn{E[S^2] = \sigma^2} and
#' \eqn{\mathrm{Var}(S^2) = \sigma^4 [2/(n-1) + \gamma_2/n]}.
#' * `"chisq2"` (default): \eqn{S^2 \approx \sigma^2 \chi^2_\nu / \nu} with
#'   \eqn{\nu = 2/[2/(n-1) + \gamma_2/n]} (Box / Satterthwaite).
#' * `"chisq3"`: \eqn{S^2 \approx a + b\chi^2_\nu}, matching the first three
#'   cumulants of \eqn{S^2}.
#' * `"cf"`: Cornish-Fisher expansion of \eqn{S^2} using its third cumulant
#'   \eqn{\kappa_3(S^2) = \kappa_6/n^2 + 12\kappa_2\kappa_4/(n(n-1)) +
#'   4(n-2)\kappa_3^2/(n(n-1)^2) + 8\kappa_2^3/(n-1)^2}.
#' * `"exact"`: the exact distribution by dynamic programming (see
#'   [s2_exact()]); small `n` only.
#'
#' @section Accuracy:
#' The approximations are least accurate in the lower tail of a skewed score
#' distribution at small `n`. There the `"chisq3"` and `"cf"` lower bounds
#' can be anti-conservative: at \eqn{m = 63}, `mean = 2`, `n = 30`,
#' \eqn{\rho = 0.1} (unconditional), about 2.6% of samples fall below a
#' nominal 0.5% lower bound. For \eqn{n \le 20}, use `method = "exact"`; above that, `"chisq2"`
#' with `condition_on_mean = TRUE`, whose upper bound is reliable. See
#' `paper/sd-plausibility-bounds/` in the package sources.
#'
#' @section Conditioning on the mean:
#' A reported SD comes with a reported mean. When the distribution is
#' skewed, \eqn{S^2} and \eqn{\bar X} are correlated
#' (\eqn{\mathrm{Cov}(\bar X, S^2) = \kappa_3/n}), so samples whose mean is
#' near `mean` have a narrower SD distribution than samples in general.
#' `condition_on_mean = TRUE` (the default) asks the question a check of a
#' reported (mean, SD) pair asks: `"chisq2"` uses the conditional variance
#' \eqn{\sigma^4[2/(n-1) + (\gamma_2 - \gamma_1^2)/n]}, and `"exact"`
#' conditions on the sums whose mean rounds to `mean` (at `digits_mean`
#' decimals). `"chisq3"` and `"cf"` need `condition_on_mean = FALSE`.
#'
#' @section Heterogeneity:
#' With `tau > 0`, \eqn{\rho} is drawn from \eqn{N(\rho, \tau^2)} truncated to
#' \eqn{[0, 1)}, represented by `K` equal-weight quantile nodes. The bound is
#' the quantile of the mixture: the per-node CDFs are averaged, then
#' inverted. When `rho` and `tau` come from [rho_corpus_fit()], consider
#' adding the prediction's squared standard error to `tau^2`.
#'
#' @param l,u Scale limits (scalars). `u - l` must be a whole number: the
#'   model is for integer scores. Rescale mean-scored scales to sum scores.
#' @param n Sample size(s), whole numbers \eqn{\ge 2}.
#' @param mean Reported mean(s).
#' @param sd Reported SD(s), used to set `rho` by the method of moments.
#'   This is circular for checking that same SD and is meant for
#'   illustration. Give `sd` or `rho`, not both.
#' @param ... Must be empty. Arguments after the dots must be named in full.
#' @param rho Intra-class correlation of the data-generating process, e.g.
#'   a corpus prediction from [rho_corpus_fit()]. Values below 0 are
#'   under-dispersed relative to the binomial and are set to 0 (with a note);
#'   values of 1 or more give the endpoint-only limit.
#' @param tau Between-study SD of `rho`, \eqn{\ge 0}. Default 0.
#' @param p Tail probability in **each** tail, in (0, 0.5); several values
#'   give several rows. Default 0.01.
#' @param method `"chisq2"`, `"chisq3"`, `"cf"` or `"exact"`.
#' @param condition_on_mean Logical; see the section above. Default `TRUE`.
#' @param digits_mean Decimal places of the reported mean, if known. Used for
#'   the hard bounds (minimum and maximum over every sum that rounds to
#'   `mean`) and for exact conditioning. `NULL` or `NA` treats `mean` as
#'   exact.
#' @param K Number of quadrature nodes for the \eqn{\rho} mixture.
#'
#' @return A tibble with one row per (`mean`, `n`, `sd`, `rho`, `tau`,
#'   `digits_mean`) element (recycled together; each of length 1 or a common
#'   length) and `p` value (crossed): the inputs, `mu`, `rho`, `sigma`
#'   (population SD at `rho`), `hard_min` and `hard_max` (integer data),
#'   `lower` and `upper` (the plausibility bounds, clipped to the hard
#'   bounds), `clipped` (which side, if any, was clipped) and `note`. A
#'   missing input or a mean outside the scale gives `NA` bounds with a note.
#' @seealso [sd_min_integer()], [sd_max_structure_s()] for the hard bounds;
#'   [rho_from_summary()], [rho_corpus_fit()] for estimating `rho`.
#' @examples
#' # BDI-II (0-63), M = 25, SD = 8 sets rho by the method of moments
#' sd_plausibility(0, 63, n = c(20, 50, 100, 500), mean = 25, sd = 8,
#'                 p = c(0.01, 0.001))
#'
#' # a corpus prediction with heterogeneity
#' sd_plausibility(0, 63, n = 50, mean = 25, rho = 0.05, tau = 0.02)
#' @export
sd_plausibility <- function(
  l,
  u,
  n,
  mean,
  sd = NULL,
  ...,
  rho = NULL,
  tau = 0,
  p = 0.01,
  method = c("chisq2", "chisq3", "cf", "exact"),
  condition_on_mean = TRUE,
  digits_mean = NULL,
  K = 41
) {
  grim_aggregate_check_dots(..., .renamed = sd_renamed)
  method <- rlang::arg_match(method)
  m <- plausibility_check_scale(l, u)
  if (is.null(rho) == is.null(sd)) {
    cli::cli_abort(
      "Supply exactly one of {.arg rho} (e.g. a corpus prediction) and
      {.arg sd} (method of moments; circular for checking that SD)."
    )
  }
  if (
    !is.logical(condition_on_mean) ||
      length(condition_on_mean) != 1 ||
      is.na(condition_on_mean)
  ) {
    cli::cli_abort("{.arg condition_on_mean} must be TRUE or FALSE.")
  }
  if (condition_on_mean && method %in% c("chisq3", "cf")) {
    cli::cli_abort(
      "Methods {.val chisq3} and {.val cf} cannot condition on the mean: set
      {.code condition_on_mean = FALSE}, or use {.val chisq2} or
      {.val exact}."
    )
  }
  plausibility_check_n(n)
  if (!is.numeric(p) || length(p) == 0 || anyNA(p) || any(p <= 0 | p >= 0.5)) {
    cli::cli_abort("{.arg p} must contain tail probabilities in (0, 0.5).")
  }
  if (!is.numeric(tau) || any(tau < 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg tau} must be non-negative.")
  }
  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K < 1 || K != round(K)) {
    cli::cli_abort("{.arg K} must be a whole number >= 1.")
  }
  if (!is.null(digits_mean)) {
    digits_known <- digits_mean[!is.na(digits_mean)]
    if (
      !is.numeric(digits_mean) ||
        any(digits_known < 0) ||
        any(digits_known != round(digits_known))
    ) {
      cli::cli_abort(
        "{.arg digits_mean} must contain whole numbers >= 0 (or NA)."
      )
    }
  }
  args <- list(
    mean = mean,
    n = n,
    sd = sd,
    rho = rho,
    tau = tau,
    digits_mean = digits_mean
  )
  args <- purrr::discard(args, is.null)
  lengths_all <- lengths(args)
  size <- max(lengths_all)
  if (any(!lengths_all %in% c(1L, size))) {
    cli::cli_abort(
      "{.arg mean}, {.arg n}, {.arg sd}, {.arg rho}, {.arg tau} and
      {.arg digits_mean} must each have length 1 or a common length."
    )
  }
  mean <- rep_len(mean, size)
  n <- rep_len(n, size)
  tau <- rep_len(tau, size)
  rho_in <- if (is.null(rho)) rep(NA_real_, size) else rep_len(rho, size)
  sd_in <- if (is.null(sd)) rep(NA_real_, size) else rep_len(sd, size)
  digits_in <- if (is.null(digits_mean)) {
    rep(NA_real_, size)
  } else {
    rep_len(digits_mean, size)
  }

  rows <- vector("list", size)
  for (i in seq_len(size)) {
    mu <- (mean[i] - l) / m
    row <- list(
      l = l,
      u = u,
      mean = mean[i],
      n = n[i],
      p = p,
      method = method,
      conditional = condition_on_mean,
      mu = mu,
      rho = rho_in[i],
      tau = tau[i],
      sigma = NA_real_,
      hard_min = NA_real_,
      hard_max = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      clipped = NA_character_,
      note = NA_character_
    )
    note <- character(0)
    if (
      anyNA(c(mean[i], n[i], tau[i])) ||
        (is.na(rho_in[i]) && is.na(sd_in[i]))
    ) {
      row$note <- "missing input"
    } else if (mean[i] < l || mean[i] > u) {
      row$note <- "mean outside [l, u]"
    } else if (mean[i] == l || mean[i] == u) {
      row[c("sigma", "hard_min", "hard_max", "lower", "upper")] <- 0
      row$clipped <- "none"
      row$note <- "mean at a scale limit: every score equals it"
    } else {
      if (m == 1) {
        # Bernoulli(mu) whatever rho: the pmf ignores it
        note <- "u - l = 1: scores are binary, so rho is not identified"
        rho_i <- 0
        tau_i <- 0
        row$rho <- NA_real_
      } else {
        rho_i <- if (!is.na(rho_in[i])) {
          rho_in[i]
        } else {
          (sd_in[i]^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
        }
        tau_i <- tau[i]
        if (rho_i < 0) {
          note <- sprintf(
            "rho = %.4g < 0 is under-dispersed; bounds use rho = 0",
            rho_i
          )
          rho_i <- 0
        }
        if (rho_i >= 1) {
          note <- "rho >= 1: all mass at the scale endpoints"
          rho_i <- 1
        }
        row$rho <- rho_i
      }
      digits_i <- if (is.na(digits_in[i])) NULL else digits_in[i]
      hard <- hard_bounds(l, u, mean[i], n[i], m, digits_i)
      bounds <- sd_plausibility_row(
        l,
        u,
        mean[i],
        n[i],
        p,
        rho_i,
        tau_i,
        method,
        condition_on_mean,
        digits_i,
        K
      )
      clip_lower <- bounds$lower < hard[1] - TOLERANCE_TIGHT
      clip_upper <- bounds$upper > hard[2] + TOLERANCE_TIGHT
      row$sigma <- sqrt(cumulants_from_pmf(bb_pmf(m, mu, rho_i))$k2)
      row$hard_min <- hard[1]
      row$hard_max <- hard[2]
      row$lower <- pmax(bounds$lower, hard[1])
      row$upper <- pmin(bounds$upper, hard[2])
      # nested if_else() rather than case_when(): an NA comparison stays NA
      row$clipped <- dplyr::if_else(
        clip_lower & clip_upper,
        "lower+upper",
        dplyr::if_else(
          clip_lower,
          "lower",
          dplyr::if_else(clip_upper, "upper", "none")
        )
      )
      clip_note <- dplyr::if_else(
        row$clipped == "none",
        NA_character_,
        paste(
          "clipped to the hard",
          sub("+", " and ", row$clipped, fixed = TRUE),
          "bound"
        )
      )
      row$note <- purrr::map_chr(clip_note, function(clip_note_current) {
        all_notes <- c(note, if (!is.na(clip_note_current)) clip_note_current)
        if (length(all_notes)) {
          paste(all_notes, collapse = "; ")
        } else {
          NA_character_
        }
      })
    }
    rows[[i]] <- tibble::as_tibble(row)
  }
  dplyr::bind_rows(rows)
}

# Internal: old argument and column names, for clear errors.
sd_renamed <- c(
  L = "l",
  U = "u",
  M = "mean",
  N = "n",
  SD = "sd",
  digits = "digits_mean"
)

# Internal: validate scalar scale limits; return m = u - l.
plausibility_check_scale <- function(l, u) {
  if (
    !is.numeric(l) ||
      !is.numeric(u) ||
      length(l) != 1 ||
      length(u) != 1 ||
      anyNA(c(l, u)) ||
      !(u > l) ||
      !is.finite(u - l)
  ) {
    cli::cli_abort(
      "{.arg l} and {.arg u} must be single finite numbers with u > l.",
      call = NULL
    )
  }
  m <- u - l
  if (abs(m - round(m)) > TOLERANCE) {
    cli::cli_abort(
      "{.code u - l} must be a whole number: the beta-binomial model is for
      integer scores. For a mean-scored scale, convert to sum scores (multiply
      l, u, mean and sd by the number of items).",
      call = NULL
    )
  }
  round(m)
}

# ---- estimating rho from summary statistics ----------------------------------

#' Method-of-moments rho from summary statistics
#'
#' \eqn{\hat\mu = (mean - l)/m} and \eqn{\hat\rho = [sd^2/(m\hat\mu(1 -
#' \hat\mu)) - 1]/(m - 1)}. \eqn{\hat\rho} is linear in \eqn{S^2}, so it is
#' unbiased for \eqn{\rho} (given \eqn{\mu}), with sampling variance
#' \eqn{\mathrm{Var}(S^2)/[m(m - 1)\mu(1 - \mu)]^2}. Negative
#' \eqn{\hat\rho} are kept, not truncated.
#'
#' @param l,u Scale limits (scalars); `u - l` a whole number \eqn{\ge 2}
#'   (on a two-point scale \eqn{\rho} is not identified).
#' @param n,mean,sd Reported sample size, mean and SD (vectors, each of
#'   length 1 or a common length).
#' @param ... Must be empty.
#' @param rho_var Optional \eqn{\rho} at which to evaluate the kurtosis in
#'   the variance (default: \eqn{\max(0, \hat\rho)}). [rho_corpus_fit()]
#'   passes the shrunken corpus prediction here.
#' @return A tibble with `mu_hat`, `rho_hat` and `var_rho`.
#' @note The variance treats \eqn{\mu} as known. Its sampling error adds a
#'   term of order \eqn{1/n} that is small next to \eqn{\mathrm{Var}(S^2)}
#'   except near the scale limits.
#' @examples
#' rho_from_summary(0, 63, n = c(50, 120), mean = c(25, 12), sd = c(8, 7))
#' @export
rho_from_summary <- function(l, u, n, mean, sd, ..., rho_var = NULL) {
  grim_aggregate_check_dots(..., .renamed = sd_renamed)
  m <- plausibility_check_scale(l, u)
  if (m < 2) {
    cli::cli_abort(
      "{.code u - l} must be at least 2: on a two-point scale scores are
      binary and rho is not identified."
    )
  }
  lengths_all <- lengths(list(n, mean, sd))
  size <- max(lengths_all)
  if (any(!lengths_all %in% c(1L, size))) {
    cli::cli_abort(
      "{.arg n}, {.arg mean} and {.arg sd} must each have length 1 or a common
      length."
    )
  }
  n <- rep_len(n, size)
  mean <- rep_len(mean, size)
  sd <- rep_len(sd, size)
  mu <- (mean - l) / m
  rho_hat <- (sd^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
  rho_kurtosis <- if (is.null(rho_var)) {
    pmax(0, rho_hat)
  } else {
    rep_len(rho_var, size)
  }
  var_rho <- purrr::map_dbl(seq_len(size), function(i) {
    if (
      anyNA(c(mu[i], rho_kurtosis[i], n[i])) ||
        mu[i] <= 0 ||
        mu[i] >= 1 ||
        n[i] < 2
    ) {
      return(NA_real_)
    }
    rho_clamped <- min(max(rho_kurtosis[i], 0), 1)
    cumulants <- cumulants_from_pmf(bb_pmf(m, mu[i], rho_clamped))
    s2_cumulants(cumulants, n[i])$k2 / (m * (m - 1) * mu[i] * (1 - mu[i]))^2
  })
  tibble::tibble(mu_hat = mu, rho_hat = rho_hat, var_rho = var_rho)
}

# Internal: rho_from_summary() row by row, for corpora mixing scales.
rho_corpus_rows <- function(corpus, rho_var = rep(NA_real_, nrow(corpus))) {
  seq_len(nrow(corpus)) |>
    purrr::map(function(i) {
      rho_from_summary(
        corpus$l[i],
        corpus$u[i],
        n = corpus$n[i],
        mean = corpus$mean[i],
        sd = corpus$sd[i],
        rho_var = if (!is.na(rho_var[i])) rho_var[i]
      )
    }) |>
    dplyr::bind_rows()
}

#' Corpus meta-regression of rho on mu
#'
#' Fits a random-effects meta-regression of per-study method-of-moments
#' \eqn{\hat\rho_i} on \eqn{\hat\mu_i}, with a natural spline in
#' \eqn{\hat\mu}, from summary statistics alone (no raw data). The sampling
#' variance of \eqn{\hat\rho_i} depends on the kurtosis at \eqn{(\mu_i,
#' \rho_i)}; it is first evaluated at \eqn{\max(0, \hat\rho_i)}, then (with
#' `iterate = TRUE`) once more at the fitted corpus prediction, and the model
#' is refitted.
#'
#' @param corpus A data frame with columns `l`, `u`, `n`, `mean` and `sd`,
#'   one row per study. Scales may differ between rows.
#' @param df Degrees of freedom of the natural spline. Default 3.
#' @param iterate Logical; re-evaluate the variances at the prediction once.
#' @param ... Further arguments passed to [metafor::rma()].
#' @return An object of class `rho_corpus`: the `fit` ([metafor::rma()]),
#'   the spline `basis`, the per-study `data` (a tibble with `mu_hat`,
#'   `rho_hat`, `var_rho`) and `tau` (\eqn{\hat\tau}). Use [predict()] on it.
#' @seealso [predict.rho_corpus()], [sd_plausibility()]
#' @export
rho_corpus_fit <- function(corpus, df = 3, iterate = TRUE, ...) {
  rlang::check_installed("metafor", reason = "to fit the corpus model.")
  need <- c("l", "u", "n", "mean", "sd")
  if (!is.data.frame(corpus) || !all(need %in% names(corpus))) {
    hint <- if (
      is.data.frame(corpus) &&
        all(names(sd_renamed)[1:5] %in% names(corpus))
    ) {
      " Columns `L`, `U`, `M`, `SD` and `N` were renamed."
    } else {
      ""
    }
    cli::cli_abort(
      "{.arg corpus} must have columns {.arg l}, {.arg u}, {.arg n},
      {.arg mean} and {.arg sd}.{hint}"
    )
  }
  studies <- corpus |>
    tibble::as_tibble() |>
    dplyr::bind_cols(rho_corpus_rows(corpus))
  studies <- studies[
    stats::complete.cases(studies[c("mu_hat", "rho_hat", "var_rho")]),
  ]
  basis <- splines::ns(studies$mu_hat, df = df)
  fit_once <- function(var_rho) {
    metafor::rma(yi = studies$rho_hat, vi = var_rho, mods = ~basis, ...)
  }
  fit <- fit_once(studies$var_rho)
  if (iterate) {
    prediction <- as.numeric(stats::predict(fit)$pred)
    studies$var_rho <- rho_corpus_rows(studies, rho_var = prediction)$var_rho
    fit <- fit_once(studies$var_rho)
  }
  structure(
    list(
      fit = fit,
      basis = basis,
      data = studies,
      tau = sqrt(fit$tau2)
    ),
    class = "rho_corpus"
  )
}

#' Predict rho and tau from a corpus fit
#'
#' @param object A `rho_corpus` object from [rho_corpus_fit()].
#' @param mu Mean proportions \eqn{(mean - l)/(u - l)} at which to predict.
#' @param ... Unused.
#' @return A tibble with `mu`, `rho` (predicted mean), `se` (its standard
#'   error), `tau` (between-study SD) and `tau_total` \eqn{=\sqrt{\tau^2 +
#'   se^2}}, the spread to pass as `tau` to [sd_plausibility()] when the
#'   prediction's own uncertainty should count.
#' @export
predict.rho_corpus <- function(object, mu, ...) {
  mods <- stats::predict(object$basis, mu)
  prediction <- stats::predict(object$fit, newmods = mods)
  se <- as.numeric(prediction$se)
  tibble::tibble(
    mu = mu,
    rho = as.numeric(prediction$pred),
    se = se,
    tau = object$tau,
    tau_total = sqrt(object$tau^2 + se^2)
  )
}

#' @export
print.rho_corpus <- function(x, ...) {
  cat(
    "Corpus meta-regression of rho on mu (",
    nrow(x$data),
    " studies, ",
    "natural spline df = ",
    ncol(x$basis),
    ")\n",
    sep = ""
  )
  cat("  tau (between-study SD of rho): ", signif(x$tau, 3), "\n", sep = "")
  invisible(x)
}
