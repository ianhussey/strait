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
.bb_pmf <- function(m, mu, rho) {
  x <- 0:m
  if (rho <= 1e-10) return(stats::dbinom(x, m, mu))
  if (rho >= 1) return(c(1 - mu, rep(0, m - 1), mu))
  phi <- 1 / rho - 1
  a <- mu * phi
  b <- (1 - mu) * phi
  p <- exp(lchoose(m, x) + lbeta(x + a, m - x + b) - lbeta(a, b))
  p / sum(p)
}

# Internal: population cumulants of X from the pmf.
.cumulants_from_pmf <- function(pr) {
  x <- seq_along(pr) - 1
  mean <- sum(pr * x)
  cm <- vapply(2:6, function(r) sum(pr * (x - mean)^r), numeric(1))
  m2 <- cm[1]; m3 <- cm[2]; m4 <- cm[3]; m5 <- cm[4]; m6 <- cm[5]
  k4 <- m4 - 3 * m2^2
  list(
    mean = mean, k2 = m2, k3 = m3, k4 = k4,
    k5 = m5 - 10 * m3 * m2,
    k6 = m6 - 15 * m4 * m2 - 10 * m3^2 + 30 * m2^3,
    g1 = m3 / m2^1.5, g2 = k4 / m2^2
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
  tibble::as_tibble(.cumulants_from_pmf(.bb_pmf(m, mu, rho)))
}

# Internal: cumulants of the unbiased sample variance S^2 of N iid draws.
# conditional = TRUE gives the leading-order moments of S^2 given that the
# sample mean equals the population mean, removing the part of Var(S^2) that
# Cov(xbar, S^2) = kappa_3 / N explains.
.s2_cumulants <- function(cu, N, conditional = FALSE) {
  k2 <- cu$k2; k3 <- cu$k3; k4 <- cu$k4; k6 <- cu$k6
  v <- k2^2 * (2 / (N - 1) + cu$g2 / N)
  if (conditional) v <- v - k3^2 / (N * k2)
  t3 <- k6 / N^2 + 12 * k2 * k4 / (N * (N - 1)) +
    4 * (N - 2) * k3^2 / (N * (N - 1)^2) + 8 * k2^3 / (N - 1)^2
  list(k1 = k2, k2 = v, k3 = t3)
}

# ---- approximate CDFs of S ---------------------------------------------------

# Internal: CDF of S at s for one rho, under an approximation of S^2.
.cdf_s_approx <- function(s, sc, method) {
  s2 <- s^2
  switch(method,
    chisq2 = {
      nu <- 2 * sc$k1^2 / sc$k2
      stats::pchisq(s2 * nu / sc$k1, nu)
    },
    chisq3 = {
      if (sc$k3 <= 0) return(.cdf_s_approx(s, sc, "chisq2"))
      b <- sc$k3 / (4 * sc$k2)
      nu <- 8 * sc$k2^3 / sc$k3^2
      a <- sc$k1 - b * nu
      stats::pchisq(pmax(0, (s2 - a) / b), nu)
    },
    cf = {
      # invert the Cornish-Fisher quantile k1 + sqrt(k2) (z + g/6 (z^2 - 1))
      # for z; the + root is the monotone branch
      w <- (s2 - sc$k1) / sqrt(sc$k2)
      h <- sc$k3 / sc$k2^1.5 / 6
      if (abs(h) < 1e-8) return(stats::pnorm(w))
      disc <- 1 + 4 * h * (h + w)
      z <- ifelse(disc < 0, if (h > 0) -Inf else Inf,
                  (-1 + sqrt(pmax(disc, 0))) / (2 * h))
      stats::pnorm(z)
    },
    stop("unknown method: ", method)
  )
}

# ---- exact distribution by dynamic programming -------------------------------

# Internal: minimum and maximum of sum x^2 over n integers in [0, m] summing
# to s.
.q_min <- function(n, s) {
  if (n == 0) return(0)
  a <- s %/% n; r <- s %% n
  (n - r) * a^2 + r * (a + 1)^2
}
.q_max <- function(n, s, m) {
  k <- s %/% m; rem <- s %% m
  k * m^2 + rem^2
}

# Internal: exact joint distribution of (sum X, sum X^2) for N iid draws from
# the pmf `pr` on 0:m, by dynamic programming over observations. State: for
# each sum s, a vector over q in [q_min(n, s), q_max(n, s)]. Values with
# p(x) < prune and sums whose total mass is below prune are dropped; the
# discarded mass is returned as `lost`. `keep` lists the sample sizes whose
# joint distribution is returned (all <= N), so one run serves several N.
.sq_joint <- function(pr, N, keep = N, prune = 1e-17, max_cells = 4e7) {
  m <- length(pr) - 1
  cells <- sum(vapply(0:(N * m), function(s)
    .q_max(N, s, m) - .q_min(N, s) + 1, numeric(1)))
  if (cells > max_cells) {
    stop(
      "The exact method needs about ",
      format(cells, big.mark = ","),
      " cells for n = ",
      N,
      " on a ",
      m,
      "-step scale, more than ",
      "max_cells = ",
      format(max_cells, big.mark = ","),
      ". Use an ",
      "approximate method for this n.",
      call. = FALSE
    )
  }
  xs <- which(pr >= prune) - 1
  state <- list(1)                     # n = 0: s = 0, q = 0
  lost <- 0
  out <- list()
  for (n in 0:(N - 1)) {
    smax_new <- (n + 1) * m
    new <- vector("list", smax_new + 1)
    for (s2 in 0:smax_new) {
      base <- .q_min(n + 1, s2)
      v <- NULL
      for (x in xs) {
        s <- s2 - x
        if (s < 0 || s > n * m) next
        o <- state[[s + 1]]
        if (is.null(o)) next
        if (is.null(v)) v <- numeric(.q_max(n + 1, s2, m) - base + 1)
        idx <- .q_min(n, s) + x^2 - base + seq_along(o)
        v[idx] <- v[idx] + pr[x + 1] * o
      }
      if (!is.null(v)) {
        tot <- sum(v)
        if (tot < prune) lost <- lost + tot else new[[s2 + 1]] <- v
      }
    }
    state <- new
    if ((n + 1) %in% keep) {
      out[[as.character(n + 1)]] <- list(state = state, N = n + 1, m = m)
    }
  }
  structure(out, lost = lost + (1 - sum(pr[xs + 1])) * N)
}

# Internal: distribution of S^2 from a joint state, optionally restricted to
# the sums in `sums` (conditioning on the sample mean). Returns a data frame
# (s2, prob) sorted by s2, with prob normalised within the condition.
.s2_from_joint <- function(j, sums = NULL) {
  N <- j$N; st <- j$state
  ss <- which(!vapply(st, is.null, logical(1))) - 1
  if (!is.null(sums)) ss <- intersect(ss, sums)
  if (length(ss) == 0) {
    stop("No attainable sample has a mean in the conditioning set.",
         call. = FALSE)
  }
  num <- list(); pr <- list()
  for (s in ss) {
    v <- st[[s + 1]]
    nz <- which(v > 0)
    q <- .q_min(N, s) + nz - 1
    num[[length(num) + 1]] <- N * q - s^2   # N (N - 1) S^2, an integer
    pr[[length(pr) + 1]] <- v[nz]
  }
  agg <- rowsum(unlist(pr), unlist(num), reorder = TRUE)
  k <- as.numeric(rownames(agg))
  data.frame(s2 = k / (N * (N - 1)), prob = agg[, 1] / sum(agg[, 1]))
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
  .sdp_check_n(n, allow_na = FALSE)
  j <- .sq_joint(.bb_pmf(m, mu, rho), max(n), keep = n)
  out <- lapply(j, function(x) tibble::as_tibble(.s2_from_joint(x, sums)))
  attr(out, "lost") <- attr(j, "lost")
  out
}

# Internal: stop unless n holds whole numbers >= 2 (NA allowed if allow_na).
.sdp_check_n <- function(n, allow_na = TRUE) {
  v <- if (allow_na) n[!is.na(n)] else n
  if (
    !is.numeric(n) ||
      anyNA(v) ||
      any(!is.finite(v)) ||
      any(v < 2) ||
      any(v != round(v))
  ) {
    stop("`n` must contain finite whole numbers >= 2.", call. = FALSE)
  }
}
# ---- hard bounds for a reported mean -----------------------------------------

# Internal: sums t (on the 0..m scale) whose mean l + t / n rounds to `mean`
# at `digits` decimals (half-width tolerance), or the exact sum when digits
# is NULL and the mean is GRIM-consistent. `digits` is a single value.
.consistent_sums <- function(l, mean, n, m, digits = NULL) {
  if (is.null(digits)) {
    t <- round(n * (mean - l))
    if (abs(t - n * (mean - l)) > 1e-7) {
      return(integer(0))
    }
    return(as.integer(t))
  }
  h <- 0.5 * 10^-digits + 1e-9
  lo <- ceiling(n * (mean - l - h))
  hi <- floor(n * (mean - l + h))
  if (hi < lo) {
    return(integer(0))
  }
  t <- lo:hi
  t[t >= 0 & t <= n * m]
}

# Internal: hard min / max SD for integer data on [l, u] with reported mean.
.hard_bounds <- function(l, u, mean, n, m, digits) {
  t <- .consistent_sums(l, mean, n, m, digits)
  if (length(t) == 0) {
    # no integer sample has this mean: fall back to the quasi-integer floor,
    # which is a valid floor at every real mean
    return(c(
      sd_min_quasi_integer(mean - l, n),
      sd_max_structure_s(mean, n, l, u)
    ))
  }
  mt <- l + t / n
  c(min(sd_min_integer(mt - l, n)), max(sd_max_structure_s(mt, n, l, u)))
}

# ---- the main function -------------------------------------------------------

# Internal: rho nodes of N(rho, tau^2) truncated to [0, 1), equal weights.
.rho_nodes <- function(rho, tau, K) {
  if (tau <= 0) return(rho)
  pa <- stats::pnorm(0, rho, tau)
  pb <- stats::pnorm(1, rho, tau)
  u <- (seq_len(K) - 0.5) / K
  pmin(stats::qnorm(pa + (pb - pa) * u, rho, tau), 1 - 1e-9)
}

# Internal: lower and upper quantiles of S for one row, at every tail
# probability in `p`. The distribution (exact) or the cumulants (approximate)
# are computed once per rho node and shared by all `p`.
.sd_plaus_row <- function(
  l,
  u,
  mean,
  n,
  p,
  rho,
  tau,
  method,
  conditional,
                          digits, K) {
  m <- u - l
  mu <- (mean - l) / m
  rhos <- .rho_nodes(rho, tau, K)
  if (method == "exact") {
    sums <- if (conditional) {
      t <- .consistent_sums(l, mean, n, m, digits)
      if (length(t) == 0) {
        stop(
          "Conditioning on the mean with the exact method needs a ",
          "GRIM-consistent mean, or `digits_mean` to define the sums ",
          "that round to it.",
          call. = FALSE
        )
      }
      t
    }
    d <- do.call(
      rbind,
      lapply(rhos, function(r) {
        .s2_from_joint(.sq_joint(.bb_pmf(m, mu, r), n)[[1]], sums)
      })
    )
    agg <- rowsum(
      d$prob / length(rhos),
      round(d$s2 * n * (n - 1)),
      reorder = TRUE
    )
    s <- sqrt(as.numeric(rownames(agg)) / (n * (n - 1)))
    cdf <- cumsum(agg[, 1])
    q <- function(pp) s[which(cdf >= pp - 1e-12)[1]]
  } else {
    scs <- lapply(rhos, function(r) {
      .s2_cumulants(.cumulants_from_pmf(.bb_pmf(m, mu, r)), n, conditional)
  })
  cdf <- function(s) mean(vapply(scs, function(sc) .cdf_s_approx(s, sc, method),
                                 numeric(1)))
  hi <- 2 * sqrt(max(vapply(scs, `[[`, numeric(1), "k1"))) + 1
  q <- function(pp) {
    if (cdf(0) >= pp) return(0)
    stats::uniroot(function(s) cdf(s) - pp, c(0, hi), extendInt = "upX",
        tol = 1e-10
      )$root
    }
  }
  list(lower = vapply(p, q, numeric(1)), upper = vapply(1 - p, q, numeric(1)))
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
  .ga_check_dots(..., .renamed = .sd_renamed)
  method <- match.arg(method)
  m <- .sdp_check_scale(l, u)
  if (is.null(rho) == is.null(sd)) {
    stop(
      "Supply exactly one of `rho` (e.g. a corpus prediction) and `sd` ",
      "(method of moments; circular for checking that SD).",
      call. = FALSE
    )
  }
  if (
    !is.logical(condition_on_mean) ||
      length(condition_on_mean) != 1 ||
      is.na(condition_on_mean)
  ) {
    stop("`condition_on_mean` must be TRUE or FALSE.", call. = FALSE)
  }
  if (condition_on_mean && method %in% c("chisq3", "cf")) {
    stop(
      "Methods \"chisq3\" and \"cf\" cannot condition on the mean: set ",
      "`condition_on_mean = FALSE`, or use \"chisq2\" or \"exact\".",
      call. = FALSE
    )
  }
  .sdp_check_n(n)
  if (!is.numeric(p) || length(p) == 0 || anyNA(p) || any(p <= 0 | p >= 0.5)) {
    stop("`p` must contain tail probabilities in (0, 0.5).", call. = FALSE)
  }
  if (!is.numeric(tau) || any(tau < 0, na.rm = TRUE)) {
    stop("`tau` must be non-negative.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1 || is.na(K) || K < 1 || K != round(K)) {
    stop("`K` must be a whole number >= 1.", call. = FALSE)
  }
  if (!is.null(digits_mean)) {
    dg <- digits_mean[!is.na(digits_mean)]
    if (!is.numeric(digits_mean) || any(dg < 0) || any(dg != round(dg))) {
      stop(
        "`digits_mean` must contain whole numbers >= 0 (or NA).",
        call. = FALSE
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
  args <- args[!vapply(args, is.null, logical(1))]
  lens <- lengths(args)
  len <- max(lens)
  if (any(!lens %in% c(1L, len))) {
    stop(
      "`mean`, `n`, `sd`, `rho`, `tau` and `digits_mean` must each have ",
      "length 1 or a common length.",
      call. = FALSE
    )
  }
  mean <- rep_len(mean, len)
  n <- rep_len(n, len)
  tau <- rep_len(tau, len)
  rho_in <- if (is.null(rho)) rep(NA_real_, len) else rep_len(rho, len)
  sd_in <- if (is.null(sd)) rep(NA_real_, len) else rep_len(sd, len)
  digits_in <- if (is.null(digits_mean)) {
    rep(NA_real_, len)
  } else {
    rep_len(digits_mean, len)
  }

  rows <- vector("list", len)
  for (i in seq_len(len)) {
    mu <- (mean[i] - l) / m
    res <- list(
      l = l,
      u = u,
      mean = mean[i],
      n = n[i],
      p = p,
      method = method,
      conditional = condition_on_mean,
      mu = mu,
      rho = rho_in[i],
                  tau = tau[i], sigma = NA_real_, hard_min = NA_real_,
                  hard_max = NA_real_, lower = NA_real_, upper = NA_real_,
                  clipped = NA_character_, note = NA_character_)
      note <- character(0)
    if (
      anyNA(c(mean[i], n[i], tau[i])) ||
        (is.na(rho_in[i]) && is.na(sd_in[i]))
    ) {
        res$note <- "missing input"
    } else if (mean[i] < l || mean[i] > u) {
      res$note <- "mean outside [l, u]"
    } else if (mean[i] == l || mean[i] == u) {
        res[c("sigma", "hard_min", "hard_max", "lower", "upper")] <- 0
      res$clipped <- "none"
      res$note <- "mean at a scale limit: every score equals it"
    } else {
      if (m == 1) {
        # Bernoulli(mu) whatever rho: the pmf ignores it
        note <- "u - l = 1: scores are binary, so rho is not identified"
        rho_i <- 0
        tau_i <- 0
        res$rho <- NA_real_
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
        res$rho <- rho_i
      }
      digits_i <- if (is.na(digits_in[i])) NULL else digits_in[i]
      hb <- .hard_bounds(l, u, mean[i], n[i], m, digits_i)
      q <- .sd_plaus_row(
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
      lo_clip <- q$lower < hb[1] - 1e-12
      hi_clip <- q$upper > hb[2] + 1e-12
      res$sigma <- sqrt(.cumulants_from_pmf(.bb_pmf(m, mu, rho_i))$k2)
        res$hard_min <- hb[1]; res$hard_max <- hb[2]
      res$lower <- pmax(q$lower, hb[1])
      res$upper <- pmin(q$upper, hb[2])
      res$clipped <- ifelse(
        lo_clip & hi_clip,
        "lower+upper",
        ifelse(lo_clip, "lower", ifelse(hi_clip, "upper", "none"))
      )
      clip_note <- ifelse(
        res$clipped == "none",
        NA,
        paste(
          "clipped to the hard",
          sub("+", " and ", res$clipped, fixed = TRUE),
          "bound"
        )
      )
      res$note <- vapply(
        clip_note,
        function(cn) {
          all_notes <- c(note, if (!is.na(cn)) cn)
          if (length(all_notes)) {
            paste(all_notes, collapse = "; ")
          } else {
            NA_character_
          }
        },
        character(1),
        USE.NAMES = FALSE
      )
    }
    rows[[i]] <- tibble::as_tibble(res)
  }
  do.call(rbind, rows)
}

# Internal: old argument and column names, for clear errors.
.sd_renamed <- c(
  L = "l",
  U = "u",
  M = "mean",
  N = "n",
  SD = "sd",
  digits = "digits_mean"
)

# Internal: validate scalar scale limits; return m = u - l.
.sdp_check_scale <- function(l, u) {
  if (
    !is.numeric(l) ||
      !is.numeric(u) ||
      length(l) != 1 ||
      length(u) != 1 ||
      anyNA(c(l, u)) ||
      !(u > l) ||
      !is.finite(u - l)
  ) {
    stop("`l` and `u` must be single finite numbers with u > l.", call. = FALSE)
  }
  m <- u - l
  if (abs(m - round(m)) > 1e-9) {
    stop(
      "`u - l` must be a whole number: the beta-binomial model is for ",
      "integer scores. For a mean-scored scale, convert to sum scores ",
      "(multiply l, u, mean and sd by the number of items).",
      call. = FALSE
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
  .ga_check_dots(..., .renamed = .sd_renamed)
  m <- .sdp_check_scale(l, u)
  if (m < 2) {
    stop(
      "`u - l` must be at least 2: on a two-point scale scores are ",
      "binary and rho is not identified.",
      call. = FALSE
    )
  }
  lens <- lengths(list(n, mean, sd))
  len <- max(lens)
  if (any(!lens %in% c(1L, len))) {
    stop(
      "`n`, `mean` and `sd` must each have length 1 or a common length.",
      call. = FALSE
    )
  }
  n <- rep_len(n, len)
  mean <- rep_len(mean, len)
  sd <- rep_len(sd, len)
  mu <- (mean - l) / m
  rho_hat <- (sd^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
  rv <- if (is.null(rho_var)) pmax(0, rho_hat) else rep_len(rho_var, len)
  var_rho <- vapply(seq_len(len), function(i) {
      if (
        anyNA(c(mu[i], rv[i], n[i])) || mu[i] <= 0 || mu[i] >= 1 || n[i] < 2
      ) {
      return(NA_real_)
    }
      cu <- .cumulants_from_pmf(.bb_pmf(m, mu[i], min(max(rv[i], 0), 1)))
      .s2_cumulants(cu, n[i])$k2 / (m * (m - 1) * mu[i] * (1 - mu[i]))^2
  }, numeric(1))
  tibble::tibble(mu_hat = mu, rho_hat = rho_hat, var_rho = var_rho)
}

# Internal: rho_from_summary() row by row, for corpora mixing scales.
.rho_corpus_rows <- function(d, rho_var = rep(NA_real_, nrow(d))) {
  do.call(
    rbind,
    lapply(seq_len(nrow(d)), function(i) {
      rho_from_summary(
        d$l[i],
        d$u[i],
        n = d$n[i],
        mean = d$mean[i],
        sd = d$sd[i],
        rho_var = if (!is.na(rho_var[i])) rho_var[i]
      )
    })
  )
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
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("rho_corpus_fit() needs the metafor package.", call. = FALSE)
  }
  need <- c("l", "u", "n", "mean", "sd")
  if (!is.data.frame(corpus) || !all(need %in% names(corpus))) {
    hint <- if (
      is.data.frame(corpus) &&
        all(names(.sd_renamed)[1:5] %in% names(corpus))
    ) {
      " Columns `L`, `U`, `M`, `SD` and `N` were renamed."
    } else {
      ""
    }
    stop(
      "`corpus` must have columns `l`, `u`, `n`, `mean` and `sd`.",
      hint,
      call. = FALSE
    )
  }
  d <- cbind(corpus, .rho_corpus_rows(corpus))
  d <- d[stats::complete.cases(d[c("mu_hat", "rho_hat", "var_rho")]), ]
  basis <- splines::ns(d$mu_hat, df = df)
  fit_once <- function(v) {
    metafor::rma(yi = d$rho_hat, vi = v, mods = ~ basis, ...)
  }
  fit <- fit_once(d$var_rho)
  if (iterate) {
    pred <- as.numeric(stats::predict(fit)$pred)
    d$var_rho <- .rho_corpus_rows(d, rho_var = pred)$var_rho
    fit <- fit_once(d$var_rho)
  }
  structure(list(fit = fit, basis = basis, data = tibble::as_tibble(d),
                 tau = sqrt(fit$tau2)),
            class = "rho_corpus")
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
  X <- stats::predict(object$basis, mu)
  pr <- stats::predict(object$fit, newmods = X)
  tibble::tibble(mu = mu, rho = as.numeric(pr$pred), se = as.numeric(pr$se),
                 tau = object$tau,
                 tau_total = sqrt(object$tau^2 + as.numeric(pr$se)^2))
}

#' @export
print.rho_corpus <- function(x, ...) {
  cat("Corpus meta-regression of rho on mu (", nrow(x$data), " studies, ",
      "natural spline df = ", ncol(x$basis), ")\n", sep = "")
  cat("  tau (between-study SD of rho): ", signif(x$tau, 3), "\n", sep = "")
  invisible(x)
}
