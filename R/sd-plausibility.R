# Plausibility bounds for the sample SD of bounded integer scores.
#
# Hard bounds (sd_min_quasi_integer(), sd_max_structure_s()) say which SDs are
# arithmetically possible. Plausibility bounds say which SDs a realistic
# data-generating process would produce: here, scores X in {0, ..., m},
# m = U - L, drawn iid from BetaBinomial(m, mu, rho), with mu = alpha /
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
#' @param m Whole number of scale steps, `U - L`.
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
    stop("The exact method needs about ", format(cells, big.mark = ","),
         " cells for N = ", N, " on a ", m, "-step scale, more than ",
         "max_cells = ", format(max_cells, big.mark = ","), ". Use an ",
         "approximate method for this N.", call. = FALSE)
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
#' The exact distribution of the unbiased sample variance \eqn{S^2} of \eqn{N}
#' iid BetaBinomial\eqn{(m, \mu, \rho)} scores, by dynamic programming over
#' the joint state \eqn{(\sum X, \sum X^2)}. Feasible because \eqn{X} is a
#' bounded integer; cost grows roughly as \eqn{N^2 m^3}, so this is for small
#' \eqn{N} (about \eqn{N \le 25} on a 0-63 scale).
#'
#' @inheritParams bb_cumulants
#' @param N Sample size(s), each \eqn{\ge 2}. A vector is computed in one
#'   pass up to `max(N)`.
#' @param sums Optional integer vector of attainable sums \eqn{\sum X}
#'   (on the \eqn{0..m} scale) to condition on, i.e. conditioning on the
#'   sample mean. Applied to every `N`.
#' @return A named list (by `N`) of data frames with columns `s2` and
#'   `prob`, sorted by `s2`. The attribute `lost` is the probability mass
#'   discarded by pruning (values and sums below \eqn{10^{-17}}).
#' @examples
#' d <- s2_exact(m = 4, mu = 0.4, rho = 0.1, N = 5)[["5"]]
#' sum(d$prob)
#' sum(d$s2 * d$prob)   # equals the population variance
#' @export
s2_exact <- function(m, mu, rho, N, sums = NULL) {
  j <- .sq_joint(.bb_pmf(m, mu, rho), max(N), keep = N)
  out <- lapply(j, .s2_from_joint, sums = sums)
  attr(out, "lost") <- attr(j, "lost")
  out
}

# ---- hard bounds for a reported mean -----------------------------------------

# Internal: sums t (on the 0..m scale) whose mean L + t / N rounds to M at
# `digits` decimals (half-width tolerance), or the exact sum when digits is
# NULL and M is GRIM-consistent.
.consistent_sums <- function(L, M, N, m, digits = NULL) {
  if (is.null(digits)) {
    t <- round(N * (M - L))
    if (abs(t - N * (M - L)) > 1e-7) return(integer(0))
    return(as.integer(t))
  }
  h <- 0.5 * 10^-digits + 1e-9
  t <- ceiling(N * (M - L - h)):floor(N * (M - L + h))
  t[t >= 0 & t <= N * m]
}

# Internal: hard min / max SD for integer data on [L, U] with reported mean M.
.hard_bounds <- function(L, U, M, N, m, digits) {
  t <- .consistent_sums(L, M, N, m, digits)
  if (length(t) == 0) {
    # no integer sample has this mean: fall back to the quasi-integer floor,
    # which is a valid floor at every real mean
    return(c(sd_min_quasi_integer(M - L, N), sd_max_structure_s(M, N, L, U)))
  }
  mt <- L + t / N
  c(min(sd_min_integer(mt - L, N)), max(sd_max_structure_s(mt, N, L, U)))
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

# Internal: lower and upper quantiles of S for one row.
.sd_plaus_row <- function(L, U, M, N, p, rho, tau, method, conditional,
                          digits, K) {
  m <- U - L
  mu <- (M - L) / m
  rhos <- .rho_nodes(rho, tau, K)
  if (method == "exact") {
    sums <- if (conditional) {
      t <- .consistent_sums(L, M, N, m, digits)
      if (length(t) == 0) {
        stop("Conditioning on the mean with the exact method needs a ",
             "GRIM-consistent M, or `digits` to define the sums that round ",
             "to it.", call. = FALSE)
      }
      t
    } else NULL
    dists <- lapply(rhos, function(r) s2_exact(m, mu, r, N, sums)[[1]])
    d <- do.call(rbind, dists)
    d$prob <- d$prob / length(rhos)
    agg <- rowsum(d$prob, round(d$s2 * N * (N - 1)), reorder = TRUE)
    s <- sqrt(as.numeric(rownames(agg)) / (N * (N - 1)))
    cdf <- cumsum(agg[, 1])
    q <- function(pp) s[which(cdf >= pp - 1e-12)[1]]
    return(c(q(p), q(1 - p)))
  }
  scs <- lapply(rhos, function(r) {
    .s2_cumulants(.cumulants_from_pmf(.bb_pmf(m, mu, r)), N, conditional)
  })
  cdf <- function(s) mean(vapply(scs, function(sc) .cdf_s_approx(s, sc, method),
                                 numeric(1)))
  hi <- 2 * sqrt(max(vapply(scs, `[[`, numeric(1), "k1"))) + 1
  q <- function(pp) {
    if (cdf(0) >= pp) return(0)
    stats::uniroot(function(s) cdf(s) - pp, c(0, hi), extendInt = "upX",
                   tol = 1e-10)$root
  }
  c(q(p), q(1 - p))
}

#' Plausibility bounds for a reported SD on a bounded integer scale
#'
#' Quantiles of the sampling distribution of the sample SD under a
#' beta-binomial data-generating process: the range of SDs that samples of
#' size `N` would produce at mean `M` if scores were
#' BetaBinomial\eqn{(U - L, \mu, \rho)} with \eqn{\mu = (M - L)/(U - L)}.
#' Returned next to the hard (arithmetic) bounds for integer data.
#'
#' @section Model:
#' \eqn{X = \text{score} - L \in \{0, \dots, m\}}, \eqn{m = U - L}, with
#' \eqn{E[X] = m\mu} and \eqn{\sigma^2 = m\mu(1 - \mu)[1 + (m - 1)\rho]}.
#' \eqn{\rho = 0} is the binomial, the least dispersed member of the family.
#' All moments are computed from the pmf.
#'
#' @section Methods:
#' For iid data with finite fourth moment, \eqn{E[S^2] = \sigma^2} and
#' \eqn{\mathrm{Var}(S^2) = \sigma^4 [2/(N-1) + \gamma_2/N]}.
#' * `"chisq2"` (default): \eqn{S^2 \approx \sigma^2 \chi^2_\nu / \nu} with
#'   \eqn{\nu = 2/[2/(N-1) + \gamma_2/N]} (Box / Satterthwaite).
#' * `"chisq3"`: \eqn{S^2 \approx a + b\chi^2_\nu}, matching the first three
#'   cumulants of \eqn{S^2}.
#' * `"cf"`: Cornish-Fisher expansion of \eqn{S^2} using its third cumulant
#'   \eqn{\kappa_3(S^2) = \kappa_6/N^2 + 12\kappa_2\kappa_4/(N(N-1)) +
#'   4(N-2)\kappa_3^2/(N(N-1)^2) + 8\kappa_2^3/(N-1)^2}.
#' * `"exact"`: the exact distribution by dynamic programming (see
#'   [s2_exact()]); small `N` only.
#'
#' @section Conditioning on the mean:
#' A reported SD comes with a reported mean. When the distribution is
#' skewed, \eqn{S^2} and \eqn{\bar X} are correlated
#' (\eqn{\mathrm{Cov}(\bar X, S^2) = \kappa_3/N}), so samples whose mean is
#' near \eqn{M} have a narrower SD distribution than samples in general. With
#' `condition_on_mean = TRUE`, `"chisq2"` uses the conditional variance
#' \eqn{\sigma^4[2/(N-1) + (\gamma_2 - \gamma_1^2)/N]}, and `"exact"`
#' conditions on the sums whose mean rounds to `M` (at `digits` decimals).
#' `"chisq3"` and `"cf"` do not support conditioning.
#'
#' @section Heterogeneity:
#' With `tau > 0`, \eqn{\rho} is drawn from \eqn{N(\rho, \tau^2)} truncated to
#' \eqn{[0, 1)}, represented by `K` equal-weight quantile nodes. The bound is
#' the quantile of the mixture: the per-node CDFs are averaged, then
#' inverted. Averaging per-node quantiles would be wrong. When `rho` and `tau`
#' come from [rho_corpus_fit()], consider adding the prediction's squared
#' standard error to `tau^2`.
#'
#' @param L,U Scale limits (scalars). `U - L` must be a whole number: the
#'   model is for integer scores. Rescale mean-scored scales to sum scores.
#' @param M Reported mean(s).
#' @param N Sample size(s), each \eqn{\ge 2}.
#' @param p Tail probability in **each** tail; several values give several
#'   rows. Default 0.01.
#' @param rho Intra-class correlation of the data-generating process, e.g.
#'   a corpus prediction from [rho_corpus_fit()]. Values below 0 are
#'   under-dispersed relative to the binomial and are set to 0 (with a note);
#'   values of 1 or more give the endpoint-only limit.
#' @param tau Between-study SD of `rho`. Default 0.
#' @param sd Reported SD(s), used to set `rho` by the method of moments when
#'   `rho` is `NULL`. This is circular for checking that same SD and is
#'   meant for illustration.
#' @param method `"chisq2"`, `"chisq3"`, `"cf"` or `"exact"`.
#' @param condition_on_mean Logical; see the section above. Default `FALSE`.
#' @param digits Decimal places of the reported mean, if known. Used for the
#'   hard bounds (minimum and maximum over every sum that rounds to `M`) and
#'   for exact conditioning. `NULL` treats `M` as exact.
#' @param K Number of quadrature nodes for the \eqn{\rho} mixture.
#'
#' @return A tibble with one row per (`M`, `N`, `rho`, `tau`, `sd`) element
#'   (recycled together) and `p` value (crossed): the inputs, `mu`, `rho`,
#'   `sigma` (population SD at `rho`), `hard_min` and `hard_max` (integer
#'   data), `lower` and `upper` (the plausibility bounds, clipped to the hard
#'   bounds), `clipped` (which side, if any, was clipped) and `note`.
#' @seealso [sd_min_integer()], [sd_max_structure_s()] for the hard bounds;
#'   [rho_from_summary()], [rho_corpus_fit()] for estimating `rho`.
#' @examples
#' # BDI-II (0-63), M = 25, SD = 8 sets rho by the method of moments
#' sd_plausibility(0, 63, M = 25, N = c(20, 50, 100, 500), p = c(0.01, 0.001),
#'                 sd = 8)
#'
#' # a corpus prediction with heterogeneity
#' sd_plausibility(0, 63, M = 25, N = 50, rho = 0.05, tau = 0.02)
#' @export
sd_plausibility <- function(L, U, M, N, p = 0.01, rho = NULL, tau = 0,
                            sd = NULL,
                            method = c("chisq2", "chisq3", "cf", "exact"),
                            condition_on_mean = FALSE, digits = NULL,
                            K = 41) {
  method <- match.arg(method)
  if (length(L) != 1 || length(U) != 1 || anyNA(c(L, U)) || !(U > L)) {
    stop("`L` and `U` must be single numbers with U > L.", call. = FALSE)
  }
  m <- U - L
  if (abs(m - round(m)) > 1e-9) {
    stop("`U - L` must be a whole number: the beta-binomial model is for ",
         "integer scores. For a mean-scored scale, convert to sum scores ",
         "(multiply L, U, M and SD by the number of items).", call. = FALSE)
  }
  m <- round(m)
  if (m < 1) stop("`U - L` must be at least 1.", call. = FALSE)
  if (is.null(rho) && is.null(sd)) {
    stop("Supply `rho` (e.g. a corpus prediction) or `sd` (method of ",
         "moments; circular for checking that SD).", call. = FALSE)
  }
  if (condition_on_mean && method %in% c("chisq3", "cf")) {
    stop("`condition_on_mean` is available for methods \"chisq2\" and ",
         "\"exact\".", call. = FALSE)
  }
  len <- max(length(M), length(N), length(rho), length(tau), length(sd))
  M <- rep_len(M, len); N <- rep_len(N, len); tau <- rep_len(tau, len)
  rho_in <- if (is.null(rho)) rep(NA_real_, len) else rep_len(rho, len)
  sd_in <- if (is.null(sd)) rep(NA_real_, len) else rep_len(sd, len)

  rows <- vector("list", len * length(p))
  r <- 0
  for (i in seq_len(len)) {
    for (pp in p) {
      r <- r + 1
      mu <- (M[i] - L) / m
      rho_i <- if (!is.na(rho_in[i])) rho_in[i] else
        (sd_in[i]^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
      res <- list(L = L, U = U, M = M[i], N = N[i], p = pp, method = method,
                  conditional = condition_on_mean, mu = mu, rho = rho_i,
                  tau = tau[i], sigma = NA_real_, hard_min = NA_real_,
                  hard_max = NA_real_, lower = NA_real_, upper = NA_real_,
                  clipped = NA_character_, note = NA_character_)
      note <- character(0)
      if (anyNA(c(M[i], N[i], rho_i, tau[i], pp))) {
        res$note <- "missing input"
      } else if (N[i] < 2 || N[i] != round(N[i])) {
        res$note <- "N must be a whole number >= 2"
      } else if (!(pp > 0 && pp < 0.5)) {
        res$note <- "p must be in (0, 0.5)"
      } else if (M[i] < L || M[i] > U) {
        res$note <- "M outside [L, U]"
      } else if (M[i] == L || M[i] == U) {
        res[c("sigma", "hard_min", "hard_max", "lower", "upper")] <- 0
        res$note <- "M at a scale limit: every score equals it"
      } else {
        if (rho_i < 0) {
          note <- c(note, sprintf(
            "rho = %.4g < 0 is under-dispersed; bounds use rho = 0", rho_i))
          rho_i <- 0
        }
        if (rho_i >= 1) {
          note <- c(note, "rho >= 1: all mass at the scale endpoints")
          rho_i <- 1
        }
        hb <- .hard_bounds(L, U, M[i], N[i], m, digits)
        q <- .sd_plaus_row(L, U, M[i], N[i], pp, rho_i, tau[i], method,
                           condition_on_mean, digits, K)
        clipped <- c(if (q[1] < hb[1] - 1e-12) "lower",
                     if (q[2] > hb[2] + 1e-12) "upper")
        if (length(clipped)) {
          note <- c(note, paste("clipped to the hard",
                                paste(clipped, collapse = " and "),
                                "bound"))
        }
        res$rho <- rho_i
        res$sigma <- sqrt(.cumulants_from_pmf(.bb_pmf(m, mu, rho_i))$k2)
        res$hard_min <- hb[1]; res$hard_max <- hb[2]
        res$lower <- max(q[1], hb[1]); res$upper <- min(q[2], hb[2])
        res$clipped <- if (length(clipped)) paste(clipped, collapse = "+") else "none"
        res$note <- if (length(note)) paste(note, collapse = "; ") else NA_character_
      }
      rows[[r]] <- tibble::as_tibble(res)
    }
  }
  do.call(rbind, rows)
}

# ---- estimating rho from summary statistics ----------------------------------

#' Method-of-moments rho from summary statistics
#'
#' \eqn{\hat\mu = (M - L)/m} and \eqn{\hat\rho = [SD^2/(m\hat\mu(1 -
#' \hat\mu)) - 1]/(m - 1)}. \eqn{\hat\rho} is linear in \eqn{S^2}, so it is
#' unbiased for \eqn{\rho} (given \eqn{\mu}), with sampling variance
#' \eqn{\mathrm{Var}(S^2)/[m(m - 1)\mu(1 - \mu)]^2}. Negative
#' \eqn{\hat\rho} are kept, not truncated.
#'
#' @param L,U Scale limits; `U - L` a whole number. Scalars or vectors.
#' @param M,SD,N Reported mean, SD and sample size (vectors).
#' @param rho_var Optional \eqn{\rho} at which to evaluate the kurtosis in
#'   the variance (default: \eqn{\max(0, \hat\rho)}). [rho_corpus_fit()]
#'   passes the shrunken corpus prediction here.
#' @return A tibble with `mu_hat`, `rho_hat` and `var_rho`.
#' @note The variance treats \eqn{\mu} as known. Its sampling error adds a
#'   term of order \eqn{1/N} that is small next to \eqn{\mathrm{Var}(S^2)}
#'   except near the scale limits.
#' @examples
#' rho_from_summary(0, 63, M = c(25, 12), SD = c(8, 7), N = c(50, 120))
#' @export
rho_from_summary <- function(L, U, M, SD, N, rho_var = NULL) {
  len <- max(length(L), length(U), length(M), length(SD), length(N))
  L <- rep_len(L, len); U <- rep_len(U, len); M <- rep_len(M, len)
  SD <- rep_len(SD, len); N <- rep_len(N, len)
  m <- U - L
  if (any(abs(m - round(m)) > 1e-9, na.rm = TRUE)) {
    stop("`U - L` must be a whole number.", call. = FALSE)
  }
  m <- round(m)
  mu <- (M - L) / m
  rho_hat <- (SD^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
  rv <- if (is.null(rho_var)) pmax(0, rho_hat) else rep_len(rho_var, len)
  var_rho <- vapply(seq_len(len), function(i) {
    if (anyNA(c(mu[i], rv[i], N[i])) || mu[i] <= 0 || mu[i] >= 1 || N[i] < 2) {
      return(NA_real_)
    }
    cu <- .cumulants_from_pmf(.bb_pmf(m[i], mu[i], min(max(rv[i], 0), 1)))
    .s2_cumulants(cu, N[i])$k2 / (m[i] * (m[i] - 1) * mu[i] * (1 - mu[i]))^2
  }, numeric(1))
  tibble::tibble(mu_hat = mu, rho_hat = rho_hat, var_rho = var_rho)
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
#' @param corpus A data frame with columns `L`, `U`, `M`, `SD` and `N`.
#' @param df Degrees of freedom of the natural spline. Default 3.
#' @param iterate Logical; re-evaluate the variances at the prediction once.
#' @param ... Further arguments passed to [metafor::rma()].
#' @return An object of class `rho_corpus`: the `fit` ([metafor::rma()]),
#'   the spline `basis`, the per-study `data` (with `mu_hat`, `rho_hat`,
#'   `var_rho`) and `tau` (\eqn{\hat\tau}). Use [predict()] on it.
#' @seealso [predict.rho_corpus()], [sd_plausibility()]
#' @export
rho_corpus_fit <- function(corpus, df = 3, iterate = TRUE, ...) {
  if (!requireNamespace("metafor", quietly = TRUE)) {
    stop("rho_corpus_fit() needs the metafor package.", call. = FALSE)
  }
  need <- c("L", "U", "M", "SD", "N")
  if (!is.data.frame(corpus) || !all(need %in% names(corpus))) {
    stop("`corpus` must have columns L, U, M, SD and N.", call. = FALSE)
  }
  d <- cbind(corpus, rho_from_summary(corpus$L, corpus$U, corpus$M,
                                      corpus$SD, corpus$N))
  d <- d[stats::complete.cases(d[c("mu_hat", "rho_hat", "var_rho")]), ]
  basis <- splines::ns(d$mu_hat, df = df)
  fit_once <- function(v) {
    metafor::rma(yi = d$rho_hat, vi = v, mods = ~ basis, ...)
  }
  fit <- fit_once(d$var_rho)
  if (iterate) {
    pred <- as.numeric(stats::predict(fit)$pred)
    d$var_rho <- rho_from_summary(d$L, d$U, d$M, d$SD, d$N,
                                  rho_var = pred)$var_rho
    fit <- fit_once(d$var_rho)
  }
  structure(list(fit = fit, basis = basis, data = tibble::as_tibble(d),
                 tau = sqrt(fit$tau2)),
            class = "rho_corpus")
}

#' Predict rho and tau from a corpus fit
#'
#' @param object A `rho_corpus` object from [rho_corpus_fit()].
#' @param mu Mean proportions \eqn{(M - L)/(U - L)} at which to predict.
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
