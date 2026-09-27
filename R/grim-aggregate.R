# Pooling many GRIM tests from one paper (or one group of values) into a
# small set of quantities, each tagged with the assumptions it rests on.
#
# The design principle is an assumption ladder. Step 0-2 are deterministic
# counts; step 3 is a partial-identification bound that needs only A0; steps
# 4-6 compare the observed failures with what uniform-random values would
# produce; steps 7-8 fit a mixture model; step 9 compares against a
# user-supplied authentic baseline. The user stops at whichever rung they are
# willing to defend. See paper/grim-aggregation/ for the derivations.
#
# Assumption labels:
#   A0  no authentic failures: correct values never fail GRIM
#   U   uniform counterfactual: a non-genuine value is equally likely to be
#       any value on the reporting grid (or distributed as `g`)
#   I   independence: values pass or fail independently
#   M   mixture: each value is non-genuine with probability lambda
#   B   authentic baseline: an authentic value fails with probability eps_i

# Grid-consistency vectors, keyed by "digits|n|items|rounding". Each entry is
# the logical vector c(j / 10^d; n, items, d) for j = 0, ..., 10^d - 1. GRIM
# consistency is periodic with period 1 in the mean, so one unit interval is
# exact. The cache persists for the session: pass rates are a property of the
# design, not of the reported value.
.grim_grid_cache <- new.env(parent = emptyenv())

# Internal: scrutiny::grim_map() consistency across both interfaces, for
# means given as character strings that all share `digits` decimal places.
#
# scrutiny <= 0.6.1 reads the decimal places off a character x. Later builds
# reject a character x outright and take a numeric x with a scalar digits_x
# argument (see R/scrutiny-compat.R). Returns a logical vector in row order.
.grim_map_consistency <- function(x, n, items, digits, rounding) {
  df <- tibble::tibble(x = x, n = n, items = items)
  if (.takes_digits(scrutiny::grim_map)) {
    df$x <- as.numeric(df$x)
    out <- scrutiny::grim_map(df, digits_x = digits, rounding = rounding)
  } else {
    out <- scrutiny::grim_map(df, rounding = rounding)
  }
  as.logical(out$consistency)
}

# Internal: observed GRIM consistency for rows with mixed decimal places.
.grim_observed <- function(x, n, items, digits, rounding) {
  out <- logical(length(x))
  for (d in unique(digits)) {
    i <- which(digits == d)
    out[i] <- .grim_map_consistency(x[i], n[i], items[i], d, rounding)
  }
  out
}

# Internal: the grid-consistency vector for every unique (digits, n, items),
# computed in one grim_map() call per decimal count and cached.
.grim_grid <- function(digits, n, items, rounding) {
  key <- paste(digits, n, items, rounding, sep = "|")
  todo <- unique(data.frame(key, digits, n, items)[
    !vapply(key, exists, logical(1), envir = .grim_grid_cache,
            inherits = FALSE), ])
  for (d in unique(todo$digits)) {
    t <- todo[todo$digits == d, ]
    grid <- formatC((0:(10^d - 1)) / 10^d, format = "f", digits = d)
    len <- length(grid)
    cons <- .grim_map_consistency(
      x = rep(grid, nrow(t)), n = rep(t$n, each = len),
      items = rep(t$items, each = len), digits = d, rounding = rounding
    )
    for (r in seq_len(nrow(t))) {
      assign(t$key[r], cons[(r - 1) * len + seq_len(len)],
             envir = .grim_grid_cache)
    }
  }
  lapply(key, get, envir = .grim_grid_cache, inherits = FALSE)
}

# Internal: grid weights for decimal count d under counterfactual g.
# NULL is uniform; a function maps final digits 0:9 to weights; a numeric
# vector gives one weight per grid value j / 10^d, j = 0, ..., 10^d - 1.
.grid_weights <- function(g, d) {
  len <- 10^d
  if (is.null(g)) {
    w <- rep(1, len)
  } else if (is.function(g)) {
    wd <- g(0:9)
    if (!is.numeric(wd) || length(wd) != 10L) {
      stop("`g` as a function must return 10 weights, one per final digit ",
           "0-9, when called on 0:9.", call. = FALSE)
    }
    w <- wd[(0:(len - 1)) %% 10 + 1]
  } else if (is.numeric(g)) {
    if (length(g) != len) {
      stop("`g` as a numeric vector must have one weight per grid value, ",
           "i.e. length 10^d = ", len, " for values reported to ", d,
           " decimals (got length ", length(g), "). Use a function of the ",
           "final digit when decimal places vary.", call. = FALSE)
    }
    w <- g
  } else {
    stop("`g` must be NULL, a function of the final digit, or a numeric ",
         "vector of grid weights.", call. = FALSE)
  }
  if (anyNA(w) || any(!is.finite(w)) || any(w < 0) || sum(w) <= 0) {
    stop("`g` weights must be finite, non-negative and not all zero.",
         call. = FALSE)
  }
  w / sum(w)
}

# Internal: Poisson-binomial pmf of the number of successes among independent
# Bernoulli(prob_i) trials, by dynamic programming. Element j + 1 is P(F = j).
.poisbinom_pmf <- function(prob) {
  P <- 1
  for (p in prob) P <- c(P * (1 - p), 0) + c(0, P * p)
  P
}

# Internal: tails P(F' <= f) and P(F' >= f) of a Poisson binomial.
.poisbinom_tails <- function(prob, f) {
  P <- .poisbinom_pmf(prob)
  c(lower = min(1, sum(P[seq_len(f + 1)])),
    upper = min(1, sum(P[(f + 1):length(P)])))
}

# Internal: mixture log-likelihood and its maximiser on [0, 1].
.lambda_loglik <- function(lambda, f, q) {
  # sum over failures only: 0 * log(0) is NaN in R, but the term is absent
  sum(log(lambda * q[f == 1])) + sum(log1p(-lambda * q[f == 0]))
}

.lambda_mle <- function(f, q) {
  F <- sum(f)
  if (F == 0) return(0)
  if (F == length(f)) return(1)
  score <- function(l) F / l - sum((1 - f) * q / (1 - l * q))
  # the log-likelihood is concave, so a non-negative score at 1 means the
  # maximum on [0, 1] is at the boundary
  if (score(1) >= 0) return(1)
  stats::uniroot(score, c(1e-12, 1), tol = 1e-12)$root
}

# Internal: likelihood support interval {lambda : l(lambda) >= l(mle) - log c}.
.lambda_support <- function(f, q, mle, c) {
  ll <- function(l) .lambda_loglik(l, f, q)
  thr <- ll(mle) - log(c)
  h <- function(l) ll(l) - thr
  lower <- if (mle == 0 || h(1e-12) >= 0) 0 else
    stats::uniroot(h, c(1e-12, mle), tol = 1e-12)$root
  upper <- if (mle == 1 || h(1) >= 0) 1 else
    stats::uniroot(h, c(mle, 1), tol = 1e-12)$root
  c(lower, upper)
}

# Internal: confidence interval for lambda by inverting Poisson-binomial tails
# with rates lambda * q_i.
.lambda_ci <- function(f, q, conf_level) {
  F <- sum(f)
  a2 <- (1 - conf_level) / 2
  ge <- function(l) .poisbinom_tails(l * q, F)[["upper"]] - a2
  le <- function(l) .poisbinom_tails(l * q, F)[["lower"]] - a2
  lower <- if (F == 0) 0 else if (ge(1) < 0) 1 else
    stats::uniroot(ge, c(0, 1), tol = 1e-10)$root
  upper <- if (le(1) > 0) 1 else
    stats::uniroot(le, c(0, 1), tol = 1e-10)$root
  c(lower, upper)
}

# Internal: validate input and derive digits and items.
.grim_agg_input <- function(data, items, percent) {
  if (isTRUE(percent)) {
    stop("Percentages are out of scope for grim_aggregate(): the pass-rate ",
         "grid assumes means, not percentages.", call. = FALSE)
  }
  if (!is.data.frame(data) || !all(c("x", "n") %in% names(data))) {
    stop("`data` must be a data frame with columns `x` (reported means, as ",
         "character) and `n` (sample sizes).", call. = FALSE)
  }
  if (!is.character(data$x)) {
    stop("`x` must be character, so that trailing zeros define the ",
         "reporting precision (\"5.10\", not 5.1).", call. = FALSE)
  }
  if (anyNA(data$x) || anyNA(data$n)) {
    stop("`x` and `n` must not contain missing values.", call. = FALSE)
  }
  n <- data$n
  if (!is.numeric(n) || any(n < 1) || any(n != round(n))) {
    stop("`n` must contain positive whole numbers.", call. = FALSE)
  }
  m <- if ("items" %in% names(data)) data$items else
    rep_len(items, nrow(data))
  if (!is.numeric(m) || any(m < 1) || any(m != round(m))) {
    stop("`items` must contain positive whole numbers.", call. = FALSE)
  }
  d <- scrutiny::decimal_places(data$x)
  if (any(d > 6)) {
    stop("Values reported to more than 6 decimal places are not supported ",
         "(the pass-rate grid has 10^d points).", call. = FALSE)
  }
  list(x = data$x, n = as.integer(n), m = as.integer(m), d = as.integer(d))
}

# Internal: resolve eps to one baseline failure probability per input row.
.grim_agg_eps <- function(eps, q, n, d) {
  if (is.function(eps)) eps <- eps(q = q, n = n, d = d)
  if (!is.numeric(eps) || !(length(eps) %in% c(1L, length(q)))) {
    stop("`eps` must be a single probability, one per row of `data`, or a ",
         "function(q, n, d) returning one of those.", call. = FALSE)
  }
  eps <- rep_len(eps, length(q))
  if (anyNA(eps) || any(eps < 0 | eps > 1)) {
    stop("`eps` must lie in [0, 1].", call. = FALSE)
  }
  eps
}

# Internal: the summary tibble rows, in ladder order.
.grim_agg_row <- function(step, quantity, value = NA_real_, lower = NA_real_,
                          upper = NA_real_, assumptions, interpretation) {
  tibble::tibble(step = step, quantity = quantity, value = value,
                 lower = lower, upper = upper, assumptions = assumptions,
                 interpretation = interpretation)
}

#' Pool many GRIM tests into quantities tagged with their assumptions
#'
#' `grim_aggregate()` GRIM-tests every reported mean in `data` and pools the
#' results into an *assumption ladder*: a short list of quantities, each
#' reported next to the assumptions it rests on, so that you can stop at
#' whichever rung you are willing to defend. The deterministic and
#' assumption-light quantities are the headline; the model-based ones are
#' secondary.
#'
#' @section Assumptions:
#' * **A0** (no authentic failures): correct values never fail GRIM. Holds
#'   only if the rounding convention, items and extraction are right.
#' * **U** (uniform counterfactual): a non-genuine value is equally likely to
#'   be any value on the reporting grid (or distributed as `g`), so it fails
#'   with probability \eqn{q_i}.
#' * **I** (independence): values pass or fail independently. Violated by a
#'   shared wrong N, or by copied or derived values.
#' * **M** (mixture): each value is independently non-genuine with
#'   probability \eqn{\lambda}, otherwise correct.
#' * **B** (authentic baseline): an authentic value fails with probability
#'   \eqn{\varepsilon_i}, supplied by you through `eps`.
#'
#' @section The ladder:
#' 0. **GRIM result** per value, \eqn{c_i}, from [scrutiny::grim_map()];
#'    \eqn{f_i = 1 - c_i}. Assumptions: none.
#' 1. **Pass rate** \eqn{p_i}: the share of decimal endings on the reporting
#'    grid that are consistent with \eqn{N_i m_i}, found by enumerating one
#'    unit interval (exact, since consistency has period 1). \eqn{q_i = 1 -
#'    p_i}. A property of the design, not of the reported value. Assumptions:
#'    none.
#' 2. **k and F**: \eqn{k} values with \eqn{q_i > 0} are informative; the
#'    rest cannot fail and are excluded from everything below (their count is
#'    reported). A value that fails although \eqn{q_i = 0} is impossible — an
#'    extraction or reporting error — and is warned about, flagged and
#'    excluded. \eqn{F} is the number of failures among informative values.
#'    Assumptions: none.
#' 3. **Partial identification**: every failure is non-genuine and any pass
#'    could be either, so the share \eqn{\theta} of non-genuine informative
#'    values satisfies \eqn{\theta \in [F/k, 1]}. The upper bound is always
#'    1. Assumptions: A0.
#' 4. **Expected failures** \eqn{E = \sum q_i} if every value were random: how
#'    much information the set carries. Assumptions: A0, U.
#' 5. **Observed-to-expected ratio** \eqn{O/E = F/E}. Needs no independence;
#'    not capped. 0: no failures; below 1: fewer failures than random values
#'    produce; 1: as many; above 1: more than random values can produce, i.e.
#'    the uniform model misfits (shared wrong N, wrong rounding or items,
#'    extraction errors). With small E it moves in steps of 1/E, so always
#'    read it next to E. Assumptions: A0, U.
#' 6. **Distribution of F under all-random**: Poisson-binomial\eqn{(q_1,
#'    \dots, q_k)}. `p_fewer_than_random` \eqn{= P(F' \le F)} is evidence
#'    against "every value random"; `p_more_than_random` \eqn{= P(F' \ge F)}
#'    is a misfit diagnostic. These are randomisation-test statements about a
#'    specified mechanism, not about a population. Assumptions: A0, U, I.
#' 7. **Mixture rate** \eqn{\lambda}, with \eqn{P(f_i = 1) = \lambda q_i}.
#'    The moment estimator \eqn{\tilde\lambda = F/E} is unbiased under M but
#'    unbounded. The MLE \eqn{\hat\lambda} is restricted to \eqn{[0, 1]}; when
#'    all \eqn{q_i} are equal it is \eqn{\min(1, O/E)}. The cap biases
#'    \eqn{\hat\lambda} downward when \eqn{k q} is small. Assumptions: A0, U,
#'    I, M.
#' 8. **Intervals for** \eqn{\lambda}: likelihood support intervals
#'    \eqn{\{\lambda : \ell(\lambda) \ge \ell(\hat\lambda) - \log c\}} for
#'    each `support_c` (Royall's 8 and 32); optionally (`ci = TRUE`) a
#'    confidence interval from inverting Poisson-binomial tails. \eqn{\lambda}
#'    is a *process* parameter: suited to pooled analyses across papers, not
#'    recommended as the headline for a single paper. Assumptions: A0, U, I,
#'    M.
#' 9. **Authentic baseline** (only with `eps`): \eqn{E_A = \sum
#'    \varepsilon_i}, the standardized inconsistency ratio \eqn{SIR = F /
#'    E_A}, and `p_vs_baseline` \eqn{= P(F' \ge F)} under
#'    Poisson-binomial\eqn{(\varepsilon_1, \dots, \varepsilon_k)}. If
#'    \eqn{\varepsilon_i = \tau q_i} then \eqn{SIR = (O/E)/\tau}. Assumptions:
#'    B (+ I for the tail probability). No \eqn{\varepsilon} enters steps
#'    0-8.
#'
#' @param data A data frame with columns `x` (reported means, **character**,
#'   so that trailing zeros define the reporting precision) and `n` (sample
#'   sizes), and optionally `items`.
#' @param items Number of items per mean, used when `data` has no `items`
#'   column. Default 1.
#' @param rounding Rounding rule, passed to [scrutiny::grim_map()] and used
#'   identically for observed consistency and for pass rates. Default
#'   `"up_or_down"`.
#' @param g Counterfactual distribution of non-genuine values over the
#'   reporting grid, for sensitivity analyses. `NULL` (default) is uniform.
#'   Either a **function of the final digit**, called once on `0:9` and
#'   returning 10 non-negative weights (grid value \eqn{j/10^d} gets the
#'   weight of its final digit \eqn{j \bmod 10}); or a **numeric vector of
#'   length \eqn{10^d}** giving one weight per grid value \eqn{j/10^d}, \eqn{j
#'   = 0, \dots, 10^d - 1} (only when every value has the same \eqn{d}).
#'   Weights are normalised to sum to 1.
#' @param eps Optional authentic-baseline failure probabilities for step 9: a
#'   single probability, one per row of `data`, or a `function(q, n, d)`
#'   returning either.
#' @param support_c Likelihood ratios defining the support intervals for
#'   \eqn{\lambda}. Default `c(8, 32)`.
#' @param ci Logical; also report a confidence interval for \eqn{\lambda}.
#'   Default `FALSE`.
#' @param conf_level Confidence level for that interval. Default 0.90.
#' @param percent Must be `FALSE`: percentages are out of scope.
#'
#' @return An object of class `grim_aggregate`: a list with
#'   * `summary`: a long tibble, one row per quantity in ladder order, with
#'     columns `step`, `quantity`, `value` (point quantities), `lower` and
#'     `upper` (intervals), `assumptions` and `interpretation`;
#'   * `data`: the input with added columns `digits`, `consistency`, `f`,
#'     `p_pass`, `q`, `informative` and `flag_impossible`;
#'   * `settings`: the rounding rule and counterfactual used.
#' @seealso [grim_aggregate_by()] for one summary per group,
#'   [grim_sensitivity()] for a sweep over counterfactuals.
#' @examples
#' dat <- data.frame(
#'   x = c("3.41", "2.87", "4.12", "3.05", "2.66", "3.90"),
#'   n = c(28, 31, 45, 22, 37, 19)
#' )
#' res <- grim_aggregate(dat)
#' res
#' res$summary
#' @export
grim_aggregate <- function(data, items = 1, rounding = "up_or_down",
                           g = NULL, eps = NULL, support_c = c(8, 32),
                           ci = FALSE, conf_level = 0.90, percent = FALSE) {
  inp <- .grim_agg_input(data, items, percent)
  x <- inp$x; n <- inp$n; m <- inp$m; d <- inp$d

  # 0. GRIM result per value
  consistency <- .grim_observed(x, n, m, d, rounding)
  f <- as.integer(!consistency)

  # 1. pass rates, under g and (for the impossibility flag) uniform
  grids <- .grim_grid(d, n, m, rounding)
  p_unif <- vapply(grids, mean, numeric(1))
  p_pass <- vapply(seq_along(grids), function(i) {
    sum(.grid_weights(g, d[i]) * grids[[i]])
  }, numeric(1))
  q <- pmax(0, 1 - p_pass)

  # 2. informative values; impossible failures
  tol <- 1e-12
  flag_impossible <- f == 1L & p_unif >= 1 - tol
  if (any(flag_impossible)) {
    warning(sum(flag_impossible), " value(s) fail GRIM although every ",
            "ending on the grid is consistent with their N (rows ",
            paste(which(flag_impossible), collapse = ", "), "). This ",
            "indicates an extraction or reporting error; they are flagged ",
            "and excluded.", call. = FALSE)
  }
  bad_g <- f == 1L & !flag_impossible & q <= tol
  if (any(bad_g)) {
    stop("`g` gives zero weight to every inconsistent ending for failing ",
         "row(s) ", paste(which(bad_g), collapse = ", "), ", so the ",
         "counterfactual cannot produce the observed failure.", call. = FALSE)
  }
  informative <- q > tol & !flag_impossible

  out_data <- tibble::as_tibble(data)
  out_data$digits <- d
  out_data$consistency <- consistency
  out_data$f <- f
  out_data$p_pass <- p_pass
  out_data$q <- q
  out_data$informative <- informative
  out_data$flag_impossible <- flag_impossible

  fi <- f[informative]; qi <- q[informative]
  k <- length(fi); F <- sum(fi)

  rows <- list(
    .grim_agg_row(0, "n_values", nrow(data), assumptions = "none",
                  interpretation = "values GRIM-tested"),
    .grim_agg_row(1, "mean_q", mean(q), assumptions = "none",
                  interpretation = "mean fail rate q over tested values"),
    .grim_agg_row(2, "n_uninformative", sum(!informative & !flag_impossible),
                  assumptions = "none",
                  interpretation = "values that cannot fail (q = 0); excluded"),
    .grim_agg_row(2, "n_flagged_impossible", sum(flag_impossible),
                  assumptions = "none",
                  interpretation = "failures with q = 0: extraction errors; excluded"),
    .grim_agg_row(2, "k", k, assumptions = "none",
                  interpretation = "informative values"),
    .grim_agg_row(2, "F", F, assumptions = "none",
                  interpretation = "GRIM failures among informative values")
  )

  if (k == 0) {
    summary <- do.call(rbind, rows)
    return(structure(
      list(summary = summary, data = out_data,
           settings = list(rounding = rounding, g = g, eps = !is.null(eps))),
      class = "grim_aggregate"))
  }

  E <- sum(qi)
  tails <- .poisbinom_tails(qi, F)
  mle <- .lambda_mle(fi, qi)

  rows <- c(rows, list(
    .grim_agg_row(3, "theta", lower = F / k, upper = 1, assumptions = "A0",
                  interpretation = "share of informative values that are non-genuine"),
    .grim_agg_row(4, "E", E, assumptions = "A0, U",
                  interpretation = "failures expected if every value were random"),
    .grim_agg_row(5, "O_E", F / E, assumptions = "A0, U",
                  interpretation = "0 none; <1 fewer than random; 1 as random; >1 misfit"),
    .grim_agg_row(6, "var_F", sum(qi * (1 - qi)), assumptions = "A0, U, I",
                  interpretation = "variance of F if every value were random"),
    .grim_agg_row(6, "p_fewer_than_random", tails[["lower"]],
                  assumptions = "A0, U, I",
                  interpretation = "P(F' <= F): evidence against all-random"),
    .grim_agg_row(6, "p_more_than_random", tails[["upper"]],
                  assumptions = "A0, U, I",
                  interpretation = "P(F' >= F): misfit of the random model"),
    .grim_agg_row(7, "lambda_moment", F / E, assumptions = "A0, U, I, M",
                  interpretation = "non-genuine rate, moment estimate (unbounded)"),
    .grim_agg_row(7, "lambda_mle", mle, assumptions = "A0, U, I, M",
                  interpretation = "non-genuine rate, MLE on [0, 1]")
  ))
  for (cc in support_c) {
    si <- .lambda_support(fi, qi, mle, cc)
    rows <- c(rows, list(.grim_agg_row(
      8, paste0("lambda_support_", cc), lower = si[1], upper = si[2],
      assumptions = "A0, U, I, M",
      interpretation = paste0("1/", cc, " likelihood support interval"))))
  }
  if (isTRUE(ci)) {
    lc <- .lambda_ci(fi, qi, conf_level)
    rows <- c(rows, list(.grim_agg_row(
      8, paste0("lambda_ci_", round(100 * conf_level)), lower = lc[1],
      upper = lc[2], assumptions = "A0, U, I, M",
      interpretation = paste0(round(100 * conf_level),
                              "% confidence interval (tail inversion)"))))
  }

  if (!is.null(eps)) {
    ei <- .grim_agg_eps(eps, q, n, d)[informative]
    E_A <- sum(ei)
    rows <- c(rows, list(
      .grim_agg_row(9, "E_A", E_A, assumptions = "B",
                    interpretation = "failures expected from authentic values"),
      .grim_agg_row(9, "SIR", F / E_A, assumptions = "B",
                    interpretation = "standardized inconsistency ratio F / E_A"),
      .grim_agg_row(9, "p_vs_baseline", .poisbinom_tails(ei, F)[["upper"]],
                    assumptions = "B, I",
                    interpretation = "P(F' >= F) if every value were authentic")
    ))
    out_data$eps <- .grim_agg_eps(eps, q, n, d)
  }

  structure(
    list(summary = do.call(rbind, rows), data = out_data,
         settings = list(rounding = rounding, g = g, eps = !is.null(eps))),
    class = "grim_aggregate")
}

#' @export
print.grim_aggregate <- function(x, digits = 3, ...) {
  s <- x$summary
  cf <- if (is.null(x$settings$g)) "uniform" else "custom g"
  cat("GRIM aggregation  (rounding = ", x$settings$rounding,
      ", counterfactual = ", cf, ")\n", sep = "")
  fmt <- function(v) formatC(v, digits = digits, format = "fg", flag = "#")
  shown <- ifelse(
    is.na(s$value),
    paste0("[", fmt(s$lower), ", ", fmt(s$upper), "]"),
    ifelse(s$quantity %in% c("n_values", "n_uninformative",
                             "n_flagged_impossible", "k", "F"),
           as.character(s$value), fmt(s$value)))
  shown <- gsub("\\s+", "", shown)
  tab <- data.frame(step = s$step, quantity = s$quantity, value = shown,
                    assumptions = s$assumptions,
                    interpretation = s$interpretation)
  w_q <- max(nchar(tab$quantity)); w_v <- max(nchar(tab$value))
  w_a <- max(nchar(tab$assumptions))
  last <- -1
  for (i in seq_len(nrow(tab))) {
    step <- if (tab$step[i] != last) format(tab$step[i], width = 2) else "  "
    last <- tab$step[i]
    cat(step, "  ", formatC(tab$quantity[i], width = -w_q), "  ",
        formatC(tab$value[i], width = w_v), "  ",
        formatC(tab$assumptions[i], width = -w_a), "  ",
        tab$interpretation[i], "\n", sep = "")
  }
  if (any(s$step == 8)) {
    cat("Note: lambda is a process parameter; prefer steps 0-6 as the ",
        "headline for a single paper.\n", sep = "")
  }
  invisible(x)
}

# Internal: one-row wide version of a summary: points as `quantity`,
# intervals as `quantity_lower` / `quantity_upper`.
.grim_agg_wide <- function(s) {
  out <- list()
  for (i in seq_len(nrow(s))) {
    if (is.na(s$value[i]) && !is.na(s$lower[i])) {
      out[[paste0(s$quantity[i], "_lower")]] <- s$lower[i]
      out[[paste0(s$quantity[i], "_upper")]] <- s$upper[i]
    } else {
      out[[s$quantity[i]]] <- s$value[i]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

#' GRIM aggregation by group
#'
#' Runs [grim_aggregate()] separately within each group of `data` (e.g. one
#' group per paper) and returns one row per group.
#'
#' @param data A data frame as for [grim_aggregate()].
#' @param by Character vector of grouping column names.
#' @param ... Further arguments passed to [grim_aggregate()].
#' @return A wide tibble with the grouping columns followed by one column per
#'   point quantity and `_lower` / `_upper` columns per interval. Groups with
#'   no informative values have `NA` for everything past `F`.
#' @seealso [grim_aggregate()]
#' @examples
#' dat <- data.frame(
#'   paper = rep(c("A", "B"), each = 4),
#'   x = c("3.41", "2.87", "4.12", "3.05", "2.60", "3.93", "1.14", "2.29"),
#'   n = c(28, 31, 45, 22, 37, 19, 26, 33)
#' )
#' grim_aggregate_by(dat, by = "paper")
#' @export
grim_aggregate_by <- function(data, by, ...) {
  if (!all(by %in% names(data))) {
    stop("`by` must name columns of `data`.", call. = FALSE)
  }
  key <- interaction(data[by], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), key)
  wide <- lapply(groups, function(i) {
    res <- grim_aggregate(data[i, , drop = FALSE], ...)
    cbind(data[i[1], by, drop = FALSE], .grim_agg_wide(res$summary))
  })
  # groups with k = 0 lack the later columns
  cols <- unique(unlist(lapply(wide, names)))
  wide <- lapply(wide, function(w) {
    w[setdiff(cols, names(w))] <- NA_real_
    w[cols]
  })
  out <- do.call(rbind, wide)
  rownames(out) <- NULL
  tibble::as_tibble(out)
}

#' GRIM aggregation across counterfactual digit distributions
#'
#' Runs [grim_aggregate()] under each digit distribution in `g_list` (plus
#' the uniform one) so you can see whether conclusions hold across
#' counterfactuals, e.g. avoidance of, or heaping on, final digits 0 and 5.
#'
#' @param data A data frame as for [grim_aggregate()].
#' @param g_list A named list of counterfactuals, each in any form accepted
#'   by the `g` argument of [grim_aggregate()]. A `uniform` entry (`NULL`) is
#'   added first unless one is already present.
#' @param ... Further arguments passed to [grim_aggregate()].
#' @return A tibble with one row per counterfactual and columns
#'   `counterfactual`, `k`, `F`, `E`, `O_E`, `p_fewer_than_random` and
#'   `p_more_than_random`.
#' @seealso [grim_aggregate()]
#' @examples
#' dat <- data.frame(
#'   x = c("3.41", "2.87", "4.12", "3.05", "2.66", "3.90"),
#'   n = c(28, 31, 45, 22, 37, 19)
#' )
#' grim_sensitivity(dat, list(
#'   avoid_0_5 = function(d) ifelse(d %in% c(0, 5), 0.5, 1),
#'   heap_0_5  = function(d) ifelse(d %in% c(0, 5), 3, 1)
#' ))
#' @export
grim_sensitivity <- function(data, g_list, ...) {
  if (!is.list(g_list) || is.null(names(g_list)) || any(names(g_list) == "")) {
    stop("`g_list` must be a named list of counterfactuals.", call. = FALSE)
  }
  if (!"uniform" %in% names(g_list)) {
    g_list <- c(list(uniform = NULL), g_list)
  }
  keep <- c("k", "F", "E", "O_E", "p_fewer_than_random", "p_more_than_random")
  rows <- lapply(names(g_list), function(nm) {
    s <- grim_aggregate(data, g = g_list[[nm]], ...)$summary
    v <- stats::setNames(s$value[match(keep, s$quantity)], keep)
    tibble::as_tibble(c(list(counterfactual = nm), as.list(v)))
  })
  do.call(rbind, rows)
}
