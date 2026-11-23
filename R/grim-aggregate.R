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

# Grid-consistency vectors, keyed by "digits|n|items|rounding|threshold|
# symmetric". Each entry is the logical vector c(j / 10^d; n, items, d) for
# j = 0, ..., 10^d - 1. GRIM consistency is periodic with period 1 in the
# mean, so one unit interval is exact. The cache persists for the session:
# pass rates are a property of the design, not of the reported value.
grim_pass_grid_cache <- new.env(parent = emptyenv())

# Internal: scrutiny::grim_map() consistency for numeric means, one call for
# all rows. `set` holds rounding, threshold and symmetric. Returns a logical
# vector in row order.
grim_consistency <- function(x, n, items, digits, set) {
  df <- tibble::tibble(x = x, n = n, items_averaged = items, digits_x = digits)
  out <- scrutiny::grim_map(
    df,
    data_type = "integer",
    rounding = set$rounding,
    threshold = set$threshold,
    symmetric = set$symmetric
  )
  as.logical(out$consistency)
}

# Internal: observed verdicts; a separate binding so tests can stub them
# without touching the pass-rate grid.
grim_observed <- function(...) grim_consistency(...)

# Internal: the grid-consistency vector for every unique (digits, n, items),
# computed in one grim_map() call and cached.
grim_pass_grid <- function(digits, n, items, set) {
  key <- paste(
    digits,
    n,
    items,
    set$rounding,
    set$threshold,
    set$symmetric,
    sep = "|"
  )
  cached <- purrr::map_lgl(key, function(key_current) {
    exists(key_current, envir = grim_pass_grid_cache, inherits = FALSE)
  })
  todo <- unique(tibble::tibble(key, digits, n, items)[!cached, ])
  if (nrow(todo) > 0) {
    grid_size <- 10^todo$digits
    todo_row <- rep(seq_len(nrow(todo)), grid_size)
    j <- sequence(grid_size) - 1
    consistency <- grim_consistency(
      j / grid_size[todo_row],
      todo$n[todo_row],
      todo$items[todo_row],
      todo$digits[todo_row],
      set
    )
    for (i in seq_len(nrow(todo))) {
      assign(
        todo$key[i],
        consistency[todo_row == i],
        envir = grim_pass_grid_cache
      )
    }
  }
  purrr::map(key, function(key_current) {
    get(key_current, envir = grim_pass_grid_cache, inherits = FALSE)
  })
}

# Internal: grid weights for decimal count d under counterfactual g.
# NULL is uniform; a function maps final digits 0:9 to weights; a numeric
# vector gives one weight per grid value j / 10^d, j = 0, ..., 10^d - 1.
grid_weights <- function(g, d) {
  grid_size <- 10^d
  if (is.null(g)) {
    weights <- rep(1, grid_size)
  } else if (is.function(g)) {
    digit_weights <- g(0:9)
    if (!is.numeric(digit_weights) || length(digit_weights) != 10L) {
      cli::cli_abort(
        "{.arg g} as a function must return 10 weights, one per final digit
        0-9, when called on 0:9.",
        call = NULL
      )
    }
    weights <- digit_weights[(0:(grid_size - 1)) %% 10 + 1]
  } else if (is.numeric(g)) {
    if (length(g) != grid_size) {
      cli::cli_abort(
        "{.arg g} as a numeric vector must have one weight per grid value,
        i.e. length 10^d = {grid_size} for values reported to {d} decimals
        (got length {length(g)}). Use a function of the final digit when
        decimal places vary.",
        call = NULL
      )
    }
    weights <- g
  } else {
    cli::cli_abort(
      "{.arg g} must be NULL, a function of the final digit, or a numeric
      vector of grid weights.",
      call = NULL
    )
  }
  if (
    anyNA(weights) ||
      any(!is.finite(weights)) ||
      any(weights < 0) ||
      sum(weights) <= 0
  ) {
    cli::cli_abort(
      "{.arg g} weights must be finite, non-negative and not all zero.",
      call = NULL
    )
  }
  weights / sum(weights)
}

# Internal: Poisson-binomial pmf of the number of successes among independent
# Bernoulli(prob_i) trials, by dynamic programming. Element j + 1 is P(F = j).
poisbinom_pmf <- function(prob) {
  P <- 1
  for (p in prob) {
    P <- c(P * (1 - p), 0) + c(0, P * p)
  }
  P
}

# Internal: tails P(F' <= f) and P(F' >= f) of a Poisson binomial.
poisbinom_tails <- function(prob, f) {
  P <- poisbinom_pmf(prob)
  c(
    lower = min(1, sum(P[seq_len(f + 1)])),
    upper = min(1, sum(P[(f + 1):length(P)]))
  )
}

# Internal: mixture log-likelihood and its maximiser on [0, 1].
lambda_loglik <- function(lambda, f, q) {
  # sum over failures only: 0 * log(0) is NaN in R, but the term is absent
  sum(log(lambda * q[f == 1])) + sum(log1p(-lambda * q[f == 0]))
}

lambda_mle <- function(f, q) {
  F <- sum(f)
  if (F == 0) {
    return(0)
  }
  if (F == length(f)) {
    return(1)
  }
  score <- function(l) F / l - sum((1 - f) * q / (1 - l * q))
  # the log-likelihood is concave, so a non-negative score at 1 means the
  # maximum on [0, 1] is at the boundary
  if (score(1) >= 0) {
    return(1)
  }
  stats::uniroot(score, c(1e-12, 1), tol = 1e-12)$root
}

# Internal: likelihood support interval {lambda : l(lambda) >= l(mle) - log c}.
lambda_support <- function(f, q, mle, c) {
  threshold <- lambda_loglik(mle, f, q) - log(c)
  excess <- function(l) lambda_loglik(l, f, q) - threshold
  lower <- if (mle == 0 || excess(1e-12) >= 0) {
    0
  } else {
    stats::uniroot(excess, c(1e-12, mle), tol = 1e-12)$root
  }
  upper <- if (mle == 1 || excess(1) >= 0) {
    1
  } else {
    stats::uniroot(excess, c(mle, 1), tol = 1e-12)$root
  }
  c(lower, upper)
}

# Internal: confidence interval for lambda by inverting Poisson-binomial tails
# with rates lambda * q_i.
lambda_ci <- function(f, q, conf_level) {
  F <- sum(f)
  alpha_half <- (1 - conf_level) / 2
  upper_gap <- function(l) poisbinom_tails(l * q, F)[["upper"]] - alpha_half
  lower_gap <- function(l) poisbinom_tails(l * q, F)[["lower"]] - alpha_half
  if (F > 0 && upper_gap(1) < 0) {
    cli::cli_warn(
      "Even lambda = 1 makes {F} or more failures less likely than
      (1 - conf_level) / 2: the mixture model is rejected, so the confidence
      interval is NA."
    )
    return(c(NA_real_, NA_real_))
  }
  lower <- if (F == 0) {
    0
  } else {
    stats::uniroot(upper_gap, c(0, 1), tol = 1e-10)$root
  }
  upper <- if (lower_gap(1) > 0) {
    1
  } else {
    stats::uniroot(lower_gap, c(0, 1), tol = 1e-10)$root
  }
  c(lower, upper)
}

# Internal: one per-row helper from either a `data` column or an argument,
# never both, never neither. `old` names a column renamed in this version.
grim_aggregate_helper <- function(data, name, arg, old = NULL) {
  in_data <- name %in% names(data)
  if (in_data && !is.null(arg)) {
    cli::cli_abort(
      "{.arg {name}} is both a column of {.arg data} and an argument. Give it
      one way only.",
      call = NULL
    )
  }
  if (!in_data && is.null(arg)) {
    hint <- if (!is.null(old) && old %in% names(data)) {
      cli::format_inline(
        " {.arg data} has a column {.arg {old}}: rename it to {.arg {name}}."
      )
    } else {
      ""
    }
    cli::cli_abort(
      "{.arg {name}} is required, as a column of {.arg data} or as an
      argument.{hint}",
      call = NULL
    )
  }
  values <- if (in_data) data[[name]] else arg
  if (!(length(values) %in% c(1L, nrow(data)))) {
    cli::cli_abort(
      "{.arg {name}} must have length 1 or one value per row of {.arg data}.",
      call = NULL
    )
  }
  rep_len(values, nrow(data))
}

# Internal: stop unless `values` holds whole numbers >= `min`, none missing.
grim_aggregate_check_whole <- function(values, name, min) {
  if (
    !is.numeric(values) ||
      anyNA(values) ||
      any(!is.finite(values)) ||
      any(values < min) ||
      any(values != round(values))
  ) {
    cli::cli_abort(
      "{.arg {name}} must contain whole numbers >= {min}, with no missing or
      infinite values.",
      call = NULL
    )
  }
}

# Internal: validate input; return x, n, items and digits, one per row.
grim_aggregate_input <- function(data, digits_x, items_averaged) {
  if (!is.data.frame(data) || !all(c("x", "n") %in% names(data))) {
    cli::cli_abort(
      "{.arg data} must be a data frame with columns {.arg x} (reported means)
      and {.arg n} (sample sizes).",
      call = NULL
    )
  }
  if (nrow(data) == 0) {
    cli::cli_abort("{.arg data} has no rows.", call = NULL)
  }
  if (!is.numeric(data$x)) {
    cli::cli_abort(
      "{.arg x} must be numeric, with the reported decimal places in
      {.arg digits_x} (a column or an argument). For means stored as strings,
      {.fn scrutiny::decimal_places} counts them.",
      call = NULL
    )
  }
  if (anyNA(data$x) || any(!is.finite(data$x))) {
    bad_rows <- which(!is.finite(data$x))
    cli::cli_abort(
      "{.arg x} must not contain missing or infinite values (row(s)
      {bad_rows}).",
      call = NULL
    )
  }
  d <- grim_aggregate_helper(data, "digits_x", digits_x)
  m <- grim_aggregate_helper(
    data,
    "items_averaged",
    items_averaged,
    old = "items"
  )
  grim_aggregate_check_whole(data$n, "n", 1)
  grim_aggregate_check_whole(d, "digits_x", 0)
  grim_aggregate_check_whole(m, "items_averaged", 1)
  if (any(d > 6)) {
    cli::cli_abort(
      "Values reported to more than 6 decimal places are not supported (the
      pass-rate grid has 10^d points).",
      call = NULL
    )
  }
  list(x = data$x, n = as.integer(data$n), m = as.integer(m), d = as.integer(d))
}

# Internal: resolve eps to one baseline failure probability per input row.
grim_aggregate_eps <- function(eps, q, n, d) {
  if (is.function(eps)) {
    eps <- eps(q = q, n = n, d = d)
  }
  if (!is.numeric(eps) || !(length(eps) %in% c(1L, length(q)))) {
    cli::cli_abort(
      "{.arg eps} must be a single probability, one per row of {.arg data}, or
      a {.code function(q, n, d)} returning one of those.",
      call = NULL
    )
  }
  eps <- rep_len(eps, length(q))
  if (anyNA(eps) || any(eps < 0 | eps > 1)) {
    cli::cli_abort("{.arg eps} must lie in [0, 1].", call = NULL)
  }
  eps
}

# Internal: the summary tibble rows, in ladder order.
grim_aggregate_row <- function(
  step,
  quantity,
  value = NA_real_,
  lower = NA_real_,
  upper = NA_real_,
  assumptions,
  interpretation
) {
  tibble::tibble(
    step = step,
    quantity = quantity,
    value = value,
    lower = lower,
    upper = upper,
    assumptions = assumptions,
    interpretation = interpretation
  )
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
#' @param data A data frame with columns `x` (reported means, numeric) and
#'   `n` (sample sizes), and optionally `digits_x` and `items_averaged`.
#' @param digits_x Decimal places of the reported means, as reported
#'   (`5.10` has 2). Required: give it here or as a column of `data`, not
#'   both. A single value or one per row.
#' @param ... Must be empty. Arguments after the dots must be named in full.
#' @param items_averaged Number of items averaged per mean. Required, as an
#'   argument or a column of `data`, like `digits_x`.
#' @param rounding,threshold,symmetric Rounding settings, passed to
#'   [scrutiny::grim_map()] and used identically for observed consistency
#'   and for pass rates. Defaults as in scrutiny.
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
#' @param support_c Likelihood ratios (each > 1) defining the support
#'   intervals for \eqn{\lambda}. Default `c(8, 32)`.
#' @param ci Logical; also report a confidence interval for \eqn{\lambda}.
#'   Default `FALSE`. If even \eqn{\lambda = 1} makes the observed `F`
#'   implausibly large, the model is rejected: the interval is `NA`, with a
#'   warning.
#' @param conf_level Confidence level for that interval, in (0, 1). Default
#'   0.90.
#'
#' @return An object of class `grim_aggregate`: a list with
#'   * `summary`: a long tibble, one row per quantity in ladder order, with
#'     columns `step`, `quantity`, `value` (point quantities), `lower` and
#'     `upper` (intervals), `assumptions` and `interpretation`;
#'   * `data`: the input as a tibble with columns `digits_x`,
#'     `items_averaged`, `consistency`, `f`, `p_pass`, `q`, `informative`,
#'     `flag_impossible` and, with `eps`, `eps`;
#'   * `settings`: the rounding settings and counterfactual used.
#' @seealso [grim_aggregate_by()] for one summary per group,
#'   [grim_sensitivity()] for a sweep over counterfactuals.
#' @examples
#' dat <- data.frame(
#'   x = c(3.41, 2.87, 4.12, 3.05, 2.66, 3.90),
#'   n = c(28, 31, 45, 22, 37, 19)
#' )
#' res <- grim_aggregate(dat, digits_x = 2, items_averaged = 1)
#' res
#' res$summary
#' @export
grim_aggregate <- function(
  data,
  digits_x = NULL,
  ...,
  items_averaged = NULL,
  rounding = "up_or_down",
  threshold = 5,
  symmetric = FALSE,
  g = NULL,
  eps = NULL,
  support_c = c(8, 32),
  ci = FALSE,
  conf_level = 0.90
) {
  grim_aggregate_check_dots(...)
  if (
    !is.numeric(support_c) ||
      length(support_c) == 0 ||
      anyNA(support_c) ||
      any(support_c <= 1) ||
      any(!is.finite(support_c))
  ) {
    cli::cli_abort("{.arg support_c} must be finite likelihood ratios > 1.")
  }
  if (
    !is.numeric(conf_level) ||
      length(conf_level) != 1 ||
      is.na(conf_level) ||
      conf_level <= 0 ||
      conf_level >= 1
  ) {
    cli::cli_abort("{.arg conf_level} must be a single number in (0, 1).")
  }
  input <- grim_aggregate_input(data, digits_x, items_averaged)
  x <- input$x
  n <- input$n
  m <- input$m
  d <- input$d
  set <- list(rounding = rounding, threshold = threshold, symmetric = symmetric)

  # 0. GRIM result per value
  consistency <- grim_observed(x, n, m, d, set)
  f <- as.integer(!consistency)

  # 1. pass rates, under g and (for the impossibility flag) uniform
  grids <- grim_pass_grid(d, n, m, set)
  p_unif <- purrr::map_dbl(grids, mean)
  p_pass <- purrr::map2_dbl(grids, d, function(grid, d_current) {
    sum(grid_weights(g, d_current) * grid)
  })
  # snap floating-point residue (1 - 1 = 1.1e-16) to exact 0
  tol <- 1e-12
  q <- 1 - p_pass
  q[q < tol] <- 0

  # 2. informative values; impossible failures
  flag_impossible <- f == 1L & p_unif >= 1 - tol
  if (any(flag_impossible)) {
    impossible_rows <- which(flag_impossible)
    cli::cli_warn(
      "{length(impossible_rows)} value(s) fail GRIM although every ending on
      the grid is consistent with their N (rows {impossible_rows}). This
      indicates an extraction or reporting error; they are flagged and
      excluded."
    )
  }
  bad_g <- f == 1L & !flag_impossible & q == 0
  if (any(bad_g)) {
    cli::cli_abort(
      "{.arg g} gives zero weight to every inconsistent ending for failing
      row(s) {which(bad_g)}, so the counterfactual cannot produce the observed
      failure."
    )
  }
  informative <- q > 0 & !flag_impossible
  eps_all <- if (!is.null(eps)) grim_aggregate_eps(eps, q, n, d)

  out_data <- tibble::as_tibble(data)
  out_data$digits_x <- d
  out_data$items_averaged <- m
  out_data$consistency <- consistency
  out_data$f <- f
  out_data$p_pass <- p_pass
  out_data$q <- q
  out_data$informative <- informative
  out_data$flag_impossible <- flag_impossible
  if (!is.null(eps)) {
    out_data$eps <- eps_all
  }

  f_informative <- f[informative]
  q_informative <- q[informative]
  k <- length(f_informative)
  F <- sum(f_informative)

  rows <- list(
    grim_aggregate_row(
      0,
      "n_values",
      nrow(data),
      assumptions = "none",
      interpretation = "values GRIM-tested"
    ),
    grim_aggregate_row(
      1,
      "mean_q",
      mean(q),
      assumptions = "none",
      interpretation = "mean fail rate q over tested values"
    ),
    grim_aggregate_row(
      2,
      "n_uninformative",
      sum(!informative & !flag_impossible),
      assumptions = "none",
      interpretation = "values that cannot fail (q = 0); excluded"
    ),
    grim_aggregate_row(
      2,
      "n_flagged_impossible",
      sum(flag_impossible),
      assumptions = "none",
      interpretation = "failures with q = 0: extraction errors; excluded"
    ),
    grim_aggregate_row(
      2,
      "k",
      k,
      assumptions = "none",
      interpretation = "informative values"
    ),
    grim_aggregate_row(
      2,
      "F",
      F,
      assumptions = "none",
      interpretation = "GRIM failures among informative values"
    )
  )
  settings <- c(set, list(g = g, eps = !is.null(eps)))

  if (k == 0) {
    return(structure(
      list(
        summary = dplyr::bind_rows(rows),
        data = out_data,
        settings = settings
      ),
      class = "grim_aggregate"
    ))
  }

  E <- sum(q_informative)
  tails <- poisbinom_tails(q_informative, F)
  mle <- lambda_mle(f_informative, q_informative)

  rows <- c(
    rows,
    list(
      grim_aggregate_row(
        3,
        "theta",
        lower = F / k,
        upper = 1,
        assumptions = "A0",
        interpretation = "share of informative values that are non-genuine"
      ),
      grim_aggregate_row(
        4,
        "E",
        E,
        assumptions = "A0, U",
        interpretation = "failures expected if every value were random"
      ),
      grim_aggregate_row(
        5,
        "O_E",
        F / E,
        assumptions = "A0, U",
        interpretation = "0 none; <1 fewer than random; 1 as random; >1 misfit"
      ),
      grim_aggregate_row(
        6,
        "var_F",
        sum(q_informative * (1 - q_informative)),
        assumptions = "A0, U, I",
        interpretation = "variance of F if every value were random"
      ),
      grim_aggregate_row(
        6,
        "p_fewer_than_random",
        tails[["lower"]],
        assumptions = "A0, U, I",
        interpretation = "P(F' <= F): evidence against all-random"
      ),
      grim_aggregate_row(
        6,
        "p_more_than_random",
        tails[["upper"]],
        assumptions = "A0, U, I",
        interpretation = "P(F' >= F): misfit of the random model"
      ),
      grim_aggregate_row(
        7,
        "lambda_moment",
        F / E,
        assumptions = "A0, U, I, M",
        interpretation = "non-genuine rate, moment estimate (unbounded)"
      ),
      grim_aggregate_row(
        7,
        "lambda_mle",
        mle,
        assumptions = "A0, U, I, M",
        interpretation = "non-genuine rate, MLE on [0, 1]"
      )
    )
  )
  for (c_current in support_c) {
    interval <- lambda_support(f_informative, q_informative, mle, c_current)
    rows <- c(
      rows,
      list(grim_aggregate_row(
        8,
        paste0("lambda_support_", c_current),
        lower = interval[1],
        upper = interval[2],
        assumptions = "A0, U, I, M",
        interpretation = paste0("1/", c_current, " likelihood support interval")
      ))
    )
  }
  if (isTRUE(ci)) {
    interval <- lambda_ci(f_informative, q_informative, conf_level)
    percent <- format(100 * conf_level)
    rows <- c(
      rows,
      list(grim_aggregate_row(
        8,
        paste0("lambda_ci_", percent),
        lower = interval[1],
        upper = interval[2],
        assumptions = "A0, U, I, M",
        interpretation = paste0(
          percent,
          "% confidence interval (tail inversion)"
        )
      ))
    )
  }

  if (!is.null(eps)) {
    eps_informative <- eps_all[informative]
    E_A <- sum(eps_informative)
    rows <- c(
      rows,
      list(
        grim_aggregate_row(
          9,
          "E_A",
          E_A,
          assumptions = "B",
          interpretation = "failures expected from authentic values"
        ),
        grim_aggregate_row(
          9,
          "SIR",
          F / E_A,
          assumptions = "B",
          interpretation = "standardized inconsistency ratio F / E_A"
        ),
        grim_aggregate_row(
          9,
          "p_vs_baseline",
          poisbinom_tails(eps_informative, F)[["upper"]],
          assumptions = "B, I",
          interpretation = "P(F' >= F) if every value were authentic"
        )
      )
    )
  }

  structure(
    list(
      summary = dplyr::bind_rows(rows),
      data = out_data,
      settings = settings
    ),
    class = "grim_aggregate"
  )
}

# Internal: stop if `...` is not empty, pointing renamed arguments to their
# new names.
grim_aggregate_check_dots <- function(
  ...,
  .renamed = c(items = "items_averaged", percent = "")
) {
  if (...length() == 0) {
    return(invisible())
  }
  arg_names <- names(list(...))
  if (is.null(arg_names)) {
    arg_names <- rep("", ...length())
  }
  hint <- dplyr::if_else(
    arg_names %in% names(.renamed),
    dplyr::if_else(
      .renamed[arg_names] == "",
      paste0(" `", arg_names, "` was removed."),
      paste0(" `", arg_names, "` was renamed to `", .renamed[arg_names], "`.")
    ),
    ""
  )
  arg_list <- paste0("`", arg_names, "`", collapse = ", ")
  hint <- paste(hint, collapse = "")
  cli::cli_abort(
    "Unknown or unnamed argument(s) in `...`: {arg_list}.{hint}",
    call = rlang::caller_env()
  )
}
#' @export
print.grim_aggregate <- function(x, digits = 3, ...) {
  summary <- x$summary
  counterfactual <- if (is.null(x$settings$g)) "uniform" else "custom g"
  cat(
    "GRIM aggregation  (rounding = ",
    x$settings$rounding,
    ", counterfactual = ",
    counterfactual,
    ")\n",
    sep = ""
  )
  format_value <- function(value) {
    formatC(value, digits = digits, format = "fg", flag = "#")
  }
  shown <- dplyr::if_else(
    is.na(summary$value),
    paste0(
      "[",
      format_value(summary$lower),
      ", ",
      format_value(summary$upper),
      "]"
    ),
    dplyr::if_else(
      summary$quantity %in%
        c("n_values", "n_uninformative", "n_flagged_impossible", "k", "F"),
      as.character(summary$value),
      format_value(summary$value)
    )
  )
  shown <- gsub("\\s+", "", shown)
  width_quantity <- max(nchar(summary$quantity))
  width_value <- max(nchar(shown))
  width_assumptions <- max(nchar(summary$assumptions))
  last <- -1
  for (i in seq_len(nrow(summary))) {
    step <- if (summary$step[i] != last) {
      format(summary$step[i], width = 2)
    } else {
      "  "
    }
    last <- summary$step[i]
    cat(
      step,
      "  ",
      formatC(summary$quantity[i], width = -width_quantity),
      "  ",
      formatC(shown[i], width = width_value),
      "  ",
      formatC(summary$assumptions[i], width = -width_assumptions),
      "  ",
      summary$interpretation[i],
      "\n",
      sep = ""
    )
  }
  if (any(summary$step == 8)) {
    cat(
      "Note: lambda is a process parameter; prefer steps 0-6 as the ",
      "headline for a single paper.\n",
      sep = ""
    )
  }
  invisible(x)
}

# Internal: one-row wide version of a summary: points as `quantity`,
# intervals as `quantity_lower` / `quantity_upper`.
grim_aggregate_wide <- function(summary) {
  out <- list()
  interval <- grepl("^(theta|lambda_support_|lambda_ci_)", summary$quantity)
  for (i in seq_len(nrow(summary))) {
    if (interval[i]) {
      out[[paste0(summary$quantity[i], "_lower")]] <- summary$lower[i]
      out[[paste0(summary$quantity[i], "_upper")]] <- summary$upper[i]
    } else {
      out[[summary$quantity[i]]] <- summary$value[i]
    }
  }
  tibble::as_tibble(out)
}

#' GRIM aggregation by group
#'
#' Runs [grim_aggregate()] separately within each group of `data` (e.g. one
#' group per paper) and returns one row per group. A missing `by` value forms
#' its own group.
#'
#' @param data A data frame as for [grim_aggregate()]. Give per-row
#'   `digits_x` and `items_averaged` as columns.
#' @param by Character vector of grouping column names. They must not clash
#'   with the output columns (such as `k`, `F`, `E` or `O_E`).
#' @param ... Further arguments passed to [grim_aggregate()].
#' @return A wide tibble with the grouping columns followed by one column per
#'   point quantity and `_lower` / `_upper` columns per interval. Groups with
#'   no informative values have `NA` for everything past `F`.
#' @seealso [grim_aggregate()]
#' @examples
#' dat <- data.frame(
#'   paper = rep(c("A", "B"), each = 4),
#'   x = c(3.41, 2.87, 4.12, 3.05, 2.60, 3.93, 1.14, 2.29),
#'   n = c(28, 31, 45, 22, 37, 19, 26, 33)
#' )
#' grim_aggregate_by(dat, by = "paper", digits_x = 2, items_averaged = 1)
#' @export
grim_aggregate_by <- function(data, by, ...) {
  if (!is.data.frame(data) || !is.character(by) || !all(by %in% names(data))) {
    cli::cli_abort("{.arg by} must name columns of {.arg data}.")
  }
  # sort by the grouping columns, NA last, and keep NA as its own group
  ordering <- do.call(order, unname(as.list(data[by])))
  key_parts <- purrr::map(data[by], function(column) {
    paste0(is.na(column), column)
  })
  key <- do.call(paste, c(key_parts, sep = "\r"))[ordering]
  # unnamed, so that purrr's error context names the group's index, not its key
  groups <- unname(split(ordering, factor(key, levels = unique(key))))
  # bind_rows() fills the later columns that groups with k = 0 lack with NA
  groups |>
    purrr::map(function(rows) {
      summary <- grim_aggregate(data[rows, , drop = FALSE], ...)$summary
      wide <- grim_aggregate_wide(summary)
      clash <- intersect(by, names(wide))
      if (length(clash)) {
        cli::cli_abort(
          "{.arg by} column(s) {.arg {clash}} clash with output columns of the
          same name. Rename them.",
          call = NULL
        )
      }
      dplyr::bind_cols(tibble::as_tibble(data[rows[1], by, drop = FALSE]), wide)
    }) |>
    dplyr::bind_rows()
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
#'   x = c(3.41, 2.87, 4.12, 3.05, 2.66, 3.90),
#'   n = c(28, 31, 45, 22, 37, 19)
#' )
#' grim_sensitivity(dat, list(
#'   avoid_0_5 = function(d) ifelse(d %in% c(0, 5), 0.5, 1),
#'   heap_0_5  = function(d) ifelse(d %in% c(0, 5), 3, 1)
#' ), digits_x = 2, items_averaged = 1)
#' @export
grim_sensitivity <- function(data, g_list, ...) {
  if (!is.list(g_list) || is.null(names(g_list)) || any(names(g_list) == "")) {
    cli::cli_abort("{.arg g_list} must be a named list of counterfactuals.")
  }
  if (!"uniform" %in% names(g_list)) {
    g_list <- c(list(uniform = NULL), g_list)
  }
  keep <- c("k", "F", "E", "O_E", "p_fewer_than_random", "p_more_than_random")
  names(g_list) |>
    purrr::map(function(name) {
      summary <- grim_aggregate(data, g = g_list[[name]], ...)$summary
      values <- stats::setNames(
        summary$value[match(keep, summary$quantity)],
        keep
      )
      tibble::as_tibble(c(list(counterfactual = name), as.list(values)))
    }) |>
    dplyr::bind_rows()
}
