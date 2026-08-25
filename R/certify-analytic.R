# Analytic certification for brimmest(): decide one (S, Q-window) target by
# the structure of the attainable spectrum, with no state-space sweep and no
# node budget.
#
# The set A(S; W, n) = { sum(y^2) : y_1..y_n integers in [0, W], sum(y) = S }
# has a specific shape (derived and validated in dev/analytic_certification/):
# all integers of S's parity between the clustered minimum and the Structure-S
# maximum are attainable EXCEPT for gaps confined to a thin self-similar zone
# below the maximum. Membership of a whole Q-window is decided exactly by the
# *demotion decomposition*: every configuration has some number t = a - j of
# values at the scale ceiling, with the remaining j*W + b units re-partitioned
# below it. That case split is complete, so recursing on it enumerates exactly
# the attainable deficits from the maximum -- cheap where the gaps live
# (the thin top zone), and short-circuited elsewhere.
#
# Exactness of the verdict does not rest on the structure theorem; the
# decomposition alone guarantees it. The theorem is why the recursion is fast:
# impossible windows only exist in the top zone, where the tree is small.
#
# The public entry point is .certify_analytic(), a drop-in sibling of
# .certify_fast() consuming the same .target_states() rows, hence inheriting
# the endpoint-inclusion semantics of every rounding rule. Unlike
# .certify_fast() it never returns NA: there is no budget to exhaust and
# therefore no corridor-DP fallback.
#
# The two engines are cross-checked against each other in
# tests/testthat/test-certify-analytic.R, and brimmest(engine = "both") does
# the same on live data, erroring on any disagreement.

# Internal: Structure-S maximum with the parts-count check: the largest
# sum of squares of n integers in [0, cap] summing to S, or NA when S exceeds
# n * cap (no such sample exists). The count-free .q_max_int() is not enough
# inside the recursion, where the remaining slots genuinely run out.
.q_max_cnt <- function(S, cap, n) {
  if (cap <= 0) return(if (S == 0) 0 else NA_real_)
  a <- S %/% cap
  b <- S - a * cap
  if (a > n || (a == n && b > 0)) return(NA_real_)
  a * cap * cap + b * b
}

# Internal: is any deficit in [dlo, dhi] attainable? Deficits are measured
# down from .q_max_cnt(S, cap, n). Exact by completeness of the demotion
# decomposition: j demoted ceiling-values leave S_j = j*cap + b units to
# place below the ceiling among n_j = n - (a - j) slots, contributing an
# offset deficit d_off plus, recursively, a deficit of the capped subproblem.
.hit_deficit <- function(S, cap, n, dlo, dhi) {
  if (S == 0) return(dlo <= 0 && 0 <= dhi)
  if (cap <= 0 || n <= 0 || dhi < 0) return(FALSE)
  if (cap > S) cap <- S
  Qx <- .q_max_cnt(S, cap, n)
  if (is.na(Qx)) return(FALSE)
  if (dlo > Qx - .q_min_int(S, n)) return(FALSE)
  if (dlo <= 0 && 0 <= dhi) return(TRUE)
  a <- S %/% cap
  b <- S - a * cap
  for (j in 0:a) {
    S_j <- j * cap + b
    n_j <- n - (a - j)
    if (n_j <= 0) next
    qx <- .q_max_cnt(S_j, cap - 1, n_j)
    if (is.na(qx)) next
    d_off <- j * cap * cap + b * b - qx
    if (d_off > dhi) next
    if (.hit_deficit(S_j, cap - 1, n_j, dlo - d_off, dhi - d_off)) return(TRUE)
  }
  FALSE
}

# Internal: exact verdict for one integer Q-window [k_lo, k_hi] (as produced
# by .k_window(): sandwich-clipped, parity-aligned). Reflects to the lower
# half of the scale, then decides.
#
# Windows reaching deep below the maximum are settled TRUE almost immediately
# by the constructive search (the dense zone is where its greedy first descent
# shines); the demotion recursion is the decider everywhere else, and the
# final authority when the search's small budget runs out. The budget only
# affects which code settles a tuple, never the verdict.
.analytic_window <- function(W, n, S, k_lo, k_hi, witness_budget = 5000) {
  if (2 * S > n * W) {
    shift <- n * W * W - 2 * W * S
    S <- n * W - S
    k_lo <- k_lo + shift
    k_hi <- k_hi + shift
  }
  Qx <- .q_max_cnt(S, W, n)
  if (is.na(Qx)) return(FALSE)
  Qn <- .q_min_int(S, n)
  k_lo <- max(k_lo, Qn)
  k_hi <- min(k_hi, Qx)
  if ((k_lo - S) %% 2 != 0) k_lo <- k_lo + 1
  if ((k_hi - S) %% 2 != 0) k_hi <- k_hi - 1
  if (k_lo > k_hi) return(FALSE)
  if (k_hi == Qx || k_lo == Qn) return(TRUE)
  if (k_lo - Qn < Qx - k_hi) {
    got <- .witness_search(W, n, S, k_lo, k_hi, budget = witness_budget)
    if (isTRUE(got$possible)) return(TRUE)
    if (isFALSE(got$possible)) return(FALSE)  # tree exhausted: a proof
  }
  .hit_deficit(S, W, n, Qx - k_hi, Qx - k_lo)
}

# Internal: the analytic engine for one reported tuple under one rounding
# rule. Same contract as .certify_fast() -- takes .target_states() rows --
# but total: always TRUE or FALSE, never NA.
.certify_analytic <- function(W, n, tg) {
  if (is.null(tg) || !nrow(tg)) return(FALSE)
  for (i in seq_len(nrow(tg))) {
    kw <- .k_window(tg$S[i], tg$Q_lo[i], tg$Q_hi[i],
                    tg$lo_incl[i], tg$hi_incl[i], n, W)
    if (is.null(kw)) next
    if (.analytic_window(W, n, tg$S[i], kw[1], kw[2])) return(TRUE)
  }
  FALSE
}
