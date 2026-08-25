# Prototype: analytic certification of reported (mean, sd) tuples, R translation of
# the band-recursion certifier (see analytic-certification.md and certifier.py).
#
# .ach_member(S, Q, cap, n): exact membership of Q in the spectrum
#   { sum(y^2) : y_1..y_n integers in [0, cap], sum(y) = S }.
#
# certify_analytic(l, u, n, mean, sd, digits): possible/impossible verdict for a
# reported tuple under the union of half_up/half_down rounding (closed window),
# matching brimmest()'s default. No lattice, no corridor DP, no search budget:
# memoized band recursion, guaranteed small because spectral gaps exist only in the
# self-similar top zone.
#
# Run validate.R in this folder for the agreement test against brimmest().

.qmin_i <- function(S, n) {
  m <- S %/% n
  r <- S - n * m
  (n - r) * m * m + r * (m + 1) * (m + 1)
}

.qmax_i <- function(S, cap, n) {
  if (cap == 0) return(if (S == 0) 0 else NA_real_)
  a <- S %/% cap
  b <- S - a * cap
  if (a > n || (a == n && b > 0)) return(NA_real_)
  a * cap * cap + b * b
}

.ach_member <- function(S, Q, cap, n, memo = new.env(parent = emptyenv())) {
  if (S == 0) return(Q == 0)
  if (n <= 0 || cap <= 0) return(FALSE)
  if (cap > S) cap <- S
  key <- paste(S, Q, cap, n)
  hit <- memo[[key]]
  if (!is.null(hit)) return(hit)
  lo <- .qmin_i(S, n)
  hi <- .qmax_i(S, cap, n)
  if (is.na(hi) || Q < lo || Q > hi || (Q - S) %% 2 != 0) {
    memo[[key]] <- FALSE
    return(FALSE)
  }
  if (Q == lo || Q == hi) {
    memo[[key]] <- TRUE
    return(TRUE)
  }
  res <- FALSE
  p_hi <- min(cap, floor(sqrt(Q) + 1e-9))
  p_lo <- ceiling(S / n)
  for (p1 in seq(p_hi, p_lo, by = -1)) {
    if (.ach_member(S - p1, Q - p1 * p1, p1, n - 1, memo)) { res <- TRUE; break }
  }
  memo[[key]] <- res
  res
}

.ach_member_top <- function(S, Q, W, n, memo = new.env(parent = emptyenv())) {
  if (2 * S > n * W) {
    Q <- Q + n * W * W - 2 * W * S
    S <- n * W - S
  }
  .ach_member(S, Q, W, n, memo)
}

#' Certify one or more reported (mean, sd) tuples analytically.
#' Verdicts match brimmest(l, u, n, mean, sd, digits) with the default
#' rounding = c("half_up", "half_down") (whose union is a closed interval).
certify_analytic <- function(l, u, n, mean, sd, digits) {
  W <- u - l
  h <- 0.5 * 10^-digits
  tol <- 1e-9
  nn <- max(length(mean), length(sd))
  mean <- rep(mean, length.out = nn)
  sd <- rep(sd, length.out = nn)
  memo <- new.env(parent = emptyenv())   # shared across tuples of one design
  possible <- logical(nn)
  for (i in seq_len(nn)) {
    my <- mean[i] - l
    S_lo <- ceiling((my - h) * n - tol)
    S_hi <- floor((my + h) * n + tol)
    ok <- FALSE
    for (S in seq(max(0, S_lo), min(n * W, S_hi))) {
      if (S_lo > S_hi) break
      q_lo <- (n - 1) * max(0, sd[i] - h)^2 + S * S / n
      q_hi <- (n - 1) * (sd[i] + h)^2 + S * S / n
      Q_lo <- max(ceiling(q_lo - tol), .qmin_i(S, n))
      Q_hi <- min(floor(q_hi + tol), .qmax_i(S, W, n))
      if (is.na(Q_hi) || Q_lo > Q_hi) next
      if ((Q_lo - S) %% 2 != 0) Q_lo <- Q_lo + 1
      if (Q_lo > Q_hi) next
      for (Q in seq(Q_lo, Q_hi, by = 2)) {
        if (.ach_member_top(S, Q, W, n, memo)) { ok <- TRUE; break }
      }
      if (ok) break
    }
    possible[i] <- ok
  }
  data.frame(mean = mean, sd = sd, possible = possible)
}

# ---- Window-level certification (v2) -----------------------------------------
# Decide a whole Q-window at once: does [Qlo, Qhi] contain an achievable Q?
# .hit(S, cap, n, dlo, dhi): exists a partition of S (<= n parts, parts <= cap)
# whose deficit from qmax(S, cap, n) lies in [dlo, dhi]. Demotion decomposition:
# every configuration has some number t = a - j of parts equal to cap and a
# remainder re-partitioned below the cap; complete, hence exact. Impossible
# windows exist only in the shallow top zone (structure theorem), so exhausting
# the tree there is cheap; possible windows hit a witness branch quickly.
.hit <- function(S, cap, n, dlo, dhi) {
  if (S == 0) return(dlo <= 0 && 0 <= dhi)
  if (cap == 0 || n <= 0 || dhi < 0) return(FALSE)
  if (cap > S) cap <- S
  Qx <- .qmax_i(S, cap, n)
  if (is.na(Qx)) return(FALSE)
  if (dlo > Qx - .qmin_i(S, n)) return(FALSE)
  if (dlo <= 0 && 0 <= dhi) return(TRUE)
  a <- S %/% cap
  b <- S - a * cap
  for (j in 0:a) {
    S_j <- j * cap + b
    n_j <- n - (a - j)
    if (n_j <= 0) next
    qx <- .qmax_i(S_j, cap - 1, n_j)
    if (is.na(qx)) next
    d_off <- j * cap * cap + b * b - qx
    if (d_off > dhi) next
    if (.hit(S_j, cap - 1, n_j, dlo - d_off, dhi - d_off)) return(TRUE)
  }
  FALSE
}

.window_possible <- function(S, Qlo, Qhi, W, n) {
  if (2 * S > n * W) {
    shift <- n * W * W - 2 * W * S
    S <- n * W - S
    Qlo <- Qlo + shift; Qhi <- Qhi + shift
  }
  Qx <- .qmax_i(S, W, n)
  if (is.na(Qx)) return(FALSE)
  Qlo <- max(Qlo, .qmin_i(S, n)); Qhi <- min(Qhi, Qx)
  if ((Qlo - S) %% 2 != 0) Qlo <- Qlo + 1
  if ((Qhi - S) %% 2 != 0) Qhi <- Qhi - 1
  if (Qlo > Qhi) return(FALSE)
  .hit(S, W, n, Qx - Qhi, Qx - Qlo)
}

#' v2: window-level analytic certification (same verdicts as certify_analytic,
#' but one call per (S, window) instead of one membership test per candidate Q).
certify_analytic2 <- function(l, u, n, mean, sd, digits) {
  W <- u - l
  h <- 0.5 * 10^-digits
  tol <- 1e-9
  nn <- max(length(mean), length(sd))
  mean <- rep(mean, length.out = nn)
  sd <- rep(sd, length.out = nn)
  possible <- logical(nn)
  for (i in seq_len(nn)) {
    my <- mean[i] - l
    S_lo <- max(0, ceiling((my - h) * n - tol))
    S_hi <- min(n * W, floor((my + h) * n + tol))
    ok <- FALSE
    if (S_lo <= S_hi) for (S in S_lo:S_hi) {
      q_lo <- (n - 1) * max(0, sd[i] - h)^2 + S * S / n
      q_hi <- (n - 1) * (sd[i] + h)^2 + S * S / n
      if (.window_possible(S, ceiling(q_lo - tol), floor(q_hi + tol), W, n)) {
        ok <- TRUE; break
      }
    }
    possible[i] <- ok
  }
  data.frame(mean = mean, sd = sd, possible = possible)
}

# ---- Greedy window witness (fast path for 'possible' verdicts) ---------------
# Pure descent, no backtracking: place the largest part that keeps the remaining
# (sum, Q-window) feasible against the sandwich + parity. Succeeds on virtually
# every achievable window (dense zone); returns FALSE on failure, in which case
# the exact .hit() decides. Exactness is unaffected -- this only short-circuits
# TRUE verdicts.
.greedy_window <- function(S, cap, n, lo, hi) {
  repeat {
    if (S == 0) return(lo <= 0 && 0 <= hi)
    if (n <= 0) return(FALSE)
    if (cap > S) cap <- S
    p_hi <- min(cap, floor(sqrt(hi) + 1e-9))
    p_lo <- ceiling(S / n)
    if (p_hi < p_lo) return(FALSE)
    stepped <- FALSE
    # candidate parts: near the natural size hi/S (near-uniform configurations,
    # right for windows deep below Qmax) and near the cap (right for top-zone
    # windows), largest first within each family
    v <- max(p_lo, min(p_hi, floor(hi / max(S, 1))))
    cand <- unique(pmin(p_hi, pmax(p_lo, c(v + 2, v + 1, v, seq(p_hi, max(p_lo, p_hi - 12), by = -1), v + 4, v + 8))))
    for (p1 in cand) {
      S2 <- S - p1
      if (S2 == 0) {
        if (lo <= p1 * p1 && p1 * p1 <= hi) return(TRUE)
        next
      }
      if (n - 1 <= 0) next
      lo2 <- max(lo - p1 * p1, .qmin_i(S2, n - 1))
      hi2 <- min(hi - p1 * p1, .qmax_i(S2, min(p1, S2), n - 1))
      if (is.na(hi2)) next
      if ((lo2 - S2) %% 2 != 0) lo2 <- lo2 + 1
      if (lo2 > hi2) next
      S <- S2; cap <- p1; n <- n - 1; lo <- lo2; hi <- hi2
      stepped <- TRUE
      break
    }
    if (!stepped) return(FALSE)
  }
}
# use the greedy fast path inside .window_possible
.window_possible <- function(S, Qlo, Qhi, W, n) {
  if (2 * S > n * W) {
    shift <- n * W * W - 2 * W * S
    S <- n * W - S
    Qlo <- Qlo + shift; Qhi <- Qhi + shift
  }
  Qx <- .qmax_i(S, W, n)
  if (is.na(Qx)) return(FALSE)
  Qlo <- max(Qlo, .qmin_i(S, n)); Qhi <- min(Qhi, Qx)
  if ((Qlo - S) %% 2 != 0) Qlo <- Qlo + 1
  if ((Qhi - S) %% 2 != 0) Qhi <- Qhi - 1
  if (Qlo > Qhi) return(FALSE)
  # bottom-half windows sit in the dense zone: the lowest parity-consistent Q is
  # nearly always achievable, and direct membership finds its witness in O(n).
  if (Qlo - .qmin_i(S, n) < Qx - Qhi &&
      .ach_member(S, Qlo, W, n, new.env(parent = emptyenv()))) return(TRUE)
  if (.greedy_window(S, W, n, Qlo, Qhi)) return(TRUE)
  .hit(S, W, n, Qx - Qhi, Qx - Qlo)
}
