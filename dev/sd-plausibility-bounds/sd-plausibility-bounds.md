# Task: analytic plausibility bounds for the sample SD of bounded questionnaire scores

## Context

I do forensic meta-science. A common check is whether a reported (M, SD) pair for a bounded 
scale (e.g. BDI-II, range 0–63) is plausible given the sample size N. Hard bounds exist 
(Bhatia–Davis upper bound, integer-scale lower bound), but I want *plausibility* bounds: 
the range of sample SDs that would be observed at a given M and N if the data came from a 
realistic data-generating process.

I previously did this by simulation: fit beta-binomial regressions to a corpus of raw BDI-II 
data to get predicted μ and dispersion at each mean, sampled from that beta-binomial at 
various N, and took the 1% / 0.1% quantiles of the sample SD across iterations. That works 
but is slow and has Monte Carlo error in the tails. I want an analytic implementation, 
validated against simulation.

## Model

Scores X ∈ {0, ..., m} where m = U − L (63 for BDI-II). X ~ BetaBinomial(m, α, β), 
parameterised as μ = α/(α+β) ∈ (0,1) and ρ = 1/(α+β+1) ∈ [0,1), or equivalently 
φ = α+β = 1/ρ − 1.

- E[X] = mμ
- σ² = mμ(1−μ)[1 + (m−1)ρ]
- ρ = 0 is binomial (floor of the family); ρ → 1 puts all mass at the endpoints

Population moments should be computed exactly from the pmf (64 terms), not from closed 
forms: log pmf = lchoose(m,x) + lbeta(x+α, m−x+β) − lbeta(α,β).

## Method of moments from summary statistics

Given (L, U, M, SD): μ̂ = (M−L)/m, ρ̂ = [SD² / (mμ̂(1−μ̂)) − 1] / (m−1).
Work on the ρ scale: ρ̂ is linear in S², so it's unbiased and its sampling variance is 
Var(S²) / [m(m−1)μ(1−μ)]². Do not truncate negative ρ̂.

## Analytic sampling distribution of S²

For any iid DGP with finite fourth moment:
- E[S²] = σ²
- Var(S²) = σ⁴ [2/(N−1) + γ₂/N], where γ₂ = μ₄/σ⁴ − 3 is the excess kurtosis

Two-moment approximation (Box/Satterthwaite): S² ≈ σ² χ²_ν/ν with 
ν = 2 / [2/(N−1) + γ₂/N]. Normal data gives ν = N−1 as a check.
Then q_p(S) = σ √(χ²_{ν,p} / ν).

Third cumulant of S², for a Cornish–Fisher or saddlepoint refinement of the extreme tails 
(in terms of population cumulants κ₂..κ₆):
κ₃(S²) = κ₆/N² + 12κ₂κ₄/(N(N−1)) + 4(N−2)κ₃²/(N(N−1)²) + 8κ₂³/(N−1)²
Check: normal case gives 8σ⁶/(N−1)², matching χ²_{N−1}.

## Heterogeneity

If ρ comes from a corpus meta-regression (ρ̂_i on μ̂_i, weights from Var(ρ̂_i), via 
metafor::rma), there is a between-study SD τ on the ρ scale. The plausibility bound should 
be the quantile of the *mixture* over ρ ~ N(ρ, τ²) truncated to [0,1): compute the mixture 
CDF of S as the average of the per-ρ CDFs, then invert. Do not average per-ρ quantiles.

## Existing R implementation (starting point)

```r
sd_bounds <- function(L, U, M, N, p = 0.01, rho = NULL, tau = 0, sd = NULL, K = 41) {
  m  <- U - L
  mu <- (M - L) / m
  stopifnot(mu > 0, mu < 1, N > 1)
  hard_max <- sqrt(N / (N - 1) * (M - L) * (U - M))   # Bhatia–Davis, sample-corrected
  hard_min <- 0                                        # continuous case; integer data is tighter
  if (is.null(rho)) {
    if (is.null(sd)) stop("supply rho (corpus prediction) or sd (MoM, circular)")
    rho <- (sd^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
    tau <- 0
  }
  bb_moments <- function(rho) {
    x <- 0:m
    if (rho <= 1e-10) {
      pr <- dbinom(x, m, mu)
    } else {
      phi <- 1 / rho - 1; a <- mu * phi; b <- (1 - mu) * phi
      pr <- exp(lchoose(m, x) + lbeta(x + a, m - x + b) - lbeta(a, b))
    }
    mean <- sum(pr * x)
    s2 <- sum(pr * (x - mean)^2)
    m4 <- sum(pr * (x - mean)^4)
    list(s2 = s2, g2 = m4 / s2^2 - 3)
  }
  eff_df <- function(g2) 2 / (2 / (N - 1) + g2 / N)
  if (tau > 0) {
    q <- (seq_len(K) - 0.5) / K
    rhos <- pmin(pmax(qnorm(q, rho, tau), 0), 1 - 1e-6)
  } else {
    rhos <- rho
  }
  mom <- lapply(rhos, bb_moments)
  s2  <- sapply(mom, `[[`, "s2")
  nu  <- sapply(mom, function(z) eff_df(z$g2))
  cdf <- function(s) mean(pchisq(s^2 * nu / s2, nu))
  qsd <- function(pp) uniroot(function(s) cdf(s) - pp, c(1e-6, hard_max * 1.5))$root
  data.frame(mu = mu, rho = rho, tau = tau, N = N, p = p,
             hard_min = hard_min, hard_max = hard_max,
             model_lower = qsd(p), model_upper = qsd(1 - p))
}
```

Reference values it should reproduce (L=0, U=63, M=25, SD=8 → μ=0.397, ρ=0.052, γ₂=−0.242):

| N    | 1% lower | 1% upper | 0.1% lower | 0.1% upper |
| ---- | -------- | -------- | ---------- | ---------- |
| 20   | 5.24     | 10.86    | 4.47       | 11.90      |
| 50   | 6.26     | 9.78     | 5.74       | 10.40      |
| 100  | 6.77     | 9.25     | 6.39       | 9.68       |
| 500  | 7.45     | 8.55     | 7.27       | 8.74       |

## What I want you to do

1. Review and make suggestions on this plan, and the underlying logic of it.

2. Turn this into a small, well-documented R file with 
   testthat tests. Handle edge cases: μ near 0 or 1, ρ ≤ 0 or ≥ 1, 
   N = 2, non-integer scale bounds, missing inputs.

3. Validate the two-moment approximation against Monte Carlo (≥ 10^5 replicates) across a 
   grid of μ ∈ {0.1, 0.3, 0.5}, ρ ∈ {0.02, 0.05, 0.15}, N ∈ {10, 20, 50, 200}. Report where 
   the 1% and 0.1% quantiles diverge by more than, say, 2% relative error.

4. Where it diverges (expected: small N, extreme tails), implement and compare:
   a. Cornish–Fisher correction using κ₃(S²) above (all cumulants from the pmf)
   b. Exact distribution of S² for small N via dynamic programming over the joint 
      (ΣX, ΣX²) state space — feasible because X is a bounded integer. Use this as ground 
      truth where N ≤ ~25.

5. Known defect to fix: at very small N the χ² lower quantile can fall below the binomial 
   floor σ_binom = √(mμ(1−μ)), which is impossible under the model. Report when this 
   happens and clip or replace with the exact method.

6. Replace hard_min = 0 with the correct integer-scale minimum SD for a given (M, N, L, U) — 
   the minimum sample variance over integer vectors of length N with mean M. Implement 
   directly (it's a small integer optimisation: values as tight as possible around M).

7. Add a function that takes a corpus data frame with columns (L, U, M, SD, N) — no raw 
   data — computes per-study (μ̂, ρ̂, Var(ρ̂)), fits a weighted meta-regression of ρ̂ on μ̂ 
   (metafor::rma with a natural spline on μ̂), and returns predicted ρ and τ at any μ. 
   Var(ρ̂) depends on γ₂ at (μ̂, ρ̂); do one iteration re-evaluating it at the shrunken 
   corpus prediction.

8. Write a short vignette / README explaining the method with the formulas above and 
   the caveat that summary statistics cannot verify the beta-binomial's kurtosis, on 
   which the bounds mainly depend.

Use base R plus extraDistr (optional, for cross-checking the pmf), metafor, and testthat. 
Keep functions vectorised over N and M where natural. Tidyverse style is fine.
