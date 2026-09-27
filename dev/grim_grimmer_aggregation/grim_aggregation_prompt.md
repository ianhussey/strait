# Task: implement GRIM aggregation quantities in R, each tagged with its assumptions

## Context

I'm a forensic meta-scientist. I want R functions that pool the results of multiple GRIM tests from one paper (or one group of values) into a small set of quantities. The design principle is an **assumption ladder**: each quantity is computed and reported alongside the assumptions it rests on, so that the user can stop at whichever level of assumption they are willing to defend. The deterministic and assumption-light quantities are the headline; model-based quantities are secondary.

Use `scrutiny` for the GRIM test itself (`grim_map()`), base R for everything else, and return tibbles. Means must be supplied as **character** so trailing zeros define the reporting precision; error if they are not.

If the working directory is an R package, add the code under `R/` with roxygen docs and tests under `tests/testthat/`. Otherwise create a single `grim_aggregate.R` plus a `test_grim_aggregate.R`.

## Notation

Value i has reported mean xᵢ (character) with dᵢ decimals, sample size Nᵢ, and number of items mᵢ (default 1). The rounding rule defaults to `"up_or_down"` and must be used identically everywhere below (both for observed consistency and for pass rates).

## Assumptions (labels used throughout)

- **A0 (no authentic failures):** correct values never fail GRIM. Holds only if the rounding convention, items and extraction are right. There is no ε in the core model.
- **U (uniform counterfactual):** a non-genuine value is equally likely to be any value on the reporting grid, so it fails with probability qᵢ.
- **I (independence):** values pass or fail independently. Violated by shared wrong N, copied or derived values.
- **M (mixture):** each value is independently non-genuine with probability λ, otherwise correct.
- **B (authentic baseline):** an authentic value fails with probability εᵢ, supplied by the user (e.g. from an empirically calibrated reference corpus).

## Quantities, in order of increasing assumptions

### 0. GRIM result per value — assumptions: none

cᵢ = 1 if some integer S exists with round(S / (Nᵢmᵢ), dᵢ) = xᵢ under the rounding rule, else 0. Failure indicator fᵢ = 1 − cᵢ. Compute with `scrutiny::grim_map()`.

### 1. Pass rate pᵢ and fail rate qᵢ — assumptions: none

pᵢ = (1 / 10^dᵢ) × Σⱼ c(j / 10^dᵢ; Nᵢ, mᵢ, dᵢ), for j = 0, …, 10^dᵢ − 1

qᵢ = 1 − pᵢ

This is the share of possible decimal endings consistent with Nᵢ. GRIM consistency is periodic with period 1, so enumerating one unit interval is exact. pᵢ is a property of (N, m, d, rounding), not of the reported value. Cache by unique (d, N, m, rounding). Percentages are out of scope: error clearly if flagged.

**Generalisation (implement as an option):** allow a digit distribution g over the grid instead of uniform: pᵢ(g) = Σⱼ g(j) · c(j; Nᵢ, mᵢ, dᵢ). Uniform is g(j) = 1/10^d. This supports sensitivity analyses over non-uniform counterfactuals (e.g. avoidance of 0 and 5 in the final digit, or heaping on 0 and 5). Accept g either as a vector of length 10^d or as a function returning weights for final digits; document which.

### 2. k and F — assumptions: none

k = #{i : qᵢ > 0}, the number of informative values. Values with qᵢ = 0 cannot fail and carry no information; exclude them from everything below but report how many there were.

A value that fails while qᵢ = 0 is impossible and indicates an extraction or reporting error: warn, flag the row, and exclude it.

F = Σ fᵢ over informative values.

### 3. Partial identification bound on θ — assumptions: A0

θ is the proportion of informative values that are non-genuine (error or fabrication). Every failure is non-genuine; any pass could be either. So F ≤ R ≤ k and

θ ∈ [F/k, 1]

No probability involved. The upper bound is always 1; that is the point, not a bug.

### 4. Expected failures under the random counterfactual — assumptions: A0, U

E = Σ qᵢ

The number of failures expected if every value were uniform-random. It measures how much information the set carries. Always report it next to any ratio below.

### 5. Observed-to-expected ratio (random-reference failure ratio) — assumptions: A0, U

O/E = F / E

Needs no independence (expectations add regardless). Not capped. Interpretation: 0 = no failures; < 1 = fewer failures than random values would produce; 1 = as many as random; > 1 = more than random values can produce, i.e. the uniform model misfits (shared wrong N, wrong pᵢ from rounding or items, extraction errors). With small E the ratio moves in steps of 1/E; do not suppress it, but E must be reported with it.

### 6. Distribution of F under all-random — assumptions: A0, U, I

F ~ PoissonBinomial(q₁, …, qₖ), mean E, variance Σ qᵢ(1 − qᵢ).

Implement the pmf by dynamic programming (no dependency): start with P = (1); for each qᵢ, P ← (P·(1 − qᵢ), 0) + (0, P·qᵢ).

Report two one-sided tail probabilities:

- `p_fewer_than_random` = P(F′ ≤ F): more passes than random values would produce; evidence against "every value random".
- `p_more_than_random` = P(F′ ≥ F): more failures than random values would produce; a misfit diagnostic.

These are statements about a specified mechanism (randomisation-test logic), not about a population.

### 7. Mixture model λ — assumptions: A0, U, I, M

P(fᵢ = 1) = λqᵢ

ℓ(λ) = Σ [fᵢ log(λqᵢ) + (1 − fᵢ) log(1 − λqᵢ)]

- Moment estimator λ̃ = F / E = O/E: unbiased under M, unbounded.
- MLE λ̂ on [0, 1]: λ̂ = 0 if F = 0; λ̂ = 1 if F = k; otherwise maximise ℓ numerically and check the boundary at 1. When all qᵢ are equal, λ̂ = min(1, O/E). Note in the docs that the cap causes downward bias when k·q is small.

### 8. Intervals for λ — assumptions: A0, U, I, M

- **Support interval (default):** {λ : ℓ(λ) ≥ ℓ(λ̂) − log c}, with c = 8 and c = 32 (Royall). Find endpoints with `uniroot()`, clamped to [0, 1]. When F = 0 the lower end is 0.
- **Confidence interval (optional, off by default):** invert the Poisson-binomial tails with rates λqᵢ. Lower bound solves P(F′ ≥ F | λ) = α/2 (0 if F = 0; 1 if even λ = 1 is too small); upper bound solves P(F′ ≤ F | λ) = α/2 (1 if P(F′ ≤ F | 1) > α/2). Default conf level 0.90.

Document that λ is a process parameter: appropriate for pooled analyses across papers, not recommended as the headline for a single paper.

### 9. Authentic-baseline comparison (SIR) — assumptions: B (+ I for the tail)

Only computed if the user supplies εᵢ (a vector, or a function of qᵢ / N / d).

- E_A = Σ εᵢ
- SIR = F / E_A (standardized inconsistency ratio)
- `p_vs_baseline` = P(F′ ≥ F) under PoissonBinomial(ε₁, …, εₖ)

Note in the docs that if εᵢ = τqᵢ, SIR = (O/E)/τ.

## Output

Main function, e.g. `grim_aggregate(data, items = 1, rounding = "up_or_down", g = NULL, eps = NULL, support_c = c(8, 32), ci = FALSE, conf_level = 0.90)`, returning a list with:

1. `summary`: a **long** tibble with one row per quantity and columns `quantity`, `value` (or `lower`/`upper` for intervals), `assumptions` (e.g. `"A0, U"`), and `interpretation` (one short phrase). Order the rows by the assumption ladder above.
2. `data`: the input with added columns `consistency`, `f`, `p_pass`, `q`, `informative`, `flag_impossible`.

Also provide:

- `grim_aggregate_by(data, by, ...)`: one summary per group, returned as a wide tibble (one row per group).
- `grim_sensitivity(data, g_list, ...)`: runs `grim_aggregate()` over a named list of digit distributions (including uniform) and returns E, O/E and the tail probabilities for each, so the user can see whether conclusions hold across counterfactuals.
- A `print` method that shows the ladder compactly.

## Tests (testthat)

Check these known values (2 decimals, items = 1, `up_or_down`):

- p_pass at N = 30 is 0.30; at N = 50 is 0.50.
- p_pass at N = 40 and N = 60 is 0.60, because S/N lands on .x25/.x75 ties that both round directions accept. If scrutiny disagrees, report the discrepancy rather than silently changing the expectation.
- p_pass = 1 for N ≥ 100 at 2 decimals.

Also test:

- F = 0 gives θ lower bound 0, O/E 0, λ̂ 0, support interval lower end 0.
- All informative values failing gives λ̂ = 1.
- A failing value with q = 0 is flagged, warned about and excluded.
- The Poisson-binomial pmf sums to 1 and matches `dbinom()` when all qᵢ are equal.
- Simulation check: with all values drawn uniformly on the grid (true λ = 1), the mean of O/E over many replicates is ≈ 1; with true λ = 0.5, ≈ 0.5.
- Non-character `x` errors.

## Things not to do

- Do not introduce ε into the core model (sections 0–8).
- Do not cap O/E.
- Do not label per-paper λ intervals as the main result.
- Do not assume independence anywhere it is not listed.

Before writing code, check the current signatures of `scrutiny::grim_map()` and `scrutiny::decimal_places()` in the installed version, and tell me if anything differs from how they are used here.
