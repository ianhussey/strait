# An analytic characterization of integer-achievable (mean, SD) tuples

*Why BRIMMER-consistent tuples without integer solutions exist, exactly where they live, and how to certify
achievability/non-achievability analytically — without CLOSURE-style enumeration, Bernoulli/exhaustive
search, or the brimmest() state-space DP.*

Working notes, 2026-08-25. All claims below are validated empirically; validation scripts and
counts are in this folder (`certifier.py`, `phantom_map.py`, `validate.R`).

---

## 1. Setup and notation

Shift to the integer scale `y = x − l`, so responses are integers `y_i ∈ [0, W]`, `W = u − l`,
sample size `n`. A reported (mean, sd) pins down

- `S = Σ y_i` — one (or a few) integer sums, via the mean's rounding interval (GRIM),
- a window `[Q_lo, Q_hi]` for `Q = Σ y_i²`, via the SD's rounding interval and
  `Q = (n−1)·sd² + S²/n`.

The certification question is therefore purely arithmetic:

> **For which triples (S, Q) does there exist a multiset of n integers in [0, W] with
> Σy = S and Σy² = Q?**

Call the set of achievable Q for a given (S, W, n) the **spectrum** `A(S; W, n)`.
Everything the GRIM family does is a necessary condition on membership of the reported
window in this spectrum. A tuple is a *phantom* when it passes every closed-form screen
(GRIM, SD bounds, GRIMMER/parity) yet its Q-window contains no member of the spectrum.

There is a reflection symmetry `y → W − y` mapping `(S, Q) → (nW − S, Q + nW² − 2WS)`,
so WLOG `S ≤ nW/2` (left half of the umbrella).

## 2. The three closed-form necessary conditions (what BRIMMER already knows)

1. **Parity.** `y² ≡ y (mod 2)`, hence `Q ≡ S (mod 2)` for *every* sample. This is the
   arithmetic core of GRIMMER.
2. **Lower bound** (clustered configuration): with `S = mn + r`,
   `Q_min = (n−r)m² + r(m+1)²`.
3. **Upper bound** (Structure-S configuration): with `S = aW + b`, `0 ≤ b < W`,
   `Q_max = aW² + b²` — `a` values at the ceiling, one remainder, rest at the floor.

These are necessary, not sufficient. The whole question is: **which parity-consistent
integers between Q_min and Q_max are missing from the spectrum?**

## 3. A useful reformulation: triangular numbers and Gauss

Let `T = (Q − S)/2`. Then

```
T = Σ y_i(y_i − 1)/2 = Σ C(y_i, 2)
```

— a sum of *triangular numbers*, one per response, where values 0 and 1 contribute
nothing. Since `8·C(p,2) + 1 = (2p−1)²`, a representation of T by k parts `p_i ≥ 2` is
the same as writing `8T + k` as a sum of k odd squares `(2p_i − 1)²`.

Achievability of (S, Q) is then *exactly*: T can be written as `Σ C(p_i, 2)` with parts
`p_i ∈ [2, W]`, subject to two size/count constraints — `Σ p_i ≤ S` (the leftover
`S − Σp_i` is filled with 1s) and `Σ p_i − k ≥ S − n` (there must be enough response
slots for those 1s).

By Gauss's Eureka theorem, **every** nonnegative integer is a sum of three triangular
numbers (`8T + 3 ≡ 3 (mod 8)` is always a sum of three odd squares). So *without* the
size constraints, every parity-consistent Q would be achievable. **Phantoms exist purely
because, near the top of the umbrella, the size constraint `Σ p_i ≤ S` forces the
representation to use one huge triangular number plus very little else — and huge
triangular numbers are sparse.** That is the entire phenomenon in one sentence.

## 4. Structure theorem for the spectrum

Empirically (exhaustive computation over all designs `W ≤ 8, n ≤ 12`, all S and Q —
859,672 membership queries — plus spot checks at W = 63) the spectrum has this exact
structure:

**(i) Gaps are confined to a top zone.** There is a *dense ceiling* `Q_d(S, W, n)` such
that every parity-consistent Q in `[Q_min, Q_d]` is achievable, and all missing values
lie in `(Q_d, Q_max)`. There are no gaps anywhere near the bottom or middle of the
umbrella (for n ≥ 3; n ≤ 2 is a trivial closed form). Intuition: wherever a
configuration has two equal parts with headroom, the move `(x, x) → (x+1, x−1)` changes
Q by exactly +2, so the spectrum climbs in steps of 2 until configurations are forced
to be near-extremal.

**(ii) The top zone is self-similar ("bands").** Writing a top-zone Q by its largest
part `p₁ = S − j`:

```
Q ∈ A(S; W, n)  ⇔  ∃ j :  Q − p₁² ∈ A(j; min(p₁, W), n−1)
```

For the cap-free corner case (`S ≤ W`), band j occupies
`[(S−j)² + j , (S−j)² + j²]`: two-part configurations `(S−j, j)` sit at band tops,
`(S−j, 1, 1, …)` at band bottoms. Bands for small j are **disjoint** — that is the gap
structure — and they merge into a continuum once `j² + 3j + 4 ≥ 2S + D*(j)`, i.e.
around `j ≈ √(2S)`. The gap set right below `Q_max` is a *translated, shrunken copy of
the whole problem at scale j* — the spectrum is self-similar, with recursion depth
`O(log log S)`.

This band description was verified exactly: over all corner cases `S ≤ 140`, the
recursive characterization reproduced the DP truth with **0 mismatches in 149,540
top-zone queries** (and 0 mismatches on full-range checks for `S ≤ 40`).

**(iii) Depth of the gap zone.** For the corner problem the deepest gap sits at deficit
`D* = Q_max − Q_d ≈ c·S^{3/2}` with c ≈ 3.0–3.4 in the tested range `S ≤ 160`
(band-merging heuristics give `Θ(S^{3/2})` with leading constant ≥ 2√2 ≈ 2.83). For the
capped interior (`S = aW + b`, `a ≥ 2`), the gap zone depth depends only on `(W, b)`
and is `O(W^{3/2})` — shallow — which is why phantoms hug the top edge tightly at
interior means and dig much deeper at the corners.

**(iv) First-gap catalogue (the practically important closed forms).** Achievable
deficits `D = Q_max − Q` near the top are generated by demoting `j` ceiling-parts and
re-partitioning `jW + b` below the cap. The smallest achievable positive deficits, in
closed form (interior, `a ≥ 1`):

| case | smallest positive achievable deficit `D_min` | source |
|---|---|---|
| `b = 0` | `2(W − 1)` | `(W) → (W−1, 1)` |
| `b = 1` | `2(W − 2)` | `(W, 1) → (W−1, 2)` |
| `2 ≤ b ≤ W−2` | `min( 2(b−1), 2(W−b−1) )` | split `b` / demote a ceiling part |
| `b = W−1` | `2(W − 2)` | split `b` |
| corner `a = 0` (i.e. `S < W`) | `2(S − 1)` | `(S) → (S−1, 1)` |

Every parity-consistent Q with `0 < Q_max − Q < D_min` is non-achievable — these
intervals are the single biggest source of phantoms, and the `b ∈ {0, 1}` rows explain
the phantom streaks at (near-)integer means right at the top-centre of the umbrella.
Deeper deficits follow the same generative rule: `D` is achievable iff
`D = jW² + b² − q` for some `j ≥ 0` and `q ∈ A(jW + b; W−1, n−a+j)` — a finite union
of translated corner spectra, each a union of explicit quadratic families
(`D = 2iS' − i² − q'` etc.), evaluable in closed form to any desired depth.
Checked against the exact spectra: e.g. `W=6, S=15` gives `𝔻 ∩ [0,8] = {0, 4, 6}` —
predicting exactly the observed missing values Q ∈ {79, 73} — and similarly across the
full W ≤ 8 sweep.

## 5. Where phantoms live on the umbrella (geography, now exact)

A reported tuple's Q-window has width `ΔQ ≈ 2(n−1)·sd·10^(−d)` (d = reported SD
decimals). A phantom occurs when that window (after parity filtering) fits inside a
spectral gap. Hence:

- **Phantoms ring the entire top edge**, not just the corners. At d = 2 they form an
  almost continuous band hugging `sd_max(mean)` all the way around — including the
  centre. Computed examples (all confirmed impossible by `brimmest()`):
  - 1–7, n = 25, d = 2: mean 4.00, sd ∈ {2.94, 2.99}; mean 3.00, sd 2.80 — dead centre.
  - 1–7, n = 50, d = 2: mean 4.00, sd ∈ {2.99, 3.01, 3.02}.
- **The "corners only" impression is a resolution artifact.** It is accurate for coarse
  reporting (d = 1) and small n, for two compounding reasons: (a) at the corners the
  relevant scale is `S` itself, so the first gap `2(S−1)` and the `O(S^{3/2})` gap-zone
  depth are enormous relative to the scale, letting phantoms sit several reporting
  steps below the top edge; (b) at interior means the gap zone is only `O(W^{3/2})`
  deep and the window `ΔQ` is at its widest (sd is maximal), so only the one or two
  reporting steps immediately under `sd_max` can be phantoms — and at d = 1 those
  steps usually straddle an achievable value. At d = 2 they frequently don't.
- **Streak structure at special means.** At `b = 0` and `b = 1` (integer and
  near-integer means) the first gap `2(W−1)` / `2(W−2)` is maximal for the interior,
  producing the visible phantom streaks at integer means near the top-centre.
- The bottom and middle of the umbrella can never host phantoms (structure theorem (i)):
  any GRIMMER-consistent window there contains an achievable Q.

`phantom_map.py` reproduces the full maps; counts for example designs:
1–5/n9/d1: 26 phantom tuples; 1–7/n9/d1: 50; 1–7/n25/d2: 392; 0–10/n15/d2: 1218;
1–7/n50/d2: 591.

## 6. Is there a *fully* analytic certificate?

Three honest answers, in increasing strength of claim:

1. **A single closed-form formula (GRIM-style) for all cases: no, and one shouldn't be
   expected.** Membership deep in the top zone encodes conditions of the form "is
   `8T + 2` a sum of two odd squares with constrained roots" — Fermat/Jacobi-type
   conditions that depend on the factorization of the quantity, not on any smooth
   formula in (S, Q, W, n). (Gauss rescues k = 3, but the size constraint near the top
   edge forces k ≤ 2 on the large parts — exactly where two-square arithmetic bites.)
   The numbers involved are tiny (≤ ~2S), so this is no computational obstacle, but it
   rules out a one-line algebraic test.

2. **Closed-form gap catalogue for the practically reachable depth: yes.** For any
   reporting precision, only gaps wider than `ΔQ ≈ 2(n−1)·sd·10^(−d)` can produce
   phantoms. These live at bounded deficit, where the catalogue in §4(iv) is an explicit
   finite union of quadratic sequences — evaluable in O(1) per tuple. This yields an
   *analytic phantom test* in exactly the sense of Allard's A-GRIMMER: piecewise closed
   form, no search. For most designs, enumerating the catalogue to deficit
   `D ≤ ΔQ + D_min`-scale settles every tuple the closed-form screens leave open.

3. **An exact, provably terminating analytic certificate for all cases: yes — the band
   recursion.** Membership `Q ∈ A(S; W, n)` is decided by recursion on the largest part
   with the parity/sandwich prune and memoization (`certifier.py: member()`). Because
   gaps exist only in the self-similar top zone, the recursion never explores more than
   the band structure: worst case observed across all stress tests is ~26k integer
   operations (scale 0–200, n = 500, adversarial top-zone impossible queries — a design
   whose brimmest() lattice would need ~10¹⁰ cells), typical cases are O(n) or O(1).
   It certifies impossibility *without* exhausting any search space in the CLOSURE
   sense: what it exhausts is the (provably small) band tree.

Validation of route 3 (this is the headline number): **0 mismatches in 859,672
exhaustive membership checks** (all `W ≤ 8`, `2 ≤ n ≤ 12`, all S, all Q), 0 mismatches
in 10,210 randomized checks at `W = 63, n = 6`, and — at the reported-tuple level —
**0 disagreements with `brimmest()` across 17,580 grid tuples** (full reporting grids
for 1–5/n9/d1, 1–7/n9/d1, 1–7/n11/d1 and a subsampled 0–4/n15/d2 grid; `validate.R`),
including the canonical blind-spot tuple 1–5, n=9, mean 1.3, sd 0.9.

## 7. Relation to brimmest() and integration notes

`brimmest()`'s layer-2 constructive search (`R/certify-sandwich.R`) is already
recognizably the band recursion — largest-part-descending with the sandwich re-applied
per node. The analytic results above upgrade it from heuristic-with-budget to
certified-fast:

- **Memoize failed states** `(S_left, Q_left, cap, n_left)`. The current search
  re-explores equivalent subtrees; with memoization the search *is* the band recursion,
  and the structure theorem guarantees it terminates after the band tree — so the
  `search_budget`/NA path and the corridor-DP fallback become unreachable in practice.
  This removes the last non-analytic component from single-tuple certification.
- **O(1) phantom screen first**: apply the §4(iv) first-gap catalogue before any
  search. For windows entirely inside `(Q_max − D_min, Q_max)` (once parity-filtered),
  impossibility is settled in closed form; for windows reaching `Q ≤ Q_max − D*`-ish
  territory, possibility is settled in closed form. Only windows overlapping the thin
  band structure need the recursion at all.
- The umbrella-plot phantom overlay can be *predicted* rather than enumerated: the
  quadratic families of §4(iv) drawn directly in (mean, sd) coordinates.

Prototypes in this folder: `certifier.py` (reference implementation + exhaustive
validator), `phantom_map.py` (tuple-level phantom maps), `certify_analytic.R`
(R translation + agreement test against `brimmest()`).

## 8. Benchmark vs brimmest() (2026-08-25)

Two analytic implementations were benchmarked (`bench.R`, `bench2.R`, `bench3.R`):
v1 = per-Q membership recursion (`certify_analytic`), v2 = window-level deficit
decomposition (`certify_analytic2`, the recommended form — one `.hit()` call decides a
whole Q-window at once via the demotion decomposition).

**Correctness.** Zero disagreements with `brimmest()` across every comparison run —
~35,000 reported-tuple verdicts over 10 designs spanning 1–5/n=9/d=1 to 0–63/n=200/d=2
and 1–7/n=2000/d=2, edge-hugging (phantom-rich) and interior tuples alike. Since the
two methods share no code or algorithmic idea, this is also an independent validation
that brimmest() is producing correct verdicts on those designs. (`window.py`
additionally validates the window-level algorithm against brute-force DP:
180,180 checks, 0 mismatches.) In 360 probed edge tuples on wide/large designs,
brimmest()'s corridor-DP fallback was never triggered — its layer-2 search settled
everything within budget — so the corridor path itself was not stress-exercised here.

**Speed** (per-tuple, R prototype vs brimmest, half_up∪half_down):

| design | brimmest med / max | analytic v2 med / max |
|---|---|---|
| 0–63, n=50, d=2 (edge tuples) | 0.5 / 50 ms | 0.08 / 4.4 ms |
| 0–63, n=200, d=2 (edge) | 2.7 / 10 ms | 0.09 / 0.2 ms |
| 1–7, n=2000, d=2 (edge) | 19.8 / 65 ms | 0.13 / 0.3 ms |
| small-design full grids | ~equal (lattice cache amortizes) | ~equal |

v2 is 6–150× faster on the hard designs with far tighter tails, and its cost is
essentially independent of n and of the lattice size (no cache, no per-design state).
On tiny designs brimmest's cached lattice is already effectively free, so there is
nothing to win. The v1 prototype showed that naive per-Q membership has bad worst
cases in R (a 7.9 s outlier on a wide window at 0–63/n=200 that v2 settles in 0.1 ms);
window-level certification is the right formulation.

### 8b. Wide-scale / awkward-n verification (sum-score designs)

Because the earlier sweeps used narrow scales and round n, a second round targeted
sum-score-style designs — e.g. a 50-item 1–7 scale summed (l=50, u=350, W=300), plus
W=173 and deliberately weird bounds 13–530 (W=517) — with non-round n ∈ {29, 37, 61,
83, 217}, oversampling the regions where phantoms should live (top edge, corners,
integer-mean streaks) and should not (interior, bottom). Results (`wide.py`/`bench4.R`):

- **Ground truth**: brute-force DP at W=300, n=4 — 69,954 three-way checks (DP vs both
  analytic algorithms), 0 mismatches.
- **Cross-validation** of the two algorithmically independent analytic certifiers
  (largest-part recursion vs demotion decomposition) at W ∈ {173, 300, 517},
  n ∈ {29, 37, 61, 83, 217}: 187,918 top-zone-oversampled point and window queries,
  0 mismatches.
- **vs brimmest()**: 1,200 region-stratified tuples across 4 wide designs (d=1 and
  d=2): 0 disagreements, 0 brimmest errors, corridor-DP fallback never needed.
- **Geography confirmed at scale**: interior and bottom controls produced 0 impossible
  verdicts (as the structure theorem requires); corners were phantom-rich at d=2
  (49–52 of 60 tuples impossible); at d=1 with n=83 phantoms all but vanish
  (ΔQ ≈ 2(n−1)·sd·10⁻¹ swamps the gaps) — 1 corner survivor.
- **Speed** (after adding two exactness-preserving fast paths to the R prototype — a
  natural-part-size greedy descent for mid-depth windows and a direct-membership probe
  for bottom-half windows): v2 medians 0.2–3.5 ms across all regions and designs,
  5–15× faster than brimmest on the phantom-critical edge/corner regions,
  comparable elsewhere.

The wide-scale round also earned its keep by exposing two bugs in the *newly added
speed paths* of the R prototype (a greedy candidate-order flaw and a division-by-zero
at the terminal step) — both caused fallthrough/halt, never a wrong verdict; the exact
core (`.hit()`, `member`) has produced zero incorrect verdicts in every test to date.

### 8c. d = 3 and the consolidated validation suite

All validations were consolidated into self-contained scripts
(`validate_01…validate_06`, see `README.md`), and extended to **three reported
decimals** — the sharpest test, since a d = 3 window typically pins at most one or two
integer Q values, stressing the endpoint/parity semantics and the gap structure
directly. Results, all clean:

- `validate_04` (vs brute-force ground truth, d = 1/2/3, full top-edge strips +
  random controls): 147,671 tuples, **0 mismatches**.
- `validate_05` (vs `brimmest()`, small designs, d = 1 exhaustive grids and d = 2/3
  strips): 280,801 tuples, **0 disagreements**.
- `validate_06` (vs `brimmest()`, wide sum-score designs, d = 1/2/3, region-stratified):
  2,100 tuples, **0 disagreements, 0 errors**.

Two d = 3 phenomena worth naming so the region counts are not misread:

1. **Phantom counts barely grow from d = 2 to d = 3** (e.g. 1–7/n=25: 399 → 390;
   1–5/n=9: 40 → 39 in the sampled strips). At d = 3 most extra impossible tuples
   contain *no parity-consistent integer Q at all*, so GRIMMER's parity screen already
   catches them; the residual blind spot (pass-everything-yet-impossible) is governed
   by the gap structure, which does not change with reporting precision.
2. **At d = 3, impossible verdicts appear in interior/bottom regions on wide scales**
   (e.g. 50–350, n=37: 27 of 40 bottom tuples impossible). These are empty-parity-window
   failures of the same GRIMMER class, not gap-phantoms — validation 04's phantom
   accounting confirms the dense zone still contains **zero** tuples that pass the
   closed-form screens yet lack integer solutions, exactly as the structure theorem
   requires.

### 8d. Integration into brimmest() (branch `analytic-engine`)

The engine now lives in the package: `R/certify-analytic.R` implements
`.hit_deficit()` / `.analytic_window()` / `.certify_analytic()` consuming
`.target_states()` rows (so all seven rounding rules inherit the existing
endpoint-inclusion semantics), and `brimmest()` gained
`engine = c("analytic", "search", "both")` — default `"analytic"`, with the
original search/corridor machinery retained verbatim as `"search"` and a
tuple-by-tuple cross-check as `"both"` (also checked against the lattice when
the lattice route is taken). `brimmest_multiple()` passes `engine` through.
Tests: `tests/testthat/test-certify-analytic.R`; full suite 2,459 pass.

**The cross-check paid for itself on its first extended run** (`validate_07`):
it exposed a pre-existing bug in the lattice route's forward rounding — an
attainable mean of exactly 113/25 = 4.52 scales to 451.99999999999994, which
`trunc` filed under 4.51, so `brimmest(1, 7, 25, mean = 4.52, sd = 2.95,
digits = 2, rounding = "trunc")` wrongly returned impossible for a report the
sample {7×14, 5, 1×10} produces exactly. Directional rules were affected at
grid boundaries, half rules at half-boundaries. Fixed by `.round_scaled()`
(snap within 1e-9 of a half-grid boundary, round in scaled integer space);
regression-tested. After the fix, validation 07 — 4,000 tuples × 3 rounding
rules (half_up, native, trunc) across 10 design/digit combos incl. d = 3 and
wide sum-score scales — completes with zero engine disagreements.

## 9. Open items

- A fully rigorous proof of the density theorem (§4(i)) with a sharp explicit
  `Q_d` bound. The band-overlap induction gives the `Θ(S^{3/2})` shape and is airtight
  level-by-level; turning it into a clean closed-form constant (`D* ≤ 4S^{3/2}` for all
  S, say) needs a careful treatment of the capped/count-limited cases. All empirical
  evidence (S ≤ 160 exhaustive) is consistent.
- The `n_items > 1` (`quasiinteger`/mean-scored) generalization: the same theory
  applies on the `mg`-refined lattice, but the parity condition generalizes to a
  congruence mod `2·mg²`-type conditions worth deriving explicitly.
- Possible publication angle: the spectrum characterization appears to be new even as
  pure combinatorics (the partition-spectrum literature has counting results, e.g.
  OEIS A069999, but no achievability characterization); framed for forensic
  metascience it is "analytic CLOSURE certification".
