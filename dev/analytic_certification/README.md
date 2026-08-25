# Analytic certification: theory, implementations, validation suite

See `analytic-certification.md` for the mathematics (spectrum structure theorem,
gap catalogue, phantom geography) and benchmark narrative.

## Implementations

| file | what |
|---|---|
| `certifier.py` | membership certifier: `member_top(S, Q, W, n)` — largest-part band recursion with memoization; running it directly executes validation 01a |
| `window.py` | window certifier: `window_possible(S, Qlo, Qhi, W, n)` — demotion-decomposition `hit()`; algorithmically independent of `certifier.py`; running it directly executes a self-check |
| `certify_analytic.R` | R prototype: `.ach_member()`, `.hit()`, `.window_possible()`, and the reported-tuple wrappers `certify_analytic()` (v1, per-Q) and `certify_analytic2()` (v2, window-level — the recommended form) |
| `phantom_map.py` | enumerates all phantom tuples of a design and where they sit on the umbrella |

## Validation suite (all expected to print 0 mismatches / 0 disagreements)

| script | validates | against | last run |
|---|---|---|---|
| `validate_01_membership.py` | membership certifier, all designs W≤8 n≤12 (every S, Q) + W=63 randomized | brute-force DP | 869,882 checks, 0 mismatches |
| `validate_02_windows.py` | window certifier, point + random windows, W≤8 n≤12 | brute-force DP | 180,180 checks, 0 mismatches |
| `validate_03_wide.py` | wide scales / awkward n: W=300 n=4 three-way vs DP; cross-validation of the two independent algorithms at W∈{173,300,517}, n∈{29,37,61,83,217} | DP + each other | 257,872 checks, 0 mismatches |
| `validate_04_tuples_dp.py` | reported-tuple layer at d=1, 2, **3**: full top-edge strip + random controls | brute-force DP | 147,671 tuples, 0 mismatches |
| `validate_05_grids.R` | reported-tuple layer at d=1 (exhaustive grids), d=2, **3** (top-edge strip + random), small designs | `brimmest()` | 280,801 tuples, 0 disagreements |
| `validate_06_wide.R` | wide sum-score designs (50–350, 13–530), awkward n, d=1, 2, **3**, region-stratified (edge/corner/interior/bottom) | `brimmest()` | 2,100 tuples, 0 disagreements, 0 errors |
| `validate_07_engine_both.R` | integrated `brimmest(engine = "both")`: analytic vs search/lattice inside the package, rules half_up/native/trunc, d=1/2/3, small + wide designs | itself (errors on disagreement) | 12,000 verdicts, 0 disagreements (after fixing the lattice trunc-boundary bug it exposed) |

Benchmarks (speed, not correctness): `bench.R`, `bench2.R`, `bench3.R`, `bench4.R`.
Superseded by the numbered suite but kept runnable: `validate.R`.

Run everything:

```
python3 validate_01_membership.py
python3 validate_02_windows.py
python3 -u validate_03_wide.py          # ~10-30 min
python3 -u validate_04_tuples_dp.py
Rscript dev/analytic_certification/validate_05_grids.R   # from package root
Rscript dev/analytic_certification/validate_06_wide.R    # from package root
```

## Notes on d = 3

At three reported decimals the sd window pins at most one or two integer Q values,
so most reporting-grid tuples fail GRIMMER's parity screen outright and the
*phantom* count (tuples passing every closed-form screen yet impossible) stays
comparable to d = 2 while total impossibles explode. d = 3 is the sharpest test of
endpoint/parity semantics and of the gap structure — validated clean against both
brute force (04) and `brimmest()` (05, 06).

## Known coverage limits

- Wide scales (W ≥ 100) with 5 ≤ n ≤ 28 are covered by two-algorithm cross-validation
  only (brute-force DP is intractable there beyond n = 4).
- `n_items > 1` (quasi-integer granularity) not yet treated.
