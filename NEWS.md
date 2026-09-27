# strait 0.4.12

## New features

* `grim_aggregate()` pools many GRIM tests from one paper (or one group of
  values) into an assumption ladder: counts and a partial-identification
  bound on the non-genuine share that need no model; expected failures,
  the observed-to-expected ratio and exact Poisson-binomial tail
  probabilities under a uniform-random counterfactual; a mixture-model
  non-genuine rate with likelihood support intervals; and an optional
  comparison against a user-supplied authentic baseline. Every quantity is
  reported with the assumptions it rests on. Pass rates are exact, from
  enumerating the reporting grid with the same scrutiny call used for the
  observed values, so they include the rounding ties the closed-form
  approximation misses (at two decimals, p = 0.60 at N = 40, not 0.40).
  `grim_aggregate_by()` gives one row per group, and `grim_sensitivity()`
  reruns the ladder under non-uniform digit distributions. The logic and
  mathematics are laid out in `paper/grim-aggregation/`.

# strait 0.4.11

## New features

* `brimmest()` gains an `engine` argument. The new default, `"analytic"`,
  certifies single reports against the structure of the attainable
  sum-of-squares spectrum (the demotion decomposition in
  `R/certify-analytic.R`): exact in both directions, with no node budget and
  no corridor-sweep fallback, and typically 5-150x faster than the search on
  wide or large-n designs. `"search"` keeps the original constructive-search /
  corridor-DP machinery unchanged, and `"both"` runs the two independent
  methods on every tuple and errors on any disagreement -- a built-in
  cross-check (against the full lattice too, when the lattice route is
  taken). All seven rounding rules are supported identically by every engine.
  `brimmest_multiple()` accepts `engine` as a call-level constant.

  The underlying characterization of which GRIMMER-consistent tuples lack
  integer solutions, with its validation suite (about 1.5 million checks
  against brute force, the search engine, and the lattice, at 1-3 reported
  decimals), is developed in `dev/analytic_certification/`.

## Bug fixes

* The lattice route declared some reports impossible that real samples
  produce exactly, under the directional rounding rules. An attainable mean
  of exactly 113/25 = 4.52 scales to the double 451.99999999999994, which
  `trunc` files under 4.51, so `brimmest(1, 7, 25, mean = 4.52, sd = 2.95,
  digits = 2, rounding = "trunc")` returned FALSE for a report that the
  sample {7 x 14, 5, 1 x 10} yields exactly. `trunc`, `floor`, `ceiling` and
  `anti_trunc` were affected at grid boundaries, the half rules at
  half-boundaries. The lattice's forward rounding now snaps scaled values to
  the boundary when within 1e-9 of one -- the same tolerance the targeted
  route has always used -- and rounds in scaled integer space. Found by the
  new `engine = "both"` cross-check on its first extended run.

# strait 0.4.10

## Bug fixes

* `plot_sd_bounds(shade = "outside")` drew the feasible region as thousands of
  one-column slivers instead of a single band, for any design with a dense kink
  lattice — roughly `n * (u - l) > 333` at the default `by`. At `l = 1, u = 7,
  n = 60` it produced 921 rings where there is one region; at `n = 200`, 2001.
  The verdicts were never affected, only the figure, but the figure asserted
  that the means between the slivers admitted no SD at all, which is the one
  thing that shading convention must not say.

  The cause was that `plot_sd_bounds()` inferred the mean-grid spacing
  `band_polygon()` needs — to tell a sampling gap from a genuine one — as
  `median(diff(curve$mean))`. `sd_bounds_curve()`'s grid is deliberately not
  uniform: on top of `seq(l, u, by)` it places every kink of the
  `1/(n * n_items)` lattice and a pair of neighbours `1e-9` away, so each kink
  contributes three closely spaced points. Once the kinks outnumber the uniform
  grid the median falls *below* the plain grid spacing, and every ordinary
  interval is then read as a gap.

  `sd_bounds_curve()` now records the spacing it used in a `"step"` attribute
  and `plot_sd_bounds()` reads it, the same way `plot_sd_region()` has always
  passed the step it knows. A curve built by hand, without the attribute, still
  falls back to the median, which is correct on the uniform grid such a curve
  should be on.

# strait 0.4.9

## New features

* `brimmest_multiple()` applies `brimmest()` across a data frame, as
  `brimmer_multiple()` does for the closed-form screen. Columns matching
  `brimmest()` arguments are read per row and anything in `...` is broadcast.
  Rows are grouped by design and each group certified in one call, because one
  attainable lattice serves every report of a design and is built once —
  certifying a 200-row table on a single design takes about 2 ms. Repeated
  `(mean, sd)` pairs within a design are computed once. `rounding`,
  `max_cells` and `search_budget` may be given only as constants: `rounding`
  is a set of admitted rules rather than one value per row, so reading it per
  row would silently mangle it.

* `brimmest()` gains a fast certification path that decides most single
  reports by arithmetic and a constructive search, instead of sweeping a state
  space sized by the design. The verdict is unchanged — the new route is
  verified to agree with the existing ones cell for cell — but the cost is
  now set by the report rather than by the scale width and sample size.

  Two layers do the work. A closed-form sandwich screen bounds the achievable
  sum of squares at a candidate sample sum between the clustered and
  Structure-S configurations, and adds the parity condition `Q = S (mod 2)`
  that holds for every sample because `y^2 = y (mod 2)`. A window meeting no
  admissible integer between those bounds is impossible on arithmetic alone.
  What survives goes to a search over non-increasing samples — the partitions
  of the sum — re-applying the sandwich at every step to what is left to
  place. Completing a sample proves possibility and produces a witness;
  exhausting the search proves impossibility. Reflecting `y -> W - y` when the
  sum sits above the scale midpoint means the search always builds from the
  nearer wall, where the tree is shallow, which is exactly where the
  closed-form screens leak.

  The search runs depth-first on an explicit stack rather than by recursion,
  so sample size is not a constraint on it: a report at `n = 20000` is
  certified, with a witness, in under half a second.

* `brimmest()` gains a `search_budget` argument (default `2e5` nodes) bounding
  that search. Reports it does not settle in budget fall through to the
  corridor dynamic program as before, so the argument trades route for route
  and never changes a verdict.

## Performance

* A single report on a 0-63 inventory at `n = 50` — a design the full lattice
  refuses outright — falls from 78 s to 0.63 ms. A 0-100 scale at `n = 100`,
  previously beyond either sweeping route, certifies in 1.5 ms. The blind-spot
  cell of the documentation (mean 1.3, SD 0.9, `n = 9`, 1-5) falls from 493 to
  41 microseconds, deciding both rounding rules in nine node expansions. Every
  figure here is produced by the benchmark table in
  `validation/validate_sd_bounds_functions.qmd`, medians of repeated runs over
  both rounding rules.

* The search itself costs about 6 microseconds per node, down from 22. Almost
  all of that was R's own overhead rather than arithmetic — the cost was flat
  in scale width, which is the signature. Children are now ordered lazily:
  only the best one is wanted on the way down, and the way down is nearly
  always the whole search, so `which.min()` serves the first visit and
  `order()` is paid only by a frame that is returned to. `order()` inspects
  its arguments through `match.arg()` and two `vapply()` passes before sorting
  anything, which at these lengths costs several times the sort. `ifelse()`,
  `pmax()`/`pmin()` and a general-purpose `.q_max_int()` call were likewise
  replaced with masked arithmetic in the one function that runs per node, and
  two conditions in its filter dropped as provably redundant given the
  candidate range. No verdict changes: the ordering only decides which witness
  is found first.
* Grid-sized workloads are unaffected: `brimmest()` still routes a whole
  reporting grid to the cached lattice, which remains much cheaper per cell
  when the design is being enumerated anyway.

## Validation

* The new path was checked exhaustively against enumeration for every sample
  sum and every narrow window at `n = 2..4` on scales up to `W = 5`, in both
  directions, with every returned witness verified to be a real sample hitting
  its window.
* Against the lattice route it agreed on all 2,719,364 reporting-grid cells
  spanning six designs, two reporting precisions and two rounding rules, and
  was never cut short by its node budget.

# strait 0.4.8

## Breaking changes

* `sd_bounds_check()` is renamed `brimmer()` and `sd_bounds_check_multiple()`
  to `brimmer_multiple()`, naming the SD-side bounds test by analogy with
  GRIM / GRIMMER. No deprecated aliases are kept: the package has not been
  released under the old names.
* `brimmer()` no longer requires a reported `sd`. Supplying only a mean runs
  the mean-side tests and returns `NA` in the SD columns.
* `brimmer()` gains an `in_scale_range` column and no longer reports a
  granularity failure as a bounds failure. Previously a GRIM-impossible mean
  under `Z = "integer"` set `failed_tests` to `"feasibility,grim"`, because
  `sd_bounds()` folds GRIM into `feasible` — correctly, since with no
  GRIM-consistent mean no bounds exist, but that made an out-of-range mean and
  a granular-impossible one indistinguishable. `in_scale_range` is now
  computed independently from `feasible_mean_band()` and the mean's rounding
  interval, so the two defects are named separately, and `feasibility` is
  reserved as a residual for infeasibility no other test accounts for.

## New features

* `band_polygon()` turns a `(mean, lo, hi)` band into closed rings, and the
  outside-shading branches of `plot_sd_bounds()` and `plot_sd_region()` now
  knock those rings out of a shaded panel rather than assembling the shading
  around them. This fixes a real defect: some constraint sets leave stretches
  where no sample exists and `lo`/`hi` are `NA` — 174 of 3001 means for
  `sd_region_data(0, 3, 7, rule = "alpha", n_items = 2, alpha = 0.70)` — and a
  ribbon draws nothing across them. Assembled shading therefore left those
  means unshaded, asserting that any SD at all was possible there. Rings leave
  them shaded by construction, and close themselves, so a region ending
  mid-scale gets a vertical edge instead of two curves trailing off.
* `plot_sd_region()` gains `shade`, defaulting to `"outside"`, matching
  `plot_sd_bounds()`. `shade = "inside"` fills the band as before.
* `sd_delta()` exposes the count-parity correction relating Muilwijk's
  mean-conditional ceiling to the sharp one. The package computed the
  corrected maximum but never the correction, which is interpretable on its
  own as the fraction of Muilwijk's ceiling that is reachable: 6/7 at the
  midpoint of a 1-5 scale at `n = 7`, about 0.11 at a mean of 1.07. It also
  pins two independent implementations against each other, since
  `sd_max_structure_s()^2` must equal `sd_max_muilwijk()^2 * sd_delta()`.
* `plot_umbrella()` gains `style`, defaulting to `"points"`: the panel is
  shaded and only the consistent tuples are drawn. The previous `"tiles"` view
  is retained and is still the right choice for methods exposition, but it was
  the wrong default for reading an answer off the page — at `n = 14` on a 1-7
  scale at two decimal places it drew 19,683 cells to show 3,270 consistent tuples, spent
  two thirds of its ink and its most saturated colour on impossible ones, and
  produced interference banding that hid the vertical striping. The points style also
  accepts an already-filtered `sd_region_data(rule = "integer")` lattice.
* `rule = "muilwijk"` is an alias for `rule = "mean"` in `plot_sd_region()` and
  `sd_region_data()`, since that rule was named after neither its constraint
  set nor its author.


* `brimmest()` now certifies designs the full lattice refuses. A reported
  tuple pins the sum to a few candidate values and the sum of squares to a
  narrow window, so only states that can still reach one of those targets are
  visited. That corridor is roughly a tenth of the full state table, and a hit
  can be declared at any layer because scores at the scale minimum pad a short
  sample out to `n`. A 0-63 inventory at `n = 50` needs a 3151 x 24801 table
  and previously failed with a size error; it now returns a verdict. On large designs where
  both routes work the targeted one is about twice as fast. Small designs
  still use the lattice, which is cheaper there and is now cached per
  `(l, u, n)` so repeated calls pay for it once. Verdicts are unchanged:
  equivalence with lattice membership is asserted over 1.75 million cells
  spanning six designs, two precisions and two rounding rules.

* `plot_sd_bounds()` gains `shade`, defaulting to `"outside"`: the infeasible
  region is shaded and the feasible band left clear, so the native scale now
  reads the same way as `plot_sd_bounds_pomp(reference = "sharp")` — shaded
  means unreachable. `shade = "inside"` restores the previous filled band.
* `plot_sd_bounds()` gains `expand`, the padding around the plotted region as
  a proportion of the scale width `u - l` rather than a fixed number of SD
  units, so the margin looks the same on a 1-5 scale and a 0-100 one. The
  limits always stretch to include `points`, so an out-of-bounds report cannot
  be clipped out of view by the finite limits that outside-shading requires.

* `brimmest()` completes the `brim()` / `brimmer()` / `brimmest()` family with
  an exact possible / impossible certificate for reported `(mean, sd)` tuples
  on a bounded integer scale. The superlative marks the substantive
  difference: `brim()` and `brimmer()`, like GRIM and GRIMMER, test conditions
  that are necessary but not sufficient, so passing them proves nothing;
  `brimmest()` is necessary *and* sufficient, and so can certify a report as
  having an integer solution or as having none. It closes the residual blind
  spot the closed-form screen admits. It enumerates the attainable lattice by
  dynamic programming over the reachable (sum, sum of squares) states and tests
  membership after rounding, so no dataset is reconstructed — the same
  certificate `unsum::closure_generate()` provides, obtained analytically. It
  is vectorised over tuples, since one lattice serves a whole design. Verified
  cell-for-cell against CLOSURE across six designs (~5,700 tuples, zero
  disagreements) at roughly 700-1000x the speed. The certification document in
  `validation/` covers the comparison; `unsum` is now suggested.
* `brim()` is the mean-side bounds test: is the reported mean attainable at
  all, given the scale limits, `n` and any attained extremes? It reports the
  feasible mean band (`band_lo`, `band_hi`) alongside the verdict. Under
  `Z = "integer"` it also runs GRIM (deferred to `scrutiny`) but never
  GRIMMER, so `brimmer()` is nested on `brim()` exactly as GRIMMER is on GRIM.

## Compatibility

* The package now works against both generations of the `scrutiny` GRIM /
  GRIMMER interface. CRAN `scrutiny` (0.6.1 and earlier) takes the reported
  mean and SD as zero-padded strings; later versions take numbers plus their
  decimal places and reject strings outright, so no single call satisfies
  both. `R/scrutiny-compat.R` detects which interface is installed, per
  function, and dispatches accordingly. Without this the package failed on
  CRAN `scrutiny` with `unused argument (digits_x = ...)`.

## Documentation

* Every exported function now carries a running example.
* New validation document `validation/certification.qmd`, comparing CLOSURE
  reconstruction with analytic certification and reporting the measured cost
  of each.
* `Language` is now `en-GB`, matching the prose, with `inst/WORDLIST`
  regenerated; `spelling::spell_check_package()` reports nothing.
* The copyright year is recorded as a range (2024-2026) in `LICENSE` and
  `LICENSE.md`, which had drifted apart in format.

## Known issues

* **CRAN `scrutiny` (0.6.1) returns wrong GRIMMER verdicts.** Its GRIMMER
  test 3 flags attainable values as inconsistent
  (<https://github.com/lhdjung/scrutiny/issues/80>), so `brimmer()`,
  `umbrella_data()` and the lattice rules inherit those false flags when run
  against it — about 2% of grid cells in the designs tested, always in the
  direction of wrongly rejecting a legitimate report. `brimmest()` is
  unaffected, and disproves the false flags constructively. Until `scrutiny`
  0.6.2 reaches CRAN, prefer the development version of `scrutiny`, or read
  `grimmer` verdicts alongside `brimmest()`. GRIMMER on 0.6.1 is also roughly
  25x slower, which dominates check time.

# strait 0.4.4

* The package was renamed from `tides` to `strait`, to avoid a collision with
  the existing CRAN package `Tides`. The exported API was unchanged; only the
  package name, its vignette (`vignette("strait")`) and the repository URL
  differed. Entries below describe the package under its former name, and
  historical function names (`tides()`, `tides_df()`, `plot_tides()`,
  `plot_tides_relative()`) are left as they were written.

# strait 0.4.0

A ground-up rewrite of the engine around the closed-form standard-deviation
bounds derived in the STRAIT article. This is a breaking change: the API is new
throughout.

## Breaking changes

* Removed `tides()`, `tides_df()`, `umbrella()`, `approximate_sd_bounds()`,
  `plot_tides()` and `plot_tides_relative()`, and the internal
  `grimmer_consistent()` scrutiny shim.
* `sd_bounds()` has a new signature and semantics. It now takes logical scale
  limits `l`, `u` (formerly `min`, `max`), optional attained extremes `a`, `b`,
  a granularity argument `Z` (`"continuous"`, `"integer"`, `"quasiinteger"`),
  `scoring` (`"singleitem"`, `"sumscored"`, `"meanscored"`) with `n_items`, and
  an optional Cronbach's `alpha`. It returns bounds in closed form (no
  constructed distributions) together with the binding rule for each bound.
* `plot_umbrella()` now takes the output of `umbrella_data()`.

## New features

* Bounds are now sharp under many more constraint sets: attained observed
  extremes (`a`, `b`), the quasi-integer floor (defined at every mean, so bound
  curves are hole-free), the sharp Structure-S mean-conditional ceiling, and a
  reported Cronbach's alpha for sum- or mean-scored composites.
* `n_items` greater than 1 is now fully supported via the granularity grid
  (previously forced to 1).
* Rounded and truncated reported values are handled explicitly
  (`unround_interval()`, and a `rounding` argument with the usual scrutiny-style
  options); GRIM and GRIMMER verdicts are deferred to `scrutiny`.
* `sd_bounds_check()` gives a single consistent/inconsistent verdict, the
  failing tests, and percent-of-maximum-possible transforms (`pomp_mean`,
  `pomp_sd_parity`, `pomp_sd_sharp`). `sd_bounds_check_multiple()` applies it
  across a data frame, de-duplicating repeated constraint tuples.
* `sd_bounds_curve()` and `umbrella_data()` build the plotting data;
  `plot_sd_bounds()`, `plot_sd_bounds_pomp()` (native and POMP scales) and
  `plot_umbrella()` visualise it.
* The alpha-conditional floor is sharpened for strictly integer composites via
  the composite's Gini mean difference (`sd_min_alpha_gini()`, STRAIT article
  Theorem H5): a reported positive alpha now yields a strictly positive minimum
  SD even at whole-number sum-score means, where the earlier amplified floor
  vanished. The exact envelope is used when the composite range is small enough
  to enumerate; otherwise the proven amplified floor is retained (flagged in
  `note`).
* The single-purpose bound primitives (e.g. `sd_max_structure_s()`,
  `sd_min_quasi_integer()`) are exported and documented.
* `sd_bounds_sample()` constructs a sample that attains the maximum or minimum
  SD bound (the constructive companion to the closed-form `sd_bounds()`),
  replacing the old constructive `return_distributions` output.
* `v_max_alpha()` (the Appendix H allocation maximum of the item-variance sum)
  is now exported for reuse.
* `plot_sd_region()` and `sd_region_data()` draw the feasible SD region for any
  one constraint set in the nested framework, reproducing the panels of the
  article's nested-constraints figure from a single entry point: `"range"`,
  `"range_n"`, `"mean"`, `"mean_naive_floor"`, `"pesant_regin"`, `"mestdagh"`,
  `"quasi"`, `"alpha"`, plus two strictly integer views, `"integer"` and
  `"integer_alpha"`, which plot the lattice of GRIM- and GRIMMER-attainable
  reported tuples (no band exists at means integer data cannot produce). Each
  panel can carry the sharp alpha-free quasi-integer band as a dashed reference.
* `sd_max_muilwijk()` (the smooth Muilwijk / Bhatia-Davis arch) is exported.
* `umbrella_data()` gains an `alpha` argument, so the grid can be restricted to
  the jointly GRIM-, GRIMMER- and alpha-consistent tuples.

## Dependencies

* `Imports` trimmed to `ggplot2` and `scrutiny` (dropping `dplyr`, `tidyr`,
  `purrr`, `tibble`, `janitor`, `forcats`, `scales`, `rlang`).

# strait 0.3.2

This release prepares the package for CRAN and fixes several bugs that affected
the installed package.

## Bug fixes

* Core functions (`tides()`, `tides_df()`, `umbrella()`, `approximate_sd_bounds()`)
  no longer error when the package is installed: functions that were used but
  never imported (e.g. `dplyr::case_when()`, `dplyr::mutate_if()`,
  `dplyr::arrange()`) are now imported or fully qualified.
* `umbrella()` and the vignette now work across `scrutiny` versions. The
  GRIM/GRIMMER interface changed across `scrutiny` releases (older versions take
  zero-padded strings; newer ones take numbers plus their decimal places); the
  installed interface is now detected at run time (`grimmer_consistent()`).
* `plot_umbrella()` no longer warns about a missing `digits` column; it gains a
  `digits` argument (defaulting to a `digits` column when present, otherwise 2).

## New features and changes

* `tides_df()` now falls back to the `tides()` defaults when the optional
  columns `n_items`, `digits`, `calculate_min_sd` or `method` are absent, so an
  input data frame need not carry them.

## Documentation

* Rewrote the README with badges and problem / method / installation / usage /
  API / limitations / citation sections, and expanded the `DESCRIPTION`.
* Added a testthat test suite, a GitHub Actions `R-CMD-check` workflow,
  `cran-comments.md`, `inst/WORDLIST`, and `LICENSE.md`.
* Moved the package-level documentation to `R/tides-package.R` (from `zzz.R`),
  made all Rd content ASCII, and made examples runnable (`\donttest` where heavy).

## Known limitations

* `n_items` greater than 1 is not yet supported and is currently forced to 1
  internally in `tides()`, `tides_df()` and `umbrella()`.
