# strait 0.5.0

## New features

* `umbrella_plot(style = "contour")` draws the outline of the consistent
  tuples rather than the tuples, and `umbrella_contour()` returns that outline
  as data. `umbrella_plot(shade = "none")` leaves the panel unshaded.

* `implied_range()` finds the range a consistent report forces the data to
  cover: the smallest observation is at most `min_at_most` and the largest at
  least `max_at_least`. An SD near the ceiling of a 1-7 scale, for instance,
  forces both a 1 and a 7 into the data. Each bound is where `brimmer()` stops
  calling the report consistent as one scale limit moves inward, so the two
  always agree.

## Breaking changes

* `Z` is now `granularity`, and `brimmer_multiple()` and
  `brimmest_multiple()` are now `brimmer_map()` and `brimmest_map()`, after
  scrutiny's mappers such as `grim_map()`.

* Arguments, rounding methods and return values follow scrutiny:
  * Decimal places are `digits_mean` and `digits_sd`, as `digits_x` and
    `digits_sd` in scrutiny, instead of `mean_digits` and `sd_digits`, in
    `sd_bounds()`, `brim()`, `brimmer()`, `brimmer_map()`,
    `implied_range()`, `brimmest()` and `brimmest_map()`. `brimmest()`
    no longer takes `digits` for both; it requires `digits_mean` and
    `digits_sd`.
  * `rounding` takes scrutiny's names everywhere. `brimmest()` and
    `brimmest_map()` take `"up_or_down"` (the new default, standing for
    `"up"` and `"down"`), `"ceiling_or_floor"`, `"up"`, `"down"`, `"even"`,
    `"ceiling"`, `"floor"`, `"trunc"` and `"anti_trunc"`, instead of
    `"half_up"`, `"half_down"` and `"native"`; `sd_region_data()` takes one of
    these, with `"up"` as the default. The old names are an error that names
    the new one.
  * Reported values are unrounded by `roundwork::unround()`, so strait and
    scrutiny agree on every method. `unround_interval()` is removed: use
    `roundwork::unround()`. `infer_digits()` is removed: use
    `scrutiny::decimal_places()`.
  * Every function that returns a table returns a tibble.
  * `umbrella_data()` takes the scale first, like every other function:
    `umbrella_data(l, u, n, ...)` instead of `umbrella_data(n, l, u, ...)`.

* As in scrutiny, a missing input gives a missing verdict. `brimmer()`,
  `brim()`, `brimmer_map()` and `implied_range()` return `NA` with a note
  naming the input, and `brimmest()` returns `possible = NA` for a missing mean
  or SD. Before, a missing SD passed `brimmer()` as consistent, or failed it as
  GRIMMER-inconsistent under `granularity = "integer"`.

* Inputs are validated with clear errors: `n` must be one whole number of at
  least 2, scale limits and reported values single finite numbers, SDs
  non-negative, and a reported value must not have more decimal places than
  its `digits_*` states. `brimmest()` now checks this on its targeted route
  too, where a non-integer scale width, a non-integer `n` or a negative SD used
  to be accepted.

* `sd_plausibility()`, `rho_from_summary()`, `s2_exact()` and
  `rho_corpus_fit()` use strait's argument names:
  `sd_plausibility(l, u, n, mean, sd = NULL, ..., digits_mean = NULL)`,
  `rho_from_summary(l, u, n, mean, sd, ...)` and `s2_exact(m, mu, rho, n,
  sums)`. `rho_corpus_fit()` expects corpus columns `l`, `u`, `n`, `mean` and
  `sd`, and scales may differ between rows. Positional calls in the old
  `(L, U, M, N)` order now swap `n` and `mean`. `s2_exact()` returns tibbles.

* `sd_plausibility()` conditions on the reported mean by default
  (`condition_on_mean = TRUE`), as recommended for checking a reported
  (mean, SD) pair. `"chisq3"` and `"cf"` cannot condition, so they need
  `condition_on_mean = FALSE`.

* `grim_aggregate()`, `grim_aggregate_by()` and `grim_sensitivity()` follow
  scrutiny 1.0: `x` must be numeric, with `digits_x` given as an argument or a
  `data` column (not both); a character `x` is an error that points to
  `scrutiny::decimal_places()`. `items` is now `items_averaged` (argument and
  column) and has no default. `percent` is removed. New `threshold` and
  `symmetric` arguments are passed to scrutiny with `rounding`. The `$data`
  column `digits` is now `digits_x`.

* The plotting functions are named after the data they draw, as the data
  builders are: `plot_sd_bounds()` is now `sd_bounds_plot()`,
  `plot_sd_bounds_pomp()` is `sd_bounds_pomp_plot()`, `plot_sd_region()` is
  `sd_region_plot()`, and `plot_umbrella()` is `umbrella_plot()`. Each takes
  what it draws as its first argument, `data`.

* `sd_region_plot()` takes the output of `sd_region_data()` instead of computing
  it from `l`, `u`, `n`, and `rule`. The dashed reference band, which it used to
  add by default, is now explicit: `reference = sd_region_data(l, u, n, rule =
  "quasi")`.

* `umbrella_plot()` takes `data` and `bounds` instead of `umbrella` and
  `curve`. Given a band, with `min_sd` and `max_sd` (such as
  `sd_bounds_curve()` output), instead of a lattice, it fills the umbrella
  itself, in the new `fill` colour; `sd_region_data(rule = "mean_naive_floor")`
  gives a smooth one. `bounds` are drawn as dashed lines in `reference_colour`
  under both lattice styles.

* `sd_bounds_plot()` and `sd_bounds_pomp_plot()` take the reported values as
  `reports` instead of `points`. `sd_bounds_plot()` also takes a band from
  `sd_region_data()`.

* Bands from `sd_region_data()`, and the input of `band_polygon()`, name their
  columns `min_sd` and `max_sd` instead of `lo` and `hi`, as in
  `sd_bounds_curve()` and `umbrella_data()`. `sd_region_data()` records its grid
  spacing in a `"step"` attribute.

* strait requires scrutiny 1.0.0 or later, and imports roundwork, to which
  scrutiny's rounding functions moved. The run-time switch between scrutiny's
  old and new interfaces is gone.

## Bug fixes

* `brimmest()` gave different verdicts on its two routes, so `engine = "both"`
  stopped with "This is a bug":
  * Under banker's rounding (`"even"`, formerly `"native"`), the targeted route
    let a tie round either way while the lattice route rounded it to even. Both
    now let a tie round either way, as `roundwork::unround()` does, because
    software does not reliably round decimal ties to even.
  * Under `"trunc"`, a reported 0 was unrounded to `[0, unit)`, though
    truncation sends all of `(-unit, unit)` to 0, so a real sample with a small
    negative mean was called impossible. `"anti_trunc"` had the mirror error and
    admitted impossible reports.

* `brimmest()`'s lattice route rounds the attainable means and SDs exactly, in
  integer arithmetic. It used to snap values within 1e-9 of a rounding tie onto
  it, which misfiles SDs that lie that close to a tie without being on it, as
  large samples produce. Its cache keeps at most four designs.

* `brimmer_map()` reads factor columns as strings, and passes a constant
  given as `NULL`, such as `rounding = NULL`, instead of replacing it with
  `brimmer()`'s default.

* `implied_range()` moves a limit along the `1 / n_items` grid of mean-scored
  data rather than in whole numbers, which gave looser bounds. Its help page no
  longer calls the results exact: they are the range `brimmer()` can prove,
  and real data may be forced to cover more.

* `sd_bounds()` accepts limits on the `1 / n_items` grid of mean-scored data
  under `granularity = "integer"` and `"quasiinteger"`, not only whole numbers.

* `umbrella_data()` tests each SD with its own rounding interval under the
  chosen `rounding`, rather than a symmetric one.

* `sd_max_span_n()` is vectorised over `n`.

* `grim_aggregate()`:
  * Missing values in `x`, `n`, `digits_x` or `items_averaged`, an infinite
    `n`, and an empty `data` are errors. Before, a malformed mean such as
    `"3,14"` silently passed GRIM.
  * If even lambda = 1 cannot explain the observed failures, the confidence
    interval for lambda is `NA` with a warning that the model is rejected,
    instead of collapsing to `[1, 1]`.
  * `conf_level` must lie in (0, 1) and `support_c` above 1. The interval label
    keeps decimals: `lambda_ci_97.5`, not `lambda_ci_98`.
  * `q` and `mean_q` are exactly 0 where every ending passes. A function `eps`
    is called once, and the `eps` column is present even when no value is
    informative.

* `grim_aggregate_by()` keeps rows whose grouping value is `NA` as their own
  group instead of dropping them, and errors if a `by` column clashes with an
  output column such as `F` or `k`.

* `sd_plausibility()`:
  * `digits_mean` is recycled with `mean`. Before, every row after the first
    got wrong hard bounds.
  * A mean at a scale limit given with `sd` returns zero bounds, as with
    `rho`, instead of "missing input".
  * On a two-point scale (`u - l == 1`), `rho` is `NA` with a note that it is
    not identified; `rho_from_summary()` errors there.
  * Giving both `rho` and `sd` is an error rather than ignoring `sd`.
  * `method = "exact"` computes the distribution once per row and node and
    reuses it for every `p`.
  * The help page gains an Accuracy section: with strong skew, the `"chisq3"`
    and `"cf"` lower bounds can be anti-conservative.

* The plots use scrutiny's colours for consistent and inconsistent values,
  which colour-blind readers can tell apart.

* Under `scoring = "sumscored"`, the GRIM and GRIMMER verdicts of `sd_bounds()`,
  `brimmer()`, and `umbrella_data()` treated each sum score as the average of
  `n_items` items, which made them too lenient. A sum score is one whole number
  per person.

* The plot functions no longer draw across a gap in a band, at means where no
  SD is possible. `sd_bounds_plot(shade = "inside")` and the parity band of
  `sd_bounds_pomp_plot()` filled it in, and their outlines, like the dashed
  reference lines of `sd_region_plot()` and `umbrella_plot()`, joined its two
  sides with a straight line.

* `sd_max_structure_s()` lost precision next to a scale limit, where it
  computed a tiny sum of squares as the difference of large ones. Within about
  1e-8 of a limit it returned 0, below the quasi-integer floor, so `sd_bounds()`
  called such a mean infeasible; `sd_bounds_curve()` samples two of them, next
  to the kinks at `l` and `u`.

* `sd_region_plot()` draws only the consistent tuples of a lattice with a
  `consistent` column, such as `umbrella_data()` output, as `umbrella_plot()`
  does, rather than every cell of the grid.

# strait 0.4.13

## New features

* `sd_plausibility()` gives plausibility bounds for a reported SD on a
  bounded integer scale: quantiles of the sample SD's sampling distribution
  under a beta-binomial data-generating process, next to the hard integer
  bounds (to which they are clipped, with a note). Four methods: a
  two-moment scaled chi-square (default), a three-moment shifted
  chi-square, Cornish-Fisher, and the exact distribution by dynamic
  programming over (sum X, sum X^2) for small N. `condition_on_mean = TRUE`
  conditions on the sample mean -- which is what checking a reported
  (M, SD) pair asks -- using Var(S^2 | mean) = sigma^4 [2/(N-1) +
  (gamma_2 - gamma_1^2)/N], or exactly. Between-study heterogeneity in rho
  enters as a mixture over a truncated normal, averaged on the CDF scale.
* Supporting functions: `bb_cumulants()`, `s2_exact()`,
  `rho_from_summary()` (method-of-moments rho and its sampling variance),
  and `rho_corpus_fit()` with a `predict()` method (a metafor
  meta-regression of rho on mu from summary statistics alone, returning
  rho, its SE, tau and tau_total). Method, validation and caveats are in
  `paper/sd-plausibility-bounds/`.

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
