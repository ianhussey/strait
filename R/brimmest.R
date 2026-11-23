# Exact certification of reported (mean, SD) tuples against the attainable
# lattice of strictly integer data. See R/sd_region_plot.R for the dynamic
# program that enumerates the lattice.

# The rounding methods brimmest() certifies against. A certificate needs a
# definite direction, so the compound methods stand for their two parts.
BRIMMEST_ROUNDING <- c(
  "up_or_down",
  "ceiling_or_floor",
  "up",
  "down",
  "even",
  "ceiling",
  "floor",
  "trunc",
  "anti_trunc"
)

expand_rounding <- function(rounding) {
  rounding |>
    purrr::map(function(method) {
      switch(
        method,
        up_or_down = c("up", "down"),
        ceiling_or_floor = c("ceiling", "floor"),
        method
      )
    }) |>
    unlist() |>
    unique()
}

# Internal: forward-round lattice values to whole numbers of reporting-grid
# steps, exactly. A value y on that scale is given as the whole numbers whose
# ratio it is, (num / den for a mean, sqrt(num / den) for an SD), so every
# decision -- below, on, or above a half step -- is integer arithmetic.
# Scaling the double mean instead misfiles values on a boundary: 113/25 = 4.52
# scales to 451.99999999999994, which trunc() files under 4.51; and snapping to
# the nearest half step within a tolerance misfiles SDs that lie that close to
# one without being on it, which large samples produce.
#
# Returns `lo` and `hi`, the step counts with ties resolved down and up. They
# differ only under "even", which, like roundwork::unround(), lets a tie go
# either way: software does not reliably round decimal ties to even.
# Directions follow roundwork: "up" takes ties toward +Inf, "trunc" toward 0.
round_exact <- function(num, den, rounding, root = FALSE) {
  if (root) {
    whole <- floor(sqrt(num / den))
    whole <- whole + ((whole + 1)^2 * den <= num) - (whole^2 * den > num)
    on_int <- whole^2 * den == num
    half <- sign(4 * num - (2 * whole + 1)^2 * den)
  } else {
    whole <- floor(num / den)
    rem <- num - whole * den
    whole <- whole + (rem >= den) - (rem < 0)
    rem <- num - whole * den
    on_int <- rem == 0
    half <- sign(2 * rem - den)
  }
  nonneg <- num >= 0
  up <- whole + (half >= 0)
  down <- whole + (half > 0)
  steps <- switch(
    rounding,
    up = up,
    down = down,
    even = return(list(lo = down, hi = up)),
    ceiling = whole + !on_int,
    floor = whole,
    trunc = dplyr::if_else(nonneg, whole, whole + !on_int),
    anti_trunc = dplyr::if_else(nonneg, whole + !on_int, whole),
    cli::cli_abort("unknown rounding method: {.val {rounding}}")
  )
  list(lo = steps, hi = steps)
}

#' BRIMMEST: certify whether a reported (mean, SD) is attainable
#'
#' Bounds-Related Inconsistency of Means and Errors, Settled Test. An exact
#' possible / impossible certificate for reported summary statistics on a
#' bounded integer scale, obtained without reconstructing any dataset.
#'
#' `brimmest()` completes the family, and the superlative is the point.
#' [brim()] and [brimmer()] — like GRIM and GRIMMER before them — test
#' conditions that are *necessary but not sufficient*: failing one proves a
#' report impossible, but passing them all proves nothing. A handful of tuples
#' clear the SD bounds, GRIM and GRIMMER together and still correspond to no
#' integer dataset.
#'
#' `brimmest()` is *necessary and sufficient*. It enumerates the exact
#' attainable `(mean, sd)` lattice for the design by dynamic programming over
#' the reachable `(sum, sum of squares)` states, rounds that lattice to the
#' reporting precision, and asks whether the reported tuple is a member. The
#' verdict is therefore exact in both directions: a hit means some integer
#' sample on `[l, u]` of size `n` rounds to exactly this report, so it is
#' **possible**; a miss means no such sample exists, so it is **impossible** —
#' a proof, not a heuristic flag. Nothing in the GRIM family can make the
#' second claim.
#'
#' The three tests are nested, cheapest first:
#'
#' | test | asks | verdict |
#' | --- | --- | --- |
#' | [brim()] | is the reported mean attainable? | necessary only |
#' | [brimmer()] | and is the reported SD attainable with it? | necessary only |
#' | `brimmest()` | is the pair jointly attainable by real integer data? | necessary and sufficient |
#'
#' This is the same question the CLOSURE algorithm answers (see
#' `unsum::closure_generate()`), reached analytically rather than by search —
#' an *analytic CLOSURE certification*. The trade is that no witness datasets
#' are produced: the lattice records which `(mean, sd)` pairs are reachable,
#' not how. In exchange the cost depends only on `l`, `u` and `n`, not on how
#' many datasets happen to satisfy the constraints. The certification
#' validation document in `validation/` reports the comparison against CLOSURE
#' in full (it ships with the package's source repository, not the installed
#' package).
#'
#' @section How the work is done:
#' Two routes reach the same verdict, chosen automatically by design size.
#'
#' For a small design the whole attainable lattice is enumerated once and the
#' reports matched against it, so a vector of tuples costs barely more than a
#' single one. The lattice is cached per `(l, u, n)`, so repeated calls on one
#' design pay for it once.
#'
#' For a large design asked about a handful of reports, enumerating everything
#' to answer a membership question is wasteful — and on a wide scale it is
#' impossible, since the state table outgrows `max_cells`. A reported tuple
#' pins both coordinates almost completely: the mean's rounding interval
#' admits only a few integer sums, and each of those pins the sum of squares
#' to a narrow window. That reduces certification to a question about one
#' sample sum at a time, which is answered in three escalating steps.
#'
#' A closed-form sandwich comes first. For a given sum the achievable sums of
#' squares lie between the clustered configuration and the Structure-S one,
#' and all share the sum's parity; a window meeting no admissible integer
#' between them is impossible on arithmetic alone. What survives is decided by
#' the engine chosen via `engine`. The default `"analytic"` engine settles the
#' window against the structure of the attainable spectrum directly (the
#' demotion decomposition in `R/certify-analytic.R`): exact in both
#' directions, with no node budget and no fallback, so a verdict is always
#' produced. The `"search"` engine is the original machinery — a constructive
#' search over non-increasing samples with the sandwich re-applied at every
#' step, escalating to the corridor sweep (a dynamic program confined to the
#' states that can still reach a target) when `search_budget` runs out.
#' `"both"` runs the two independently and stops with an error if they ever
#' disagree, which is the recommended setting for validation work.
#'
#' The practical effect is that designs the lattice refuses outright become
#' answerable — a 0-63 inventory at `n = 50` needs a 3151 x 24801 table — and
#' that the overwhelming majority of single reports are settled without any
#' sweep at all.
#'
#' @section Rounding and the direction of proof:
#' A hit certifies possibility under whichever rounding rule produced it. A
#' miss certifies impossibility only *relative to the rules in `rounding`*,
#' because a report unreachable under one convention may be reachable under
#' another. The default, `"up_or_down"`, takes the union of `"up"` and `"down"`,
#' which spans the usual ambiguity about how a source rounded its halves.
#' Narrow it only when the source's convention is actually known: certifying
#' against one rule when the paper used another manufactures a false
#' impossibility. Rounding the lattice forward needs a definite direction, so
#' the compound methods `"up_or_down"` and `"ceiling_or_floor"` stand for
#' their two parts, and `"even"` lets a tie round either way, as
#' [roundwork::unround()] does, since software does not reliably round decimal
#' ties to even.
#'
#' @param l,u Numeric scalars, the scale limits. Integer-valued (in mean-score
#'   units when `n_items > 1` and `scoring` is `"meanscored"`).
#' @param n Integer scalar, sample size.
#' @param mean,sd Numeric vectors of reported means and SDs, recycled against
#'   each other. All are certified against one lattice. A missing value gives
#'   a missing verdict.
#' @param digits_mean,digits_sd Integer scalars, the reported decimal places
#'   of the mean and the SD, as in scrutiny.
#' @param rounding Character vector of rounding methods to admit, named as in
#'   scrutiny: `"up_or_down"` (default), `"ceiling_or_floor"`, `"up"`,
#'   `"down"`, `"even"`, `"ceiling"`, `"floor"`, `"trunc"`, `"anti_trunc"`. A
#'   tuple is possible if it is reachable under any of them.
#' @param scoring "singleitem" (default), "sumscored", or "meanscored".
#' @param n_items Positive whole number of response items (default 1).
#' @param max_cells Guard on the dynamic program's state space (default
#'   `2e7`); the enumeration errors rather than exhausting memory above it.
#' @param search_budget Node budget for the constructive search (default
#'   `2e5`). Raising it lets the search settle more reports on its own;
#'   lowering it hands them to the corridor sweep sooner. The verdict is the
#'   same either way — only which route produces it changes. Used by the
#'   `"search"` and `"both"` engines; the `"analytic"` engine has no budget.
#' @param engine How single reports are certified: `"analytic"` (default) for
#'   the exact spectrum-structure certification, `"search"` for the original
#'   constructive search with corridor-sweep fallback, or `"both"` to run the
#'   two independent methods and error on any disagreement (a built-in
#'   cross-check; also applied against the lattice when the lattice route is
#'   taken). All rounding rules are supported identically by every engine.
#' @return A tibble with one row per reported tuple: `mean`, `sd`,
#'   `possible` (logical), and `rules` — the rounding methods under which the
#'   tuple is reachable, comma-separated and `""` when none.
#' @seealso [brim()] and [brimmer()] for the closed-form screens to run
#'   first, and [sd_region_data()] for the lattice itself.
#' @examples
#' # a report a real 1-5 scale sample can produce
#' brimmest(l = 1, u = 5, n = 9, mean = 3.0, sd = 1.0, digits_mean = 1,
#'          digits_sd = 1)
#'
#' # inside the SD bounds and passing GRIM and GRIMMER, yet no integer
#' # sample produces it: the residual blind spot of the closed-form screen
#' brimmest(l = 1, u = 5, n = 9, mean = 1.3, sd = 0.9, digits_mean = 1,
#'          digits_sd = 1)
#'
#' # one lattice certifies many reports at once
#' brimmest(l = 1, u = 5, n = 9, digits_mean = 1, digits_sd = 1,
#'          mean = c(3.0, 1.3, 2.5), sd = c(1.0, 0.9, 1.2))
#' @export
brimmest <- function(
  l,
  u,
  n,
  mean,
  sd,
  digits_mean,
  digits_sd,
  rounding = "up_or_down",
  scoring = c("singleitem", "sumscored", "meanscored"),
  n_items = 1,
  max_cells = 2e7,
  search_budget = 2e5,
  engine = c("analytic", "search", "both")
) {
  scoring <- rlang::arg_match(scoring)
  engine <- rlang::arg_match(engine)
  rounding <- rounding |>
    check_rounding(BRIMMEST_ROUNDING, several = TRUE) |>
    expand_rounding()
  if (missing(digits_mean) || missing(digits_sd)) {
    cli::cli_abort(
      "reported decimal places are required: give {.arg digits_mean} and {.arg digits_sd}"
    )
  }
  check_scalar(digits_mean, "digits_mean", whole = TRUE, min = 0)
  check_scalar(digits_sd, "digits_sd", whole = TRUE, min = 0)
  check_scalar(l, "l")
  check_scalar(u, "u")
  check_scalar(n, "n", whole = TRUE, min = 2)
  check_scalar(n_items, "n_items", whole = TRUE, min = 1)

  k <- as.integer(round(n_items))
  geometry <- scoring_geometry(scoring, k, l, u)
  mg <- geometry$mg
  # both routes work on the integer grid y = mg * (x - l) in 0..W
  if (
    abs(mg * l - round(mg * l)) > TOLERANCE ||
      abs(mg * u - round(mg * u)) > TOLERANCE ||
      u <= l
  ) {
    cli::cli_abort(
      "{.arg l} and {.arg u} must lie on the response grid (whole numbers, or multiples of 1/{.arg n_items} for mean-scored data), with {.code u > l}"
    )
  }

  if (!is.numeric(mean) || !is.numeric(sd)) {
    cli::cli_abort("{.arg mean} and {.arg sd} must be numeric")
  }
  n_reports <- max(length(mean), length(sd))
  mean <- rep(mean, length.out = n_reports)
  sd <- rep(sd, length.out = n_reports)
  if (any(sd < 0, na.rm = TRUE)) {
    cli::cli_abort("{.arg sd} must not be negative")
  }

  # scrutiny's NA rule: a missing mean or SD gets a missing verdict
  ok <- !is.na(mean) & !is.na(sd)
  if (!all(ok)) {
    out <- tibble::tibble(
      mean = mean,
      sd = sd,
      possible = NA,
      rules = NA_character_
    )
    if (any(ok)) {
      known <- brimmest(
        l,
        u,
        n,
        mean[ok],
        sd[ok],
        digits_mean,
        digits_sd,
        rounding = rounding,
        scoring = scoring,
        n_items = n_items,
        max_cells = max_cells,
        search_budget = search_budget,
        engine = engine
      )
      out$possible[ok] <- known$possible
      out$rules[ok] <- known$rules
    }
    return(out)
  }
  # reject values with more decimals than stated, as both routes assume
  unround_reported(mean, digits_mean, "up", "mean", "digits_mean")
  unround_reported(sd, digits_sd, "up", "sd", "digits_sd")

  # Two routes to the same verdict, chosen by how much work each implies.
  #
  # Targeted (see R/attainable-target.R): reachability of just the states the
  # report pins down. Cost is per tuple, so it wins for a handful of reports,
  # and it is the only route on wide scales, where the full lattice exceeds
  # its own size guard.
  #
  # Lattice: enumerate every attainable pair once and match against it. Cost
  # is per design regardless of how many tuples are asked about, so it wins
  # once there are enough of them, as when brimmest() is handed a whole grid.
  # A small lattice is so cheap to build, and now cached, that it beats the
  # targeted route even for one tuple; the targeted route earns its overhead
  # only once the lattice is genuinely large, and is the sole option once the
  # lattice will not fit at all.
  W <- as.integer(round(mg * (u - l)))
  cells <- lattice_cells(l, u, n, mg)
  affordable <- !is.na(cells) && cells <= max_cells
  use_target <- !affordable ||
    (cells > 1e6 && n_reports * length(rounding) <= 64)

  targets_of <- function(i, rule) {
    target_states(l, u, n, mg, mean[i], sd[i], digits_mean, digits_sd, rule)
  }
  abort_disagreement <- function(i, rule, verdict_analytic, other, verdict) {
    cli::cli_abort(c(
      "engine disagreement at mean = {mean[i]}, sd = {sd[i]}, rounding = {.val {rule}}: analytic says {verdict_analytic}, {other} says {verdict}.",
      "i" = "This is a bug; please report it with these inputs."
    ))
  }

  # One target certified under one rounding rule, by the requested engine.
  # "analytic" is total (never NA); "search" keeps the original escalation
  # (constructive search, then corridor sweep); "both" runs the two
  # independent methods and refuses to answer if they disagree.
  certify_one <- function(i, rule) {
    targets <- targets_of(i, rule)
    verdict_analytic <- if (engine != "search") {
      certify_analytic(W, n, targets)
    } else {
      NA
    }
    verdict_search <- NA
    if (engine != "analytic") {
      verdict_search <- certify_fast(W, n, targets, budget = search_budget)
      if (is.na(verdict_search)) {
        verdict_search <- attainable_target(
          W,
          n,
          targets,
          max_cells = max_cells
        )
      }
      if (is.na(verdict_search) && engine == "search") {
        cli::cli_abort(c(
          "this design is too large to certify at the requested precision",
          "i" = "Raise {.arg max_cells}, or report to fewer decimal places."
        ))
      }
    }
    if (
      engine == "both" &&
        !is.na(verdict_search) &&
        verdict_analytic != verdict_search
    ) {
      abort_disagreement(i, rule, verdict_analytic, "search", verdict_search)
    }
    if (engine == "search") verdict_search else verdict_analytic
  }

  if (use_target) {
    hits <- rounding |>
      purrr::map(function(rule) {
        purrr::map_lgl(seq_len(n_reports), function(i) certify_one(i, rule))
      }) |>
      unlist()
  } else {
    # one lattice per design, reused across every reported tuple and cached so
    # that repeated calls on the same design pay for it once
    lattice <- lattice_cached(l, u, n, mg, max_cells)
    # Both sides are on the 10^-digits reporting grid, so compare them as
    # whole numbers of grid steps rather than as strings: exact, and it keeps
    # a multi-million-row lattice cheap to match against. Integers also
    # sidestep the negative zero some rounding rules return at 0, which
    # compares equal but formats differently.
    report_mean_steps <- round(mean * 10^digits_mean)
    report_sd_steps <- round(sd * 10^digits_sd)
    # the lattice's mean and SD on the reporting grid, as whole-number ratios:
    # mean * 10^digits_mean = (l * n * mg + S) * 10^digits_mean / (n * mg),
    # and (sd * 10^digits_sd)^2 = nss * 10^(2 * digits_sd) /
    # (n * (n - 1) * mg^2)
    m_num <- (round(l * mg) * n + lattice$S) * 10^digits_mean
    m_den <- n * mg
    s_num <- lattice$nss * 10^(2 * digits_sd)
    s_den <- n * (n - 1) * mg^2
    if (max(4 * s_num, abs(m_num)) >= 2^53) {
      cli::cli_abort(c(
        "the lattice is too fine to round exactly at this precision",
        "i" = "Report to fewer decimal places."
      ))
    }
    hits <- rounding |>
      purrr::map(function(rule) {
        mean_steps <- round_exact(m_num, m_den, rule)
        sd_steps <- round_exact(s_num, s_den, rule, root = TRUE)
        span <- max(c(sd_steps$hi, report_sd_steps), 0) + 1 # same packing
        key <- report_mean_steps * span + report_sd_steps
        # under "even" a tie may round either way, in the mean and SD alike
        key %in%
          (mean_steps$lo * span + sd_steps$lo) |
          key %in% (mean_steps$lo * span + sd_steps$hi) |
          key %in% (mean_steps$hi * span + sd_steps$lo) |
          key %in% (mean_steps$hi * span + sd_steps$hi)
      }) |>
      unlist()
  }
  hit_mat <- matrix(hits, n_reports, length(rounding))

  if (!use_target && engine == "both") {
    # cross-check the lattice verdicts against the analytic engine, tuple by
    # tuple -- two entirely independent computations of the same answer
    for (j in seq_along(rounding)) {
      for (i in seq_len(n_reports)) {
        verdict_analytic <- certify_analytic(W, n, targets_of(i, rounding[j]))
        if (verdict_analytic != hit_mat[i, j]) {
          abort_disagreement(
            i,
            rounding[j],
            verdict_analytic,
            "lattice",
            hit_mat[i, j]
          )
        }
      }
    }
  }

  tibble::tibble(
    mean = mean,
    sd = sd,
    possible = rowSums(hit_mat) > 0,
    rules = purrr::map_chr(seq_len(n_reports), function(i) {
      paste(rounding[hit_mat[i, ]], collapse = ",")
    })
  )
}

# ---- Batch certification -----------------------------------------------------

#' Certify many reported (mean, SD) rows exactly
#'
#' Applies [brimmest()] to each row of a data frame. Columns of `data` whose
#' names match [brimmest()] arguments are taken per row; any argument given in
#' `...` is a constant broadcast to every row. Supplying one name both ways is
#' an error. The counterpart of [brimmer_map()] for the exact test.
#'
#' Rows are **grouped by design** and each group certified in a single
#' [brimmest()] call, which is the whole point of the function: one attainable
#' lattice serves every report of a design and is built once, so certifying a
#' 200-row table on one design costs barely more than certifying one row of it.
#' Repeated `(mean, sd)` pairs within a design are computed once and reused.
#' Row order is preserved.
#'
#' Recognised per row (column or constant): `l`, `u`, `n`, `mean`, `sd`,
#' `digits_mean`, `digits_sd`, `scoring`, `n_items`.
#'
#' `rounding`, `max_cells`, `search_budget` and `engine` may be given only as
#' constants, never as columns. `rounding` is a *set* of admitted rules rather
#' than one value per row, and the others are controls on the call; passing
#' any of them as a column is an error rather than a silent per-row reading.
#'
#' @param data A data frame, one reported statistic set per row.
#' @param ... Constant arguments applied to all rows.
#' @param include_inputs If TRUE (default), returns `data` column-bound to the
#'   results; if FALSE, only the result columns (same row order), for drop-in
#'   use inside a `dplyr::mutate()` / `purrr` pipeline.
#' @return A tibble with `possible` (logical) and `rules` (the rounding
#'   rules under which the row is reachable, comma-separated and `""` when
#'   none), optionally with `data` prepended. Unlike [brimmest()], the reported
#'   `mean` and `sd` are not echoed back: they are already columns of `data`
#'   when supplied that way, and repeating them would collide.
#' @seealso [brimmest()] for one design at a time, [brimmer_map()] for the
#'   closed-form screen to run first.
#' @examples
#' reports <- data.frame(mean = c(3.0, 1.3, 2.5),
#'                       sd   = c(1.0, 0.9, 1.2))
#'
#' # one design, one lattice. Rows 2 and 3 are both impossible, for different
#' # reasons: row 3 already fails GRIM and GRIMMER, while row 2 clears every
#' # closed-form screen and is caught only here.
#' brimmest_map(reports, l = 1, u = 5, n = 9, digits_mean = 1,
#'                   digits_sd = 1)
#'
#' # designs may vary per row; rows sharing one are certified together
#' mixed <- data.frame(mean = c(3.0, 1.3, 4.0),
#'                     sd   = c(1.0, 0.9, 1.5),
#'                     n    = c(9, 9, 12))
#' brimmest_map(mixed, l = 1, u = 5, digits_mean = 1, digits_sd = 1,
#'                   include_inputs = FALSE)
#' @export
brimmest_map <- function(data, ..., include_inputs = TRUE) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame")
  }
  row_args <- c(
    "l",
    "u",
    "n",
    "mean",
    "sd",
    "digits_mean",
    "digits_sd",
    "scoring",
    "n_items"
  )
  call_args <- c("rounding", "max_cells", "search_budget", "engine")
  constants <- list(...)
  unknown <- setdiff(names(constants), c(row_args, call_args))
  if (length(unknown)) {
    cli::cli_abort("unknown constant argument{?s}: {.arg {unknown}}")
  }
  clash <- intersect(call_args, names(data))
  if (length(clash)) {
    cli::cli_abort(
      "{.arg {clash}} appl{?ies/y} to the whole call, not to one row; supply {?it/them} as a constant rather than a column"
    )
  }

  n_rows <- nrow(data)
  if (!n_rows) {
    cli::cli_abort("{.arg data} has no rows")
  }
  resolve <- function(name) {
    in_data <- name %in% names(data)
    in_constants <- name %in% names(constants) && !is.null(constants[[name]])
    if (in_data && in_constants) {
      cli::cli_abort(
        "{.arg {name}} supplied as both a column and a constant"
      )
    }
    values <- if (in_data) {
      data[[name]]
    } else if (in_constants) {
      rep(constants[[name]], length.out = n_rows)
    }
    if (is.factor(values)) as.character(values) else values
  }
  cols <- row_args |>
    rlang::set_names() |>
    purrr::map(resolve)
  present <- names(purrr::discard(cols, is.null))
  for (name in c("l", "u", "n", "mean", "sd", "digits_mean", "digits_sd")) {
    if (!(name %in% present)) {
      cli::cli_abort(
        "{.arg {name}} is required (as a column of data or a constant)"
      )
    }
  }

  passthru <- constants[intersect(call_args, names(constants))]

  # Group by design, not by row. brimmest() amortises one lattice across every
  # tuple it is handed, and it also uses the number of tuples to choose its
  # route, so handing it a whole group at once is both cheaper and better
  # routed than calling it per row would be.
  design <- setdiff(present, c("mean", "sd"))
  design_key <- cols[design] |>
    purrr::map(format, nsmall = 6, trim = TRUE) |>
    purrr::reduce(paste, sep = "\r")

  possible <- rep(NA, n_rows)
  rules <- rep(NA_character_, n_rows)
  for (key_current in unique(design_key)) {
    rows <- which(design_key == key_current)
    args <- purrr::map(cols[design], function(col) col[rows[1]])
    means <- cols$mean[rows]
    sds <- cols$sd[rows]
    tuple_key <- paste(
      format(means, nsmall = 6, trim = TRUE),
      format(sds, nsmall = 6, trim = TRUE),
      sep = "\r"
    )
    first <- !duplicated(tuple_key)
    certified <- do.call(
      brimmest,
      c(args, list(mean = means[first], sd = sds[first]), passthru)
    )
    back <- match(tuple_key, tuple_key[first])
    possible[rows] <- certified$possible[back]
    rules[rows] <- certified$rules[back]
  }

  out <- tibble::tibble(possible = possible, rules = rules)
  if (include_inputs) bind_inputs(data, out) else out
}
