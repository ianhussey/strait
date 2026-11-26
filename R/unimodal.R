# The unimodal ceiling: the largest sample SD of integer data whose frequency
# distribution is unimodal. The floor needs no counterpart, because the
# smallest SD at every mean is attained by a unimodal sample: its
# observations lie on the two integers bracketing the mean.

# Guard on the recursion's tables, in cells summed over all blocks. A design
# above it would take minutes, so it is refused.
MAX_UNIMODAL_CELLS <- 5e8

# Internal: the recursion's results, keyed by n and the number of steps
# between the limits. sd_bounds_curve() and the rounding envelope ask for many
# means at the same n and scale, and one table answers all of them.
unimodal_cache <- new.env(parent = emptyenv())

# Internal: the largest sum of squared positions among the unimodal samples of
# n observations on the positions 0, ..., top, for every sum of positions
# 0, ..., n * top. Element Y + 1 holds the value for the sum Y.
#
# Every unimodal sample is a stack of nested blocks of neighbouring positions,
# each block placing one observation on each of its positions. For the block
# a..b, a table holds the largest sum of squares among the stacks whose blocks
# all contain a..b, with a row for each total number of positions nu and a
# column for each sum of positions tau. A stack either has a..b as its
# smallest block, or all its blocks contain a - 1..b or a..b + 1, so each
# table follows from the tables of the two blocks one position longer.
unimodal_max_sum_of_squares <- function(n, top) {
  key <- paste(n, top)
  if (!is.null(unimodal_cache[[key]])) {
    return(unimodal_cache[[key]])
  }

  cells <- (top + 1) * (top + 2) / 2 * (n + 1) * (n * top + 1)
  if (cells > MAX_UNIMODAL_CELLS) {
    cli::cli_abort(paste(
      "the unimodal ceiling for n = {n} on {top + 1} scale values is too",
      "large to compute ({format(cells, big.mark = ',')} table cells)"
    ))
  }

  # Only the empty stack.
  empty_stack <- matrix(-Inf, n + 1, n * top + 1)
  empty_stack[1, 1] <- 0

  tables <- list()
  for (len in (top + 1):1) {
    for (a in 0:(top + 1 - len)) {
      b <- a + len - 1
      block_sum <- sum(a:b)
      block_squares <- sum((a:b)^2)

      # Stacks without the block a..b: all blocks contain a - 1 or b + 1.
      current <- list(
        tables[[paste(a - 1, b)]],
        tables[[paste(a, b + 1)]]
      ) |>
        purrr::compact() |>
        purrr::reduce(pmax, .init = empty_stack)

      # Stacks with the block a..b: one such block on top of a stack with len
      # fewer positions.
      cols <- seq_len(max(n * top + 1 - block_sum, 0))
      for (row in seq_len(max(n + 1 - len, 0))) {
        current[row + len, cols + block_sum] <- pmax(
          current[row + len, cols + block_sum],
          current[row, cols] + block_squares
        )
      }

      tables[[paste(a, b)]] <- current
    }
  }

  result <- purrr::reduce(tables[paste(0:top, 0:top)], pmax)[n + 1, ]
  unimodal_cache[[key]] <- result
  result
}

# Internal: the largest sum of squared deviations of a unimodal sample of n
# observations on the positions 0, ..., top, over all means. With
# W = min(top + 1, n) and B = c * W * (W - 1), it is the largest value of
# (2W - 1) * B / 6 - B^2 / (4n) over the whole numbers 1 <= c <= n / W,
# attained by c observations on each of the W - 1 positions below the top
# and the remaining n - c (W - 1) on the top.
unimodal_ssd_max_any_mean <- function(n, top) {
  w <- min(top + 1, n)
  count <- seq_len(n %/% w)
  b <- count * w * (w - 1)

  max((2 * w - 1) * b / 6 - b^2 / (4 * n))
}

#' Maximum SD of unimodal integer data
#'
#' The sharp ceiling on the sample SD of `n` whole numbers in
#' `[lower, upper]` whose frequency distribution is unimodal: the frequencies
#' of the values `lower, lower + 1, ..., upper` rise, or stay level, up to a
#' mode and then fall, or stay level, where a value that no observation takes
#' counts as frequency 0. Unimodality lowers the ceiling of
#' [sd_max_structure_s()], which puts the observations at the two limits.
#'
#' Given the mean, the ceiling is computed exactly by a recursion over the
#' layers of the frequency distribution: every unimodal sample is a stack of
#' nested runs of neighbouring values, each adding one observation to each of
#' its values. Without the mean, it is a closed-form maximum over all means.
#' The recursion's cost grows with `n^2` and the cube of the
#' number of scale values; it answers every mean of a design at once, and the
#' answer is cached for the session. Defined only at GRIM-consistent means; the
#' caller is responsible for that check (see [is_grim_consistent()]), and the
#' result is `NA` at any other mean.
#'
#' For data on a finer grid, such as means of several integer items, rescale
#' the data to whole numbers first, as [sd_bounds()] does.
#'
#' @param mean Numeric vector, mean(s) in `[lower, upper]`, or NULL for the
#'   maximum over all means.
#' @param n Integer scalar, sample size, `n >= 2`.
#' @param lower,upper Whole numbers, the limits of the scale.
#' @return Numeric vector, the maximum sample SD at each mean (one value if
#'   `mean` is NULL).
#' @examples
#' # 30 single-item responses on a 1-5 scale with a mean of 4.5
#' sd_max_unimodal(mean = 4.5, n = 30, lower = 1, upper = 5)
#'
#' # without unimodality, the observations would sit at the two limits
#' sd_max_structure_s(mean = 4.5, n = 30, lower = 1, upper = 5)
#'
#' # over all means, at n = 10
#' sd_max_unimodal(n = 10, lower = 1, upper = 5)
#' @export
sd_max_unimodal <- function(mean = NULL, n, lower, upper) {
  top <- upper - lower
  if (is.null(mean)) {
    return(sd_from_ss(unimodal_ssd_max_any_mean(n, top), n))
  }

  sums <- n * (mean - lower)
  whole <- abs(sums - round(sums)) < TOLERANCE &
    sums > -TOLERANCE &
    sums < n * top + TOLERANCE
  y <- round(sums[whole])
  table <- unimodal_max_sum_of_squares(n, top)

  out <- rep(NA_real_, length(mean))
  out[whole] <- sd_from_ss(pmax(table[y + 1] - y^2 / n, 0), n)
  out
}
