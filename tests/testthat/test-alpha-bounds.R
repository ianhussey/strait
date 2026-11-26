# The bounds on the item-variance sum of a k-item integer battery under the
# quasi-integer constraint against brute force over allocations of the item
# means. Each allocation is evaluated with the single-item quasi-integer
# maximum and minimum, not with the composite formulas, so the search is
# independent of the results it checks.

item_l <- 1
item_u <- 7

# The single-item quasi-integer maximum and minimum population variances.
item_var_max <- function(mean, n) {
  sd_max_structure_s(mean, n, item_l, item_u)^2 * (n - 1) / n
}

item_var_min <- function(mean, n) {
  sd_min_quasi_integer(mean, n)^2 * (n - 1) / n
}

clamp <- function(x) pmin(pmax(x, item_l), item_u)

test_that("v_max_alpha() is a valid and attained ceiling", {
  skip_on_cran() # brute force over item-mean allocations
  for (design in list(c(5, 2), c(10, 2), c(7, 3))) {
    n <- design[1]
    k <- design[2]
    sweeps <- 200
    sum_means <- k * item_l + k * (item_u - item_l) * (1:sweeps) / (sweeps + 1)

    for (mean_sum in sum_means) {
      ceiling <- v_max_alpha(mean_sum, k, n, item_l, item_u)

      # Validity: no allocation of the item means exceeds the ceiling.
      if (k == 2) {
        m1 <- item_l + (item_u - item_l) * (0:4000) / 4000
        m2 <- mean_sum - m1
        ok <- m2 >= item_l - 1e-12 & m2 <= item_u + 1e-12
        brute <- max(item_var_max(m1[ok], n) + item_var_max(clamp(m2[ok]), n))
      } else {
        grid <- item_l + (item_u - item_l) * (0:200) / 200
        allocation <- expand.grid(m1 = grid, m2 = grid)
        allocation$m3 <- mean_sum - allocation$m1 - allocation$m2
        allocation <- allocation[
          allocation$m3 >= item_l - 1e-12 & allocation$m3 <= item_u + 1e-12,
        ]
        brute <- max(
          item_var_max(allocation$m1, n) +
            item_var_max(allocation$m2, n) +
            item_var_max(clamp(allocation$m3), n)
        )
      }
      expect_lte(brute, ceiling + 1e-7)

      # Sharpness: k - 1 two-point items at the means bracketing
      # mean_sum / k and one remainder item reach the ceiling.
      h <- (item_u - item_l) / n
      t <- (mean_sum - k * item_l) / h
      j <- pmin(floor(t / k + 1e-12), n - 1)
      m <- pmin(floor(t - k * j + 1e-12), k - 1)
      low <- item_l + j * h
      high <- item_l + (j + 1) * h
      remainder <- mean_sum - (k - 1 - m) * low - m * high
      attained <- (k - 1 - m) * item_var_max(low, n) +
        m * item_var_max(high, n) +
        item_var_max(clamp(remainder), n)
      expect_equal(attained, ceiling, tolerance = 1e-8)
    }
  }
})

test_that("the sum-score quasi-integer floor is the exact minimum", {
  skip_on_cran() # brute force over item-mean allocations
  for (design in list(c(5, 2), c(10, 2), c(7, 3))) {
    n <- design[1]
    k <- design[2]
    sweeps <- 200
    sum_means <- k * item_l + k * (item_u - item_l) * (1:sweeps) / (sweeps + 1)
    # k - 1 strictly integer items (means on the 1/n grid) and one item that
    # carries the free observation
    grid <- item_l + (0:((item_u - item_l) * n)) / n

    for (mean_sum in sum_means) {
      floor_sum <- item_var_min(mean_sum, n)
      if (k == 2) {
        m2 <- mean_sum - grid
        ok <- m2 >= item_l - 1e-12 & m2 <= item_u + 1e-12
        d1 <- grid[ok] - floor(grid[ok])
        brute <- min(d1 * (1 - d1) + item_var_min(m2[ok], n))
      } else {
        allocation <- expand.grid(m1 = grid, m2 = grid)
        allocation$m3 <- mean_sum - allocation$m1 - allocation$m2
        allocation <- allocation[
          allocation$m3 >= item_l - 1e-12 & allocation$m3 <= item_u + 1e-12,
        ]
        d1 <- allocation$m1 - floor(allocation$m1)
        d2 <- allocation$m2 - floor(allocation$m2)
        brute <- min(
          d1 * (1 - d1) + d2 * (1 - d2) + item_var_min(allocation$m3, n)
        )
      }
      expect_equal(brute, floor_sum, tolerance = 1e-9)
    }
  }
})
