# Extracted from test-sd_region_plot.R:217

# prequel ----------------------------------------------------------------------
l <- 1
u <- 5
n <- 7
m <- seq(l, u, by = 0.002)
bes <- sqrt(n / (n - 1))

# test -------------------------------------------------------------------------
for (cfg in list(c(4, 1, 5), c(5, 0, 3), c(7, 0, 6), c(6, 0, 2))) {
    n <- cfg[1]; l <- cfg[2]; u <- cfg[3]
    d <- strait:::.attainable_lattice(l, u, n, 1)
    g <- multisets(l:u, n)
    b <- unique(round(cbind(rowMeans(g), apply(g, 1, stats::sd)), 10))
    b <- b[order(b[, 1], b[, 2]), , drop = FALSE]
    a <- unique(round(as.matrix(d[c("mean", "sd")]), 10))
    a <- a[order(a[, 1], a[, 2]), , drop = FALSE]
    expect_equal(unname(a), unname(b), tolerance = 1e-9)
  }
expect_equal(nrow(strait:::.attainable_lattice(0, 3, 7, 1)), 100L)
expect_equal(nrow(strait:::.attainable_lattice(0, 6, 23, 1)), 8634L)
d <- strait:::.attainable_lattice(0, 3, 6, 2)
expect_true(all(abs(d$mean * 6 * 2 - round(d$mean * 6 * 2)) < 1e-9))
expect_error(strait:::.attainable_lattice(1, 7, 61, 5), "too large")
expect_lt(strait:::.lattice_cells(0, 6, 300, 1), 2e7)
expect_gt((300 * 6 + 1) * (300 * 6^2 + 1), 2e7)
