# Validation 05: certify_analytic2() vs brimmest() on small designs,
# at d = 1, 2 and 3 reported decimals.
#
# d = 1: exhaustive over the full (mean, sd) reporting grid.
# d = 2, 3: the full top-edge strip (every grid mean, sd within 0.3 of sd_max)
#           plus random tuples anywhere -- top edge is where all phantoms live,
#           random controls cover the dense zone.
#
# Expected output: 0 disagreements.
# Run from the package root: Rscript dev/analytic_certification/validate_05_grids.R
suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

sd_max_of <- function(S, W, n) {
  a <- S %/% W; b <- S - a * W
  sqrt((a * W^2 + b^2 - S^2 / n) / (n - 1))
}

full_grid <- function(l, u, n, d) {
  W <- u - l
  means <- round(seq(l, u, by = 10^-d), d)
  sds <- round(seq(0, (W / 2) * sqrt(n / (n - 1)) * 1.05, by = 10^-d), d)
  expand.grid(mean = means, sd = sds)
}

strip_plus_random <- function(l, u, n, d, strip = 0.3, n_random = 1500, seed = 9) {
  set.seed(seed)
  W <- u - l
  step <- 10^-d
  rows <- list()
  for (m in round(seq(l, u, by = max(step, (u - l) / 2000)), d)) {
    S0 <- round((m - l) * n)
    if (S0 < 1 || S0 > n * W - 1) next
    sm <- sd_max_of(S0, W, n)
    ks <- seq(max(0, floor((sm - strip) / step)), floor(sm / step) + 1)
    if (length(ks) > 40) ks <- sort(sample(ks, 40))
    for (k in ks) rows[[length(rows) + 1]] <- c(m, round(k * step, d))
  }
  for (i in seq_len(n_random)) {
    m <- round(runif(1, l, u), d)
    S0 <- min(max(round((m - l) * n), 1), n * W - 1)
    rows[[length(rows) + 1]] <- c(m, round(runif(1, 0, sd_max_of(S0, W, n)), d))
  }
  g <- do.call(rbind, rows)
  data.frame(mean = g[, 1], sd = g[, 2])
}

check <- function(g, l, u, n, d) {
  ref <- brimmest(l, u, n, mean = g$mean, sd = g$sd, digits = d)
  new <- certify_analytic2(l, u, n, g$mean, g$sd, d)
  dis <- sum(ref$possible != new$possible)
  cat(sprintf("[05] %d-%d n=%d d=%d: %d tuples, %d impossible, %d disagreements\n",
              l, u, n, d, nrow(g), sum(!new$possible), dis))
  if (dis > 0) {
    bad <- which(ref$possible != new$possible)
    print(cbind(g[head(bad, 12), ], brimmest = ref$possible[head(bad, 12)],
                analytic = new$possible[head(bad, 12)]))
  }
  dis
}

tot <- 0
tot <- tot + check(full_grid(1, 5, 9, 1), 1, 5, 9, 1)
tot <- tot + check(full_grid(1, 7, 9, 1), 1, 7, 9, 1)
tot <- tot + check(full_grid(1, 7, 11, 1), 1, 7, 11, 1)
tot <- tot + check(strip_plus_random(1, 5, 9, 2), 1, 5, 9, 2)
tot <- tot + check(strip_plus_random(1, 5, 9, 3), 1, 5, 9, 3)
tot <- tot + check(strip_plus_random(1, 7, 25, 2), 1, 7, 25, 2)
tot <- tot + check(strip_plus_random(1, 7, 25, 3), 1, 7, 25, 3)
tot <- tot + check(strip_plus_random(0, 10, 15, 3), 0, 10, 15, 3)
cat(sprintf("[05] TOTAL disagreements: %d\n", tot))
quit(status = as.integer(tot > 0))
