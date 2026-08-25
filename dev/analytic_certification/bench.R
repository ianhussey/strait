# Benchmark: certify_analytic() (band recursion, plain R prototype) vs brimmest()
# (lattice / targeted routes). Correctness agreement + wall-clock.
suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

mk_tuples <- function(l, u, n, d, n_edge = 100, n_mid = 100, seed = 1) {
  set.seed(seed)
  W <- u - l
  means <- round(runif(n_edge + n_mid, l, u), d)
  sd_max_of <- function(m) {
    S <- round((m - l) * n); a <- S %/% W; b <- S - a * W
    sqrt((a * W^2 + b^2 - S^2 / n) / (n - 1))
  }
  sds <- numeric(n_edge + n_mid)
  for (i in seq_len(n_edge)) {          # hug the top edge (hard cases, phantoms)
    sm <- sd_max_of(means[i])
    sds[i] <- round(max(0, sm - runif(1, 0, 3) * 10^-d), d)
  }
  for (i in (n_edge + 1):(n_edge + n_mid)) {   # interior (easy cases)
    sm <- sd_max_of(means[i])
    sds[i] <- round(runif(1, 0.1, max(0.11, sm * 0.85)), d)
  }
  data.frame(mean = means, sd = sds)
}

bench_design <- function(l, u, n, d, n_edge = 100, n_mid = 100) {
  tp <- mk_tuples(l, u, n, d, n_edge, n_mid)
  cat(sprintf("\n=== design %d-%d, n=%d, digits=%d (%d tuples: %d edge, %d interior) ===\n",
              l, u, n, d, nrow(tp), n_edge, n_mid))
  # brimmest: cold single, then batch
  t0 <- Sys.time()
  r1 <- brimmest(l, u, n, mean = tp$mean[1], sd = tp$sd[1], digits = d)
  t_cold <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  rb <- brimmest(l, u, n, mean = tp$mean, sd = tp$sd, digits = d)
  t_batch <- as.numeric(Sys.time() - t0, units = "secs")
  # analytic: single, then batch (no cross-call cache at all)
  t0 <- Sys.time()
  a1 <- certify_analytic(l, u, n, tp$mean[1], tp$sd[1], d)
  ta_cold <- as.numeric(Sys.time() - t0, units = "secs")
  t0 <- Sys.time()
  ab <- certify_analytic(l, u, n, tp$mean, tp$sd, d)
  ta_batch <- as.numeric(Sys.time() - t0, units = "secs")
  dis <- sum(rb$possible != ab$possible)
  cat(sprintf("  brimmest : first tuple %7.3fs | batch of %d %8.3fs (%.1f ms/tuple warm)\n",
              t_cold, nrow(tp), t_batch, 1000 * t_batch / nrow(tp)))
  cat(sprintf("  analytic : first tuple %7.3fs | batch of %d %8.3fs (%.1f ms/tuple)\n",
              ta_cold, nrow(tp), ta_batch, 1000 * ta_batch / nrow(tp)))
  cat(sprintf("  verdicts : %d impossible of %d | disagreements: %d\n",
              sum(!ab$possible), nrow(tp), dis))
  if (dis > 0) {
    bad <- which(rb$possible != ab$possible)
    print(cbind(tp[head(bad, 12), ], brimmest = rb$possible[head(bad, 12)],
                analytic = ab$possible[head(bad, 12)]))
  }
  invisible(dis)
}

tot <- 0
tot <- tot + bench_design(1, 5, 9, 1)          # tiny: lattice route
tot <- tot + bench_design(1, 7, 25, 2)         # small: lattice route
tot <- tot + bench_design(1, 9, 100, 2)        # medium-n
tot <- tot + bench_design(0, 63, 50, 2)        # wide scale: targeted route (doc example)
cat(sprintf("\nTOTAL disagreements: %d\n", tot))
