suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))
set.seed(7)

probe <- function(l, u, n, d, n_tuples = 120, label = "") {
  W <- u - l
  sd_max_of <- function(m) {
    S <- round((m - l) * n); a <- S %/% W; b <- S - a * W
    sqrt((a * W^2 + b^2 - S^2/n) / (n - 1))
  }
  means <- round(runif(n_tuples, l + 0.02*(u-l), u - 0.02*(u-l)), d)
  sds <- vapply(means, function(m) round(max(0, sd_max_of(m) - runif(1,0,2)*10^-d), d), 0)
  cat(sprintf("\n=== %s: %d-%d n=%d d=%d, %d edge tuples ===\n", label, l, u, n, d, n_tuples))
  # how often does brimmest's fast path fail to settle (-> corridor DP)?
  na_ct <- 0; t_b <- numeric(n_tuples); t_a <- numeric(n_tuples); agree <- logical(n_tuples)
  mg <- 1L
  for (i in seq_len(n_tuples)) {
    for (rr in c("half_up","half_down")) {
      tg <- strait:::.target_states(l, u, n, mg, means[i], sds[i], d, d, rr)
      got <- strait:::.certify_fast(W, n, tg, budget = 2e5)
      if (is.na(got)) { na_ct <- na_ct + 1; break }
    }
    t0 <- Sys.time(); rb <- brimmest(l, u, n, mean = means[i], sd = sds[i], digits = d)
    t_b[i] <- as.numeric(Sys.time()-t0, units="secs")
    t0 <- Sys.time(); ra <- certify_analytic(l, u, n, means[i], sds[i], d)
    t_a[i] <- as.numeric(Sys.time()-t0, units="secs")
    agree[i] <- rb$possible == ra$possible
  }
  cat(sprintf("  corridor-DP fallback needed: %d/%d tuples\n", na_ct, n_tuples))
  cat(sprintf("  per-tuple  brimmest: median %.1f ms, max %.0f ms | analytic: median %.1f ms, max %.0f ms\n",
              1000*median(t_b), 1000*max(t_b), 1000*median(t_a), 1000*max(t_a)))
  cat(sprintf("  disagreements: %d\n", sum(!agree)))
  wb <- which.max(t_b); wa <- which.max(t_a)
  cat(sprintf("  slowest for brimmest: mean=%.2f sd=%.2f (%.0f ms; analytic %.0f ms)\n",
              means[wb], sds[wb], 1000*t_b[wb], 1000*t_a[wb]))
  cat(sprintf("  slowest for analytic: mean=%.2f sd=%.2f (%.0f ms; brimmest %.0f ms)\n",
              means[wa], sds[wa], 1000*t_a[wa], 1000*t_b[wa]))
  invisible(sum(!agree))
}

tot <- 0
tot <- tot + probe(0, 63, 50, 2, label = "wide scale")
tot <- tot + probe(0, 63, 200, 2, label = "wide scale, larger n")
tot <- tot + probe(1, 7, 2000, 2, label = "huge n")
cat(sprintf("\nTOTAL disagreements: %d\n", tot))
