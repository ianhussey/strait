# Validate the approximate SD plausibility bounds against ground truth.
#
# Grid: m = 63 (BDI-II), mu in {0.1, 0.3, 0.5}, rho in {0.02, 0.05, 0.15},
# N in {10, 20, 50, 200}, tail p in {0.01, 0.001}.
# Ground truth: the exact DP distribution of S^2 where N <= 20; Monte Carlo
# with 10^6 replicates everywhere (sampled as multinomial counts, so the cost
# does not grow with N). Also: the conditional-on-mean distribution (exact,
# N <= 20), and how often a lower quantile falls below the binomial SD.
#
# Run from the package root:  Rscript dev/sd-plausibility-bounds/validate.R
# Writes dev/sd-plausibility-bounds/validation.rds (~10 minutes).

pkgload::load_all(quiet = TRUE)
set.seed(20260927)

m <- 63
grid <- expand.grid(mu = c(0.1, 0.3, 0.5), rho = c(0.02, 0.05, 0.15))
Ns <- c(10, 20, 50, 200)
ps <- c(0.01, 0.001)
probs <- c(ps, 1 - ps)          # lower tails, then upper tails
reps <- 1e6

mc_s <- function(pr, N, reps, chunk = 1e5) {
  x <- 0:m
  out <- numeric(0)
  for (i in seq_len(reps / chunk)) {
    cnt <- stats::rmultinom(chunk, N, pr)
    s <- as.numeric(crossprod(x, cnt))
    q <- as.numeric(crossprod(x^2, cnt))
    out <- c(out, sqrt(pmax(0, (q - s^2 / N) / (N - 1))))
  }
  out
}
# order-statistic 95% interval for an MC quantile, as a relative half-width
mc_rel_se <- function(s_sorted, pp) {
  n <- length(s_sorted)
  k <- n * pp; h <- 1.96 * sqrt(n * pp * (1 - pp))
  lo <- s_sorted[max(1, floor(k - h))]; hi <- s_sorted[min(n, ceiling(k + h))]
  (hi - lo) / 2 / s_sorted[round(k)]
}
exact_q <- function(d, pp) sqrt(d$s2[which(cumsum(d$prob) >= pp - 1e-12)[1]])

rows <- list()
for (g in seq_len(nrow(grid))) {
  mu <- grid$mu[g]; rho <- grid$rho[g]
  pr <- strait:::.bb_pmf(m, mu, rho)
  sigma_binom <- sqrt(m * mu * (1 - mu))
  t0 <- Sys.time()
  ex <- s2_exact(m, mu, rho, c(10, 20))
  exc <- lapply(c(10, 20), function(N) {
    s2_exact(m, mu, rho, N, sums = round(N * m * mu))[[1]]
  })
  names(exc) <- c("10", "20")
  message(sprintf("mu = %.1f rho = %.2f: exact in %.0f s", mu, rho,
                  as.numeric(Sys.time() - t0, units = "secs")))
  for (N in Ns) {
    s <- sort(mc_s(pr, N, reps))
    for (pp in probs) {
      approx <- vapply(c("chisq2", "chisq3", "cf"), function(meth) {
        q <- strait:::.sd_plaus_row(0, m, m * mu, N, min(pp, 1 - pp), rho, 0,
                                    meth, FALSE, NULL, 1)
        if (pp < 0.5) q[1] else q[2]
      }, numeric(1))
      cond <- strait:::.sd_plaus_row(0, m, m * mu, N, min(pp, 1 - pp), rho, 0,
                                     "chisq2", TRUE, NULL, 1)
      has_exact <- N <= 20
      rows[[length(rows) + 1]] <- data.frame(
        mu = mu, rho = rho, N = N, prob = pp,
        tail = if (pp < 0.5) "lower" else "upper",
        p = min(pp, 1 - pp),
        mc = unname(stats::quantile(s, pp, type = 1)),
        mc_rel_ci = mc_rel_se(s, pp),
        exact = if (has_exact) exact_q(ex[[as.character(N)]], pp) else NA,
        exact_cond = if (has_exact) exact_q(exc[[as.character(N)]], pp) else NA,
        chisq2 = approx[["chisq2"]], chisq3 = approx[["chisq3"]],
        cf = approx[["cf"]], chisq2_cond = cond[if (pp < 0.5) 1 else 2],
        sigma_binom = sigma_binom,
        p_below_binom_exact = if (has_exact) {
          d <- ex[[as.character(N)]]; sum(d$prob[sqrt(d$s2) < sigma_binom])
        } else NA,
        p_below_binom_mc = mean(s < sigma_binom)
      )
    }
  }
}
res <- do.call(rbind, rows)
res$truth <- ifelse(is.na(res$exact), res$mc, res$exact)
for (meth in c("chisq2", "chisq3", "cf")) {
  res[[paste0("rel_", meth)]] <- res[[meth]] / res$truth - 1
}
res$rel_chisq2_cond <- res$chisq2_cond / res$exact_cond - 1
res$rel_chisq2_vs_cond <- res$chisq2 / res$exact_cond - 1
res$rel_exact_vs_mc <- res$exact / res$mc - 1
saveRDS(res, "dev/sd-plausibility-bounds/validation.rds")
message("done")
