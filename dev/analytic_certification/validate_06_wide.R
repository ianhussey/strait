# Validation 06: certify_analytic2() vs brimmest() on wide sum-score designs
# with awkward n, at d = 1, 2 and 3 -- region-stratified tuples:
#   edge     sd hugging sd_max across the whole mean range (incl. b = 0/1 streaks)
#   corner   mean within ~1.5 scale units of l or u, sd near max
#   interior mid-range mean, mid-range sd    (structure theorem: never phantoms)
#   bottom   sd near the minimum             (structure theorem: never phantoms)
#
# Expected output: 0 disagreements, 0 brimmest errors.
# Run from the package root: Rscript dev/analytic_certification/validate_06_wide.R
suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

run <- function(l, u, n, d, N_edge = 120, N_corner = 60, N_int = 80, N_bottom = 40,
                seed = 101) {
  W <- u - l
  set.seed(seed)
  sd_max_of <- function(S) { a <- S %/% W; b <- S - a * W
                             sqrt((a * W^2 + b^2 - S^2 / n) / (n - 1)) }
  mk <- function(S) round(l + S / n, d)
  tuples <- list()
  Ss <- c(sample(seq(2, n * W - 2), N_edge / 2),
          sapply(sample(seq(1, n - 1), N_edge / 4, replace = TRUE),
                 function(k) k * W + sample(0:1, 1)),
          round(seq(0.02, 0.98, length.out = N_edge / 4) * n * W))
  for (S in Ss)
    tuples[[length(tuples) + 1]] <-
      c(mk(S), round(max(0, sd_max_of(S) - runif(1, 0, 2) * 10^-d), d))
  for (i in seq_len(N_corner)) {
    S <- sample(2:round(1.5 * n), 1)
    if (runif(1) < .5) S <- n * W - S
    tuples[[length(tuples) + 1]] <-
      c(mk(S), round(max(0, sd_max_of(S) - runif(1, 0, 3) * 10^-d), d))
  }
  for (i in seq_len(N_int)) {
    S <- sample(seq(round(0.15 * n * W), round(0.85 * n * W)), 1)
    tuples[[length(tuples) + 1]] <- c(mk(S), round(runif(1, .15, .8) * sd_max_of(S), d))
  }
  for (i in seq_len(N_bottom)) {
    S <- sample(seq(round(0.1 * n * W), round(0.9 * n * W)), 1)
    tuples[[length(tuples) + 1]] <- c(mk(S), round(runif(1, .02, .1) * sd_max_of(S), d))
  }
  tp <- do.call(rbind, tuples)
  region <- c(rep("edge", N_edge), rep("corner", N_corner),
              rep("interior", N_int), rep("bottom", N_bottom))
  rb <- ra <- rep(NA, nrow(tp)); nerr <- 0
  for (i in seq_len(nrow(tp))) {
    rb[i] <- tryCatch(
      brimmest(l, u, n, mean = tp[i, 1], sd = tp[i, 2], digits = d)$possible,
      error = function(e) { nerr <<- nerr + 1; NA })
    ra[i] <- certify_analytic2(l, u, n, tp[i, 1], tp[i, 2], d)$possible
  }
  cmp <- !is.na(rb)
  dis <- sum(rb[cmp] != ra[cmp])
  cat(sprintf("[06] %d-%d (W=%d) n=%d d=%d: %d tuples | brimmest errors %d | disagreements %d\n",
              l, u, W, n, d, nrow(tp), nerr, dis))
  for (rg in unique(region)) {
    ix <- region == rg
    cat(sprintf("     %-8s: %3d tuples, %3d impossible\n", rg, sum(ix), sum(!ra[ix])))
  }
  if (dis > 0) {
    bad <- which(cmp & rb != ra)
    print(cbind(as.data.frame(tp[head(bad, 12), , drop = FALSE]),
                region = region[head(bad, 12)],
                brimmest = rb[head(bad, 12)], analytic = ra[head(bad, 12)]))
  }
  dis
}

tot <- 0
tot <- tot + run(50, 350, 37, 2)
tot <- tot + run(50, 350, 37, 3)
tot <- tot + run(50, 350, 83, 1)
tot <- tot + run(50, 350, 83, 2)
tot <- tot + run(50, 350, 83, 3)
tot <- tot + run(13, 530, 61, 2)
tot <- tot + run(13, 530, 61, 3)
cat(sprintf("[06] TOTAL disagreements: %d\n", tot))
quit(status = as.integer(tot > 0))
