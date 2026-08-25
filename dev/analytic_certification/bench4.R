# Wide-scale sum-score designs (e.g. 50-item 1-7 scale summed: l=50, u=350),
# awkward n. Oversample where phantoms should live (top edge, corners,
# integer-mean streaks) and where they should not (interior, bottom).
suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

run <- function(l, u, n, d, N_edge = 120, N_corner = 60, N_int = 80, N_bottom = 40) {
  W <- u - l
  set.seed(101)
  sd_max_of <- function(S) { a <- S %/% W; b <- S - a*W; sqrt((a*W^2 + b^2 - S^2/n)/(n-1)) }
  mk <- function(S) round(l + S/n, d)
  tuples <- list()
  # top edge across the full mean range (incl. exact-integer-mean streaks b=0/1)
  Ss <- c(sample(seq(2, n*W-2), N_edge/2),
          sapply(sample(seq(1, n-1), N_edge/4, replace=TRUE), function(k) k*W + sample(0:1,1)),
          round(seq(0.02, 0.98, length.out = N_edge/4) * n*W))
  for (S in Ss) {
    sm <- sd_max_of(S)
    tuples[[length(tuples)+1]] <- c(mk(S), round(max(0, sm - runif(1, 0, 2)*10^-d), d))
  }
  # corners: mean within ~1.5 scale units of l or u, sd near max
  for (i in seq_len(N_corner)) {
    S <- sample(2:round(1.5*n), 1)
    if (runif(1) < .5) S <- n*W - S
    sm <- sd_max_of(S)
    tuples[[length(tuples)+1]] <- c(mk(S), round(max(0, sm - runif(1, 0, 3)*10^-d), d))
  }
  # interior controls (should never be phantoms)
  for (i in seq_len(N_int)) {
    S <- sample(seq(round(0.15*n*W), round(0.85*n*W)), 1)
    sm <- sd_max_of(S)
    tuples[[length(tuples)+1]] <- c(mk(S), round(runif(1, 0.15, 0.8) * sm, d))
  }
  # bottom controls (sd near minimum)
  for (i in seq_len(N_bottom)) {
    S <- sample(seq(round(0.1*n*W), round(0.9*n*W)), 1)
    tuples[[length(tuples)+1]] <- c(mk(S), round(runif(1, 0.02, 0.1) * sd_max_of(S), d))
  }
  tp <- do.call(rbind, tuples)
  region <- c(rep("edge", N_edge), rep("corner", N_corner),
              rep("interior", N_int), rep("bottom", N_bottom))
  t_b <- t_a <- rep(NA_real_, nrow(tp)); rb <- ra <- rep(NA, nrow(tp)); err <- character(nrow(tp))
  for (i in seq_len(nrow(tp))) {
    t0 <- Sys.time()
    r <- tryCatch(brimmest(l, u, n, mean = tp[i,1], sd = tp[i,2], digits = d)$possible,
                  error = function(e) { err[i] <<- conditionMessage(e); NA })
    t_b[i] <- as.numeric(Sys.time()-t0, units="secs"); rb[i] <- r
    t0 <- Sys.time()
    ra[i] <- certify_analytic2(l, u, n, tp[i,1], tp[i,2], d)$possible
    t_a[i] <- as.numeric(Sys.time()-t0, units="secs")
  }
  cat(sprintf("\n=== %d-%d (W=%d), n=%d, d=%d: %d tuples ===\n", l, u, u-l, n, d, nrow(tp)))
  cmp <- !is.na(rb)
  cat(sprintf("  brimmest errors/NA: %d | comparable: %d | disagreements: %d\n",
              sum(!cmp), sum(cmp), sum(rb[cmp] != ra[cmp])))
  if (any(cmp & rb != ra)) {
    bad <- which(cmp & rb != ra)
    print(cbind(as.data.frame(tp[head(bad,10),,drop=FALSE]), region=region[head(bad,10)],
                brimmest=rb[head(bad,10)], analytic=ra[head(bad,10)]))
  }
  for (rg in unique(region)) {
    ix <- region == rg
    cat(sprintf("  %-8s: %3d tuples, %3d impossible | brimmest med %6.1f max %7.0f ms | v2 med %5.2f max %6.1f ms\n",
        rg, sum(ix), sum(!ra[ix]), 1000*median(t_b[ix], na.rm=TRUE), 1000*max(t_b[ix], na.rm=TRUE),
        1000*median(t_a[ix]), 1000*max(t_a[ix])))
  }
  invisible(sum(cmp & rb != ra))
}

tot <- 0
tot <- tot + run(50, 350, 37, 2)
tot <- tot + run(50, 350, 83, 2)
tot <- tot + run(50, 350, 83, 1)
tot <- tot + run(13, 530, 61, 2)   # deliberately weird bounds and n
cat(sprintf("\nTOTAL disagreements: %d\n", tot))
