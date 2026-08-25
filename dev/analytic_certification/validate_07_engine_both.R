# Validation 07: brimmest(engine = "both") sweep after integration.
# Every tuple is certified independently by the analytic engine AND the
# search/lattice machinery inside brimmest() itself, which errors on any
# disagreement -- so a silent completion IS the validation. Covers d = 1/2/3,
# small + wide designs, phantom-rich top strips, and the rounding rules
# half_up, native (bankers) and trunc.
#
# Run from the package root: Rscript dev/analytic_certification/validate_07_engine_both.R
suppressMessages(devtools::load_all(quiet = TRUE))
rules <- c("half_up", "native", "trunc")

sweep <- function(l, u, n, d, N = 400, seed = 13) {
  set.seed(seed)
  W <- u - l
  sd_max_of <- function(S) { a <- S %/% W; b <- S - a * W
                             sqrt((a * W^2 + b^2 - S^2 / n) / (n - 1)) }
  S <- sample(seq(2, n * W - 2), N, replace = TRUE)
  means <- round(l + S / n, d)
  # half the tuples hug the top edge (phantom-rich), half anywhere
  sds <- ifelse(seq_len(N) %% 2 == 0,
                pmax(0, vapply(S, sd_max_of, 0) - runif(N, 0, 2) * 10^-d),
                runif(N, 0, vapply(S, sd_max_of, 0)))
  sds <- round(sds, d)
  r <- brimmest(l, u, n, mean = means, sd = sds, digits = d,
                rounding = rules, engine = "both")
  cat(sprintf("[07] %d-%d n=%d d=%d: %d tuples x %d rules cross-checked, %d impossible -- OK\n",
              l, u, n, d, N, length(rules), sum(!r$possible)))
}

sweep(1, 5, 9, 1)
sweep(1, 5, 9, 3)
sweep(1, 7, 25, 2)
sweep(1, 7, 25, 3)
sweep(0, 10, 15, 3)
sweep(0, 63, 50, 2)      # targeted route: analytic vs search
sweep(0, 63, 50, 3)
sweep(50, 350, 83, 2)    # wide sum-score scale
sweep(50, 350, 83, 3)
sweep(13, 530, 61, 3)
cat("[07] all sweeps completed with zero engine disagreements\n")
