# Agreement test: certify_analytic() vs brimmest() over full reporting grids.
# Run from the package root: Rscript dev/analytic_certification/validate.R
suppressMessages(devtools::load_all(quiet = TRUE))
source(file.path("dev", "analytic_certification", "certify_analytic.R"))

check_design <- function(l, u, n, digits, sd_step_limit = NULL) {
  W <- u - l
  means <- round(seq(l, u, by = 10^-digits), digits)
  sd_max <- (W / 2) * sqrt(n / (n - 1)) * 1.05
  sds <- round(seq(0, sd_max, by = 10^-digits), digits)
  grid <- expand.grid(mean = means, sd = sds)
  if (!is.null(sd_step_limit)) grid <- grid[seq(1, nrow(grid), by = sd_step_limit), ]
  ref <- brimmest(l = l, u = u, n = n, mean = grid$mean, sd = grid$sd, digits = digits)
  new <- certify_analytic(l, u, n, grid$mean, grid$sd, digits)
  bad <- which(ref$possible != new$possible)
  cat(sprintf("design %d-%d n=%d d=%d: %d tuples, %d disagreements\n",
              l, u, n, digits, nrow(grid), length(bad)))
  if (length(bad)) print(cbind(grid[head(bad, 20), ], ref = ref$possible[head(bad, 20)],
                               new = new$possible[head(bad, 20)]))
  invisible(length(bad))
}

total <- 0
total <- total + check_design(1, 5, 9, 1)
total <- total + check_design(1, 7, 9, 1)
total <- total + check_design(1, 7, 11, 1)
total <- total + check_design(0, 4, 15, 2, sd_step_limit = 7)
cat(sprintf("TOTAL disagreements: %d\n", total))
