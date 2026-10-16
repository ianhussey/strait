# Quasi-integer composite bounds, observed (attained) extremes, AND plausibility
#
# Copy of app_quasi_integer_alpha_bounds_and_observed_range.R that adds
# PLAUSIBILITY bounds for the SD (strait::sd_plausibility(); method and
# validation in paper/sd-plausibility-bounds/). Everything below the next
# paragraph is unchanged from that app, except where marked PLAUSIBILITY.
#
# POSSIBLE vs PLAUSIBLE. The envelope says which SDs are arithmetically
# possible. The plausibility bands say which SDs a realistic data-generating
# process would produce. The sum score is modelled as BetaBinomial(m, mu, rho)
# on [k*min, k*max], m = k*(max - min). The user supplies an assumed POPULATION
# mean and SD, which fix rho by the method of moments, and tau, the
# between-study SD of rho. The bands hold rho fixed and let the mean vary along
# the x axis: at each mean they are quantiles of the sample SD among samples of
# size n whose sample mean is that mean (condition_on_mean = TRUE), mixed over
# rho ~ N(rho, tau^2) truncated to [0, 1). They are clipped to the possibility
# envelope. The reported mean is compared with the population mean through the
# sampling range of a mean of n draws, M0 +/- z * SD0 / sqrt(n) (no
# between-study heterogeneity in the mean). The beta-binomial models the sum
# score directly: it does not use alpha, the item structure or the observed
# extremes, which constrain only the possibility envelope.
#
# The curves use the conditional two-moment approximation. At the reported
# mean, the verdict uses the EXACT conditional distribution (dynamic
# programming over the sum and sum of squares) when that is cheap enough,
# because the approximation is poor in the lower tail at small n.
#
# Needs this repository's strait (loaded from source, two levels up) or an
# installed strait >= 0.4.13.
#
# ---- original header ----
#
# Copy of app_quasi_integer_alpha_bounds.R that adds the option to condition on
# the OBSERVED minimum and maximum of the composite, in addition to the logical
# item limits. Companion to master.qmd (@sec-alpha-integer, Appendix H) and to
# the working notes in bounds_observed_extremes.qmd.
#
# THE DISTINCTION. The logical limits l, u say only that every observation lies
# somewhere in [l, u]; nothing need sit at either end, so the floor is driven
# purely by integer granularity and hugs zero. Reported observed extremes A, B
# say more: at least one observation IS A and at least one IS B (constraint E).
# The sample is then forced to span B - A, which creates a substantial floor,
# and it restricts which means are reachable at all.
#
# WHAT E BUYS (sum-score units; p = S_bar - A, q = B - S_bar, R = B - A):
#
#   feasible mean band   A + R/n <= S_bar <= B - R/n
#                        outside it, no sample with both extremes attained has
#                        that mean, at any SD -> the mean itself is impossible
#
#   floor, continuous    SS = p^2 + q^2 + (p-q)^2/(n-2)
#                           = R^2/2 + (2n/(n-2))(S_bar - (A+B)/2)^2
#                        a parabola whose vertex is the Nagy (1918) / Thomson
#                        (1955) value R^2/2; strictly sharper than the
#                        Laguerre-Samuelson corollary max(p,q)^2
#
#   floor, integer       + (n-2) d(1-d),  d = frac(m), m = (n S_bar - A - B)/(n-2)
#   floor, quasi-integer - g(1-g),        g = frac(n S_bar)   (same GRIM quantity
#                        as master, since A + B is an integer)
#
#   ceiling              UNCHANGED by attainment inside the band: the attaining
#                        configuration of the mean-conditional maximum already
#                        occupies both extremes. So the ceiling gain comes only
#                        from the narrower range, via Structure S on [A, B].
#
# COMBINING WITH ALPHA. The alpha-amplified floor and the attained floor are
# different mechanisms (item-variance granularity vs. forced span of the sum
# score), and the alpha-ceiling and the observed-range ceiling likewise. Each is
# separately a valid bound, so we take
#     floor   = max(alpha-amplified quasi floor, attained floor)
#     ceiling = min(alpha-ceiling, observed-range Structure S ceiling)
# which is valid but not guaranteed sharp: a configuration attaining one need
# not attain the other. Curves are labelled accordingly in the app.
#
# All formulas for the attained floor were brute-force validated in
# bounds_observed_extremes.qmd (exhaustive enumeration + constrained optimization).
#
# Run with:  shiny::runApp("app_quasi_integer_alpha_bounds_observed_range_and_plausibility.R")
# (needs shiny + ggplot2, and pkgload to load strait from source)

library(shiny)
library(ggplot2)

# PLAUSIBILITY: strait from this repository when run from inside it, without
# attaching (the app defines its own bound primitives, some with names strait
# also exports); otherwise an installed strait.
local({
  # walk up from the working directory to the strait package root
  is_root <- function(d) {
    f <- file.path(d, "DESCRIPTION")
    file.exists(f) && identical(readLines(f, 1), "Package: strait")
  }
  root <- getwd()
  while (!is_root(root) && dirname(root) != root) root <- dirname(root)
  if (is_root(root) && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, attach = FALSE, quiet = TRUE)
  } else {
    requireNamespace("strait")
  }
  if (!"sd_plausibility" %in% getNamespaceExports("strait")) {
    stop("This app needs strait >= 0.4.13 (sd_plausibility()).")
  }
})

INK  <- "#1c1917"; MUTE <- "#78716c"; GREY <- "#57534e"
BLUE <- "#1d4ed8"; RED  <- "#dc2626"; TEAL <- "#0d9488"; AMBER <- "#d97706"
VIO  <- "#7c3aed"
TOL  <- 1e-9
PLAUS <- "#0e7490"; PLAUS_FILL <- "#cffafe"   # PLAUSIBILITY bands
EXACT_BUDGET <- 3e6   # max DP cells x rho nodes for the exact check

# ---- bound primitives (population variance of one item / of the sum score) ----

# per-item quasi-integer MAX population variance summed over k items (Theorem H3)
v_max_alpha <- function(muS, k, n, a, b) {
  h  <- (b - a) / n
  T  <- (muS - k * a) / h
  x  <- T / k
  j  <- pmin(floor(x + 1e-12), n - 1)
  th <- x - j
  phi <- T - floor(T + 1e-12)
  phi <- ifelse(phi > 1 - 1e-9, 0, phi)
  A <- function(jj) h^2 * jj * (n - jj)
  pmax(0, k * ((1 - th) * A(j) + th * A(j + 1)) - h^2 * (n - 1) * phi * (1 - phi))
}

# quasi-integer MIN sample SD of the sum score (sec-quasi-floor), range-free
sd_min_quasi <- function(muS, n) {
  d <- muS - floor(muS)
  g <- n * muS - floor(n * muS)
  sqrt(pmax(0, (n * d * (1 - d) - g * (1 - g)) / (n - 1)))
}

# ---- observed-extremes primitives (bounds_observed_extremes.qmd) -------------

# Structure S maximum applied to the sum score on the OBSERVED range [A, B].
# Unchanged by attainment inside the feasible band (the attaining configuration
# already places observations at both extremes).
sd_max_quasi_obs <- function(muS, n, A, B) {
  R <- B - A
  S_excess <- n * (muS - A)
  n_u <- pmin(floor(S_excess / R + 1e-9), n - 1)
  x_r <- A + (S_excess - n_u * R)
  n_l <- n - n_u - 1
  sqrt(pmax(0, (n_l * A^2 + x_r^2 + n_u * B^2 - n * muS^2) / (n - 1)))
}

# Minimum sample SD with both observed extremes attained. quasi = TRUE applies
# the Z_{n-1} relaxation to the interior so the curve is defined at every mean.
sd_min_attained <- function(muS, n, A, B, quasi = TRUE) {
  R <- B - A
  if (n <= 2) {
    return(ifelse(abs(muS - (A + B) / 2) <= 1e-9, sqrt(R^2 / 2 / (n - 1)), NA_real_))
  }
  k2 <- n - 2
  p <- muS - A
  q <- B - muS
  cont <- p^2 + q^2 + (p - q)^2 / k2          # continuous attained floor
  T <- n * muS - A - B                        # interior sum
  m <- T / k2                                 # interior mean
  d <- m - floor(m)
  g <- T - floor(T)                           # = frac(n*muS) when A + B integer
  gran <- k2 * d * (1 - d) - if (quasi) g * (1 - g) else 0
  sqrt(pmax(0, (cont + pmax(0, gran)) / (n - 1)))
}

# ---- PLAUSIBILITY primitives ----------------------------------------------

# Method-of-moments rho from a population mean and SD (sum-score units).
rho_mom <- function(M0, SD0, L, U) {
  m <- U - L; mu <- (M0 - L) / m
  (SD0^2 / (m * mu * (1 - mu)) - 1) / (m - 1)
}

# Number of DP cells the exact method needs at sample size n on 0..m.
exact_cells <- function(n, m) {
  s <- 0:(n * m)
  a <- s %/% n; r <- s %% n
  qmin <- (n - r) * a^2 + r * (a + 1)^2
  qmax <- (s %/% m) * m^2 + (s %% m)^2
  sum(qmax - qmin + 1)
}

# Exact conditional plausibility bounds at one reported mean (sum-score
# units): the distribution of S over samples whose mean rounds to muS within
# +/- h, mixed over rho nodes. Returns lower/upper for each p, or NULL.
exact_bounds <- function(L, U, muS, n, rho, tau, h, ps) {
  m <- U - L; mu <- (muS - L) / m
  t <- ceiling(n * (muS - L - h) - 1e-9):floor(n * (muS - L + h) + 1e-9)
  t <- t[t >= 0 & t <= n * m]
  if (length(t) == 0) return(NULL)
  K <- if (tau > 0) 9 else 1
  if (exact_cells(n, m) * K > EXACT_BUDGET) return(NULL)
  rhos <- strait:::.rho_nodes(max(0, rho), tau, K)
  d <- do.call(rbind, lapply(rhos, function(r) {
    strait::s2_exact(m, mu, min(r, 1), n, sums = t)[[1]]
  }))
  agg <- rowsum(d$prob / length(rhos), round(d$s2 * n * (n - 1)), reorder = TRUE)
  s <- sqrt(as.numeric(rownames(agg)) / (n * (n - 1)))
  cdf <- cumsum(agg[, 1])
  q <- function(pp) s[which(cdf >= pp - 1e-12)[1]]
  data.frame(p = ps, lower = vapply(ps, q, 0), upper = vapply(1 - ps, q, 0))
}

ui <- fluidPage(
  tags$head(tags$style(HTML(sprintf("
    body { background:#fafaf9; color:%s; font-family: ui-sans-serif, system-ui, sans-serif; }
    .title-x { font-family: Georgia, 'Times New Roman', serif; font-size:22px; letter-spacing:-0.01em; }
    .lede { font-size:13px; color:#57534e; max-width:860px; }
    .stat-card { border:1px solid #e7e5e4; background:#fff; border-radius:8px; padding:8px 10px; }
    .stat-label { font-size:11px; color:%s; } .stat-val { font-size:17px; font-family: ui-monospace, Menlo, monospace; }
    .verdict { font-size:15px; font-weight:600; padding:8px 0; }
    .form-group { margin-bottom:8px; }
  ", INK, MUTE)))),

  div(class = "title-x", HTML("Composite SD bounds: possible (logical limits, observed extremes) and plausible")),
  div(class = "lede", style = "margin:4px 0 12px;",
      HTML("The possibility envelope for a composite of <i>k</i> integer items. <b>Logical limits alone</b>
            (faint curves) place every observation in [min, max] but require nothing at either end, so the
            floor is driven only by integer granularity and hugs zero. <b>Observed extremes</b> (bold curves)
            additionally require that at least one observation equals the reported minimum and one the
            reported maximum &mdash; which forces the sample to span that range, lifts the floor sharply, and
            restricts which means are reachable at all. The x axis always spans the <i>logical</i> range so
            the two regimes can be compared directly. <b>Plausibility bands</b> (teal) add a different
            question: not whether an SD is possible, but whether a realistic population with the assumed
            mean, SD and between-study heterogeneity would produce it in a sample of this size.")),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      numericInput("lo",    "logical item min", 1, step = 1),
      numericInput("hi",    "logical item max", 7, step = 1),
      numericInput("k",     "k (items)", 3, min = 1, step = 1),
      numericInput("n",     "n (sample size)", 5, min = 2, step = 1),
      numericInput("mean",  "reported mean", 12, step = 0.1),
      numericInput("sd",    "reported SD", 3.0, min = 0, step = 0.01),
      numericInput("alpha", HTML("Cronbach&rsquo;s &alpha;"), 0.80, min = -0.99, max = 0.99, step = 0.1),
      radioButtons("scale", "reported mean/SD are:",
                   c("sum score" = "sum", "mean score (item average)" = "mean"),
                   selected = "sum"),
      tags$hr(style = "margin:10px 0;"),
      checkboxInput("use_obs", HTML("<b>also condition on observed min/max</b>"), FALSE),
      conditionalPanel(
        "input.use_obs == true",
        numericInput("obs_lo", "observed minimum", 6, step = 1),
        numericInput("obs_hi", "observed maximum", 18, step = 1),
        div(class = "lede", style = "font-size:11px; margin-top:4px;",
            HTML("In the <b>same units</b> as the reported mean/SD above, and taken to be
                  <i>attained</i>: at least one observation sits at each. Requires
                  <i>n</i>&nbsp;&ge;&nbsp;2."))
      ),
      tags$hr(style = "margin:10px 0;"),
      checkboxInput("use_plaus", HTML("<b>show plausibility bounds</b>"), TRUE),
      conditionalPanel(
        "input.use_plaus == true",
        numericInput("pop_mean", "assumed population mean", 12, step = 0.1),
        numericInput("pop_sd",   "assumed population SD", 3.5, min = 0, step = 0.01),
        numericInput("tau", HTML("&tau; (between-study SD of &rho;)"), 0.02,
                     min = 0, step = 0.005),
        numericInput("mean_dp", "reported mean: decimal places", 2, min = 0,
                     max = 6, step = 1),
        checkboxInput("use_exact",
                      "exact distribution at the reported mean (when feasible)", TRUE),
        div(class = "lede", style = "font-size:11px; margin-top:4px;",
            HTML("Population mean and SD are in the <b>same units</b> as the reported values.
                  They fix the beta-binomial dispersion &rho; of the sum score; &tau; spreads
                  &rho; across studies. Bands: 1% (dark) and 0.1% (light) in each tail."))
      ),
      div(class = "lede", style = "font-size:11px; margin-top:8px;",
          HTML("With <i>k</i> items each on [min, max], the sum score ranges
                [<i>k</i>&middot;min, <i>k</i>&middot;max]; the mean score is the sum &divide; <i>k</i>."))
    ),
    mainPanel(
      width = 9,
      uiOutput("verdict"),
      uiOutput("plaus_verdict"),
      uiOutput("stats"),
      uiOutput("plaus_stats"),
      plotOutput("envelope_plot", height = "460px"),
      uiOutput("obs_note"),
      uiOutput("plaus_note")
    )
  )
)

server <- function(input, output, session) {

  base <- reactive({
    lo <- input$lo; hi <- input$hi; k <- input$k; n <- input$n
    a <- input$alpha; meanscore <- input$scale == "mean"
    ok <- all(vapply(list(lo, hi, k, n, input$mean, input$sd, a), is.finite, logical(1))) &&
          hi > lo && k >= 1 && n >= 2 && a > -1 && a < 1 && input$sd >= 0
    validate(need(ok, "check inputs: need max > min, k >= 1, n >= 2, -1 < alpha < 1, SD >= 0"))

    conv <- if (meanscore) k else 1        # sum-score = reported * conv
    muS  <- input$mean * conv
    sS   <- input$sd   * conv
    L <- k * lo; U <- k * hi
    validate(need(muS > L + TOL & muS < U - TOL,
                  sprintf("mean must lie strictly inside the %s range [%g, %g]",
                          if (meanscore) "mean-score" else "sum-score",
                          if (meanscore) lo else L, if (meanscore) hi else U)))
    cc  <- (k - 1) / k                     # c = (k-1)/k; k = 1 gives cc = 0
    denom <- 1 - cc * a
    validate(need(denom > TOL, "alpha too high for this k: 1 - c*alpha must be positive"))

    # observed extremes, converted to sum-score units
    use_obs <- isTRUE(input$use_obs)
    A <- B <- NA_real_
    if (use_obs) {
      validate(need(all(vapply(list(input$obs_lo, input$obs_hi), is.finite, logical(1))),
                    "enter both observed minimum and maximum"))
      A <- input$obs_lo * conv
      B <- input$obs_hi * conv
      validate(need(B > A + TOL, "observed maximum must exceed observed minimum"))
      validate(need(A >= L - TOL & B <= U + TOL,
                    "observed extremes must lie within the logical range"))
    }

    list(lo = lo, hi = hi, k = k, n = n, alpha = a, muS = muS, sS = sS,
         L = L, U = U, cc = cc, denom = denom, meanscore = meanscore,
         use_obs = use_obs, A = A, B = B,
         disp = if (meanscore) 1 / k else 1,   # multiply sum-score values by this for display
         bes = sqrt(n / (n - 1)))
  })

  # feasible mean band under attainment (sum-score units); NULL if not in use
  band <- reactive({
    b <- base()
    if (!b$use_obs) return(NULL)
    R <- b$B - b$A
    c(lo = b$A + R / b$n, hi = b$B - R / b$n)
  })

  # all curves, in sum-score units
  bounds_at <- function(b, muS) {
    mu <- muS / b$k
    ceil_any   <- b$bes * sqrt(pmax(0, (b$U - muS) * (muS - b$L)))                # alpha = 1 outer
    ceil_cont  <- b$bes * sqrt(pmax(0, b$k * (b$hi - mu) * (mu - b$lo)) / b$denom)
    ceil_quasi <- b$bes * sqrt(v_max_alpha(muS, b$k, b$n, b$lo, b$hi) / b$denom)
    floor_free <- sd_min_quasi(muS, b$n)
    floor_amp  <- floor_free / sqrt(b$denom)

    out <- list(ceil_any = ceil_any, ceil_cont = ceil_cont, ceil_quasi = ceil_quasi,
                floor_free = floor_free, floor_amp = floor_amp,
                ceil_obs = NA_real_, floor_att = NA_real_,
                ceil_comb = ceil_quasi, floor_comb = floor_amp, in_band = TRUE)

    if (b$use_obs) {
      R <- b$B - b$A
      inb <- muS >= b$A + R / b$n - 1e-9 & muS <= b$B - R / b$n + 1e-9
      co <- sd_max_quasi_obs(muS, b$n, b$A, b$B)
      fa <- sd_min_attained(muS, b$n, b$A, b$B, quasi = TRUE)
      co[!inb] <- NA_real_; fa[!inb] <- NA_real_
      out$ceil_obs   <- co
      out$floor_att  <- fa
      out$ceil_comb  <- pmin(ceil_quasi, co)
      out$floor_comb <- pmax(floor_amp, fa)
      out$in_band    <- inb
    }
    out
  }

  # ---- PLAUSIBILITY ---------------------------------------------------------

  # population model; depends only on the design and the population inputs,
  # so moving the reported point does not recompute the curves
  plaus_base <- reactive({
    if (!isTRUE(input$use_plaus)) return(NULL)
    lo <- input$lo; hi <- input$hi; k <- input$k; n <- input$n
    vals <- list(lo, hi, k, n, input$pop_mean, input$pop_sd, input$tau)
    validate(need(all(vapply(vals, is.finite, logical(1))),
                  "plausibility: enter the population mean, SD and tau"))
    validate(need(hi > lo && k >= 1 && n >= 2, ""))
    conv <- if (input$scale == "mean") k else 1
    L <- k * lo; U <- k * hi; m <- U - L
    validate(need(abs(m - round(m)) < 1e-9,
                  "plausibility: k * (max - min) must be a whole number (integer scores)"))
    M0 <- input$pop_mean * conv; SD0 <- input$pop_sd * conv
    validate(need(M0 > L + TOL && M0 < U - TOL,
                  "plausibility: the population mean must lie strictly inside the scale"))
    validate(need(SD0 > 0 && input$tau >= 0,
                  "plausibility: population SD must be > 0 and tau >= 0"))
    rho_raw <- rho_mom(M0, SD0, L, U)
    note <- character(0)
    rho <- rho_raw
    if (rho < 0) {
      note <- c(note, sprintf(paste(
        "The population SD is below the binomial SD at this mean (rho = %.3f):",
        "under-dispersed data are outside the beta-binomial family, so the bands use rho = 0."),
        rho_raw))
      rho <- 0
    }
    if (rho >= 1) {
      note <- c(note, "The population SD is at or above the endpoint-only maximum (rho >= 1).")
      rho <- 1
    }
    z <- qnorm(c(0.99, 0.999))
    list(L = L, U = U, m = round(m), n = n, conv = conv, M0 = M0, SD0 = SD0,
         rho_raw = rho_raw, rho = rho, tau = input$tau, note = note,
         mean_range = rbind(M0 - z * SD0 / sqrt(n), M0 + z * SD0 / sqrt(n)))
  })

  # bands across the mean (sum-score units), conditional two-moment method
  plaus_curves <- reactive({
    pb <- plaus_base()
    if (is.null(pb)) return(NULL)
    xs <- seq(pb$L, pb$U, length.out = 123)
    xs <- xs[-c(1, length(xs))]
    r <- strait::sd_plausibility(pb$L, pb$U, M = xs, N = pb$n,
                                 p = c(0.01, 0.001), rho = pb$rho, tau = pb$tau,
                                 method = "chisq2", condition_on_mean = TRUE)
    data.frame(muS = xs,
               lo1  = r$lower[r$p == 0.01],  hi1  = r$upper[r$p == 0.01],
               lo01 = r$lower[r$p == 0.001], hi01 = r$upper[r$p == 0.001])
  })

  # bands at the reported mean: exact when feasible, else two-moment
  plaus_at <- reactive({
    pb <- plaus_base()
    if (is.null(pb)) return(NULL)
    b <- base(); v <- bounds_at(b, b$muS)
    ex <- NULL
    if (isTRUE(input$use_exact) && is.finite(input$mean_dp)) {
      h <- 0.5 * 10^-input$mean_dp * pb$conv
      ex <- tryCatch(exact_bounds(pb$L, pb$U, b$muS, pb$n, pb$rho, pb$tau, h,
                                  c(0.01, 0.001)),
                     error = function(e) NULL)
    }
    if (!is.null(ex)) {
      out <- ex; method <- "exact"
    } else {
      r <- strait::sd_plausibility(pb$L, pb$U, M = b$muS, N = pb$n,
                                   p = c(0.01, 0.001), rho = pb$rho,
                                   tau = pb$tau, method = "chisq2",
                                   condition_on_mean = TRUE)
      out <- data.frame(p = r$p, lower = r$lower, upper = r$upper)
      method <- "two-moment"
    }
    # clip to the possibility envelope, which knows about alpha and the
    # observed extremes; the beta-binomial does not
    fl <- if (is.finite(v$floor_comb)) v$floor_comb else 0
    ce <- if (is.finite(v$ceil_comb)) v$ceil_comb else Inf
    out$lower <- pmin(pmax(out$lower, fl), ce)
    out$upper <- pmax(pmin(out$upper, ce), fl)
    list(b = out, method = method,
         exact_wanted = isTRUE(input$use_exact) && is.null(ex))
  })

  output$plaus_verdict <- renderUI({
    pa <- plaus_at()
    if (is.null(pa)) return(NULL)
    pb <- plaus_base(); b <- base(); s <- b$sS; bb <- pa$b
    in1  <- s >= bb$lower[bb$p == 0.01]  - TOL & s <= bb$upper[bb$p == 0.01]  + TOL
    in01 <- s >= bb$lower[bb$p == 0.001] - TOL & s <= bb$upper[bb$p == 0.001] + TOL
    side <- if (s < bb$lower[bb$p == 0.01]) "small" else "large"
    sd_msg <- if (!identical(verdict()[[1]], TEAL)) {
      list(GREY, "SD not assessed: it is already outside the possible range.")
    } else if (!in01) {
      list(AMBER, sprintf("SD implausibly %s: outside the 0.1%% plausibility band at this mean.", side))
    } else if (!in1) {
      list(AMBER, sprintf("SD unusually %s: outside the 1%% band, inside the 0.1%% band.", side))
    } else {
      list(PLAUS, "SD plausible: inside the 1% band at this mean.")
    }
    mr <- pb$mean_range
    m_msg <- if (b$muS < mr[1, 2] || b$muS > mr[2, 2]) {
      "Mean far from the assumed population mean: outside the 0.1% sampling range for this n."
    } else if (b$muS < mr[1, 1] || b$muS > mr[2, 1]) {
      "Mean unusually far from the assumed population mean: outside the 1% sampling range."
    } else {
      "Mean within the 1% sampling range of the assumed population mean."
    }
    div(class = "verdict", style = paste0("color:", sd_msg[[1]], "; padding-top:0;"),
        HTML(paste0("Plausibility: ", sd_msg[[2]],
                    "<br><span style='font-weight:400; color:", GREY, "'>",
                    m_msg, "</span>")))
  })

  output$plaus_stats <- renderUI({
    pa <- plaus_at()
    if (is.null(pa)) return(NULL)
    pb <- plaus_base(); d <- base()$disp; bb <- pa$b
    f <- function(x) sprintf("%.3f", x * d)
    card <- function(lab, val) div(class = "stat-card",
      style = paste0("border-color:", PLAUS_FILL, ";"),
      div(class = "stat-label", HTML(lab)), div(class = "stat-val", val))
    mr <- pb$mean_range
    cards <- list(
      card("assumed &rho; / &tau;",
           sprintf("%.4f / %.3f", pb$rho_raw, pb$tau)),
      card(paste0("plausible SD, 1% (", pa$method, ")"),
           paste0(f(bb$lower[bb$p == 0.01]), " - ", f(bb$upper[bb$p == 0.01]))),
      card(paste0("plausible SD, 0.1% (", pa$method, ")"),
           paste0(f(bb$lower[bb$p == 0.001]), " - ", f(bb$upper[bb$p == 0.001]))),
      card("sampling range of the mean, 0.1%",
           paste0(f(mr[1, 2]), " - ", f(mr[2, 2]))))
    div(style = "display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin:0 0 12px;", cards)
  })

  output$plaus_note <- renderUI({
    pb <- plaus_base()
    if (is.null(pb)) return(NULL)
    pa <- plaus_at()
    extra <- c(pb$note,
               if (isTRUE(pa$exact_wanted))
                 paste("The exact distribution was too costly (or no integer sum has a mean",
                       "that rounds to the reported mean) at this n and scale length, so the",
                       "verdict uses the two-moment approximation.",
                       if (pb$n <= 50) "That approximation is unreliable in the lower tail at small n." else ""))
    div(class = "lede", style = "font-size:11px; margin-top:8px;",
        HTML(paste0(
          "<b>Plausibility bands</b> (teal): quantiles of the sample SD among samples of this <i>n</i> whose
           mean equals the value on the x axis, if the sum score were beta-binomial with the dispersion
           &rho; implied by the assumed population mean and SD, mixed over &rho; ~ N(&rho;, &tau;&sup2;).
           They hold &rho; fixed as the mean varies, are clipped to the possibility envelope, and ignore
           &alpha; and the observed extremes. Dotted teal lines: the 0.1% sampling range of the mean
           around the assumed population mean (marked &#9651;). Error bar at the reported mean: its 1%
           band from the ", pa$method, " distribution.",
          if (length(extra)) paste0("<br><span style='color:", AMBER, "'>",
                                    paste(extra, collapse = " "), "</span>") else "")))
  })

  verdict <- reactive({
    b <- base(); v <- bounds_at(b, b$muS); s <- b$sS
    if (b$use_obs && !isTRUE(v$in_band)) {
      bd <- band(); d <- b$disp
      return(list(RED, sprintf(
        "Mean outside the feasible band [%.3f, %.3f]: with the observed minimum and maximum both attained, no sample of this n has this mean, at any SD.",
        bd["lo"] * d, bd["hi"] * d)))
    }
    if (s > v$ceil_any + TOL)
      list(RED, "Above the maximum for any bounded data of this size: impossible whatever the internal consistency.")
    else if (b$use_obs && is.finite(v$ceil_obs) && s > v$ceil_obs + TOL)
      list(RED, "Above the maximum possible given the observed range: impossible for data spanning only [observed min, observed max].")
    else if (s > v$ceil_quasi + TOL)
      list(AMBER, "Above the α-ceiling: impossible given the reported α (would require higher internal consistency).")
    else if (b$use_obs && is.finite(v$floor_att) && s < v$floor_att - TOL)
      list(RED, "Below the attained-extremes floor: a sample containing both the reported minimum and maximum must be at least this dispersed.")
    else if (s < v$floor_amp - TOL)
      list(RED, "Below the α-amplified integer floor: inconsistent with the reported α for integer items.")
    else
      list(TEAL, if (b$use_obs)
        "Within the combined bounds: achievable at this α, n, item structure, and observed range."
        else "Within the quasi-integer α-bounds: achievable at this α, n, and item structure.")
  })

  output$verdict <- renderUI({
    m <- verdict()
    div(class = "verdict", style = paste0("color:", m[[1]]), m[[2]])
  })

  output$stats <- renderUI({
    b <- base(); v <- bounds_at(b, b$muS); d <- b$disp
    f <- function(x) ifelse(is.finite(x), sprintf("%.3f", x * d), "--")
    card <- function(lab, val) div(class = "stat-card",
      div(class = "stat-label", HTML(lab)), div(class = "stat-val", val))
    amp <- 1 / sqrt(b$denom)
    cards <- list(
      card("reported SD",              f(b$sS)),
      card("&alpha;-ceiling (quasi)",  f(v$ceil_quasi)),
      card("&alpha;-amplified floor",  f(v$floor_amp)),
      card("&alpha;-free floor",       f(v$floor_free)))
    if (b$use_obs) {
      bd <- band()
      cards <- c(cards, list(
        card("observed-range ceiling", f(v$ceil_obs)),
        card("<b>attained floor</b>",  f(v$floor_att)),
        card("combined floor / ceiling",
             paste0(f(v$floor_comb), " / ", f(v$ceil_comb))),
        card("feasible mean band",
             sprintf("%.2f - %.2f", bd["lo"] * d, bd["hi"] * d))))
    } else {
      cards <- c(cards, list(
        card("max for any data",         f(v$ceil_any)),
        card("&alpha;-ceiling (cont.)",  f(v$ceil_cont)),
        card("amplification 1/&radic;(1&minus;c&alpha;)", sprintf("%.3f", amp)),
        card("units", if (b$meanscore) "mean score" else "sum score")))
    }
    div(style = "display:grid; grid-template-columns:repeat(4,1fr); gap:10px; margin:0 0 12px;", cards)
  })

  output$obs_note <- renderUI({
    b <- base()
    if (!b$use_obs) return(NULL)
    div(class = "lede", style = "font-size:11px; margin-top:8px;",
        HTML("The combined floor is the <i>larger</i> of the &alpha;-amplified and attained floors, and the
              combined ceiling the <i>smaller</i> of the &alpha; and observed-range ceilings. Each is
              separately valid, so their intersection is valid, but it is not guaranteed sharp: a
              configuration attaining one need not attain the other. Outside the shaded feasible band the
              mean itself is unreachable with both extremes attained."))
  })

  output$envelope_plot <- renderPlot({
    b <- base(); d <- b$disp

    # x axis ALWAYS spans the logical range
    mus <- seq(b$L, b$U, length.out = 3001)
    cv  <- bounds_at(b, mus)
    df  <- data.frame(x = mus * d,
                      ceil_any   = cv$ceil_any   * d,
                      ceil_quasi = cv$ceil_quasi * d,
                      floor_free = cv$floor_free * d,
                      floor_amp  = cv$floor_amp  * d,
                      ceil_obs   = cv$ceil_obs   * d,
                      floor_att  = cv$floor_att  * d,
                      ceil_comb  = cv$ceil_comb  * d,
                      floor_comb = cv$floor_comb * d)

    vv <- verdict(); pcol <- vv[[1]]
    ymax <- max(c(df$ceil_any, b$sS * d), na.rm = TRUE) * 1.08
    xr <- b$muS * d
    xmid <- (b$L + b$U) / 2

    p <- ggplot(df, aes(x))

    if (b$use_obs) {
      bd <- band()
      # feasible region under attainment
      p <- p +
        geom_ribbon(aes(ymin = floor_comb, ymax = ceil_comb), fill = "#ede9fe",
                    na.rm = TRUE)
    }

    # PLAUSIBILITY bands, clipped to the possibility envelope, drawn under the
    # envelope curves
    pc <- plaus_curves()
    if (!is.null(pc)) {
      env <- bounds_at(b, pc$muS)
      fl <- ifelse(is.finite(env$floor_comb), env$floor_comb, NA_real_)
      ce <- ifelse(is.finite(env$ceil_comb),  env$ceil_comb,  NA_real_)
      clip <- function(y) pmin(pmax(y, fl), ce)
      pl <- data.frame(x = pc$muS * d,
                       lo1 = clip(pc$lo1) * d,  hi1 = clip(pc$hi1) * d,
                       lo01 = clip(pc$lo01) * d, hi01 = clip(pc$hi01) * d)
      pb <- plaus_base()
      p <- p +
        geom_ribbon(data = pl, aes(x, ymin = lo01, ymax = hi01), inherit.aes = FALSE,
                    fill = PLAUS_FILL, alpha = 0.8, na.rm = TRUE) +
        geom_ribbon(data = pl, aes(x, ymin = lo1, ymax = hi1), inherit.aes = FALSE,
                    fill = PLAUS, alpha = 0.22, na.rm = TRUE) +
        geom_vline(xintercept = pb$mean_range[, 2] * d, colour = PLAUS,
                   linetype = "dotted", linewidth = 0.5) +
        annotate("point", x = pb$M0 * d, y = 0, shape = 2, size = 3, colour = PLAUS)
      pa <- plaus_at()
      b1 <- pa$b[pa$b$p == 0.01, ]
      p <- p + annotate("errorbar", x = xr, ymin = b1$lower * d, ymax = b1$upper * d,
                        width = (b$U - b$L) * d / 60, colour = PLAUS, linewidth = 0.8)
    }

    if (b$use_obs) {
      p <- p +
        geom_vline(xintercept = c(bd["lo"], bd["hi"]) * d,
                   colour = VIO, linetype = "dashed", linewidth = 0.4) +
        geom_vline(xintercept = c(b$A, b$B) * d,
                   colour = "grey65", linetype = "dotted", linewidth = 0.4)
      # faint logical-limits reference
      p <- p +
        geom_line(aes(y = ceil_quasi), colour = "grey72", linewidth = 0.5) +
        geom_line(aes(y = floor_amp),  colour = "grey72", linewidth = 0.5) +
        # observed-extremes curves
        geom_line(aes(y = ceil_comb),  colour = VIO, linewidth = 1.0, na.rm = TRUE) +
        geom_line(aes(y = floor_comb), colour = VIO, linewidth = 1.0, na.rm = TRUE) +
        geom_line(aes(y = floor_att),  colour = INK, linetype = "dashed",
                  linewidth = 0.45, na.rm = TRUE)
      lab_x <- mean(c(bd["lo"], bd["hi"])) * d
      lab_y <- suppressWarnings(max(df$floor_comb, na.rm = TRUE))
      if (is.finite(lab_y))
        p <- p + annotate("text", x = lab_x, y = lab_y * 1.12,
                          label = "attained-extremes floor", size = 3.2, colour = VIO)
    } else {
      p <- p +
        geom_line(aes(y = ceil_any),   colour = "grey70", linewidth = 0.5) +
        geom_line(aes(y = ceil_quasi), colour = INK, linewidth = 0.9) +
        geom_line(aes(y = floor_free), colour = GREY, linetype = "dotted", linewidth = 0.6) +
        geom_line(aes(y = floor_amp),  colour = INK, linewidth = 0.9) +
        annotate("text", x = xmid * d, y = bounds_at(b, xmid)$ceil_any * d,
                 label = "max for any data (α = 1)", vjust = -0.6, size = 3.2, colour = "grey55") +
        annotate("text", x = xmid * d, y = bounds_at(b, xmid)$ceil_quasi * d,
                 label = "α-ceiling (quasi-integer)", vjust = 1.6, size = 3.2, colour = INK)
    }

    p +
      geom_point(data = data.frame(x = xr, y = b$sS * d), aes(x, y),
                 colour = pcol, size = 3.6) +
      scale_x_continuous(if (b$meanscore) "mean score" else "sum-score mean",
                         limits = c(b$L, b$U) * d) +
      scale_y_continuous("standard deviation", limits = c(0, ymax)) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(),
            axis.title = element_text(size = 11, colour = MUTE))
  })
}

shinyApp(ui, server)
