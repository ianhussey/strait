# Quasi-integer composite bounds WITH optional observed (attained) extremes
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
# Run with:  shiny::runApp("app_observed_extremes_bounds.R")   (needs shiny + ggplot2)

library(shiny)
library(ggplot2)

INK  <- "#1c1917"; MUTE <- "#78716c"; GREY <- "#57534e"
BLUE <- "#1d4ed8"; RED  <- "#dc2626"; TEAL <- "#0d9488"; AMBER <- "#d97706"
VIO  <- "#7c3aed"
TOL  <- 1e-9

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

  div(class = "title-x", HTML("Composite SD bounds: logical limits vs. observed (attained) extremes")),
  div(class = "lede", style = "margin:4px 0 12px;",
      HTML("The possibility envelope for a composite of <i>k</i> integer items. <b>Logical limits alone</b>
            (faint curves) place every observation in [min, max] but require nothing at either end, so the
            floor is driven only by integer granularity and hugs zero. <b>Observed extremes</b> (bold curves)
            additionally require that at least one observation equals the reported minimum and one the
            reported maximum &mdash; which forces the sample to span that range, lifts the floor sharply, and
            restricts which means are reachable at all. The x axis always spans the <i>logical</i> range so
            the two regimes can be compared directly.")),

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
      div(class = "lede", style = "font-size:11px; margin-top:8px;",
          HTML("With <i>k</i> items each on [min, max], the sum score ranges
                [<i>k</i>&middot;min, <i>k</i>&middot;max]; the mean score is the sum &divide; <i>k</i>."))
    ),
    mainPanel(
      width = 9,
      uiOutput("verdict"),
      uiOutput("stats"),
      plotOutput("envelope_plot", height = "460px"),
      uiOutput("obs_note")
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
                    na.rm = TRUE) +
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
