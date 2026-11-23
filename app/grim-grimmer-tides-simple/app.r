# app.R
# Simplified Shiny app: single output pane with table, plot, and HTML report download.

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(scrutiny)
  library(tides)
  library(rmarkdown)
})

# ---------------- UI ----------------

ui <- fluidPage(
  titlePanel("GRIM+GRIMMER+TIDES Consistency Checker"),
  sidebarLayout(
    sidebarPanel(
      numericInput("digits", "Reported decimal places (digits)", value = 2, min = 0, step = 1),
      numericInput("n_items", "Number of items averaged over (use 1 for sum scores, single-item measures, and measures that are implicitly single-item like 'age in years')", value = 1, min = 1, step = 1),
      #tags$hr(),
      numericInput("mean", "Reported mean", value = 4.00, step = 0.01),
      numericInput("sd",   "Reported SD",   value = 2.15, min = 0, step = 0.01),
      numericInput("n",    "Reported N", value = 14, min = 1, step = 1),
      #tags$hr(),
      numericInput("min",  "Scale minimum (optional)", value = 1, step = 1),
      numericInput("max",  "Scale maximum (optional)", value = 7, step = 1),
      checkboxInput("show_umbrella", "Overlay umbrella bounds on TIDES plot", value = FALSE),
      actionButton("run", "Run checks", class = "btn-primary"),
      br(),br(),
      downloadButton("download_report", "Download report", class = "btn-success"),
      width = 3
    ),
    mainPanel(
      tableOutput("res_table"),
      br(),
      plotOutput("res_plot", height = 520)
    )
  )
)

# ---------------- Helper functions ----------------

compute_with_tides <- function(mean, sd, n, min, max, digits, n_items) {
  dat <- tibble(
    mean  = mean,
    sd    = sd,
    n     = n,
    min   = min,
    max   = max,
    digits = digits,
    n_items = n_items,
    calculate_min_sd = TRUE
  )
  
  dat %>%
    mutate(method = "approximate", verbose = FALSE) %>%
    mutate(
      tides = pmap_dfr(
        list(
          mean = mean, sd = sd, n = n, min = min, max = max,
          n_items = n_items, digits = digits, method = method, verbose = FALSE
        ),
        tides
      )
    ) %>%
    unnest(tides) %>%
    mutate(
      mean_str = restore_zeros(as.character(mean), width = digits),
      sd_str   = restore_zeros(as.character(sd),   width = digits)
    ) %>%
    mutate(grim = pmap(list(x = mean_str, n = n, items = n_items), scrutiny::grim)) %>%
    unnest(grim) %>%
    rename(grim_only = grim) %>%
    mutate(
      grimmer = pmap(
        list(x = mean_str, sd = sd_str, n = n, items = n_items),
        show_reason = TRUE,
        scrutiny::grimmer
      )
    ) %>%
    unnest_wider(grimmer, names_sep = "_") %>%
    rename(
      grim_grimmer   = grimmer_1,
      grimmer_reason = grimmer_2
    ) %>%
    mutate(
      grimmer_only = !grepl("GRIMMER inconsistent", grimmer_reason)
    )
}

compute_without_tides <- function(mean, sd, n, digits, n_items) {
  dat <- tibble(
    mean = mean,
    sd   = sd,
    n    = n,
    digits = digits,
    n_items = n_items
  )
  
  dat %>%
    mutate(
      mean_str = restore_zeros(as.character(mean), width = digits),
      sd_str   = restore_zeros(as.character(sd),   width = digits)
    ) %>%
    mutate(grim = pmap(list(x = mean_str, n = n, items = n_items), scrutiny::grim)) %>%
    unnest(grim) %>%
    rename(grim_only = grim) %>%
    mutate(
      grimmer = pmap(
        list(x = mean_str, sd = sd_str, n = n, items = n_items),
        show_reason = TRUE,
        scrutiny::grimmer
      )
    ) %>%
    unnest_wider(grimmer, names_sep = "_") %>%
    rename(
      grim_grimmer   = grimmer_1,
      grimmer_reason = grimmer_2
    ) %>%
    mutate(
      grimmer_only = !grepl("GRIMMER inconsistent", grimmer_reason)
    )
}

make_summary_table <- function(res, include_tides) {
  if (isTRUE(include_tides) && "tides_consistent" %in% names(res)) {
    res %>%
      select(
        grim_only,
        grimmer_only,
        tides_only = tides_consistent,
        grim_grimmer
      ) %>%
      mutate(all_pass = (grim_grimmer + tides_only) >= 2) %>%
      pivot_longer(everything(), names_to = "Test", values_to = "Result")
  } else {
    res %>%
      select(
        grim_only,
        grimmer_only,
        grim_grimmer
      ) %>%
      pivot_longer(everything(), names_to = "Test", values_to = "Result")
  }
}

# ---------------- Server ----------------

server <- function(input, output, session) {
  
  # Compute results when "Run checks" is clicked
  res_reactive <- eventReactive(input$run, {
    validate(
      need(input$n >= 1, "n must be >= 1"),
      need(input$sd >= 0, "SD must be >= 0")
    )
    
    use_tides <- !is.na(input$min) && !is.na(input$max) && input$n_items == 1
    
    res <- if (use_tides) {
      compute_with_tides(
        mean   = input$mean,
        sd     = input$sd,
        n      = input$n,
        min    = input$min,
        max    = input$max,
        digits = input$digits,
        n_items = input$n_items
      )
    } else {
      compute_without_tides(
        mean   = input$mean,
        sd     = input$sd,
        n      = input$n,
        digits = input$digits,
        n_items = input$n_items
      )
    }
    
    list(res = res, use_tides = use_tides)
  })
  
  # Table
  output$res_table <- renderTable({
    rr <- res_reactive()
    req(rr)
    make_summary_table(rr$res, include_tides = rr$use_tides)
  })
  
  # Plot
  output$res_plot <- renderPlot({
    rr <- res_reactive()
    req(rr, rr$use_tides)
    
    res <- rr$res
    gp  <- plot_tides(res)
    
    if (isTRUE(input$show_umbrella)) {
      umb <- tides::umbrella(
        n       = unique(res$n),
        min     = unique(res$min),
        max     = unique(res$max),
        n_items = 1,
        digits  = unique(res$digits)
      )
      gp <- gp +
        geom_tile(
          data = umb,
          aes(x = mean, y = sd),
          width = 0.02,
          height = 0.02,
          alpha = 0.8
        )
    }
    
    gp
  })
  
  # Report download (uses external report.Rmd)
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("consistency-report-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".html")
    },
    content = function(file) {
      td  <- tempdir()
      rmd <- file.path(td, "report.Rmd")
      # assumes report.Rmd is in the same directory as app.R
      file.copy("report.Rmd", rmd, overwrite = TRUE)
      
      rmarkdown::render(
        input       = rmd,
        output_file = file,
        params      = list(
          mean          = input$mean,
          sd            = input$sd,
          n             = input$n,
          digits        = input$digits,
          n_items       = input$n_items,
          min           = if (!is.na(input$min)) input$min else NULL,
          max           = if (!is.na(input$max)) input$max else NULL,
          show_umbrella = input$show_umbrella
        ),
        envir = new.env(parent = globalenv())
      )
    }
  )
}

shinyApp(ui, server)