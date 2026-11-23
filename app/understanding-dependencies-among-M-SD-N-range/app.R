library(shiny)
library(tidyverse)

# Helper
round_to_decimal <- function(number, decimal_places) {
  format_string <- paste0("%.", decimal_places, "f")
  sprintf(format_string, number)
}

# Presets
presets <- list(
  "None (manual entry)" = NULL,
  "Uniform distribution (1–5)" = list(min = 1, max = 5, counts = c(`1` = 20, `2` = 20, `3` = 20, `4` = 20, `5` = 20)),
  "Approximately normal (1–7)" = list(min = 1, max = 7, counts = c(`1` =  0, `2` =  5, `3` = 25, `4` = 40, `5` = 25, `6` =  5, `7` =  0)),
  "Approximately negative-binomial (1–7)" = list(min = 1, max = 7, counts = c(`1` = 50, `2` = 28, `3` = 14, `4` = 6, `5` =  2, `6` =  0, `7` =  0))
)

# UI
ui <- fluidPage(
  titlePanel("Understanding Dependencies Among Mean, SD, N, and Range"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("preset", "Choose a Preset:", choices = names(presets)),
      h4("Scale"),
      numericInput("min_score", "Minimum Score:", 1),
      numericInput("max_score", "Maximum Score:", 7),
      h4("Data"),
      uiOutput("count_inputs")
    ),
    
    mainPanel(
      plotOutput("histPlot")
    )
  )
)

# Server
server <- function(input, output, session) {
  preset_counts <- reactiveVal(NULL)
  
  # When preset changes, update min/max and store counts for use in UI generation
  observeEvent(input$preset, {
    selected <- presets[[input$preset]]
    if (!is.null(selected)) {
      updateNumericInput(session, "min_score", value = selected$min)
      updateNumericInput(session, "max_score", value = selected$max)
      preset_counts(selected$counts)
    } else {
      preset_counts(NULL)
    }
  })
  
  # Generate count inputs based on current min/max, using preset values if available
  output$count_inputs <- renderUI({
    req(input$min_score, input$max_score)
    min_val <- input$min_score
    max_val <- input$max_score
    if (min_val > max_val) return(h5("Minimum score must be less than or equal to maximum score."))
    
    preset <- preset_counts()
    
    lapply(min_val:max_val, function(i) {
      numericInput(
        inputId = paste0("n_", i),
        label = paste0("Count of values of ", i, ":"),
        value = if (!is.null(preset) && as.character(i) %in% names(preset)) preset[[as.character(i)]] else 0,
        min = 0
      )
    })
  })
  
  # Generate plot
  output$histPlot <- renderPlot({
    req(input$min_score, input$max_score)
    min_val <- input$min_score
    max_val <- input$max_score
    if (min_val > max_val) return(NULL)
    
    scores <- min_val:max_val
    counts <- sapply(scores, function(i) input[[paste0("n_", i)]])
    
    # Handle NULLs or missing values safely
    counts[is.na(counts)] <- 0
    
    if (any(counts < 0) || length(counts) != length(scores)) return(NULL)
    
    dat <- tibble(score = rep(scores, times = counts))
    if (nrow(dat) == 0) return(NULL)
    
    M <- round_to_decimal(mean(dat$score), 2)
    SD <- round_to_decimal(sd(dat$score), 2)
    N <- nrow(dat)
    
    ggplot(dat, aes(score)) +
      geom_histogram(binwidth = 1, boundary = min_val - 0.5, fill = "steelblue", color = "black") +
      scale_x_continuous(
        labels = scores,
        breaks = scores,
        limits = c(min_val - 0.5, max_val + 0.5)
      ) +
      theme_linedraw(base_size = 24) +
      theme(panel.grid.minor = element_blank()) +
      ylab("Frequency") +
      xlab("Score") +
      ggtitle(paste0("N = ", N, "\nM = ", M, "\nSD = ", SD, "\nScale range = ", max_val-min_val))
  })
}

# Run app
shinyApp(ui = ui, server = server)