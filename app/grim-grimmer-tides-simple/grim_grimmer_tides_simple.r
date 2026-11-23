
# dependencies
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
#library(janitor)
library(scrutiny)
library(tides)

# shiny inputs
digits <- 2 # integer, default to 2
n_items <- 1L # integer, default to 1
mean <- 1.25 # numeric, default to empty
sd <- 1.00 # numeric, default to empty
n <- 17 # integer, default to empty
min <- 1L # integer, default to empty
max <- 7L # integer, default to empty
umbrella <- FALSE # logical, default to FALSE

# # shiny inputs
# digits <- 1 # integer, default to 2
# n_items <- 1L # integer, default to 1
# mean <- 58.7 # numeric, default to empty
# sd <- 10.8 # numeric, default to empty
# n <- 396 # integer, default to empty
# min <- 18 # integer, default to empty
# max <- 65 # integer, default to empty
# umbrella <- FALSE # logical, default to FALSE
# 
# digits <- 1 # integer, default to 2
# n_items <- 1L # integer, default to 1
# mean <- 57.8 # numeric, default to empty
# sd <- 10.7 # numeric, default to empty
# n <- 136 # integer, default to empty
# min <- 18 # integer, default to empty
# max <- 65 # integer, default to empty
# umbrella <- FALSE # logical, default to FALSE
# 
# digits <- 1 # integer, default to 2
# n_items <- 1L # integer, default to 1
# mean <- 59.2 # numeric, default to empty
# sd <- 10.9 # numeric, default to empty
# n <- 260 # integer, default to empty
# min <- 18 # integer, default to empty
# max <- 65 # integer, default to empty
# umbrella <- FALSE # logical, default to FALSE

# shiny logic
if(!is.null(min) & !is.null(max) & n_items == 1){
  
  dat <- data.frame(mean = mean,
                    sd = sd,
                    n = n,
                    min = min,
                    max = max,
                    digits  = digits, 
                    n_items = n_items,
                    calculate_min_sd = TRUE)
  
  res <- dat |>
    mutate(method = "approximate",
           verbose = FALSE) |>
    # tides
    mutate(tides = pmap_dfr(.l = list(mean = mean,
                                      sd = sd,
                                      n = n, 
                                      min = min, 
                                      max = max,
                                      n_items = n_items,
                                      digits = digits,
                                      method = method,
                                      verbose = FALSE),
                            .f = tides)) |> 
    unnest(tides) |>
    # grim/mer prep
    mutate(mean_str = as.character(mean),
           sd_str = as.character(sd),
           mean_str = restore_zeros(mean_str, width = 2),
           sd_str = restore_zeros(sd_str, width = 2)) |>
    # grim
    mutate(grim = pmap(.l = list(x = mean_str, 
                                 n = n,
                                 items = n_items), 
                       .f = scrutiny::grim)) |>
    unnest(grim) |>
    rename(grim_only = grim) |>
    # grimmer
    mutate(grimmer = pmap(.l = list(x = mean_str, 
                                    sd = sd_str,
                                    n = n,
                                    items = n_items),
                          show_reason = TRUE,
                          .f = scrutiny::grimmer)) |>
    #unnest(grimmer) |>
    unnest_wider(grimmer, names_sep = "_") |>
    rename(grim_grimmer = grimmer_1,
           grimmer_reason = grimmer_2) |>
    mutate(result = case_when(grimmer_reason == "Passed all" ~ "Consistent",
                              grimmer_reason == "GRIM inconsistent" ~ "GRIM inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 1)" ~ "GRIMMER inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 2)" ~ "GRIMMER inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 3)" ~ "GRIMMER inconsistent")) |>
    mutate(grimmer_only = case_when(grimmer_reason == "GRIM inconsistent" ~ TRUE,
                                    grimmer_reason == "Passed all" ~ TRUE,
                                    grimmer_reason == "GRIMMER inconsistent (test 1)" ~ FALSE,
                                    grimmer_reason == "GRIMMER inconsistent (test 2)" ~ FALSE,
                                    grimmer_reason == "GRIMMER inconsistent (test 3)" ~ FALSE)) 
  
  if(umbrella){
    dat_umbrella <- umbrella(n = dat$n, min = dat$min, max = dat$max, n_items = 1, digits = dat$digits) 
    
    res_plot <- plot_tides(res) +
      geom_tile(data = dat_umbrella, aes(x = mean, y = sd), width = 0.02, height = 0.02, alpha = 0.8)
  } else {
    res_plot <- plot_tides(res)
  }
  
  res_table <- res |>
    select(
      #relative_location,
      #relative_dispersion,
      grim_only,
      grimmer_only,
      #sd_range_calculable,
      #mean_inside_range,
      #sd_inside_range,
      tides_only = tides_consistent,
      grim_grimmer
    ) |>
    # master variable for passing all three
    mutate(grim_grimmer_tides = if_else(grim_grimmer + tides_only < 2, FALSE, TRUE)) |>
    pivot_longer(cols = everything(),
                 names_to = "Test",
                 values_to = "Result")
  
  # return:
  res_table
  res_plot
  
} else {
  
  dat <- data.frame(mean = mean,
                    sd = sd,
                    n = n,
                    min = min,
                    max = max,
                    digits  = digits, 
                    n_items = n_items,
                    calculate_min_sd = TRUE)
  
  res <- dat |>
    # grim/mer prep
    mutate(mean_str = as.character(mean),
           sd_str = as.character(sd),
           mean_str = restore_zeros(mean_str, width = 2),
           sd_str = restore_zeros(sd_str, width = 2)) |>
    # grim
    mutate(grim = pmap(.l = list(x = mean_str, 
                                 n = n,
                                 items = n_items), 
                       .f = scrutiny::grim)) |>
    unnest(grim) |>
    rename(grim_only = grim) |>
    # grimmer
    mutate(grimmer = pmap(.l = list(x = mean_str, 
                                    sd = sd_str,
                                    n = n,
                                    items = n_items),
                          show_reason = TRUE,
                          .f = scrutiny::grimmer)) |>
    #unnest(grimmer) |>
    unnest_wider(grimmer, names_sep = "_") |>
    rename(grim_grimmer = grimmer_1,
           grimmer_reason = grimmer_2) |>
    mutate(result = case_when(grimmer_reason == "Passed all" ~ "Consistent",
                              grimmer_reason == "GRIM inconsistent" ~ "GRIM inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 1)" ~ "GRIMMER inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 2)" ~ "GRIMMER inconsistent",
                              grimmer_reason == "GRIMMER inconsistent (test 3)" ~ "GRIMMER inconsistent")) |>
    mutate(grimmer_only = case_when(grimmer_reason == "GRIM inconsistent" ~ TRUE,
                                    grimmer_reason == "Passed all" ~ TRUE,
                                    grimmer_reason == "GRIMMER inconsistent (test 1)" ~ FALSE,
                                    grimmer_reason == "GRIMMER inconsistent (test 2)" ~ FALSE,
                                    grimmer_reason == "GRIMMER inconsistent (test 3)" ~ FALSE)) 
  
  res_table <- res |>
    select(
      #relative_location,
      #relative_dispersion,
      grim_only,
      grimmer_only,
      #sd_range_calculable,
      #mean_inside_range,
      #sd_inside_range,
      grim_grimmer
    ) |>
    # master variable for passing all three
    pivot_longer(cols = everything(),
                 names_to = "Test",
                 values_to = "Result")
  
  # return:
  res_table

}
