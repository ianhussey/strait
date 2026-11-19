#' strait: Truncation-Induced Dependency among Summary Statistics
#'
#' A forensic meta-science / trustworthiness-assessment toolkit for checking
#' whether reported means, standard deviations and sample sizes measured on a
#' bounded (truncated) scale are mutually consistent. When a measure has a known
#' minimum and maximum, the mean it can take constrains the standard deviation
#' that is arithmetically possible: the two summary statistics are dependent.
#' \code{sd_bounds()} computes the smallest and largest sample SD consistent
#' with a chosen set of constraints (scale limits or attained extremes, sample
#' size, mean, granularity, and a reported Cronbach's alpha), in closed form.
#' \code{brimmer()} turns those bounds into a report-level consistency
#' verdict and its percent-of-maximum-possible (POMP) transforms;
#' \code{brimmer_multiple()} applies it across a data frame. Companion
#' functions trace the feasible envelope (\code{sd_bounds_curve()}), build the
#' jointly GRIM-, GRIMMER- and bounds-consistent grid (\code{umbrella_data()}),
#' and visualise both on the native and POMP scales (\code{sd_bounds_plot()},
#' \code{sd_bounds_pomp_plot()}, \code{umbrella_plot()}).
#'
#' @seealso The package README and vignette (\code{vignette("strait")}) for a
#'   worked example and the method background.
#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"
