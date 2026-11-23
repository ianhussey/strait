# Build-time constants. Named here so that each is one decision rather than
# dozens of literals. Exported functions spell their defaults out instead, so
# that the values show in their usage.

# Slack on floating-point comparisons of quantities that are whole numbers,
# grid steps or bounds by construction. Any legitimate difference in this
# package is at least 1/n, orders of magnitude larger.
TOLERANCE <- 1e-9

# The same idea one step tighter, for two quantities computed the same way,
# which differ only by accumulated rounding error.
TOLERANCE_TIGHT <- 1e-12

# Guard on the attainable lattice's state table, in cells. A design above it
# takes brimmest()'s targeted route instead.
MAX_LATTICE_CELLS <- 2e7

# Node budget for the constructive witness search. Across 2,719,364 reporting
# grid cells and six designs it was never reached, so it is a safety valve
# rather than a tuning knob.
SEARCH_BUDGET <- 2e5

# Enumeration budget for the exact alpha-conditional Gini envelope, in integer
# sum-score profiles.
MAX_GINI_PROFILES <- 5e5

# Default density of a mean grid, as divisions of the scale width.
MEAN_GRID_DIVISIONS <- 1000

# Points in the dense part of the rounding-envelope candidate grid, before the
# analytic breakpoints are added to it.
ENVELOPE_GRID_POINTS <- 401

# Consistent and inconsistent values in plots, as in scrutiny.
COLOUR_CONSISTENT <- "#0047ff"
COLOUR_INCONSISTENT <- "#ff2a00"

# The rlang data pronoun, bound here so that tidy-eval code can name columns as
# `.data$x` without an @importFrom directive. A data mask shadows this binding
# with its own pronoun; the binding exists so that R CMD check sees a
# definition.
.data <- rlang::.data
