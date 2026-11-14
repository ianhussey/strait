# This script draws man/figures/logo.png; run from the package root. It needs
# strait itself, and the text is set in Space Grotesk, an open font. Before
# drawing, install Space Grotesk, or the text falls back to a generic sans.
#
# Two umbrella plots face each other across the hex like two shores, and the
# word sits in the water between them: a strait. Both plots are the same
# 5-point scale, one turned 180 degrees from the other, and the pair is tilted
# so that their x-axes run parallel to the upper left and lower right edges.

water <- "#1d5a8e"
border <- "#0f3a5f"
land <- "white"

# The umbrella, smooth: both edges are variances of two-point distributions.
# The dome takes the two points to be the scale limits (the Bhatia-Davis
# inequality), each arch the whole numbers on either side of the mean. Here for
# n = 200 responses on a 1-5 scale.
plot <- strait::plot_sd_region(
  l = 1,
  u = 5,
  n = 200,
  rule = "mean_naive_floor",
  by = 0.0005,
  reference = FALSE,
  shade = "inside",
  fill = land,
  line_colour = NA, # no outline: the fill alone is the region
  expand = 0
)
plot <- plot +
  ggplot2::theme_void() +
  ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

# The x-axes sit d units from the center, toward the lower right and upper left
# edges, so the domes point inward and leave a gap of 2 * (d - h) for the word.
# The hex's apothem is 251.1; the border takes 8 of that, and a 10-unit margin
# keeps the umbrellas off it. Each x-axis is as wide as the hex allows at its
# distance d, which is also the widest point of each umbrella, so they cannot
# touch the border anywhere.
d <- 220
h <- 165
inner <- 251.1 - 8 - 10
w <- 2 * tan(pi / 6) * (2 * inner - d)

# The plot as a transparent PNG, w by h logo units and rendered at four times
# that size so it stays sharp in the final 2x logo.
png_path <- tempfile(fileext = ".png")
ggplot2::ggsave(
  png_path,
  plot,
  width = 4 * w,
  height = 4 * h,
  units = "px",
  dpi = 300,
  bg = "transparent",
  device = ragg::agg_png
)
href <- paste0("data:image/png;base64,", base64enc::base64encode(png_path))

# One umbrella, its x-axis centered on (x, y) and turned by angle degrees. In
# its own frame the dome points up, away from the x-axis.
umbrella <- function(x, y, angle) {
  sprintf(
    paste(
      '<g transform="translate(%.1f,%.1f) rotate(%g)">',
      '<image href="%s" x="%g" y="%g" width="%g" height="%g"',
      'preserveAspectRatio="none"/></g>'
    ),
    x,
    y,
    angle,
    href,
    -w / 2,
    -h,
    w,
    h
  )
}

dx <- d * sin(pi / 6)
dy <- d * cos(pi / 6)

font_size <- 84

pts <- function(x, y) paste(sprintf("%.1f,%.1f", x, y), collapse = " ")
hexs <- pts(
  c(260, 511.1, 511.1, 260, 8.9, 8.9),
  c(10, 155, 445, 590, 445, 155)
)

svg <- c(
  paste(
    '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600"',
    'width="520" height="600">'
  ),
  sprintf('<clipPath id="hex"><polygon points="%s"/></clipPath>', hexs),
  '<g clip-path="url(#hex)">',
  sprintf('<rect width="520" height="600" fill="%s"/>', water),
  umbrella(260 + dx, 300 + dy, -30),
  umbrella(260 - dx, 300 - dy, 150),
  # the baseline sits below the center by half the height of the letters, so
  # the word is centered in the gap
  sprintf(
    paste(
      '<text x="260" y="%.1f" transform="rotate(-30 260 300)"',
      'text-anchor="middle" font-family="Space Grotesk, sans-serif"',
      'font-weight="bold" font-size="%g" fill="%s">strait</text>'
    ),
    300 + 0.349 * font_size,
    font_size,
    land
  ),
  "</g>",
  sprintf(
    paste(
      '<polygon points="%s" fill="none" stroke="%s" stroke-width="16"',
      'stroke-linejoin="round"/>'
    ),
    hexs,
    border
  ),
  "</svg>"
)

# Write logo to disk
svg_path <- tempfile(fileext = ".svg")
writeLines(svg, svg_path)
rsvg::rsvg_png(svg_path, "man/figures/logo.png", width = 1040)
