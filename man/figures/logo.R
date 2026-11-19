# The hex logo: the umbrella plot of a 1-7 scale, drawn as the arches it is made
# of (after Heathers 2017), and mirrored. A sample that uses only the values a
# and b has an SD that traces an arch from a to b as its mean moves between
# them, sd_max_muilwijk(mean, n, a, b). The ceiling is the arch between the
# extremes, the floor (sd_min_integer()) the arches between neighbouring values,
# and the arches between all other pairs lie inside. All are computed with
# strait's own functions. Drawn with one SD unit as long as one mean unit
# divided by the Bessel factor, every arch is a semicircle.
#
# Being circles, they take their size from the hex: each ceiling is tangent to
# the two edges that make the hex's point above (or below) it, less a small
# margin. The two umbrellas stand apart across the hex's waist, so that each
# floor arch faces its mirror image across a narrow channel, the strait, and
# each such open pair holds one letter of the name.
#
# Heathers 2017:
# https://jamesheathers.medium.com/sprite-interlude-the-umbrella-graph-connecting-grim-and-sprite-also-brunch-sucks-1266c629c974
#
# Run from the package root: Rscript man/figures/logo.R
pkgload::load_all(quiet = TRUE)
out <- "man/figures/logo.png"

navy <- "#173d63"
border <- "#0d2740"
paper <- "#f6f1e7"
family <- "Space Grotesk"

# rsvg falls back to another font without a word, so check it is installed
stopifnot(family %in% systemfonts::system_fonts()$family)

# Draw the faint arches inside each canopy?
inner_arches <- TRUE

# Scale points
k <- 7

cx <- 260
cy <- 300 # the hex's centre, on the waist
gap <- 24 # half the width of the channel between the umbrellas
base <- cy - gap # the upper umbrella's rib tips

# Distance from the centre of the ceiling, (cx, base), to the upper edges
# (which meet at (260, 10) at 30 degrees to the horizontal), less the half
# border and a margin.
R <- (base - 10) * cos(pi / 6) - 8 - 9
x <- cx - R + (seq_len(k) - 1) * 2 * R / (k - 1)

# From (mean, SD) to the drawing. The sample size only sets the Bessel factor,
# which scales every arch alike, so the SD axis is scaled to make the ceiling's
# peak R high, as the ceiling is 2R wide.
n <- 30
X <- function(mean) x[1] + (mean - 1) * 2 * R / (k - 1)
px_per_sd <- R / sd_max_muilwijk((1 + k) / 2, n, lower = 1, upper = k)
Y <- function(sd) base - sd * px_per_sd
polyline <- function(mean, sd) {
  paste0("M", paste(sprintf("%.2f,%.2f", X(mean), Y(sd)), collapse = "L"))
}
means <- function(from, to) seq(from, to, length.out = 120 * (to - from) + 1)

m <- means(1, k)
canopy <- paste0(
  polyline(
    c(m, rev(m)),
    c(sd_max_muilwijk(m, n, lower = 1, upper = k), rev(sd_min_integer(m, n)))
  ),
  "Z"
)

pairs <- t(combn(k, 2))
inner <- pairs[pairs[, 2] - pairs[, 1] > 1 & pairs[, 2] - pairs[, 1] < k - 1, ]

arches <- apply(inner, 1, function(p) {
  m <- means(p[1], p[2])
  sprintf(
    '<path d="%s"/>',
    polyline(m, sd_max_muilwijk(m, n, lower = p[1], upper = p[2]))
  )
})

# Drawn in line: the inner arches faint, the outline of the canopy (ceiling
# and floor) bold, so the umbrella stands out from the arches inside it.
half <- c(
  if (inner_arches) {
    c(
      sprintf(
        '<g fill="none" stroke="%s" stroke-width="2" opacity="0.55">',
        paper
      ),
      arches,
      "</g>"
    )
  },
  sprintf(
    '<path d="%s" fill="none" stroke="%s" stroke-width="3.5" stroke-linejoin="round"/>',
    canopy,
    paper
  )
)

# One letter between each pair of facing floor arches, on a common baseline
# that puts the middle of the x-height on the waist. Each letter is centred on
# its arches by its ink, not its advance.
size <- 68
font <- systemfonts::match_fonts(family, weight = "medium")
chars <- strsplit("strait", "")[[1]]

g <- systemfonts::glyph_info(
  chars,
  path = font$path,
  index = font$index,
  size = size
)

ink_mid <- g$x_bearing + g$width / 2
x_mid <- g$y_bearing[g$glyph == "s"] / 2
centres <- (x[-1] + x[-k]) / 2

letters <- sprintf(
  '<text x="%.2f" y="%.2f">%s</text>',
  centres - ink_mid,
  cy + x_mid,
  chars
)

hex <- "260,10 511.1,155 511.1,445 260,590 8.9,445 8.9,155"

svg <- c(
  '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 520 600" width="520" height="600">',
  sprintf('<polygon points="%s" fill="%s"/>', hex, navy),
  half,
  sprintf('<g transform="matrix(1 0 0 -1 0 %g)">', 2 * cy),
  half,
  "</g>",
  sprintf(
    '<g font-family="%s" font-weight="500" font-size="%g" fill="%s">',
    family,
    size,
    paper
  ),
  letters,
  "</g>",
  sprintf(
    '<polygon points="%s" fill="none" stroke="%s" stroke-width="16" stroke-linejoin="round"/>',
    hex,
    border
  ),
  "</svg>"
)

svg_path <- tempfile(fileext = ".svg")
writeLines(svg, svg_path)
rsvg::rsvg_png(svg_path, out, width = 520)
