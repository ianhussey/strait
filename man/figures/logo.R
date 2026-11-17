# This script draws man/figures/logo.png; run from the package root. It needs
# strait itself, and the text is set in Fraunces, an open variable font. Before
# drawing, install its variable version (with the weight axis), or systemfonts
# falls back to a generic sans.
#
# Two umbrella plots face each other across the hex like two shores, and the
# word sits in the water between them: a strait. Both plots are the same
# 5-point scale, one turned 180 degrees from the other, and the pair is tilted
# so that their x-axes run parallel to the upper left and lower right edges.
# Like a sea chart, each shore is sandy land ringed by shallow water that
# deepens in bands, and the letters of the word grow as the channel widens.

water <- "#1d5a8e"
border <- "#0f3a5f"
land <- "#f2e2b8"
shallows <- c("#3f7fb3", "#6aa3cc", "#9cc7e3") # deepest first
ink <- "white"

# The umbrella, smooth: both edges are variances of two-point distributions.
# The dome takes the two points to be the scale limits (the Bhatia-Davis
# inequality), each arch the whole numbers on either side of the mean. Here for
# n = 200 responses on a 1-5 scale.
reg <- strait::sd_region_data(
  l = 1,
  u = 5,
  n = 200,
  rule = "mean_naive_floor",
  by = 0.0005
)
reg <- reg[!is.na(reg$max_sd), ]

# The hex's apothem is 251.1; the border takes 8 of that, and a 10-unit margin
# keeps the shores off it. The shallows ring each umbrella b units deep, so the
# umbrella itself must keep b units further in: its x-axis sits d units from
# the center, and is as wide as that shrunken hex allows at that distance.
band <- 7
b <- band * length(shallows)
inner <- 251.1 - 8 - 10 - b
d <- inner
h <- 140
w <- 2 * tan(pi / 6) * (2 * inner - d)

# One umbrella in its own frame: x-axis on y = 0, centered, dome toward -y.
ux <- (reg$mean - 1) / 4 * w - w / 2
uy <- -c(reg$max_sd, rev(reg$min_sd)) / max(reg$max_sd) * h
outline <- paste0(
  "M",
  paste(sprintf("%.2f,%.2f", c(ux, rev(ux)), uy), collapse = "L"),
  "Z"
)

# Each band is the outline stroked ever narrower, lightest last, so that the
# land on top leaves rings of shallows that round off the coast.
shore <- function(x, y, angle) {
  widths <- 2 * band * rev(seq_along(shallows))
  c(
    sprintf('<g transform="translate(%.1f,%.1f) rotate(%g)">', x, y, angle),
    sprintf(
      '<path d="%s" fill="%s" stroke="%s" stroke-width="%g" stroke-linejoin="round"/>',
      outline,
      shallows,
      shallows,
      widths
    ),
    sprintf('<path d="%s" fill="%s"/>', outline, land),
    "</g>"
  )
}

# The word lives in the frame of the strait: centered on the hex, x along the
# channel, y across it (down). There the lower dome is at y = d - height, the
# upper one its mirror image. At each x, g is how far the word may reach from
# the middle line: clear of the shallows by a gap, and inside the hex.
gap <- 10
reach <- b + gap
dome_x <- ux
dome_y <- d + uy[seq_along(ux)]
xs <- seq(-280, 280, by = 0.5)
near <- outer(xs, dome_x, function(x, px) reach^2 - (px - x)^2)
below <- matrix(dome_y, length(xs), length(dome_y), byrow = TRUE) -
  sqrt(pmax(near, 0))
below[near < 0] <- Inf
# in this frame the hex has corners at either end of the channel
ap <- 251.1 - 8 - 10
g <- pmin(
  apply(below, 1, min),
  rev(apply(below, 1, min)),
  ap,
  (ap / cos(pi / 6) - abs(xs)) * tan(pi / 3)
)

# The word, letter by letter, from a variable font, so that each letter can
# have its own size and its own weight. Outlines are at size 100, with the
# baseline at y = 0 and y pointing up.
family <- "Fraunces"
weights <- c(100, 900) # the font's weight axis
weight <- 500 # for the letters in the middle
track <- 4 # extra space after each letter, at size 100
font <- systemfonts::match_fonts(family)
chars <- strsplit("strait", "")[[1]]
letter <- function(ch, wt) {
  v <- systemfonts::font_variation(weight = wt)
  gi <- systemfonts::glyph_info(
    ch,
    path = font$path,
    index = font$index,
    size = 100,
    variation = v
  )
  o <- systemfonts::glyph_outline(
    gi$index,
    font$path,
    font$index,
    size = 100,
    variation = v
  )
  list(outline = o, advance = gi$x_advance)
}

# How thick a letter's lines are on the whole: twice its area over the length
# of its outline, as for a long thin strip.
thickness <- function(o) {
  parts <- split(o, o$contour)
  area <- vapply(
    parts,
    function(p) sum(p$x * c(p$y[-1], p$y[1]) - c(p$x[-1], p$x[1]) * p$y) / 2,
    0
  )
  len <- vapply(
    parts,
    function(p) sum(sqrt(diff(c(p$x, p$x[1]))^2 + diff(c(p$y, p$y[1]))^2)),
    0
  )
  2 * abs(sum(area)) / sum(len)
}
grid <- seq(weights[1], weights[2], by = 10)
thick <- sapply(grid, function(wt) {
  vapply(chars, function(ch) thickness(letter(ch, wt)$outline), 0)
})

# Each letter gets its own size, and the middle of its x-height (the height of
# "r" and "a") sits on the middle line of the channel, so that the bodies of
# all letters line up and only the "t" and the "i" reach out above. The sizes
# follow the channel's width: each is the largest at which a letter the
# word's full height would fit all along that letter's width, and the letter
# itself fits too, so the letters keep their shape and grow as the channel
# widens. The "ra" sits in the middle, the gap between them on the narrowest
# point, and both take the size that fits there. Past the umbrellas' feet the
# channel opens into the corners of the hex, so the letters stop growing at
# cap times that size.
#
# A bigger letter would have thicker lines, so each takes the weight at which
# its lines come out as thick as those of "ra" at `weight`, on average. Sizes,
# weights and places all depend on each other; a few rounds settle them.
cap <- 1.5
pair <- 3:4 # "ra"
n <- length(chars)
s_mid <- 1
size <- rep(1, n)
for (k in 1:30) {
  target <- s_mid * mean(thick[pair, grid == weight])
  wt <- vapply(
    seq_len(n),
    function(i) approx(thick[i, ], grid, target / size[i])$y,
    0
  )
  if (anyNA(wt)) {
    stop("the weight axis is too short for these sizes; change `weight`")
  }
  ls <- Map(letter, chars, wt)
  ink_lo <- vapply(ls, function(l) min(l$outline$x), 0)
  ink_hi <- vapply(ls, function(l) max(l$outline$x), 0)
  bottom <- vapply(ls, function(l) min(l$outline$y), 0)
  top <- vapply(ls, function(l) max(l$outline$y), 0)
  step <- (vapply(ls, function(l) l$advance, 0) + track) * size
  left <- cumsum(step) - step
  left <- left -
    (left[3] + ink_hi[3] * size[3] + left[4] + ink_lo[4] * size[4]) / 2
  half <- max(top - bottom) / 2
  xmid <- (max(top[pair]) + min(bottom[pair])) / 2
  size <- vapply(
    seq_len(n),
    function(i) {
      # and the letter itself, above and below the line, must fit too
      fits <- pmin(
        g / half,
        g / (top[i] - xmid),
        g / (xmid - bottom[i]),
        cap * s_mid
      )
      min(fits[
        xs >= left[i] + ink_lo[i] * size[i] &
          xs <= left[i] + ink_hi[i] * size[i]
      ])
    },
    0
  )
  s_mid <- min(size[pair])
  size[pair] <- s_mid
  if (!all(is.finite(size) & size > 0)) {
    stop("the word is longer than the channel; lower `track` or `cap`")
  }
}
glyphs <- do.call(
  rbind,
  lapply(seq_along(ls), function(i) {
    cbind(ls[[i]]$outline[c("contour", "x", "y")], glyph = i)
  })
)
x_at <- left[glyphs$glyph] + glyphs$x * size[glyphs$glyph]
y_at <- -(glyphs$y - xmid) * size[glyphs$glyph]

# No letter may come nearer to either umbrella than the shallows plus the gap.
shores <- cbind(c(dome_x, -dome_x), c(dome_y, -dome_y))
nearest <- min(sqrt(
  outer(x_at, shores[, 1], "-")^2 +
    outer(y_at, shores[, 2], "-")^2
))
stopifnot(nearest >= reach - 0.5)
word <- vapply(
  split(seq_along(x_at), list(glyphs$glyph, glyphs$contour), drop = TRUE),
  function(i) {
    paste0(
      "M",
      paste(sprintf("%.2f,%.2f", x_at[i], y_at[i]), collapse = "L"),
      "Z"
    )
  },
  ""
)

dx <- d * sin(pi / 6)
dy <- d * cos(pi / 6)

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
  shore(260 + dx, 300 + dy, -30),
  shore(260 - dx, 300 - dy, 150),
  sprintf(
    '<path transform="translate(260,300) rotate(-30)" fill="%s" d="%s"/>',
    ink,
    paste(word, collapse = "")
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
