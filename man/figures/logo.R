args <- commandArgs(TRUE)
out <- if (length(args)) args[1] else "man/figures/logo.png"
env <- function(name, default) Sys.getenv(name, default)

water <- "#1d5a8e"
border <- "#0f3a5f"
land <- env("LAND", "white")
coast <- "#8dbfe2"
grid_col <- "#8dbfe2"
ink <- env("INK", "#a9d0ec")

# Change to adjust the size of the compass roses
rose_fill_value <- 0.65

reg <- strait::sd_region_data(
  l = 1,
  u = 5,
  n = 200,
  rule = "mean_naive_floor",
  by = 0.0005
)
reg <- reg[!is.na(reg$max_sd), ]

# The water inside the border is a hex with apothem A; its corner is R from
# the center. Each umbrella's x-axis lies d from the center, parallel to the
# upper left or lower right edge, and spans that edge's full length, so that
# the grid lines at means 1 and 5 run into the corners at its ends.
A <- 251.1 - 8
R <- A / cos(pi / 6)
d <- 205
h <- 160
w <- 2 * A * tan(pi / 6)
line <- 2.5 # the light-blue coastline, outside the land

ux <- (reg$mean - 1) / 4 * w - w / 2
top_y <- -reg$max_sd / max(reg$max_sd) * h
low_y <- -reg$min_sd / max(reg$max_sd) * h
path <- function(x, y) {
  paste0("M", paste(sprintf("%.2f,%.2f", x, y), collapse = "L"), "Z")
}
umbrella <- path(c(ux, rev(ux)), c(top_y, rev(low_y)))

# A small compass rose in the middle of each umbrella, north toward the
# dome's peak: an eight-pointed two-tone star on one ring, in the coast's
# light blue.
rose_ink <- env("ROSE_INK", coast)
edge_x <- c(ux, ux)
edge_y <- c(top_y, low_y)
cys <- seq(-h + 10, -10, by = 0.5)
room <- vapply(cys, function(cy) min(sqrt(edge_x^2 + (edge_y - cy)^2)), 0)
rose_y <- cys[which.max(room)]
rose_r <- max(room) / 1.42 * as.numeric(env("ROSE_FILL", rose_fill_value))
polar <- function(r, deg) {
  cbind(r * sin(deg * pi / 180), -r * cos(deg * pi / 180))
}
poly <- function(p, fill, stroke = "none", sw = 0) {
  sprintf(
    '<polygon points="%s" fill="%s" stroke="%s" stroke-width="%.2f" stroke-linejoin="round"/>',
    paste(sprintf("%.2f,%.2f", p[, 1], p[, 2]), collapse = " "),
    fill,
    stroke,
    sw
  )
}
circle <- function(r, sw) {
  sprintf(
    '<circle r="%.2f" fill="none" stroke="%s" stroke-width="%.2f"/>',
    r,
    rose_ink,
    sw
  )
}
star_point <- function(deg, tip, base, half) {
  o <- c(0, 0)
  t <- polar(tip, deg)
  l <- polar(base, deg - half)
  r <- polar(base, deg + half)
  c(
    poly(rbind(o, t, l), land, rose_ink, 0.6),
    poly(rbind(o, t, r), rose_ink)
  )
}
rose <- function(r) {
  c(
    sprintf('<g transform="translate(0,%.2f)">', rose_y),
    circle(0.72 * r, 1.1),
    unlist(lapply(
      c(45, 135, 225, 315),
      star_point,
      tip = 0.75 * r,
      base = 0.2 * r,
      half = 45
    )),
    unlist(lapply(
      c(0, 90, 180, 270),
      star_point,
      tip = 1.1 * r,
      base = 0.2 * r,
      half = 45
    )),
    "</g>"
  )
}
shore <- function(x, y, angle) {
  c(
    sprintf('<g transform="translate(%.1f,%.1f) rotate(%g)">', x, y, angle),
    sprintf(
      '<path d="%s" fill="%s" stroke="%s" stroke-width="%g" stroke-linejoin="round"/>',
      umbrella,
      coast,
      coast,
      2 * line
    ),
    sprintf('<path d="%s" fill="%s"/>', umbrella, land),
    if (Sys.getenv("ROSE", "1") == "1") rose(rose_r),
    "</g>"
  )
}

# The strait's frame: centered on the hex, x along the channel, y across it
# (down). The lower dome is at y = d + top_y, the upper one its point mirror.
# gl(x) is how far below the middle line a letter may reach at x: clear of the
# coastline by a gap, and inside the hex.
gap <- 7
reach <- line + gap
xs <- seq(-300, 300, by = 0.5)
near <- outer(xs, ux, function(x, px) reach^2 - (px - x)^2)
below <- matrix(d + top_y, length(xs), length(ux), byrow = TRUE) -
  sqrt(pmax(near, 0))
below[near < 0] <- Inf
ap <- A - 10
gl <- pmin(
  apply(below, 1, min),
  ap,
  (ap / cos(pi / 6) - abs(xs)) * tan(pi / 3)
)

# The chart's grid, which is also the plots' grid: a line across the channel
# at each whole-number mean, and lines along it at SD 0, 1 and 2 of either
# plot. Each runs from edge to edge of the water, its dashes stretched a
# little so that it starts and ends on a full dash.
dashed <- function(x1, y1, x2, y2, on = 7, off = 5) {
  len <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  k <- len / (round((len - on) / (on + off)) * (on + off) + on)
  sprintf(
    '<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke-dasharray="%.3f %.3f"/>',
    x1,
    y1,
    x2,
    y2,
    on * k,
    off * k
  )
}
across <- (0:4) / 4 * w - w / 2
across_end <- pmin(A, (R - abs(across)) * tan(pi / 3))
along <- c(-1, 1) %o% (d - c(0, 0.5, 1) * h)
along_end <- R - abs(along) / tan(pi / 3)
graticule <- c(
  sprintf(
    paste(
      '<g transform="translate(260,300) rotate(-30)" stroke="%s"',
      'stroke-width="1.3" opacity="0.45">'
    ),
    grid_col
  ),
  dashed(across, -across_end, across, across_end),
  dashed(-along_end, along, along_end, along),
  "</g>"
)

# The word, in the italic that sea charts use for the names of waters. The
# letters come in pairs about the middle, "ra", then "ti", then "st"; each pair
# takes the largest size at which it clears both shores and the hex, at most
# `jump` times the size of the pair inside it. The middle of the x-height sits
# on the middle line of the channel, so the word widens evenly toward both
# ends. The whole word is then set at `scale` of that size, to leave water
# around it.
family <- env("FONT", "STIX Two Text")
jump <- 1.35
scale <- 0.72
track <- 14
font <- systemfonts::match_fonts(family, italic = TRUE)
chars <- strsplit("strait", "")[[1]]
ls <- lapply(chars, function(ch) {
  gi <- systemfonts::glyph_info(
    ch,
    path = font$path,
    index = font$index,
    size = 100
  )
  list(
    o = systemfonts::glyph_outline(gi$index, font$path, font$index, size = 100),
    adv = gi$x_advance
  )
})
glyphs <- do.call(
  rbind,
  lapply(seq_along(ls), function(i) {
    cbind(ls[[i]]$o[c("contour", "x", "y")], glyph = i)
  })
)
ink_lo <- vapply(ls, function(l) min(l$o$x), 0)
ink_hi <- vapply(ls, function(l) max(l$o$x), 0)
adv <- vapply(ls, function(l) l$adv, 0)
xmid <- max(ls[[3]]$o$y, ls[[4]]$o$y) / 2 # middle of the x-height

ring <- c(3, 2, 1, 1, 2, 3)
place <- function(size) {
  step <- (adv + track) * size
  left <- cumsum(step) - step
  left <- left -
    (left[3] + ink_hi[3] * size[3] + left[4] + ink_lo[4] * size[4]) / 2
  list(
    x = left[glyphs$glyph] + glyphs$x * size[glyphs$glyph],
    y = -(glyphs$y - xmid) * size[glyphs$glyph]
  )
}
fits <- function(size) {
  p <- place(size)
  all(p$y <= approx(xs, gl, p$x)$y & -p$y <= approx(xs, gl, -p$x)$y)
}
size <- rep(0, 6)
for (k in 1:3) {
  lo <- if (k == 1) 0.1 else max(size)
  hi <- if (k == 1) 5 else jump * max(size)
  for (it in 1:40) {
    m <- (lo + hi) / 2
    size[ring >= k] <- m
    if (fits(size)) lo <- m else hi <- m
  }
  size[ring >= k] <- lo
}
word <- place(scale * size)
letters <- vapply(
  split(seq_along(word$x), list(glyphs$glyph, glyphs$contour), drop = TRUE),
  function(i) path(word$x[i], word$y[i]),
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
  graticule,
  shore(260 + dx, 300 + dy, -30),
  shore(260 - dx, 300 - dy, 150),
  sprintf(
    '<path transform="translate(260,300) rotate(-30)" fill="%s" d="%s"/>',
    ink,
    paste(letters, collapse = "")
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
svg_path <- tempfile(fileext = ".svg")
writeLines(svg, svg_path)
rsvg::rsvg_png(svg_path, out, width = as.numeric(env("PX", 520)))
