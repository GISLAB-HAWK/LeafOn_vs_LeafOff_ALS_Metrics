# Figure: one random example raster cell per forest type as 3D point cloud
# (fixed oblique view) with its vertical return profile.
# Rows: coniferous / deciduous; columns: leaf-on / leaf-off
# Requires: lidR, ggplot2, patchwork, dplyr
#
# The 3D view is a fixed parallel projection drawn with ggplot2, so every panel
# is shown from exactly the same angle and the figure can be saved as vector PDF.

library(lidR)
library(ggplot2)
library(patchwork)
library(dplyr)

# ---- 1. Inputs ----------------------------------------------------------------
samples_file <- "R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS_Metrics/data/metadata/sample_selection_n400_balanced.csv"
ctg_on  <- readLAScatalog("R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS_Metrics/data/processed_data/pc_leafon_2023/ppm20")
ctg_off <- readLAScatalog("R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS_Metrics/data/processed_data/pc_leafoff_2024/ppm20")

cell_size <- 10   # m
zmax      <- 30   # common height axis
bin_width <- 1    # profile bin width (m)

# View: azimuth (rotation around the vertical axis) and elevation (tilt), degrees
view_theta <- 35
view_phi   <- 20

colour_by_height <- TRUE   # TRUE: light (ground) -> dark (top), helps depth perception

# ---- 2. Draw one random valid cell per forest type ---------------------------
set.seed(36)
cells <- read.csv(samples_file) |>
  filter(status == "valid") |>
  group_by(species) |>
  slice_sample(n = 1) |>
  ungroup() |>
  mutate(type = ifelse(species == "coniferous", "Coniferous", "Deciduous"))
print(cells)

# x/y are assumed to be cell CENTRES
get_points <- function(ctg, cell, condition) {
  las <- clip_rectangle(ctg, cell$x - cell_size / 2, cell$y - cell_size / 2,
                        cell$x + cell_size / 2, cell$y + cell_size / 2)
  data.frame(x         = las$X - (cell$x - cell_size / 2),   # local coords 0..10 m
             y         = las$Y - (cell$y - cell_size / 2),
             Z         = las$HAG,
             ground    = las$Classification == 2L,
             type      = cell$type,
             condition = condition)
}

pts <- bind_rows(lapply(seq_len(nrow(cells)), function(i) bind_rows(
  get_points(ctg_on,  cells[i, ], "Leaf-on"),
  get_points(ctg_off, cells[i, ], "Leaf-off")
)))

prof <- pts |>
  filter(!ground) |>
  mutate(bin = floor(Z / bin_width) * bin_width + bin_width / 2) |>
  count(type, condition, bin) |>
  group_by(type, condition) |>
  mutate(pct = 100 * n / sum(n)) |>
  ungroup()

# ---- 3. Oblique projection ------------------------------------------------------
# Rotate around the cell centre by theta, tilt by phi.
# u = horizontal screen position, v = vertical screen position, depth = distance
# into the screen (larger = further back; drawn first so near points lie on top)
project <- function(x, y, z, theta = view_theta, phi = view_phi) {
  th <- theta * pi / 180; ph <- phi * pi / 180
  xc <- x - cell_size / 2; yc <- y - cell_size / 2
  xr <- xc * cos(th) - yc * sin(th)
  yr <- xc * sin(th) + yc * cos(th)
  data.frame(u = xr, v = yr * sin(ph) + z * cos(ph), depth = yr)
}

# Bounding box: 4 ground corners, 4 top corners, 12 edges
corners <- expand.grid(x = c(0, cell_size), y = c(0, cell_size))
corner_order <- c(1, 2, 4, 3, 1)   # walk around the square
box_edges <- bind_rows(
  # ground and top squares
  bind_rows(lapply(c(0, zmax), function(z) {
    p <- project(corners$x[corner_order], corners$y[corner_order], z)
    data.frame(u = head(p$u, -1), v = head(p$v, -1), u2 = p$u[-1], v2 = p$v[-1])
  })),
  # vertical edges
  bind_rows(lapply(seq_len(4), function(k) {
    p0 <- project(corners$x[k], corners$y[k], 0)
    p1 <- project(corners$x[k], corners$y[k], zmax)
    data.frame(u = p0$u, v = p0$v, u2 = p1$u, v2 = p1$v)
  }))
)

# Reference = right-most vertical edge. The profile is drawn right of it in the
# same coordinate system, so its heights line up exactly with the box.
cp    <- project(corners$x, corners$y, 0)
k_ref <- which.max(cp$u)
u_ref <- cp$u[k_ref]               # screen position of that edge
v0    <- cp$v[k_ref]               # screen height of z = 0 on that edge
zv    <- function(z) v0 + z * cos(view_phi * pi / 180)   # height -> screen v

prof_gap   <- 2.6                  # space between box and profile (for labels)
prof_width <- 6                    # screen width of the profile at its maximum
u_ax       <- u_ref + prof_gap     # profile baseline / height axis

tick_z <- seq(0, zmax, 10)
ticks  <- data.frame(z = tick_z, v = zv(tick_z))

# Common limits for all panels -> identical scale everywhere
lim_u <- c(min(c(box_edges$u, box_edges$u2)) - 1.5, u_ax + prof_width + 0.8)
lim_v <- c(min(c(box_edges$v, box_edges$v2)) - 3.2, max(c(box_edges$v, box_edges$v2)) + 1.5)

# ---- 4. Styling -----------------------------------------------------------------
col_pts  <- "#2a78d6"
col_prof <- "#9a988f"
col_ink  <- "#0b0b0b"
col_box  <- "#9a988f"

theme_fig <- theme_classic(base_size = 8) +
  theme(panel.grid.major.y = element_line(colour = "#e4e3df", linewidth = 0.3),
        axis.line  = element_line(colour = "#52514e", linewidth = 0.3),
        axis.ticks = element_line(colour = "#52514e", linewidth = 0.3),
        axis.text  = element_text(colour = "#52514e"),
        plot.margin = margin(2, 2, 2, 2))

# ---- 5. One panel = 3D cloud + profile in one coordinate system ---------------
make_panel <- function(ft, cond, tag, show_x) {
  other   <- setdiff(c("Leaf-on", "Leaf-off"), cond)
  xmax_prof <- 1.08 * max(filter(prof, type == ft)$pct)   # same scale within a row
  sc <- prof_width / xmax_prof                            # % -> screen units
  
  p_this <- filter(prof, type == ft, condition == cond) |>
    mutate(ymin = zv(bin - 0.45 * bin_width), ymax = zv(bin + 0.45 * bin_width),
           xmin = u_ax, xmax = u_ax + pct * sc)
  # outline of the other condition as explicit path (all bins, empty bins = 0)
  all_bins <- seq(bin_width / 2, zmax, by = bin_width)
  po <- filter(prof, type == ft, condition == other)
  pct_o <- po$pct[match(all_bins, po$bin)]; pct_o[is.na(pct_o)] <- 0
  p_other <- data.frame(
    u = u_ax + rep(pct_o, each = 2) * sc,
    v = zv(as.vector(rbind(all_bins - bin_width / 2, all_bins + bin_width / 2))))
  
  pct_ticks <- data.frame(p = c(0, 10, 20)) |> filter(p <= xmax_prof) |>
    mutate(u = u_ax + p * sc)
  v_base <- zv(0) - 0.6                                       # % axis below the profile
  
  d <- filter(pts, type == ft, condition == cond)
  d <- cbind(d, project(d$x, d$y, d$Z)) |> arrange(desc(depth))
  
  ggplot() +
    # box
    geom_segment(data = box_edges, aes(u, v, xend = u2, yend = v2),
                 colour = col_box, linewidth = 0.25) +
    # points (far first)
    { if (colour_by_height)
      geom_point(data = d, aes(u, v, colour = Z), size = 0.8, stroke = 0, alpha = 0.8)
      else
        geom_point(data = d, aes(u, v), colour = col_pts, size = 0.8, stroke = 0, alpha = 0.8) } +
    scale_colour_gradient(low = "#b7d3f2", high = "#0d3d7a", limits = c(0, zmax),
                          guide = "none") +
    # height grid + axis shared by box and profile
    geom_segment(data = ticks, aes(x = u_ax, xend = u_ax + prof_width, y = v, yend = v),
                 colour = "#e4e3df", linewidth = 0.25) +
    geom_segment(aes(x = u_ax, xend = u_ax, y = zv(0), yend = zv(zmax)),
                 colour = "#52514e", linewidth = 0.3) +
    geom_text(data = ticks, aes(x = u_ax - 0.3, y = v, label = z),
              hjust = 1, size = 2.3, colour = "#52514e") +
    # profile: this condition (bars) and other condition (outline)
    geom_rect(data = p_this, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = col_prof) +
    geom_path(data = p_other, aes(u, v), colour = col_ink, linewidth = 0.3) +
    # % axis
    geom_segment(aes(x = u_ax, xend = u_ax + prof_width, y = v_base, yend = v_base),
                 colour = "#52514e", linewidth = 0.3) +
    geom_segment(data = pct_ticks, aes(x = u, xend = u, y = v_base, yend = v_base - 0.4),
                 colour = "#52514e", linewidth = 0.3) +
    { if (show_x) list(
      geom_text(data = pct_ticks, aes(u, v_base - 0.6, label = p),
                vjust = 1, size = 2.3, colour = "#52514e"),
      annotate("text", x = u_ax + prof_width / 2, y = v_base - 2.2,
               label = "Returns (%)", vjust = 1, size = 2.6)) } +
    annotate("text", x = u_ax + 0.4, y = zv(zmax) + 0.3, label = "Height (m)",
             hjust = 0, vjust = 0, size = 2.3, colour = "#52514e") +
    annotate("text", x = lim_u[1], y = lim_v[2], label = tag,
             hjust = 0, vjust = 1, fontface = "bold", size = 2.8) +
    coord_fixed(xlim = lim_u, ylim = lim_v, expand = FALSE, clip = "off") +
    theme_void()
}

# ---- 6. Assemble 2 x 2 ----------------------------------------------------------
lab <- function(txt, rot = 0) wrap_elements(full = grid::textGrob(
  txt, rot = rot, gp = grid::gpar(fontsize = 10, fontface = "bold")))

a <- make_panel("Coniferous", "Leaf-on",  "(a)", FALSE)
b <- make_panel("Coniferous", "Leaf-off", "(b)", FALSE)
c <- make_panel("Deciduous",  "Leaf-on",  "(c)", TRUE)
d <- make_panel("Deciduous",  "Leaf-off", "(d)", TRUE)

fig <- wrap_plots(
  plot_spacer(),         lab("Leaf-on"), lab("Leaf-off"),
  lab("Coniferous", 90), a,              b,
  lab("Deciduous", 90),  c,              d,
  ncol = 3, widths = c(0.06, 1, 1), heights = c(0.05, 1, 1)
)

fig

ggsave("R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS_Metrics/output/figures/Figure2_pointcloud_3d.png", fig, width = 150, height = 170,
       units = "mm", dpi = 600)

ggsave("Figure2_pointcloud_3d.pdf", fig, width = 180, height = 190, units = "mm")