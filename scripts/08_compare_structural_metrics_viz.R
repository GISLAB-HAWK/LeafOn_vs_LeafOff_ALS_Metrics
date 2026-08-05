#-------------------------------------------------------------------------------
# Name:         compare_structural_metrics_viz.R
# Description:  Figures for the structural metric comparison. Reads the prepared
#               pixel tables written by compare_structural_metrics.R, so the
#               rasters do not have to be extracted again. Builds one panel per
#               metric: kernel density per season with a leaf-on vs leaf-off
#               scatter plot as inset, split over two pages of four panels.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#
# Requires:     compare_structural_metrics.R must have been run for the same
#               PPM setting.
#-------------------------------------------------------------------------------

source('src/setup.R')


# 01 - configuration
#-------------------

# Must match the settings used in compare_structural_metrics.R
PPM <- "ppm20"

# Prepared pixel tables produced by compare_structural_metrics.R
long_file <- file.path(output_dir, "stats", paste0("df_long_structural_", PPM, ".rds"))
wide_file <- file.path(output_dir, "stats", paste0("df_wide_structural_", PPM, ".rds"))

fig_dir <- file.path(output_dir, "figures")

# Figure export settings
FIG_DPI    <- 350
FIG_UNITS  <- "cm"
FIG_WIDTH  <- 14
FIG_HEIGHT <- 10
N_COL      <- 2    # panel columns per page

# Metric whose density peak sits on the left, so its inset and the season
# legend move to the right half of the panel
VAR_INSET_RIGHT <- "Rumple"

# Season aesthetics
fc_season <- c("leaf on" = "#afafaf", "leaf off" = "#e5e6eb")

# Species aesthetics
sc_cols   <- c(coniferous = "grey30", deciduous = "grey60")
sc_shapes <- c(coniferous = 17,       deciduous = 16)
sc_labs   <- c(coniferous = "Coniferous", deciduous = "Deciduous")

# Axis labels per metric
y_labels <- c(
  ENL0D         = "Effective number of layers (ENL0)",
  ENL1D         = "Effective number of layers (ENL1)",
  ENL2D         = "Effective number of layers (ENL2)",
  Box_dimension = "Db",
  VCI           = "VCI",
  CHM_mean      = "Mean Canopy Height",
  Canopy_cover  = "Canopy Cover",
  Rumple        = "Rumple Index"
)


# 02 - load prepared data
#------------------------

df_long <- readRDS(long_file)
df_wide <- readRDS(wide_file)

message("loaded ", nrow(df_wide), " rows from: ", wide_file)

# Variable order comes from the factor levels stored in the RDS file,
# so it stays in sync with compare_structural_metrics.R without duplicating
# the list.
lvl <- levels(df_wide$variable)

stopifnot(length(lvl) == 8)

# Figures are split into two pages of four panels each
vars_1 <- lvl[1:4]
vars_2 <- lvl[5:8]

# Panel letters a) b) c) ... continuous over both pages
panel_labels <- setNames(paste0(letters[seq_along(lvl)], ")"), lvl)

# Normalized pairs for the scatter insets, species must be known.
# Both seasons share one range per metric, so the 1:1 line stays interpretable.
df_norm <- df_wide %>%
  filter(!is.na(species)) %>%
  group_by(variable) %>%
  mutate(
    .lo = min(c(leaf_on, leaf_off), na.rm = TRUE),
    .hi = max(c(leaf_on, leaf_off), na.rm = TRUE),
    leaf_on_n  = (leaf_on  - .lo) / (.hi - .lo),
    leaf_off_n = (leaf_off - .lo) / (.hi - .lo)
  ) %>%
  ungroup() %>%
  dplyr::select(-.lo, -.hi)


# 03 - plotting setup
#--------------------

density_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title       = element_text(size = 6, margin = margin(t = 1)),
    axis.text.x      = element_text(size = 6, margin = margin(t = 1)),
    axis.text.y      = element_text(size = 6, margin = margin(t = 8)),
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.2),
    panel.grid       = element_blank(),
    plot.title       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    plot.margin      = margin(5, 5, 5, 5),
    legend.position  = "none"
  )

scatter_theme <- theme_classic(base_size = 6) +
  theme(
    axis.title      = element_text(size = 4, colour = "black"),
    axis.text       = element_text(size = 4, colour = "black"),
    axis.ticks      = element_line(linewidth = 0.1, colour = "black"),
    axis.line       = element_line(linewidth = 0.1, colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    legend.position = "none",
    strip.text      = element_blank(),
    panel.spacing   = unit(0, "lines"),
    panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.4)
  )

# Season legend: one tile per season
p_leg_dens <- ggplot(
  data.frame(season = names(fc_season), label = c("Leaf-on", "Leaf-off")),
  aes(x = 1, y = label, fill = season)
) +
  geom_tile(width = 0.1, height = 0.3, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = fc_season) +
  theme_void() +
  theme(legend.position = "none",
        plot.margin     = margin(2, 2, 2, 2),
        axis.text.y     = element_text(size = 6))

# Species legend: label, point and line per species
leg_df <- data.frame(
  species = rep(c("coniferous", "deciduous"), each = 3),
  type    = rep(c("label", "point", "line"), times = 2),
  x       = c(1, 2.3, 2.3, 3, 4.3, 4.3),
  y       = 1
)

p_leg_species <- ggplot() +
  geom_text(
    data = subset(leg_df, type == "label"),
    aes(x = x, y = y, label = sc_labs[species]),
    fontface = "bold", size = 1.2, hjust = 0.1
  ) +
  geom_point(
    data = subset(leg_df, type == "point"),
    aes(x = x, y = y, shape = species, colour = species), size = 0.3
  ) +
  geom_segment(
    data = subset(leg_df, type == "line"),
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y, colour = species),
    linewidth = 0.2
  ) +
  scale_colour_manual(values = sc_cols, guide = "none") +
  scale_shape_manual(values = sc_shapes, guide = "none") +
  coord_cartesian(xlim = c(0.8, 10), ylim = c(0.9, 1.1), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")


# 04 - panels
#------------

#' Density plot of one metric with a scatter plot inset
#'
#' @param v Metric name, must be present in df_long and df_norm.
#' @param inset_right Put the scatter inset on the right half of the panel and
#'   add the season legend. Used for metrics whose density peak sits on the
#'   left, see VAR_INSET_RIGHT.
#'
#' @return patchwork object.
make_panel <- function(v, inset_right = FALSE) {
  
  p_density <- ggplot(filter(df_long, variable == v), aes(x = value, fill = season)) +
    geom_density(alpha = 0.7, linewidth = 0.3) +
    geom_text(
      data = data.frame(label = panel_labels[[v]]),
      aes(label = label), x = -Inf, y = Inf,
      hjust = -0.4, vjust = 1.4, size = 2, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = fc_season) +
    labs(x = y_labels[v], y = "Density (KDE)", title = NULL) +
    density_theme
  
  p_scatter <- ggplot(filter(df_norm, variable == v),
                      aes(x = leaf_on_n, y = leaf_off_n, colour = species)) +
    geom_point(aes(shape = species), alpha = 0.4, size = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.45) +
    # 1:1 line as reference for a perfect season agreement
    geom_abline(intercept = 0, slope = 1, colour = "red",
                linetype = "dashed", linewidth = 0.2) +
    ggpubr::stat_cor(
      aes(label = after_stat(r.label)),
      method = "spearman", cor.coef.name = "rho", size = 1.5,
      show.legend = FALSE, na.rm = TRUE, geom = "text",
      label.x.npc = 0.65, label.y.npc = 0.25, lineheight = 0.7
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0.2, 0.8, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0.2, 0.8, 0.2)) +
    scale_colour_manual(values = sc_cols) +
    scale_shape_manual(values = sc_shapes) +
    labs(x = "Leaf-On Value", y = "Leaf-Off Value") +
    scatter_theme
  
  if (inset_right) {
    p_density +
      inset_element(p_scatter,     left = 0.51, bottom = 0.5, right = 0.99, top = 1.00,
                    align_to = "panel") +
      inset_element(p_leg_species, left = 0.61, bottom = 1.0, right = 1.51, top = 1.05,
                    align_to = "panel") +
      inset_element(p_leg_dens,    left = 0.70, bottom = 0.2, right = 0.93, top = 0.45,
                    align_to = "panel")
  } else {
    p_density +
      inset_element(p_scatter,     left = 0.01, bottom = 0.5, right = 0.47, top = 1.00,
                    align_to = "panel") +
      inset_element(p_leg_species, left = 0.11, bottom = 1.0, right = 0.99, top = 1.05,
                    align_to = "panel")
  }
}

#' One figure page for a subset of metrics
#'
#' @param vars_subset Character vector of metric names (max. 4).
#'
#' @return patchwork object with N_COL columns.
make_page <- function(vars_subset) {
  panels <- lapply(vars_subset, function(v) {
    make_panel(v, inset_right = v %in% VAR_INSET_RIGHT)
  })
  wrap_plots(panels, ncol = N_COL)
}


# 05 - assemble
#--------------

p1_final <- make_page(vars_1)
p2_final <- make_page(vars_2)

print(p1_final)
print(p2_final)


#--- export figures ------------------------------------------------------------

ggsave(
  file.path(fig_dir, paste0("density_structural_metrics_", PPM, "_part1.png")),
  p1_final, units = FIG_UNITS, dpi = FIG_DPI, width = FIG_WIDTH, height = FIG_HEIGHT
)

ggsave(
  file.path(fig_dir, paste0("density_structural_metrics_", PPM, "_part2.png")),
  p2_final, units = FIG_UNITS, dpi = FIG_DPI, width = FIG_WIDTH, height = FIG_HEIGHT
)

message("figures written to: ", fig_dir)
