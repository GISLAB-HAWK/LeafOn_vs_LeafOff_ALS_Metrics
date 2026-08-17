#-------------------------------------------------------------------------------
# Name:         compare_metrics_viz.R
# Description:  Figures for the leaf-on vs leaf-off metric comparison.
#               Reads the prepared pixel table written by compare_metrics.R,
#               so the rasters do not have to be extracted again.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#
# Requires:     compare_metrics.R must have been run for the same PPM settings.
#-------------------------------------------------------------------------------

source('src/setup.R')


# 01 - configuration
#-------------------

# Must match the settings used in compare_metrics.R
PPM       <- "ppm20"

# Prepared pixel table produced by compare_metrics.R
pix_file <- file.path(output_dir, "stats", paste0("df_wide_", PPM, ".rds"))


# Figure export settings
FIG_DPI   <- 350
FIG_UNITS <- "cm"


# 02 - load prepared data
#------------------------

df_wide <- readRDS(pix_file)

message("loaded ", nrow(df_wide), " rows from: ", pix_file)

# Variable order comes from the factor levels stored in the RDS file,
# so it stays in sync with compare_metrics.R without duplicating the list.
lvl <- levels(df_wide$variable)

stopifnot(length(lvl) == 24)

# Figures are split into two pages of 12 panels each
vars_1 <- lvl[1:12]
vars_2 <- lvl[13:24]


# 03 - plotting setup
#--------------------

# Shared theme for all figures
custom_theme <- theme_classic(base_size = 12) +
  theme(
    axis.title       = element_text(size = 8, color = "black"),
    axis.text        = element_text(size = 6, color = "black"),
    axis.ticks       = element_line(color = "black"),
    axis.line        = element_line(color = "black"),
    panel.grid       = element_blank(),
    panel.spacing    = unit(0, "lines"),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    plot.title       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    legend.position  = "none"
  )

# Shared aesthetics for species
sc_cols   <- c(coniferous = "grey30", decidious = "grey60")
sc_labs   <- c(coniferous = "Coniferous", decidious = "Deciduous")
sc_shapes <- c(coniferous = 16, decidious = 17)

# Panel letters a) b) c) ... in the order defined by lvl
panel_labels <- df_wide %>%
  distinct(variable) %>%
  arrange(variable) %>%
  mutate(label = paste0(letters[seq_along(variable)], ")"))


# 04 - scatterplots: leaf-on vs leaf-off
#---------------------------------------

# Min-max normalisation per variable across both seasons
df_norm <- df_wide %>%
  filter(!is.na(species)) %>%
  group_by(variable) %>%
  mutate(
    leaf_on_n  = (leaf_on  - min(c(leaf_on, leaf_off), na.rm = TRUE)) /
      (max(c(leaf_on, leaf_off), na.rm = TRUE) -
         min(c(leaf_on, leaf_off), na.rm = TRUE)),
    leaf_off_n = (leaf_off - min(c(leaf_on, leaf_off), na.rm = TRUE)) /
      (max(c(leaf_on, leaf_off), na.rm = TRUE) -
         min(c(leaf_on, leaf_off), na.rm = TRUE))
  ) %>%
  ungroup()

#' Scatterplot panel page for a subset of variables
#'
#' @param vars_subset Character vector of variable names (max. 12).
#'
#' @return A ggplot object with 4 x 3 facets.
make_scatter <- function(vars_subset) {
  df_sub  <- df_norm      %>% filter(variable %in% vars_subset)
  lab_sub <- panel_labels %>% filter(variable %in% vars_subset)
  
  ggplot(df_sub, aes(x = leaf_on_n, y = leaf_off_n, colour = species)) +
    geom_point(alpha = 0.4, size = 0.85, aes(shape = species)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.65,
                linetype = "solid", aes(color = species)) +
    geom_abline(intercept = 0, slope = 1, color = "red",
                linetype = "dashed", linewidth = 0.4) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0.2, 0.8, 0.2)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0.2, 0.8, 0.2)) +
    facet_wrap(~ variable, scales = "fixed", ncol = 3, nrow = 4) +
    geom_text(
      data = lab_sub,
      aes(label = label),
      x = -Inf, y = Inf,
      hjust = -0.4, vjust = 1.4,
      size = 3,
      inherit.aes = FALSE
    ) +
    ggpubr::stat_cor(
      aes(label = ..r.label..),
      method = "spearman",
      cor.coef.name = "rho",
      size = 2.5,
      show.legend = FALSE,
      na.rm = TRUE,
      geom = "text",
      label.x.npc = 0.7,
      label.y.npc = 0.2,
      lineheight = 0.7
    ) +
    labs(x = "Leaf-On Value", y = "Leaf-Off Value") +
    scale_colour_manual(values = sc_cols) +
    scale_shape_manual(values = sc_shapes) +
    custom_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 40, l = 5, unit = "pt"))
}

# Manual legend: label, point and line per species
leg_df <- data.frame(
  species = c("coniferous", "coniferous", "coniferous",
              "decidious", "decidious", "decidious"),
  type    = c("label", "point", "line",
              "label", "point", "line"),
  x       = c(1, 1, 1, 2, 2, 2),
  y       = c(2.1, 1.5, 1, 2.1, 1.5, 1)
)

p_leg <- ggplot() +
  geom_text(
    data = subset(leg_df, type == "label"),
    aes(x = x, y = y, label = sc_labs[species]),
    fontface = "bold",
    size = 2
  ) +
  geom_point(
    data = subset(leg_df, type == "point"),
    aes(x = x, y = y, shape = species, colour = species),
    size = 2
  ) +
  geom_segment(
    data = subset(leg_df, type == "line"),
    aes(x = x - 0.2, xend = x + 0.2, y = y, yend = y, colour = species),
    linewidth = 0.7
  ) +
  scale_colour_manual(values = sc_cols, guide = "none") +
  scale_shape_manual(values = sc_shapes, guide = "none") +
  coord_cartesian(xlim = c(0.6, 2.4), ylim = c(0.6, 3.4), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

s1_final <- make_scatter(vars_1) +
  inset_element(p_leg, left = 0.78, bottom = 0.0,
                right = 0.95, top = 0.10, align_to = "full")

s2_final <- make_scatter(vars_2) +
  inset_element(p_leg, left = 0.78, bottom = 0.0,
                right = 0.95, top = 0.10, align_to = "full")

print(s1_final)
print(s2_final)

ggsave(file.path(output_dir, "figures",
                 paste0("pix_scatterplot_", PPM, "_part1.png")),
       s1_final, units = FIG_UNITS, dpi = FIG_DPI, width = 14, height = 18)

ggsave(file.path(output_dir, "figures",
                 paste0("pix_scatterplot_", PPM, "_part2.png")),
       s2_final, units = FIG_UNITS, dpi = FIG_DPI, width = 14, height = 18)


# 05 - difference distributions
#------------------------------

# Legend built from a dummy normal density, one facet per species
x_leg <- seq(-3, 3, length.out = 300)
dens_leg <- rbind(
  data.frame(species = "coniferous", x = x_leg, y = dnorm(x_leg)),
  data.frame(species = "deciduous",  x = x_leg, y = dnorm(x_leg))
)

p_leg_dens <- ggplot(dens_leg, aes(x = x, y = y, fill = species)) +
  geom_area(alpha = 0.7, colour = "black", linewidth = 0.25) +
  facet_wrap(~ species, nrow = 1, labeller = labeller(species = sc_labs)) +
  scale_fill_manual(values = sc_cols) +
  theme_void() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(size = 6, face = "bold"),
    strip.background = element_blank(),
    plot.margin      = margin(1, 1, 1, 1)
  )

#' Density panel page for a subset of variables
#'
#' @param vars_subset Character vector of variable names (max. 12).
#'
#' @return A ggplot object with 4 x 3 facets.
make_density <- function(vars_subset) {
  df_sub  <- df_wide      %>% filter(variable %in% vars_subset)
  lab_sub <- panel_labels %>% filter(variable %in% vars_subset)
  
  ggplot(df_sub, aes(x = diff, fill = species)) +
    geom_density(alpha = 0.7, linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    facet_wrap(~ variable, scales = "free", ncol = 3) +
    geom_text(
      data = lab_sub,
      aes(label = label),
      x = -Inf, y = Inf,
      hjust = -0.4, vjust = 1.4,
      size = 3,
      inherit.aes = FALSE
    ) +
    labs(y = "Density (KDE)", x = "Difference (leaf_on − leaf_off)") +
    scale_fill_manual(values = sc_cols) +
    custom_theme +
    theme(plot.margin = margin(t = 5, r = 5, b = 45, l = 5, unit = "pt"),
          panel.border = element_blank())
}

d1_final <- make_density(vars_1) +
  inset_element(p_leg_dens, left = 0.75, bottom = 0.02,
                right = 0.93, top = 0.09, align_to = "full")

d2_final <- make_density(vars_2) +
  inset_element(p_leg_dens, left = 0.75, bottom = 0.02,
                right = 0.93, top = 0.09, align_to = "full")

print(d1_final)
print(d2_final)

ggsave(file.path(output_dir, "figures",
                 paste0("pix_diff_distribution_", PPM, "_part1.png")),
       d1_final, units = FIG_UNITS, dpi = FIG_DPI, width = 16, height = 18)

ggsave(file.path(output_dir, "figures",
                 paste0("pix_diff_distribution_", PPM, "_part2.png")),
       d2_final, units = FIG_UNITS, dpi = FIG_DPI, width = 16, height = 18)

message("plots written to: ",  file.path(output_dir, "figures"))
