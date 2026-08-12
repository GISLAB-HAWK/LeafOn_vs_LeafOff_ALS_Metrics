#-------------------------------------------------------------------------------
# Name:         pc_summary_viz.R
# Description:  Visualize the point cloud level comparison of the harmonized
#               laz files (see pc_summary.R). Produces a point count
#               table, a return statistics table with conditional return
#               probabilities and mirrored height histograms for leaf-on and
#               leaf-off conditions.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R')


# 01 - configuration
#-------------------

# Point density of the underlying point clouds.
# Used for both input paths and output file names.
PPM <- "ppm20"

# lasR summaries produced by compare_pointcloud.R
inputs <- list(
  leaf_on  = paste0("lon23_",  PPM, "_pc_description.rds"),
  leaf_off = paste0("loff24_", PPM, "_pc_description.rds")
)


# Point classes in the harmonized clouds
CLASS_NONGROUND    <- "1"
CLASS_GROUND <- "2"

# Return statistics: returns above this number are pooled into one row,
# because their individual counts are negligible.
RETURN_GROUP_FROM <- 7

# Height histogram
Z_MAX      <- 40   # highest normalized height shown (m)
ZBIN       <- 2     # bin width of the Z histogram (m), see pc_summary.R
PCT_MAX    <- 10    # y axis limit (%)
LOESS_SPAN <- 0.65

# Shared aesthetics: fill, linetype and labels per dataset
fc          <- c(lon23 = "#afafaf", loff24 = "#e5e6eb")
lt          <- c(lon23 = "solid",   loff24 = "dashed")
dataset_lab <- c(lon23 = "Leaf-on", loff24 = "Leaf-off")


# 02 - data preparation
#----------------------

lon23  <- readRDS(file.path(stats_dir, inputs$leaf_on))
loff24 <- readRDS(file.path(stats_dir, inputs$leaf_off))

#' Convert a named lasR histogram into a long data frame
#'
#' @param x Named numeric vector, names are the bin values.
#' @param dataset Dataset label added as a column.
#'
#' @return Data frame with columns bin, count and dataset.
as_hist_df <- function(x, dataset) {
  data.frame(
    bin     = as.numeric(names(x)),
    count   = as.numeric(x),
    dataset = dataset
  )
}

# Point counts per dataset and class
point_count <- data.frame(
  dataset = c("lon23", "loff24"),
  n_points = c(
    lon23$npoints,
    loff24$npoints
  ),
  non_ground_points = c(
    lon23$npoints_per_class[CLASS_NONGROUND],
    loff24$npoints_per_class[CLASS_NONGROUND]
  ),
  ground_points = c(
    lon23$npoints_per_class[CLASS_GROUND],
    loff24$npoints_per_class[CLASS_GROUND]
  )
) %>%
  mutate(
    canopy_to_ground_ratio = non_ground_points / ground_points,
    vegetation_fraction    = non_ground_points /
      (non_ground_points + ground_points) * 100
  )

print(point_count)


# 03 - return statistics
#-----------------------
# Points per return, long format
ppr <- bind_rows(
  as_hist_df(lon23$npoints_per_return,  "lon23"),
  as_hist_df(loff24$npoints_per_return, "loff24")
) %>%
  rename(return_no = bin, npoints = count) %>%
  filter(return_no >= 1) %>%
  arrange(dataset, return_no)

#' Ordinal label for a return number (1 -> "1st", 4 -> "4th")
ordinal <- function(n) {
  suffix <- ifelse(n == 1, "st", ifelse(n == 2, "nd", ifelse(n == 3, "rd", "th")))
  paste0(n, suffix)
}

# Label of the pooled row, e.g. "7-10th"
group_lab <- paste0(RETURN_GROUP_FROM, "-", ordinal(max(ppr$return_no)))

# Conditional probability that a pulse produces a further return:
# P(n+1 | n) = count(n+1) / count(n). Undefined for the pooled returns.
ppr_stats <- ppr %>%
  group_by(dataset) %>%
  mutate(
    percent = 100 * npoints / sum(npoints),
    p_next  = if_else(return_no < RETURN_GROUP_FROM - 1,
                      lead(npoints) / npoints,
                      NA_real_)
  ) %>%
  ungroup() %>%
  mutate(
    return_grp = if_else(return_no < RETURN_GROUP_FROM,
                         ordinal(return_no), group_lab)
  )

# Pool the high returns, keep p_next of the single returns
ppr_tbl <- ppr_stats %>%
  group_by(dataset, return_grp) %>%
  dplyr::summarise(
    npoints = sum(npoints),
    percent = sum(percent),
    p_next  = first(p_next),
    .groups = "drop"
  ) %>%
  mutate(return_grp = factor(
    return_grp,
    levels = c(ordinal(seq_len(RETURN_GROUP_FROM - 1)), group_lab)
  )) %>%
  arrange(return_grp, dataset)

# Wide layout as in the manuscript: one column block per season
ppr_wide <- ppr_tbl %>%
  pivot_wider(
    names_from  = dataset,
    values_from = c(npoints, percent, p_next)
  ) %>%
  dplyr::select(return_grp,
                npoints_lon23,  percent_lon23,  p_next_lon23,
                npoints_loff24, percent_loff24, p_next_loff24)

print(ppr_wide, n = Inf)


# 04 - height distribution
#-------------------------

# Height distribution, as share of all points within a dataset
z_all <- bind_rows(
  as_hist_df(lon23$z_histogram,  "lon23"),
  as_hist_df(loff24$z_histogram, "loff24")
) %>%
  rename(z = bin) %>%
  filter(z > 0, z <= Z_MAX) %>%
  group_by(dataset) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

# Both panels share this theme; the legend is drawn separately (see below)
custom_theme <- theme_classic(base_size = 14) +
  theme(
    axis.title       = element_text(size = 34),
    axis.text.x      = element_text(size = 30, margin = margin(t = 10)),
    axis.text.y      = element_text(size = 30, margin = margin(r = 10)),
    panel.grid       = element_blank(),
    plot.title       = element_blank(),
    strip.background = element_blank(),
    strip.text       = element_blank(),
    plot.margin      = margin(t = 10, r = 30, b = 10, l = 10),
    legend.position  = "none"
  )

#' Height histogram of one dataset with the loess curves of both datasets
#'
#' @param focal Dataset shown as bars, either "lon23" or "loff24".
#'
#' @return ggplot object with height on the y axis (flipped coordinates).
plot_z_hist <- function(focal) {
  ggplot() +
    geom_col(
      data = filter(z_all, dataset == focal),
      aes(x = z, y = percent, fill = dataset),
      width = ZBIN, alpha = 0.8, colour = "black", linewidth = 0.1
    ) +
    # Both curves in both panels, so the seasons can be compared directly
    geom_smooth(
      data = z_all,
      aes(x = z, y = percent, linetype = dataset, group = dataset),
      method = "loess", formula = y ~ x, se = FALSE,
      colour = "black", linewidth = 1, span = LOESS_SPAN
    ) +
    scale_fill_manual(values = fc) +
    scale_linetype_manual(values = lt) +
    scale_x_continuous(breaks = seq(0, Z_MAX, by = 5)) +
    scale_y_continuous(
      limits = c(0, PCT_MAX),
      breaks = seq(0, PCT_MAX, by = 2.5)
    ) +
    coord_flip() +
    labs(x = "Height (m)", y = "Frequency (%)") +
    custom_theme
}

# Manual legend: fill box and linetype per dataset, label above both
leg_df <- data.frame(
  dataset = c("loff24", "loff24", "lon23", "lon23"),
  type    = c("box", "line", "box", "line"),
  x       = c(1, 1, 2, 2),
  y       = c(2, 1, 2, 1)
)

p_leg <- ggplot() +
  # colored boxes
  geom_rect(
    data = subset(leg_df, type == "box"),
    aes(xmin = x - 0.32, xmax = x + 0.32,
        ymin = y - 0.52, ymax = y + 0.02, fill = dataset),
    colour = "black", linewidth = 0.3
  ) +
  # line types
  geom_segment(
    data = subset(leg_df, type == "line"),
    aes(x = x - 0.32, xend = x + 0.32, y = y, yend = y, linetype = dataset),
    colour = "black", linewidth = 1
  ) +
  # labels above box and line
  geom_text(
    data = data.frame(x = c(1, 2), y = c(2.35, 2.35),
                      dataset = c("loff24", "lon23")),
    aes(x = x, y = y, label = dataset_lab[dataset]),
    size = 8, fontface = "bold"
  ) +
  scale_fill_manual(values = fc) +
  scale_linetype_manual(values = lt) +
  coord_cartesian(xlim = c(0.6, 2.4), ylim = c(0.6, 2.6), expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

# Leaf-off left, leaf-on right, shared legend as inset
p_hist <- (plot_z_hist("loff24") + plot_z_hist("lon23")) +
  inset_element(p_leg, left = 0.48, bottom = 0.76, right = 0.92, top = 1.0)

p_hist


#--- export tables and figures -------------------------------------------------

# directories for the output 
stats_dir <- file.path(output_dir, "stats")
fig_dir <- file.path(output_dir, "figures")

# Point counts per class
write.csv(
  point_count,
  file.path(stats_dir, paste0("pc_point_count_", PPM, ".csv")),
  row.names = FALSE
)

# Return statistics: formatted for the manuscript ...
write.csv(
  ppr_wide,
  file.path(stats_dir, paste0("pc_returns_", PPM, ".csv")),
  row.names = FALSE
)

# ... and unrounded for reuse in R
saveRDS(
  ppr_wide,
  file.path(stats_dir, paste0("pc_returns_", PPM, ".rds"))
)


# save figure 
ggsave(
  file.path(fig_dir, paste0("pc_z_hist_", PPM, ".png")),
  p_hist, dpi = 500, width = 16, height = 10
)

message("tables written to: ", stats_dir)
message("figures written to: ", fig_dir)
