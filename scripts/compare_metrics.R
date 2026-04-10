#-------------------------------------------------------------------------------
# Name:         compare_metrics.R
# Description:  compare Metrics for the two rasters between leaf-on and leaf-off 
#               conditions.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

# source setup script
source('src/setup.R', local = TRUE)

#### Data preparation ####

# read raster files  
loff_r <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling24_loff_ppm20_lt20_indices.tiff')
               ) 
lon_r <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling23_lon_ppm20_lt20_indices.tiff')
  )
 
pix_valid <-  read.csv(
  file.path(output_dir, 'sample_selection.csv')
) %>% filter(status == "valid")

pts <- vect(pix_valid, geom = c("x", "y"), crs = crs(lon_r))

#  Extract corresponding values from both rasters

lon_vals  <- terra::extract(lon_r, pts)
loff_vals <- terra::extract(loff_r, pts)

#lon_vals  <- cbind(pts, lon_vals)
#loff_vals <- cbind(pts, loff_vals)

# extract data
lon_df <- as.data.frame(lon_vals,xy = TRUE, na.rm = T) 
loff_df <- as.data.frame(loff_vals, xy = TRUE, na.rm = T) 

# combine datasets and bring into long format 
df_long <- bind_rows(
  lon_df  %>% mutate(season = "leaf on"),
  loff_df %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(season, ts, ID), names_to = "variable", values_to = "value") %>%
  mutate(species = case_when(
    ts == "1" ~ "decidious",
    ts == "2" ~ "coniferous",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-ts)

# summarize the data
df_summary <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value))

print(df_summary)

# brint data to wide format 
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value) %>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`) %>%
  mutate(diff = leaf_on - leaf_off)
  

# summarize the data
df_summary <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(mean_lon = mean(leaf_on), sd_lon = sd(leaf_on), mean_loff = mean(leaf_off), sd_loff = sd(leaf_off)) %>%
  mutate(diff = mean_lon-mean_loff)

print(df_summary)

#### wilcoxon-test ####

# first, testing for normal distribution 
shapiro <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(
    shapiro = list(shapiro.test(diff)), ### REVISION NEEDED! testing on difference (lon-loff) not the real data 
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = shapiro$statistic,
    p.value   = shapiro$p.value,
    method    = shapiro$method
  ) %>%
  dplyr::select(-shapiro, -method) %>%
  dplyr::mutate(signif = ifelse(p.value < 0.05, "***", ""))

print(shapiro, n = 34)
any(shapiro$p.value>0.05)

# data not normally distributed. using non-parametric wilcoxon-test instead of simple t-test

wilcox_test <- df_wide %>%
  #filter(species == "decidious") %>%
  #filter(species == "coniferous") %>%
  tidyr::drop_na(leaf_on, leaf_off) %>%   # remove NAs
  group_by(variable) %>%
  dplyr::summarise(
    n_pairs   = n(),
    mean_diff = mean(leaf_on - leaf_off),
    
    # run Wilcoxon once
    test      = list(wilcox.test(
      leaf_on, leaf_off,
      paired = TRUE,
      exact  = FALSE
    )),
    
    W         = test[[1]]$statistic,
    p.value   = test[[1]]$p.value,
    
    # rank-biserial correlation
    r_rb      = (2 * W) / (n_pairs * (n_pairs + 1)) - 1,
    
    .groups = "drop"
  ) %>%
  mutate(
    signif = ifelse(p.value < 0.05, "***", ""),
  )
print(wilcox_test)

any(wilcox_test$p.value > 0.05)


# change order of variables for later plotting
lvl = c("BE_H_MAX","BE_H_P90","BE_H_P80","BE_H_P50","BE_H_P20",
        "BE_H_P10", "BE_PR_02","BE_PR_10", "BE_RD_02", "BE_RD_10", 
        "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
        "point_density", "pulse_returns_mean")

df_wide$variable <- factor(df_wide$variable, 
                           levels = lvl )  
df_long$variable <- factor(df_long$variable, 
                           levels = lvl)  


## plot results 
df_norm <- df_wide %>% 
  filter(!is.na(species)) %>% # keeping only samples where species is known 
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



panel_labels <- df_wide %>%
  distinct(variable) %>%
  arrange(variable) %>%              # ensures stable order
  mutate(
    label = paste0(letters[seq_along(variable)], ")")
  )

letter_labeller <- function(x) {
  setNames(paste0(letters[seq_along(x)], ")"), x)
}

# aesthetics for species
sc_cols <- c(coniferous = "grey30", decidious = "grey60")
sc_shapes <- c(coniferous = 17, decidious = 16)
sc_labs <- c(coniferous = "Coniferous", decidious = "Deciduous")   # or whatever you want


s <- ggplot(df_norm, aes(x = leaf_on_n, y = leaf_off_n, colour = species)) +
  geom_point(alpha = 0.4, size = 0.85,aes(shape = species)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.65, linetype = "solid", aes(color = species)) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed",
    linewidth = 0.4
  ) + 
  scale_x_continuous(limits = c(0, 1),breaks = seq(0.2,0.8,0.2)) +
  scale_y_continuous(limits = c(0, 1),breaks = seq(0.2,0.8,0.2)) +
  facet_wrap(~ variable, labeller = as_labeller(letter_labeller), scales = "fixed", ncol = 3, nrow = 8) +
  geom_text(
    data = panel_labels,
    aes(label = label),
    x = -Inf, y = Inf,                # top-left corner
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
  theme_classic(base_size = 12) +
  labs(
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  ) +
  scale_colour_manual(values = sc_cols) +
  scale_shape_manual(values = sc_shapes) +
  theme(
    axis.title = element_text(size = 8, color = "black"),
    axis.text = element_text(size = 6, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    legend.position = "none",
    strip.text = element_blank(), 
    panel.spacing = unit(0, "lines"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )
s


## create custom legend
leg_df <- data.frame(
  species = c("coniferous","coniferous","coniferous",
              "decidious","decidious","decidious"),
  type    = c("label","point","line",
              "label","point","line"),
  x       = c(1,1,1, 2,2,2),
  y       = c(3,2,1, 3,2,1)
)

p_leg <- ggplot() +
  # text labels (top)
  geom_text(
    data = subset(leg_df, type == "label"),
    aes(x = x, y = y, label = sc_labs[species]),
    fontface = "bold",
    size = 2
  ) +
  # points (middle)
  geom_point(
    data = subset(leg_df, type == "point"),
    aes(x = x, y = y, shape = species, colour = species),
    size = 2
  ) +
  # solid lines (bottom)
  geom_segment(
    data = subset(leg_df, type == "line"),
    aes(x = x - 0.2, xend = x + 0.2,
        y = y, yend = y, colour = species),
    linewidth = 0.7
  ) +
  scale_colour_manual(values = sc_cols, guide = "none") +
  scale_shape_manual(values = sc_shapes, guide = "none") +
  coord_cartesian(xlim = c(0.6, 2.4),
                  ylim = c(0.6, 3.4),
                  expand = FALSE) +
  theme_void() +
  theme(legend.position = "none")

s <- s +
  inset_element(
    p_leg,
    left   = 0.75,
    bottom = 0.05,
    right  = 0.95,
    top    = 0.1
  )

print(s)

ggsave(paste0(output_dir,"/pix_scatterplot_ppm20_lt20.png"),s, units = "cm", dpi = 350, width = 14, height = 20)

#### density plot ####

## define custom theme for plotting
custom_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title = element_text(size = 12),
    axis.text.x = element_text(size = 6, margin = margin(t = 2)),
    axis.text.y = element_text(size = 6, margin = margin(t = 10)),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    #plot.margin = margin(t = 10, r = 30, b = 10, l = 10),
    legend.position = "none"
  )

## define linetype, color and legend labels
lt <- c(lon23 = "solid", loff24 = "dashed")
fc <- c(lon23 = "#1b9e77", loff24 = "#d95f02")
labs <- c(lon23 = "Leaf-on", loff24 = "Leaf-off")

fc_season <- c("leaf on" = "#1b9e77",
               "leaf off" = "#d95f02")
labs_season <- c("leaf on" = "Leaf-on",
                 "leaf off" = "Leaf-off")


d <- ggplot(df_long, aes(x = value, fill = season)) +
  geom_density(alpha = 0.7, linewidth = 0.4) +
  facet_wrap(~ variable, scales = "free", ncol = 3) +
  geom_text(
    data = panel_labels,
    aes(label = label),
    x = -Inf, y = Inf,                # top-left corner
    hjust = -0.4, vjust = 1.4,
    size = 4,
    inherit.aes = FALSE
  ) + 
  labs(y = "Density (KDE)", x = "Metric Value") + 
  scale_fill_manual(
    values = fc_season) + 
  custom_theme


# fake density shape
x <- seq(-3, 3, length.out = 300)
dens_leg <- rbind(
  data.frame(season = "leaf on",  x = x, y = dnorm(x)),
  data.frame(season = "leaf off", x = x, y = dnorm(x))
)

p_leg_dens <- ggplot(dens_leg, aes(x = x, y = y, fill = season)) +
  geom_area(alpha = 0.7, colour = "black", linewidth = 0.25) +
  facet_wrap(~season, nrow = 1, labeller = labeller(season = labs_season)) +
  scale_fill_manual(values = fc_season) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )

d <- d +
  inset_element(
    p_leg_dens,
    left = 0.72, bottom = 0.001,
    right = 0.95, top = 0.12
  )

print(d)

ggsave(paste0(output_dir,"/pix_densityplot_ppm20_lt20.png"),d, units = "cm", dpi = 350, width = 14, height = 20)

##############################################################################
