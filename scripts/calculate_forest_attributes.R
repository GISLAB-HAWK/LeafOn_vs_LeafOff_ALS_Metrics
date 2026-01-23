#-------------------------------------------------------------------------------
# Name:         calculate_forest_attributes.R
# Description:  Script calculates metrics in forest inventory plots using 
#               the point clouds from leaf-on and leaf-off season.
#               Plot coordinates are available in two variants: 
#               RTK-corrected and non-RTK (uncorrected).
#               The RSDB (Remote Sensing Database) R-package is used for this.
#               For further information see
#               https://github.com/environmentalinformatics-marburg/rsdb-data and 
#               https://environmentalinformatics-marburg.github.io/rsdb/docs/r_package_installation/ 
#               The calculated metrics are joined with the forest inventory data.
# Author:       Florian Franz, Svenja Dobelmann
# Contact:      florian.franz@nw-fva.de
#               svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)


# 01: setup the connection to RSDB
#-------------------------------------------------------------------------------
if(!require("remotes")) install.packages("remotes")

# install RSDB package and automatically install updated versions
remotes::install_github("environmentalinformatics-marburg/rsdb/r-package")
# In some cases a restart of R is needed to work with a updated version of RSDB package (in RStudio - Session - Terminate R).

# logging into the server
fileName <- r'{C:\rsdb\rsdb_login.txt}' # file containing your username and pw in the form "username:password"
userpwd <- readChar(fileName, file.info(fileName)$size) #  read account from file
remotesensing <-RSDB::RemoteSensing$new("https://gislab.hawk.de",userpwd)



# 02: metrics calculation
#-------------------------------------------------------------------------------

# list point cloud layers
remotesensing$pointclouds

# load the point cloud 
pointcloud_lon <- remotesensing$pointcloud('solling23_lon_ppm20')
pointcloud_loff <- remotesensing$pointcloud('solling24_loff_ppm20')

# list POI layers
remotesensing$poi_groups

# get all POIs of POI layer
# once with RTK-based coordinates and once without RTK-based coordinates 
pois_rtk <- remotesensing$poi_group('inv_attr_bi_plots_solling_rtk')

# create spatial object with point geometry
pois_rtk_sf <- sf::st_as_sf(pois_rtk, coords = c('x', 'y'), crs = 25832)

# buffer the points
pois_rtk_sf_buffered <- sf::st_buffer(pois_rtk_sf, dist = 13)

# convert to sp format
areas_rtk_sp <- as(pois_rtk_sf_buffered, 'Spatial')

# extract the polygons objects
polygons_rtk <- areas_rtk_sp@polygons

# name them using kspnr
names(polygons_rtk) <- pois_rtk_sf_buffered$name

# set metrics to calculate
metrics <- c(
  'chm_height_mean',
  'ENL0',
  'chm_surface_area',
  'dtm_surface_area',
  'BE_FHD',
  'vegetation_coverage_05m_CHM',
  'vegetation_point_density'
)

# calculate indices (on the RSDB server)
# leaf-on and leaf-of
pc_lon_metrics_rtk <- pointcloud_lon$indices(areas = polygons_rtk, functions = metrics)
pc_loff_metrics_rtk <- pointcloud_loff$indices(areas = polygons_rtk, functions = metrics)

# rename id column to kspnr
pc_metrics_list <- list(
  pc_lon_metrics_rtk = pc_lon_metrics_rtk,
  pc_loff_metrics_rtk = pc_loff_metrics_rtk
)

for (name in names(pc_metrics_list)) {
  names(pc_metrics_list[[name]])[names(pc_metrics_list[[name]]) == 'name'] <- 'kspnr'
}

# unpack back to individual variables
list2env(pc_metrics_list, envir = .GlobalEnv)


# write to disk

for (name in names(pc_metrics_list)) {
  write.csv2(
    pc_metrics_list[[name]],
    file.path(processed_data_dir, 'metrics', 'plt_level', paste0(name, '.csv'))
  )
}



# 03: Data Analysis
#-------------------------------------------------------------------------------

#### Data preparation ####
# read the leaf-off data 
loff_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','pc_loff_metrics_rtk.csv'), 
  header = T, sep = ";", dec = ","
)

# read the leaf-on data 
lon_df <-read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','pc_lon_metrics_rtk.csv'), 
  header = T, sep = ";", dec = ","
)

## calculate rumple index (RI)
lon_df <- lon_df %>%
  mutate(RI = chm_surface_area / dtm_surface_area )


loff_df <- loff_df %>%
  mutate(RI = chm_surface_area / dtm_surface_area )

# get tree species information from plot file 
ts <-read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','vol_stp_GR.csv')
)


ts <- ts %>%
  mutate(dominant_species = case_when(
    dominant_species == "LB" ~ "decidious",
    dominant_species == "NB" ~ "coniferous",
    TRUE ~ NA_character_
  ))

# filter complete cases 
keep <- intersect(lon_df$kspnr, loff_df$kspnr)

# combine datasets in bring into long format 
df_long <- bind_rows(
  lon_df  %>% filter(kspnr %in% keep) %>% mutate(season = "leaf on") %>% select(-X),
  loff_df %>% filter(kspnr %in% keep) %>% mutate(season = "leaf off")%>% select(-X)
) %>%
  pivot_longer(cols = -c(kspnr, season), names_to = "variable", values_to = "value") %>%
  left_join(ts %>% select(kspnr, dominant_species),
            by = c("kspnr" = "kspnr"))  %>%
  rename(species = dominant_species) %>%
  filter(variable %in% c("chm_height_mean", "ENL0", "FHD", "RI", "vegetation_coverage_05m_CHM"))


# summarize the data
df_summary <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value))

print(df_summary)

# brint data to wide format 
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value)%>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`)

# summarize the data
df_summary <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(mean_lon = mean(leaf_on), sd_lon = sd(leaf_on), mean_loff = mean(leaf_off), sd_loff = sd(leaf_off)) %>%
  mutate(diff = mean_lon-mean_loff)

print(df_summary)

#### wilcoxon-test ####
# first, testing for normal distribution 
shapiro <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(
    shapiro = list(shapiro.test(value)), 
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = shapiro$statistic,
    p.value   = shapiro$p.value,
    method    = shapiro$method
  ) %>%
  select(-shapiro, -method) %>%
  dplyr::mutate(signif = ifelse(p.value < 0.05, "***", ""))

print(shapiro, n = 34)

which(shapiro$p.value > 0.05) # lon skewness and SD normally distributed, rest non-normally

# data not normally distributed. using non-parametric wilcoxon-test instead of simple t-test
wilcox_test <- df_wide %>%
  tidyr::drop_na(leaf_on, leaf_off) %>%   # remove NAs
  group_by(variable) %>%
  dplyr::summarise(
    n_pairs  = n(),
    mean_diff = mean(leaf_on - leaf_off),  # mean difference (on - off)
    p.value  = wilcox.test(leaf_on, leaf_off, paired = TRUE)$p.value,
    .groups = "drop"
  )  %>%
  dplyr::mutate(signif = ifelse(p.value < 0.05, "***", ""))

print(wilcox_test)

which(wilcox_test$p.value > 0.05) # significant differences for all pairs, exept BE_RD_10

#### scatterplot ####

## normalize data 
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


## define custom theme for plotting
custom_theme <- theme_classic(base_size = 20) +
  theme(
    axis.title.x = element_text(size = 6, margin = margin(t = 1)),
    axis.title.y = element_text(size = 6, margin = margin(t = 1)),
    axis.text.x = element_text(size = 6, margin = margin(t = 1)),
    axis.text.y = element_text(size = 6, margin = margin(t = 8)),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.2),
    panel.grid = element_blank(),
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    legend.position = "none"
  )


fc_season <- c("leaf on" = "#1b9e77",
               "leaf off" = "#d95f02")
labs_season <- c("leaf on" = "Leaf-on",
                 "leaf off" = "Leaf-off")

leg_season <- data.frame(
  season = c("leaf on", "leaf off"),
  label  = c("Leaf-on", "Leaf-off")
)

p_leg_dens <- ggplot(leg_season, aes(x = 1, y = label, fill = season)) +
  geom_tile(width = 0.1, height = 0.3, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = fc_season) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2),
    axis.text.y = element_text(size = 6)
  )


vars <- unique(df_long$variable)
var_ur <- tail(vars, 1)   # last variable

y_labels <- c(
  chm_height_mean = "Canopy Height",
  ENL0 = "Effective number of layers (ENL0)",
  vegetation_coverage_05m_CHM = "Canopy Cover",
  RI = "Rumple Index (RI)"
)

# aesthetics for species
sc_cols <- c(coniferous = "grey30", decidious = "grey60")
sc_shapes <- c(coniferous = 17, decidious = 16)
sc_labs <- c(coniferous = "Coniferous", decidious = "Deciduous")   # or whatever you want


## create custom legend
leg_df <- data.frame(
  species = c("coniferous","coniferous","coniferous",
              "decidious","decidious","decidious"),
  type    = c("label","point","line",
              "label","point","line"),
  x       = c(1,2.3,2.3, 3,4.3,4.3),
  y       = c(1,1,1, 1,1,1)
)

p_scatter_legend <- ggplot() +
  # text labels (top)
  geom_text(
    data = subset(leg_df, type == "label"),
    aes(x = x, y = y, label = sc_labs[species]),
    fontface = "bold",
    size = 1.2,
    hjust = 0.1
    ) +
  # points (middle)
  geom_point(
    data = subset(leg_df, type == "point"),
    aes(x = x, y = y, shape = species, colour = species),
    size = 0.3
  ) +
  # solid lines (bottom)
  geom_segment(
    data = subset(leg_df, type == "line"),
    aes(x = x - 0.2, xend = x + 0.2,
        y = y, yend = y, colour = species),
    linewidth = 0.2
  ) +
  scale_colour_manual(values = sc_cols, guide = "none") +
  scale_shape_manual(values = sc_shapes, guide = "none") +
  coord_cartesian(
    xlim = c(0.8, 10),
    ylim = c(0.9, 1.1),
    expand = FALSE
  ) + 
  theme_void() +
  theme(legend.position = "none")

# function for plotting density plot with scatterplot inset
make_panel <- function(v){
  
  dsub_long <- df_long %>% filter(variable == v)
  dsub_wide <- df_norm %>% filter(variable == v)
  
  
  p_density <- ggplot(dsub_long, aes(x = value, fill = season)) +
    geom_density(alpha = 0.7, linewidth = 0.3) + 
    labs(y = "Density (KDE)", x = y_labels[v], title = NULL) + 
    scale_fill_manual(
      values = fc_season) + 
    custom_theme 
  
  p_scatter <- ggplot(dsub_wide, aes(x = leaf_on_n, y = leaf_off_n, colour = species)) +
    geom_point(alpha = 0.4, size = 0.3,aes(shape = species)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.45, linetype = "solid", aes(color = species)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "red",
      linetype = "dashed",
      linewidth = 0.2
    ) + 
    scale_x_continuous(limits = c(0, 1),breaks = seq(0.2,0.8,0.2)) +
    scale_y_continuous(limits = c(0, 1),breaks = seq(0.2,0.8,0.2)) +
    ggpubr::stat_cor(
      aes(label = ..r.label..),
      method = "spearman",
      cor.coef.name = "rho",
      size = 1.5,
      show.legend = FALSE,
      na.rm = TRUE,
      geom = "text",
      label.x.npc = 0.65,
      label.y.npc = 0.25,
      lineheight = 0.7
    ) +
    theme_classic(base_size = 6) +
    labs(
      x = "Leaf-On Value",
      y = "Leaf-Off Value"
    ) +
    scale_colour_manual(values = sc_cols) +
    scale_shape_manual(values = sc_shapes) +
    theme(
      axis.title = element_text(size = 4, color = "black"),
      axis.text = element_text(size = 4, color = "black"),
      axis.ticks = element_line(linewidth = 0.1, color = "black"),
      axis.line = element_line(linewidth = 0.1, color = "black"),
      plot.background  = element_rect(fill = NA, color = NA),
      legend.position = "none",
      strip.text = element_blank(), 
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.4
      )
    )
  # conditional inset position
  if (v == var_ur) {
    p_density +
      inset_element(
        p_scatter,
        left = 0.51, bottom = 0.42,
        right = 0.99, top = 0.95,
        align_to = "panel"
      ) +
      inset_element(
        p_scatter_legend,
        left = 0.61, bottom = 0.93,
        right = 1.51, top = 1,
        align_to = "panel"
      ) +
      inset_element(
        p_leg_dens,
        left = 0.7, bottom = 0.2,
        right = 0.93, top = 0.45,
        align_to = "panel"
      ) 
  } else {
    p_density +
      inset_element(
        p_scatter,
        left = 0.01, bottom = 0.42,
        right = 0.47, top = 0.95,
        align_to = "panel"
      ) +
      inset_element(
        p_scatter_legend,
        left = 0.11, bottom = 0.93,
        right = 0.99, top = 1,
        align_to = "panel"
      ) 
  }
}

# Build all variable panels and arrange like your facet (ncol=3)
plot <- wrap_plots(lapply(vars, make_panel), ncol = 2)

plot
ggsave(paste0(output_dir,"/BI_densityplot_ppm20_metrics2.png"), plot ,  units = "cm", dpi = 350, width = 14, height = 10)

#################################### END OF SCRIPT #############################