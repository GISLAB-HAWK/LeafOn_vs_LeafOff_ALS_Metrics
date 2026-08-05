#-------------------------------------------------------------------------------
# Name:         calculate_metrics.R
# Description:  Script calculates metrics in forest inventory plots using 
#               the point clouds from leaf-on and leaf-off season.
#               The point clouds used were pulse-based thinned to a 
#               pulse-density of 20 p/m² to ensure consistency.
#               The RSDB (Remote Sensing Database) R-package is used for
#               calculating the metrics. For further information see:
#               https://github.com/environmentalinformatics-marburg/rsdb-data and 
#               https://environmentalinformatics-marburg.github.io/rsdb/docs/r_package_installation/ 
#               The calculated metrics are joined with forest inventory data.
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
# In some cases a restart of R is needed to work with a updated version
# of RSDB package (in RStudio - Session - Terminate R).

# logging into the server
# file containing username and pw in the form "username:password"
fileName <- r'{C:\Users\ffranz\rsdb\rsdb_login.txt}'
userpwd <- readChar(fileName, file.info(fileName)$size)
# read account from file
remotesensing <-RSDB::RemoteSensing$new("https://gislab.hawk.de",userpwd)



# 02: metrics calculation
#-------------------------------------------------------------------------------

# list point cloud layers
remotesensing$pointclouds

# load the point cloud 
pointcloud_lon <- remotesensing$pointcloud('solling23_lon_ppm20_lt20_debugged')
pointcloud_loff <- remotesensing$pointcloud('solling24_loff_ppm20_lt20_debugged')

# list POI layers
remotesensing$poi_groups

# get all POIs of POI layer
pois <- remotesensing$poi_group('inv_attr_plots')

# create spatial object with point geometry
pois_sf <- sf::st_as_sf(pois, coords = c('x', 'y'), crs = 25832)

# buffer the points
pois_sf_buffered <- sf::st_buffer(pois_sf, dist = 13)

# convert to sp format
areas_sp <- as(pois_sf_buffered, 'Spatial')

# extract the polygons objects
polygons <- areas_sp@polygons

# name them using kspnr
names(polygons) <- pois_sf_buffered$name

# set metrics to calculate
metrics <- c(
  'BE_H_MEAN',
  'BE_H_MAX',
  'BE_H_P5',
  'BE_H_P10',
  'BE_H_P20',
  'BE_H_P25',
  'BE_H_P30',
  'BE_H_P40',
  'BE_H_P50',
  'BE_H_P60',
  'BE_H_P70',
  'BE_H_P75',
  'BE_H_P80',
  'BE_H_P90',
  'BE_H_P95',
  'BE_H_P99',
  'BE_H_VAR',
  'BE_H_SD',
  'BE_H_KURTOSIS',
  'BE_H_SKEW',
  'BE_PR_02',
  'BE_PR_10',
  'BE_RD_02',
  'BE_RD_10',
  'point_density',
  'pulse_returns_mean'
)

# calculate indices (on the RSDB server)
# leaf-on and leaf-off
pc_lon_metrics <- pointcloud_lon$indices(areas = polygons, functions = metrics)
pc_loff_metrics <- pointcloud_loff$indices(areas = polygons, functions = metrics)

# rename id column ('name') to kspnr
names(pc_lon_metrics)[names(pc_lon_metrics) == 'name'] <- 'kspnr'
names(pc_loff_metrics)[names(pc_loff_metrics) == 'name'] <- 'kspnr'



# 03: join calculated metrics with forest inventory data
#-------------------------------------------------------------------------------

# output paths for the joined plot metrics (leaf-on and leaf-off)
plot_metrics_lon_file <- file.path(
  processed_data_dir, 'metrics', 'plot_metrics_lon.gpkg'
)
plot_metrics_loff_file <- file.path(
  processed_data_dir, 'metrics', 'plot_metrics_loff.gpkg'
)

if (file.exists(plot_metrics_lon_file) & file.exists(plot_metrics_loff_file)) {

  # load existing joined plot metrics
  cat('Loading existing joined plot metrics from disk...\n')
  plot_metrics_lon <- sf::st_read(plot_metrics_lon_file, quiet = T)
  plot_metrics_loff <- sf::st_read(plot_metrics_loff_file, quiet = T)

} else {

  cat('No joined plot metrics found, joining metrics with inventory data...\n')

  # read forest inventory plots
  # (already clipped to the AOI and filtered)
  inv_attr_plots <- sf::st_read(
    file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots.gpkg')
  )

  # join metrics with inventory data
  plot_metrics_lon <- inv_attr_plots %>%
    dplyr::left_join(
      pc_lon_metrics %>% dplyr::mutate(kspnr = as.integer(kspnr)),
      by = 'kspnr'
    )

  plot_metrics_loff <- inv_attr_plots %>%
    dplyr::left_join(
      pc_loff_metrics %>% dplyr::mutate(kspnr = as.integer(kspnr)),
      by = 'kspnr'
    )

  # write to disk
  sf::st_write(plot_metrics_lon, plot_metrics_lon_file, delete_dsn = T)
  sf::st_write(plot_metrics_loff, plot_metrics_loff_file, delete_dsn = T)

}

head(plot_metrics_lon)
head(plot_metrics_loff)



# 04: correlation analysis of metrics with forest inventory attributes
#-------------------------------------------------------------------------------

# response variables (forest inventory attributes)
response_vars <- c('agb_ha')
# response_vars <- c('agb_ha', 'total_vol_ha', 'merch_vol_ha',
#                    'tree_density', 'basal_area_ha', 'dg')

# metric (predictor) columns
metric_cols <- names(plot_metrics_lon)
metric_cols <- metric_cols[
  grepl('^BE_|^point_density$|^pulse_returns_mean$', metric_cols)
]

# leaf-on and leaf-off datasets (without geometry)
metrics_data <- list(
  lon = sf::st_drop_geometry(plot_metrics_lon),
  loff = sf::st_drop_geometry(plot_metrics_loff)
)

# leaf types to analyse ('all' = deciduous and coniferous combined)
leaf_types <- c('all', 'deciduous', 'coniferous')

# compute Pearson correlations between each metric and each response variable,
# separately for leaf-on/leaf-off and for all/deciduous/coniferous plots
correlation_results <- do.call(rbind, lapply(names(metrics_data), function(leaf_condition) {

  df_full <- metrics_data[[leaf_condition]]

  do.call(rbind, lapply(leaf_types, function(lt) {

    # subset by dominant leaf type ('all' keeps every plot)
    df <- if (lt == 'all') df_full else df_full[df_full$dominant_leaf_type == lt, ]

    do.call(rbind, lapply(response_vars, function(resp) {
      data.frame(
        leaf_condition = leaf_condition,
        leaf_type = lt,
        response_var = resp,
        metric = metric_cols,
        correlation = sapply(metric_cols, function(m) {
          stats::cor(df[[m]], df[[resp]], use = 'complete.obs', method = 'pearson')
        }),
        row.names = NULL
      )
    }))
  }))
}))

# display the correlation table (strongest correlations first)
correlation_results %>%
  dplyr::arrange(response_var, leaf_type, leaf_condition, dplyr::desc(abs(correlation))) %>%
  knitr::kable(digits = 2)

# store the correlation table
write.csv(
  correlation_results,
  file.path(output_dir, 'metric_inventory_correlations.csv'),
  row.names = F
)

# plot: leaf-on vs leaf-off correlation with AGB per metric,
# faceted by leaf type (all / deciduous / coniferous)
corr_agb <- correlation_results %>%
  dplyr::filter(response_var == 'agb_ha') %>%
  dplyr::mutate(
    leaf_condition = ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
    leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
    metric = factor(metric, levels = rev(metric_cols))
  )

corr_plot <- ggplot(corr_agb, aes(x = metric, y = correlation, fill = leaf_condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, colour = 'grey40') +
  facet_wrap(~ leaf_type, nrow = 1) +
  scale_fill_manual(values = c('leaf-on' = '#009E73', 'leaf-off' = '#E69F00'),
                    name = '') +
  labs(x = '', y = 'Pearson correlation with AGB') +
  coord_flip() +
  theme_bw() +
  theme(panel.grid = element_blank())

print(corr_plot)

# save the plot
ggplot2::ggsave(
  filename = file.path(output_dir, 'metric_agb_correlations.pdf'),
  plot = corr_plot,
  width = 10,
  height = 6
)


































