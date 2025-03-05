#--------------------------------------------------------------------------
# Name:         forest_metrics.R
# Description:  Script calculates several metrics derived from ALS-based
#               point clouds. The point clouds are already filtered and 
#               normalized. The metrics are calculated in terrestrial
#               forest inventory plots (Betriebsinventur Niedersachsen).
#               LAScatalog is used to process several point clouds.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#--------------------------------------------------------------------------



# 01 - file path definitions
#----------------------------

# define raw data directory
raw_data_dir <- 'data/raw/'

# define processed data directory
processed_data_dir <- 'data/processed/'

# define output directory
output_dir <- 'output/'



# 02 - set file paths to input files
#-------------------------------------

# input path to point clouds
path_pc_leafoff <- file.path(raw_data_dir, 'pc_LeafOff_2024')
path_pc_leafon <- file.path(raw_data_dir, 'pc_LeafOn_2023')

# input path to forest inventory data
path_forest_inventory <- file.path(raw_data_dir, 'forest_inventory')



# 03 - data reading
#-------------------------------------

# read point clouds with LAScatalog
pc_ctg_leafoff <- lidR::readLAScatalog(path_pc_leafoff)
pc_ctg_leafon <- lidR::readLAScatalog(path_pc_leafon)

pc_ctg_leafoff
pc_ctg_leafon

# quick plot
lidR::plot(pc_ctg_leafoff)
lidR::plot(pc_ctg_leafon)

# read forest inventory data (BI plots)
# contains timber volume per sample plot
bi_plots <- sf::st_read(file.path(path_forest_inventory, 'vol_stp.gpkg'))

bi_plots
str(bi_plots)



# 04 - data preparation
#-------------------------------------

# check if CRS in both ALS datasets is the same
lidR::crs(pc_ctg_leafoff)
lidR::crs(pc_ctg_leafon)

# assign CRS from pc_ctg_leafoff to pc_ctg_leafon
# ETRS89 / UTM zone 32N
lidR::crs(pc_ctg_leafon) <- lidR::crs(pc_ctg_leafoff)

# reproject BI plots to the CRS of the point clouds
# DHDN / 3-degree Gauss-Kruger zone 3 --> ETRS89 / UTM zone 32N
bi_plots <- sf::st_transform(bi_plots, sf::st_crs(25832))

# crop BI plots to the area covered by leaf-on and leaf-off
bi_plots_aoi_leafon <- sf::st_crop(bi_plots, pc_ctg_leafon)
bi_plots_aoi_leafoff <- sf::st_crop(bi_plots, pc_ctg_leafoff)

# remove outliers (two exceptionally high volumes)
boxplot(bi_plots_aoi_leafon$vol_ha)
boxplot(bi_plots_aoi_leafoff$vol_ha)

outliers_leafon_ds <- boxplot.stats(bi_plots_aoi_leafon$vol_ha)$out
outliers_leafoff_ds <- boxplot.stats(bi_plots_aoi_leafoff$vol_ha)$out

bi_plots_aoi_leafon <- bi_plots_aoi_leafon[
  !bi_plots_aoi_leafon$vol_ha %in% outliers_leafon_ds, 
  ] 
bi_plots_aoi_leafoff <- bi_plots_aoi_leafoff[
  !bi_plots_aoi_leafoff$vol_ha %in% outliers_leafoff_ds, 
  ] 

boxplot(bi_plots_aoi_leafon$vol_ha)
boxplot(bi_plots_aoi_leafoff$vol_ha)

# visualize locations of the BI plots
lidR::plot(pc_ctg_leafon, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(bi_plots_aoi_leafon, col.regions = 'black', cex = 5)

lidR::plot(pc_ctg_leafoff, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(bi_plots_aoi_leafoff, col.regions = 'black', cex = 5)



# 04 - calculation of metrics
#--------------------------------------------------------

# calculate predefined metrics for each plot (radius = 13 m) 
# within the normalized point clouds (leaf-on and leaf-off)
# the .stdmetrics function implemented in the lidR package is used
# non-canopy elements (e.g. stones, shrubs --> points below 2 m) are ignored
# save data frame with the plots and calculated metrics
# if the data frame with the metrics already exists, read it
if (!file.exists(file.path(processed_data_dir, 'plot_metrics_aoi_leafon.RDS'))) {
  
  lidR::opt_filter(pc_ctg_leafon) <- '-drop_z_below 2'
  
  plot_metrics_aoi_leafon <- lidR::plot_metrics(
    pc_ctg_leafon, lidR::.stdmetrics,
    bi_plots_aoi_leafon, radius = 13)
  
  saveRDS(
    plot_metrics_aoi_leafon, 
    file.path(processed_data_dir, 'plot_metrics_aoi_leafon.RDS')
    )
  
} else {
  
  plot_metrics_aoi_leafon <- readRDS(
    file.path(processed_data_dir, 'plot_metrics_aoi_leafon.RDS')
    )
  
}

if (!file.exists(file.path(processed_data_dir, 'plot_metrics_aoi_leafoff.RDS'))) {
  
  lidR::opt_filter(pc_ctg_leafoff) <- '-drop_z_below 2'
  
  plot_metrics_aoi_leafoff <- lidR::plot_metrics(
    pc_ctg_leafoff, lidR::.stdmetrics,
    bi_plots_aoi_leafoff, radius = 13)
  
  saveRDS(
    plot_metrics_aoi_leafoff, 
    file.path(processed_data_dir, 'plot_metrics_aoi_leafoff.RDS'))
  
} else {
  
  plot_metrics_aoi_leafoff <- readRDS(
    file.path(processed_data_dir, 'plot_metrics_aoi_leafoff.RDS'))
  
}
























