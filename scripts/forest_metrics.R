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

# crop BI plots to the area covered by leaf-on
bi_plots_aoi_leafon <- sf::st_crop(bi_plots, pc_ctg_leafon) 

# visualize locations of the BI plots
lidR::plot(pc_ctg_leafon, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(bi_plots_aoi_leafon, col.regions = 'black', cex = 5)



# 04 - calculation of metrics
#--------------------------------------------------------

# source function for metrics calculation
source('src/calc_metrics.R', local = T)


# calculate the predefined metrics for each plot (radius = 13 m) 
# within the normalized point cloud
# 2 m height threshold according to literature
# save data frame with the plots and calculated metrics
# if the data frame with the metrics already exists, read it
if (!file.exists(paste0(processed_data_dir, 'plot_metrics_pc_solling.RDS'))) {
  
  lidR::opt_filter(ndsm_pc_ctg) <- '-drop_z_below 2'
  
  plot_metrics <- lidR::plot_metrics(ndsm_pc_ctg, ~calc_metrics(Z, R, B),
                                     bi_plots_projected, radius = 13)
  
  saveRDS(plot_metrics, file = paste0(processed_data_dir, 'plot_metrics_pc_solling.RDS'))
  
} else {
  
  plot_metrics <- readRDS(paste0(processed_data_dir, 'plot_metrics_pc_solling.RDS'))
  
}
























