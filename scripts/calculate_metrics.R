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

# read forest inventory plots
# (already clipped to the AOI and filtered)
inv_attr_plots <- sf::st_read(
  file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots.gpkg')
)

# join metrics with inventory data
inv_attr_plots_metrics_lon <- inv_attr_plots %>%
  dplyr::left_join(
    pc_lon_metrics %>% dplyr::mutate(kspnr = as.integer(kspnr)),
    by = 'kspnr'
  )

inv_attr_plots_metrics_loff <- inv_attr_plots %>%
  dplyr::left_join(
    pc_loff_metrics %>% dplyr::mutate(kspnr = as.integer(kspnr)),
    by = 'kspnr'
  )

head(inv_attr_plots_metrics_lon)
head(inv_attr_plots_metrics_loff)

# write to disk
sf::st_write(
  inv_attr_plots_metrics_lon,
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon.gpkg')
)
sf::st_write(
  inv_attr_plots_metrics_loff,
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff.gpkg')
)


































