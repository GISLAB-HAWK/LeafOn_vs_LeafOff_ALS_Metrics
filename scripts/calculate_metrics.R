#-------------------------------------------------------------------------------
# Name:         calculate_metrics.R
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
fileName <- r'{C:\Users\ffranz\rsdb\rsdb_login.txt}' # file containing your username and pw in the form "username:password"
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
pois_non_rtk <- remotesensing$poi_group('inv_attr_bi_plots_solling_non_rtk')

# create spatial object with point geometry
pois_rtk_sf <- sf::st_as_sf(pois_rtk, coords = c('x', 'y'), crs = 25832)
pois_non_rtk_sf <- sf::st_as_sf(pois_non_rtk, coords = c('x', 'y'), crs = 25832)

# buffer the points
pois_rtk_sf_buffered <- sf::st_buffer(pois_rtk_sf, dist = 13)
pois_non_rtk_sf_buffered <- sf::st_buffer(pois_non_rtk_sf, dist = 13)

# convert to sp format
areas_rtk_sp <- as(pois_rtk_sf_buffered, 'Spatial')
areas_non_rtk_sp <- as(pois_non_rtk_sf_buffered, 'Spatial')

# extract the polygons objects
polygons_rtk <- areas_rtk_sp@polygons
polygons_non_rtk <- areas_non_rtk_sp@polygons

# name them using kspnr
names(polygons_rtk) <- pois_rtk_sf_buffered$name
names(polygons_non_rtk) <- pois_non_rtk_sf_buffered$name 

# set metrics to calculate
metrics <- c(
  'BE_H_MEAN',
  'BE_H_MAX',
  'BE_H_P10',
  'BE_H_P20',
  'BE_H_P30',
  'BE_H_P40',
  'BE_H_P50',
  'BE_H_P60',
  'BE_H_P70',
  'BE_H_P80',
  'BE_H_P90',
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
# leaf-on and leaf-of, RTK and non-RTK
pc_lon_metrics_rtk <- pointcloud_lon$indices(areas = polygons_rtk, functions = metrics)
pc_loff_metrics_rtk <- pointcloud_loff$indices(areas = polygons_rtk, functions = metrics)
pc_lon_metrics_non_rtk <- pointcloud_lon$indices(areas = polygons_non_rtk, functions = metrics)
pc_loff_metrics_non_rtk <- pointcloud_loff$indices(areas = polygons_non_rtk, functions = metrics)

# rename id column to kspnr
pc_metrics_list <- list(
  pc_lon_metrics_rtk = pc_lon_metrics_rtk,
  pc_loff_metrics_rtk = pc_loff_metrics_rtk,
  pc_lon_metrics_non_rtk = pc_lon_metrics_non_rtk,
  pc_loff_metrics_non_rtk = pc_loff_metrics_non_rtk
)

for (name in names(pc_metrics_list)) {
  names(pc_metrics_list[[name]])[names(pc_metrics_list[[name]]) == 'name'] <- 'kspnr'
}

# unpack back to individual variables
list2env(pc_metrics_list, envir = .GlobalEnv)



# 03: join calculated metrics with forest inventory data
#-------------------------------------------------------------------------------

# read forest inventory plots 
# (already clipped to the AOI and filtered, see script inv_attr_plots.R)
inv_attr_plots_rtk <- sf::st_read(
  file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots_rtk.gpkg')
)
inv_attr_plots_non_rtk <- sf::st_read(
  file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots_non_rtk.gpkg')
)

# join metrics with inventory data
inv_attr_plots_metrics_lon_rtk <- inv_attr_plots_rtk %>%
  dplyr::left_join(
    pc_lon_metrics_rtk %>% dplyr::mutate(kspnr = as.integer(kspnr)), 
    by = 'kspnr'
  )

inv_attr_plots_metrics_loff_rtk <- inv_attr_plots_rtk %>%
  dplyr::left_join(
    pc_loff_metrics_rtk %>% dplyr::mutate(kspnr = as.integer(kspnr)), 
    by = 'kspnr'
  )

inv_attr_plots_metrics_lon_non_rtk <- inv_attr_plots_non_rtk %>%
  dplyr::left_join(
    pc_lon_metrics_non_rtk %>% dplyr::mutate(kspnr = as.integer(kspnr)), 
    by = 'kspnr'
  )

inv_attr_plots_metrics_loff_non_rtk <- inv_attr_plots_non_rtk %>%
  dplyr::left_join(
    pc_loff_metrics_non_rtk %>% dplyr::mutate(kspnr = as.integer(kspnr)), 
    by = 'kspnr'
  )

head(inv_attr_plots_metrics_lon_rtk)
head(inv_attr_plots_metrics_loff_rtk)
head(inv_attr_plots_metrics_lon_non_rtk)
head(inv_attr_plots_metrics_loff_non_rtk)

# write to disk
output_list <- list(
  plot_metrics_lon_rtk = inv_attr_plots_metrics_lon_rtk,
  plot_metrics_loff_rtk = inv_attr_plots_metrics_loff_rtk,
  plot_metrics_lon_non_rtk = inv_attr_plots_metrics_lon_non_rtk,
  plot_metrics_loff_non_rtk = inv_attr_plots_metrics_loff_non_rtk
)

for (name in names(output_list)) {
  sf::st_write(
    output_list[[name]],
    file.path(processed_data_dir, 'metrics', paste0(name, '.gpkg'))
  )
}


































