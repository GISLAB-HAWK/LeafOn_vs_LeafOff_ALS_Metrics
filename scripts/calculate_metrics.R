#-------------------------------------------------------------------------------
# Name:         calculate_metrics.R
# Description:  Script calculates metrics for raster cells with the Area Based approach
#               using the point clouds from leaf-on and leaf-off season and for 
#               different point densities.
#               The RSDB (Remote Sensing Database) R-package is used for this.
#               For further information see
#               https://github.com/environmentalinformatics-marburg/rsdb-data and 
#               https://environmentalinformatics-marburg.github.io/rsdb/docs/r_package_installation/ 
#               The pointclouds and a output raster dummy need to be uploaded to the RSDB before
#               running this script!
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#               
#-------------------------------------------------------------------------------

# source setup script
source('src/setup.R', local = TRUE)

# 01: setup the connection to RSDB
#-------------------------------------------------------------------------------
#if(!require("remotes")) install.packages("remotes")

# install RSDB package and automatically install updated versions
#remotes::install_github("environmentalinformatics-marburg/rsdb/r-package")
# In some cases a restart of R is needed to work with a updated version of RSDB package (in RStudio - Session - Terminate R).

# logging into the server
fileName <- r'{/home/sdobelma/.config/rsdb_credentials.txt}' # file containing your username and pw in the form "username:password"
userpwd <- trimws(readChar(fileName, file.info(fileName)$size)) #  read account from file
remotesensing <-RSDB::RemoteSensing$new("https://gislab.hawk.de",userpwd)


# 02: metrics calculation
#-------------------------------------------------------------------------------

# list point cloud layers
remotesensing$pointclouds
remotesensing$rasterdbs

# set of metrics to calculate
metrics <- c(
  'point_density',
  'BE_H_KURTOSIS',
  'BE_H_SKEW',
  'BE_H_VAR',
  'BE_H_SD',
  'BE_H_MAX',
  'BE_H_MEAN',
  'BE_H_P10',
  'BE_H_P20',
  'BE_H_P50',
  'BE_H_P80',
  'BE_H_P90',
  'BE_H_P95',
  'BE_PR_02',
  'BE_PR_05',
  'BE_PR_10',
  'BE_PR_20',
  'BE_PR_30',
  'BE_RD_02',
  'BE_RD_05',
  'BE_RD_10',
  'BE_RD_20',
  'BE_RD_30',
  'pulse_returns_mean'
)

# forest attributes
forest_attributes <- c(
  'chm_height_mean',
  'ENL0',
  'chm_surface_area',
  'dtm_surface_area',
  'BE_FHD',
  'vegetation_coverage_05m_CHM',
  'vegetation_point_density'
)

## define seasons and densities to build pointcloud names 
seasons <- c("solling23_lon", "solling24_loff")
ppms <- c("ppm20", "ppm10", "ppm4")

## function for the task  
submit_index_task <- function(season,
                              ppm,
                              metric_set = c("indices", "forest_attributes"),
                              scanangle_filter = TRUE,
                              scanangle = "lt20",
                              debugged = TRUE,
                              task_pointcloud = "index_raster") {
  
  metric_set <- match.arg(metric_set)
  selected_metrics <- if (metric_set == "indices") metrics else forest_attributes
  
  # build pointcloud name
  pointcloud_name <- paste(season, ppm, sep = "_")
  
  if (scanangle_filter) {
    pointcloud_name <- paste(pointcloud_name, scanangle, sep = "_")
  }
  
  if (debugged) {
    pointcloud_name <- paste(pointcloud_name, "debugged", sep = "_")
  }
  
  # build rasterdb name
  rasterdb_name <- paste(pointcloud_name, metric_set, sep = "_")
  
  # step 1: rebuild CLMS layer
  remotesensing$submit_task(task = list(
    task_rasterdb = "rebuild",
    rasterdb = "CLMS_DLT2023_clipped"
  ))
  
  # step 2: rename to rasterdb_name
  remotesensing$submit_task(task = list(
    task_rasterdb = "rename",
    rasterdb = "CLMS_DLT2023_clipped_rebuild",
    new_name = rasterdb_name
  ))
  
  # warte bis RasterDB existiert
  start <- Sys.time()
  repeat {
    existing <- remotesensing$rasterdbs$name
    if (rasterdb_name %in% existing) break
    if (as.numeric(Sys.time() - start) > 60) stop("Timeout: RasterDB not found: ", rasterdb_name)
    Sys.sleep(3)
  }
  cat("RasterDB ready:", rasterdb_name, "\n")
  
  # step 3: submit index task
  result <- remotesensing$submit_task(task = list(
    task_pointcloud = task_pointcloud,
    pointcloud = pointcloud_name,
    rasterdb = rasterdb_name,
    indices = selected_metrics
  ))
  
  return(result)
}

## save results in a list
results <- list()

for (season in seasons) {
  for (ppm in ppms) {
    name <- paste(season, ppm, sep = "_")
    
    results[[name]] <- submit_index_task(
      season = season,
      ppm = ppm,
      metric_set = "indices", # "indices",  # oder "forest_attributes"
      scanangle_filter = TRUE,
      scanangle = "lt20"
    )
  }
}

################################################################################