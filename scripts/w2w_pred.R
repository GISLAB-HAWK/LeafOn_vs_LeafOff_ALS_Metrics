#-----------------------------------------------------------------------------
# Name:         w2w_pred.R
# Description:  Script predicts the growing stock (m³/ha) based on previously
#               derived ALS-based metrics in forest inventory plots.
#               A random forest trained before with ALS leaf-on data is used.
#               The metrics derived on plot-level before are calculated on
#               pixel-level to generate wall-to-wall raster of these metrics. 
#               This is done twice, with leaf-on and leaf-off data.
#               Finally, wall-to-wall predictions of the growing stock
#               are generated, again twice with the pixel-metrics based on
#               leaf-on data and with pixel-metrics based on leaf-off data.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-----------------------------------------------------------------------------



# 01 - file path definitions
#----------------------------

# define raw data directory
raw_data_dir <- 'data/raw/'

# define processed data directory
processed_data_dir <- 'data/processed/'

# define output directory
output_dir <- 'output/'



# 02 - data reading
#-------------------------------------

# input path to point clouds
path_pc_leafoff <- file.path(raw_data_dir, 'pc_LeafOff_2024')
path_pc_leafon <- file.path(raw_data_dir, 'pc_LeafOn_2023')

# read point clouds with LAScatalog
pc_ctg_leafoff <- lidR::readLAScatalog(path_pc_leafoff)
pc_ctg_leafon <- lidR::readLAScatalog(path_pc_leafon)

pc_ctg_leafoff
pc_ctg_leafon

# assign CRS from pc_ctg_leafoff to pc_ctg_leafon
# ETRS89 / UTM zone 32N
lidR::crs(pc_ctg_leafon) <- lidR::crs(pc_ctg_leafoff)



# 03 - wall-to-wall modeling
#-------------------------------------

# 1. calculate forest metrics for the entire collection of files
# (normalized point clouds in LAScatalog)
# output resolution of the metrics = 20 m (13 m plot radius = 531 m²)
# this is done twice, for leaf-on and for leaf-off data
if (!file.exists(file.path(output_dir, 'metrics_w2w_leafon.tif'))) {
  
  lidR::opt_filter(pc_ctg_leafon) <- '-drop_z_below 2'
  
  metrics_w2w_leafon <- lidR::pixel_metrics(
    pc_ctg_leafon, .stdmetrics,
    res = 20,
    pkg = 'terra'
  )
  
  terra::writeRaster(
    metrics_w2w_leafon,
    file.path(output_dir, 'metrics_w2w_leafon.tif'),
    overwrite = T)
  
} else {
  
  metrics_w2w_leafon <- terra::rast(
    file.path(output_dir, 'metrics_w2w_leafon.tif')
    )
  
}

if (!file.exists(file.path(output_dir, 'metrics_w2w_leafoff.tif'))) {
  
  lidR::opt_filter(pc_ctg_leafoff) <- '-drop_z_below 2'
  
  metrics_w2w_leafoff <- lidR::pixel_metrics(
    pc_ctg_leafoff, .stdmetrics,
    res = 20,
    pkg = 'terra'
  )
  
  terra::writeRaster(
    metrics_w2w_leafoff,
    file.path(output_dir, 'metrics_w2w_leafoff.tif'),
    overwrite = T)
  
} else {
  
  metrics_w2w_leafoff <- terra::rast(
    file.path(output_dir, 'metrics_w2w_leafoff.tif')
  )
  
}

# 2. predict growing stock for the whole area using random forest model
# trained in script model_train.R with leaf-on data
# two predictions:
# one using the metrics calculated with leaf-on data,
# and the other one using the metrics calculated with leaf-off data

# leaf-on model --> predicted on leaf-on data
if (!file.exists(file.path(output_dir, 'vol_ha_pred_leafon.tif'))) {
  
  if (!exists('metrics_w2w_leafon.tif')) {
    metrics_w2w_leafon <- terra::rast(
      file.path(output_dir, 'metrics_w2w_leafon.tif')
      )
  }
  
  rf_model_leafon <- readRDS(file.path(output_dir, 'rf_model_leafon.RDS'))
  
  vol_ha_pred_leafon <- terra::predict(
    metrics_w2w_leafon,
    rf_model_leafon,
    na.rm = T
    )
  
  terra::writeRaster(
    vol_ha_pred_leafon,
    file.path(output_dir, 'vol_ha_pred_leafon.tif'),
    overwrite = T
    )
  
} else {
  
  vol_ha_pred_leafon <- terra::rast(
    file.path(output_dir, 'vol_ha_pred_leafon.tif')
    )
  
}

# leaf-off model --> predicted on leaf-off data
if (!file.exists(file.path(output_dir, 'vol_ha_pred_leafoff.tif'))) {
  
  if (!exists('metrics_w2w_leafoff.tif')) {
    metrics_w2w_leafoff <- terra::rast(
      file.path(output_dir, 'metrics_w2w_leafoff.tif')
    )
  }
  
  rf_model_leafoff <- readRDS(file.path(output_dir, 'rf_model_leafoff.RDS'))
  
  vol_ha_pred_leafoff <- terra::predict(
    metrics_w2w_leafoff,
    rf_model_leafoff,
    na.rm = T
  )
  
  terra::writeRaster(
    vol_ha_pred_leafoff,
    file.path(output_dir, 'vol_ha_pred_leafoff.tif'),
    overwrite = T
  )
  
} else {
  
  vol_ha_pred_leafoff <- terra::rast(
    file.path(output_dir, 'vol_ha_pred_leafoff.tif')
  )
  
}

# leaf-on model --> predicted on leaf-off data
if (!file.exists(file.path(output_dir, 'vol_ha_pred_leafon_leafoff.tif'))) {
  
  if (!exists('metrics_w2w_leafoff.tif')) {
    metrics_w2w_leafoff <- terra::rast(
      file.path(output_dir, 'metrics_w2w_leafoff.tif')
    )
  }
  
  rf_model_leafon <- readRDS(file.path(output_dir, 'rf_model_leafon.RDS'))
  
  vol_ha_pred_leafon_leafoff <- terra::predict(
    metrics_w2w_leafoff,
    rf_model_leafon,
    na.rm = T
  )
  
  terra::writeRaster(
    vol_ha_pred_leafon_leafoff,
    file.path(output_dir, 'vol_ha_pred_leafon_leafoff.tif'),
    overwrite = T
  )
  
} else {
  
  vol_ha_pred_leafon_leafoff <- terra::rast(
    file.path(output_dir, 'vol_ha_pred_leafon_leafoff.tif')
  )
  
}

# leaf-off model --> predicted on leaf-on data
if (!file.exists(file.path(output_dir, 'vol_ha_pred_leafoff_leafon.tif'))) {
  
  if (!exists('metrics_w2w_leafon.tif')) {
    metrics_w2w_leafon <- terra::rast(
      file.path(output_dir, 'metrics_w2w_leafon.tif')
    )
  }
  
  rf_model_leafoff <- readRDS(file.path(output_dir, 'rf_model_leafoff.RDS'))
  
  vol_ha_pred_leafoff_leafon <- terra::predict(
    metrics_w2w_leafon,
    rf_model_leafoff,
    na.rm = T
  )
  
  terra::writeRaster(
    vol_ha_pred_leafoff_leafon,
    file.path(output_dir, 'vol_ha_pred_leafoff_leafon.tif'),
    overwrite = T
  )
  
} else {
  
  vol_ha_pred_leafoff_leafon <- terra::rast(
    file.path(output_dir, 'vol_ha_pred_leafoff_leafon.tif')
  )
  
}

# crop leaf-on predictions to leaf-off prediction
vol_ha_pred_leafon_cropped <- terra::crop(
  vol_ha_pred_leafon,
  vol_ha_pred_leafoff,
  mask = T
)

vol_ha_pred_leafoff_leafon_cropped <- terra::crop(
  vol_ha_pred_leafoff_leafon,
  vol_ha_pred_leafon_leafoff,
  mask = T
)

# write to disk
terra::writeRaster(
  vol_ha_pred_leafon_cropped,
  file.path(output_dir, 'vol_ha_pred_leafon_cropped.tif')
  )

terra::writeRaster(
  vol_ha_pred_leafoff_leafon_cropped,
  file.path(output_dir, 'vol_ha_pred_leafoff_leafon_cropped.tif')
)

# quick visualization
par(mfrow = c(1,2))
terra::plot(
  vol_ha_pred_leafon_cropped,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-on model - leaf-on prediction'
)

terra::plot(
  vol_ha_pred_leafoff,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-off model - leaf-off prediction'
)

par(mfrow = c(1,2))
terra::plot(
  vol_ha_pred_leafoff_leafon_cropped,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-off model - leaf-on prediction'
)

terra::plot(
  vol_ha_pred_leafon_leafoff,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-on model - leaf-off prediction'
)

par(mfrow = c(1,4))
terra::plot(
  vol_ha_pred_leafon_cropped,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-on model - leaf-on prediction'
)

terra::plot(
  vol_ha_pred_leafoff,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-off model - leaf-off prediction'
)

terra::plot(
  vol_ha_pred_leafoff_leafon_cropped,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-off model - leaf-on prediction'
)

terra::plot(
  vol_ha_pred_leafon_leafoff,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-on model - leaf-off prediction'
)




















