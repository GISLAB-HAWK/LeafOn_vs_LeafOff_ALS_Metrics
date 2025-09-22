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



# source setup script
source('src/setup.R', local = TRUE)


# 01 - data reading
#-------------------------------------------------------------------------------

# input paths
metrics_path <- file.path(processed_data_dir, 'metrics')

# read w2w metrics (leaf- on & leaf-off)
metrics_w2w_leafon <- terra::rast(
  file.path(metrics_path, 'Solling23_leafon_indices_harm.tiff')
)

metrics_w2w_leafoff <- terra::rast(
  file.path(metrics_path, 'Solling24_leafoff_indices_harm.tiff')
)

metrics_w2w_leafon
metrics_w2w_leafoff
metrics_w2w_leafon@pntr$names
metrics_w2w_leafoff@pntr$names



# 02 - wall-to-wall modeling
#-------------------------------------

# path to the models
model_path <- file.path(processed_data_dir, 'models')

# output path for the predictions
vol_ha_pred_path <- file.path(processed_data_dir, 'predictions')

# leaf-on model --> predicted on leaf-on data
if (!file.exists(file.path(vol_ha_pred_path, 'vol_ha_pred_leafon.tif'))) {
  
  # deciduous model
  if (metrics_w2w_leafon$band1 == 1) {
    
    # read leaf-on deciduous model
    rf_model_lon_decid <- readRDS(
      file.path(model_path, 'ffs_rf_model_leafon_deciduous_filtered')
      )
    
    # predict for deciduous dominated pixels
    vol_ha_pred_leafon <- terra::predict(
      metrics_w2w_leafon,
      rf_model_lon_decid,
      na.rm = T
    )
  
  # coniferous model
  } else {
      
    # read leaf-on coniferous model
    rf_model_lon_conif <- readRDS(
      file.path(model_path, 'ffs_rf_model_leafon_coniferous_filtered')
    )
    
    # predict for coniferous dominated pixels
    vol_ha_pred_leafon <- terra::predict(
      metrics_w2w_leafon,
      rf_model_lon_conif,
      na.rm = T
    )
  }
}
  
### workaround because of different variable names ###
# Get the current names
current_names <- names(metrics_w2w_leafon)
# Remove "lon_" prefix if it exists
new_names <- gsub("^lon_", "", current_names)
# Rename the raster bands
names(metrics_w2w_leafon) <- new_names
###



# leaf-on model --> predicted on leaf-on data
if (!file.exists(file.path(vol_ha_pred_path, 'vol_ha_pred_leafon.tif'))) {
  
  # read deciduous and coniferous models
  rf_model_lon_decid <- readRDS(
    file.path(model_path, 'ffs_rf_model_leafon_deciduous_filtered.RDS')
  )
  
  rf_model_lon_conif <- readRDS(
    file.path(model_path, 'ffs_rf_model_leafon_coniferous_filtered.RDS')
  )
  
  # create masks for deciduous and coniferous pixels
  deciduous_mask <- metrics_w2w_leafon$band1 == 1
  coniferous_mask <- metrics_w2w_leafon$band1 != 1
  
  # predict for deciduous pixels
  vol_ha_pred_lon_decid <- terra::predict(
    metrics_w2w_leafon,
    rf_model_lon_decid,
    na.rm = T
  )
  
  # predict for coniferous pixels
  vol_ha_pred_lon_conif <- terra::predict(
    metrics_w2w_leafon,
    rf_model_lon_conif,
    na.rm = T
  )
  
  # combine predictions based on tree species
  vol_ha_pred_lon <- terra::ifel(
    deciduous_mask,
    vol_ha_pred_lon_decid,
    vol_ha_pred_lon_conif
  )
  
  # clean up temporary rasters
  rm(vol_ha_pred_lon_decid, vol_ha_pred_lon_conif, deciduous_mask, coniferous_mask)
  
}















# 02 - wall-to-wall modeling
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




















