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
#-------------------------------------------------------------------------------

# path to the models
model_path <- file.path(processed_data_dir, 'models')

# output path for the predictions
vol_ha_pred_path <- file.path(processed_data_dir, 'predictions')

################################################################################

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
  
  # write final prediction raster
  terra::writeRaster(
    vol_ha_pred_lon,
    file.path(vol_ha_pred_path, 'vol_ha_pred_leafon.tif'),
    overwrite = T
  )
  
  # clean up temporary rasters
  rm(vol_ha_pred_lon_decid, vol_ha_pred_lon_conif, deciduous_mask, coniferous_mask)
  
}

################################################################################

# leaf-off model --> predicted on leaf-off data
if (!file.exists(file.path(vol_ha_pred_path, 'vol_ha_pred_leafoff.tif'))) {
  
  # read deciduous and coniferous models
  rf_model_loff_decid <- readRDS(
    file.path(model_path, 'ffs_rf_model_leafoff_deciduous_filtered.RDS')
  )
  
  rf_model_loff_conif <- readRDS(
    file.path(model_path, 'ffs_rf_model_leafoff_coniferous_filtered.RDS')
  )
  
  # create masks for deciduous and coniferous pixels
  deciduous_mask <- metrics_w2w_leafoff$band1 == 1
  coniferous_mask <- metrics_w2w_leafoff$band1 != 1
  
  # predict for deciduous pixels
  vol_ha_pred_loff_decid <- terra::predict(
    metrics_w2w_leafoff,
    rf_model_loff_decid,
    na.rm = T
  )
  
  # predict for coniferous pixels
  vol_ha_pred_loff_conif <- terra::predict(
    metrics_w2w_leafoff,
    rf_model_loff_conif,
    na.rm = T
  )
  
  # combine predictions based on tree species
  vol_ha_pred_loff <- terra::ifel(
    deciduous_mask,
    vol_ha_pred_loff_decid,
    vol_ha_pred_loff_conif
  )
  
  # write final prediction raster
  terra::writeRaster(
    vol_ha_pred_loff,
    file.path(vol_ha_pred_path, 'vol_ha_pred_leafoff.tif'),
    overwrite = T
  )
  
  # clean up temporary rasters
  rm(vol_ha_pred_loff_decid, vol_ha_pred_loff_conif, deciduous_mask, coniferous_mask)
  
}

# crop leaf-on predictions to leaf-off prediction
vol_ha_pred_lon_cropped <- terra::crop(
  vol_ha_pred_lon,
  vol_ha_pred_loff,
  mask = T
)

terra::writeRaster(
  vol_ha_pred_lon_cropped,
  file.path(vol_ha_pred_path, 'vol_ha_pred_leafon.tif'),
  overwrite = T
)

# quick visualization
par(mfrow = c(1,2))
terra::plot(
  vol_ha_pred_lon_cropped,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-on model - leaf-on prediction'
)

terra::plot(
  vol_ha_pred_loff,
  col = grDevices::hcl.colors(
    n = 50, palette = 'YlGn', rev = T
  ),
  main = 'leaf-off model - leaf-off prediction'
)

################################################################################


# 03 - validation
#-------------------------------------------------------------------------------

# read test plots
test_lon_deciduous <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon_deciduous_filtered.RDS')
)
test_lon_coniferous <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon_coniferous_filtered.RDS')
)
test_loff_deciduous <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff_deciduous_filtered.RDS')
)
test_loff_coniferous <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff_coniferous_filtered.RDS')
)

# combine test datasets for leaf-on and leaf-off
test_lon <- rbind(test_lon_deciduous, test_lon_coniferous)
test_lon <- sf::st_drop_geometry(test_lon)
test_loff <- rbind(test_loff_deciduous, test_loff_coniferous)
test_loff <- sf::st_drop_geometry(test_loff)

################################################################################
# leaf-on --> leaf-on
################################################################################

# initialize vectors to store predictions
pred_values_lon <- numeric(nrow(test_lon))

# predict for deciduous plots
deciduous_mask_lon <- test_lon$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_lon) > 0) {
  pred_values_lon[deciduous_mask_lon] <- stats::predict(
    rf_model_lon_decid, 
    test_lon[deciduous_mask_lon, ]
  )
}

# predict for coniferous plots
coniferous_mask_lon <- test_lon$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_lon) > 0) {
  pred_values_lon[coniferous_mask_lon] <- stats::predict(
    rf_model_lon_conif, 
    test_lon[coniferous_mask_lon, ]
  )
}

# create validation dataframe
val_df_lon <- data.frame(
  kspnr = test_lon$kspnr,
  observed = test_lon$vol_ha,
  predicted = pred_values_lon,
  dominant_leaf_type = test_lon$dominant_leaf_type
)

# remove any NA values
# val_df <- val_df_lon[!is.na(val_df_lon$predicted) & !is.na(val_df_lon$observed), ]

# calculate validation metrics
rmse_lon <- round(sqrt(mean((val_df_lon$predicted - val_df_lon$observed)^2, na.rm = T)), 2)
rel_rmse_lon <- round((rmse_lon / mean(val_df_lon$observed)) * 100, 2)
bias_lon <- round(mean(val_df_lon$predicted - val_df_lon$observed), 2)
rel_bias_lon <- round((bias_lon / mean(val_df_lon$observed)) * 100, 2)
cat("Overall RMSE leaf-on:", rmse_lon, "\n")
cat("Overall relative RMSE leaf-on:", rel_rmse_lon, "%\n")
cat("Overall bias leaf-on:", bias_lon, "\n")
cat("Overall relative bias leaf-on:", rel_bias_lon, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_lon$dominant_leaf_type)) {
  subset_data <- val_df_lon[val_df_lon$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-on: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-on: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-on: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-on: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-off --> leaf-off
################################################################################

# initialize vectors to store predictions
pred_values_loff <- numeric(nrow(test_loff))

# predict for deciduous plots
deciduous_mask_loff <- test_loff$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_loff) > 0) {
  pred_values_loff[deciduous_mask_loff] <- stats::predict(
    rf_model_loff_decid, 
    test_loff[deciduous_mask_loff, ]
  )
}

# predict for coniferous plots
coniferous_mask_loff <- test_loff$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_loff) > 0) {
  pred_values_loff[coniferous_mask_loff] <- stats::predict(
    rf_model_loff_conif, 
    test_loff[coniferous_mask_loff, ]
  )
}

# create validation dataframe
val_df_loff <- data.frame(
  kspnr = test_loff$kspnr,
  observed = test_loff$vol_ha,
  predicted = pred_values_loff,
  dominant_leaf_type = test_loff$dominant_leaf_type
)

# remove any NA values
#val_df_loff <- val_df_loff[!is.na(val_df_loff$predicted) & !is.na(val_df_loff$observed), ]

# calculate validation metrics
rmse_loff <- round(sqrt(mean((val_df_loff$predicted - val_df_loff$observed)^2, na.rm = T)), 2)
rel_rmse_loff <- round((rmse_loff / mean(val_df_loff$observed)) * 100, 2)
bias_loff <- round(mean(val_df_loff$predicted - val_df_loff$observed), 2)
rel_bias_loff <- round((bias_loff / mean(val_df_loff$observed)) * 100, 2)
cat("Overall RMSE leaf-off:", rmse_loff, "\n")
cat("Overall relative RMSE leaf-off:", rel_rmse_loff, "%\n")
cat("Overall bias leaf-off:", bias_loff, "\n")
cat("Overall relative bias leaf-off:", rel_bias_loff, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_loff$dominant_leaf_type)) {
  subset_data <- val_df_loff[val_df_loff$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-off: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-off: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-off: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-off: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-on --> leaf-off
################################################################################

# initialize vectors to store predictions
pred_values_lon_on_loff <- numeric(nrow(test_loff))

# predict for deciduous plots using leaf-on deciduous model
deciduous_mask_loff <- test_loff$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_loff) > 0) {
  pred_values_lon_on_loff[deciduous_mask_loff] <- stats::predict(
    rf_model_lon_decid, 
    test_loff[deciduous_mask_loff, ]
  )
}

# predict for coniferous plots using leaf-on coniferous model
coniferous_mask_loff <- test_loff$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_loff) > 0) {
  pred_values_lon_on_loff[coniferous_mask_loff] <- stats::predict(
    rf_model_lon_conif, 
    test_loff[coniferous_mask_loff, ]
  )
}

# create validation dataframe
val_df_lon_on_loff <- data.frame(
  kspnr = test_loff$kspnr,
  observed = test_loff$vol_ha,
  predicted = pred_values_lon_on_loff,
  dominant_leaf_type = test_loff$dominant_leaf_type
)

# calculate validation metrics
rmse_lon_on_loff <- round(sqrt(mean((val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed)^2, na.rm = T)), 2)
rel_rmse_lon_on_loff <- round((rmse_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)
bias_lon_on_loff <- round(mean(val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed), 2)
rel_bias_lon_on_loff <- round((bias_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)

cat("Overall RMSE leaf-on model on leaf-off data:", rmse_lon_on_loff, "\n")
cat("Overall relative RMSE leaf-on model on leaf-off data:", rel_rmse_lon_on_loff, "%\n")
cat("Overall bias leaf-on model on leaf-off data:", bias_lon_on_loff, "\n")
cat("Overall relative bias leaf-on model on leaf-off data:", rel_bias_lon_on_loff, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_lon_on_loff$dominant_leaf_type)) {
  subset_data <- val_df_lon_on_loff[val_df_lon_on_loff$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-on model on leaf-off data: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-on model on leaf-off data: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-on model on leaf-off data: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-on model on leaf-off data: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-off --> leaf-on
################################################################################

# initialize vectors to store predictions
pred_values_loff_on_lon <- numeric(nrow(test_lon))

# predict for deciduous plots using leaf-off deciduous model
deciduous_mask_lon <- test_lon$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_lon) > 0) {
  pred_values_loff_on_lon[deciduous_mask_lon] <- stats::predict(
    rf_model_loff_decid, 
    test_lon[deciduous_mask_lon, ]
  )
}

# predict for coniferous plots using leaf-off coniferous model
coniferous_mask_lon <- test_lon$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_lon) > 0) {
  pred_values_loff_on_lon[coniferous_mask_lon] <- stats::predict(
    rf_model_loff_conif, 
    test_lon[coniferous_mask_lon, ]
  )
}

# create validation dataframe
val_df_loff_on_lon <- data.frame(
  kspnr = test_lon$kspnr,
  observed = test_lon$vol_ha,
  predicted = pred_values_loff_on_lon,
  dominant_leaf_type = test_lon$dominant_leaf_type
)

# calculate validation metrics
rmse_loff_on_lon <- round(sqrt(mean((val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed)^2, na.rm = T)), 2)
rel_rmse_loff_on_lon <- round((rmse_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)
bias_loff_on_lon <- round(mean(val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed), 2)
rel_bias_loff_on_lon <- round((bias_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)

cat("Overall RMSE leaf-off model on leaf-on data:", rmse_loff_on_lon, "\n")
cat("Overall relative RMSE leaf-off model on leaf-on data:", rel_rmse_loff_on_lon, "%\n")
cat("Overall bias leaf-off model on leaf-on data:", bias_loff_on_lon, "\n")
cat("Overall relative bias leaf-off model on leaf-on data:", rel_bias_loff_on_lon, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_loff_on_lon$dominant_leaf_type)) {
  subset_data <- val_df_loff_on_lon[val_df_loff_on_lon$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-off model on leaf-on data: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-off model on leaf-on data: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-off model on leaf-on data: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-off model on leaf-on data: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# final data frame with error metrics
################################################################################

calc_metrics <- function(df, dataset_name, leaf_type_col = "dominant_leaf_type") {
  
  # overall metrics
  overall <- data.frame(
    dataset = dataset_name,
    leaf_type = "overall",
    RMSE = round(sqrt(mean((df$predicted - df$observed)^2, na.rm = TRUE)), 2),
    Rel_RMSE = round((sqrt(mean((df$predicted - df$observed)^2, na.rm = TRUE)) / mean(df$observed, na.rm = TRUE)) * 100, 2),
    Bias = round(mean(df$predicted - df$observed, na.rm = TRUE), 2),
    Rel_Bias = round((mean(df$predicted - df$observed, na.rm = TRUE) / mean(df$observed, na.rm = TRUE)) * 100, 2)
  )
  
  # metrics by leaf type (only if leaf_type_col is not NULL)
  if (!is.null(leaf_type_col) && leaf_type_col %in% colnames(df)) {
    by_leaf <- do.call(rbind, lapply(unique(df[[leaf_type_col]]), function(lt) {
      sub <- df[df[[leaf_type_col]] == lt, ]
      data.frame(
        dataset = dataset_name,
        leaf_type = lt,
        RMSE = round(sqrt(mean((sub$predicted - sub$observed)^2, na.rm = TRUE)), 2),
        Rel_RMSE = round((sqrt(mean((sub$predicted - sub$observed)^2, na.rm = TRUE)) / mean(sub$observed, na.rm = TRUE)) * 100, 2),
        Bias = round(mean(sub$predicted - sub$observed, na.rm = TRUE), 2),
        Rel_Bias = round((mean(sub$predicted - sub$observed, na.rm = TRUE) / mean(sub$observed, na.rm = TRUE)) * 100, 2)
      )
    }))
    rbind(overall, by_leaf)
  } else {
    overall
  }
}

metrics_all <- rbind(
  calc_metrics(val_df_lon, "Leaf-on model on leaf-on data"),
  calc_metrics(val_df_loff, "Leaf-off model on leaf-off data"),
  calc_metrics(val_df_lon_on_loff, "Leaf-on model on leaf-off data"),
  calc_metrics(val_df_loff_on_lon, "Leaf-off model on leaf-on data")
)

################################################################################
# plots predicted vs. observed growing stock
################################################################################

plot_leafon <- ggplot(val_df_lon, aes(x=observed, y=predicted)) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red') #+
  #ggtitle('leaf-on model | leaf-on data')

plot_leafoff <- ggplot(val_df_loff, aes(x=observed, y=predicted)) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red') #+
  #ggtitle('leaf-off model | leaf-off data')

plot_leafon_leafoff <- ggplot(val_df_lon_on_loff, aes(x=observed, y=predicted)) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red')# +
  #ggtitle('leaf-on model | leaf-off data')

plot_leafoff_leafon <- ggplot(val_df_loff_on_lon, aes(x=observed, y=predicted)) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14)
  ) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red') #+
  #ggtitle('leaf-off model | leaf-on data')

ggpubr::ggarrange(
  plot_leafon, plot_leafoff,
  plot_leafon_leafoff, plot_leafoff_leafon,
  ncol = 2, nrow = 2)


################################################################################


# 04 - validation only with remeasured plots
#-------------------------------------------------------------------------------

################################################################################
# leaf-on --> leaf-on
################################################################################

# filter only remeasured plots
test_lon_re <- subset(test_lon, remeasured == "yes")

# initialize vectors to store predictions
pred_values_lon <- numeric(nrow(test_lon_re))

# predict for deciduous plots
deciduous_mask_lon <- test_lon_re$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_lon) > 0) {
  pred_values_lon[deciduous_mask_lon] <- stats::predict(
    rf_model_lon_decid, 
    test_lon_re[deciduous_mask_lon, ]
  )
}

# predict for coniferous plots
coniferous_mask_lon <- test_lon_re$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_lon) > 0) {
  pred_values_lon[coniferous_mask_lon] <- stats::predict(
    rf_model_lon_conif, 
    test_lon_re[coniferous_mask_lon, ]
  )
}

# create validation dataframe
val_df_lon <- data.frame(
  kspnr = test_lon_re$kspnr,
  observed = test_lon_re$vol_ha,
  predicted = pred_values_lon,
  dominant_leaf_type = test_lon_re$dominant_leaf_type
)

# calculate validation metrics
rmse_lon <- round(sqrt(mean((val_df_lon$predicted - val_df_lon$observed)^2, na.rm = T)), 2)
rel_rmse_lon <- round((rmse_lon / mean(val_df_lon$observed)) * 100, 2)
bias_lon <- round(mean(val_df_lon$predicted - val_df_lon$observed), 2)
rel_bias_lon <- round((bias_lon / mean(val_df_lon$observed)) * 100, 2)
cat("Overall RMSE leaf-on:", rmse_lon, "\n")
cat("Overall relative RMSE leaf-on:", rel_rmse_lon, "%\n")
cat("Overall bias leaf-on:", bias_lon, "\n")
cat("Overall relative bias leaf-on:", rel_bias_lon, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_lon$dominant_leaf_type)) {
  subset_data <- val_df_lon[val_df_lon$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-on: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-on: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-on: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-on: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-off --> leaf-off
################################################################################

# filter only remeasured plots
test_loff_re <- subset(test_loff, remeasured == "yes")

# initialize vectors to store predictions
pred_values_loff <- numeric(nrow(test_loff_re))

# predict for deciduous plots
deciduous_mask_loff <- test_loff_re$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_loff) > 0) {
  pred_values_loff[deciduous_mask_loff] <- stats::predict(
    rf_model_loff_decid, 
    test_loff_re[deciduous_mask_loff, ]
  )
}

# predict for coniferous plots
coniferous_mask_loff <- test_loff_re$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_loff) > 0) {
  pred_values_loff[coniferous_mask_loff] <- stats::predict(
    rf_model_loff_conif, 
    test_loff_re[coniferous_mask_loff, ]
  )
}

# create validation dataframe
val_df_loff <- data.frame(
  kspnr = test_loff_re$kspnr,
  observed = test_loff_re$vol_ha,
  predicted = pred_values_loff,
  dominant_leaf_type = test_loff_re$dominant_leaf_type
)

# calculate validation metrics
rmse_loff <- round(sqrt(mean((val_df_loff$predicted - val_df_loff$observed)^2, na.rm = T)), 2)
rel_rmse_loff <- round((rmse_loff / mean(val_df_loff$observed)) * 100, 2)
bias_loff <- round(mean(val_df_loff$predicted - val_df_loff$observed), 2)
rel_bias_loff <- round((bias_loff / mean(val_df_loff$observed)) * 100, 2)
cat("Overall RMSE leaf-off:", rmse_loff, "\n")
cat("Overall relative RMSE leaf-off:", rel_rmse_loff, "%\n")
cat("Overall bias leaf-off:", bias_loff, "\n")
cat("Overall relative bias leaf-off:", rel_bias_loff, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_loff$dominant_leaf_type)) {
  subset_data <- val_df_loff[val_df_loff$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-off: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-off: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-off: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-off: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-on --> leaf-off
################################################################################

# initialize vectors to store predictions
pred_values_lon_on_loff <- numeric(nrow(test_loff_re))

# predict for deciduous plots using leaf-on deciduous model
deciduous_mask_loff <- test_loff_re$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_loff) > 0) {
  pred_values_lon_on_loff[deciduous_mask_loff] <- stats::predict(
    rf_model_lon_decid, 
    test_loff_re[deciduous_mask_loff, ]
  )
}

# predict for coniferous plots using leaf-on coniferous model
coniferous_mask_loff <- test_loff_re$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_loff) > 0) {
  pred_values_lon_on_loff[coniferous_mask_loff] <- stats::predict(
    rf_model_lon_conif, 
    test_loff_re[coniferous_mask_loff, ]
  )
}

# create validation dataframe
val_df_lon_on_loff <- data.frame(
  kspnr = test_loff_re$kspnr,
  observed = test_loff_re$vol_ha,
  predicted = pred_values_lon_on_loff,
  dominant_leaf_type = test_loff_re$dominant_leaf_type
)

# calculate validation metrics
rmse_lon_on_loff <- round(sqrt(mean((val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed)^2, na.rm = T)), 2)
rel_rmse_lon_on_loff <- round((rmse_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)
bias_lon_on_loff <- round(mean(val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed), 2)
rel_bias_lon_on_loff <- round((bias_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)

cat("Overall RMSE leaf-on model on leaf-off data:", rmse_lon_on_loff, "\n")
cat("Overall relative RMSE leaf-on model on leaf-off data:", rel_rmse_lon_on_loff, "%\n")
cat("Overall bias leaf-on model on leaf-off data:", bias_lon_on_loff, "\n")
cat("Overall relative bias leaf-on model on leaf-off data:", rel_bias_lon_on_loff, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_lon_on_loff$dominant_leaf_type)) {
  subset_data <- val_df_lon_on_loff[val_df_lon_on_loff$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-on model on leaf-off data: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-on model on leaf-off data: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-on model on leaf-off data: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-on model on leaf-off data: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# leaf-off --> leaf-on
################################################################################

# initialize vectors to store predictions
pred_values_loff_on_lon <- numeric(nrow(test_lon_re))

# predict for deciduous plots using leaf-off deciduous model
deciduous_mask_lon <- test_lon_re$dominant_leaf_type == 'deciduous'
if (sum(deciduous_mask_lon) > 0) {
  pred_values_loff_on_lon[deciduous_mask_lon] <- stats::predict(
    rf_model_loff_decid, 
    test_lon_re[deciduous_mask_lon, ]
  )
}

# predict for coniferous plots using leaf-off coniferous model
coniferous_mask_lon <- test_lon_re$dominant_leaf_type == 'coniferous'
if (sum(coniferous_mask_lon) > 0) {
  pred_values_loff_on_lon[coniferous_mask_lon] <- stats::predict(
    rf_model_loff_conif, 
    test_lon_re[coniferous_mask_lon, ]
  )
}

# create validation dataframe
val_df_loff_on_lon <- data.frame(
  kspnr = test_lon_re$kspnr,
  observed = test_lon_re$vol_ha,
  predicted = pred_values_loff_on_lon,
  dominant_leaf_type = test_lon_re$dominant_leaf_type
)

# calculate validation metrics
rmse_loff_on_lon <- round(sqrt(mean((val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed)^2, na.rm = T)), 2)
rel_rmse_loff_on_lon <- round((rmse_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)
bias_loff_on_lon <- round(mean(val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed), 2)
rel_bias_loff_on_lon <- round((bias_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)

cat("Overall RMSE leaf-off model on leaf-on data:", rmse_loff_on_lon, "\n")
cat("Overall relative RMSE leaf-off model on leaf-on data:", rel_rmse_loff_on_lon, "%\n")
cat("Overall bias leaf-off model on leaf-on data:", bias_loff_on_lon, "\n")
cat("Overall relative bias leaf-off model on leaf-on data:", rel_bias_loff_on_lon, "%\n")

# calculate metrics separately for each dominant leaf type
for (leaf_type in unique(val_df_loff_on_lon$dominant_leaf_type)) {
  subset_data <- val_df_loff_on_lon[val_df_loff_on_lon$dominant_leaf_type == leaf_type, ]
  
  if (nrow(subset_data) > 0) {
    rmse_subset <- round(sqrt(mean((subset_data$predicted - subset_data$observed)^2, na.rm = T)), 2)
    rel_rmse_subset <- round((rmse_subset / mean(subset_data$observed)) * 100, 2)
    bias_subset <- round(mean(subset_data$predicted - subset_data$observed), 2)
    rel_bias_subset <- round((bias_subset / mean(subset_data$observed)) * 100, 2)
    
    cat(paste0(leaf_type, " RMSE leaf-off model on leaf-on data: "), rmse_subset, "\n")
    cat(paste0(leaf_type, " relative RMSE leaf-off model on leaf-on data: "), rel_rmse_subset, "%\n")
    cat(paste0(leaf_type, " bias leaf-off model on leaf-on data: "), bias_subset, "\n")
    cat(paste0(leaf_type, " relative bias leaf-off model on leaf-on data: "), rel_bias_subset, "%\n")
  }
}

################################################################################
# final data frame with error metrics
################################################################################

metrics_remeasured <- rbind(
  calc_metrics(val_df_lon, "Leaf-on model on leaf-on data", leaf_type_col = NULL),
  calc_metrics(val_df_loff, "Leaf-off model on leaf-off data", leaf_type_col = NULL),
  calc_metrics(val_df_lon_on_loff, "Leaf-on model on leaf-off data", leaf_type_col = NULL),
  calc_metrics(val_df_loff_on_lon, "Leaf-off model on leaf-on data", leaf_type_col = NULL)
)


#############################################################################

# 05 - validation using the model only trained with RTK plots
#-------------------------------------------------------------------------------

# read test plots
test_lon <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon_filtered_rtk_only.RDS')
)
test_loff <- readRDS( 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff_filtered_rtk_only.RDS')
)
test_lon <- sf::st_drop_geometry(test_lon)
test_loff <- sf::st_drop_geometry(test_loff)

################################################################################
# leaf-on --> leaf-on
################################################################################

# predict for all plots with the leaf-on model
pred_values_lon <- stats::predict(rf_model_lon, test_lon)

# create validation dataframe
val_df_lon <- data.frame(
  kspnr = test_lon$kspnr,
  observed = test_lon$vol_ha,
  predicted = pred_values_lon
)

# calculate validation metrics
rmse_lon <- round(sqrt(mean((val_df_lon$predicted - val_df_lon$observed)^2, na.rm = TRUE)), 2)
rel_rmse_lon <- round((rmse_lon / mean(val_df_lon$observed)) * 100, 2)
bias_lon <- round(mean(val_df_lon$predicted - val_df_lon$observed, na.rm = TRUE), 2)
rel_bias_lon <- round((bias_lon / mean(val_df_lon$observed)) * 100, 2)
cat("Overall RMSE leaf-on:", rmse_lon, "\n")
cat("Overall relative RMSE leaf-on:", rel_rmse_lon, "%\n")
cat("Overall bias leaf-on:", bias_lon, "\n")
cat("Overall relative bias leaf-on:", rel_bias_lon, "%\n")

################################################################################
# leaf-off --> leaf-off
################################################################################

pred_values_loff <- stats::predict(rf_model_loff, test_loff)

val_df_loff <- data.frame(
  kspnr = test_loff$kspnr,
  observed = test_loff$vol_ha,
  predicted = pred_values_loff
)

rmse_loff <- round(sqrt(mean((val_df_loff$predicted - val_df_loff$observed)^2, na.rm = TRUE)), 2)
rel_rmse_loff <- round((rmse_loff / mean(val_df_loff$observed)) * 100, 2)
bias_loff <- round(mean(val_df_loff$predicted - val_df_loff$observed, na.rm = TRUE), 2)
rel_bias_loff <- round((bias_loff / mean(val_df_loff$observed)) * 100, 2)
cat("Overall RMSE leaf-off:", rmse_loff, "\n")
cat("Overall relative RMSE leaf-off:", rel_rmse_loff, "%\n")
cat("Overall bias leaf-off:", bias_loff, "\n")
cat("Overall relative bias leaf-off:", rel_bias_loff, "%\n")

################################################################################
# leaf-on --> leaf-off
################################################################################

pred_values_lon_on_loff <- stats::predict(rf_model_lon, test_loff)

val_df_lon_on_loff <- data.frame(
  kspnr = test_loff$kspnr,
  observed = test_loff$vol_ha,
  predicted = pred_values_lon_on_loff
)

rmse_lon_on_loff <- round(sqrt(mean((val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed)^2, na.rm = TRUE)), 2)
rel_rmse_lon_on_loff <- round((rmse_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)
bias_lon_on_loff <- round(mean(val_df_lon_on_loff$predicted - val_df_lon_on_loff$observed, na.rm = TRUE), 2)
rel_bias_lon_on_loff <- round((bias_lon_on_loff / mean(val_df_lon_on_loff$observed)) * 100, 2)
cat("Overall RMSE leaf-on model on leaf-off data:", rmse_lon_on_loff, "\n")
cat("Overall relative RMSE leaf-on model on leaf-off data:", rel_rmse_lon_on_loff, "%\n")
cat("Overall bias leaf-on model on leaf-off data:", bias_lon_on_loff, "\n")
cat("Overall relative bias leaf-on model on leaf-off data:", rel_bias_lon_on_loff, "%\n")

################################################################################
# leaf-off --> leaf-on
################################################################################

pred_values_loff_on_lon <- stats::predict(rf_model_loff, test_lon)

val_df_loff_on_lon <- data.frame(
  kspnr = test_lon$kspnr,
  observed = test_lon$vol_ha,
  predicted = pred_values_loff_on_lon
)

rmse_loff_on_lon <- round(sqrt(mean((val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed)^2, na.rm = TRUE)), 2)
rel_rmse_loff_on_lon <- round((rmse_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)
bias_loff_on_lon <- round(mean(val_df_loff_on_lon$predicted - val_df_loff_on_lon$observed, na.rm = TRUE), 2)
rel_bias_loff_on_lon <- round((bias_loff_on_lon / mean(val_df_loff_on_lon$observed)) * 100, 2)
cat("Overall RMSE leaf-off model on leaf-on data:", rmse_loff_on_lon, "\n")
cat("Overall relative RMSE leaf-off model on leaf-on data:", rel_rmse_loff_on_lon, "%\n")
cat("Overall bias leaf-off model on leaf-on data:", bias_loff_on_lon, "\n")
cat("Overall relative bias leaf-off model on leaf-on data:", rel_bias_loff_on_lon, "%\n")

################################################################################
# final data frame with error metrics
################################################################################

metrics_all <- rbind(
  calc_metrics(val_df_lon, "Leaf-on model on leaf-on data", leaf_type_col = NULL),
  calc_metrics(val_df_loff, "Leaf-off model on leaf-off data", leaf_type_col = NULL),
  calc_metrics(val_df_lon_on_loff, "Leaf-on model on leaf-off data", leaf_type_col = NULL),
  calc_metrics(val_df_loff_on_lon, "Leaf-off model on leaf-on data", leaf_type_col = NULL)
)

print(metrics_all)

################################################################################




























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




















