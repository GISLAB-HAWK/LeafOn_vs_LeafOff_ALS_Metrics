#-------------------------------------------------------------------------------
# Name:         model_comp_remeasured_plots.R
# Description:  Script compares predictions of growing stock volume (GSV) (m³/ha),
#               once using a model trained on metrics calculated in plots 
#               measured with Real-Time Kinematic (RTK) and once using a model
#               trained on metrics not calculated from these RTK plots.
#               Two models are trained,
#               one using the the metrics calculated from leaf-on dataset,
#               and one using the metrics calculated from leaf-off dataset.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read data with inventory plots (BI plots RTK- and non RTK -based) 
# and calculated metrics (leaf-on and leaf-off)
plot_metrics <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'vol_stp_metrics.gpkg')
)

plot_metrics_rtk <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'vol_stp_metrics_rtk.gpkg')
)

head(plot_metrics)
str(plot_metrics)
head(plot_metrics_rtk)
str(plot_metrics_rtk)



# 02: data preparation
#-------------------------------------------------------------------------------

# remove NA (empty plots)
plot_metrics <- na.omit(plot_metrics)
plot_metrics_rtk <- na.omit(plot_metrics_rtk)

# filter plots to only include those which were remeasured
plot_metrics_rtk <- plot_metrics_rtk[plot_metrics_rtk$remeasured == 'yes', ]
remeasured_kspnr <- plot_metrics_rtk$kspnr
plot_metrics <- plot_metrics[plot_metrics$kspnr %in% remeasured_kspnr, ]

# split data for leaf-on and leaf-off
id_cols <- c('key', 'kspnr', 'abt', 'vol_ha', 'ts')

# select identifier columns + leaf-on metrics
plot_metrics_lon <- dplyr::select(
  plot_metrics,
  dplyr::any_of(id_cols),
  dplyr::starts_with('lon_')
)
plot_metrics_rtk_lon <- dplyr::select(
  plot_metrics_rtk,
  dplyr::any_of(id_cols),
  dplyr::starts_with('lon_')
)

# select identifier columns + leaf-off metrics
plot_metrics_loff <- dplyr::select(
  plot_metrics,
  dplyr::any_of(id_cols),
  dplyr::starts_with('loff_')
)
plot_metrics_rtk_loff <- dplyr::select(
  plot_metrics_rtk,
  dplyr::any_of(id_cols),
  dplyr::starts_with('loff_')
)

# remove prefix lon/loff
names(plot_metrics_lon) <- sub('^lon_', '', names(plot_metrics_lon))
names(plot_metrics_rtk_lon) <- sub('^lon_', '', names(plot_metrics_rtk_lon))
names(plot_metrics_loff) <- sub('^loff_', '', names(plot_metrics_loff))
names(plot_metrics_rtk_loff) <- sub('^loff_', '', names(plot_metrics_rtk_loff))



# 03 - model preparation
#-------------------------------------------------------------------------------

# define predictors and response
set.seed(11)
predictors_lon <- sf::st_drop_geometry(plot_metrics_lon[,4:length(plot_metrics_lon)])
predictors_rtk_lon <- sf::st_drop_geometry(plot_metrics_rtk_lon[,5:length(plot_metrics_rtk_lon)])
predictors_loff <- sf::st_drop_geometry(plot_metrics_loff[,4:length(plot_metrics_loff)])
predictors_rtk_loff <- sf::st_drop_geometry(plot_metrics_rtk_loff[,5:length(plot_metrics_rtk_loff)])
response <- sf::st_drop_geometry(plot_metrics_lon[,'vol_ha'])

# initialize cross validation - using k-fold instead of LOOCV to avoid issues
ctrl <- caret::trainControl(
  method = 'cv',
  number = 5,
  savePredictions = T,
  allowParallel = T
  )

# create grid for tuning features
tgrid <- expand.grid(
  mtry = 1:length(predictors_lon)
)



# 04 - model training
#-------------------------------------------------------------------------------

# train random forest models
# implementing forward feature selection

### leaf-on not using RTK ###
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_lon <- CAST::ffs(
  predictors_lon,
  response$vol_ha,
  method = 'rf',
  trControl = ctrl,
 tuneGrid = tgrid,
  ntree = 100,
  importance = F,
  seed = 999
)

parallel::stopCluster(cl)

### leaf-on using RTK ###
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_rtk_lon <- CAST::ffs(
  predictors_rtk_lon,
  response$vol_ha,
  method = 'rf',
  trControl = ctrl,
  tuneGrid = tgrid,
  ntree = 100,
  importance = F,
  seed = 999
)

parallel::stopCluster(cl)

### leaf-off not using RTK ###
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_loff <- CAST::ffs(
  predictors_loff,
  response$vol_ha,
  method = 'rf',
  trControl = ctrl,
  tuneGrid = tgrid,
  ntree = 100,
  importance = F,
  seed = 999
)

parallel::stopCluster(cl)

### leaf-off using RTK ###
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_rtk_loff <- CAST::ffs(
  predictors_rtk_loff,
  response$vol_ha,
  method = 'rf',
  trControl = ctrl,
  tuneGrid = tgrid,
  ntree = 100,
  importance = F,
  seed = 999
)

parallel::stopCluster(cl)

# table with results
rbind(
  data.frame(outcome='GSV', model='leaf-on no RTK',
             t(as.data.frame(CAST::global_validation(ffs_rf_model_lon)))),
  data.frame(outcome='GSV', model='leaf-on with RTK',
             t(as.data.frame(CAST::global_validation(ffs_rf_model_rtk_lon)))),
  data.frame(outcome='GSV', model='leaf-off no RTK',
             t(as.data.frame(CAST::global_validation(ffs_rf_model_loff)))),
  data.frame(outcome='GSV', model='leaf-off with RTK',
             t(as.data.frame(CAST::global_validation(ffs_rf_model_rtk_loff))))
) |> 
  knitr::kable(digits=2, row.names = F)



# 05 - plotting predicted vs observed
#-------------------------------------------------------------------------------

# extract cross-validation predictions from the trained models
cv_pred_lon <- ffs_rf_model_lon$pred
cv_pred_rtk_lon <- ffs_rf_model_rtk_lon$pred
cv_pred_loff <- ffs_rf_model_loff$pred
cv_pred_rtk_loff <- ffs_rf_model_rtk_loff$pred

# the cv predictions contain multiple rows per observation due to resampling
# get final predictions for each observation
final_pred_lon <- cv_pred_lon %>%
  dplyr::filter(
    mtry == ffs_rf_model_lon$bestTune$mtry) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_rtk_lon <- cv_pred_rtk_lon %>%
  dplyr::filter(
    mtry == ffs_rf_model_rtk_lon$bestTune$mtry) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_loff <- cv_pred_loff %>%
  dplyr::filter(
    mtry == ffs_rf_model_loff$bestTune$mtry) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_rtk_loff <- cv_pred_rtk_loff %>%
  dplyr::filter(
    mtry == ffs_rf_model_rtk_loff$bestTune$mtry) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

# extract validation metrics for annotations
val_lon <- CAST::global_validation(ffs_rf_model_lon)
val_rtk_lon <- CAST::global_validation(ffs_rf_model_rtk_lon)
val_loff <- CAST::global_validation(ffs_rf_model_loff)
val_rtk_loff <- CAST::global_validation(ffs_rf_model_rtk_loff)

# create individual plots for each model
p1 <- ggplot(final_pred_lon, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-on (no RTK)') +
  annotate('text', x = 50, y = 750, 
           label = paste0('RMSE = ', round(val_lon[1], 1), '\nR² = ', round(val_lon[2], 2)), 
           hjust = 0, vjust = 1, size = 3.5, fontface = 'bold')

p2 <- ggplot(final_pred_rtk_lon, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-on (with RTK)') +
  annotate('text', x = 50, y = 750, 
           label = paste0('RMSE = ', round(val_rtk_lon[1], 1), '\nR² = ', round(val_rtk_lon[2], 2)), 
           hjust = 0, vjust = 1, size = 3.5, fontface = 'bold')

p3 <- ggplot(final_pred_loff, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-off (no RTK)') +
  annotate('text', x = 50, y = 750, 
           label = paste0('RMSE = ', round(val_loff[1], 1), '\nR² = ', round(val_loff[2], 2)), 
           hjust = 0, vjust = 1, size = 3.5, fontface = 'bold')

p4 <- ggplot(final_pred_rtk_loff, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,800), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-off (with RTK)') +
  annotate('text', x = 50, y = 750, 
           label = paste0('RMSE = ', round(val_rtk_loff[1], 1), '\nR² = ', round(val_rtk_loff[2], 2)), 
           hjust = 0, vjust = 1, size = 3.5, fontface = 'bold')

combined_plot <- gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
print(combined_plot)



# 06 - coordinate deviation analysis
#-------------------------------------------------------------------------------

# extract coordinates from both datasets
coords_no_rtk <- sf::st_coordinates(plot_metrics)
coords_rtk <- sf::st_coordinates(plot_metrics_rtk)

# create dataframes with coordinates and plot identifiers
coord_df_no_rtk <- data.frame(
  kspnr = plot_metrics$kspnr,
  x_no_rtk = coords_no_rtk[,1],
  y_no_rtk = coords_no_rtk[,2]
)

coord_df_rtk <- data.frame(
  kspnr = plot_metrics_rtk$kspnr,
  x_rtk = coords_rtk[,1],
  y_rtk = coords_rtk[,2]
)

# merge datasets by plot identifier to ensure matching plots
coord_deviation_df <- merge(coord_df_no_rtk, coord_df_rtk, by = 'kspnr')

# calculate coordinate deviations
coord_deviation_df$x_deviation <- coord_deviation_df$x_rtk - coord_deviation_df$x_no_rtk
coord_deviation_df$y_deviation <- coord_deviation_df$y_rtk - coord_deviation_df$y_no_rtk

# calculate euclidean distance deviation
coord_deviation_df$euclidean_deviation <- sqrt(
  coord_deviation_df$x_deviation^2 + coord_deviation_df$y_deviation^2
)

# summary statistics
cat('Coordinate deviation summary:\n')
cat('X deviation (m): mean =', round(mean(coord_deviation_df$x_deviation), 2), 
    ', sd =', round(sd(coord_deviation_df$x_deviation), 2), '\n')
cat('Y deviation (m): mean =', round(mean(coord_deviation_df$y_deviation), 2), 
    ', sd =', round(sd(coord_deviation_df$y_deviation), 2), '\n')
cat('Euclidean deviation (m): mean =', round(mean(coord_deviation_df$euclidean_deviation), 2), 
    ', sd =', round(sd(coord_deviation_df$euclidean_deviation), 2), '\n')
cat('Max euclidean deviation (m):', round(max(coord_deviation_df$euclidean_deviation), 2), '\n')

# create scatter plot of coordinate deviations
p_coord_dev <- ggplot(coord_deviation_df, aes(x = x_deviation, y = y_deviation)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'red', linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'red', linewidth = 0.8) +
  coord_fixed(ratio = 1) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = 'bold'),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  xlab('X coordinate deviation (RTK - no RTK) [m]') +
  ylab('Y coordinate deviation (RTK - no RTK) [m]') +
  ggtitle('Coordinate Deviations: RTK vs non-RTK Plot Positions') +
  annotate('text', x = max(coord_deviation_df$x_deviation), 
           y = max(coord_deviation_df$y_deviation), 
           label = paste0('n = ', nrow(coord_deviation_df), ' plots\n',
                         'Mean euclidean dev. = ', 
                         round(mean(coord_deviation_df$euclidean_deviation), 2), ' m\n',
                         'Max euclidean dev. = ', 
                         round(max(coord_deviation_df$euclidean_deviation), 2), ' m'), 
           hjust = 1, vjust = 1, size = 3.5, 
           bbox = list(boxstyle = 'round,pad=0.3', facecolor = 'white', alpha = 0.8))

print(p_coord_dev)








































