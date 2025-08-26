#-------------------------------------------------------------------------------
# Name:         model_train.R
# Description:  Script trains random forest models for predicting
#               growing stock volume (GSV) (m³/ha).
#               ALS-based metrics previously derived in forest inventory plots
#               (BI plots) are used as predictors. Two models are trained,
#               one using the the metrics calculated from leaf-on dataset,
#               and one using the metrics calculated from leaf-off dataset.
#               Leave-Location-Out cross-validation (LLO CV) is used as a 
#               spatial cross validation method.
#               80% of each dataset is used for training the models.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read data with inventory plots (BI) 
# and calculated metrics (leaf-on and leaf-off) and extracted tree species
plot_metrics <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'vol_stp_metrics_rtk.gpkg')
  )

head(plot_metrics)
str(plot_metrics)

# read extent of the leaf-off dataset
ext_loff <- sf::st_read(
  file.path(processed_data_dir, 'pc_leafoff_2024', 'leafoff.vpc')
  )

# plot data
ggplot() +
  geom_sf(data = ext_loff, fill = "grey", alpha = 0.1) +
  geom_sf(data = plot_metrics, aes(col = vol_ha)) +
  scale_colour_distiller(palette = "YlGn", direction = 1) +
  theme_bw() +
  labs(col = "") +
  ggtitle("Growing Stock Volume (m³/ha)")



# 02: data preparation
#-------------------------------------------------------------------------------

# remove NA (empty plots)
plot_metrics <- na.omit(plot_metrics)

# split data for leaf-on and leaf-off
id_cols <- c('key', 'kspnr', 'abt', 'vol_ha', 'ts')

# select identifier columns + leaf-on metrics
plot_metrics_lon <- dplyr::select(
  plot_metrics,
  dplyr::any_of(id_cols),
  dplyr::starts_with('lon_')
)

# select identifier columns + leaf-off metrics
plot_metrics_loff <- dplyr::select(
  plot_metrics,
  dplyr::any_of(id_cols),
  dplyr::starts_with('loff_')
)

# df versions
plot_metrics_lon_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_lon))
plot_metrics_loff_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_loff))



# 03 - model preparation
#-------------------------------------------------------------------------------

################################################################################
#                           START DRAFT - cv testing                           #

# random 5-fold cross-validation (cv)
fold5 <- caret::createFolds(1:nrow(plot_metrics_loff), k = 5, returnTrain = F)

predcond_cv_fold5 <- CAST::geodist(
  x = plot_metrics_loff,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = fold5
  )

# leave-one-out cv (loocv)
loo_cv <- caret::createFolds(1:nrow(plot_metrics_loff), k = nrow(plot_metrics_loff), returnTrain = F)

predcond_loo_cv <- CAST::geodist(
  x = plot_metrics_loff,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = loo_cv
)

# nearest neighbour distance matching (nndm) 
# leave-one-out (loo) cv
nndm <- CAST::nndm(
  tpoints = plot_metrics_lon,
  modeldomain = sf::st_transform(sf::st_zm(ext_loff), sf::st_crs(plot_metrics_lon)),
  samplesize = 1000
)

predcond_cv_nndm <- CAST::geodist(
  x = plot_metrics_lon,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = nndm$indx_test
  )

# leave-location-out (llo) cv
# test different k-folds
k_values <- c(3, 5, 7, 10)
indices_llo_cv_results <- list()

for(k in k_values) {
  spatial_folds <- CAST::CreateSpacetimeFolds(
    x = plot_metrics_loff,
    spacevar = 'abt',
    k = k,
    seed = 11
  )
  
  indices_llo_cv_results[[paste0('k', k)]] <- CAST::geodist(
    x = plot_metrics_loff,
    modeldomain = sf::st_zm(ext_loff),
    cvfolds = spatial_folds$indexOut
  )
}

sapply(indices_llo_cv_results, function(x) attr(x, 'W_CV'))

plot_metrics_lon$abt <- as.character(plot_metrics_lon$abt)

indices_llo_cv <- CAST::CreateSpacetimeFolds(
  x = plot_metrics_lon,
  spacevar = 'abt',
  k = 10,
  seed = 22
)

predcond_llo_cv <- CAST::geodist(
  x = plot_metrics_lon,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = indices_llo_cv$indexOut
)

# plot ECDF function
plot(predcond_loo_cv, stat = 'ecdf')
plot(predcond_cv_fold5, stat = 'ecdf')
plot(nndm, type = 'simple')
plot(predcond_llo_cv, stat = 'ecdf')

# plot density function
plot(predcond_cv_fold5) + scale_x_log10(labels = round)
plot(predcond_cv_nndm) + scale_x_log10(labels = round)
plot(llo_cv) + scale_x_log10(labels = round)

# fit models and estimate their performance
# loo cv
vol_ha_lon_loo_ctrl <- caret::trainControl(method = 'LOOCV', savePredictions = T)
vol_ha_lon_loo_mod <- caret::train(
  plot_metrics_loff_df[,5:length(plot_metrics_loff_df)],
  plot_metrics_loff_df[,'vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_loo_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_loo_mod)

# 5-fold cv
vol_ha_lon_fold5_ctrl <- caret::trainControl(
  method = 'cv', number = 5, savePredictions = T
  )
vol_ha_lon_fold5_mod <- caret::train(
  plot_metrics_loff_df[,5:length(plot_metrics_loff_df)],
  plot_metrics_loff_df[,'vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_fold5_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_fold5_mod)

# nndm loo cv
vol_ha_lon_nndm_ctrl <- caret::trainControl(
  method = 'cv',
  index = nndm$indx_train,
  indexOut = nndm$indx_test,
  savePredictions = T
)
vol_ha_lon_nndm_mod <- caret::train(
  plot_metrics_lon_df[,6:length(plot_metrics_lon_df)],
  plot_metrics_lon_df[,'vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_nndm_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_nndm_mod)

# llo cv
vol_ha_lon_llo_ctrl <- caret::trainControl(
  method = 'cv',
  index = indices_llo_cv$index,
  savePredictions = T
)
vol_ha_lon_llo_mod <- caret::train(
  plot_metrics_lon_df[,6:length(plot_metrics_lon_df)],
  plot_metrics_lon_df[,'vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_llo_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_llo_mod)

# table with results
rbind(
  #data.frame(outcome="GSV", validation="LOO CV",
  #           t(as.data.frame(CAST::global_validation(vol_ha_lon_loo_mod)))),
  #data.frame(outcome="GSV", validation="5-fold CV",
  #           t(as.data.frame(CAST::global_validation(vol_ha_lon_fold5_mod)))),
  data.frame(outcome="GSV", validation="NNDM LOO CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_nndm_mod)))),
  data.frame(outcome="GSV", validation="LLO CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_llo_mod))))
) |> 
  knitr::kable(digits=2, row.names = F)

#                           END DRAFT - cv testing                             #
################################################################################
################################################################################
################################################################################
#                         START DRAFT - model training                         #

# split data into training and testing
set.seed(11)
trainIndex_lon <- caret::createDataPartition(
  plot_metrics_lon$vol_ha,
  p = 0.8, list = F)

train_lon <- plot_metrics_lon[trainIndex_lon,]
test_lon <- plot_metrics_lon[-trainIndex_lon,]

set.seed(11)
trainIndex_loff <- caret::createDataPartition(
  plot_metrics_loff$vol_ha,
  p = 0.8, list = F)

train_loff <- plot_metrics_loff[trainIndex_loff,]
test_loff <- plot_metrics_loff[-trainIndex_loff,]

names(train_lon) <- sub('^lon_', '', names(train_lon))
names(test_lon) <- sub('^lon_', '', names(test_lon))
names(train_loff) <- sub('^loff_', '', names(train_loff))
names(test_loff) <- sub('^loff_', '', names(test_loff))

saveRDS(
  train_lon,
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafon.RDS')
)
saveRDS(
  test_lon, 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon.RDS')
)
saveRDS(
  train_loff, 
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafoff.RDS')
)
saveRDS(
  test_loff,
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff.RDS')
)

# define predictors and response
predictors_lon <- sf::st_drop_geometry(train_lon[,5:length(train_lon)])
predictors_loff <- sf::st_drop_geometry(train_loff[,5:length(train_loff)])
response <- sf::st_drop_geometry(train_lon[,'vol_ha'])

# when including tree species
# convert to factor
predictors_lon$ts <- as.factor(predictors_lon$ts)
predictors_loff$ts <- as.factor(predictors_loff$ts)
levels(predictors_lon$ts)
levels(predictors_loff$ts)

# option with grouping into deciduous and coniferous
# 0,1,2,3 --> coniferous (1)
# 5,6,7,8 --> deciduous (2)
# 9 --> other (3)
# 666 --> canopy cover loss (0)
train_lon <- train_lon %>%
  dplyr::mutate(ts_gr = dplyr::case_when(
    ts %in% c(0, 1, 2, 3) ~ 1,
    ts %in% c(5, 6, 7, 8) ~ 2,
    ts == 9 ~ 3,
    ts == 666 ~ 0,
  ))
train_loff <- train_loff %>%
  dplyr::mutate(ts_gr = dplyr::case_when(
    ts %in% c(0, 1, 2, 3) ~ 1,
    ts %in% c(5, 6, 7, 8) ~ 2,
    ts == 9 ~ 3,
    ts == 666 ~ 0,
  ))

predictors_lon <- sf::st_drop_geometry(train_lon[,6:length(train_lon)])
predictors_loff <- sf::st_drop_geometry(train_loff[,6:length(train_loff)])
response <- sf::st_drop_geometry(train_lon[,'vol_ha'])

predictors_lon$ts_gr <- as.factor(predictors_lon$ts_gr)
predictors_loff$ts_gr <- as.factor(predictors_loff$ts_gr)
levels(predictors_lon$ts_gr)
levels(predictors_loff$ts_gr)

# initialize nearest neighbour distance matching
# leave-one-out cross-validation (NNDM LOO CV)
nndm <- CAST::nndm(
  tpoints = train_lon,
  modeldomain = sf::st_transform(sf::st_zm(ext_loff), sf::st_crs(train_lon)),
  samplesize = 1000
)

# control parameters for the train function
ctrl <- caret::trainControl(
  method = 'cv',
  index = nndm$indx_train,
  indexOut = nndm$indx_test,
  savePredictions = T,
  allowParallel = T
)

# initialize leave-location-out cross-validation (LLO CV)
train_lon$abt <- as.character(train_lon$abt)

indices <- CAST::CreateSpacetimeFolds(
  train_lon,
  spacevar = 'abt',
  k = 10,
  seed = 9999 
)

# control parameters for the train function
ctrl <- caret::trainControl(
  method = 'cv',
  index = indices$index,
  savePredictions = T,
  allowParallel = T
)

# create grid for tuning features
tgrid <- expand.grid(
  mtry = 1:length(predictors_lon),
  splitrule = c('variance', 'extratrees', 'maxstat'),
  min.node.size = c(5,10,15,20)
)

# create parallel cluster to increase computing speed
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# train random forest model
set.seed(11)
rf_model_lon <- caret::train(
  predictors_lon,
  response$vol_ha,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation'
)

set.seed(11)
rf_model_loff <- caret::train(
  predictors_loff,
  response$vol_ha,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation'
)

# stop parallel cluster
parallel::stopCluster(cl)

# get model performance and information
CAST::global_validation(rf_model_lon)
CAST::global_validation(rf_model_loff)
print(rf_model_lon)
print(rf_model_loff)

# correlation plots for lon and loff predictors
cor_lon <- stats::cor(predictors_lon, method = 'pearson')
cor_loff <- stats::cor(predictors_loff, method = 'pearson')
corrplot::corrplot(cor_lon, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)
corrplot::corrplot(cor_loff, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)

# variable importance
plot(caret::varImp(rf_model_lon))
plot(caret::varImp(rf_model_loff))

# train random forest model
# implementing forward feature selection
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_lon <- CAST::ffs(
  predictors_lon,
  response$vol_ha,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_loff <- CAST::ffs(
  predictors_loff,
  response$vol_ha,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

# save trained models
saveRDS(
  ffs_rf_model_lon,
  file.path(processed_data_dir, 'models' , 'ffs_rf_model_leafon.RDS')
  )
saveRDS(
  ffs_rf_model_loff,
  file.path(processed_data_dir, 'models', 'ffs_rf_model_leafoff.RDS')
  )

# get model performance and information
ffs_rf_model_lon
CAST::global_validation(ffs_rf_model_lon)
ffs_rf_model_loff
CAST::global_validation(ffs_rf_model_loff)

# plot model performance
plot(ffs_rf_model_lon)
plot(ffs_rf_model_loff)

# extract cross-validation predictions from the trained models
cv_pred_lon <- ffs_rf_model_lon$pred
cv_pred_loff <- ffs_rf_model_loff$pred

# the cv predictions contain multiple rows per observation due to resampling
# get final predictions for each observation
final_pred_lon <- cv_pred_lon %>%
  dplyr::filter(
    mtry == ffs_rf_model_lon$bestTune$mtry,
    splitrule == ffs_rf_model_lon$bestTune$splitrule,
    min.node.size == ffs_rf_model_lon$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_loff <- cv_pred_loff %>%
  dplyr::filter(
    mtry == ffs_rf_model_loff$bestTune$mtry,
    splitrule == ffs_rf_model_loff$bestTune$splitrule,
    min.node.size == ffs_rf_model_loff$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

# link it to the original geometries (BI plots used for training)
final_pred_lon_sf <- train_lon[final_pred_lon$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_lon$pred,
    obs = final_pred_lon$obs,
  )
sf::st_write(
  final_pred_lon_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_lon.gpkg')
)

final_pred_loff_sf <- train_loff[final_pred_loff$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_loff$pred,
    obs = final_pred_loff$obs
  )
sf::st_write(
  final_pred_loff_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_loff.gpkg')
)

# plot predicted vs. observed GSV with CV predictions
ggplot(final_pred_lon, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-on')

ggplot(final_pred_loff, aes(x=obs, y=pred)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, linewidth=1, color='red') +
  ggtitle('leaf-off')

#                           END DRAFT - model training                         #
################################################################################
























