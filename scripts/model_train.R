#-------------------------------------------------------------------------------
# Name:         model_train.R
# Description:  Script trains random forest models for predicting
#               growing stock volume (GSV) (m³/ha).
#               ALS-based metrics previously derived in forest inventory plots
#               (BI plots) are used as predictors. Two models are trained,
#               one using the the metrics calculated from leaf-on dataset,
#               and one using the metrics calculated from leaf-off dataset.
#               Nearest Neighbour Distance Matching (NNDM) 
#               Leave-One-Out Cross Validation (LOO CV) is used as a 
#               spatial map validation method.
#               The models are trained separately for deciduous and coniferous
#               dominated plots.
#               80% of each dataset is used for training the models.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read data with forest inventory plots (BI) 
# and calculated metrics (leaf-on and leaf-off)
plot_metrics_lon <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon.gpkg')
)

plot_metrics_loff <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff.gpkg')
)

head(plot_metrics_lon)
head(plot_metrics_lon)
str(plot_metrics_lon)
str(plot_metrics_loff)

# read extent of the leaf-off dataset
ext_loff <- sf::st_read(
  file.path(processed_data_dir, 'pc_leafoff_2024', 'leafoff_ppm4.vpc')
  )

# plot data
#ggplot() +
#  geom_sf(data = ext_loff, fill = "grey", alpha = 0.1) +
#  geom_sf(data = plot_metrics_loff, aes(col = total_vol_ha)) +
#  scale_colour_distiller(palette = "YlGn", direction = 1) +
#  theme_bw() +
#  labs(col = "") +
#  ggtitle("Growing Stock Volume (m³/ha)")



# 02: data preparation
#-------------------------------------------------------------------------------

# filter plots by dominant leaf type
plots_lon_deciduous <- plot_metrics_lon %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')

plots_lon_coniferous <- plot_metrics_lon %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

plots_loff_deciduous <- plot_metrics_loff %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')

plots_loff_coniferous <- plot_metrics_loff %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# df versions
plots_lon_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_deciduous))
plots_lon_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_coniferous))
plots_loff_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_deciduous))
plots_loff_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_coniferous))



# 03 - model preparation
#-------------------------------------------------------------------------------

################################################################################
#                           START DRAFT - cv testing                           #

# leave-one-out cv (loocv)
loo_cv <- caret::createFolds(
  1:nrow(plots_lon_coniferous), k = nrow(plots_lon_coniferous), returnTrain = F
)

predcond_loo_cv <- CAST::geodist(
  x = plots_lon_coniferous,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = loo_cv
)

# random 5-fold cross-validation (cv)
fold5 <- caret::createFolds(1:nrow(plots_lon_coniferous), k = 5, returnTrain = F)

predcond_cv_fold5 <- CAST::geodist(
  x = plots_lon_coniferous,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = fold5
  )

# nearest neighbour distance matching (nndm) 
# leave-one-out (loo) cv
nndm <- CAST::nndm(
  tpoints = plots_lon_coniferous,
  #modeldomain = sf::st_transform(
  #  sf::st_zm(ext_loff),
  #  sf::st_crs(plots_lon_coniferous)),
  predpoints = pixel_centroids,
  #samplesize = 250
)

# k-nearest neighbour distance matching (knndm)
knndm <- CAST::knndm(
  tpoints = plots_lon_coniferous,
  predpoints = pixel_centroids,
  k = 5,
  clustering = 'kmeans',
)

# plot ECDF function
plot(predcond_loo_cv, stat = 'ecdf')
plot(predcond_cv_fold5, stat = 'ecdf')
plot(nndm, type = 'simple')
plot(knndm, type = 'simple')

# plot showing how nndm works
# cv iteration with the most excluded plots
id_plot <- which.max(sapply(nndm$indx_exclude, length))
lon_coniferous_plot <- plots_lon_coniferous
lon_coniferous_plot$set <- ""
lon_coniferous_plot$set[nndm$indx_train[[id_plot]]] <- 'train'
lon_coniferous_plot$set[nndm$indx_exclude[[id_plot]]] <- 'exclude'
lon_coniferous_plot$set[nndm$indx_test[[id_plot]]] <- 'test'
lon_coniferous_plot <- lon_coniferous_plot[order(lon_coniferous_plot$set),]

ggplot() +
  geom_sf(data = ext_loff, fill = 'grey', alpha = 0.1) +
  geom_sf(data = lon_coniferous_plot, aes(col = set)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw()

# plot showing how knndm works
ggplot() +
  geom_sf(data = ext_loff, fill = "grey", alpha = 0.1) +
  geom_sf(data = plots_lon_coniferous, aes(col = as.factor(knndm$clusters))) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() + theme(legend.position = 'none')

# highlight the clusters in knndm
ggplot() +
  geom_sf(data = ext_loff, fill = "grey", alpha = 0.1) +
  geom_sf(data = plots_lon_coniferous, 
          aes(col = as.factor(knndm$clusters)), alpha = 0.3) +
  geom_sf(data = plots_lon_coniferous[knndm$clusters == 5, ], 
          col = "black", size = 2) +
  theme_bw() + theme(legend.position = 'none')

# fit models and estimate their performance
# loo cv
vol_ha_lon_loo_ctrl <- caret::trainControl(method = 'LOOCV', savePredictions = T)
vol_ha_lon_loo_mod <- caret::train(
  plots_lon_coniferous_df[,12:length(plots_lon_coniferous_df)],
  plots_lon_coniferous_df[,'total_vol_ha'],
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
  plots_lon_coniferous_df[,12:length(plots_lon_coniferous_df)],
  plots_lon_coniferous_df[,'total_vol_ha'],
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
  plots_lon_coniferous_df[,12:length(plots_lon_coniferous_df)],
  plots_lon_coniferous_df[,'total_vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_nndm_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_nndm_mod)

# knndm cv
vol_ha_lon_knndm_ctrl <- caret::trainControl(
  method = 'cv',
  index = knndm$indx_train,
  indexOut = knndm$indx_test,
  savePredictions = T
)
vol_ha_lon_knndm_mod <- caret::train(
  plots_lon_coniferous_df[,12:length(plots_lon_coniferous_df)],
  plots_lon_coniferous_df[,'total_vol_ha'],
  method = 'rf',
  importance = F,
  trControl = vol_ha_lon_knndm_ctrl,
  ntree = 100,
  tuneLength = 1
)
CAST::global_validation(vol_ha_lon_knndm_mod)

# table with results
rbind(
  data.frame(outcome="GSV", validation="LOO CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_loo_mod)))),
  data.frame(outcome="GSV", validation="5-fold CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_fold5_mod)))),
  data.frame(outcome="GSV", validation="NNDM LOO CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_nndm_mod)))),
  data.frame(outcome="GSV", validation="kNNDM LOO CV",
             t(as.data.frame(CAST::global_validation(vol_ha_lon_knndm_mod))))
) |> 
  knitr::kable(digits=2, row.names = F)



#                           END DRAFT - cv testing                             #
################################################################################
################################################################################
################################################################################
#                         START DRAFT - model training                         #

# split data into training and testing

# deciduous
set.seed(11)
trainIndex_lon_deciduous <- caret::createDataPartition(
  plots_lon_deciduous$vol_ha,
  p = 0.8, list = F)

train_lon_deciduous <- plots_lon_deciduous[trainIndex_lon_deciduous,]
test_lon_deciduous <- plots_lon_deciduous[-trainIndex_lon_deciduous,]

set.seed(11)
trainIndex_loff_deciduous <- caret::createDataPartition(
  plots_loff_deciduous$vol_ha,
  p = 0.8, list = F)

train_loff_deciduous <- plots_loff_deciduous[trainIndex_loff_deciduous,]
test_loff_deciduous <- plots_loff_deciduous[-trainIndex_loff_deciduous,]

names(train_lon_deciduous) <- sub('^lon_', '', names(train_lon_deciduous))
names(test_lon_deciduous) <- sub('^lon_', '', names(test_lon_deciduous))
names(train_loff_deciduous) <- sub('^loff_', '', names(train_loff_deciduous))
names(test_loff_deciduous) <- sub('^loff_', '', names(test_loff_deciduous))

saveRDS(
  train_lon_deciduous,
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafon_deciduous_filtered.RDS')
)
saveRDS(
  test_lon_deciduous, 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon_deciduous_filtered.RDS')
)
saveRDS(
  train_loff_deciduous, 
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafoff_deciduous_filtered.RDS')
)
saveRDS(
  test_loff_deciduous,
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff_deciduous_filtered.RDS')
)

# coniferous
set.seed(11)
trainIndex_lon_coniferous <- caret::createDataPartition(
  plots_lon_coniferous$vol_ha,
  p = 0.8, list = F)

train_lon_coniferous <- plots_lon_coniferous[trainIndex_lon_coniferous,]
test_lon_coniferous <- plots_lon_coniferous[-trainIndex_lon_coniferous,]

set.seed(11)
trainIndex_loff_coniferous <- caret::createDataPartition(
  plots_loff_coniferous$vol_ha,
  p = 0.8, list = F)

train_loff_coniferous <- plots_loff_coniferous[trainIndex_loff_coniferous,]
test_loff_coniferous <- plots_loff_coniferous[-trainIndex_loff_coniferous,]

names(train_lon_coniferous) <- sub('^lon_', '', names(train_lon_coniferous))
names(test_lon_coniferous) <- sub('^lon_', '', names(test_lon_coniferous))
names(train_loff_coniferous) <- sub('^loff_', '', names(train_loff_coniferous))
names(test_loff_coniferous) <- sub('^loff_', '', names(test_loff_coniferous))

saveRDS(
  train_lon_coniferous,
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafon_coniferous_filtered.RDS')
)
saveRDS(
  test_lon_coniferous, 
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafon_coniferous_filtered.RDS')
)
saveRDS(
  train_loff_coniferous, 
  file.path(processed_data_dir, 'train_test_ds', 'train_ds_leafoff_coniferous_filtered.RDS')
)
saveRDS(
  test_loff_coniferous,
  file.path(processed_data_dir, 'train_test_ds', 'test_ds_leafoff_coniferous_filtered.RDS')
)

# define predictors and response
predictors_lon_deciduous <- sf::st_drop_geometry(
  train_lon_deciduous[,10:length(train_lon_deciduous)]
  )
predictors_loff_deciduous <- sf::st_drop_geometry(
  train_loff_deciduous[,10:length(train_loff_deciduous)]
  )
predictors_lon_coniferous <- sf::st_drop_geometry(
  train_lon_coniferous[,10:length(train_lon_coniferous)]
)
predictors_loff_coniferous <- sf::st_drop_geometry(
  train_loff_coniferous[,10:length(train_loff_coniferous)]
)
response_deciduous <- sf::st_drop_geometry(train_lon_deciduous[,'vol_ha'])
response_coniferous <- sf::st_drop_geometry(train_lon_coniferous[,'vol_ha'])

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
nndm_deciduous <- CAST::nndm(
  tpoints = train_lon_deciduous,
  modeldomain = sf::st_transform(sf::st_zm(ext_loff), sf::st_crs(train_lon_deciduous)),
  samplesize = 500
)
nndm_coniferous <- CAST::nndm(
  tpoints = train_lon_coniferous,
  modeldomain = sf::st_transform(sf::st_zm(ext_loff), sf::st_crs(train_lon_coniferous)),
  samplesize = 500
)

# control parameters for the train function
ctrl_deciduous <- caret::trainControl(
  method = 'cv',
  index = nndm_deciduous$indx_train,
  indexOut = nndm_deciduous$indx_test,
  savePredictions = T,
  allowParallel = T
)
ctrl_coniferous <- caret::trainControl(
  method = 'cv',
  index = nndm_coniferous$indx_train,
  indexOut = nndm_coniferous$indx_test,
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
  mtry = 1:length(predictors_lon_deciduous),
  splitrule = c('variance', 'extratrees', 'maxstat'),
  min.node.size = c(5,10,15,20)
)

# create parallel cluster to increase computing speed
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# train random forest model
set.seed(11)
rf_model_lon_deciduous <- caret::train(
  predictors_lon_deciduous,
  response$vol_ha,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation'
)

set.seed(11)
rf_model_loff_deciduous <- caret::train(
  predictors_loff_deciduous,
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
CAST::global_validation(rf_model_lon_deciduous)
CAST::global_validation(rf_model_loff_deciduous)
print(rf_model_lon)
print(rf_model_loff)

# correlation plots for lon and loff predictors
cor_lon <- stats::cor(predictors_lon_deciduous, method = 'pearson')
cor_loff <- stats::cor(predictors_loff_deciduous, method = 'pearson')
corrplot::corrplot(cor_lon, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)
corrplot::corrplot(cor_loff, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)

# variable importance
plot(caret::varImp(rf_model_lon_deciduous))
plot(caret::varImp(rf_model_loff_deciduous))

# train random forest model
# implementing forward feature selection

# deciduous
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_lon_deciduous <- CAST::ffs(
  predictors_lon_deciduous,
  response_deciduous$vol_ha,
  method = 'ranger',
  trControl = ctrl_deciduous,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_loff_deciduous <- CAST::ffs(
  predictors_loff_deciduous,
  response_deciduous$vol_ha,
  method = 'ranger',
  trControl = ctrl_deciduous,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

# save trained models
saveRDS(
  ffs_rf_model_lon_deciduous,
  file.path(processed_data_dir, 'models' , 'ffs_rf_model_leafon_deciduous_filtered.RDS')
  )
saveRDS(
  ffs_rf_model_loff_deciduous,
  file.path(processed_data_dir, 'models', 'ffs_rf_model_leafoff_deciduous_filtered.RDS')
  )

# coniferous
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_lon_coniferous <- CAST::ffs(
  predictors_lon_coniferous,
  response_coniferous$vol_ha,
  method = 'ranger',
  trControl = ctrl_coniferous,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

ffs_rf_model_loff_coniferous <- CAST::ffs(
  predictors_loff_coniferous,
  response_coniferous$vol_ha,
  method = 'ranger',
  trControl = ctrl_coniferous,
  tuneGrid = tgrid,
  num.trees = 100,
  importance = 'permutation',
  seed = 999
)

parallel::stopCluster(cl)

# save trained models
saveRDS(
  ffs_rf_model_lon_coniferous,
  file.path(processed_data_dir, 'models' , 'ffs_rf_model_leafon_coniferous_filtered.RDS')
)
saveRDS(
  ffs_rf_model_loff_coniferous,
  file.path(processed_data_dir, 'models', 'ffs_rf_model_leafoff_coniferous_filtered.RDS')
)

# get model performance and information
ffs_rf_model_lon_deciduous
CAST::global_validation(ffs_rf_model_lon_deciduous)
ffs_rf_model_loff_deciduous
CAST::global_validation(ffs_rf_model_loff_deciduous)

ffs_rf_model_lon_coniferous
CAST::global_validation(ffs_rf_model_lon_coniferous)
ffs_rf_model_loff_coniferous
CAST::global_validation(ffs_rf_model_loff_coniferous)

# plot model performance
plot(ffs_rf_model_lon)
plot(ffs_rf_model_loff)

# extract cross-validation predictions from the trained models
cv_pred_lon_deciduous <- ffs_rf_model_lon_deciduous$pred
cv_pred_loff_deciduous <- ffs_rf_model_loff_deciduous$pred
cv_pred_lon_coniferous <- ffs_rf_model_lon_coniferous$pred
cv_pred_loff_coniferous <- ffs_rf_model_loff_coniferous$pred

# the cv predictions contain multiple rows per observation due to resampling
# get final predictions for each observation

# deciduous
final_pred_lon_deciduous <- cv_pred_lon_deciduous %>%
  dplyr::filter(
    mtry == ffs_rf_model_lon_deciduous$bestTune$mtry,
    splitrule == ffs_rf_model_lon_deciduous$bestTune$splitrule,
    min.node.size == ffs_rf_model_lon_deciduous$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_loff_deciduous <- cv_pred_loff_deciduous %>%
  dplyr::filter(
    mtry == ffs_rf_model_loff_deciduous$bestTune$mtry,
    splitrule == ffs_rf_model_loff_deciduous$bestTune$splitrule,
    min.node.size == ffs_rf_model_loff_deciduous$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

# link it to the original geometries (BI plots used for training)
final_pred_lon_deciduous_sf <- train_lon_deciduous[final_pred_lon_deciduous$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_lon_deciduous$pred,
    obs = final_pred_lon_deciduous$obs,
  )
sf::st_write(
  final_pred_lon_deciduous_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_lon_deciduous_filtered.gpkg')
)

final_pred_loff_deciduous_sf <- train_loff_deciduous[final_pred_loff_deciduous$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_loff_deciduous$pred,
    obs = final_pred_loff_deciduous$obs
  )
sf::st_write(
  final_pred_loff_deciduous_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_loff_deciduous_filtered.gpkg')
)

# coniferous
final_pred_lon_coniferous <- cv_pred_lon_coniferous %>%
  dplyr::filter(
    mtry == ffs_rf_model_lon_coniferous$bestTune$mtry,
    splitrule == ffs_rf_model_lon_coniferous$bestTune$splitrule,
    min.node.size == ffs_rf_model_lon_coniferous$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

final_pred_loff_coniferous <- cv_pred_loff_coniferous %>%
  dplyr::filter(
    mtry == ffs_rf_model_loff_coniferous$bestTune$mtry,
    splitrule == ffs_rf_model_loff_coniferous$bestTune$splitrule,
    min.node.size == ffs_rf_model_loff_coniferous$bestTune$min.node.size) %>%
  dplyr::group_by(rowIndex) %>%
  dplyr::summarise(pred = mean(pred), obs = first(obs), .groups = 'drop')

# link it to the original geometries (BI plots used for training)
final_pred_lon_coniferous_sf <- train_lon_coniferous[final_pred_lon_coniferous$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_lon_coniferous$pred,
    obs = final_pred_lon_coniferous$obs,
  )
sf::st_write(
  final_pred_lon_coniferous_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_lon_coniferous_filtered.gpkg')
)

final_pred_loff_coniferous_sf <- train_loff_coniferous[final_pred_loff_coniferous$rowIndex, ] %>%
  dplyr::select(key, kspnr, abt) %>%
  dplyr::mutate(
    pred = final_pred_loff_coniferous$pred,
    obs = final_pred_loff_coniferous$obs
  )
sf::st_write(
  final_pred_loff_coniferous_sf,
  file.path(processed_data_dir, 'predictions', 'pred_obsv_train_loff_coniferous_filtered.gpkg')
)

# plot predicted vs. observed GSV with CV predictions

# deciduous
ggplot(final_pred_lon_deciduous, aes(x=obs, y=pred)) +
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

ggplot(final_pred_loff_deciduous, aes(x=obs, y=pred)) +
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

final_pred_lon_loff_deciduous <- dplyr::bind_rows(
  final_pred_lon_deciduous %>% dplyr::mutate(condition = 'leaf-on'),
  final_pred_loff_deciduous %>% dplyr::mutate(condition = 'leaf-off')
)

stats_summary_deciduous <- final_pred_lon_loff_deciduous %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    rrmse = (sqrt(mean((pred - obs)^2, na.rm = T)) / mean(obs, na.rm = T)),
    r_squared = stats::cor(obs, pred, use = 'complete.obs')^2,
    .groups = 'drop'
  )

text_labels_deciduous <- paste0(
  'leaf-on: rRMSE = ', round(stats_summary_deciduous$rrmse[stats_summary_deciduous$condition == 'leaf-on'], 2),
  ', R² = ', round(stats_summary_deciduous$r_squared[stats_summary_deciduous$condition == 'leaf-on'], 2), '\n',
  'leaf-off: rRMSE = ', round(stats_summary_deciduous$rrmse[stats_summary_deciduous$condition == 'leaf-off'], 2),
  ', R² = ', round(stats_summary_deciduous$r_squared[stats_summary_deciduous$condition == 'leaf-off'], 2)
)

ggplot(final_pred_lon_loff_deciduous, aes(x = obs, y = pred, color = condition)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200)) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = c('leaf-on' = 'gray60', 'leaf-off' = 'gray30'),
                     name = 'Condition') +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, color = 'red') +
  annotate('text', x = 475, y = 10, label = text_labels_deciduous, 
           hjust = 0, vjust = 0, size = 3.5, color = 'black')

# coniferous
ggplot(final_pred_lon_coniferous, aes(x=obs, y=pred)) +
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

ggplot(final_pred_loff_coniferous, aes(x=obs, y=pred)) +
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

final_pred_lon_loff_coniferous <- dplyr::bind_rows(
  final_pred_lon_coniferous %>% dplyr::mutate(condition = 'leaf-on'),
  final_pred_loff_coniferous %>% dplyr::mutate(condition = 'leaf-off')
)

stats_summary_coniferous <- final_pred_lon_loff_coniferous %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(
    rmse = sqrt(mean((pred - obs)^2, na.rm = T)),
    rrmse = (sqrt(mean((pred - obs)^2, na.rm = T)) / mean(obs, na.rm = T)),
    r_squared = stats::cor(obs, pred, use = 'complete.obs')^2,
    .groups = 'drop'
  )

text_labels_coniferous <- paste0(
  'leaf-on: rRMSE = ', round(stats_summary_coniferous$rrmse[stats_summary_coniferous$condition == 'leaf-on'], 2),
  ', R² = ', round(stats_summary_coniferous$r_squared[stats_summary_coniferous$condition == 'leaf-on'], 2), '\n',
  'leaf-off: rRMSE = ', round(stats_summary_coniferous$rrmse[stats_summary_coniferous$condition == 'leaf-off'], 2),
  ', R² = ', round(stats_summary_coniferous$r_squared[stats_summary_coniferous$condition == 'leaf-off'], 2)
)

ggplot(final_pred_lon_loff_coniferous, aes(x = obs, y = pred, color = condition)) +
  geom_point() +
  xlab(expression(paste('observed GSV [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted GSV [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top") +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200)) +
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200)) +
  scale_color_manual(values = c('leaf-on' = 'gray60', 'leaf-off' = 'gray30'),
                     name = 'Condition') +
  geom_abline(slope = 1, intercept = 0, linewidth = 1, color = 'red') +
  annotate('text', x = 475, y = 10, label = text_labels_coniferous, 
           hjust = 0, vjust = 0, size = 3.5, color = 'black')

#                           END DRAFT - model training                         #
################################################################################
























