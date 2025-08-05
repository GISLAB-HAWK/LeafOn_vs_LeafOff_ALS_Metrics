#-------------------------------------------------------------------------------
# Name:         model_train.R
# Description:  Script trains a random forest to model the growing stock (m³/ha).
#               ALS-based metrics previously derived in forest inventory plots
#               are used as predictors. Two models are trained, one using the
#               the metrics calculated from leaf-on dataset and one using
#               the metrics calculated from leaf-off dataset.
#               Leave-Location-Out cross-validation (LLO CV) is used as a 
#               spatial cross validation method.
#               80% of the data is used for training, 20% for testing.
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
  file.path(processed_data_dir, 'metrics', 'vol_stp_species_metrics.gpkg')
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

# HERE OR LATER?
# remove NA (empty plots)
plot_metrics <- na.omit(plot_metrics)

# split data for leaf-on and leaf-off
id_cols <- c('key', 'kspnr', 'vol_ha', 'ts')

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


# HERE OR BEFORE?
# remove NA (empty plots)
plot_metrics_lon <- na.omit(plot_metrics_lon)
plot_metrics_loff <- na.omit(plot_metrics_loff)



# 03 - model preparation
#-------------------------------------------------------------------------------

################################################################################
#                               NEW - DRAFT                                    #

# split data into training and testing
set.seed(11)

trainIndex_lon <- caret::createDataPartition(
  plot_metrics_lon$vol_ha,
  p = 0.8, list = F)

train_lon <- plot_metrics_lon[trainIndex_lon,]
test_lon <- plot_metrics_lon[-trainIndex_lon,]

trainIndex_loff <- caret::createDataPartition(
  plot_metrics_loff$vol_ha,
  p = 0.8, list = F)

train_loff <- plot_metrics_loff[trainIndex_loff,]
test_loff <- plot_metrics_loff[-trainIndex_loff,]

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

# df versions
train_lon_df <- as.data.frame(sf::st_drop_geometry(train_lon))
test_lon_df <- as.data.frame(sf::st_drop_geometry(test_lon))
train_loff_df <- as.data.frame(sf::st_drop_geometry(train_loff))
test_loff_df <- as.data.frame(sf::st_drop_geometry(test_loff))

# random 5-fold cross-validation (cv)
fold5 <- caret::createFolds(1:nrow(train_lon), k = 5, returnTrain = F)

# explore geographic predictive conditions
predcond_cv_fold5 <- CAST::geodist(
  x = train_lon,
  modeldomain = sf::st_zm(ext_loff),
  cvfolds = fold5
  )

predcond_cv_fold5 <- CAST::geodist(
  x = train_lon,
  modeldomain = sf::st_zm(ext_loff),
  cvtrain = trainIndex_lon,
  testdata = test_lon
)

# nearest neighbour distance matching (nndm)
nndm <- CAST::nndm(
  train_lon,
  modeldomain = sf::st_transform(sf::st_zm(ext_loff), sf::st_crs(train_lon)),
  samplesize = 1000
)

# plot ECDF function
plot(predcond_cv_fold5, stat = 'ecdf')
plot(nndm, type = 'simple')




################################################################################








# split data into training and testing
plot_metrics_leafon_df <- as.data.frame(plot_metrics_leafon)
rownames(plot_metrics_leafon_df) <- 1:nrow(plot_metrics_leafon_df)
plot_metrics_leafon_df <- plot_metrics_leafon_df[1:length(plot_metrics_leafon_df)-1]

plot_metrics_leafoff_df <- as.data.frame(plot_metrics_leafoff)
rownames(plot_metrics_leafoff_df) <- 1:nrow(plot_metrics_leafoff_df)
plot_metrics_leafoff_df <- plot_metrics_leafoff_df[1:length(plot_metrics_leafoff_df)-1]

set.seed(11)

trainIndex_leafon <- caret::createDataPartition(
  plot_metrics_leafon_df$vol_ha,
  p = 0.8, list = F)

train_leafon <- plot_metrics_leafon_df[trainIndex_leafon,]
test_leafon <- plot_metrics_leafon_df[-trainIndex_leafon,]

trainIndex_leafoff <- caret::createDataPartition(
  plot_metrics_leafoff_df$vol_ha,
  p = 0.8, list = F)

train_leafoff <- plot_metrics_leafoff_df[trainIndex_leafoff,]
test_leafoff <- plot_metrics_leafoff_df[-trainIndex_leafoff,]

saveRDS(train_leafon, file.path(processed_data_dir, 'train_ds_leafon.RDS'))
saveRDS(test_leafon, file.path(processed_data_dir, 'test_ds_leafon.RDS'))
saveRDS(train_leafoff, file.path(processed_data_dir, 'train_ds_leafoff.RDS'))
saveRDS(test_leafoff, file.path(processed_data_dir, 'test_ds_leafoff.RDS'))

# define predictors and response
predictors_leafon <- train_leafon[,7:length(train_leafon)]
response_leafon <- train_leafon[,'vol_ha']
predictors_leafoff <- train_leafoff[,7:length(train_leafoff)]
response_leafoff <- train_leafoff[,'vol_ha']

# initialize leave-location-out cross-validation (LLO CV)
# this requires a spatial units variable as character
# in this case, the 'kspnr' variable is used for this
train_leafon$kspnr <- as.character(train_leafon$kspnr)
train_leafoff$kspnr <- as.character(train_leafoff$kspnr)

indices_leafon <- CAST::CreateSpacetimeFolds(
  train_leafon,
  spacevar = 'kspnr',
  k = 10
  )

indices_leafoff <- CAST::CreateSpacetimeFolds(
  train_leafoff,
  spacevar = 'kspnr',
  k = 10
)

# control parameters for the train function
ctrl_leafon <- caret::trainControl(
  method = 'cv',
  index = indices_leafon$index,
  savePredictions = T,
  allowParallel = T
  )

ctrl_leafoff <- caret::trainControl(
  method = 'cv',
  index = indices_leafoff$index,
  savePredictions = T,
  allowParallel = T
)

# create grid for tuning features
tgrid_leafon <- expand.grid(
  mtry = 1:length(predictors_leafon),
  splitrule = 'variance',
  min.node.size = c(10,20,30,40,50)
 )

tgrid_leafoff <- expand.grid(
  mtry = 1:length(predictors_leafoff),
  splitrule = 'variance',
  min.node.size = c(10,20,30,40,50)
)



# 04 - model training
#-------------------------------------

# create parallel cluster to increase computing speed
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# train random forest model
rf_model_leafon <- caret::train(
  predictors_leafon,
  response_leafon,
  method = 'ranger',
  trControl = ctrl_leafon,
  tuneGrid = tgrid_leafon,
  num.trees = 500,
  importance = 'permutation'
)

rf_model_leafoff <- caret::train(
  predictors_leafoff,
  response_leafoff,
  method = 'ranger',
  trControl = ctrl_leafoff,
  tuneGrid = tgrid_leafoff,
  num.trees = 500,
  importance = 'permutation'
)

# stop parallel cluster
parallel::stopCluster(cl)

# save trained model
saveRDS(rf_model_leafon, file.path(output_dir, 'rf_model_leafon.RDS'))
saveRDS(rf_model_leafoff, file.path(output_dir, 'rf_model_leafoff.RDS'))

print(rf_model_leafon)
summary(rf_model_leafon)
print(rf_model_leafoff)
summary(rf_model_leafoff)

# plot predicted vs. observed growing stock
library(ggplot2)
library(ggpubr)

plot_leafon <- ggplot(train_leafon, aes(x=vol_ha, y=stats::predict(rf_model_leafon))) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red') +
  ggtitle('leaf-on')

plot_leafoff <- ggplot(train_leafoff, aes(x=vol_ha, y=stats::predict(rf_model_leafoff))) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red') +
  ggtitle('leaf-off')

ggpubr::ggarrange(plot_leafon, plot_leafoff, ncol = 2, nrow = 1)
























