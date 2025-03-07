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

# read data frame with inventory plots 
# and calculated metrics (leaf-on and leaf-off)
# --> see script forest_metrics.R
plot_metrics_leafon <- readRDS(
  file.path(processed_data_dir, 'plot_metrics_leafon.RDS')
  )

plot_metrics_leafoff <- readRDS(
  file.path(processed_data_dir, 'plot_metrics_leafoff.RDS')
)

head(plot_metrics_leafon)
str(plot_metrics_leafon)
head(plot_metrics_leafoff)
str(plot_metrics_leafoff)

# remove NA (empty plots)
plot_metrics_leafon <- na.omit(plot_metrics_leafon)
plot_metrics_leafoff <- na.omit(plot_metrics_leafoff)



# 03 - model preparation
#-------------------------------------

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
























