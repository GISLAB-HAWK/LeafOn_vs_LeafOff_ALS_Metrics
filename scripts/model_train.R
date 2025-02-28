#-------------------------------------------------------------------------------
# Name:         model_train.R
# Description:  Script trains a random forest to model the growing stock (m³/ha).
#               ALS-based metrics previously derived in forest inventory plots
#               are used as predictors.
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
# and calculated metrics
# --> see script forest_metrics.R
plot_metrics_aoi_leafon <- readRDS(
  file.path(processed_data_dir, 'plot_metrics_aoi_leafon.RDS')
  )

head(plot_metrics_aoi_leafon)
str(plot_metrics_aoi_leafon)

# remove NA
plot_metrics_aoi_leafon <- na.omit(plot_metrics_aoi_leafon)



# 03 - model preparation
#-------------------------------------

# split data into training and testing
plot_metrics_aoi_leafon_df <- as.data.frame(plot_metrics_aoi_leafon)
rownames(plot_metrics_aoi_leafon_df) <- 1:nrow(plot_metrics_aoi_leafon_df)
plot_metrics_aoi_leafon_df <- plot_metrics_aoi_leafon_df[1:length(plot_metrics_aoi_leafon_df)-1]

set.seed(11)

trainIndex <- caret::createDataPartition(
  plot_metrics_aoi_leafon_df$vol_ha,
  p = 0.8, list = F)

train <- plot_metrics_aoi_leafon_df[trainIndex,]
test <- plot_metrics_aoi_leafon_df[-trainIndex,]

saveRDS(train, file.path(processed_data_dir, 'train_ds_leafon.RDS'))
saveRDS(test, file.path(processed_data_dir, 'test_ds_leafon.RDS'))

# define predictors and response
predictors <- train[,7:length(train)]
response <- train[,'vol_ha']

# initialize leave-location-out cross-validation (LLO CV)
# this requires a spatial units variable as character
# in this case, the 'kspnr' variable is used for this
train$kspnr <- as.character(train$kspnr)

indices <- CAST::CreateSpacetimeFolds(
  train,
  spacevar = 'kspnr',
  k = 10
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
  mtry = 1:length(predictors),
  splitrule = 'variance',
  min.node.size = c(10,20,30,40,50)
 )

#tgrid <- expand.grid(
#  mtry = 1:10,
#  splitrule = 'variance',
#  min.node.size = c(10,20,30,40,50)
#)



# 04 - model training
#-------------------------------------

# create parallel cluster to increase computing speed
n_cores <- parallel::detectCores() - 2 
cl <- parallel::makeCluster(n_cores)
doParallel::registerDoParallel(cl)

# train random forest model
rf_model <- caret::train(
  predictors,
  response,
  method = 'ranger',
  trControl = ctrl,
  tuneGrid = tgrid,
  num.trees = 500,
  importance = 'permutation'
)

# stop parallel cluster
parallel::stopCluster(cl)

# save trained model
saveRDS(rf_model, file.path(output_dir, 'rf_model.RDS'))

print(rf_model)
summary(rf_model)

# plot predicted vs. observed growing stock
library(ggplot2)
ggplot(train, aes(x=vol_ha, y=stats::predict(rf_model))) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,1200), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red')
























