#-----------------------------------------------------------------------------
# Name:         val.R
# Description:  Script validates the predicted growing stock (m³/ha).
#               Validation is done on plot-level. Independent test plots
#               not used during model training are used.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-----------------------------------------------------------------------------



# 01 - file path definitions
#----------------------------

# define processed data directory
processed_data_dir <- 'data/processed/'

# define output directory
output_dir <- 'output/'



# 02 - data reading
#-------------------------------------

# trained random forest model
rf_model <- readRDS(file.path(output_dir, 'rf_model.RDS'))

# test dataset leaf-on (20% of the plots)
test_ds_leafon <- readRDS(file.path(processed_data_dir, 'test_ds_leafon.RDS'))



# 03 - validation: plot-level
#-------------------------------------

# predict on unseen test dataset
pred_test_leafon <- stats::predict(rf_model, test_ds_leafon)

# create data frame with predicted and observed growing stock
val_df <- data.frame(
  kspnr = dplyr::pull(test_ds_leafon, 'kspnr'),
  observed = dplyr::pull(test_ds_leafon, 'vol_ha'),
  predicted = pred_test_leafon)

# calculate error metrics (RMSE, rel. RMSE, bias, R-squared)
rmse <- round(sqrt(
  mean((val_df$predicted - val_df$observed)^2, na.rm = T)), 2
  )

rrmse <- round((rmse / mean(val_df$observed)) * 100, 2)

bias <- mean(val_df$predicted - val_df$observed)

ssr <- sum((val_df$predicted - mean(val_df$observed))^2)
sst <- sum((val_df$observed - mean(val_df$observed))^2)
rsquared <- 1 - (ssr / sst)

cat('RMSE:', rmse, '\n')
cat('rRMSE:', rrmse, '\n')
cat('Bias:', bias, '\n')
cat('R-squared:', rsquared, '\n')

# predicted vs. observed grwoing stock
library(ggplot2)
ggplot(val_df, aes(x=observed, y=predicted)) +
  geom_point() +
  xlab(expression(paste('observed growing stock [', m^3, ha^-1, ']', sep = ''))) +
  ylab(expression(paste('predicted growing stock [', m^3, ha^-1, ']', sep = ''))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  scale_y_continuous(limits=c(0,700), breaks=seq(0,1500, by=300)) +
  geom_abline(slope=1, intercept=0, size=1, color='red')















