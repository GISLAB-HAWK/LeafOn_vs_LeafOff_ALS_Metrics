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

# read pixel centroids from w2w-metrics raster
pixel_centroids <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'pixel_centroids_dom_leaf_type.gpkg')
  )

# split pixel centroids based on dominant leaf type
pixel_centroids_deciduous <- pixel_centroids %>% 
  dplyr::filter(dominant_leaf_type1 == 1)
pixel_centroids_coniferous <- pixel_centroids %>% 
  dplyr::filter(dominant_leaf_type1 == 2)

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

# define predictors (ALS metrics)
predictors_lon_deciduous <- sf::st_drop_geometry(
  plots_lon_deciduous[,12:length(plots_lon_deciduous)]
  )
predictors_loff_deciduous <- sf::st_drop_geometry(
  plots_loff_deciduous[,12:length(plots_loff_deciduous)]
  )
predictors_lon_coniferous <- sf::st_drop_geometry(
  plots_lon_coniferous[,12:length(plots_lon_coniferous)]
)
predictors_loff_coniferous <- sf::st_drop_geometry(
  plots_loff_coniferous[,12:length(plots_loff_coniferous)]
)

# initialize nearest neighbour distance matching
# leave-one-out cross-validation (NNDM LOO CV)
nndm_deciduous <- CAST::nndm(
  tpoints = plots_lon_deciduous,
  predpoints = pixel_centroids_deciduous
)
nndm_coniferous <- CAST::nndm(
  tpoints = plots_lon_coniferous,
  predpoints = pixel_centroids_coniferous
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

# create grid for tuning features
tgrid <- expand.grid(
  mtry = 1:3,
  splitrule = c('variance', 'extratrees', 'maxstat'),
  min.node.size = c(5, 10, 15, 20)
)

# correlation plots for lon and loff predictors
cor_lon_deciduous <- stats::cor(predictors_lon_deciduous, method = 'pearson')
cor_loff_deciduous <- stats::cor(predictors_loff_deciduous, method = 'pearson')
cor_lon_coniferous <- stats::cor(predictors_lon_coniferous, method = 'pearson')
cor_loff_coniferous <- stats::cor(predictors_loff_coniferous, method = 'pearson')
corrplot::corrplot(cor_lon_deciduous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)
corrplot::corrplot(cor_loff_deciduous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)
corrplot::corrplot(cor_lon_coniferous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)
corrplot::corrplot(cor_loff_coniferous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA)

# train random forest models
# implementing forward feature selection
# for all response variables, leaf types, and leaf conditions

# define response variables with their names (for file naming)
response_vars <- list(
  list(name = 'total_vol_ha', col = 'total_vol_ha'),
  list(name = 'merch_vol_ha', col = 'merch_vol_ha'),
  list(name = 'agb_ha', col = 'agb_ha'),
  list(name = 'tree_density', col = 'tree_density'),
  list(name = 'basal_area_ha', col = 'basal_area_ha'),
  list(name = 'dg', col = 'dg')
)

# define leaf types and conditions
leaf_types <- c('deciduous', 'coniferous')
leaf_conditions <- c('lon', 'loff')

# create list to store all models
ffs_models <- list()

# check which models already exist
models_to_train <- list()
for (resp in response_vars) {
  for (leaf_type in leaf_types) {
    for (leaf_cond in leaf_conditions) {
      
      model_name <- paste0('ffs_rf_', resp$name, '_', leaf_cond, '_', leaf_type)
      file_name <- paste0('ffs_rf_', resp$name, '_leaf', 
                          ifelse(leaf_cond == 'lon', 'on', 'off'), 
                          '_', leaf_type, '.RDS')
      file_path <- file.path(processed_data_dir, 'models', file_name)
      
      if (file.exists(file_path)) {
        # load existing model
        message(paste0('Loading existing model: ', file_name))
        ffs_models[[model_name]] <- readRDS(file_path)
      } else {
        # mark for training
        models_to_train[[model_name]] <- list(
          resp = resp,
          leaf_type = leaf_type,
          leaf_cond = leaf_cond,
          file_name = file_name
        )
      }
    }
  }
}

# train only models that don't exist yet
if (length(models_to_train) > 0) {
  
  message(paste0('\n--- ', length(models_to_train), ' models to train ---\n'))
  
  # setup parallel processing
  n_cores <- parallel::detectCores() - 2 
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  for (model_name in names(models_to_train)) {
    
    model_info <- models_to_train[[model_name]]
    resp <- model_info$resp
    leaf_type <- model_info$leaf_type
    leaf_cond <- model_info$leaf_cond
    file_name <- model_info$file_name
    
    message(paste0('\n--- Training model: ', model_name, ' ---\n'))
    
    # select predictors based on leaf type and condition
    if (leaf_type == 'deciduous' && leaf_cond == 'lon') {
      predictors <- predictors_lon_deciduous
      ctrl <- ctrl_deciduous
      response_data <- plots_lon_deciduous
    } else if (leaf_type == 'deciduous' && leaf_cond == 'loff') {
      predictors <- predictors_loff_deciduous
      ctrl <- ctrl_deciduous
      response_data <- plots_loff_deciduous
    } else if (leaf_type == 'coniferous' && leaf_cond == 'lon') {
      predictors <- predictors_lon_coniferous
      ctrl <- ctrl_coniferous
      response_data <- plots_lon_coniferous
    } else {
      predictors <- predictors_loff_coniferous
      ctrl <- ctrl_coniferous
      response_data <- plots_loff_coniferous
    }
    
    # get response variable
    response <- sf::st_drop_geometry(response_data[[resp$col]])
    
    # train model with forward feature selection
    ffs_model <- CAST::ffs(
      predictors,
      response,
      method = 'ranger',
      trControl = ctrl,
      tuneGrid = tgrid,
      num.trees = 100,
      importance = 'permutation',
      seed = 999
    )
    
    # store model in list
    ffs_models[[model_name]] <- ffs_model
    
    # save model to file
    saveRDS(
      ffs_model,
      file.path(processed_data_dir, 'models', file_name)
    )
    
    message(paste0('Model saved: ', file_name))
  }
  
  # stop parallel processing
  parallel::stopCluster(cl)
  
} else {
  message('\n--- All models already exist, no training needed ---\n')
}



# validation
#-------------------------------------------------------------------------------

# print summary of all trained models
message('\n--- Summary of all trained models ---\n')
for (model_name in names(ffs_models)) {
  message(paste0('\n', model_name, ':'))
  print(CAST::global_validation(ffs_models[[model_name]]))
}

# extract cross-validation predictions from the trained models
# create list to store all CV predictions
cv_predictions <- list()

# loop over all combinations 
# (forest inventory attribute, leaf-type, leaf-condition)
for (resp in response_vars) {
  for (leaf_type in leaf_types) {
    for (leaf_cond in leaf_conditions) {
      
      # create model name
      model_name <- paste0('ffs_rf_', resp$name, '_', leaf_cond, '_', leaf_type)
      
      # check if model exists
      if (!model_name %in% names(ffs_models)) {
        message(paste0('Model not found: ', model_name, ' - skipping'))
        next
      }
      
      # get model and extract CV predictions
      model <- ffs_models[[model_name]]
      cv_pred <- model$pred
      
      # select appropriate plot data for linking geometries
      if (leaf_type == 'deciduous' && leaf_cond == 'lon') {
        plot_data <- plots_lon_deciduous
      } else if (leaf_type == 'deciduous' && leaf_cond == 'loff') {
        plot_data <- plots_loff_deciduous
      } else if (leaf_type == 'coniferous' && leaf_cond == 'lon') {
        plot_data <- plots_lon_coniferous
      } else {
        plot_data <- plots_loff_coniferous
      }
      
      # link to the original geometries (BI plots used for training)
      cv_pred_sf <- plot_data[cv_pred$rowIndex, ] %>%
        dplyr::select(key, kspnr, abt) %>%
        dplyr::mutate(
          pred = cv_pred$pred,
          obs = cv_pred$obs
        )
      
      # store spatial predictions
      pred_name <- paste0('pred_', resp$name, '_', leaf_cond, '_', leaf_type)
      cv_predictions[[pred_name]] <- cv_pred_sf
      
      # save to file
      file_name <- paste0('pred_obsv_', resp$name, '_leaf',
                          ifelse(leaf_cond == 'lon', 'on', 'off'),
                          '_', leaf_type, '.gpkg')
      sf::st_write(
        cv_pred_sf,
        file.path(processed_data_dir, 'predictions', file_name),
        delete_dsn = T
      )
      
      message(paste0('CV predictions saved: ', file_name))
    }
  }
}

# test where differences between pred and obsv are highest using RMSE

# If you have names for your models
model_names <- names(cv_predictions)

# Loop through all 24 models
results_list <- lapply(1:length(cv_predictions), function(i) {
  
  # Merge current model predictions with bi_plots_rtk
  merged_data <- cv_predictions[[i]] %>%
    dplyr::left_join(sf::st_drop_geometry(bi_plots_rtk), by = "kspnr")
  
  # Calculate squared error and absolute error for relative RMSE
  merged_data$squared_error <- (merged_data$pred - merged_data$obs)^2
  merged_data$abs_error <- abs(merged_data$pred - merged_data$obs)
  
  # Summarize by estimated status
  summary_by_estimated <- merged_data %>%
    dplyr::group_by(estimated) %>%
    dplyr::summarise(
      n = dplyr::n(),
      RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
      mean_obs = mean(obs, na.rm = TRUE),
      rel_RMSE = sqrt(mean(squared_error, na.rm = TRUE)) / mean(obs, na.rm = TRUE) * 100,
      MAE = mean(abs_error, na.rm = TRUE),
      max_error = max(abs_error, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Summarize for all plots combined
  summary_all <- merged_data %>%
    dplyr::summarise(
      estimated = "all",
      n = dplyr::n(),
      RMSE = sqrt(mean(squared_error, na.rm = TRUE)),
      mean_obs = mean(obs, na.rm = TRUE),
      rel_RMSE = sqrt(mean(squared_error, na.rm = TRUE)) / mean(obs, na.rm = TRUE) * 100,
      MAE = mean(abs_error, na.rm = TRUE),
      max_error = max(abs_error, na.rm = TRUE)
    )
  
  # Combine both summaries
  summary_stats <- dplyr::bind_rows(summary_by_estimated, summary_all) %>%
    dplyr::mutate(model_id = i,
                  model_name = ifelse(exists("model_names"), model_names[i], paste0("model_", i)))
  
  return(summary_stats)
})

# Combine all results
all_results <- dplyr::bind_rows(results_list)

# Extract response variable from model_name
all_results <- all_results %>%
  dplyr::mutate(
    # Extract the response variable (everything between "pred_" and "_lon" or "_loff")
    response_var = gsub("pred_(.*?)_(lon|loff)_.*", "\\1", model_name),
    # Extract leaf condition (lon/loff)
    leaf_condition = ifelse(grepl("_lon_", model_name), "lon", "loff"),
    # Extract leaf type (deciduous/coniferous)
    leaf_type = ifelse(grepl("deciduous", model_name), "deciduous", "coniferous")
  )

# Aggregate by response variable and leaf type (averaged over leaf condition)
summary_by_response_leaf <- all_results %>%
  dplyr::group_by(response_var, leaf_type, estimated) %>%
  dplyr::summarise(
    n_models = dplyr::n(),
    total_plots = sum(n, na.rm = TRUE),
    avg_RMSE = mean(RMSE, na.rm = TRUE),
    avg_rel_RMSE = mean(rel_RMSE, na.rm = TRUE),
    avg_MAE = mean(MAE, na.rm = TRUE),
    avg_max_error = mean(max_error, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_by_response_leaf, n = Inf)
summary_by_response_leaf %>%
  knitr::kable(digits = 2)

# Summary separated by leaf condition (no aggregation)
summary_by_response_leaf_cond <- all_results %>%
  dplyr::group_by(response_var, leaf_type, leaf_condition, estimated) %>%
  dplyr::summarise(
    n = sum(n, na.rm = TRUE),
    RMSE = mean(RMSE, na.rm = TRUE),
    mean_obs = mean(mean_obs, na.rm = TRUE),
    rel_RMSE = mean(rel_RMSE, na.rm = TRUE),
    MAE = mean(MAE, na.rm = TRUE),
    max_error = mean(max_error, na.rm = TRUE),
    .groups = 'drop'
  )

print(summary_by_response_leaf_cond, n = Inf)
summary_by_response_leaf_cond %>%
  knitr::kable(digits = 2)

# prepare data for plotting
all_results_plot <- all_results %>%
  dplyr::mutate(
    plot_type = dplyr::case_when(
      estimated == "all" ~ "All plots",
      estimated == "yes" ~ "RTK estimated",
      estimated == "no" ~ "RTK remeasured",
      is.na(estimated) ~ "non-RTK"
    ),
    plot_type = factor(plot_type, levels = c("All plots", "RTK remeasured", "non-RTK", "RTK estimated")),
    leaf_condition = ifelse(leaf_condition == "lon", "leaf-on", "leaf-off")
  )

# plot 1: relative RMSE by response variable and leaf type (aggregated over leaf condition)
ggplot(all_results_plot, aes(x = response_var, y = rel_RMSE, fill = plot_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ leaf_type) +
  labs(title = "Relative RMSE (%) by Response Variable and Leaf Type",
       x = "Response Variable", 
       y = "Relative RMSE (%)",
       fill = "Plot Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot 2: relative RMSE for leaf-on only
ggplot(all_results_plot %>% dplyr::filter(leaf_condition == "leaf-on"), 
       aes(x = response_var, y = rel_RMSE, fill = plot_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ leaf_type) +
  labs(title = "Relative RMSE (%) - Leaf-on",
       x = "Response Variable", 
       y = "Relative RMSE (%)",
       fill = "Plot Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot 3: relative RMSE for leaf-off only
ggplot(all_results_plot %>% dplyr::filter(leaf_condition == "leaf-off"), 
       aes(x = response_var, y = rel_RMSE, fill = plot_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ leaf_type) +
  labs(title = "Relative RMSE (%) - Leaf-off",
       x = "Response Variable", 
       y = "Relative RMSE (%)",
       fill = "Plot Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##########################################################

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
























