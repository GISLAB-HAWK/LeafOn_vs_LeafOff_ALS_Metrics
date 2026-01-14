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
# and calculated metrics (leaf-on and leaf-off, RTK and non-RTK)
plot_metrics_lon_rtk <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon_rtk.gpkg')
)

plot_metrics_loff_rtk <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff_rtk.gpkg')
)

plot_metrics_lon_non_rtk <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon_non_rtk.gpkg')
)

plot_metrics_loff_non_rtk <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff_non_rtk.gpkg')
)

head(plot_metrics_lon_rtk)
head(plot_metrics_loff_rtk)
head(plot_metrics_lon_non_rtk)
head(plot_metrics_loff_non_rtk)
str(plot_metrics_lon_rtk)
str(plot_metrics_loff_rtk)

# read extent of the leaf-off dataset
ext_loff <- sf::st_read(
  file.path(processed_data_dir, 'pc_leafoff_2024', 'leafoff_ppm4.vpc')
  )

# read pixel centroids from w2w-metrics raster
pixel_centroids <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'pixel_centroids.gpkg')
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

# filter plots by dominant leaf type (leaf-on RTK)
plots_lon_rtk_deciduous <- plot_metrics_lon_rtk %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_lon_rtk_coniferous <- plot_metrics_lon_rtk %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# filter plots by dominant leaf type (leaf-off RTK)
plots_loff_rtk_deciduous <- plot_metrics_loff_rtk %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_loff_rtk_coniferous <- plot_metrics_loff_rtk %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# filter plots by dominant leaf type (leaf-on non-RTK)
plots_lon_non_rtk_deciduous <- plot_metrics_lon_non_rtk %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_lon_non_rtk_coniferous <- plot_metrics_lon_non_rtk %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# filter plots by dominant leaf type (leaf-off non-RTK)
plots_loff_non_rtk_deciduous <- plot_metrics_loff_non_rtk %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_loff_non_rtk_coniferous <- plot_metrics_loff_non_rtk %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# df versions (leaf-on RTK)
plots_lon_rtk_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_lon_rtk))
plots_lon_rtk_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_rtk_deciduous))
plots_lon_rtk_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_rtk_coniferous))

# df versions (leaf-off RTK)
plots_loff_rtk_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_loff_rtk))
plots_loff_rtk_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_rtk_deciduous))
plots_loff_rtk_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_rtk_coniferous))

# df versions (leaf-on non-RTK)
plots_lon_non_rtk_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_lon_non_rtk))
plots_lon_non_rtk_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_non_rtk_deciduous))
plots_lon_non_rtk_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_non_rtk_coniferous))

# df versions (leaf-off non-RTK)
plots_loff_non_rtk_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_loff_non_rtk))
plots_loff_non_rtk_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_non_rtk_deciduous))
plots_loff_non_rtk_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_non_rtk_coniferous))



# 03: CV testing
#-------------------------------------------------------------------------------

# define data combinations for CV comparison
# (RTK vs non-RTK) x (deciduous vs coniferous vs all)
# leaf-on vs leaf-off not differentiated as plot positions are the same

data_list <- list(
  rtk_all = list(data = plot_metrics_lon_rtk, predpoints = pixel_centroids),
  rtk_deciduous = list(data = plots_lon_rtk_deciduous, predpoints = pixel_centroids_deciduous),
  rtk_coniferous = list(data = plots_lon_rtk_coniferous, predpoints = pixel_centroids_coniferous),
  non_rtk_all = list(data = plot_metrics_lon_non_rtk, predpoints = pixel_centroids),
  non_rtk_deciduous = list(data = plots_lon_non_rtk_deciduous, predpoints = pixel_centroids_deciduous),
  non_rtk_coniferous = list(data = plots_lon_non_rtk_coniferous, predpoints = pixel_centroids_coniferous)
)

# function to compute CV methods for the different datasets
compute_cv_geodist <- function(data, predpoints, modeldomain) {
  n <- nrow(data)
  
  # LOO CV
  loo_cv <- caret::createFolds(1:n, k = n, returnTrain = F)
  geodist_loo <- CAST::geodist(
    x = data,
    modeldomain = sf::st_zm(modeldomain),
    cvfolds = loo_cv
  )
  
  # 5-fold CV
  fold5 <- caret::createFolds(1:n, k = 5, returnTrain = F)
  geodist_fold5 <- CAST::geodist(
    x = data,
    modeldomain = sf::st_zm(modeldomain),
    cvfolds = fold5
  )
  
  # NNDM LOO CV
  nndm_result <- CAST::nndm(
    tpoints = data,
    predpoints = predpoints
  )
  
  list(
    loo_cv = loo_cv,
    geodist_loo = geodist_loo,
    fold5 = fold5,
    geodist_fold5 = geodist_fold5,
    nndm = nndm_result
  )
}

# compute CV for all datasets
cv_results <- lapply(data_list, function(item) {
  compute_cv_geodist(item$data, item$predpoints, ext_loff)
})

# plot ECDF functions for all combinations
# LOO CV ECDF plot
loo_plots <- lapply(names(cv_results), function(name) {
  plot(cv_results[[name]]$geodist_loo, stat = 'ecdf') +
    ggplot2::ggtitle(paste0('LOO CV - ', name)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
})
names(loo_plots) <- names(cv_results)

# 5-fold CV ECDF plot
fold5_plots <- lapply(names(cv_results), function(name) {
  plot(cv_results[[name]]$geodist_fold5, stat = 'ecdf') +
    ggplot2::ggtitle(paste0('5-fold CV - ', name)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
})
names(fold5_plots) <- names(cv_results)

# NNDM ECDF plot
nndm_plots <- lapply(names(cv_results), function(name) {
  plot(cv_results[[name]]$nndm, type = 'simple') +
    ggplot2::ggtitle(paste0('NNDM - ', name)) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 10))
})
names(nndm_plots) <- names(cv_results)

# combine ECDF plots
gridExtra::grid.arrange(grobs = loo_plots, ncol = 3, 
                        top = 'LOO CV - ECDF Comparison')

gridExtra::grid.arrange(grobs = fold5_plots, ncol = 3, 
                        top = '5-fold CV - ECDF Comparison')

gridExtra::grid.arrange(grobs = nndm_plots, ncol = 3, 
                        top = 'NNDM - ECDF Comparison')

# plot showing how NNDM works - RTK coniferous example
# cv iteration with the most excluded plots
nndm_rtk_coniferous <- cv_results$rtk_coniferous$nndm
id_plot_coniferous <- which.max(sapply(nndm_rtk_coniferous$indx_exclude, length))

rtk_coniferous_plot <- plots_lon_rtk_coniferous
rtk_coniferous_plot$set <- ""
rtk_coniferous_plot$set[nndm_rtk_coniferous$indx_train[[id_plot_coniferous]]] <- 'train'
rtk_coniferous_plot$set[nndm_rtk_coniferous$indx_exclude[[id_plot_coniferous]]] <- 'exclude'
rtk_coniferous_plot$set[nndm_rtk_coniferous$indx_test[[id_plot_coniferous]]] <- 'test'
rtk_coniferous_plot <- rtk_coniferous_plot[order(rtk_coniferous_plot$set),]

ggplot() +
  geom_sf(data = ext_loff, fill = 'grey', alpha = 0.1) +
  geom_sf(data = rtk_coniferous_plot, aes(col = set)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  ggtitle('NNDM visualization - RTK coniferous',
          subtitle = 'CV iteration with most excluded plots')

# plot showing how NNDM works - RTK deciduous example
# cv iteration with the most excluded plots
nndm_rtk_deciduous <- cv_results$rtk_deciduous$nndm
id_plot_deciduous <- which.max(sapply(nndm_rtk_deciduous$indx_exclude, length))

rtk_deciduous_plot <- plots_lon_rtk_deciduous
rtk_deciduous_plot$set <- ""
rtk_deciduous_plot$set[nndm_rtk_deciduous$indx_train[[id_plot_deciduous]]] <- 'train'
rtk_deciduous_plot$set[nndm_rtk_deciduous$indx_exclude[[id_plot_deciduous]]] <- 'exclude'
rtk_deciduous_plot$set[nndm_rtk_deciduous$indx_test[[id_plot_deciduous]]] <- 'test'
rtk_deciduous_plot <- rtk_deciduous_plot[order(rtk_deciduous_plot$set),]

ggplot() +
  geom_sf(data = ext_loff, fill = 'grey', alpha = 0.1) +
  geom_sf(data = rtk_deciduous_plot, aes(col = set)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  ggtitle('NNDM visualization - RTK deciduous',
          subtitle = 'CV iteration with most excluded plots')

# fit models and estimate their performance (RTK and leaf-on only)
data_df_list <- list(
  rtk_all = plots_lon_rtk_all_df,
  rtk_deciduous = plots_lon_rtk_deciduous_df,
  rtk_coniferous = plots_lon_rtk_coniferous_df
)

# function to fit RF models with different CV methods
fit_cv_models <- function(data_df, nndm_result, predictor_start_col = 13) {
  
  # LOO CV
  loo_ctrl <- caret::trainControl(method = 'LOOCV', savePredictions = T)
  loo_mod <- caret::train(
    data_df[, predictor_start_col:ncol(data_df)],
    data_df[, 'total_vol_ha'],
    method = 'rf',
    importance = F,
    trControl = loo_ctrl,
    ntree = 100,
    tuneLength = 1
  )
  
  # 5-fold CV
  fold5_ctrl <- caret::trainControl(method = 'cv', number = 5, savePredictions = T)
  fold5_mod <- caret::train(
    data_df[, predictor_start_col:ncol(data_df)],
    data_df[, 'total_vol_ha'],
    method = 'rf',
    importance = F,
    trControl = fold5_ctrl,
    ntree = 100,
    tuneLength = 1
  )
  
  # NNDM LOO CV
  nndm_ctrl <- caret::trainControl(
    method = 'cv',
    index = nndm_result$indx_train,
    indexOut = nndm_result$indx_test,
    savePredictions = T
  )
  nndm_mod <- caret::train(
    data_df[, predictor_start_col:ncol(data_df)],
    data_df[, 'total_vol_ha'],
    method = 'rf',
    importance = F,
    trControl = nndm_ctrl,
    ntree = 100,
    tuneLength = 1
  )
  
  list(
    loo = loo_mod,
    fold5 = fold5_mod,
    nndm = nndm_mod
  )
}

# fit models for all datasets
cv_models <- lapply(names(data_df_list), function(name) {
  message(paste0('Fitting models for: ', name))
  fit_cv_models(
    data_df = data_df_list[[name]],
    nndm_result = cv_results[[name]]$nndm
  )
})
names(cv_models) <- names(data_df_list)

# create results table for all datasets
cv_validation_results <- do.call(rbind, lapply(names(cv_models), function(name) {
  rbind(
    data.frame(
      data = name,
      validation = 'LOO CV',
      t(as.data.frame(CAST::global_validation(cv_models[[name]]$loo)))
    ),
    data.frame(
      data = name,
      validation = '5-fold CV',
      t(as.data.frame(CAST::global_validation(cv_models[[name]]$fold5)))
    ),
    data.frame(
      data = name,
      validation = 'NNDM LOO CV',
      t(as.data.frame(CAST::global_validation(cv_models[[name]]$nndm)))
    )
  )
}))

cv_validation_results |> 
  knitr::kable(digits = 2, row.names = F)



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
# for all forest inventory attributes, leaf types, and leaf conditions

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

# test where differences between pred and obsv are highest using rel. RMSE

# get model names
model_names <- names(cv_predictions)

# loop through all 24 models
results_list <- lapply(1:length(cv_predictions), function(i) {
  
  # merge current model predictions with bi_plots_rtk
  merged_data <- cv_predictions[[i]] %>%
    dplyr::left_join(sf::st_drop_geometry(bi_plots_rtk), by = "kspnr")
  
  # calculate squared error and absolute error
  merged_data$squared_error <- (merged_data$pred - merged_data$obs)^2
  merged_data$abs_error <- abs(merged_data$pred - merged_data$obs)
  
  # summarize by estimated status
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
  
  # summarize for all plots combined
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
  
  # combine both summaries
  summary_stats <- dplyr::bind_rows(summary_by_estimated, summary_all) %>%
    dplyr::mutate(model_id = i,
                  model_name = ifelse(exists("model_names"), model_names[i], paste0("model_", i)))
  
  return(summary_stats)
})

# combine all results
all_results <- dplyr::bind_rows(results_list)

# extract response variable from model_name
all_results <- all_results %>%
  dplyr::mutate(
    # extract the response variable (everything between "pred_" and "_lon" or "_loff")
    response_var = gsub("pred_(.*?)_(lon|loff)_.*", "\\1", model_name),
    # extract leaf condition (lon/loff)
    leaf_condition = ifelse(grepl("_lon_", model_name), "lon", "loff"),
    # extract leaf type (deciduous/coniferous)
    leaf_type = ifelse(grepl("deciduous", model_name), "deciduous", "coniferous")
  )

# aggregate by response variable and leaf type (averaged over leaf condition)
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

summary_by_response_leaf %>%
  knitr::kable(digits = 2)

# summary separated by leaf condition (no aggregation)
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

summary_by_response_leaf_cond %>%
  knitr::kable(digits = 2)

# calculate absolute difference between leaf-on and leaf-off rel. RMSE
diff_lon_loff <- all_results_plot %>%
  sf::st_drop_geometry() %>%
  tidyr::pivot_wider(
    id_cols = c(response_var, leaf_type, plot_type),
    names_from = leaf_condition,
    values_from = rel_RMSE
  ) %>%
  dplyr::mutate(
    abs_diff = abs(`leaf-on` - `leaf-off`),
    better_condition = ifelse(`leaf-on` < `leaf-off`, "leaf-on", "leaf-off")
  )

diff_lon_loff %>%
  dplyr::select(response_var, leaf_type, plot_type, `leaf-on`, `leaf-off`, abs_diff, better_condition) %>%
  knitr::kable(digits = 2)

# summary: mean difference by leaf type and plot type across all response variables
diff_summary_by_leaf_type <- diff_lon_loff %>%
  dplyr::filter(plot_type %in% c("All plots", "RTK remeasured")) %>%
  dplyr::group_by(leaf_type, plot_type) %>%
  dplyr::summarise(
    mean_abs_diff = mean(abs_diff, na.rm = TRUE),
    sd_abs_diff = sd(abs_diff, na.rm = TRUE),
    n_leaf_on_better = sum(better_condition == "leaf-on"),
    n_leaf_off_better = sum(better_condition == "leaf-off"),
    .groups = 'drop'
  )

diff_summary_by_leaf_type %>%
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
# with error bars showing range between leaf-on and leaf-off
all_results_summary <- all_results_plot %>%
  dplyr::group_by(response_var, leaf_type, plot_type) %>%
  dplyr::summarise(
    mean_rel_RMSE = mean(rel_RMSE, na.rm = TRUE),
    lon_rel_RMSE = rel_RMSE[leaf_condition == "leaf-on"],
    loff_rel_RMSE = rel_RMSE[leaf_condition == "leaf-off"],
    min_rel_RMSE = min(rel_RMSE, na.rm = TRUE),
    max_rel_RMSE = max(rel_RMSE, na.rm = TRUE),
    # determine which condition is at min/max
    min_label = ifelse(lon_rel_RMSE <= loff_rel_RMSE, "leaf-on", "leaf-off"),
    max_label = ifelse(lon_rel_RMSE >= loff_rel_RMSE, "leaf-on", "leaf-off"),
    .groups = 'drop'
  )

ggplot(all_results_summary, aes(x = response_var, y = mean_rel_RMSE, fill = plot_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = min_rel_RMSE, ymax = max_rel_RMSE),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_text(aes(y = max_rel_RMSE, label = max_label), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 2.5) +
  geom_text(aes(y = min_rel_RMSE, label = min_label), 
            position = position_dodge(width = 0.8), vjust = 1.5, size = 2.5) +
  facet_wrap(~ leaf_type) +
  labs(title = "Relative RMSE (%) by Forest Invenotry Attribute and Dominant Leaf Type",
       subtitle = "Error bars show range between leaf-on and leaf-off",
       x = "", 
       y = "Relative RMSE (%)",
       fill = "Plot Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot 2: relative RMSE separated by leaf condition and leaf type (4 panels)
ggplot(all_results_plot, aes(x = response_var, y = rel_RMSE, fill = plot_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(leaf_condition ~ leaf_type) +
  labs(title = "Relative RMSE (%) by Leaf Condition and Dominant Leaf Type",
       x = "", 
       y = "Relative RMSE (%)",
       fill = "Plot Type") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot: absolute difference between leaf-on and leaf-off (All plots and RTK remeasured only)
ggplot(diff_lon_loff %>% dplyr::filter(plot_type %in% c("All plots", "RTK remeasured")), 
       aes(x = response_var, y = abs_diff, fill = better_condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(plot_type ~ leaf_type) +
  scale_fill_manual(values = c("leaf-on" = "#2E7D32", "leaf-off" = "#1565C0"),
                    name = "Better condition") +
  labs(title = "Absolute Difference in Relative RMSE between Leaf-on and Leaf-off",
       subtitle = "Bar color indicates which condition has lower (better) rRMSE",
       x = "", 
       y = "Absolute Difference in rel. RMSE (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot: compare leaf-on vs leaf-off differences between deciduous and coniferous
# for All plots and RTK remeasured
ggplot(diff_lon_loff %>% dplyr::filter(plot_type %in% c("All plots", "RTK remeasured")), 
       aes(x = leaf_type, y = abs_diff, fill = leaf_type)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = round(abs_diff, 1)), vjust = -0.5, size = 3) +
  facet_grid(plot_type ~ response_var) +
  scale_fill_manual(values = c("deciduous" = "#8BC34A", "coniferous" = "#4CAF50")) +
  labs(title = "Leaf-on vs Leaf-off Difference: Deciduous vs Coniferous",
       subtitle = "Absolute difference in rel. RMSE (%)",
       x = "", 
       y = "Absolute Difference in rel. RMSE (%)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

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
























