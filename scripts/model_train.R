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



# 04: model training
#-------------------------------------------------------------------------------

# define all training datasets
# structure: list with data (sf), predictors (df), and matching pixel_centroids
predictor_start_col <- 13

training_data <- list(
  # RTK datasets
  lon_rtk_all = list(
    data = plot_metrics_lon_rtk,
    predictors = sf::st_drop_geometry(
      plot_metrics_lon_rtk[, predictor_start_col:ncol(plot_metrics_lon_rtk)]
      ),
    predpoints = pixel_centroids
  ),
  lon_rtk_deciduous = list(
    data = plots_lon_rtk_deciduous,
    predictors = sf::st_drop_geometry(
      plots_lon_rtk_deciduous[, predictor_start_col:ncol(plots_lon_rtk_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  lon_rtk_coniferous = list(
    data = plots_lon_rtk_coniferous,
    predictors = sf::st_drop_geometry(
      plots_lon_rtk_coniferous[, predictor_start_col:ncol(plots_lon_rtk_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  ),
  loff_rtk_all = list(
    data = plot_metrics_loff_rtk,
    predictors = sf::st_drop_geometry(
      plot_metrics_loff_rtk[, predictor_start_col:ncol(plot_metrics_loff_rtk)]
      ),
    predpoints = pixel_centroids
  ),
  loff_rtk_deciduous = list(
    data = plots_loff_rtk_deciduous,
    predictors = sf::st_drop_geometry(
      plots_loff_rtk_deciduous[, predictor_start_col:ncol(plots_loff_rtk_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  loff_rtk_coniferous = list(
    data = plots_loff_rtk_coniferous,
    predictors = sf::st_drop_geometry(
      plots_loff_rtk_coniferous[, predictor_start_col:ncol(plots_loff_rtk_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  ),
  # non-RTK datasets
  lon_non_rtk_all = list(
    data = plot_metrics_lon_non_rtk,
    predictors = sf::st_drop_geometry(
      plot_metrics_lon_non_rtk[, predictor_start_col:ncol(plot_metrics_lon_non_rtk)]
      ),
    predpoints = pixel_centroids
  ),
  lon_non_rtk_deciduous = list(
    data = plots_lon_non_rtk_deciduous,
    predictors = sf::st_drop_geometry(
      plots_lon_non_rtk_deciduous[, predictor_start_col:ncol(plots_lon_non_rtk_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  lon_non_rtk_coniferous = list(
    data = plots_lon_non_rtk_coniferous,
    predictors = sf::st_drop_geometry(
      plots_lon_non_rtk_coniferous[, predictor_start_col:ncol(plots_lon_non_rtk_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  ),
  loff_non_rtk_all = list(
    data = plot_metrics_loff_non_rtk,
    predictors = sf::st_drop_geometry(
      plot_metrics_loff_non_rtk[, predictor_start_col:ncol(plot_metrics_loff_non_rtk)]
      ),
    predpoints = pixel_centroids
  ),
  loff_non_rtk_deciduous = list(
    data = plots_loff_non_rtk_deciduous,
    predictors = sf::st_drop_geometry(
      plots_loff_non_rtk_deciduous[, predictor_start_col:ncol(plots_loff_non_rtk_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  loff_non_rtk_coniferous = list(
    data = plots_loff_non_rtk_coniferous,
    predictors = sf::st_drop_geometry(
      plots_loff_non_rtk_coniferous[, predictor_start_col:ncol(plots_loff_non_rtk_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  )
)

# initialize NNDM for all datasets
message('Initializing NNDM for all datasets...')
nndm_results <- lapply(names(training_data), function(name) {
  message(paste0('  Computing NNDM for: ', name))
  CAST::nndm(
    tpoints = training_data[[name]]$data,
    predpoints = training_data[[name]]$predpoints
  )
})
names(nndm_results) <- names(training_data)

# create trainControl for each dataset
train_controls <- lapply(names(nndm_results), function(name) {
  caret::trainControl(
    method = 'cv',
    index = nndm_results[[name]]$indx_train,
    indexOut = nndm_results[[name]]$indx_test,
    savePredictions = T,
    allowParallel = T
  )
})
names(train_controls) <- names(training_data)

# create grid for tuning features
tgrid <- expand.grid(
  mtry = 1:3,
  splitrule = c('variance', 'extratrees', 'maxstat'),
  min.node.size = c(5, 10, 15, 20)
)

# correlation plots for selected datasets (RTK leaf-on as example)
cor_lon_rtk_deciduous <- stats::cor(training_data$lon_rtk_deciduous$predictors, method = 'pearson')
cor_lon_rtk_coniferous <- stats::cor(training_data$lon_rtk_coniferous$predictors, method = 'pearson')
corrplot::corrplot(cor_lon_rtk_deciduous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA,
                   title = 'Correlation - leaf-on RTK deciduous')
corrplot::corrplot(cor_lon_rtk_coniferous, method = 'color', type = 'full',
                   tl.col = 'black', tl.cex = 0.6, addCoef.col = NA,
                   title = 'Correlation - leaf-on RTK coniferous')

# train random forest models
# implementing forward feature selection
# for all forest inventory attributes across all datasets
# datasets: 2 leaf conditions (lon, loff) x 2 RTK types (rtk, non_rtk) x 
#           3 leaf types (all, deciduous, coniferous) = 12 datasets
# total models: 6 responses x 12 datasets = 72 models

# define response variables
response_vars <- list(
  list(name = 'total_vol_ha', col = 'total_vol_ha'),
  list(name = 'merch_vol_ha', col = 'merch_vol_ha'),
  list(name = 'agb_ha', col = 'agb_ha'),
  list(name = 'tree_density', col = 'tree_density'),
  list(name = 'basal_area_ha', col = 'basal_area_ha'),
  list(name = 'dg', col = 'dg')
)

# dataset names for iteration
dataset_names <- names(training_data)

# create list to store all models
ffs_models <- list()

# check which models already exist
models_to_train <- list()
for (resp in response_vars) {
  for (dataset_name in dataset_names) {
    
    model_name <- paste0('ffs_rf_', resp$name, '_', dataset_name)
    file_name <- paste0('ffs_rf_', resp$name, '_', dataset_name, '.RDS')
    file_path <- file.path(processed_data_dir, 'models', file_name)
    
    if (file.exists(file_path)) {
      # load existing model
      message(paste0('Loading existing model: ', file_name))
      ffs_models[[model_name]] <- readRDS(file_path)
    } else {
      # mark for training
      models_to_train[[model_name]] <- list(
        resp = resp,
        dataset_name = dataset_name,
        file_name = file_name
      )
    }
  }
}

# print summary
message(paste0('\n--- Total models: ', length(response_vars) * length(dataset_names), ' ---'))
message(paste0('--- Models already trained: ', length(ffs_models), ' ---'))
message(paste0('--- Models to train: ', length(models_to_train), ' ---\n'))

# train only models that don't exist yet
if (length(models_to_train) > 0) {
  
  # setup parallel processing
  n_cores <- parallel::detectCores() - 2 
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)
  
  for (model_name in names(models_to_train)) {
    
    model_info <- models_to_train[[model_name]]
    resp <- model_info$resp
    dataset_name <- model_info$dataset_name
    file_name <- model_info$file_name
    
    message(paste0('\n--- Training model: ', model_name, ' ---\n'))
    
    # get data for this dataset
    predictors <- training_data[[dataset_name]]$predictors
    response_data <- training_data[[dataset_name]]$data
    ctrl <- train_controls[[dataset_name]]
    
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

# loop over all combinations (response variable x dataset)
for (resp in response_vars) {
  for (dataset_name in dataset_names) {
    
    # create model name
    model_name <- paste0('ffs_rf_', resp$name, '_', dataset_name)
    
    # check if model exists
    if (!model_name %in% names(ffs_models)) {
      message(paste0('Model not found: ', model_name, ' - skipping'))
      next
    }
    
    # get model and extract CV predictions
    model <- ffs_models[[model_name]]
    cv_pred <- model$pred
    
    # get plot data for linking geometries
    plot_data <- training_data[[dataset_name]]$data
    
    # link to the original geometries (BI plots used for training)
    cv_pred_sf <- plot_data[cv_pred$rowIndex, ] %>%
      dplyr::select(key, kspnr, abt) %>%
      dplyr::mutate(
        pred = cv_pred$pred,
        obs = cv_pred$obs
      )
    
    # store spatial predictions
    pred_name <- paste0('pred_', resp$name, '_', dataset_name)
    cv_predictions[[pred_name]] <- cv_pred_sf
    
    # save to file
    file_name <- paste0('pred_obsv_', resp$name, '_', dataset_name, '.gpkg')
    sf::st_write(
      cv_pred_sf,
      file.path(processed_data_dir, 'predictions', file_name),
      delete_dsn = T
    )
    
    message(paste0('CV predictions saved: ', file_name))
  }
}

# compare RTK vs non-RTK model performance
# for corresponding combinations (same leaf condition and leaf type)

# get model names
model_names <- names(cv_predictions)

# calculate validation metrics for all models
results_list <- lapply(seq_along(cv_predictions), function(i) {
  
  pred_data <- cv_predictions[[i]]
  
  # calculate squared error and absolute error
  pred_data$squared_error <- (pred_data$pred - pred_data$obs)^2
  pred_data$abs_error <- abs(pred_data$pred - pred_data$obs)
  
  # calculate summary statistics
  summary_stats <- pred_data %>%
    sf::st_drop_geometry() %>%
    dplyr::summarise(
      n = dplyr::n(),
      RMSE = sqrt(mean(squared_error, na.rm = T)),
      mean_obs = mean(obs, na.rm = T),
      rel_RMSE = sqrt(mean(squared_error, na.rm = T)) / mean(obs, na.rm = T) * 100,
      MAE = mean(abs_error, na.rm = T),
      max_error = max(abs_error, na.rm = T)
    ) %>%
    dplyr::mutate(
      model_id = i,
      model_name = model_names[i]
    )
  
  return(summary_stats)
})

# combine all results
all_results <- dplyr::bind_rows(results_list)

# extract components from model_name
# format: pred_{response}_{leaf_cond}_{positioning}_{leaf_type}
all_results <- all_results %>%
  
  dplyr::mutate(
    
    # extract response variable (everything between 'pred_' and '_lon_' or '_loff_')
    response_var = gsub('pred_(.+)_(lon|loff)_.*', '\\1', model_name),
    
    # extract positioning method (rtk vs non-rtk)
    positioning = dplyr::case_when(
      grepl('_non_rtk_', model_name) ~ 'non_rtk',
      grepl('_rtk_', model_name) ~ 'rtk',
      T ~ NA_character_
    ),
    
    # extract leaf condition (lon/loff)
    leaf_condition = dplyr::case_when(
      grepl('_lon_', model_name) ~ 'lon',
      grepl('_loff_', model_name) ~ 'loff',
      T ~ NA_character_
    ),
    
    # extract leaf type (all/deciduous/coniferous)
    leaf_type = dplyr::case_when(
      grepl('_deciduous$', model_name) ~ 'deciduous',
      grepl('_coniferous$', model_name) ~ 'coniferous',
      grepl('_all$', model_name) ~ 'all',
      T ~ NA_character_
    )
  )

# display full results table
all_results %>%
  dplyr::select(response_var, positioning, leaf_condition, leaf_type, n, RMSE, rel_RMSE, MAE) %>%
  dplyr::arrange(response_var, positioning, leaf_condition, leaf_type) %>%
  knitr::kable(digits = 2)

# compare RTK vs non-RTK for corresponding combinations
pos_comparison <- all_results %>%
  dplyr::select(response_var, positioning, leaf_condition, leaf_type, n, RMSE, rel_RMSE, MAE) %>%
  tidyr::pivot_wider(
    id_cols = c(response_var, leaf_condition, leaf_type),
    names_from = positioning,
    values_from = c(n, RMSE, rel_RMSE, MAE),
    names_glue = '{positioning}_{.value}'
  ) %>%
  dplyr::mutate(
    # calculate differences (rtk - non-rtk)
    diff_RMSE = rtk_RMSE - non_rtk_RMSE,
    diff_rel_RMSE = rtk_rel_RMSE - non_rtk_rel_RMSE,
    diff_MAE = rtk_MAE - non_rtk_MAE,
    # determine which is better (lower error = better)
    better_positioning = ifelse(rtk_rel_RMSE < non_rtk_rel_RMSE, 'rtk', 'non_rtk'),
    # relative improvement (%)
    rel_improvement = (non_rtk_rel_RMSE - rtk_rel_RMSE) / non_rtk_rel_RMSE * 100
  )

# display RTK vs non-RTK comparison
pos_comparison %>%
  dplyr::select(response_var, leaf_condition, leaf_type, 
                rtk_rel_RMSE, non_rtk_rel_RMSE, diff_rel_RMSE, 
                better_positioning, rel_improvement) %>%
  dplyr::arrange(response_var, leaf_condition, leaf_type) %>%
  knitr::kable(digits = 2)

# summary: how often is RTK better than non-RTK?
pos_summary <- pos_comparison %>%
  dplyr::group_by(leaf_condition, leaf_type) %>%
  dplyr::summarise(
    n_comparisons = dplyr::n(),
    n_rtk_better = sum(better_positioning == 'rtk'),
    n_non_rtk_better = sum(better_positioning == 'non_rtk'),
    pct_rtk_better = n_rtk_better / n_comparisons * 100,
    mean_rel_improvement = mean(rel_improvement, na.rm = T),
    .groups = 'drop'
  )

pos_summary %>%
  knitr::kable(digits = 2)

# summary by response variable
pos_summary_by_response <- pos_comparison %>%
  dplyr::group_by(response_var) %>%
  dplyr::summarise(
    n_comparisons = dplyr::n(),
    n_rtk_better = sum(better_positioning == 'rtk'),
    n_non_rtk_better = sum(better_positioning == 'non_rtk'),
    pct_rtk_better = n_rtk_better / n_comparisons * 100,
    mean_rtk_rel_RMSE = mean(rtk_rel_RMSE, na.rm = T),
    mean_non_rtk_rel_RMSE = mean(non_rtk_rel_RMSE, na.rm = T),
    mean_rel_improvement = mean(rel_improvement, na.rm = T),
    .groups = 'drop'
  )

pos_summary_by_response %>%
  knitr::kable(digits = 2)

# overall summary
overall_summary <- pos_comparison %>%
  dplyr::summarise(
    total_comparisons = dplyr::n(),
    n_rtk_better = sum(better_positioning == 'rtk'),
    n_non_rtk_better = sum(better_positioning == 'non_rtk'),
    pct_rtk_better = n_rtk_better / total_comparisons * 100,
    mean_rtk_rel_RMSE = mean(rtk_rel_RMSE, na.rm = T),
    mean_non_rtk_rel_RMSE = mean(non_rtk_rel_RMSE, na.rm = T),
    mean_rel_improvement = mean(rel_improvement, na.rm = T)
  )

message('\n--- Overall RTK vs non-RTK comparison ---')
overall_summary %>%
  knitr::kable(digits = 2)

# define forest inventory attribute names (response variables)
forest_inv_names <- c(
  'agb_ha' = 'AGB',
  'total_vol_ha' = 'VOL_Tot',
  'merch_vol_ha' = 'VOL_Merch',
  'basal_area_ha' = 'BA',
  'dg' = 'QMD',
  'tree_density' = 'TD'
)

# define order for forest inventory attributes
forest_inv_order <- c('AGB', 'VOL_Tot', 'VOL_Merch', 'BA', 'QMD', 'TD')

# prepare data for plotting
data_visual <- all_results %>%
  dplyr::mutate(
    leaf_condition_label = ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
    response_var_label = factor(forest_inv_names[response_var], levels = forest_inv_order)
  ) %>%
  dplyr::group_by(response_var, response_var_label, positioning, leaf_type) %>%
  dplyr::summarise(
    mean_rel_RMSE = mean(rel_RMSE, na.rm = T),
    lon_rel_RMSE = rel_RMSE[leaf_condition == 'lon'],
    loff_rel_RMSE = rel_RMSE[leaf_condition == 'loff'],
    min_rel_RMSE = min(rel_RMSE, na.rm = T),
    max_rel_RMSE = max(rel_RMSE, na.rm = T),
    # labels for error bar ends
    min_label = ifelse(lon_rel_RMSE <= loff_rel_RMSE, 'leaf-on', 'leaf-off'),
    max_label = ifelse(lon_rel_RMSE >= loff_rel_RMSE, 'leaf-on', 'leaf-off'),
    .groups = 'drop'
  )

# plot: rel RMSE by response variable, positioning (rtk/non-rtk) as bars
# error bars show range between leaf-on and leaf-off
# faceted by leaf type (all, deciduous, coniferous)
ggplot(data_visual, aes(x = response_var_label, y = mean_rel_RMSE, fill = positioning)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = min_rel_RMSE, ymax = max_rel_RMSE),
                position = position_dodge(width = 0.8), width = 0.25) +
  geom_text(aes(y = max_rel_RMSE, label = max_label), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 2.5) +
  geom_text(aes(y = min_rel_RMSE, label = min_label), 
            position = position_dodge(width = 0.8), vjust = 1.5, size = 2.5) +
  facet_wrap(~ leaf_type, ncol = 3) +
  scale_fill_manual(values = c('rtk' = '#606060', 'non_rtk' = '#B0B0B0'),
                    labels = c('rtk' = 'RTK', 'non_rtk' = 'non-RTK')) +
  labs(x = '', 
       y = 'Relative RMSE (%)',
       fill = 'Positioning') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = 'lightgrey'),
        strip.text = element_text(face = 'bold'))

# plot: same but separated by leaf condition 
# (6 panels: 2 leaf conditions x 3 leaf types)
ggplot(all_results %>% 
         dplyr::mutate(leaf_condition_label = ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
                       response_var_label = factor(forest_inv_names[response_var], levels = forest_inv_order)), 
       aes(x = response_var_label, y = rel_RMSE, fill = positioning)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(leaf_condition_label ~ leaf_type) +
  scale_fill_manual(values = c('rtk' = '#606060', 'non_rtk' = '#B0B0B0'),
                    labels = c('rtk' = 'RTK', 'non_rtk' = 'non-RTK')) +
  labs(x = '', 
       y = 'Relative RMSE (%)',
       fill = 'Positioning') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = 'lightgrey'),
        strip.text = element_text(face = 'bold'))

# plot: difference between RTK and non-RTK (improvement)
ggplot(pos_comparison %>%
         dplyr::mutate(leaf_condition = ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
                       response_var_label = factor(forest_inv_names[response_var], levels = forest_inv_order)), 
       aes(x = response_var_label, y = rel_improvement, fill = rel_improvement > 0)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  facet_grid(leaf_condition ~ leaf_type) +
  scale_fill_manual(values = c('TRUE' = '#606060', 'FALSE' = '#B0B0B0'),
                    labels = c('TRUE' = 'RTK better', 'FALSE' = 'non-RTK better'),
                    name = '') +
  labs(x = '', 
       y = 'Relative Improvement (%)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = 'lightgrey'),
        strip.text = element_text(face = 'bold'))

##########################################################

# plot predicted vs. observed with CV predictions
# combining RTK and non-RTK in one plot
# 36 plots: 6 response vars x 2 leaf conditions x 3 leaf types

# combine all CV predictions into one data frame
all_cv_predictions <- do.call(rbind, lapply(names(cv_predictions), function(name) {
  cv_predictions[[name]] %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(model_name = name)
}))

# extract components from model_name
all_cv_predictions <- all_cv_predictions %>%
  dplyr::mutate(
    response_var = gsub('pred_(.+)_(lon|loff)_.*', '\\1', model_name),
    positioning = dplyr::case_when(
      grepl('_non_rtk_', model_name) ~ 'non_rtk',
      grepl('_rtk_', model_name) ~ 'rtk',
      T ~ NA_character_
    ),
    leaf_condition = dplyr::case_when(
      grepl('_lon_', model_name) ~ 'leaf-on',
      grepl('_loff_', model_name) ~ 'leaf-off',
      T ~ NA_character_
    ),
    leaf_type = dplyr::case_when(
      grepl('_deciduous$', model_name) ~ 'deciduous',
      grepl('_coniferous$', model_name) ~ 'coniferous',
      grepl('_all$', model_name) ~ 'all',
      T ~ NA_character_
    )
  )

# create one plot per response variable
for (resp in unique(all_cv_predictions$response_var)) {
  
  plot_data_resp <- all_cv_predictions %>%
    dplyr::filter(response_var == resp) %>%
    dplyr::mutate(
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous'))
    )
  
  if (nrow(plot_data_resp) == 0) next
  
  # get forest inventory attribute name for response variable
  forest_inv_label <- forest_inv_names[resp]
  
  # calculate axis limits for equal scales
  axis_max <- max(c(plot_data_resp$obs, plot_data_resp$pred), na.rm = T)
  axis_min <- min(c(plot_data_resp$obs, plot_data_resp$pred), na.rm = T)
  
  p <- ggplot(plot_data_resp, aes(x = obs, y = pred, color = positioning, shape = positioning)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linewidth = 1, color = 'red', linetype = 'dashed') +
    geom_smooth(method = 'lm', se = F, linewidth = 1) +
    facet_grid(leaf_condition ~ leaf_type) +
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    scale_color_manual(values = c('rtk' = '#404040', 'non_rtk' = '#909090'),
                       labels = c('rtk' = 'RTK', 'non_rtk' = 'non-RTK')) +
    scale_shape_manual(values = c('rtk' = 16, 'non_rtk' = 17),
                       labels = c('rtk' = 'RTK', 'non_rtk' = 'non-RTK')) +
    labs(title = paste0('Predicted vs. Observed: ', forest_inv_label),
         x = 'Observed', 
         y = 'Predicted',
         color = 'Positioning', 
         shape = 'Positioning') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          legend.position = 'bottom',
          strip.background = element_rect(fill = 'lightgrey'),
          strip.text = element_text(face = 'bold'))
  
  print(p)
}
























