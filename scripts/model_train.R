#-------------------------------------------------------------------------------
# Name:         model_train.R
# Description:  Script trains random forest models for predicting forest
#               inventory attributes: above-ground biomass (agb_ha, Mg/ha)
#               and stem density (tree_density, n/ha).
#               ALS-based metrics previously derived in forest inventory plots
#               (BI plots) are used as predictors, once with the base ALS
#               metrics only and once with the structural complexity metrics
#               added, to quantify whether the latter are worth including.
#               Two models are trained per combination, one using the metrics
#               calculated from the leaf-on dataset, and one using the metrics
#               calculated from the leaf-off dataset.
#               Nearest Neighbour Distance Matching (NNDM)
#               Leave-One-Out Cross Validation (LOO CV) is used as a
#               spatial map validation method.
#               The models are trained separately for deciduous and coniferous
#               dominated plots, and also one time not separating between them.
#               Cross-seasonal transferability is also assessed: each fitted
#               leaf-on model is applied to the leaf-off metrics and vice versa,
#               for all/deciduous/coniferous plots, to test whether the
#               metric-to-attribute relationships hold across seasons.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read forest inventory plots (BI) with their attributes
# and calculated ALS metrics (leaf-on and leaf-off)
plot_metrics_lon <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon_new.gpkg')
)

plot_metrics_loff<- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff_new.gpkg')
)

head(plot_metrics_lon)
head(plot_metrics_loff)
str(plot_metrics_lon)
str(plot_metrics_loff)

# read extent of the leaf-off dataset
ext_loff <- sf::st_read(
  file.path(raw_data_dir, 'pc_leafoff_2024', 'leafoff.vpc')
)

# read pixel centroids from w2w-metrics raster
pixel_centroids <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'pixel_centroids.gpkg')
  )

# split pixel centroids based on dominant leaf type
pixel_centroids_deciduous <- pixel_centroids %>% 
  dplyr::filter(leaf_type_code == 1)
pixel_centroids_coniferous <- pixel_centroids %>% 
  dplyr::filter(leaf_type_code == 2)

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

# filter plots by dominant leaf type (leaf-on)
plots_lon_deciduous <- plot_metrics_lon %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_lon_coniferous <- plot_metrics_lon %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# filter plots by dominant leaf type (leaf-off)
plots_loff_deciduous <- plot_metrics_loff %>%
  dplyr::filter(dominant_leaf_type == 'deciduous')
plots_loff_coniferous <- plot_metrics_loff %>%
  dplyr::filter(dominant_leaf_type == 'coniferous')

# df versions (leaf-on)
plots_lon_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_lon))
plots_lon_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_deciduous))
plots_lon_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_lon_coniferous))

# df versions (leaf-off)
plots_loff_all_df <- as.data.frame(sf::st_drop_geometry(plot_metrics_loff))
plots_loff_deciduous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_deciduous))
plots_loff_coniferous_df <- as.data.frame(sf::st_drop_geometry(plots_loff_coniferous))



# 03: CV testing
#-------------------------------------------------------------------------------

# define data combinations for CV comparison
# (deciduous vs coniferous vs all)
# leaf-on vs leaf-off not differentiated as plot positions are the same

data_list <- list(
  all = list(data = plot_metrics_lon, predpoints = pixel_centroids),
  deciduous = list(data = plots_lon_deciduous, predpoints = pixel_centroids_deciduous),
  coniferous = list(data = plots_lon_coniferous, predpoints = pixel_centroids_coniferous)
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

# plot showing how NNDM works - coniferous example
# cv iteration with the most excluded plots
nndm_coniferous <- cv_results$coniferous$nndm
id_plot_coniferous <- which.max(sapply(nndm_coniferous$indx_exclude, length))

coniferous_plot <- plots_lon_coniferous
coniferous_plot$set <- ""
coniferous_plot$set[nndm_coniferous$indx_train[[id_plot_coniferous]]] <- 'train'
coniferous_plot$set[nndm_coniferous$indx_exclude[[id_plot_coniferous]]] <- 'exclude'
coniferous_plot$set[nndm_coniferous$indx_test[[id_plot_coniferous]]] <- 'test'
coniferous_plot <- coniferous_plot[order(coniferous_plot$set),]

p_coniferous <- ggplot() +
  geom_sf(data = ext_loff, fill = 'grey', alpha = 0.1) +
  geom_sf(data = coniferous_plot, aes(col = set)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  ggtitle('NNDM visualization - coniferous',
          subtitle = 'CV iteration with most excluded plots')

# plot showing how NNDM works - deciduous example
# cv iteration with the most excluded plots
nndm_deciduous <- cv_results$deciduous$nndm
id_plot_deciduous <- which.max(sapply(nndm_deciduous$indx_exclude, length))

deciduous_plot <- plots_lon_deciduous
deciduous_plot$set <- ""
deciduous_plot$set[nndm_deciduous$indx_train[[id_plot_deciduous]]] <- 'train'
deciduous_plot$set[nndm_deciduous$indx_exclude[[id_plot_deciduous]]] <- 'exclude'
deciduous_plot$set[nndm_deciduous$indx_test[[id_plot_deciduous]]] <- 'test'
deciduous_plot <- deciduous_plot[order(deciduous_plot$set),]

p_deciduous <- ggplot() +
  geom_sf(data = ext_loff, fill = 'grey', alpha = 0.1) +
  geom_sf(data = deciduous_plot, aes(col = set)) +
  scale_color_brewer(palette = 'Dark2') +
  theme_bw() +
  ggtitle('NNDM visualization - deciduous',
          subtitle = 'CV iteration with most excluded plots')

# show the two plots side by side
cowplot::plot_grid(p_coniferous, p_deciduous, ncol = 2)

# fit models and estimate their performance (leaf-on only)
data_df_list <- list(
  all = plots_lon_all_df,
  deciduous = plots_lon_deciduous_df,
  coniferous = plots_lon_coniferous_df
)

# function to fit RF models with different CV methods
fit_cv_models <- function(data_df, nndm_result, predictor_start_col = 19) {
  
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



# 04: model training with the base ALS metrics
#-------------------------------------------------------------------------------

# Two predictor sets are trained:
#   - 04:  base ALS metrics, without the structural complexity metrics
#   - 04b: with the structural complexity metrics added
#          (box_dimension, vci, rumple, enl_richness, enl_shannon, enl_simpson)

# predictors start at 'mean' (first calc_metrics() column); fail fast if the
# column layout of plot_metrics_*_new.gpkg ever changes
predictor_start_col <- 19
stopifnot(names(plot_metrics_lon)[predictor_start_col] == 'mean')

# structural complexity metrics, excluded for the part 04 (base) predictor set
structural_cols <- c('box_dimension', 'vci', 'rumple',
                     'enl_richness', 'enl_shannon', 'enl_simpson')

# response variables: AGB and stem density
response_vars <- list(
  list(name = 'agb_ha',       col = 'agb_ha'),
  list(name = 'tree_density', col = 'tree_density')
)

# fail fast if a response column is missing or incomplete
for (resp in response_vars) {
  stopifnot(resp$col %in% names(plot_metrics_lon))
  if (any(is.na(sf::st_drop_geometry(plot_metrics_lon)[[resp$col]]))) {
    warning('NA values in response variable: ', resp$col)
  }
}

# the 6 datasets, built with the full predictor set (including the structural
# complexity metrics); part 04's base set drops the structural columns below.
# pixel_centroids (predpoints) are reused for the NNDM: the plot/raster
# locations are the same, only the metric calculation method differs
training_data_struct_comp <- list(
  lon_all = list(
    data = plot_metrics_lon,
    predictors = sf::st_drop_geometry(
      plot_metrics_lon[, predictor_start_col:ncol(plot_metrics_lon)]
      ),
    predpoints = pixel_centroids
  ),
  lon_deciduous = list(
    data = plots_lon_deciduous,
    predictors = sf::st_drop_geometry(
      plots_lon_deciduous[, predictor_start_col:ncol(plots_lon_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  lon_coniferous = list(
    data = plots_lon_coniferous,
    predictors = sf::st_drop_geometry(
      plots_lon_coniferous[, predictor_start_col:ncol(plots_lon_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  ),
  loff_all = list(
    data = plot_metrics_loff,
    predictors = sf::st_drop_geometry(
      plot_metrics_loff[, predictor_start_col:ncol(plot_metrics_loff)]
      ),
    predpoints = pixel_centroids
  ),
  loff_deciduous = list(
    data = plots_loff_deciduous,
    predictors = sf::st_drop_geometry(
      plots_loff_deciduous[, predictor_start_col:ncol(plots_loff_deciduous)]
      ),
    predpoints = pixel_centroids_deciduous
  ),
  loff_coniferous = list(
    data = plots_loff_coniferous,
    predictors = sf::st_drop_geometry(
      plots_loff_coniferous[, predictor_start_col:ncol(plots_loff_coniferous)]
      ),
    predpoints = pixel_centroids_coniferous
  )
)

# dataset names for iteration
dataset_names <- names(training_data_struct_comp)

# NNDM depends only on the plot geometries and predpoints, not on which
# predictor columns are used, so it is computed once here and reused for both
# the base (part 04) and the structure-inclusive (part 04b) predictor set
message('Initializing NNDM for all datasets...')
nndm_results <- lapply(names(training_data_struct_comp), function(name) {
  message(paste0('  Computing NNDM for: ', name))
  CAST::nndm(
    tpoints = training_data_struct_comp[[name]]$data,
    predpoints = training_data_struct_comp[[name]]$predpoints
  )
})
names(nndm_results) <- names(training_data_struct_comp)

train_controls <- lapply(names(nndm_results), function(name) {
  caret::trainControl(
    method = 'cv',
    index = nndm_results[[name]]$indx_train,
    indexOut = nndm_results[[name]]$indx_test,
    savePredictions = T,
    allowParallel = T
  )
})
names(train_controls) <- names(training_data_struct_comp)

# tuning parameters for the random forest model
tgrid <- expand.grid(
  mtry          = 1:3,
  splitrule     = 'variance',
  min.node.size = 5
)

# part 04 (base) predictor set: drop the structural complexity metrics
training_data_base <- lapply(training_data_struct_comp, function(d) {
  d$predictors <- d$predictors[, setdiff(names(d$predictors), structural_cols), drop = FALSE]
  d
})

# generic training loop, reused for both predictor sets
train_ffs_models <- function(training_data_set, prefix, label) {

  models <- list()

  models_to_train <- list()
  for (resp in response_vars) {
    for (dataset_name in dataset_names) {

      model_name <- paste0(prefix, resp$name, '_', dataset_name)
      file_name  <- paste0(prefix, resp$name, '_', dataset_name, '.RDS')
      file_path  <- file.path(processed_data_dir, 'models', file_name)

      if (file.exists(file_path)) {
        message(paste0('Loading existing model: ', file_name))
        models[[model_name]] <- readRDS(file_path)
      } else {
        models_to_train[[model_name]] <- list(
          resp = resp, dataset_name = dataset_name, file_name = file_name
        )
      }
    }
  }

  message(paste0('\n--- ', label, ' models to train: ', length(models_to_train), ' ---\n'))

  if (length(models_to_train) > 0) {

    n_cores <- parallel::detectCores() - 2
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    for (model_name in names(models_to_train)) {

      model_info   <- models_to_train[[model_name]]
      resp         <- model_info$resp
      dataset_name <- model_info$dataset_name
      file_name    <- model_info$file_name

      message(paste0('\n--- Training model: ', model_name, ' ---\n'))

      predictors    <- training_data_set[[dataset_name]]$predictors
      response_data <- training_data_set[[dataset_name]]$data
      ctrl          <- train_controls[[dataset_name]]
      response      <- sf::st_drop_geometry(response_data[[resp$col]])

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

      models[[model_name]] <- ffs_model
      saveRDS(ffs_model, file.path(processed_data_dir, 'models', file_name))
      message(paste0('Model saved: ', file_name))
    }

    parallel::stopCluster(cl)

  } else {
    message(paste0('\n--- All ', label, ' models already exist, no training needed ---\n'))
  }

  models
}

ffs_models_base <- train_ffs_models(
  training_data_base, 'ffs_rf_base_', 'base ALS'
  )

selected_features_base <- data.frame(
  model      = names(ffs_models_base),
  n_selected = sapply(ffs_models_base, function(m) length(m$selectedvars)),
  variables  = sapply(ffs_models_base, function(m) paste(m$selectedvars, collapse = ', ')),
  row.names  = NULL
)

selected_features_base %>%
  knitr::kable(caption = 'Variables selected by FFS per model (base ALS metrics)')

write.csv(
  selected_features_base,
  file.path(output_dir, 'ffs_selected_variables_base.csv'),
  row.names = F
)



# 04b: model training with structural complexity metrics added
#-------------------------------------------------------------------------------

# box_dimension, vci, rumple, enl_richness, enl_shannon, enl_simpson added to
# the base ALS metrics from part 04. Same data / NNDM folds / tgrid / seed as
# part 04, so this is directly comparable to it.

ffs_models_struct_comp <- train_ffs_models(
  training_data_struct_comp, 'ffs_rf_struct_comp_', 'ALS + structural complexity'
  )

selected_features_struct_comp <- data.frame(
  model      = names(ffs_models_struct_comp),
  n_selected = sapply(ffs_models_struct_comp, function(m) length(m$selectedvars)),
  variables  = sapply(ffs_models_struct_comp, function(m) paste(m$selectedvars, collapse = ', ')),
  row.names  = NULL
)

selected_features_struct_comp %>%
  knitr::kable(caption = 'Variables selected by FFS per model (ALS + structural complexity)')

write.csv(
  selected_features_struct_comp,
  file.path(output_dir, 'ffs_selected_variables_struct_comp.csv'),
  row.names = F
)

# compare the base ALS metrics (part 04) against the ALS + structural complexity
# metrics (part 04b), using the global NNDM LOO CV performance of each model
gv_table <- function(models, predictor_set, prefix) {
  do.call(rbind, lapply(names(models), function(mn) {
    gv <- CAST::global_validation(models[[mn]])
    data.frame(
      key           = sub(paste0('^', prefix), '', mn),
      predictor_set = predictor_set,
      RMSE          = gv[['RMSE']],
      R2            = gv[['Rsquared']],
      row.names     = NULL
    )
  }))
}

model_comparison <- dplyr::bind_rows(
  gv_table(ffs_models_base,        'base',        'ffs_rf_base_'),
  gv_table(ffs_models_struct_comp, 'struct_comp', 'ffs_rf_struct_comp_')
) %>%
  tidyr::pivot_wider(names_from = predictor_set, values_from = c(RMSE, R2)) %>%
  dplyr::mutate(
    delta_RMSE = RMSE_struct_comp - RMSE_base,
    delta_R2   = R2_struct_comp - R2_base
  )

model_comparison %>%
  knitr::kable(digits = 2,
               caption = 'Base ALS vs. ALS + structural complexity (NNDM LOO CV, global validation)')

write.csv(
  model_comparison,
  file.path(output_dir, 'model_comparison_base_vs_struct_comp.csv'),
  row.names = F
)



# 05: validation
#-------------------------------------------------------------------------------

# print summary of all trained models (base ALS and ALS + structural complexity)
all_trained_models <- c(ffs_models_base, ffs_models_struct_comp)
message('\n--- Summary of all trained models ---\n')
for (model_name in names(all_trained_models)) {
  message(paste0('\n', model_name, ':'))
  print(CAST::global_validation(all_trained_models[[model_name]]))
}

# extract cross-validation predictions from a set of trained models
# (cached to disk; file_tag keeps the base and struct_comp predictions separate)
extract_cv_predictions <- function(models, model_prefix, file_tag) {

  preds <- list()

  for (resp in response_vars) {
    for (dataset_name in dataset_names) {

      model_name <- paste0(model_prefix, resp$name, '_', dataset_name)
      pred_name  <- paste0('pred_', resp$name, '_', dataset_name)
      file_name  <- paste0('pred_obsv_', file_tag, resp$name, '_', dataset_name, '.gpkg')
      file_path  <- file.path(processed_data_dir, 'predictions', file_name)

      if (file.exists(file_path)) {

        message(paste0('Loading existing CV predictions: ', file_name))
        preds[[pred_name]] <- sf::st_read(file_path, quiet = T)

      } else {

        if (!model_name %in% names(models)) {
          message(paste0('Model not found: ', model_name, ' - skipping'))
          next
        }

        model   <- models[[model_name]]
        cv_pred <- model$pred

        # link CV predictions to the original plot geometries
        # (base and struct_comp share the same plot data, only the predictor
        # columns differ, so either training data list can be used here)
        plot_data  <- training_data_struct_comp[[dataset_name]]$data
        cv_pred_sf <- plot_data[cv_pred$rowIndex, ] %>%
          dplyr::select(key, kspnr, abt) %>%
          dplyr::mutate(pred = cv_pred$pred, obs = cv_pred$obs)

        preds[[pred_name]] <- cv_pred_sf
        sf::st_write(cv_pred_sf, file_path, delete_dsn = T)
        message(paste0('CV predictions saved: ', file_name))
      }
    }
  }

  preds
}

# base ALS predictions
cv_predictions_base <- extract_cv_predictions(
  ffs_models_base, 'ffs_rf_base_', 'base_'
)

# ALS + structural complexity predictions
cv_predictions_struct_comp <- extract_cv_predictions(
  ffs_models_struct_comp, 'ffs_rf_struct_comp_', 'struct_comp_'
)

# compute validation metrics for a set of CV predictions
compute_cv_metrics <- function(cv_preds, predictor_set) {
  do.call(rbind, lapply(names(cv_preds), function(nm) {
    d  <- sf::st_drop_geometry(cv_preds[[nm]])
    se <- (d$pred - d$obs)^2
    ae <- abs(d$pred - d$obs)
    data.frame(
      predictor_set = predictor_set,
      model_name    = nm,
      n             = nrow(d),
      RMSE          = sqrt(mean(se, na.rm = T)),
      R2            = 1 - sum(se, na.rm = T) / sum((d$obs - mean(d$obs, na.rm = T))^2, na.rm = T),
      rel_RMSE      = sqrt(mean(se, na.rm = T)) / mean(d$obs, na.rm = T) * 100,
      MAE           = mean(ae, na.rm = T),
      bias          = mean(d$pred - d$obs, na.rm = T),
      rel_bias      = mean(d$pred - d$obs, na.rm = T) / mean(d$obs, na.rm = T) * 100,
      row.names     = NULL
    )
  }))
}

# combine metrics of both predictor sets
model_metrics <- dplyr::bind_rows(
  compute_cv_metrics(cv_predictions_base,        'base'),
  compute_cv_metrics(cv_predictions_struct_comp, 'struct_comp')
) %>%
  dplyr::mutate(
    response_var = gsub('pred_(.+)_(lon|loff)_.*', '\\1', model_name),
    leaf_condition = dplyr::case_when(
      grepl('_lon_', model_name) ~ 'lon',
      grepl('_loff_', model_name) ~ 'loff',
      T ~ NA_character_
    ),
    leaf_type = dplyr::case_when(
      grepl('_deciduous$', model_name) ~ 'deciduous',
      grepl('_coniferous$', model_name) ~ 'coniferous',
      grepl('_all$', model_name) ~ 'all',
      T ~ NA_character_
    )
  )

# assemble the final metrics table (base vs. struct_comp side by side)
metrics_table <- model_metrics %>%
  dplyr::select(response_var, leaf_condition, leaf_type, predictor_set,
                n, R2, RMSE, rel_RMSE, MAE, bias, rel_bias) %>%
  dplyr::arrange(response_var, leaf_condition, leaf_type, predictor_set)

# display the table
metrics_table %>%
  knitr::kable(digits = 2)

# store the table
write.csv(
  metrics_table,
  file.path(output_dir, 'model_validation_metrics.csv'),
  row.names = F
)

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

# define units for forest inventory attributes as plotmath snippets
# (parsed into expressions for superscripts in the plot titles)
forest_inv_units <- c(
  'agb_ha'        = '"Mg" ~ "ha"^-1',
  'total_vol_ha'  = '"m"^3 ~ "ha"^-1',
  'merch_vol_ha'  = '"m"^3 ~ "ha"^-1',
  'basal_area_ha' = '"m"^2 ~ "ha"^-1',
  'dg'            = '"cm"',
  'tree_density'  = '"n" ~ "ha"^-1'
)

# plot: leaf-on vs. leaf-off performance, for both predictor sets
# (rows = metric, columns = predictor set, bars grouped by leaf condition)
season_cols <- c('leaf-on' = 'gray40', 'leaf-off' = 'gray70')

for (resp in unique(model_metrics$response_var)) {

  season_plot_data <- model_metrics %>%
    dplyr::filter(response_var == resp) %>%
    dplyr::select(predictor_set, leaf_condition, leaf_type, rel_RMSE, MAE, R2) %>%
    tidyr::pivot_longer(cols = c(rel_RMSE, MAE, R2),
                        names_to = 'metric', values_to = 'value') %>%
    dplyr::mutate(
      leaf_condition = factor(
        ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
        levels = c('leaf-on', 'leaf-off')
      ),
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
      predictor_set = factor(
        predictor_set,
        levels = c('base', 'struct_comp'),
        labels = c('base ALS', 'ALS + structural complexity')
      ),
      metric = factor(metric, levels = c('rel_RMSE', 'MAE', 'R2'),
                      labels = c('relative RMSE [%]', 'MAE', 'R²'))
    )

  if (nrow(season_plot_data) == 0) next

  p_season <- ggplot(season_plot_data,
                     aes(x = leaf_type, y = value, fill = leaf_condition)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf('%.2f', value)),
              position = position_dodge(width = 0.8),
              vjust = -0.4, size = 2.5) +
    facet_grid(metric ~ predictor_set, scales = 'free_y', switch = 'y') +
    scale_fill_manual(values = season_cols, name = '') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = bquote('Leaf-on vs. leaf-off:' ~ .(forest_inv_names[resp]) ~
                        '[' * .(parse(text = forest_inv_units[resp])[[1]]) * ']'),
         x = '', y = '') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          strip.background = element_rect(fill = 'lightgrey'),
          strip.text = element_text(face = 'bold'),
          strip.placement = 'outside')

  print(p_season)

  ggplot2::ggsave(
    filename = file.path(output_dir,
                         paste0('leafon_vs_leafoff_', resp, '.pdf')),
    plot = p_season, width = 9, height = 8
  )
}



# predicted vs. observed plots

# combine CV predictions of both predictor sets into one data frame
cv_predictions_df <- dplyr::bind_rows(
  do.call(rbind, lapply(names(cv_predictions_base), function(nm)
    cv_predictions_base[[nm]] %>% sf::st_drop_geometry() %>%
      dplyr::mutate(model_name = nm, predictor_set = 'base'))),
  do.call(rbind, lapply(names(cv_predictions_struct_comp), function(nm)
    cv_predictions_struct_comp[[nm]] %>% sf::st_drop_geometry() %>%
      dplyr::mutate(model_name = nm, predictor_set = 'struct_comp')))
) %>%
  dplyr::mutate(
    response_var = gsub('pred_(.+)_(lon|loff)_.*', '\\1', model_name),
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

# readable labels for the predictor sets (used as plot subtitle)
pset_labels <- c(
  'base'        = 'base ALS metrics',
  'struct_comp' = 'ALS + structural complexity metrics'
)

# one predicted vs. observed figure per response variable and predictor set
for (resp in unique(cv_predictions_df$response_var)) {
  for (pset in c('base', 'struct_comp')) {

    plot_data_resp <- cv_predictions_df %>%
      dplyr::filter(response_var == resp, predictor_set == pset) %>%
      dplyr::mutate(
        leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
        leaf_condition = factor(leaf_condition, levels = c('leaf-on', 'leaf-off'))
      )

    if (nrow(plot_data_resp) == 0) next

    forest_inv_label <- forest_inv_names[resp]
    axis_max <- max(c(plot_data_resp$obs, plot_data_resp$pred), na.rm = T)
    axis_min <- min(c(plot_data_resp$obs, plot_data_resp$pred), na.rm = T)

    # per-panel relative RMSE and R2 for annotation
    panel_metrics <- plot_data_resp %>%
      dplyr::group_by(leaf_condition, leaf_type) %>%
      dplyr::summarise(
        rel_RMSE = sqrt(mean((pred - obs)^2, na.rm = T)) / mean(obs, na.rm = T) * 100,
        R2 = 1 - sum((pred - obs)^2, na.rm = T) / sum((obs - mean(obs, na.rm = T))^2, na.rm = T),
        .groups = 'drop'
      ) %>%
      dplyr::mutate(label = paste0('rRMSE: ', round(rel_RMSE, 1), '%\n',
                                   'R² = ', round(R2, 2)))

    p <- ggplot(plot_data_resp, aes(x = pred, y = obs)) +
      geom_point(alpha = 0.7, size = 2) +
      geom_abline(slope = 1, intercept = 0, linewidth = 1,
                  color = 'red', linetype = 'dashed') +
      geom_smooth(method = 'lm', se = F, linewidth = 1, color = 'black') +
      geom_text(data = panel_metrics,
                aes(x = axis_min + (axis_max - axis_min) * 0.05,
                    y = axis_max - (axis_max - axis_min) * 0.05,
                    label = label),
                hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
      facet_grid(leaf_condition ~ leaf_type) +
      coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
      labs(title = bquote('Predicted vs. Observed:' ~ .(forest_inv_label) ~
                          '[' * .(parse(text = forest_inv_units[resp])[[1]]) * ']'),
           subtitle = pset_labels[pset],
           x = 'Predicted',
           y = 'Observed') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
            plot.subtitle = element_text(hjust = 0.5),
            strip.background = element_rect(fill = 'lightgrey'),
            strip.text = element_text(face = 'bold'),
            panel.grid = element_blank())

    print(p)

    ggplot2::ggsave(
      filename = file.path(output_dir,
        paste0('pred_obs_', resp, '_', gsub('[^A-Za-z]', '', pset), '.pdf')),
      plot = p, width = 10, height = 8
    )
  }
}

# plot: how each plot's prediction moves from leaf-on to leaf-off
# arrows pointing towards the 1:1 line mean the leaf-off model
# predicts that plot better, arrows pointing away mean it predicts it worse
for (resp in unique(cv_predictions_df$response_var)) {

  arrow_data <- cv_predictions_df %>%
    dplyr::filter(response_var == resp) %>%
    dplyr::select(predictor_set, leaf_type, kspnr, obs, leaf_condition, pred) %>%
    tidyr::pivot_wider(names_from = leaf_condition, values_from = pred) %>%
    dplyr::filter(!is.na(`leaf-on`), !is.na(`leaf-off`)) %>%
    dplyr::mutate(
      pred_lon  = `leaf-on`,
      pred_loff = `leaf-off`,
      direction = ifelse(abs(pred_loff - obs) < abs(pred_lon - obs),
                         'closer to 1:1', 'further from 1:1'),
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
      predictor_set = factor(
        predictor_set,
        levels = c('base', 'struct_comp'),
        labels = c('base ALS', 'ALS + structural complexity')
      )
    )

  if (nrow(arrow_data) == 0) next

  axis_min <- min(c(arrow_data$pred_lon, arrow_data$pred_loff, arrow_data$obs), na.rm = T)
  axis_max <- max(c(arrow_data$pred_lon, arrow_data$pred_loff, arrow_data$obs), na.rm = T)

  # share of plots that improve, per panel (for annotation)
  # the two-sided binomial test asks whether that share differs from 50%,
  # i.e. whether either season is systematically better on a per-plot basis
  arrow_summary <- arrow_data %>%
    dplyr::group_by(predictor_set, leaf_type) %>%
    dplyr::summarise(
      n_plots      = dplyr::n(),
      n_closer     = sum(direction == 'closer to 1:1'),
      share_closer = n_closer / n_plots * 100,
      p_value      = stats::binom.test(n_closer, n_plots, p = 0.5)$p.value,
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      label = paste0(round(share_closer), '% closer (n = ', n_plots, ')\n',
                     'p = ', format.pval(p_value, digits = 2, eps = 0.001))
    )

  p_arrows <- ggplot(arrow_data, aes(x = pred_lon, y = obs)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8,
                color = 'black', linetype = 'dashed') +
    geom_segment(aes(xend = pred_loff, yend = obs, colour = direction),
                 arrow = ggplot2::arrow(length = grid::unit(0.10, 'cm'),
                                        type = 'closed'),
                 linewidth = 0.4, alpha = 0.85) +
    geom_point(size = 0.7, colour = 'grey30', alpha = 0.6) +
    geom_text(data = arrow_summary,
              aes(x = axis_min + (axis_max - axis_min) * 0.05,
                  y = axis_max - (axis_max - axis_min) * 0.05,
                  label = label),
              hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
    facet_grid(predictor_set ~ leaf_type) +
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    scale_colour_manual(
      values = c('closer to 1:1' = '#009E73', 'further from 1:1' = '#D55E00'),
      name = 'leaf-off prediction'
    ) +
    labs(title = bquote('Shift from leaf-on to leaf-off:' ~ .(forest_inv_names[resp]) ~
                        '[' * .(parse(text = forest_inv_units[resp])[[1]]) * ']'),
         subtitle = 'arrow start = leaf-on prediction, arrow head = leaf-off prediction',
         x = 'Predicted',
         y = 'Observed') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.background = element_rect(fill = 'lightgrey'),
          strip.text = element_text(face = 'bold'),
          panel.grid = element_blank(),
          legend.position = 'bottom')

  print(p_arrows)

  ggplot2::ggsave(
    filename = file.path(output_dir,
                         paste0('pred_shift_lon_to_loff_', resp, '.pdf')),
    plot = p_arrows, width = 10, height = 8
  )
}

# plot: how each plot's prediction moves when the structural complexity
# metrics are added to the base ALS metrics
# same principle as the plot above, but the arrow now starts at the base ALS
# prediction and points to the ALS + structural complexity prediction. 
for (resp in unique(cv_predictions_df$response_var)) {

  arrow_data_struct <- cv_predictions_df %>%
    dplyr::filter(response_var == resp, leaf_condition == 'leaf-off') %>%
    dplyr::select(leaf_type, kspnr, obs, predictor_set, pred) %>%
    tidyr::pivot_wider(names_from = predictor_set, values_from = pred) %>%
    dplyr::filter(!is.na(base), !is.na(struct_comp)) %>%
    dplyr::mutate(
      pred_base   = base,
      pred_struct = struct_comp,
      direction = ifelse(abs(pred_struct - obs) < abs(pred_base - obs),
                         'closer to 1:1', 'further from 1:1'),
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous'))
    )

  if (nrow(arrow_data_struct) == 0) next

  axis_min <- min(c(arrow_data_struct$pred_base, arrow_data_struct$pred_struct,
                    arrow_data_struct$obs), na.rm = T)
  axis_max <- max(c(arrow_data_struct$pred_base, arrow_data_struct$pred_struct,
                    arrow_data_struct$obs), na.rm = T)

  # share of plots that improve, per panel (for annotation)
  # the two-sided binomial test asks whether that share differs from 50%,
  # i.e. whether adding the structural metrics systematically helps on a
  # per-plot basis
  # p_adj applies the Holm correction over the three leaf types of this figure
  arrow_summary_struct <- arrow_data_struct %>%
    dplyr::group_by(leaf_type) %>%
    dplyr::summarise(
      n_plots      = dplyr::n(),
      n_closer     = sum(direction == 'closer to 1:1'),
      share_closer = n_closer / n_plots * 100,
      p_value      = stats::binom.test(n_closer, n_plots, p = 0.5)$p.value,
      .groups = 'drop'
    ) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = 'holm'),
      label = paste0(round(share_closer), '% closer (n = ', n_plots, ')\n',
                     'p = ', format.pval(p_value, digits = 2, eps = 0.001),
                     '\nHolm: ', format.pval(p_adj, digits = 2, eps = 0.001))
    )

  p_arrows_struct <- ggplot(arrow_data_struct, aes(x = pred_base, y = obs)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8,
                color = 'black', linetype = 'dashed') +
    geom_segment(aes(xend = pred_struct, yend = obs, colour = direction),
                 arrow = ggplot2::arrow(length = grid::unit(0.10, 'cm'),
                                        type = 'closed'),
                 linewidth = 0.4, alpha = 0.85) +
    geom_point(size = 0.7, colour = 'grey30', alpha = 0.6) +
    geom_text(data = arrow_summary_struct,
              aes(x = axis_min + (axis_max - axis_min) * 0.05,
                  y = axis_max - (axis_max - axis_min) * 0.05,
                  label = label),
              hjust = 0, vjust = 1, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ leaf_type, nrow = 1) +
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    scale_colour_manual(
      values = c('closer to 1:1' = '#009E73', 'further from 1:1' = '#D55E00'),
      name = 'with structural complexity'
    ) +
    labs(title = bquote('Effect of adding structural complexity (leaf-off):' ~
                        .(forest_inv_names[resp]) ~
                        '[' * .(parse(text = forest_inv_units[resp])[[1]]) * ']'),
         subtitle = 'arrow start = base ALS prediction, arrow head = ALS + structural complexity prediction',
         x = 'Predicted',
         y = 'Observed') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          plot.subtitle = element_text(hjust = 0.5, size = 8),
          strip.background = element_rect(fill = 'lightgrey'),
          strip.text = element_text(face = 'bold'),
          panel.grid = element_blank(),
          legend.position = 'bottom')

  print(p_arrows_struct)

  ggplot2::ggsave(
    filename = file.path(output_dir,
                         paste0('pred_shift_base_to_struct_comp_', resp, '.pdf')),
    plot = p_arrows_struct, width = 10, height = 4.5
  )
}



# 06: model transferability
#-------------------------------------------------------------------------------

# leaf types
leaf_types <- c('all', 'deciduous', 'coniferous')

# leaf-on / leaf-off data frames (predictors + response) per leaf type
transfer_data <- list(
  all        = list(lon_df = plots_lon_all_df, 
                    loff_df = plots_loff_all_df),
  deciduous  = list(lon_df = plots_lon_deciduous_df,
                    loff_df = plots_loff_deciduous_df),
  coniferous = list(lon_df = plots_lon_coniferous_df,
                    loff_df = plots_loff_coniferous_df)
)

# store transferability predictions for all response variables and leaf types
transfer_predictions <- list()

for (resp in response_vars) {
  for (lt in leaf_types) {

    lon_model_name  <- paste0('ffs_rf_', resp$name, '_lon_', lt)
    loff_model_name <- paste0('ffs_rf_', resp$name, '_loff_', lt)

    # check that both models exist
    if (!lon_model_name %in% names(ffs_models) | !loff_model_name %in% names(ffs_models)) {
      message(paste0('Skipping ', resp$name, ' / ', lt, ': model(s) not found'))
      next
    }

    lon_model  <- ffs_models[[lon_model_name]]
    loff_model <- ffs_models[[loff_model_name]]

    # predictor data frames for this leaf type
    lon_df  <- transfer_data[[lt]]$lon_df
    loff_df <- transfer_data[[lt]]$loff_df
    preds_lon  <- lon_df[, predictor_start_col:ncol(lon_df)]
    preds_loff <- loff_df[, predictor_start_col:ncol(loff_df)]

    # 1) + 2) native performance from NNDM LOO CV predictions
    cv_lon  <- cv_predictions[[paste0('pred_', resp$name, '_lon_', lt)]]  %>% sf::st_drop_geometry()
    cv_loff <- cv_predictions[[paste0('pred_', resp$name, '_loff_', lt)]] %>% sf::st_drop_geometry()

    # 3) leaf-on model applied to leaf-off data (transfer)
    obs_loff         <- loff_df[[resp$col]]
    pred_lon_on_loff <- stats::predict(lon_model, preds_loff)

    # 4) leaf-off model applied to leaf-on data (transfer)
    obs_lon          <- lon_df[[resp$col]]
    pred_loff_on_lon <- stats::predict(loff_model, preds_lon)

    # combine the four scenarios into one data frame
    transfer_df <- rbind(
      data.frame(response = resp$name, leaf_type = lt,
                 scenario = 'lon model -> lon data',   type = 'native (CV)',
                 obs = cv_lon$obs,  pred = cv_lon$pred),
      data.frame(response = resp$name, leaf_type = lt,
                 scenario = 'loff model -> loff data', type = 'native (CV)',
                 obs = cv_loff$obs, pred = cv_loff$pred),
      data.frame(response = resp$name, leaf_type = lt,
                 scenario = 'lon model -> loff data',  type = 'transfer',
                 obs = obs_loff,    pred = pred_lon_on_loff),
      data.frame(response = resp$name, leaf_type = lt,
                 scenario = 'loff model -> lon data',  type = 'transfer',
                 obs = obs_lon,     pred = pred_loff_on_lon)
    )

    transfer_predictions[[paste0(resp$name, '_', lt)]] <- transfer_df
  }
}

# combine all prediction scenarios
scenario_predictions <- dplyr::bind_rows(transfer_predictions)

# calculate error metrics for all combinations
transfer_metrics <- scenario_predictions %>%
  dplyr::group_by(response, leaf_type, scenario, type) %>%
  dplyr::summarise(
    n = dplyr::n(),
    RMSE = sqrt(mean((pred - obs)^2, na.rm = T)),
    mean_obs = mean(obs, na.rm = T),
    rel_RMSE = sqrt(mean((pred - obs)^2, na.rm = T)) / mean(obs, na.rm = T) * 100,
    R2 = 1 - sum((pred - obs)^2, na.rm = T) / sum((obs - mean(obs, na.rm = T))^2, na.rm = T),
    MAE = mean(abs(pred - obs), na.rm = T),
    bias = mean(pred - obs, na.rm = T),
    rel_bias = mean(pred - obs, na.rm = T) / mean(obs, na.rm = T) * 100,
    .groups = 'drop'
  )

# order rows: per leaf type, native scenarios first, then transfer scenarios
transfer_metrics <- transfer_metrics %>%
  dplyr::arrange(
    response,
    factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
    factor(scenario, levels = c('lon model -> lon data', 'loff model -> loff data',
                                'lon model -> loff data', 'loff model -> lon data'))
  )

# display full results table
transfer_metrics %>%
  dplyr::select(response, leaf_type, scenario, type,
                n, RMSE, rel_RMSE, R2, MAE, bias, rel_bias) %>%
  knitr::kable(digits = 2,
               caption = 'Transferability: 4 scenarios per leaf type')

# store the metrics table
write.csv(
  transfer_metrics,
  file.path(output_dir, 'model_transferability_metrics.csv'),
  row.names = F
)

# prepare plotting data
transfer_plot_data <- transfer_metrics %>%
  dplyr::mutate(
    leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
    scenario_label = factor(scenario,
      levels = c('lon model -> lon data', 'loff model -> loff data',
                 'lon model -> loff data', 'loff model -> lon data'),
      labels = c('leaf-on -> leaf-on', 'leaf-off -> leaf-off',
                 'leaf-on -> leaf-off', 'leaf-off -> leaf-on'))
  )

scenario_cols <- c(
  'leaf-on -> leaf-on'   = '#009E73',
  'leaf-off -> leaf-off' = '#E69F00',
  'leaf-on -> leaf-off'  = '#80CEB9',
  'leaf-off -> leaf-on'  = '#F2CF80'
)

# plot: relative RMSE per leaf type and scenario
p_transfer_rmse <- ggplot(transfer_plot_data,
       aes(x = leaf_type, y = rel_RMSE, fill = scenario_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  geom_text(aes(label = round(rel_RMSE, 1)),
            position = position_dodge(width = 0.85), vjust = -0.5, size = 2.5) +
  scale_fill_manual(values = scenario_cols, name = '') +
  labs(x = '', y = 'Relative RMSE (%)') +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'bottom')

print(p_transfer_rmse)

# plot: relative bias per leaf type and scenario
p_transfer_bias <- ggplot(transfer_plot_data,
       aes(x = leaf_type, y = rel_bias, fill = scenario_label)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.75) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey40') +
  geom_text(aes(label = round(rel_bias, 1),
                vjust = ifelse(rel_bias >= 0, -0.5, 1.5)),
            position = position_dodge(width = 0.85), size = 2.5) +
  scale_fill_manual(values = scenario_cols, name = '') +
  labs(x = '', y = 'Relative Bias (%)') +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = 'bottom')

print(p_transfer_bias)

# save the two barplots
ggplot2::ggsave(
  filename = file.path(output_dir, 'transferability_rel_rmse.pdf'),
  plot = p_transfer_rmse,
  width = 9,
  height = 6
)
ggplot2::ggsave(
  filename = file.path(output_dir, 'transferability_rel_bias.pdf'),
  plot = p_transfer_bias,
  width = 9,
  height = 6
)

# plot: predicted vs observed for all scenarios,
# faceted leaf type (rows) x scenario (columns)
for (resp in unique(scenario_predictions$response)) {

  plot_data_transfer <- scenario_predictions %>%
    dplyr::filter(response == resp) %>%
    dplyr::mutate(
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
      scenario_label = factor(scenario,
        levels = c('lon model -> lon data', 'loff model -> loff data',
                   'lon model -> loff data', 'loff model -> lon data'),
        labels = c('leaf-on -> leaf-on', 'leaf-off -> leaf-off',
                   'leaf-on -> leaf-off', 'leaf-off -> leaf-on'))
    )

  if (nrow(plot_data_transfer) == 0) next

  forest_inv_label <- forest_inv_names[resp]
  axis_max <- max(c(plot_data_transfer$obs, plot_data_transfer$pred), na.rm = T)
  axis_min <- min(c(plot_data_transfer$obs, plot_data_transfer$pred), na.rm = T)

  # per-panel metrics for annotation
  metrics_for_plot <- transfer_metrics %>%
    dplyr::filter(response == resp) %>%
    dplyr::mutate(
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
      scenario_label = factor(scenario,
        levels = c('lon model -> lon data', 'loff model -> loff data',
                   'lon model -> loff data', 'loff model -> lon data'),
        labels = c('leaf-on -> leaf-on', 'leaf-off -> leaf-off',
                   'leaf-on -> leaf-off', 'leaf-off -> leaf-on')),
      label = paste0('rRMSE: ', round(rel_RMSE, 1), '%\n', 'R² = ', round(R2, 2))
    )

  p <- ggplot(plot_data_transfer, aes(x = pred, y = obs)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8,
                color = 'red', linetype = 'dashed') +
    geom_smooth(method = 'lm', se = F, linewidth = 0.8, color = 'black') +
    geom_text(data = metrics_for_plot,
              aes(x = axis_min + (axis_max - axis_min) * 0.05,
                  y = axis_max - (axis_max - axis_min) * 0.05,
                  label = label),
              hjust = 0, vjust = 1, size = 2.5) +
    facet_grid(leaf_type ~ scenario_label) +
    coord_fixed(ratio = 1, xlim = c(axis_min, axis_max), ylim = c(axis_min, axis_max)) +
    labs(title = paste0('Transferability: ', forest_inv_label),
         x = 'Predicted',
         y = 'Observed') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
          strip.background = element_rect(fill = 'lightgrey'),
          strip.text = element_text(face = 'bold'),
          panel.grid = element_blank())

  print(p)

  # save the predicted vs. observed plot for this response variable
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0('transferability_pred_obs_', resp, '.pdf')),
    plot = p,
    width = 10,
    height = 8
  )
}