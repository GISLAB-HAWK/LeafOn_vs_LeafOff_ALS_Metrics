#-------------------------------------------------------------------------------
# Name:         calculate_metrics_new.R
# Description:  Script calculates metrics in forest inventory plots using
#               the point clouds from leaf-on and leaf-off season.
#               The metrics are computed with the lidR package using
#               plot_metrics(). The metric definitions live in
#               src/calc_metrics.R, which is shared with the wall-to-wall
#               (pixel-level) metric calculation so that the predictors are
#               identical for training and prediction.
#               The point clouds are expected to be noise filtered already and
#               to carry a HAG (height above ground) attribute.
#               The calculated metrics are joined with forest inventory data.
# Author:       Florian Franz, Svenja Dobelmann
# Contact:      florian.franz@nw-fva.de
#               svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)

# source the shared metric function (also used for the wall-to-wall metrics)
source('src/calc_metrics.R', local = TRUE)



# 01: processing setup
#-------------------------------------------------------------------------------

# point cloud paths (leaf-on = 2023, leaf-off = 2024)
pc_lon_path  <- file.path(raw_data_dir, 'pc_leafon_2023')
pc_loff_path <- file.path(raw_data_dir, 'pc_leafoff_2024')

# read the point clouds as catalogs
ctg_lon  <- lidR::readLAScatalog(pc_lon_path)
ctg_loff <- lidR::readLAScatalog(pc_loff_path)

# catalog options
# select all attributes so that the HAG extra bytes are loaded
# no filter is set: the point clouds are already noise filtered
lidR::opt_select(ctg_lon)  <- '*'
lidR::opt_select(ctg_loff) <- '*'
lidR::opt_filter(ctg_lon)  <- ''
lidR::opt_filter(ctg_loff) <- ''

# read forest inventory plots (already clipped to the AOI and filtered)
inv_attr_plots <- sf::st_read(
  file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots.gpkg')
)

# sample-circle radius (m) and the resulting circle area (m²)
plot_radius <- 13
plot_area   <- pi * plot_radius^2

# lower bounds of the 1 m layers used for the pr_* / rd_* metrics
# layer 2 covers (2, 3] and layer 10 covers (10, 11],
# yielding pr_2_3, pr_10_11, rd_2_3 and rd_10_11
metric_layers <- c(2, 10)

# parallel processing
n_workers <- 32



# 02: metrics calculation with lidR
#-------------------------------------------------------------------------------

# output paths for the raw plot metrics (cached, computation is expensive)
metrics_lon_file  <- file.path(processed_data_dir, 'metrics', 'plot_metrics_lon_new.gpkg')
metrics_loff_file <- file.path(processed_data_dir, 'metrics', 'plot_metrics_loff_new.gpkg')

if (file.exists(metrics_lon_file) & file.exists(metrics_loff_file)) {

  cat('Loading existing plot metrics from disk...\n')
  pc_lon_metrics  <- sf::st_read(metrics_lon_file,  quiet = T)
  pc_loff_metrics <- sf::st_read(metrics_loff_file, quiet = T)

} else {

  future::plan(future::multisession, workers = n_workers)

  # plot_metrics() extracts the points within a circle of the given radius
  # around each plot centre and returns the plots with the metrics attached
  # X and Y are needed for the box dimension and the rumple index
  metric_func <- ~calc_metrics(
    X,
    Y,
    HAG,
    Classification,
    ReturnNumber,
    NumberOfReturns,
    area = plot_area,
    layers = metric_layers
    )

  cat('Computing leaf-on plot metrics...\n')
  pc_lon_metrics <- lidR::plot_metrics(
    las      = ctg_lon,
    func     = metric_func,
    geometry = inv_attr_plots,
    radius   = plot_radius
  )

  cat('Computing leaf-off plot metrics...\n')
  pc_loff_metrics <- lidR::plot_metrics(
    las      = ctg_loff,
    func     = metric_func,
    geometry = inv_attr_plots,
    radius   = plot_radius
  )

  sf::st_write(pc_lon_metrics,  metrics_lon_file,  delete_dsn = T)
  sf::st_write(pc_loff_metrics, metrics_loff_file, delete_dsn = T)

  future::plan(future::sequential)

}

head(pc_lon_metrics)
head(pc_loff_metrics)



# 03: correlation analysis of metrics with forest inventory attributes
#-------------------------------------------------------------------------------

# response variables (forest inventory attributes)
response_vars <- c('agb_ha', 'tree_density')
# response_vars <- c('agb_ha', 'total_vol_ha', 'merch_vol_ha',
#                    'tree_density', 'basal_area_ha', 'dg')

# leaf-on and leaf-off datasets (without geometry)
metrics_data <- list(
  lon  = sf::st_drop_geometry(pc_lon_metrics),
  loff = sf::st_drop_geometry(pc_loff_metrics)
)

# metric (predictor) columns: everything from the first calc_metrics() column
# ('mean') to the end of the table
predictor_start_col <- 19
stopifnot(names(metrics_data$lon)[predictor_start_col] == 'mean')
metric_cols <- names(metrics_data$lon)[predictor_start_col:ncol(metrics_data$lon)]

# leaf types to analyse ('all' = deciduous and coniferous combined)
leaf_types <- c('all', 'deciduous', 'coniferous')

# compute Pearson correlations between each metric and each response variable,
# separately for leaf-on/leaf-off and for all/deciduous/coniferous plots
correlation_results <- do.call(rbind, lapply(names(metrics_data), function(leaf_condition) {

  df_full <- metrics_data[[leaf_condition]]

  do.call(rbind, lapply(leaf_types, function(lt) {

    # subset by dominant leaf type ('all' keeps every plot)
    df <- if (lt == 'all') df_full else df_full[df_full$dominant_leaf_type == lt, ]

    do.call(rbind, lapply(response_vars, function(resp) {
      data.frame(
        leaf_condition = leaf_condition,
        leaf_type = lt,
        response_var = resp,
        metric = metric_cols,
        correlation = sapply(metric_cols, function(m) {
          stats::cor(df[[m]], df[[resp]], use = 'complete.obs', method = 'pearson')
        }),
        row.names = NULL
      )
    }))
  }))
}))

# display the correlation table (strongest correlations first)
correlation_results %>%
  dplyr::arrange(response_var, leaf_type, leaf_condition, dplyr::desc(abs(correlation))) %>%
  knitr::kable(digits = 2)

# store the correlation table
write.csv(
  correlation_results,
  file.path(output_dir, 'metric_inventory_correlations.csv'),
  row.names = F
)

# plot: leaf-on vs leaf-off correlation with the response variable per metric,
# faceted by leaf type (all / deciduous / coniferous)
for (resp in response_vars) {

  corr_resp <- correlation_results %>%
    dplyr::filter(response_var == resp) %>%
    dplyr::mutate(
      leaf_condition = factor(
        ifelse(leaf_condition == 'lon', 'leaf-on', 'leaf-off'),
        levels = c('leaf-on', 'leaf-off')
      ),
      leaf_type = factor(leaf_type, levels = c('all', 'deciduous', 'coniferous')),
      metric = factor(metric, levels = rev(metric_cols))
    )

  if (nrow(corr_resp) == 0) next

  corr_plot <- ggplot(corr_resp, aes(x = metric, y = correlation, fill = leaf_condition)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, colour = 'grey40') +
    facet_wrap(~ leaf_type, nrow = 1) +
    scale_fill_manual(values = c('leaf-on' = 'gray40', 'leaf-off' = 'gray70'),
                      name = '') +
    labs(x = '', y = paste0('Pearson correlation with ', resp)) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank())

  print(corr_plot)

  ggplot2::ggsave(
    filename = file.path(output_dir,
                         paste0('metric_', resp, '_correlations.pdf')),
    plot = corr_plot,
    width = 10,
    height = 8
  )
}