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