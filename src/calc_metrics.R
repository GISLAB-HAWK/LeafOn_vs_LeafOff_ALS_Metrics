#' Calculation of remote sensing-based forest metrics in ALS-based point clouds.
#'
#' This function calculates several metrics as potential explanatory variables
#' for a growing stock model. It takes the height values of a single point cloud 
#' or a point cloud catalog (LAScatalog from the lidR package) as input to 
#' calculate height, variability, and canopy cover metrics. 
#' The function is intended to be used in the plot_metrics function from the lidR
#' package to calculate the metrics for forest inventory plots stored as an sf object.
#'
#' @param z numeric. The height values of the point cloud.
#' 
#' @details
#' height metrics: mean, standard deviation, minimum, maximum,
#' percentile values (1st, 5th, 10th, 20th, 25th, 30th, 40th, 50th,
#'                    60th, 70th, 75th, 80th, 90th, 95th, 99th);
#' variability metrics: skewness, kurtosis, coefficient of variation;
#' canopy relief ratio (crr) --> https://doi.org/10.1016/j.foreco.2003.09.001;
#' canopy cover metrics: percentage of points above 3m and above mean height
#'
#' --> see different studies, e.g. https://doi.org/10.1139/cjfr-2014-0297
#'
#' @examples
#' # ------------------ #
#' # single point cloud #
#' # ------------------ #
#' 
#' las <- lidR::readLAS(path/to/las_file)
#' plots <- sf::st_read(path/to/shapefile_or_geopackage)
#' 
#' plot_metrics <- lidR::plot_metrics(las, ~calc_metrics(Z, R, B),
#'                                    plots, radius = 13)
#'
#' # ------------------- #
#' # point cloud catalog #
#' # ------------------- #
#' 
#' las_catalog <- lidR::readLAScatalog(path/to/las_file)
#' plots <- sf::st_read(path/to/shapefile_or_geopackage)
#' 
#' plot_metrics <- lidR::plot_metrics(las_catalog, ~calc_metrics(Z),
#'                                    plots, radius = 13)
#'

calc_metrics <- function(z) {
  
  probs <- c(0.01, 0.05, 0.1, 0.2, 0.25,
             0.3, 0.4, 0.5, 0.6, 0.7,
             0.75, 0.8, 0.9, 0.95, 0.99)
  
  zq <- stats::quantile(z, probs, na.rm = T)
  
  list(
    zmean = mean(z, na.rm = T), 
    zsd = sd(z, na.rm = T),
    zmin = min(z, na.rm = T),
    zmax = max(z, na.rm = T),
    zq1 = zq[1], zq5 = zq[2], zq10 = zq[3], zq20 = zq[4],
    zq25 = zq[5], zq30 = zq[6], zq40 = zq[7], zq50 = zq[8],
    zq60 = zq[9], zq70 = zq[10], zq75 = zq[11], zq80 = zq[12],
    zq90 = zq[13], zq95 = zq[14], zq99 = zq[15],
    zskew = moments::skewness(z, na.rm = T),
    zkurt = moments::kurtosis(z, na.rm = T),
    zcv = sd(z, na.rm = T) / mean(z, na.rm = T) * 100,
    zcrr = ((mean(z, na.rm = T) - min(z, na.rm = T)) / (max(z, na.rm = T) - min(z, na.rm = T))),
    pzabove3 = (sum(z > 3, na.rm = T) / length(z)) * 100,
    pzabovezmean = (sum(z > mean(z, na.rm = T), na.rm = T) / length(z)) * 100
  )
  
}