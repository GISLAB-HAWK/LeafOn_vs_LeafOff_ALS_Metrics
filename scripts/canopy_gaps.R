#-------------------------------------------------------------------------------
# Name:         canopy_gaps.R
# Description:  Detection of canopy gaps in airborne laser scanning (ALS)-based
#               canopy height models (CHM) in leaf-on and leaf-off condition.
#               The detected canopy gaps are cropped to forest inventory plots,
#               and the gap fraction per plot is calculated.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - processing setup 
#-------------------------------------------------------------------------------

# point cloud paths and acquisition years (leaf-on = 2023, leaf-off = 2024)
pc_lon_path  <- file.path(raw_data_dir, 'pc_leafon_2023')
pc_loff_path <- file.path(raw_data_dir, 'pc_leafoff_2024')

# output directories per acquisition
chm_dir_2023 <- file.path(processed_data_dir, 'chm_leafon')
chm_dir_2024 <- file.path(processed_data_dir, 'chm_leafoff')
gap_dir_2023 <- file.path(processed_data_dir, 'gap_polygons_leafon')
gap_dir_2024 <- file.path(processed_data_dir, 'gap_polygons_leafoff')

# get all LAZ files per acquisition
laz_files_2023 <- list.files(pc_lon_path,  pattern = '\\.laz$', full.names = T)
laz_files_2024 <- list.files(pc_loff_path, pattern = '\\.laz$', full.names = T)

# CHM resolution
chm_res <- 0.5

# number of cores used to compute the CHMs
n_cores <- 32



# 02 - CHM computation 
#-------------------------------------------------------------------------------

# compute spike-free CHMs (locally adaptive spike-free algorithm,
# Fisher F. J. 2024; freeze_distance = 0)
# point clouds are normalized first so that the rasterized Z is height above ground
compute_spikefree_chm <- function(laz_files, out_dir, res = chm_res) {

  read <- lasR::reader()
  norm <- lasR::normalize()
  chm  <- lasR::spikefree(
    res = res, freeze_distance = 0,
    ofile = file.path(out_dir, '*.tif')
  )

  pipeline <- read + norm + chm

  lasR::exec(
    pipeline, on = laz_files,
    with = list(ncores = n_cores, progress = T)
  )

  # rename outputs: 3dm_..._2023.tif -> chm_..._2023.tif
  produced <- list.files(out_dir, pattern = '^3dm_.*\\.tif$', full.names = T)
  if (length(produced) > 0) {
    file.rename(produced, file.path(out_dir, sub('^3dm_', 'chm_', basename(produced))))
  }
}

# compute the leaf-on (2023) and leaf-off (2024) CHMs
# (skip an acquisition if all its output tiles already exist)
run_spikefree_chm <- function(laz_files, out_dir, label) {
  expected <- file.path(
    out_dir, sub('^3dm_', 'chm_', sub('\\.laz$', '.tif', basename(laz_files)))
  )
  if (length(expected) > 0 && all(file.exists(expected))) {
    cat('All', label, 'CHM tiles already exist - skipping\n')
    return(invisible(NULL))
  }
  cat('Computing', label, 'spike-free CHMs...\n')
  compute_spikefree_chm(laz_files, out_dir)
}

run_spikefree_chm(laz_files_2023, chm_dir_2023, 'leaf-on (2023)')
run_spikefree_chm(laz_files_2024, chm_dir_2024, 'leaf-off (2024)')

# build a virtual mosaic (VRT) from the single CHM tiles,
# or load it if it already exists
chm_lon_mosaic_file  <- file.path(chm_dir_2023, 'chm_leafon_2023.vrt')
chm_loff_mosaic_file <- file.path(chm_dir_2024, 'chm_leafoff_2024.vrt')

load_or_build_chm_mosaic <- function(tile_dir, mosaic_file) {
  if (file.exists(mosaic_file)) {
    cat('Loading existing CHM mosaic:', basename(mosaic_file), '\n')
    return(terra::rast(mosaic_file))
  }
  cat('Building CHM mosaic:', basename(mosaic_file), '\n')
  tiles <- list.files(tile_dir, pattern = '\\.tif$', full.names = T)
  terra::vrt(tiles, filename = mosaic_file, overwrite = T)
}

chm_lon  <- load_or_build_chm_mosaic(chm_dir_2023, chm_lon_mosaic_file)
chm_loff <- load_or_build_chm_mosaic(chm_dir_2024, chm_loff_mosaic_file)

# crop leaf-on CHM to the leaf-off extent
chm_lon <- terra::crop(chm_lon, chm_loff, mask = T)

# quick look
par(mfrow = c(1, 2))
terra::plot(chm_lon,  main = 'leaf-on 2023')
terra::plot(chm_loff, main = 'leaf-off 2024')
par(mfrow = c(1, 1))



# 03 - automatic canopy gap detection
#-------------------------------------------------------------------------------

# source function for gap detection
source('src/detect_gaps_multi_stage.R', local = T)

# define height stages for multi-stage gap detection
stages <- list(
  list(
    gap_height_threshold = 5,
    size = c(10, 5000), 
    buffer_width = 20, 
    percentile_threshold = 10
    ),
  list(
    gap_height_threshold = 10,
    size = c(10, 5000), 
    buffer_width = 20, 
    percentile_threshold = 20
    ),
  list(
    gap_height_threshold = 15,
    size = c(10, 5000),
    buffer_width = 20,
    percentile_threshold = 30
    )
)

# final merged gap-polygon layer per acquisition (written by the function)
gap_file_lon  <- file.path(gap_dir_2023, 'gap_polys_lon.gpkg')
gap_file_loff <- file.path(gap_dir_2024, 'gap_polys_loff.gpkg')

# apply function to leaf-on and leaf-off CHMs
# (skip detection and load the gaps if they already exist)
if (file.exists(gap_file_lon)) {
  cat('Loading existing leaf-on canopy gaps...\n')
  canopy_gaps_lon <- sf::st_read(gap_file_lon)
} else {
  canopy_gaps_lon <- detect_gaps_multi_stage(
    chm = chm_lon,
    stages = stages,
    output_dir = gap_dir_2023,
    area_name = 'lon'
  )
}

if (file.exists(gap_file_loff)) {
  cat('Loading existing leaf-off canopy gaps...\n')
  canopy_gaps_loff <- sf::st_read(gap_file_loff)
} else {
  canopy_gaps_loff <- detect_gaps_multi_stage(
    chm = chm_loff,
    stages = stages,
    output_dir = gap_dir_2024,
    area_name = 'loff'
  )
}



# 04 - gap fraction per forest inventory plot
#-------------------------------------------------------------------------------

# plot metrics (BI plots: geometry + attributes + ALS metrics)
plot_metrics_lon  <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon.gpkg')
  )
plot_metrics_loff <- sf::st_read(
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff.gpkg')
  )

# sample-circle radius (m)
plot_radius <- 13

# gap fraction (clipped gap area / plot-circle area) per plot
compute_gap_fraction <- function(plots, gaps, radius = plot_radius, id_col = 'kspnr') {
  
  # plot circles (fine approximation) and clip the gaps to them
  plots_buf <- sf::st_buffer(plots, dist = radius, nQuadSegs = 90)
  inter <- sf::st_intersection(plots_buf[, id_col], gaps)
  inter$gap_area <- as.numeric(sf::st_area(inter))
  
  # gap area per plot (0 where none), then fraction of the circle area (%)
  gap_sum  <- tapply(inter$gap_area, inter[[id_col]], sum)
  gap_area <- as.numeric(gap_sum[as.character(plots[[id_col]])])
  gap_area[is.na(gap_area)] <- 0
  plots$gap_fraction <- gap_area / (pi * radius^2) * 100

  plots
}

plot_metrics_lon  <- compute_gap_fraction(plot_metrics_lon,  canopy_gaps_lon)
plot_metrics_loff <- compute_gap_fraction(plot_metrics_loff, canopy_gaps_loff)

# optional: write the plot metrics back with the new gap_fraction column
sf::st_write(
  plot_metrics_lon,
  file.path(processed_data_dir, 'metrics', 'plot_metrics_lon.gpkg'), 
  delete_dsn = T
  )

sf::st_write(
  plot_metrics_loff,
  file.path(processed_data_dir, 'metrics', 'plot_metrics_loff.gpkg'),
  delete_dsn = T
  )



# 05 - compare leaf-on vs leaf-off gap fractions
#-------------------------------------------------------------------------------

# long data frame with the gap fraction of both acquisitions
gap_df <- rbind(
  data.frame(season = 'leaf-on',  gap_fraction = plot_metrics_lon$gap_fraction),
  data.frame(season = 'leaf-off', gap_fraction = plot_metrics_loff$gap_fraction)
)
gap_df$season <- factor(gap_df$season, levels = c('leaf-on', 'leaf-off'))

# colors
season_cols <- c('leaf-on' = '#009E73', 'leaf-off' = '#E69F00')

# mean gap fraction per season
gap_means <- aggregate(gap_fraction ~ season, data = gap_df, FUN = mean)

# density plot with dotted vertical lines at the season means
gap_dens <- ggplot2::ggplot(gap_df, ggplot2::aes(gap_fraction, fill = season, color = season)) +
  ggplot2::geom_density(alpha = 0.4) +
  ggplot2::geom_vline(
    data = gap_means,
    ggplot2::aes(xintercept = gap_fraction, color = season),
    linetype = 'dotted', linewidth = 0.8, show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(values = season_cols, name = NULL) +
  ggplot2::scale_color_manual(values = season_cols, name = NULL) +
  ggplot2::labs(
    x = 'Gap fraction [%]',
    y = 'Density'
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = c(0.85, 0.85)
  )

print(gap_dens)


























