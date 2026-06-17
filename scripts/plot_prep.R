#-------------------------------------------------------------------------------
# Name:         plot_prep.R
# Description:  Preparation of RTK-remeasured forest inventory plots (BI)
#               in the Solling region. Reads and merges the 2023/24 and
#               2025/26 campaigns, applies spatial filters,
#               adds per-plot leaf-off / leaf-on CHM heights and their difference,
#               and supports an automatic and manual filter for plots
#               affected by tree fall or harvest.
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read RTK-remeasured forest inventory plot centers
# from 2023/24 and 2025/26
bi_plots_1 <- sf::st_read(
  file.path(raw_data_dir, 'forest_inventory', 'bi_center_points_solling_2023_24.gpkg')
)

# read RTK-remeasured forest inventory plot centers
bi_plots_2 <- sf::st_read(
  file.path(raw_data_dir, 'forest_inventory', 'bi_center_points_solling_2025_26.gpkg')
)

# read metadata of RTK-remeasured plots from 2023/24
measurement_details <- read.table(
  file.path(meta_data_dir, 'Neueinmessung Solling-etrs89-messdetails.txt'),
  skip = 10,
  header = F,
  fill = T,
  stringsAsFactors = F
)

colnames(measurement_details) <- c(
  'point_nr', 'antenna_height',
  'time', 'day', 'month', 'year',
  'epoch', 'hz_accuracy', 'v_accuracy', 'sv', 'rms'
)

lc_time_old <- Sys.getlocale('LC_TIME')
Sys.setlocale('LC_TIME', 'C')

measurement_details$datetime <- as.POSIXct(
  paste(
    measurement_details$year, measurement_details$month,
    measurement_details$day, measurement_details$time
  ),
  format = '%Y %b %d %H:%M:%S'
)

Sys.setlocale('LC_TIME', lc_time_old)

measurement_details <- measurement_details[
  , !(names(measurement_details) %in% c('time', 'day', 'month', 'year', 'rms'))
]

# read comments table to identify estimated plot centers 
# of RTK-remeasured plots from 2023/24
comments <- readxl::read_excel(
  file.path(meta_data_dir, '2024_01_15_Hauck_GNSS_BI-Punkte_Liste_Kommentare.xlsx')
)



# 02 - data cleaning
#-------------------------------------------------------------------------------

estimated_kspnr <- comments$Pktnr[grepl('ohne Nagel gemessen',
                                        comments$`Kommentar 1`,
                                        ignore.case = T)]

# clean up bi_plots_1:
bi_plots_1$kspnr <- bi_plots_1$KSPNR
bi_plots_1$center_point_estimated <- ifelse(
  bi_plots_1$KSPNR %in% estimated_kspnr, 'yes', 'no'
)
sf::st_geometry(bi_plots_1) <- 'geometry'

# clean up bi_plots_2:
bi_plots_2 <- bi_plots_2[!sf::st_is_empty(bi_plots_2), ]
bi_plots_2$center_point_estimated <- ifelse(grepl('e$', bi_plots_2$Name), 'yes', 'no')
bi_plots_2$kspnr <- as.integer(gsub('[a-zA-Z]$', '', bi_plots_2$Name))
bi_plots_2 <- sf::st_zm(bi_plots_2, drop = T, what = 'ZM')
bi_plots_2 <- bi_plots_2 %>%
  dplyr::arrange(kspnr, dplyr::desc(center_point_estimated)) %>%
  dplyr::filter(!duplicated(kspnr))

# harmonize geometry column naming for compatibility
bi_plots_2 <- dplyr::rename(bi_plots_2, geometry = geom)
sf::st_geometry(bi_plots_2) <- 'geometry'



# 03 - solution status
#-------------------------------------------------------------------------------

# bi_plots_1: derive solution status from hz_accuracy via measurement details
accuracy_per_point <- measurement_details %>%
  arrange(point_nr, hz_accuracy, desc(datetime)) %>%
  distinct(point_nr, .keep_all = T) %>%
  select(point_nr, hz_accuracy, v_accuracy, datetime)

bi_plots_1$point_nr <- as.numeric(substring(bi_plots_1$kspnr, 2))
bi_plots_1 <- dplyr::left_join(bi_plots_1, accuracy_per_point, by = 'point_nr')
bi_plots_1$point_nr <- NULL

bi_plots_1$solution_status <- dplyr::case_when(
  is.na(bi_plots_1$hz_accuracy)    ~ 'single',
  bi_plots_1$hz_accuracy <= 0.1    ~ 'fix',
  TRUE                             ~ 'float'
)

bi_plots_1$measurement_date <- as.Date(bi_plots_1$datetime)

# bi_plots_2: solution status is already in the data
bi_plots_2$solution_status <- tolower(bi_plots_2$Solution)
bi_plots_2$measurement_date <- as.Date(bi_plots_2$Averaging.start)


# exclude single-solution plots from campaign 1 for analysis
bi_plots_1 <- bi_plots_1 %>%
  dplyr::filter(solution_status != 'single')

# merge both RTK datasets
bi_plots <- rbind(
  bi_plots_1[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')],
  bi_plots_2[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')]
)



# 04 - spatial filtering and height difference calculation
#-------------------------------------------------------------------------------

# point cloud paths and acquisition years (leaf-on = 2023, leaf-off = 2024);
# defined here as they are also used by the inspection helpers in 04.4
pc_lon_path  <- file.path(raw_data_dir, 'pc_leafon_2023')
pc_loff_path <- file.path(raw_data_dir, 'pc_leafoff_2024')

year_lon  <- '2023'
year_loff <- '2024'

# cache: the merged plots after the extent filter (04.1) and the angular sector
# coverage filter (04.2), carrying per-plot leaf-off/leaf-on CHM heights and
# their difference (04.3)
merged_file <- file.path(
  processed_data_dir, 'forest_inventory', 'bi_center_points_merged.gpkg'
  )
chm_loff_file <- file.path(processed_data_dir, 'forest_inventory', 'chm_loff_2024.tif')
chm_lon_file  <- file.path(processed_data_dir, 'forest_inventory', 'chm_lon_2023.tif')

if (file.exists(merged_file)) {
  
  # load cached results
  cat('Loading cached filtered plots with height differences from disk...\n')
  bi_plots_height_diff <- sf::st_read(merged_file, quiet = T)
  chm_loff <- terra::rast(chm_loff_file)
  chm_lon  <- terra::rast(chm_lon_file)
  
} else {
  
  cat('No cache found, running the full spatial filtering and extraction...\n')
  
  # 04.1. extent filter
  #---------------------
  
  cat('Extent filter: keeping plots inside the leaf-off footprint...\n')
  
  # read leaf-off point cloud extent
  ext_loff <- sf::st_read(
    file.path(raw_data_dir, 'pc_leafoff_2024', 'leafoff.vpc')
  )
  
  ext_loff <- sf::st_transform(ext_loff, sf::st_crs(bi_plots))
  
  # keep only plots whose 13 m buffered footprint lies fully within the extent
  ext_geom <- sf::st_union(sf::st_buffer(ext_loff, 0.5))
  
  filter_within_extent <- function(plots, ext, radius = 13) {
    buf <- sf::st_buffer(plots, radius)
    plots[lengths(sf::st_within(buf, ext)) > 0, ]
  }
  
  n_before_extent <- nrow(bi_plots)
  #bi_plots_1 <- filter_within_extent(bi_plots_1, ext_geom)
  #bi_plots_2 <- filter_within_extent(bi_plots_2, ext_geom)
  bi_plots   <- filter_within_extent(bi_plots,   ext_geom)
  cat(sprintf('  merged: %d -> %d plots (%d dropped outside extent)\n',
              n_before_extent, nrow(bi_plots), n_before_extent - nrow(bi_plots)))
  
  # 04.2. angular sector coverage filter
  # drops plots whose point cloud has a large empty bin (data gap on one side),
  # even though the plot lies inside the nominal extent
  #-----------------------------------------------------------------------------
  
  cat('Angular sector coverage filter: scanning each plot circle...\n')
  
  laz_dir <- file.path(raw_data_dir, 'pc_leafoff_2024')
  
  laz_files_2024 <- list.files(
    laz_dir, pattern = "_2024\\.laz$", full.names = T
  )
  
  ctg <- lidR::readLAScatalog(laz_files_2024)
  
  angular_sector_coverage <- function(ctg, x, y, radius = 13, n_sectors = 8) {
    
    # clip circular plot from the catalog around (x, y) with the given radius
    las <- lidR::clip_circle(ctg, x, y, radius)
    if (lidR::is.empty(las)) {
      return(list(min_sector_frac = 0))
    }
    
    # total number of points inside the circle
    n_pts <- lidR::npoints(las)
    
    # angle of each point relative to the plot center, in radians (-pi..pi)
    ang <- atan2(las$Y - y, las$X - x)
    
    # assign each point to one of n_sectors equal-width angular bins
    sec <- cut(ang, breaks = seq(-pi, pi, length.out = n_sectors + 1),
               include.lowest = T)
    
    # count how many points fell into each sector
    counts <- as.integer(table(sec))
    
    # ratio of the emptiest sector to the expected share under uniform coverage
    # (1.0 = perfectly uniform, 0 = at least one bin is completely empty)
    min_sector_frac <- min(counts) / (n_pts / n_sectors)
    
    list(min_sector_frac = min_sector_frac)
  }
  
  filter_by_angular_coverage <- function(
    plots, ctg, radius = 13, n_sectors = 8, min_frac = 0
  )
  {
    
    # for each plot, measure how evenly points are distributed around the center
    # a value of 0 means at least one bin is empty -> data gap on one side
    plots$min_sector_frac <- vapply(seq_len(nrow(plots)), function(i) {
      co <- sf::st_coordinates(plots[i, ])
      angular_sector_coverage(ctg, co[1], co[2],
                              radius = radius,
                              n_sectors = n_sectors)$min_sector_frac
    }, numeric(1))
    
    plots[plots$min_sector_frac > min_frac, ]
  }
  
  n_before_coverage <- nrow(bi_plots)
  #bi_plots_1 <- filter_by_angular_coverage(bi_plots_1, ctg)
  #bi_plots_2 <- filter_by_angular_coverage(bi_plots_2, ctg)
  bi_plots   <- filter_by_angular_coverage(bi_plots,   ctg)
  cat(sprintf('  merged: %d -> %d plots (%d dropped for coverage gaps)\n',
              n_before_coverage, nrow(bi_plots),
              n_before_coverage - nrow(bi_plots)))
  
  # 04.3. height difference calculation
  # computes per-plot CHM height for leaf-off and leaf-on acquisitions
  # and stores their difference (leaf-off - leaf-on) as an attribute,
  # so plots affected by treefall or harvest can be identified later
  #-----------------------------------------------------------------------------
  
  extract_plot_heights <- function(
    pc_path,
    plots,
    year,
    res = 0.5,
    buffer_radius = 13)
  {
    
    # plot centre coordinates as regions of interest
    coords <- sf::st_coordinates(plots)
    xc <- coords[, 1]
    yc <- coords[, 2]
    
    # get all LAZ files from the folder, then keep only those of the
    # requested acquisition year (year suffix in the file name)
    laz_files <- list.files(pc_path, pattern = '\\.laz$', full.names = T)
    laz_files <- laz_files[grepl(year, basename(laz_files))]
    if (length(laz_files) == 0) {
      stop(sprintf('No LAZ files matching year "%s" found in %s', year, pc_path))
    }
    cat(sprintf('  %d LAZ files match year %s\n', length(laz_files), year))
    
    # reads and process only the buffered plot circles
    reader <- lasR::reader_circles(xc = xc, yc = yc, r = buffer_radius)
    chm_stage <- lasR::rasterize(res = res, operators = max(HAG))
    na_fill <- lasR::focal(chm_stage, size = 3, fun = 'mean')
    pipeline <- reader + chm_stage + na_fill
    
    # execute pipeline on all relevant LAZ files at once
    ans <- lasR::exec(
      pipeline,
      on = laz_files,
      with = list(ncores = lasR::half_cores(), progress = T)
    )
    
    # the NA-filled CHM is the last stage of the pipeline
    chm <- ans[[length(ans)]]
    
    # set CRS of CHM to match plots
    terra::crs(chm) <- sf::st_crs(plots)$wkt
    
    # extract mean height within each plot (buffered to match the circles)
    plots_buffered <- sf::st_buffer(plots, dist = buffer_radius)
    plot_heights <- exactextractr::exact_extract(
      chm,
      plots_buffered,
      fun = 'mean'
    )
    
    # convert list to vector
    if (is.list(plot_heights)) {
      plot_heights <- unlist(plot_heights)
    }
    
    list(plot_heights = plot_heights, chm = chm)
  }
  
  # extract per-plot mean heights from CHMs of both acquisitions
  cat('Height difference: extracting CHM heights for both acquisitions...\n')
  cat('Extracting heights from leaf-off point cloud...\n')
  loff_result <- extract_plot_heights(pc_loff_path, bi_plots, year = year_loff)
  
  cat('Extracting heights from leaf-on point cloud...\n')
  lon_result  <- extract_plot_heights(pc_lon_path,  bi_plots, year = year_lon)
  
  heights_loff <- loff_result$plot_heights
  heights_lon  <- lon_result$plot_heights
  chm_loff <- loff_result$chm
  chm_lon  <- lon_result$chm
  
  # build merged dataset with height columns
  bi_plots_height_diff <- bi_plots
  bi_plots_height_diff$height_loff <- heights_loff
  bi_plots_height_diff$height_lon  <- heights_lon
  bi_plots_height_diff$height_diff <-
    bi_plots_height_diff$height_loff - bi_plots_height_diff$height_lon
  
  # write the merged plot dataset (with height columns) and the CHMs to disk;
  cat('Writing merged plot dataset and CHMs to disk...\n')
  sf::st_write(bi_plots_height_diff, merged_file, delete_dsn = T)
  terra::writeRaster(chm_loff, chm_loff_file, overwrite = T)
  terra::writeRaster(chm_lon,  chm_lon_file,  overwrite = T)
  
}

# 04.4. manual filtering (merged file only)
# inspects the height-difference distribution to detect ambiguous plots,
# exports per-plot CHM crops and point cloud clips for ambiguous plots
# visual inspection of the point clouds
# confirmed tree-fall/harvest plots are filtered out
#-------------------------------------------------------------------------------

# buffer radius for plot circles
pc_buffer_radius <- 13

# summary statistics of the height differences
diff_vals <- bi_plots_height_diff$height_diff[
  !is.na(bi_plots_height_diff$height_diff)
]
cat('\nHeight difference (leaf-off - leaf-on) summary:\n')
print(summary(diff_vals))

# density plot of the height differences
density_plot <- ggplot2::ggplot(
  data.frame(height_diff = diff_vals),
  ggplot2::aes(x = height_diff)
) +
  ggplot2::geom_density(
    fill = 'grey80',
    colour = 'firebrick',
    linewidth = 0.8,
    alpha = 0.6
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    colour = 'black',
    linewidth = 0.6
  ) +
  ggplot2::geom_vline(
    xintercept = median(diff_vals),
    linewidth = 0.8,
    linetype = 'dashed'
  ) +
  ggplot2::labs(
    title = 'Distribution of per-plot height differences',
    subtitle = 'Solid line: zero    Dashed line: median',
    x = NULL,
    y = 'Density'
  ) +
  ggplot2::theme_minimal()

print(density_plot)

# automatically flag plots with a large leaf-off / leaf-on height difference
# (|height_diff| > 5 m), indicating treefall or harvest between acquisitions;
# the threshold was chosen based on the distribution above
height_diff_threshold <- 5
plots_to_remove_auto <- bi_plots_height_diff$kspnr[
  !is.na(bi_plots_height_diff$height_diff) &
    abs(bi_plots_height_diff$height_diff) > height_diff_threshold
]
cat(sprintf('%d plots exceed the height_diff threshold of %g m: %s\n',
            length(plots_to_remove_auto), height_diff_threshold,
            paste(plots_to_remove_auto, collapse = ', ')))

# ambiguous height-difference range (leaf-off - leaf-on)
# adjusted after inspecting the distribution plot above
ambig_lower <- -5
ambig_upper <- -2.5

ambiguous_plots <- bi_plots_height_diff[
  !is.na(bi_plots_height_diff$height_diff) &
    bi_plots_height_diff$height_diff > ambig_lower &
    bi_plots_height_diff$height_diff < ambig_upper,
]
cat(sprintf('%d plots in ambiguous height_diff range (%g, %g)\n',
            nrow(ambiguous_plots), ambig_lower, ambig_upper))

# output directories for per-plot CHM crops and point cloud clips
plot_chm_out_dir <- file.path(processed_data_dir, 'cropped_chms_plots')
plot_pc_out_dir  <- file.path(processed_data_dir, 'cropped_pcs_plots')

# helper: write per-plot CHM crops to disk
save_plot_chm_crops <- function(chm, plots, out_dir, name_suffix,
                                buffer_radius = pc_buffer_radius,
                                skip_existing = T) {
  
  dir.create(out_dir, recursive = T, showWarnings = F)
  plots_buffered <- sf::st_buffer(plots, dist = buffer_radius)
  n_saved <- 0
  n_skipped <- 0
  
  for (i in seq_len(nrow(plots_buffered))) {
    kspnr_val <- as.character(plots_buffered$kspnr[i])
    out_file <- file.path(out_dir, sprintf('plot_%s_%s.tif', kspnr_val, name_suffix))
    
    if (skip_existing && file.exists(out_file)) {
      n_skipped <- n_skipped + 1
      next
    }
    
    plot_vect <- terra::vect(plots_buffered[i, ])
    chm_crop <- terra::crop(chm, plot_vect)
    chm_mask <- terra::mask(chm_crop, plot_vect)
    terra::writeRaster(chm_mask, out_file, overwrite = T)
    n_saved <- n_saved + 1
  }
  
  cat(sprintf('CHM crop export (%s): saved %d, skipped existing %d (%s)\n',
              name_suffix, n_saved, n_skipped, out_dir))
}

# helper: write per-plot point cloud clips to disk
save_plot_pc_clips <- function(plots, pc_path, year, name_suffix, out_dir,
                               buffer_radius = pc_buffer_radius,
                               skip_existing = T) {
  
  dir.create(out_dir, recursive = T, showWarnings = F)
  laz_files <- list.files(pc_path, pattern = '\\.laz$', full.names = T)
  laz_files <- laz_files[grepl(year, basename(laz_files))]
  if (length(laz_files) == 0) {
    stop(sprintf('No LAZ files matching year "%s" found in %s', year, pc_path))
  }
  ctg <- lidR::readLAScatalog(laz_files, progress = F)
  n_saved <- 0
  n_skipped <- 0
  
  for (i in seq_len(nrow(plots))) {
    kspnr_val <- as.character(plots$kspnr[i])
    out_file <- file.path(out_dir, sprintf('plot_%s_%s.laz', kspnr_val, name_suffix))
    
    if (skip_existing && file.exists(out_file)) {
      n_skipped <- n_skipped + 1
      next
    }
    
    coords <- sf::st_coordinates(plots[i, ])
    las_clip <- lidR::clip_circle(ctg, coords[1, 1], coords[1, 2], buffer_radius)
    lidR::writeLAS(las_clip, out_file)
    n_saved <- n_saved + 1
  }
  
  cat(sprintf('PC clip export (%s): saved %d, skipped existing %d (%s)\n',
              name_suffix, n_saved, n_skipped, out_dir))
}

# export CHM crops and PC clips for the ambiguous plots
if (nrow(ambiguous_plots) > 0) {
  save_plot_chm_crops(chm_loff, ambiguous_plots, plot_chm_out_dir, 'leafoff')
  save_plot_chm_crops(chm_lon, ambiguous_plots, plot_chm_out_dir, 'leafon')
  save_plot_pc_clips(ambiguous_plots, pc_loff_path, year_loff, 'leafoff', plot_pc_out_dir)
  save_plot_pc_clips(ambiguous_plots, pc_lon_path, year_lon,  'leafon',  plot_pc_out_dir)
}

# helper: read the buffered plot circle from all LAZ files of one acquisition
read_plot_pc <- function(pc_path, year, xc, yc, r) {
  laz_files <- list.files(pc_path, pattern = '\\.laz$', full.names = T)
  laz_files <- laz_files[grepl(year, basename(laz_files))]
  if (length(laz_files) == 0) {
    stop(sprintf('No LAZ files matching year "%s" found in %s', year, pc_path))
  }
  ctg <- lidR::readLAScatalog(laz_files, progress = F)
  lidR::clip_circle(ctg, xc, yc, r)
}

# show both point clouds in one window
plot_sidebyside <- function(las1, las2, col_by = 'Z',
                            palette = viridis::viridis, size = 2,
                            top_down = T) {
  rgl::open3d()
  rgl::bg3d(color = 'black')
  rgl::mfrow3d(1, 2, sharedMouse = T)
  rgl::bg3d(color = 'black')
  rgl::plot3d(
    las1@data$X, las1@data$Y, las1@data$Z,
    col = palette(100)[cut(las1@data[[col_by]], breaks = 100)],
    size = size, decorate = F
  )
  if (top_down) rgl::view3d(theta = 0, phi = 0)
  rgl::next3d()
  rgl::bg3d(color = 'black')
  rgl::plot3d(
    las2@data$X, las2@data$Y, las2@data$Z,
    col = palette(100)[cut(las2@data[[col_by]], breaks = 100)],
    size = size, decorate = F
  )
  if (top_down) rgl::view3d(theta = 0, phi = 0)
}

# read + visualise one plot by its kspnr
inspect_plot_by_kspnr <- function(kspnr_val) {
  plot_idx <- which(bi_plots_height_diff$kspnr == kspnr_val)
  if (length(plot_idx) == 0) {
    stop(sprintf('No plot with kspnr = %s found.', as.character(kspnr_val)))
  }
  if (length(plot_idx) > 1) {
    warning(sprintf('Multiple plots with kspnr = %s; using the first.',
                    as.character(kspnr_val)))
    plot_idx <- plot_idx[1]
  }
  plot_geom <- bi_plots_height_diff[plot_idx, ]
  coords <- sf::st_coordinates(plot_geom)
  cat(sprintf('kspnr %s (height_diff = %.2f m)\n',
              as.character(kspnr_val), plot_geom$height_diff))
  
  las_loff <- read_plot_pc(pc_loff_path, year_loff,
                           coords[1, 1], coords[1, 2], pc_buffer_radius)
  las_lon  <- read_plot_pc(pc_lon_path,  year_lon,
                           coords[1, 1], coords[1, 2], pc_buffer_radius)
  plot_sidebyside(las_lon, las_loff)
  
  invisible(list(las_loff = las_loff, las_lon = las_lon))
}

# example: inspect a plot interactively (uncomment to run)
#inspect_plot_by_kspnr(49013)

# kspnr of plots to remove after visual inspection (fill in manually)
plots_to_remove_manual <- c(49013)

# combine the automatic height_diff threshold removals (see section above)
# with the manually confirmed plots
plots_to_remove <- union(plots_to_remove_auto, plots_to_remove_manual)

bi_plots_filtered <- bi_plots_height_diff[
  !bi_plots_height_diff$kspnr %in% plots_to_remove,
]

cat('Removed', length(plots_to_remove),
    'plots total (height_diff threshold + manual inspection)\n')
cat('Filtered merged dataset:',
    nrow(bi_plots_filtered), '/', nrow(bi_plots_height_diff), 'plots retained\n')

# write the final filtered dataset to disk (depends on the 04.4 filtering;
# the other derived outputs are written in section 04 on a full run)
sf::st_write(
  bi_plots_filtered,
  file.path(processed_data_dir, 'forest_inventory', 'bi_center_points_merged_filtered.gpkg'),
  delete_dsn = T
)
