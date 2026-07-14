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

# get all LAZ files per acquisition
laz_files_2023 <- list.files(pc_lon_path,  pattern = '\\.laz$', full.names = T)
laz_files_2024 <- list.files(pc_loff_path, pattern = '\\.laz$', full.names = T)

# CHM resolution
chm_res <- 0.5

# number of cores used to compute the CHMs
n_cores <- 32

# output directories per acquisition (one CHM tile per input LAZ tile)
chm_dir_2023 <- file.path(processed_data_dir, 'chm_leafon_2023')
chm_dir_2024 <- file.path(processed_data_dir, 'chm_leafoff_2024')



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

# quick look
par(mfrow = c(1, 2))
terra::plot(chm_lon,  main = 'leaf-on 2023')
terra::plot(chm_loff, main = 'leaf-off 2024')
par(mfrow = c(1, 1))



