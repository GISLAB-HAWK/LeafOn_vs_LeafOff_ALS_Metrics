#-------------------------------------------------------------------------------
# Name:         calculate_structural_metrics.R
# Description:  Calculating structural ALS metrics (CHM_mean, Canopy_cover,
#               Rumple, ENL0D/1D/2D, box dimension, VCI) from ALS data using the
#               lidR catalog_apply() function and LAScatalog objects, for all
#               combinations of season and point density.
#               Functions and their rationale: src/functions_structural_metrics.R
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R', local = TRUE)
source('src/struct_metrics_func.R', local = TRUE) # metric functions


# 01 - configuration
#-------------------

# seasons and point densities (points per square meter)
seasons <- list(
  solling23_lon = list(
    in_dir  = "pc_leafon_2023",
    out_dir = "metrics_leafon_2023"
  ),
  solling24_loff = list(
    in_dir  = "pc_leafoff_2024",
    out_dir = "metrics_leafoff_2024"
  )
)

ppms <- c("ppm20", "ppm10", "ppm4")

CHUNK_SIZE <- 500
CHUNK_BUFF <- 20
DROP_CLASS <- "-drop_class 7 18"
SELECT     <- "xyz1"   # HAG = extra byte 1, "xyz0" would read all extra bytes
N_WORKERS  <- 12       # limited by RAM, not by cores

# Metric parameters, used as defaults in src/functions_structural_metrics.R
ZMAX_FIX        <- 45      # VCI upper bound (fixed -> same scaling for all cells)
COVER_THRESHOLD <- 5       # canopy cover threshold [m]
SURF_RES        <- 0.5     # decimation grid for the surface metrics [m]
MIN_H           <- 0.5     # height filter for voxels / box dimension / ENL [m]
VOXEL_RES       <- 0.2     # fixed voxel size for ENL [m]
HEIGHT_BIN      <- 1       # height layer thickness for ENL [m]
BOX_MIN_SIZE    <- 0.2     # smallest box counting scale [m]
N_HALVINGS      <- 14
Z_MAX_VALID     <- 60      # upper plausibility filter for HAG (noise)
VCI_MIN_H       <- 0

# Reference raster: defines resolution, extent and raster grid.
# Native 10 m in EPSG:3035, reprojected to EPSG:25832, hence res = 9.9968 m 
raster_template <- rast(
  file.path(metadata_dir, "tree_species", "CLMS_2023_Solling_merged_25832.tif")
)

res_val   <- res(raster_template)[1]
start_val <- c(xmin(raster_template), ymin(raster_template))


# 02 - worker function
#---------------------

#' Compute structural metrics for one season / point density combination
#'
#' @param season Character. Name of the season, must be one of names(seasons).
#' @param ppm Character. Point density subfolder, e.g. "ppm20".
#' @param overwrite Logical. Recompute even if the output file already exists.
#'
#' @return Invisibly, the path to the written raster, or NULL on failure.
#'
calc_structural_task <- function(season, ppm, overwrite = FALSE) {
  
  cfg <- seasons[[season]]
  
  in_dir   <- file.path(processed_data_dir, cfg$in_dir, ppm)
  out_dir  <- file.path(processed_data_dir, cfg$out_dir, ppm)
  out_name <- paste0(season, "_", ppm, "_struct_metrics.tiff")
  out_file <- file.path(out_dir, out_name)
  
  if (file.exists(out_file) && !overwrite) {
    message("skipping (exists): ", out_name)
    return(invisible(out_file))
  }
  if (!dir.exists(in_dir)) {
    warning("input directory not found: ", in_dir)
    return(invisible(NULL))
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("processing: ", season, " / ", ppm)
  
  ctg <- readLAScatalog(in_dir)
  opt_chunk_size(ctg)   <- CHUNK_SIZE
  opt_chunk_buffer(ctg) <- CHUNK_BUFF
  opt_filter(ctg)       <- DROP_CLASS
  opt_select(ctg)       <- SELECT
  opt_output_files(ctg) <- ""   # keep results in RAM
  
  # catalog_apply() instead of pixel_metrics(): each chunk needs Z set from HAG
  # and two passes, point based and on the decimated cloud.
  m <- catalog_apply(
    ctg, chunk_fun,
    res_out = res_val, start_out = start_val, surf_res = SURF_RES,
    .options = list(
      raster_alignment = list(res = res_val, start = start_val),
      automerge        = TRUE
    )
  )
  
  m   <- align_to_template(m, raster_template)
  out <- c(m, raster_template)   # add the tree species class for the comparison
  
  writeRaster(out, out_file, overwrite = TRUE)
  message("written: ", out_file)
  
  invisible(out_file)
}


# 03 - run over all combinations
#-------------------------------

plan(multisession, workers = N_WORKERS)

results <- list()

for (season in names(seasons)) {
  for (ppm in ppms) {
    
    name <- paste(season, ppm, sep = "_")
    
    results[[name]] <- tryCatch(
      calc_structural_task(season = season, ppm = ppm),
      error = function(e) {
        warning("failed for ", name, ": ", conditionMessage(e))
        NULL
      }
    )
  }
}

plan(sequential)