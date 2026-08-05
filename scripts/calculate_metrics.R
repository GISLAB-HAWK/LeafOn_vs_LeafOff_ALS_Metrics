#-------------------------------------------------------------------------------
# Name:         calculate_metrics.R
# Description:  Calculating 24 ABA Metrics from ALS data using the lidR pixel_metrics() 
#               function and las Catalog object 
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R', local = TRUE)
source('src/be_metrics_func.R', local = TRUE) # load the funtion for calculating metrics


# 01 - configuration
#-------------------

# seasons and point densites (point per metersquare)
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

THRESHOLDS  <- c(2, 5, 10, 20, 30)   # Height thresholds for RD und PR [m]
CHUNK_SIZE  <- 500
CHUNK_BUFF  <- 20
DROP_CLASS  <- "-drop_class 7 18"


# Reference raster: defines resolution, extent and raster gitter 
raster_template <- rast(
  file.path(metadata_dir, "tree_species", "CLMS_DLT2023_clipped.tif")
)

res_val   <- res(raster_template)[1]
start_val <- c(xmin(raster_template), ymin(raster_template))


# 02 - worker function
#---------------------

#' Compute ABA metrics for one season / point density combination
#'
#' @param season Character. Name of the season, must be one of names(seasons).
#' @param ppm Character. Point density subfolder, e.g. "ppm20".
#' @param overwrite Logical. Recompute even if the output file already exists.
#'
#' @return Invisibly, the path to the written raster, or NULL on failure.
#' 
calc_metrics_task <- function(season, ppm, overwrite = FALSE) {
  
  cfg <- seasons[[season]]
  
  in_dir   <- file.path(processed_data_dir, cfg$in_dir, ppm)
  out_dir  <- file.path(processed_data_dir, cfg$out_dir, ppm)
  out_name <- paste0(season, "_", ppm, "_metrics.tiff")
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
  opt_select(ctg)       <- "*"
  
  m <- pixel_metrics(
    ctg,
    ~be_metrics(HAG, Classification, ReturnNumber, NumberOfReturns,
                thresholds = THRESHOLDS,
                cell_area  = res_val^2),
    res   = res_val,
    start = start_val
  )
  
  m   <- crop(extend(m, raster_template), raster_template)
  out <- c(m, raster_template)
  
  writeRaster(out, out_file, overwrite = TRUE)
  message("written: ", out_file)
  
  invisible(out_file)
}


# 03 - run over all combinations
#-------------------------------

n_cores <- max(1, parallel::detectCores() - 4)
plan(multisession, workers = n_cores)

results <- list()

for (season in names(seasons)) {
  for (ppm in ppms) {
    
    name <- paste(season, ppm, sep = "_")
    
    results[[name]] <- tryCatch(
      calc_metrics_task(season = season, ppm = ppm),
      error = function(e) {
        warning("failed for ", name, ": ", conditionMessage(e))
        NULL
      }
    )
  }
}

plan(sequential)

