#-------------------------------------------------------------------------------
# Name:         pc_summary.R
# Description:   leaf-on and leaf-off catalogs are clipped to their
#               common extent, height normalized and summarised with a
#               lasR pipeline. Results are stored as RDS per season.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R', local = TRUE)


# 01 - configuration
#-------------------

# Point density of the input point clouds.
# Used for both input paths and output file names.
PPM <- "ppm20"

# Bin width in meters for the Z histogram
ZBIN <- 2

# Cores used by the lasR pipeline
N_CORES <- max(1, parallel::detectCores() - 4)

# Input point clouds per season: subdirectory and output basename
inputs <- list(
  leaf_on = list(
    dir  = "pc_leafon_2023",
    file = paste0("lon23_", PPM, "_pc_description.rds")
  ),
  leaf_off = list(
    dir  = "pc_leafoff_2024",
    file = paste0("loff24_", PPM, "_pc_description.rds")
  )
)

lon_dir  <- file.path(processed_data_dir, inputs$leaf_on$dir,  PPM)
loff_dir <- file.path(processed_data_dir, inputs$leaf_off$dir, PPM)


# 02 - data preparation
#----------------------

lon_ctg  <- readLAScatalog(lon_dir)
loff_ctg <- readLAScatalog(loff_dir)

# Bring both datasets to the same extent
bb_lon  <- st_bbox(lon_ctg)
bb_loff <- st_bbox(loff_ctg)

# Overlapping bounding box of both seasons
xmin <- max(bb_loff["xmin"], bb_lon["xmin"])
ymin <- max(bb_loff["ymin"], bb_lon["ymin"])
xmax <- min(bb_loff["xmax"], bb_lon["xmax"])
ymax <- min(bb_loff["ymax"], bb_lon["ymax"])

# Stop early if the two catalogs do not overlap
stopifnot(xmin < xmax, ymin < ymax)

# Extent filter passed to the lasR reader
filter_arg <- sprintf("-keep_xy %f %f %f %f", xmin, ymin, xmax, ymax)


# 03 - lasR pipeline
#-------------------

read <- lasR::reader_las(filter = filter_arg)   # reader with extent filter
norm <- lasR::normalize()                       # height normalization
summ <- lasR::summarise(metrics = "count", zwbin = ZBIN)

# Build pipeline from stages
pipeline <- read + norm + summ


# 04 - calculate statistics
#--------------------------

ans_lon <- lasR::exec(
  pipeline,
  on       = lon_ctg,
  progress = TRUE,
  ncores   = N_CORES
)

ans_loff <- lasR::exec(
  pipeline,
  on       = loff_ctg,
  progress = TRUE,
  ncores   = N_CORES
)


#--- export results ------------------------------------------------------------

saveRDS(ans_lon, file = file.path(output_dir, 'stats',  inputs$leaf_on$file))
saveRDS(ans_loff, file = file.path(output_dir, 'stats',  inputs$leaf_off$file))

message("point cloud descriptions written to: ", file.path(output_dir, 'stats'))
