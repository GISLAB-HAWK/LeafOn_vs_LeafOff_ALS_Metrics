#-----------------------------------------------------------------------------
# Name:         clip_las2plots.R
# Description:  takes pre-processed and harmonized las data and clips it to the 
#               extent of the BI plots with a buffer of 15m. these will be used 
#               for ITC in the script ITC_BI.R
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-----------------------------------------------------------------------------

library(lidR)
library(sf)
library(parallel)
#library(future)

# 01 - file path definitions
#----------------------------
setwd("R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/")

# define BI plot directory 
plots_path <- "R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS/BI2trees/data/raw_data/bi_center_points_not_remeasured.gpkg"
buffer <- 20

# define harmonized data directory (leaf-on)
las_dir <- './03_harmonized/lon23/'

# define output directory
out_dir <- './05_plot_clips/lon23'

# create output dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 02 - data reading
#-------------------------------------
#plan(multisession)
ctg <- readLAScatalog(las_dir)

opt_progress(ctg) <- TRUE


plots <- st_read(plots_path)
stopifnot(!is.null(plots))
if (is.na(st_crs(plots))) stop("Plots have no CRS. Please define it.")
if (st_crs(plots) != st_crs(ctg)) plots <- st_transform(plots, st_crs(ctg))


# Require 'kspnr' attribute
if (!("kspnr" %in% names(plots))) stop("Attribute 'kspnr' not found in plots.")


# Buffer plots 
plots_buf <- st_buffer(st_make_valid(plots), buffer)
plots_buf$kspnr <- as.character(plots$kspnr)


# Small helper to make safe filenames from kspnr
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

# 03 - per-plot clipping
#-------------------------------------
i <-  46
## loop over plots 
for (i in seq_len(nrow(plots_buf))) {
  plot_geom <- plots_buf[i, ]
  plot_id   <- plots_buf$kspnr[i]
  
  message("processing: kspnr ", plot_id)
  
  # tiles intersecting with current plot 
  subcat <- catalog_intersect(ctg, plot_geom)
  if (length(subcat) == 0) {
    message("[", plot_id, "] no intersecting tiles — skip")
    return(invisible(NULL))}
  
  clip <- clip_roi(subcat, plot_geom)
  
  if (!length(clip)) {
    message("[", plot_id, "] nothing after clip")
    return(invisible(NULL))
  }
  
  out_file <- file.path(out_dir, paste0(sanitize(plot_id), ".laz"))
  writeLAS(clip, out_file)
  message("[", plot_id, "] → ", out_file)
  invisible(out_file)
  
  
}
