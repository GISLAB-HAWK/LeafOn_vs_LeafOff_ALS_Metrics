#-----------------------------------------------------------------------------
# Name:         ITC_BI.R
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

# define plot clip laz directory
las_dir <- './05_plot_clips/lon23/'

plots_path <- "R:/AG_Magdon/datensaetze/solling/dobelmann/Repo_LeafOn_vs_LeafOff_ALS/BI2trees/data/raw_data/bi_center_points_not_remeasured.gpkg"


# define output directory
out_dir <- 'R:/AG_Magdon/datensaetze/solling/dobelmann/ITC/'

# create output dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


#######

# Segmentation params
res_chm <- 0.5 # CHM resolution (m)
min_ht <- 2 # minimum tree height (m)
wsfun <- function(x) {        # variable window for LMF
  y <- 2.6 * (-(exp(-0.08*(x-2)) - 1)) + 3 
  # from https://r-lidar.github.io/lidRbook/itd.html
  y[x < 2] <- 3
  y[x > 20] <- 5
  return(y)
}

dal_max_cr <- 12 # max crown radius (m); set <= buffer - 1 if you used 15 m buffer
#concavity <- 2 # concave hull tightness; set NULL for convex hulls

# 02 - data reading
#-------------------------------------
#plan(multisession)
ctg <- readLAScatalog(las_dir)

opt_progress(ctg) <- TRUE

plots <- st_read(plots_path)
# Create 13 m buffers around plot centers
plots_buf13 <- st_buffer(plots, 13)


# Small helper to make safe filenames from kspnr
sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  trimws(x)
}

# 03 - per-plot segmentation
#-------------------------------------
# Per-file segmentation returning sf POLYGON crowns


seg_file <- function(chunk) {
  las <- readLAS(chunk)
  
  if (is.empty(las)) return(NULL)
  las_n <- normalize_height(las, tin())
  las_n <- filter_poi(las_n, Z >= 0)
  
  
  chm <- rasterize_canopy(las_n, res = res_chm, algorithm = dsmtin(max_edge = 8))
  kern <- matrix(1, 3, 3)
  chm_s <- terra::focal(chm, w = kern, fun = mean, na.rm = T)
  
  
  tt <- locate_trees(las_n, lmf(ws = wsfun, hmin = min_ht))
  if (is.null(tt) || nrow(tt) == 0) return(NULL)
  
  
  alg <- dalponte2016(chm = chm_s, treetops = tt, th_tree = min_ht, max_cr = dal_max_cr)
  seg <- segment_trees(las_n, alg)
  if (is.null(seg) || is.empty(seg)) return(NULL)
  
  crowns <- crown_metrics(seg,.stdtreemetrics, geom = "convex")
  return(crowns)
  
  
}


# Run over the catalog
crowns_list <- catalog_apply(ctg, seg_file)

crowns_all <- do.call(rbind, crowns_list)

# Spatial intersection: keep crowns that intersect plot buffers
crowns_clipped <- st_intersection(crowns_all, plots_buf13)

plot(crowns_clipped[3])

write_sf(crowns_clipped,"BI_crowns5.gpkg")


