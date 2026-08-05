#-------------------------------------------------------------------------------
# Name:         struct_metrics_func.R
# Description:  Functions used by calculate_structural_metrics.R,
#               compare_structural_metrics.R and its visualization script:
#               structural ALS metrics (box dimension, ENL, VCI, surface
#               metrics, catalog chunk function) and the paired tests on the
#               extracted pixels.
#
#               The metric functions take their defaults from the constants
#               defined in calculate_structural_metrics.R: BOX_MIN_SIZE,
#               N_HALVINGS, MIN_H, VOXEL_RES, HEIGHT_BIN, ZMAX_FIX, VCI_MIN_H,
#               COVER_THRESHOLD and Z_MAX_VALID.
#
#               Methodological decisions encoded here:
#               - ENL on a fixed VOXEL_RES grid anchored at 0, decoupled from
#                 the box counting scale ladder, so values stay comparable
#                 between cells and seasons.
#               - Box counting on raw coordinates, no re-origining per cell:
#                 min(X) differs between seasons and would shift the grid.
#               - Box dimension drops the coarsest scale (N = 1).
#               - Canopy cover is surface based, not the 0.2 m column method,
#                 which is point density dependent.
#               - Negative HAG is clamped to 0 instead of dropped: dropping
#                 removes more leaf-off than leaf-on points and biases exactly
#                 the factor under test. lidR::VCI() also returns NA if z < 0.
#
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


# ---- structural metrics -------------------------------------------------------

#' Count occupied voxels over a ladder of box sizes
#'
#' Scale ladder derived from the bounding box of the cell points (Seidel logic).
#' The voxel grid is anchored at 0 globally. The former ULS formulation
#' round((x + s/2)/s)*s is identical up to a cell label to floor(x/s) and gives
#' the same voxel count.
#'
#' @param xyz Numeric matrix with columns X, Y, Z.
#' @param min_size Smallest box size [m].
#' @param n_halvings Number of times the largest box size is halved.
#' @param min_height_filter Minimum height [m] of the points used.
#'
#' @return List with box sizes r and voxel counts N, or NULL if fewer than
#'   three scales remain (three, because the coarsest one is dropped later).
count_voxels <- function(xyz, min_size = BOX_MIN_SIZE, n_halvings = N_HALVINGS,
                         min_height_filter = MIN_H) {
  
  xyz <- xyz[xyz[, 3L] > min_height_filter, , drop = FALSE]
  if (nrow(xyz) < 2L) return(NULL)
  
  rngs <- c(diff(range(xyz[, 1L])),
            diff(range(xyz[, 2L])),
            diff(range(xyz[, 3L])))
  
  box_size_max <- max(round(rngs * 100) / 100)
  if (!is.finite(box_size_max) || box_size_max <= 0) return(NULL)
  
  sizes <- box_size_max / (2^(0:n_halvings))
  sizes <- sizes[sizes >= min_size]
  if (length(sizes) < 3L) return(NULL)
  
  N <- integer(length(sizes))
  for (s_i in seq_along(sizes)) {
    s <- sizes[s_i]
    N[s_i] <- uniqueN(data.table(vi = as.integer(xyz[, 1L] %/% s),
                                 vj = as.integer(xyz[, 2L] %/% s),
                                 vk = as.integer(xyz[, 3L] %/% s)))
  }
  
  list(r = sizes, N = N)
}


#' Box dimension as the OLS slope of log(N) ~ log(1/r)
#'
#' Relative vs absolute ruler scaling is irrelevant: log(1/ruler) is only an
#' x offset and does not change the slope.
#'
#' @param r Box sizes.
#' @param N Voxel counts.
#' @param drop_coarsest Drop the starting box (N = 1, point (0,0)).
#'
#' @return Slope, or NA if fewer than two usable scales remain.
calculate_box_dimension <- function(r, N, drop_coarsest = TRUE) {
  if (length(r) < 2L || length(N) < 2L) return(NA_real_)
  
  x <- log(1 / r)
  y <- log(N)
  
  if (drop_coarsest) { x <- x[-1L]; y <- y[-1L] }
  
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if (length(x) < 2L) return(NA_real_)
  
  xm <- mean(x); ym <- mean(y)
  den <- sum((x - xm)^2)
  if (den <= 0) return(NA_real_)
  
  ## OLS slope
  sum((x - xm) * (y - ym)) / den
}


#' Effective number of layers as Hill numbers of order 0, 1 and 2
#'
#'   ENL0D = number of occupied layers  (richness)
#'   ENL1D = exp(Shannon entropy)       (Hill q = 1)
#'   ENL2D = 1 / Simpson                (Hill q = 2)
#' ENL0D >= ENL1D >= ENL2D always holds and is a useful consistency check.
#'
#' Own fixed grid, decoupled from the box counting. Taking the voxels from the
#' finest box counting scale made the voxel size vary between 0.2 and 0.4 m
#' with stand height, with a jump at Zmax = 25.6 m, so ENL values were not
#' comparable between cells or seasons.
#'
#' @param xyz Numeric matrix with columns X, Y, Z.
#' @param voxel_resolution Voxel edge length [m].
#' @param height_bin_size Height layer thickness [m].
#' @param min_height_filter Minimum height [m] of the points used.
#'
#' @return Named numeric vector ENL0D, ENL1D, ENL2D.
calculate_enl <- function(xyz, voxel_resolution = VOXEL_RES,
                          height_bin_size = HEIGHT_BIN,
                          min_height_filter = MIN_H) {
  
  na_out <- c(ENL0D = NA_real_, ENL1D = NA_real_, ENL2D = NA_real_)
  
  xyz <- xyz[xyz[, 3L] > min_height_filter, , drop = FALSE]
  if (nrow(xyz) < 2L) return(na_out)
  
  s <- voxel_resolution
  
  vox <- unique(data.table(vi = as.integer(xyz[, 1L] %/% s),
                           vj = as.integer(xyz[, 2L] %/% s),
                           vk = as.integer(xyz[, 3L] %/% s)))
  
  # voxel centre height -> height layer
  hc <- as.integer(round((vox$vk * s + s / 2) / height_bin_size))
  
  n <- tabulate(hc - min(hc) + 1L)
  n <- n[n > 0L]
  if (!length(n)) return(na_out)
  
  p <- n / sum(n)
  
  c(ENL0D = as.numeric(length(n)),
    ENL1D = exp(-sum(p * log(p))),
    ENL2D = 1 / sum(p^2))
}


#' Point based metrics for one raster cell
#'
#' Box dimension and ENL are deliberately independent of each other: own scale
#' ladder vs fixed voxel grid.
#'
#' @param X,Y,Z Point coordinates, Z is height above ground.
#' @param min_size,n_halvings Box counting parameters.
#' @param voxel_resolution,height_bin_size ENL parameters.
#' @param zmax,vci_min_h VCI parameters.
#'
#' @return List of metrics plus the diagnostic layers Zmax and n_scales
#'   (for check_metrics_plausibility.R), or NULL for cells with fewer than
#'   ten points.
point_metrics_fun <- function(X, Y, Z,
                              min_size = BOX_MIN_SIZE, n_halvings = N_HALVINGS,
                              voxel_resolution = VOXEL_RES,
                              height_bin_size = HEIGHT_BIN,
                              zmax = ZMAX_FIX, vci_min_h = VCI_MIN_H) {
  
  ok <- !is.na(X) & !is.na(Y) & !is.na(Z)
  X <- X[ok]; Y <- Y[ok]; Z <- Z[ok]
  if (length(Z) < 10L) return(NULL)
  
  xyz <- cbind(X, Y, Z)
  
  # --- box dimension: bbox dependent scale ladder ---
  vox     <- count_voxels(xyz, min_size = min_size, n_halvings = n_halvings,
                          min_height_filter = MIN_H)
  box_dim <- if (is.null(vox)) NA_real_ else calculate_box_dimension(vox$r, vox$N)
  
  # --- ENL: fixed voxel grid ---
  enl <- calculate_enl(xyz, voxel_resolution = voxel_resolution,
                       height_bin_size = height_bin_size,
                       min_height_filter = MIN_H)
  
  # --- VCI: normalized Shannon entropy of the point heights, range [0, 1] ---
  z_vci <- if (vci_min_h > 0) Z[Z >= vci_min_h] else Z
  vci   <- if (length(z_vci) < 2L) NA_real_ else
    tryCatch(VCI(z_vci, zmax = Z_MAX_VALID, by = 1), error = function(e) NA_real_)
  
  list(Box_dimension = box_dim,
       ENL0D         = unname(enl[1]),
       ENL1D         = unname(enl[2]),
       ENL2D         = unname(enl[3]),
       VCI           = vci,
       Zmax          = max(Z),
       n_scales      = if (is.null(vox)) NA_real_ else as.numeric(length(vox$r)))
}


#' Surface metrics for one raster cell
#'
#' Input is the cloud decimated to SURF_RES (highest point per cell), which
#' replaces the former wall to wall p2r CHM:
#'   CHM_mean     -> mean of the surface heights
#'   Canopy_cover -> share of surface cells above the threshold, [0, 1]
#'   Rumple       -> rumple_index on the TIN of the surface points, >= 1
#'
#' Canopy cover is deliberately different from the ULS script, whose 0.2 m
#' column method with a bbox denominator is point density dependent.
#'
#' @param X,Y,Z Point coordinates of the decimated cloud.
#' @param cover_threshold Height threshold [m] for canopy cover.
#'
#' @return List with CHM_mean, Canopy_cover and Rumple, or NULL.
surface_metrics_fun <- function(X, Y, Z, cover_threshold = COVER_THRESHOLD) {
  
  ok <- !is.na(X) & !is.na(Y) & !is.na(Z)
  X <- X[ok]; Y <- Y[ok]; Z <- Z[ok]
  if (length(Z) < 4L) return(NULL)
  
  # Rumple only on a non degenerate footprint. Edge cells of the flight
  # sometimes hold a single 0.5 m column -> min == max -> qhull cannot
  # triangulate. Recentring keeps precision at UTM magnitudes (~5.5e5),
  # rumple itself is translation invariant.
  rum <- NA_real_
  if (diff(range(X)) > 1e-6 && diff(range(Y)) > 1e-6) {
    xc <- X - mean(X)
    yc <- Y - mean(Y)
    if (nrow(unique(cbind(xc, yc))) > 3L) {
      rum <- tryCatch(suppressMessages(rumple_index(xc, yc, Z)),
                      error = function(e) NA_real_)
    }
  }
  
  list(CHM_mean     = mean(Z),
       Canopy_cover = mean(Z > cover_threshold),
       Rumple       = rum)
}


#' Compute all metrics for one catalog chunk, in a single pass
#'
#' @param chunk LAScluster passed by catalog_apply().
#' @param res_out,start_out Resolution and origin of the output grid.
#' @param surf_res Decimation resolution for the surface metrics [m].
#'
#' @return SpatRaster of all metrics, cropped to the unbuffered chunk.
chunk_fun <- function(chunk, res_out, start_out, surf_res) {
  
  las <- readLAS(chunk)
  if (is.empty(las)) return(NULL)
  
  # Guard: without HAG,  las@data$Z <- NULL  would DELETE the Z column in a
  # data.table instead of emptying it -- silently and without an error.
  if (!"HAG" %in% names(las@data)) {
    message("Kein HAG-Attribut im Chunk gefunden.")
    return(NULL)
  }
  
  las@data[["Z"]] <- pmax(las@data[["HAG"]], 0)   # clamp negative HAG, do not drop
  las <- filter_poi(las, Z <= Z_MAX_VALID)
  if (is.empty(las)) return(NULL)
  
  # --- point based ---
  m_pts <- pixel_metrics(las, ~point_metrics_fun(X, Y, Z),
                         res = res_out, start = start_out)
  
  # --- surface based (decimated cloud = CHM replacement) ---
  las_surf <- decimate_points(las, highest(surf_res))
  m_surf   <- pixel_metrics(las_surf, ~surface_metrics_fun(X, Y, Z),
                            res = res_out, start = start_out)
  
  rm(las, las_surf); gc(verbose = FALSE)
  
  if (is.null(m_pts) || is.null(m_surf)) return(NULL)
  
  # align extents (same grid, but possibly different extent)
  common <- terra::union(terra::ext(m_pts), terra::ext(m_surf))
  m_pts  <- terra::extend(m_pts,  common)
  m_surf <- terra::extend(m_surf, common)
  out    <- c(m_pts, m_surf)
  
  # remove the buffer (st_bbox returns the UNbuffered chunk bbox)
  bb  <- sf::st_bbox(chunk)
  out <- terra::crop(out, terra::ext(bb[["xmin"]], bb[["xmax"]],
                                     bb[["ymin"]], bb[["ymax"]]))
  out
}


#' Align a metrics raster to the template grid
#'
#' Everything is computed with start = START_OUT and raster_alignment, so the
#' grids already match -> extend/crop instead of resample. resample(method =
#' "near") would silently shift by up to half a cell on any misalignment.
#'
#' @param r SpatRaster of the metrics.
#' @param template SpatRaster defining grid, extent and CRS.
#'
#' @return SpatRaster on the template grid.
align_to_template <- function(r, template) {
  r <- terra::extend(r, template)
  r <- terra::crop(r, template)
  terra::crs(r) <- terra::crs(template)
  r
}
