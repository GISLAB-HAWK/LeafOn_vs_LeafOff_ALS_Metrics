# calc_metrics.R

#' Count occupied voxels over a ladder of box sizes
#'
#' Box-counting on the raw coordinates. The voxel grid is anchored globally at 0
#' (floor(x / s)) and is deliberately NOT re-origined per region: min(X) comes
#' from the point cloud and differs between the seasons, which would make the
#' grid move with the season.
#' The scale ladder is derived from the bounding box of the points, halved
#' n_halvings times and cut off at min_size.
#'
#' @param xyz three-column matrix of coordinates (X, Y, Z), Z = height above ground.
#' @param min_size numeric, smallest box size (m).
#' @param n_halvings integer, number of halvings of the scale ladder.
#' @param min_height numeric, points at or below this height are dropped (m).
#'
#' @return list with the box sizes (r) and the occupied voxel counts (N),
#'         or NULL if the scale ladder is too short.
#' @author Svenja Dobelmann

count_voxels <- function(xyz, min_size = 0.2, n_halvings = 14, min_height = 0.5) {

  xyz <- xyz[xyz[, 3L] > min_height, , drop = FALSE]
  if (nrow(xyz) < 2L) return(NULL)

  rngs <- c(diff(range(xyz[, 1L])),
            diff(range(xyz[, 2L])),
            diff(range(xyz[, 3L])))

  box_size_max <- max(round(rngs * 100) / 100)
  if (!is.finite(box_size_max) || box_size_max <= 0) return(NULL)

  sizes <- box_size_max / (2^(0:n_halvings))
  sizes <- sizes[sizes >= min_size]

  # 3, because the coarsest scale is dropped when fitting the box dimension
  if (length(sizes) < 3L) return(NULL)

  N <- integer(length(sizes))
  for (i in seq_along(sizes)) {
    s <- sizes[i]
    N[i] <- data.table::uniqueN(
      data.table::data.table(vi = as.integer(xyz[, 1L] %/% s),
                             vj = as.integer(xyz[, 2L] %/% s),
                             vk = as.integer(xyz[, 3L] %/% s))
    )
  }

  list(r = sizes, N = N)
}


#' Box dimension from a box-counting result
#'
#' Slope of log(N) ~ log(1/r), fitted with ordinary least squares in closed form.
#' The coarsest scale is dropped: it is the starting box with N = 1, i.e. the
#' point (0, 0), which would otherwise anchor the fit.
#' Relative vs. absolute scaling of r is irrelevant, since log(1/r) only shifts
#' x by an offset and the OLS slope is unaffected by that.
#'
#' @param r numeric vector of box sizes.
#' @param N integer vector of occupied voxel counts.
#' @param drop_coarsest logical, drop the coarsest scale before fitting.
#'
#' @return numeric, the box dimension (slope), or NA.
#' @author Svenja Dobelmann

box_dimension <- function(r, N, drop_coarsest = TRUE) {

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

  # OLS slope
  sum((x - xm) * (y - ym)) / den
}


#' Effective number of layers (Hill numbers)
#'
#' Hill numbers of order 0/1/2 over the distribution of the occupied voxels
#' across height bins:
#'   richness -> number of occupied layers        (Hill q = 0)
#'   shannon  -> exp(Shannon entropy)             (Hill q = 1)
#'   simpson  -> 1 / Simpson index                (Hill q = 2)
#' It always holds that richness >= shannon >= simpson, which is a useful
#' consistency check.
#'
#' The voxel grid is FIXED (voxel_res) and independent of the box-counting scale
#' ladder, so that the values are comparable between regions and between seasons.
#' It is anchored globally at 0 for the same reason.
#'
#' @param xyz three-column matrix of coordinates (X, Y, Z), Z = height above ground.
#' @param voxel_res numeric, fixed voxel size (m).
#' @param height_bin numeric, thickness of the height layers (m).
#' @param min_height numeric, points at or below this height are dropped (m).
#'
#' @return named numeric vector with richness, shannon and simpson.
#' @author Svenja Dobelmann

effective_layers <- function(xyz, voxel_res = 0.2, height_bin = 1,
                             min_height = 0.5) {

  na_out <- c(richness = NA_real_, shannon = NA_real_, simpson = NA_real_)

  xyz <- xyz[xyz[, 3L] > min_height, , drop = FALSE]
  if (nrow(xyz) < 2L) return(na_out)

  s <- voxel_res

  vox <- unique(
    data.table::data.table(vi = as.integer(xyz[, 1L] %/% s),
                           vj = as.integer(xyz[, 2L] %/% s),
                           vk = as.integer(xyz[, 3L] %/% s))
  )

  # voxel Z as cell centre -> height layer
  hc <- as.integer(round((vox$vk * s + s / 2) / height_bin))

  n <- tabulate(hc - min(hc) + 1L)
  n <- n[n > 0L]
  if (!length(n)) return(na_out)

  p <- n / sum(n)

  c(richness = as.numeric(length(n)),
    shannon  = exp(-sum(p * log(p))),
    simpson  = 1 / sum(p^2))
}


#' Rumple index of the canopy surface
#'
#' Ratio of the canopy surface area to its ground area, computed on the TIN of
#' the decimated surface points (highest point per surf_res cell), which stands
#' in for a rasterised CHM. Values are >= 1.
#'
#' @param X,Y numeric vectors of coordinates.
#' @param Z numeric vector of heights above ground.
#' @param surf_res numeric, cell size of the decimation grid (m).
#'
#' @return numeric, the rumple index, or NA.
#' @author Svenja Dobelmann

rumple <- function(X, Y, Z, surf_res = 0.5) {

  if (length(Z) < 4L) return(NA_real_)

  # decimate to the highest point per surf_res cell (CHM surrogate)
  dt <- data.table::data.table(
    gx = as.integer(X %/% surf_res),
    gy = as.integer(Y %/% surf_res),
    X = X, Y = Y, Z = Z
  )
  surf <- dt[dt[, .I[which.max(Z)], by = list(gx, gy)]$V1]

  if (nrow(surf) < 4L) return(NA_real_)

  # a degenerate footprint cannot be triangulated by qhull
  if (diff(range(surf$X)) <= 1e-6 || diff(range(surf$Y)) <= 1e-6) return(NA_real_)

  # recentre: rumple is translation invariant, but qhull loses precision at
  # UTM magnitudes (~5.5e5)
  xc <- surf$X - mean(surf$X)
  yc <- surf$Y - mean(surf$Y)

  if (nrow(unique(cbind(xc, yc))) <= 3L) return(NA_real_)

  tryCatch(suppressMessages(lidR::rumple_index(xc, yc, surf$Z)),
           error = function(e) NA_real_)
}


#' Compute ALS metrics from a set of points
#'
#' The function is shared by the plot-level and the wall-to-wall (pixel-level)
#' metric calculation so that the predictors are defined identically for model
#' training and for prediction. The region area is passed via the area argument,
#' which makes the function applicable to both:
#'   plot_metrics()  -> area = pi * radius^2  (sample circle)
#'   pixel_metrics() -> area = res^2          (pixel)
#' Do not change a formula here without recomputing BOTH the plot metrics and
#' the wall-to-wall raster.
#'
#' Point sets differ per metric group:
#'   height metrics                 -> vegetation points only (Classification != 2)
#'   pr_*, rd_*, point_density      -> all points (ground included)
#'   pulse_returns_mean             -> first returns only (ReturnNumber == 1)
#'   structural metrics             -> all points, see below
#'
#' Metric definitions:
#'   mean, max, var, sd, kurtosis, skew, p5 ... p99
#'                      height distribution of the vegetation points
#'                      (kurtosis is the excess kurtosis)
#'   pr_A_B             penetration (pass through) rate of the 1 m layer
#'                      (A, B]: of the points that reached below the top of
#'                      the layer, the fraction that also passed through it,
#'                      #{z <= A} / #{z <= B}
#'   rd_A_B             return density of the 1 m layer (A, B]:
#'                      #{A < z <= B} / #{all points}
#'   point_density      points per m²
#'   pulse_returns_mean mean number of returns per laser pulse
#'   box_dimension      fractal box dimension of the point cloud, from
#'                      box-counting over a ladder of box sizes
#'   vci                vertical complexity index: Shannon entropy of the point
#'                      heights normalised to [0, 1]
#'   rumple             ratio of canopy surface area to ground area (>= 1)
#'   enl_richness       effective number of layers, Hill q = 0 (occupied layers)
#'   enl_shannon        effective number of layers, Hill q = 1
#'   enl_simpson        effective number of layers, Hill q = 2
#'
#' Structural metrics use all points above min_height, except vci, which uses
#' all points from vci_min_height upwards. Negative heights are clamped to 0
#' rather than dropped: lidR::VCI() returns NA as soon as min(z) < 0, and
#' dropping those points would remove systematically more points from leaf-off
#' (more ground penetration) than from leaf-on, i.e. a bias in exactly the
#' factor under test.
#'
#' The input point clouds are expected to be noise filtered already and to
#' carry a HAG (height above ground) attribute.
#'
#' @param X,Y numeric vectors of coordinates (needed for box_dimension and rumple).
#' @param Z numeric vector of heights above ground (HAG).
#' @param Classification integer vector of point classifications.
#' @param ReturnNumber integer vector of return numbers.
#' @param NumberOfReturns integer vector of the number of returns per pulse.
#' @param area numeric, area of the region in m² (used for point_density):
#'        the circle area for plot_metrics(), the pixel area for pixel_metrics().
#' @param layers numeric vector of the lower bounds of the 1 m layers used for
#'        the pr_* / rd_* metrics. Layer XX covers (XX, XX+1] and is named after
#'        both bounds, so layers = c(2, 10) yields pr_2_3, pr_10_11, rd_2_3 and
#'        rd_10_11.
#' @param min_height numeric, height filter for box_dimension and the effective
#'        number of layers (m).
#' @param box_min_size numeric, smallest box size for the box dimension (m).
#' @param n_halvings integer, halvings of the box-counting scale ladder.
#' @param voxel_res numeric, fixed voxel size for the effective number of layers (m).
#' @param height_bin numeric, height layer thickness for the effective number
#'        of layers (m).
#' @param surf_res numeric, cell size of the decimation grid for the rumple index (m).
#' @param vci_zmax numeric, upper height bound used to normalize the vci.
#' @param vci_min_height numeric, lower height bound for the vci. 0 keeps all
#'        returns including ground; the ground peak is markedly larger in
#'        leaf-off, so the vci difference then partly measures penetration
#'        rather than stand complexity. A value of 2 keeps vegetation returns
#'        only. Decide deliberately and state it in the methods.
#'
#' @return named list of metrics.
#' @author Svenja Dobelmann, Florian Franz

calc_metrics <- function(X, Y, Z, Classification, ReturnNumber, NumberOfReturns,
                         area, layers = c(2, 10),
                         min_height = 0.5, box_min_size = 0.2, n_halvings = 14,
                         voxel_res = 0.2, height_bin = 1, surf_res = 0.5,
                         vci_zmax = 45, vci_min_height = 0) {

  # remove NAs consistently across all vectors
  valid <- !is.na(X) & !is.na(Y) & !is.na(Z) & !is.na(Classification) &
           !is.na(ReturnNumber) & !is.na(NumberOfReturns)
  X <- X[valid]
  Y <- Y[valid]
  Z <- Z[valid]
  Classification <- Classification[valid]
  ReturnNumber <- ReturnNumber[valid]
  NumberOfReturns <- NumberOfReturns[valid]

  n_all <- length(Z)
  if (n_all == 0) return(NULL)

  # all points (ground included) vs vegetation points only
  Z_all <- Z
  Z_veg <- Z[Classification != 2]
  n_veg <- length(Z_veg)

  # percentiles computed for p*
  probs <- c(5, 10, 20, 25, 30, 40, 50, 60, 70, 75, 80, 90, 95, 99)

  # --- height metrics (vegetation points only, ground excluded) ---
  if (n_veg > 1 && !is.na(sd(Z_veg)) && sd(Z_veg) > 0) {

    h_perc <- setNames(
      as.list(as.numeric(quantile(Z_veg, probs / 100))),
      paste0('p', probs)
    )

    h_stats <- c(
      list(
        mean     = mean(Z_veg),
        max      = max(Z_veg),
        var      = var(Z_veg),
        sd       = sd(Z_veg),
        kurtosis = (sum((Z_veg - mean(Z_veg))^4) / n_veg) / sd(Z_veg)^4 - 3,
        skew     = (sum((Z_veg - mean(Z_veg))^3) / n_veg) / sd(Z_veg)^3
      ),
      h_perc
    )

  } else {

    h_names <- c('mean', 'max', 'var', 'sd', 'kurtosis', 'skew',
                 paste0('p', probs))
    h_stats <- setNames(as.list(rep(NA_real_, length(h_names))), h_names)

  }

  # --- penetration (pass through) rate of layer (lower, upper] (all points) ---
  # NA when no point reached below the top of the layer
  pr <- setNames(
    as.list(sapply(layers, function(lower) {
      upper <- lower + 1
      reach <- sum(Z_all <= upper)
      if (reach == 0) return(NA_real_)
      sum(Z_all <= lower) / reach
    })),
    paste0('pr_', layers, '_', layers + 1)
  )

  # --- return density of layer (lower, upper] (all points) ---
  rd <- setNames(
    as.list(sapply(layers, function(lower) {
      upper <- lower + 1
      sum(Z_all > lower & Z_all <= upper) / n_all
    })),
    paste0('rd_', layers, '_', layers + 1)
  )

  # --- point_density (all points) ---
  point_density <- n_all / area

  # --- pulse_returns_mean (first returns only) ---
  first_idx <- ReturnNumber == 1
  pulse_returns_mean <- if (any(first_idx)) {
    mean(NumberOfReturns[first_idx])
  } else {
    NA_real_
  }

  # --- structural metrics (all points, negative heights clamped to 0) ---
  Z_struct <- pmax(Z_all, 0)
  xyz <- cbind(X, Y, Z_struct)

  # box dimension: scale ladder derived from the bounding box of the points
  vox <- count_voxels(xyz, min_size = box_min_size, n_halvings = n_halvings,
                      min_height = min_height)
  box_dim <- if (is.null(vox)) NA_real_ else box_dimension(vox$r, vox$N)

  # effective number of layers: fixed voxel grid, independent of the ladder above
  enl <- effective_layers(xyz, voxel_res = voxel_res, height_bin = height_bin,
                          min_height = min_height)

  # vertical complexity index
  z_vci <- if (vci_min_height > 0) Z_struct[Z_struct >= vci_min_height] else Z_struct
  vci <- if (length(z_vci) < 2L) {
    NA_real_
  } else {
    tryCatch(lidR::VCI(z_vci, zmax = vci_zmax, by = 1), error = function(e) NA_real_)
  }

  # rumple index of the canopy surface
  rum <- rumple(X, Y, Z_struct, surf_res = surf_res)

  c(
    h_stats,
    pr,
    rd,
    list(
      point_density      = point_density,
      pulse_returns_mean = pulse_returns_mean,
      box_dimension      = box_dim,
      vci                = vci,
      rumple             = rum,
      enl_richness       = unname(enl[['richness']]),
      enl_shannon        = unname(enl[['shannon']]),
      enl_simpson        = unname(enl[['simpson']])
    )
  )
}
