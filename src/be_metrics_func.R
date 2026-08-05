#' Calculate Metrics from ALS catalogue 
#'
#' calculates 24 ABA metrics from LAS files 
#'
#' @param Z Height attribute of LAS file 
#' @param Classification PC Attribute indexing classification into gound and non-ground
#' @param ReturnNumber PC Attribute indexing the return number 
#' @param NumberOfReturns PC Attribute indexing the total number of returns
#' 
#' @return Raster file with multiple bands containing the metrics 
#'
#' @details
#' The function can be used within  LidR::pixel_metrics() 
#'
#' @examples
#'pixel_metrics(
#'  ctg_loff,
#'  ~be_metrics(HAG, Classification, ReturnNumber, NumberOfReturns),
#'  res = res_val,
#'  start = start_val
#')#'
 

be_metrics <- function(Z, Classification, ReturnNumber, NumberOfReturns, thresholds, cell_area) {
  
  #  --- remove NAs (over all attributes)  ---
  valid <- !is.na(Z) & !is.na(Classification) & !is.na(ReturnNumber) & !is.na(NumberOfReturns)
  Z <- Z[valid]
  Classification <- Classification[valid]
  ReturnNumber <- ReturnNumber[valid]
  NumberOfReturns <- NumberOfReturns[valid]
  
  n_all <- length(Z)
  if (n_all == 0) return(NULL)
  
  Z_all <- Z
  Z_veg <- Z[Classification != 2]
  n_veg <- length(Z_veg)
  
  # --- BE_H_* (ground excluded) ---
  if (n_veg > 1 && !is.na(sd(Z_veg)) && sd(Z_veg) > 0) {
    h_stats <- list(
      BE_H_KURTOSIS = (sum((Z_veg - mean(Z_veg))^4) / n_veg) / sd(Z_veg)^4 - 3,
      BE_H_SKEW     = (sum((Z_veg - mean(Z_veg))^3) / n_veg) / sd(Z_veg)^3,
      BE_H_VAR      = var(Z_veg),
      BE_H_SD       = sd(Z_veg),
      BE_H_MAX      = max(Z_veg),
      BE_H_MEAN     = mean(Z_veg),
      BE_H_P10      = as.numeric(quantile(Z_veg, 0.10)),
      BE_H_P20      = as.numeric(quantile(Z_veg, 0.20)),
      BE_H_P50      = as.numeric(quantile(Z_veg, 0.50)),
      BE_H_P80      = as.numeric(quantile(Z_veg, 0.80)),
      BE_H_P90      = as.numeric(quantile(Z_veg, 0.90)),
      BE_H_P95      = as.numeric(quantile(Z_veg, 0.95))
    )
  } else {
    h_stats <- setNames(as.list(rep(NA_real_, 12)),
                        c("BE_H_KURTOSIS","BE_H_SKEW","BE_H_VAR","BE_H_SD","BE_H_MAX",
                          "BE_H_MEAN","BE_H_P10","BE_H_P20","BE_H_P50","BE_H_P80",
                          "BE_H_P90","BE_H_P95"))
  }
  
  # --- PR ---
  pr <- setNames(
    as.list(sapply(thresholds, function(lo) {
      hi <- lo + 1
      denom <- sum(Z_all <= hi)
      if (is.na(denom) || denom == 0) return(NA_real_)
      sum(Z_all <= lo) / denom
    })),
    paste0("BE_PR_", sprintf("%02d", thresholds))
  )
  
  # --- RD ---
  rd <- setNames(
    as.list(sapply(thresholds, function(lo) {
      hi <- lo + 1
      sum(Z_all > lo & Z_all <= hi) / n_all
    })),
    paste0("BE_RD_", sprintf("%02d", thresholds))
  )
  
  # --- point_density, pulse_returns_mean ---
  point_density <- n_all / cell_area
  first_idx <- ReturnNumber == 1
  pulse_returns_mean <- if (any(first_idx, na.rm = TRUE)) mean(NumberOfReturns[first_idx], na.rm = TRUE) else NA_real_
  
  c(list(point_density = point_density),
    h_stats,
    pr,
    rd,
    list(pulse_returns_mean = pulse_returns_mean))
}
