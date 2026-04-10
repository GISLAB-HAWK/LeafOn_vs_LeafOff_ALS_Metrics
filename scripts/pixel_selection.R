#-------------------------------------------------------------------------------
# Name:         pixel_selection.R
# Description:  interactive script that plots the extent of randomly selected   
#               metrik pixels (10x10m) and shows the surface model in a higher 
#               resolution 0,5m. User can then set the validity of the pixel. 
#               Goal is to get a sample of 100 valid pixel localtions for the 
#               further analysis.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)


# Output
out_csv <- file.path(output_dir, "sample_selection.csv")

#### Parameter Setup ####

set.seed(123)

n_valid_target <- 100   # Ziel: 100 valide Samples
pixel_size <- 10
buffer_m <- 3
diff_min <- -15
diff_max <- 15

#### Load data ####

lon_r <- rast(
    file.path(processed_data_dir, 'metrics','pix_level','solling23_lon_ppm20_lt20_indices.tiff')
  )
  
loff_r <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling24_loff_ppm20_lt20_indices.tiff')
)

lon_dsm <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling23_lon_ppm20_lt20_dsm05.tiff')
)

loff_dsm <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling24_loff_ppm20_lt20_dsm05.tiff')
)


#### Pre-select suitable canditate pixels ####

# crop rasters to same extent 
lon_r <- crop(lon_r, loff_r)
lon_dsm <- crop(lon_dsm, ext(loff_dsm))

# mask only valid cells
# 1. masking out NA
is_na <- is.na(lon_r) | is.na(loff_r)
mask_na <- app(is_na, fun = function(x) any(x, na.rm = TRUE))   

# 2. masking out height difference >10m (assuming clearcut)
diff <- loff_r$BE_H_MAX - lon_r$BE_H_MAX
mask_diff <- diff <= -10

# 3. masking out non forested areas 
mask_nontree <- lon_r$ts == 0

# combine the three masks 
mask_combined <- mask_na | mask_diff | mask_nontree

# apply the mask 
loff_r <- mask(loff_r, mask_combined, maskvalue = TRUE)
lon_r <- mask(lon_r, mask_combined, maskvalue = TRUE)


valid_10m <- !is.na(lon_r) & !is.na(loff_r) 

cand <- as.data.frame(valid_10m, xy = TRUE, cells = TRUE, na.rm = TRUE) %>%
  dplyr::select(cell, x, y, ts ) %>%
  rename( valid = ts) %>%
  filter(valid == "TRUE")

if (nrow(cand) == 0) {
  stop("No valid 10-m-Pixel found.")
}


# random order of pixels 
cand <- cand[sample(seq_len(nrow(cand))), ]
row.names(cand) <- NULL

#### Helper functions ####

make_extent <- function(x, y, pixel_size, buffer_m) {
  h <- pixel_size / 2
  ext(
    x - h - buffer_m, x + h + buffer_m,
    y - h - buffer_m, y + h + buffer_m
  )
}

make_pixel_polygon <- function(x, y, pixel_size, crs_str = "") {
  h <- pixel_size / 2
  xy <- rbind(
    c(x-h, y-h),
    c(x+h, y-h),
    c(x+h, y+h),
    c(x-h, y+h),
    c(x-h, y-h)
  )
  vect(list(xy), type = "polygons", crs = crs_str)
}

save_results <- function(results_list, out_csv) {
  if (length(results_list) > 0) {
    df <- do.call(rbind, results_list)
    write.csv(df, out_csv, row.names = FALSE)
  }
}

# Difference-DSM
diff_dsm <- loff_dsm - lon_dsm

#### Interactive selection until 100 pixels are selected ####
results <- list()
n_valid <- 0
candidate_index <- 1
review_counter <- 0

while (n_valid < n_valid_target && candidate_index <= nrow(cand)) {
  
  review_counter <- review_counter + 1
  
  x <- cand$x[candidate_index]
  y <- cand$y[candidate_index]
  cell_id <- cand$cell[candidate_index]
  
  e <- make_extent(x, y, pixel_size, buffer_m)
  
  lon_dsm_crop <- crop(lon_dsm, e)
  loff_dsm_crop <- crop(loff_dsm, e)
  diff_crop  <- crop(diff_dsm, e)
  
  pix_poly <- make_pixel_polygon(x, y, pixel_size, crs(lon_r))
  
  # skipp empty pixels automaticcally
  if (ncell(diff_crop) == 0) {
    results[[length(results) + 1]] <- data.frame(
      review_id = review_counter,
      cell = cell_id,
      x = x,
      y = y,
      lon_r = cand$lon_r[candidate_index],
      loff_r = cand$loff_r[candidate_index],
      status = "skipped_empty"
    )
    candidate_index <- candidate_index + 1
    save_results(results, out_csv)
    next
  }
  
  # Plot
  op <- par(no.readonly = TRUE)
  par(mfrow = c(1, 3), mar = c(3, 3, 3, 5))
  
  plot(lon_dsm_crop, main = paste0("DSM Leaf-on "), range = minmax(lon_dsm_crop))
  lines(pix_poly, col = "red", lwd = 2)

  plot(loff_dsm_crop, main = "DSM Leaf-off", range = minmax(lon_dsm_crop))
  lines(pix_poly, col = "red", lwd = 2)
  
  plot(diff_crop, main = "Diff (lon-loff)", range = c(diff_min, diff_max))
  lines(pix_poly, col = "red", lwd = 2)
  
  par(op)
  
  cat("\n")
  cat("====================================================\n")
  cat("Review:", review_counter, "\n")
  cat("already valid:", n_valid, "out of", n_valid_target, "\n")
  cat("Candidate:", candidate_index, "von", nrow(cand), "\n")
  cat("Coordinate:", round(x, 2), round(y, 2), "\n")
  cat("Entry: v = valid | n = invalid | s = skip | q = quit\n")
  
  ans <- tolower(trimws(readline("Your Decision: ")))
  
  if (ans == "v") {
    status <- "valid"
    n_valid <- n_valid + 1
  } else if (ans == "n") {
    status <- "invalid"
  } else if (ans == "s") {
    status <- "skipped"
  } else if (ans == "q") {
    cat("\nAborted by user.\n")
    break
  } else {
    status <- "skipped"
    cat("invalid entry -> saved as skipped.\n")
  }
  
  results[[length(results) + 1]] <- data.frame(
    review_id = review_counter,
    cell = cell_id,
    x = x,
    y = y,
    status = status
  )
  
  save_results(results, out_csv)
  candidate_index <- candidate_index + 1
}

#### Save results ####

save_results(results, out_csv)

df_final <- do.call(rbind, results)
n_valid_final <- sum(df_final$status == "valid")

cat("\n")
cat("Finished\n")
cat("Valid Samples:", n_valid_final, "\n")
cat("saved in:", out_csv, "\n")

if (n_valid_final < n_valid_target && candidate_index > nrow(cand)) {
  cat("Note: Not enough valid candidates found to meet", n_valid_target,
      "valid samples.\n")
}



################################################################################