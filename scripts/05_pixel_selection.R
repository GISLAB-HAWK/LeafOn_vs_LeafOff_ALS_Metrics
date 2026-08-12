#-------------------------------------------------------------------------------
# Name:         pixel_selection.R
# Description:  interactive script that plots the extent of randomly selected   
#               metrik pixels (10x10m) and shows the surface model in a higher 
#               resolution 0,5m. User can then set the validity of the pixel. 
#               Goal is to get a sample of 400 valid pixel localtions for the 
#               further analysis.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)

# Output
out_csv <- file.path(metadata_dir, "sample_selection_n400_balanced.csv")

#### Parameter Setup ####

set.seed(123)

n_valid_target <- 400   
n_per_class <- n_valid_target / 2  # 200 coniferous, 200 deciduous
pixel_size <- 10
buffer_m <- 3
diff_min <- -15
diff_max <- 15

#### Load data ####

loff_r <- rast(
  file.path(
  processed_data_dir, "metrics_leafoff_2024", "ppm20", "solling24_loff_ppm20_metrics.tiff")
)

lon_r <- rast(
  file.path(
  processed_data_dir, "metrics_leafon_2023", "ppm20", "solling23_lon_ppm20_metrics.tiff")
)

loff_dsm <- rast(
  file.path(metadata_dir, "dsm", "solling23_lon_ppm20_dsm05.tiff")
)

lon_dsm <- rast(
  file.path(metadata_dir, "dsm", "solling24_loff_ppm20_dsm05.tiff")
)

## tree species data 
ts_r <- rast(
  file.path(metadata_dir, "tree_species", "CLMS_DLT2023_clipped.tif")
)


#### Pre-select suitable canditate pixels ####

# crop rasters to same extent 
lon_r <- crop(lon_r, loff_r)

# mask only valid cells
# 1. masking out NA
is_na <- is.na(lon_r) | is.na(loff_r)
mask_na <- app(is_na, fun = function(x) any(x, na.rm = TRUE))   

# 2. masking out height difference >10m (assuming clearcut)
diff <- loff_r$BE_H_MAX - lon_r$BE_H_MAX
mask_diff <- abs(diff) >= 5

# 3. masking out non forested areas 
mask_nontree <- lon_r$class_name == 0

# combine the three masks 
mask_combined <- mask_na | mask_diff | mask_nontree

# apply the mask 
loff_r <- mask(loff_r, mask_combined, maskvalue = TRUE)
lon_r <- mask(lon_r, mask_combined, maskvalue = TRUE)


valid_10m <- !is.na(lon_r) & !is.na(loff_r) 

cand <- as.data.frame(lon_r$class_name, xy = TRUE, cells = TRUE, na.rm = TRUE)

ts_vals <- terra::extract(lon_r$class_name, cand[, c("x", "y")])

cand <- cand %>%
  mutate(species = case_when(
    ts_vals$class_name == 1 ~ "deciduous",
    ts_vals$class_name == 2 ~ "coniferous",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(species)) %>%
  # random order within each species group
  group_by(species) %>%
  slice_sample(prop = 1) %>%  
  ungroup()

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
    df <- dplyr::bind_rows(results_list)
    write.csv(df, out_csv, row.names = FALSE)
  }
}

# Difference-DSM
e <- terra::intersect(ext(lon_dsm), ext(loff_dsm))
diff_dsm <- crop(loff_dsm, e) - crop(lon_dsm, e)

#### Interactive selection until n pixels are selected ####
results <- list()
candidate_index <- 1
review_counter <- 0
n_valid_deciduous   <- 0
n_valid_coniferous  <- 0

while ((n_valid_deciduous < n_per_class || n_valid_coniferous < n_per_class) && 
       candidate_index <= nrow(cand)) {
  
  current_species <- cand$species[candidate_index]
  
  # skip if this species class is already full
  if ((current_species == "deciduous"  && n_valid_deciduous  >= n_per_class) |
      (current_species == "coniferous" && n_valid_coniferous >= n_per_class)) {
    candidate_index <- candidate_index + 1
    next
  }
  
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
      cell      = cell_id,
      x         = x,
      y         = y,
      species   = current_species,   
      status    = "skipped_empty"
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
  
  cat("====================================================\n")
  cat("Review:", review_counter, "\n")
  cat("Species:", current_species, "\n")   # <- add this line
  cat("Valid deciduous: ",  n_valid_deciduous,  "/", n_per_class, "\n")
  cat("Valid coniferous:", n_valid_coniferous, "/", n_per_class, "\n")
  cat("Candidate:", candidate_index, "of", nrow(cand), "\n")
  cat("Coordinate:", round(x, 2), round(y, 2), "\n")
  cat("Entry: v = valid | n = invalid | s = skip | q = quit\n")
  
  ans <- tolower(trimws(readline("Your Decision: ")))
  
  if (ans == "v") {
    status <- "valid"
    if (current_species == "deciduous")  n_valid_deciduous  <- n_valid_deciduous  + 1
    if (current_species == "coniferous") n_valid_coniferous <- n_valid_coniferous + 1
  } else if (ans == "n") {
    status <- "invalid"
  } else if (ans == "s") {
    status <- "skipped"
  } else if (ans == "q") {
    cat("\nAborted by user.\n")
    break
  } else {
    status <- "skipped"
  }
  
  results[[length(results) + 1]] <- data.frame(
    review_id = review_counter,
    cell      = cell_id,
    x         = x,
    y         = y,
    species   = current_species,   
    status    = status
  )
  
  save_results(results, out_csv)
  candidate_index <- candidate_index + 1
}

#### Save results ####

save_results(results, out_csv)

message("pixel selection written to: ",out_csv)

################################################################################
