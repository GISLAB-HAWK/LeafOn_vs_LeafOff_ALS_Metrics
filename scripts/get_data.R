# workflow to get harmonized and thinned LASfiles
# due to differences in data between leaf-on (2023) and leaf-off (2024) campaigns it is hard to compare the data
# we harmonize data (same attributes, same data scheme), remove noise, and apply a pulse based thinning


source("config.R")

# pak::pkg_install("wiesehahn/managelidar")
# pak::pak("r-lidar/lasR@devel")
library(fs)
library(lasR)
library(lidR)


#### define functions

# function to write Virtual Point Cloud
get_vpc <- function(input, save_vpc_to = tempfile(fileext = ".vpc"), crs = 25832, use_gpstime = TRUE, absolute_path = TRUE) {
  if (!file_exists(save_vpc_to)) {
    vpc <- exec(
      set_crs(crs) + write_vpc(ofile = save_vpc_to, use_gpstime = use_gpstime, absolute_path = absolute_path),
      with = list(ncores = concurrent_files(half_cores())),
      on = input
    )
    return(vpc)
  }
}


# function to run process if files not present
process_missing_outputs <- function(input, output_dir, process_fn, ...) {
  if (length(input) == 1 && dir.exists(input)) {
    # Input is a directory
    input_files <- list.files(input, pattern = "*.laz$", full.names = TRUE)
  } else {
    # Input is a vector of files
    input_files <- input
  }
  
  output_files <- file.path(output_dir, basename(input_files))
  
  # Find missing outputs
  missing <- !file.exists(output_files)
  
  if (any(missing)) {
    dir_create(output_dir)
    missing_files <- input_files[missing]
    n_missing <- length(missing_files)
    
    for (i in seq_along(missing_files)) {
      file <- missing_files[i]
      message(sprintf("Processing file %d/%d: %s", i, n_missing, basename(file)))
      process_fn(file, ...)
    }
  }
  
  invisible(output_files)
}


# write VPC referencing intersecting raw data if not present
if(!all(file.exists(c(vpc_leafon_raw, vpc_leafoff_raw)))){
  # get intersecting tiles
  raw_files <- managelidar::get_intersection(folder_leafon_raw, folder_leafoff_raw, full.names = TRUE)
  get_vpc(input = raw_files[[1]]$filename, save_vpc_to = vpc_leafon_raw)
  get_vpc(input = raw_files[[2]]$filename, save_vpc_to = vpc_leafoff_raw)
}


#### standardize data

# function to classify ground using lastools
classify_ground <- function(las_files){
  for(file in las_files) {
    output_file <- path(folder_classified, path_file(file))
    # lastools need to be installed to run classification
    command <- paste("C:/LAStools/bin/lasground_new64.exe -step 5 -sub 4 -i", path_abs(file), "-o", path_abs(output_file)
                     #, "-non_ground_unchanged"
    )
    system(command)
  }
}

# classify files if not present
raw_files <- c(managelidar::get_names(vpc_leafon_raw, full.names = TRUE),
               managelidar::get_names(vpc_leafoff_raw, full.names = TRUE))
process_missing_outputs(raw_files, folder_classified, function(files) classify_ground(files))


# pipeline to homogenize data
# delete errorneous points, drop unnecessary attributes, classify noise, calc Height above ground, save in same format
# (requires lasR devel until now)
pipeline <- reader() + 
  set_crs(25832) +
  # some files showed invalid points with ReturnNumber and/or NumberOfReturns of 0
  delete_points(filter = "ReturnNumber == 0") +
  delete_points(filter = "NumberOfReturns == 0") +
  keep_attributes(c("Intensity", "ReturnNumber", "NumberOfReturns", "Classification", "UserData", "ScanAngle", "gpstime", "PointSourceID", "ScannerChannel", "EdgeOfFlightline", "ScanDirectionFlag")) +
  classify_with_ivf() +
  hag() +
  edit_attribute(filter = "HAG < -0.1", attribute = "Classification", value = 7) +
  edit_attribute(filter = "HAG > 63", attribute = "Classification", value = 18) +
  write_las(ofile = path("data/02_harmonized/*.laz"), pdrf = 6)

# homogenize files if not present
process_missing_outputs(folder_classified, folder_harmonized, function(files) exec(pipeline, on = files))


# print density of harmoized files
# average pulse density is above 26 pulses/m² for most leaf-on tiles and above 30 for most leaf-of tiles
# list.files(folder_harmonized, pattern = "*.laz$") |> managelidar::get_density()




###### replace homogenize function with fixed one temporarily (based on Chatgpt)

# --- 1) Keep a backup you can restore later
.orig_homogenize <- getFromNamespace("homogenize", "lidR")

# --- 2) Your replacement
homogenize_new <- function(density, res = 5, use_pulse = FALSE)
{
  assert_is_a_number(density)
  assert_all_are_positive(density)
  assert_is_a_bool(use_pulse)
  assert_is_a_number(res)
  assert_all_are_positive(res)
  
  density   <- lazyeval::uq(density)
  res       <- lazyeval::uq(res)
  use_pulse <- lazyeval::uq(use_pulse)
  
  f = function(las)
  {
    assert_is_valid_context(LIDRCONTEXTDEC, "homogenize")
    
    if (use_pulse & !"pulseID" %in% names(las))
    {
      warning("No 'pulseID' attribute found. Decimation by points is used.")
      use_pulse <- FALSE
    }
    
    pulseID <- NULL
    
    n       <- round(density*res^2)
    layout  <- raster_layout(las, res)
    
    if (use_pulse){
      
      # too many pulses in some areas
      last_returns <- las@data[ReturnNumber == NumberOfReturns]
      cells   <- get_group(layout, last_returns)
      selected_pulse_ids <- last_returns[, pulseID[.selected_pulses(pulseID, n)], by = cells]$V1
      return(las@data[, pulseID %in% selected_pulse_ids])
      
    } else
      cells   <- get_group(layout, las)
    return(las@data[, .I[.selected_pulses(1:.N, n)], by = cells]$V1)
  }
  
  f <- plugin_decimate(f)
  return(f)
}

# IMPORTANT: make the function live in lidR's namespace so it can "see" internals
ns <- asNamespace("lidR")
environment(homogenize_new) <- ns

# --- 3) Patch the namespace (affects lidR internal calls)
if (bindingIsLocked("homogenize", ns)) unlockBinding("homogenize", ns)
assignInNamespace("homogenize", homogenize_new, ns = "lidR")
lockBinding("homogenize", ns)

# --- 4) Refresh the exported binding (affects lidR::homogenize and attached package env)
pkg_env <- as.environment("package:lidR")
if (bindingIsLocked("homogenize", pkg_env)) unlockBinding("homogenize", pkg_env)
assign("homogenize", homogenize_new, envir = pkg_env)
lockBinding("homogenize", pkg_env)

# --- 5) Quick sanity checks
stopifnot(identical(getFromNamespace("homogenize","lidR"), homogenize_new))
stopifnot(identical(lidR:::homogenize, homogenize_new))
stopifnot(identical(lidR::homogenize, homogenize_new))


#### reduce point density
# function to apply pulse based point sampling on virtual point cloud
thin_by_pulses <- function(vpc, den) {
  
  ctg <- lidR::readLAScatalog(vpc)
  lidR::opt_chunk_size(ctg) <- 0
  lidR::opt_output_files(ctg) <- path("data/03_thinned/", paste0("ppm", den), "{ORIGINALFILENAME}")
  lidR::opt_laz_compression(ctg) <- TRUE
  
  sample_pulses <- function(las) {
    filtered <- las |> 
      lidR::filter_poi(Synthetic_flag == FALSE)  |> 
      lidR::filter_poi(Classification != 7L & Classification != 18L)  |> 
      lidR::retrieve_pulses() |> 
      lidR::decimate_points(lidR::homogenize(density = den, res = 2, use_pulse = TRUE))
    return(filtered)
  }
  
  lidR::catalog_map(ctg, sample_pulses)
}

# apply pulse based thinning if files not present
files_harmonized <- list.files(folder_harmonized, pattern = "*.laz$", full.names = TRUE) 
process_missing_outputs(files_harmonized, path(folder_thinned, "ppm20"), function(files) thin_by_pulses(files, den = 20))
process_missing_outputs(files_harmonized, path(folder_thinned, "ppm10"), function(files) thin_by_pulses(files, den = 10))
process_missing_outputs(files_harmonized, path(folder_thinned, "ppm4"), function(files) thin_by_pulses(files, den = 4))


#### Index files
files_thinned <- list.files(folder_thinned, pattern = "*.laz$", full.names = TRUE, recursive = TRUE) 

# create spatial index for files if not present
files_not_indexed <- files_thinned[!lasR:::is_indexed(files_thinned)]
if (length(files_not_indexed) > 0) {
  n_files <- length(files_not_indexed)
  
  for (i in seq_along(files_not_indexed)) {
    file <- files_not_indexed[i]
    message(sprintf("Indexing file %d/%d: %s", i, n_files, basename(file)))
    exec(write_lax(embedded = TRUE), on = path(file))
  }
}


#### create VPCs if not present
if(!file.exists(vpc_leafon_harmonized)){
  list.files(folder_harmonized, pattern = "*2023.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafon_harmonized)
}
if(!file.exists(vpc_leafoff_harmonized)){
  list.files(folder_harmonized, pattern = "*2024.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafoff_harmonized)
}

if(!file.exists(vpc_leafon_ppm4)){
  list.files(path(folder_thinned, "ppm4"), pattern = "*2023.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafon_ppm4)
}

if(!file.exists(vpc_leafon_ppm10)){
  list.files(path(folder_thinned, "ppm10"), pattern = "*2023.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafon_ppm10)
}

if(!file.exists(vpc_leafon_ppm20)){
  list.files(path(folder_thinned, "ppm20"), pattern = "*2023.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafon_ppm20)
}


if(!file.exists(vpc_leafoff_ppm4)){
  list.files(path(folder_thinned, "ppm4"), pattern = "*2024.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafoff_ppm4)
}

if(!file.exists(vpc_leafoff_ppm10)){
  list.files(path(folder_thinned, "ppm10"), pattern = "*2024.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafoff_ppm10)
}

if(!file.exists(vpc_leafoff_ppm20)){
  list.files(path(folder_thinned, "ppm20"), pattern = "*2024.laz$", full.names = TRUE) |> 
    get_vpc(save_vpc_to = vpc_leafoff_ppm20)
}
