# define tile of interest ("xxx_yyyy")
tile <- "548_5728"




source("scripts/functions.R")
source("config.R")
library(lidR)

tile_on <- paste0("3dm_32_", tile, "_1_ni_2023.laz")
tile_off <- paste0("3dm_32_", tile, "_1_ni_2024.laz")

file_leafon_raw <- path(folder_leafon_raw, tile_on)
file_leafon_classified <- path(folder_classified, tile_on)
file_leafon_harmonized <- path(folder_harmonized, tile_on)
file_leafon_ppm20 <- path(folder_thinned, "ppm20", tile_on)

file_leafoff_raw <- path(folder_leafoff_raw, tile_off)
file_leafoff_classified <- path(folder_classified, tile_off)
file_leafoff_harmonized <- path(folder_harmonized, tile_off)
file_leafoff_ppm20 <- path(folder_thinned, "ppm20", tile_off)


# function to read a subset (50m radius) at center of file 
read_subset <- function(file_path, r = 50, select = NULL) {
  
  if (endsWith(file_path, ".laz")) {
    # read sample subset in center
    # get center coordinate
    header <- lidR::readLASheader(file_path)
    xc <- header$`Min X` + (header$`Max X` - header$`Min X`) / 2
    yc <- header$`Min Y` + (header$`Max Y` - header$`Min Y`) / 2
    # read tiny subset
    if (!is.null(select)){
      ans <- lidR::readLAS(file_path, filter = paste("-keep_circle", xc, yc, r), select = select)
      
    } else {
      ans <- lidR::readLAS(file_path, filter = paste("-keep_circle", xc, yc, r))
    }
    return(ans)
    
  } else {
    warning("File type not supported: ", file_path)
  }
}


las_leafon_harmonized <- read_subset(file_leafon_harmonized)
# las_leafon_thinned <- read_subset(file_leafon_ppm4)
# las_leafon_thinned <- read_subset(file_leafon_ppm10)
las_leafon_thinned <- read_subset(file_leafon_ppm20)

las_leafoff_harmonized <- read_subset(file_leafoff_harmonized)
# las_leafoff_thinned <- read_subset(file_leafoff_ppm4)
# las_leafoff_thinned <- read_subset(file_leafoff_ppm10)
las_leafoff_thinned <- read_subset(file_leafoff_ppm20)


# plot thinned and un-thinned point cloud over each other
# leaf-on
rbind(add_attribute(las_leafon_thinned, 2, "source"), add_attribute(las_leafon_harmonized, 1, "source")) |> 
  plot(color = "source")
# leaf-off
rbind(add_attribute(las_leafoff_thinned, 2, "source"), add_attribute(las_leafoff_harmonized, 1, "source")) |> 
  plot(color = "source")

# plot leaf-on and leaf-off point cloud over each other
rbind(add_attribute(las_leafon_thinned, 2, "source"), add_attribute(las_leafoff_thinned, 1, "source")) |>
  plot(color = "source")


# function to plot point clouds next to each other
plot_sidebyside <- function(las1, las2, col_by = "Z"){
  rgl::open3d()
  rgl::mfrow3d(1, 2, sharedMouse = TRUE)
  rgl::plot3d(las1@data$X, las1@data$Y, las1@data$Z, col = viridis::cividis(100)[cut(las1@data[[col_by]], breaks = 100)], size = 2, decorate = FALSE)
  rgl::plot3d(las2@data$X, las2@data$Y, las2@data$Z, col = viridis::cividis(100)[cut(las2@data[[col_by]], breaks = 100)], size = 2, decorate = FALSE)
}

# leaf-off - leaf-on
plot_sidebyside(las_leafoff_thinned, las_leafon_thinned)
plot_sidebyside(las_leafoff_thinned, las_leafon_thinned, col_by = "Intensity")

# un-thinned - thinned
plot_sidebyside(las_leafoff_harmonized, las_leafoff_thinned)
plot_sidebyside(las_leafon_harmonized, las_leafon_thinned)



# function to plot pulse density
plot_pulsedensity <- function(las, resolution = 1, main = NULL) {
  pm <- pixel_metrics(
    las = las,
    func = ~list(pulse_count = length(unique(gpstime))),
    res = resolution,
    filter = ~ReturnNumber == NumberOfReturns
  )
  # Convert counts to density
  pulse_density <- pm / (resolution^2)
  terra::plot(pulse_density, main = main)
}


las_leafoff_harmonized |> plot_pulsedensity(resolution = 2, main = "Leaf-off")
las_leafoff_thinned |> plot_pulsedensity(resolution = 2, main = "Leaf-off thinned")
las_leafon_harmonized |> plot_pulsedensity(resolution = 2, main = "Leaf-on")
las_leafon_thinned |> plot_pulsedensity(resolution = 2, main = "Leaf-on thinned")


# plot height above ground of first and last returns
plot_hag_byreturn <- function(las, resolution = 1, main = NULL) {
  pm <- pixel_metrics(
    las = las,
    func = ~list(mean_hag = mean(HAG)),
    res = resolution,
    by_echo =c("first", "lastofmany")
  )
  terra::plot(pm, main = main, 
              plg = list(title = "height"), 
              legend = "first",
              range = c(min(terra::values(pm), na.rm = TRUE), max(terra::values(pm), na.rm = TRUE)))
}

las_leafoff_thinned |> plot_hag_byreturn(resolution = 2, main = c("Leaf-off - FirstReturns", "Leaf-off - LastReturns"))
las_leafon_thinned |> plot_hag_byreturn(resolution = 2, main = c("Leaf-on - FirstReturns", "Leaf-on - LastReturns"))

# print penetration depth (average HAG by return number)
las_leafoff_thinned |> print_penetration_depth()
las_leafon_thinned |> print_penetration_depth()
