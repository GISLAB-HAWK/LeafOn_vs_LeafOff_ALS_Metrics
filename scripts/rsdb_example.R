#--------------------------------------------------------------------------
# Name:         rsdb_example.R
# Description:  A script that demostrates how to use the HAWK Remote Sensing 
#               Database (RSDB) in R. To run the script you need to have access 
#               to RSDB.
#               For further information see
#               https://github.com/environmentalinformatics-marburg/rsdb-data and 
#               https://environmentalinformatics-marburg.github.io/rsdb/docs/r_package_installation/ 
# Author:       Svenja Dobelmann
#--------------------------------------------------------------------------

# install package 
# if(!require("remotes")) install.packages("remotes")
# remotes::install_github("environmentalinformatics-marburg/rsdb/r-package")

library(RSDB)
library(lidR)

# logging into the server
fileName <- "C:/rsdb/rsdb_login.txt" # file containing your username and pw in the form "username:password"
userpwd <- readChar(fileName, file.info(fileName)$size) #  read account from file
remotesensing <- RemoteSensing$new("https://gislab.hawk.de",userpwd)

# open web interface in browser (optional)
remotesensing$web()

## get names of available PointDBs
remotesensing$pointclouds

## load the pointcloud 
pointcloud <- remotesensing$pointcloud("solling_2023")

## metadata 
pointcloud$description
pointcloud$geocode
pointcloud$proj4

df <- pointcloud$points(extent)

las <- as.LAS(df)
plot(las)

##### TO BE CONTINUED!  #########################################################################

## define Tiling parameter
extent <- pointcloud$extent
tile_size= 1000
x_tiles <- seq(extent@xmin, extent@xmax, by = tile_size)
y_tiles <- seq(extent@ymin, extent@ymax, by = tile_size)
# Generate all tile combinations
tile_grid <- expand.grid(x = x_tiles[-length(x_tiles)], y = y_tiles[-length(y_tiles)])


### rasterize the tiles 
##  list to store rasters
raster_list <- list()

## loop over tiles 
for (i in c(1:nrow(tile_grid))){
  print(paste0("Processing tile:", i,"/",nrow(tile_grid), " at ", Sys.time()))
  tile <- tile_grid[i,]
  tile_extent <- raster::extent(tile$x, tile$x + tile_size, tile$y, tile$y + tile_size)
  
  tile_raster <-  pointcloud$raster(ext = tile_extent, res = 1, type = "dtm") # define type here! 
  
  # Store raster in the list with a unique name
  tile_name <- paste(tile$x, tile$y, sep = "_")
  raster_list[[tile_name]] <- tile_raster
}


# list available indices 
pointcloud$index_list$name

# calculate  index for first tile 
df <- pointcloud$indices(raster::extent(raster_list[[1]]), "BE_H_MAX")


# get raster of points per pixel at 1 meter resolution
r <- pointcloud$raster(ext, res=1,  type="point_count")
raster::plot(r)

# get raster of laser pulses per pixel at 1 meter resolution
r <- pointcloud$raster(ext, res=1,  type="pulse_count")
raster::plot(r)

# get raster of DSM at 0.5 meter resolution and fill missing pixels
r <- pointcloud$raster(ext, res=0.5,  type="dsm", fill=10)
terra::plot(r)
visualise_raster(r)

# get raster of DTM at 0.5 meter resolution and fill missing pixels
r <- pointcloud$raster(ext, res=0.5,  type="dtm", fill=10)
raster::plot(r)
visualise_raster(r)

# get raster of CHM at 0.5 meter resolution and fill missing pixels
r <- pointcloud$raster(ext, res=0.5,  type="chm", fill=10)
raster::plot(r)
visualise_raster(r)


# create polygon
library(sp)
x <- 597300.863 + c(1,2,1.5,1,1)
y <- 5671900.623 + c(1,1,2.7,2,1)
coords <- cbind(x, y)
p <- Polygon(coords=coords)
ps <- Polygons(list(p), ID=17)
sps <- SpatialPolygons(list(ps))
plot(sps)

# get points in polygon
df <- pointcloud$points(ps, columns=c("x", "y", "z", "classification", "returns"))

# calculate all indices for area of polygon
df <- pointcloud$indices(list(p1=ps), pointcloud$index_list$name)

# calcaulte for one polygon one index
pointcloud$indices(ps, "BE_H_MAX")

#get one ROI, create polygon and create bounding box
roi <- remotesensing$roi(group_name="hai", roi_name="HEW02")
p <- Polygon(coords=roi$polygon)
x <- p@coords[,1]
y <- p@coords[,2]
ext <- extent(x=min(x), xmax=max(x), ymin=min(y), ymax=max(y)) # create bounding box from polygon
pointcloud$indices(ext, "BE_H_MAX") # calculate index at ext
pointcloud$indices(p, "BE_H_MAX") # calculate index at polygon


