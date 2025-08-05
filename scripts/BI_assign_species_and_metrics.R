# load libraries
library(sf)
library(terra)
library(dplyr)

# set working directory 
setwd("~")


# import data 
BI <- read_sf("vol_stp.gpkg") # Betriebsinventur Punkte 
ts <- rast("treespecies_de_2022.tif") # Tree Species Data 

ind23 <- read.csv("BI_indices_2023.csv") # Indices for BI points 
ind24 <- read.csv("BI_indices_2024.csv")

names(ind23) <- c("kspnr", paste0("lon_", names(ind23)[-1])) # rename layer 
names(ind24) <- c("kspnr", paste0("loff_", names(ind24)[-1]))


# extract tree species for BI points 
crs(BI) == crs(ts) # check crs
BI_prj <- st_transform(BI, crs(ts)) # reproject BI points to mask them first 

ts_crop <- crop(ts, BI_prj) # crop to extent of BI points 
ts_utm <- project(ts_crop, crs(BI), method = "near") # reproject to EPSG:25832

dummy <- rast(extent = ext(ts_utm), resolution = 20, crs = crs(BI)) # create dummy raster with target resolution
ts_res <- resample(ts_utm,dummy, method = "near") # resample tree species to 20m 
#writeRaster(ts_res, "treespecies_de_2022_reprojected_resampled.tif")

BI_ts <- extract(ts_res,BI) # extract tree species per point


# combine all the information together in one layer 
BI_join <- BI %>% 
  mutate(ts = BI_ts[,2]) %>%
  mutate(kspnr = as.integer(kspnr)) %>% 
  left_join(ind23, by = "kspnr") %>% 
  left_join(ind24, by = "kspnr") %>%
  relocate(geom, .after = last_col() )

# export layer 
write_sf(BI_join, "./data/vol_stp_species_metrics.gpkg")

