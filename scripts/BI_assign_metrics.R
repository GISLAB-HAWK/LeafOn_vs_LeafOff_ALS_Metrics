# load libraries
library(sf)
library(terra)
library(dplyr)

# set working directory 
setwd("C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/")


# import data 
BI <- read_sf("data/vol_stp_filtered.gpkg") # Betriebsinventur Punkte 

ind23 <- read.csv("data/BI_indices_2023_RTK_harm.csv") # Indices for BI points 
ind24 <- read.csv("data/BI_indices_2024_RTK_harm.csv")

names(ind23) <- c("kspnr", paste0("lon_", names(ind23)[-1])) # rename layer 
names(ind24) <- c("kspnr", paste0("loff_", names(ind24)[-1]))

# combine all the information together in one layer 
BI_join <- BI %>% 
  mutate(kspnr = as.integer(kspnr)) %>% 
  left_join(ind23, by = "kspnr") %>% 
  left_join(ind24, by = "kspnr") %>%
  relocate(geom, .after = last_col() )

# export layer 
write_sf(BI_join, "data/vol_stp_filtered_metrics.gpkg")

