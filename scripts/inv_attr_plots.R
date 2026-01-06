#----------------------------------------------------------------------------------------------------------------------------------------------
# Name:         inv_attr_plots.R
# Description:  Calculation of forest attributes in inventory plots (Betriebsionventur (BI) Lower Saxony).
#               Inventory data is first pre-processed and then volume (total and merchantable)
#               and above ground biomass (AGB) are calculated for individual trees. 
#               The tree volumes and AGB are aggregated per sample plot to obtain the growing stock volume (GSV) [m³/ha] and AGB [t/ha].
#               Other attributes which are calculated per sample plot include tree density [n/ha],
#               basal area [m³/ha], and quadratic mean diameter (QMD) [cm].
#               The plots with the calculated forest attributes are clipped to the area of interest (AOI),
#               which is the area covered by both leaf-off and leaf-on airborne laser scanning (ALS) datasets.
#               Some plots were were remeasured with RTK-GNSS. Where available, the corrected coordinates of these plots are used.
#               Plots with a defined vegetation height change between leaf-on and leaf-off that cannot be attributed to seasonal differences,
#               but rather to treefall (harvest, natural disturbance), are removed.
# Author:       Christoph Fischer, Georgia Reeves, Florian Franz
# Contact:      christoph.fischer@nw-fva.de
#               georgia.reeves@nw-fva.de
#               florian.franz@nw-fva.de
#----------------------------------------------------------------------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)



# 01: data reading
#-------------------------------------------------------------------------------

# input paths
bi_path <- file.path(raw_data_dir, 'forest_inventory')
pc_loff_path <- file.path(raw_data_dir, 'pc_leafoff_2024')
pc_lon_path <- file.path(raw_data_dir, 'pc_leafon_2023')

# read forest inventory (BI) data
bi_data <- list.files(bi_path)

bi_points <- read.table(
  file.path(bi_path, 'tblDatPh2_ZE.txt'),
  header = T, sep = ';'
)

bi_trees <- read.table(
  file.path(bi_path, 'tblDatPh2_Vorr_ZE.txt'),
  header = T, sep = ';'
)

# select desired forestry offices (Solling --> Neuhaus, Dassel)
bi_points <- bi_points[
  bi_points$DatOrga_Key == '268-2022-002' | 
    bi_points$DatOrga_Key == '254-2022-002',
]

bi_trees <- bi_trees[
  bi_trees$DatOrga_Key == '268-2022-002' |
    bi_trees$DatOrga_Key == '254-2022-002',
]

head(bi_points)
head(bi_trees)

# read remeasured plots (RTK-GNSS)
bi_plots_rtk_1 <- sf::st_read(file.path(bi_path, 'bi_center_points_rtk.gpkg'))
bi_plots_rtk_2 <- sf::st_read(file.path(bi_path, 'remeasured_plots.shp'))

# clean up bi_plots_rtk_1: 
# - create kspnr column from KSPNR
# - rename geometry column to 'geometry' for consistency
bi_plots_rtk_1$kspnr <- bi_plots_rtk_1$KSPNR
sf::st_geometry(bi_plots_rtk_1) <- 'geometry'

# clean up bi_plots_rtk_2:
# - remove trailing letters (e.g. "e" or "b") from Name column to get kspnr
# - drop Z dimension from geometry
bi_plots_rtk_2$kspnr <- as.integer(gsub('[a-zA-Z]$', '', bi_plots_rtk_2$Name))
bi_plots_rtk_2 <- sf::st_zm(bi_plots_rtk_2, drop = T, what = 'ZM')

# merge both RTK datasets, keeping only geometry and kspnr
bi_plots_rtk <- rbind(
  bi_plots_rtk_1[, 'kspnr'],
  bi_plots_rtk_2[, 'kspnr']
)

head(bi_plots_rtk)

# read point clouds with LAScatalog
pc_ctg_loff <- lidR::readLAScatalog(pc_loff_path)
pc_ctg_lon <- lidR::readLAScatalog(pc_lon_path)

pc_ctg_loff
pc_ctg_lon



# 02: data preparation
#-------------------------------------------------------------------------------

# source and apply function for data formatting
source('src/format_data.R', local = TRUE)

bi_points <- format_data(bi_points)
bi_trees <- format_data(bi_trees)

head(bi_points)
str(bi_points)
head(bi_trees)
str(bi_trees)

## delete deadwood and used trees
bi_trees <- bi_trees[!bi_trees$ba %in% seq(100,800,100),]

bi_trees <- bi_trees[bi_trees$'1' < 3 & bi_trees$'2' < 3,]

bi_trees <- bi_trees[bi_trees$art != 1 & bi_trees$art != 2 & bi_trees$bhd > 0,]

# select needed columns
bi_trees <- bi_trees[,c(1:12)]

# assign tree species groups
bi_trees$bagr <- 
  ifelse(bi_trees$ba > 0 & bi_trees$ba < 200, "EI",
  ifelse(bi_trees$ba > 199 & bi_trees$ba < 300, "BU",
  ifelse(bi_trees$ba > 299 & bi_trees$ba < 400, "ALH",	
  ifelse(bi_trees$ba > 399 & bi_trees$ba < 500, "ALN",
  ifelse(bi_trees$ba > 499 & bi_trees$ba < 600, "FI",
  ifelse(bi_trees$ba > 599 & bi_trees$ba < 700, "DGL",
  ifelse(bi_trees$ba > 699 & bi_trees$ba < 800, "KI",	"LAE")
  ))))))

# DBH correction
# if not measured at 1.3 m (deviating measuring height),
# then correction to 1.3 m
bi_trees$ba1 <- bi_trees$ba

# red oak to oak, fir to spruce, hornbeam to beech
source('src/d_corr_func.R', local = TRUE)

bi_trees$ba <- input_d_korr(bi_trees$bagr)

# average DBH from 'Kreuzkluppung (bhdklup)',
# convert to cm
bi_trees$bhd <- ifelse(
  bi_trees$bhdklup > 0,
  (0.5 * (bi_trees$bhd + bi_trees$bhdklup)) / 10,
  bi_trees$bhd / 10
)

# separate data:
# trees with diameter at deviating measurement height
# and without deviating measurement height
bi_trees_2 <- bi_trees[bi_trees$bhddiff > 0,]
bi_trees <- bi_trees[bi_trees$bhddiff == 0,]

# convert diameter to DBH in case of deviating measuring height
d <- d_korr(du = bi_trees_2$bhd, abwmh = bi_trees_2$bhddiff, ba = bi_trees_2$ba)
bi_trees_2$bhd <- d

# bind the two data frames together again
# (correctly back to original number)
bi_trees <- rbind(bi_trees, bi_trees_2)
rm(bi_trees_2, d)

# add original tree species again
bi_trees$ba <- bi_trees$ba1
bi_trees$ba1 <- NULL

head(bi_trees)

# calculate number of stems per ha

# concentric sample circles:
#	r = 6 m all trees 
#	r = 13 m all trees with DBH >= 30 cm
# radius must be projected into the plane

# correct sample circle sizes with slope, calculate N_ha
# inclination from degrees in rad
bi_points$hang_rad <- (pi / 180) * bi_points$hang

bi_points_trees <- merge(
  bi_trees, bi_points[,c("key", "kspnr", "abt", "hang_rad", "rw", "hw")],
  by = c("key", "kspnr")
)

# r_plane = r_slope * cos(slope_rad)
bi_points_trees$nha <- ifelse(
  bi_points_trees$bhd < 30, 
  10000 / (pi * 6**2 * cos(bi_points_trees$hang_rad)),
  10000 / (pi * 13**2 * cos(bi_points_trees$hang_rad))
)

## add heights

# heights in m
bi_points_trees$hoehe <- bi_points_trees$hoehe / 10

# new ID consisting of key + sample point number
bi_points_trees$id2 <- paste(bi_points_trees$key, bi_points_trees$kspnr, sep = "_")

# data format for heights adding
source('src/ehk_func.R', local = TRUE)

dat <- input_ehk(
  id = bi_points_trees$id2, bnr = bi_points_trees$id,
  bs = bi_points_trees$bestschicht, bhd = bi_points_trees$bhd,
  hoe = bi_points_trees$hoehe, nha = bi_points_trees$nha,
  bagr = bi_points_trees$bagr
)

head(dat)

# assigning appropriate quantile percentages in order
# to properly remove outliers (unrealistic DBH height value pairs)
# in a statistically sound way by building a scam model with height and DBH
dat <- dat[dat$hoe > 0, ]
summary(dat)

m <- scam::scam(
  hoe ~ s(bhd, bs = 'mpi'),
  data = dat, 
  family = Gamma(link = 'log')
)

nd <- data.frame('bhd' = floor(min(dat$bhd)):ceiling(max(dat$bhd)))
nd$hoe <- predict(m, newdata = nd, type = 'response')

p <- ggplot(data = dat, aes(x = bhd, y  = hoe)) + 
  geom_point(color = rgb(.5, .5, .5, alpha = .2)) + 
  geom_line(dat = nd, color = 1, linewidth = 2)

tmp <- NULL
for (x in seq(10, 110, by = 10)) {
  nd2 <- data.frame('bhd' = x)
  nd2$hoe <- predict(m, newdata = nd2, type = 'response')
  v <- 1/m$sig2
  d <- stats::dgamma(1:60, shape = (nd2$hoe[1]^2)/v, scale = v/nd2$hoe[1])
  tmp <- rbind(
    tmp, 
    data.frame(
      'bhd' = x - (9 * d / max(d)), 
      'hoe' = 1:60, 
      'x' = x
    )
  )
}

p1 <- p + 
  geom_path(dat = tmp, aes(group = x, color = factor(x)), show.legend = F) +
  geom_vline(
    data = data.frame('bhd' = seq(10, 110, by = 10)), 
    aes(xintercept = bhd, color = factor(bhd)), show.legend = F, 
    linetype = 2
  )

v <- 1/m$sig2
dat$p <- stats::pgamma(q = dat$hoe, shape = (fitted(m)^2)/v, scale = v/fitted(m))
summary(dat$p)

dat$lab <- ''
ix_label <- sort(c(which(dat$p > .99), which(dat$p < .01)))
dat$lab[ix_label] <- paste0(round(dat$p[ix_label] * 100, 2), '%')

p2 <- ggplot(data = dat, aes(x = bhd,y  = hoe)) + 
  geom_point(color = rgb(.5, .5, .5, alpha = .2)) + 
  geom_line(dat = nd, color = 1, linewidth = 2) + 
  ggrepel::geom_text_repel(aes(label = lab))

dat$lab <- round(dat$p*100, 2)
dat_extremes <- dat[dat$lab >= 99.98 | dat$lab <= 0.02,]
dat_without_extremes <- dat[dat$lab < 99.98 & dat$lab > 0.02,]

cowplot::plot_grid(p1, p2, ncol = 2)

# re-assigning dat to its original value
dat <- input_ehk(
  id = bi_points_trees$id2, bnr = bi_points_trees$id,
  bs = bi_points_trees$bestschicht, bhd = bi_points_trees$bhd,
  hoe = bi_points_trees$hoehe, nha = bi_points_trees$nha,
  bagr = bi_points_trees$bagr
)

dat2 <- dat[dat$hoe == 0,]

dat3 <- subset(dat_without_extremes, select = -c(p, lab))
dat <- rbind(dat2, dat3)
rm(dat_without_extremes)
rm(dat2)
rm(dat3)
rm(dat_extremes)

# uniform height curve
dat2 <- ehk(dat)
plot(dat2$bhd, dat2$hoe_mod)

# merge modeled heights with original table
bi_points_trees <- merge(
  bi_points_trees,
  dat2[,c("id", "bnr", "hoe_mod")],
  by.x = c("id2","id"), 
  by.y = c("id", "bnr")
)

# remove unneeded data frames
rm(dat2, dat)



# 03: calculate individual tree volume and AGB
#-------------------------------------------------------------------------------

vol_agb <- bi_points_trees

# recode tree species before applying rBDAT
# mapping from BI tree species codes to rBDAT species codes
vol_agb$ba1 <- ifelse(vol_agb$ba == 110 | vol_agb$ba == 111 | vol_agb$ba == 112, 17,
               ifelse(vol_agb$ba == 113, 18,
               ifelse(vol_agb$ba == 211, 15,
               ifelse(vol_agb$ba == 221, 16,
               ifelse(vol_agb$ba == 311, 21,
               ifelse(vol_agb$ba == 320, 22,
               ifelse(vol_agb$ba == 321, 23,
               ifelse(vol_agb$ba == 322, 24,
               ifelse(vol_agb$ba == 323, 25,
               ifelse(vol_agb$ba == 330 | vol_agb$ba == 331 | vol_agb$ba == 332, 30,
               ifelse(vol_agb$ba == 341 | vol_agb$ba == 342, 27,
               ifelse(vol_agb$ba == 351, 31,
               ifelse(vol_agb$ba == 352 | vol_agb$ba == 442, 33,
               ifelse(vol_agb$ba == 353 | vol_agb$ba == 355, 35,
               ifelse(vol_agb$ba == 354 | vol_agb$ba == 452, 29,
               ifelse(vol_agb$ba == 357, 32,
               ifelse(vol_agb$ba == 410 | vol_agb$ba == 411 | vol_agb$ba == 412 | vol_agb$ba == 414, 26,
               ifelse(vol_agb$ba == 420 | vol_agb$ba == 421 | vol_agb$ba == 422, 28,
               ifelse(vol_agb$ba == 430 | vol_agb$ba == 431, 19,
               ifelse(vol_agb$ba == 441, 34,
               ifelse(vol_agb$ba == 451, 36,
               ifelse(vol_agb$ba == 511 | vol_agb$ba == 513 | vol_agb$ba == 551, 1,
               ifelse(vol_agb$ba == 512, 2,
               ifelse(vol_agb$ba == 521, 3,
               ifelse(vol_agb$ba == 523, 4,
               ifelse(vol_agb$ba == 541, 13,
               ifelse(vol_agb$ba == 542, 12,
               ifelse(vol_agb$ba == 611, 8,
               ifelse(vol_agb$ba == 711, 5,
               ifelse(vol_agb$ba == 712, 6,
               ifelse(vol_agb$ba == 731, 7,
               ifelse(vol_agb$ba == 810, 9,
               ifelse(vol_agb$ba == 811, 10,
               ifelse(vol_agb$ba == 812, 11, vol_agb$ba)
               )))))))))))))))))))))))))))))))))

# calculation of volume with rBDAT::getVolume
# https://gitlab.com/vochr/rbdat
print('Calculating total volume...')
vol_agb$total_vol <- rBDAT::getVolume(
  vol_agb,
  bark = T,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

print('Calculating merchantable volume...')
vol_agb$merch_vol <- rBDAT::getVolume(
  vol_agb,
  bark = F,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

# calculation of AGB with rBDAT::getBiomass
print('Calculating AGB...')
vol_agb$agb <- rBDAT::getBiomass(
  vol_agb,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

vol_agb$ba1 <- NULL

# add results to bi_points_trees
bi_points_trees$total_vol <- vol_agb$total_vol
bi_points_trees$merch_vol <- vol_agb$merch_vol
bi_points_trees$agb <- vol_agb$agb

summary(bi_points_trees$total_vol)
summary(bi_points_trees$merch_vol)
summary(bi_points_trees$agb)

par(mfrow = c(1,3))
plot(bi_points_trees$bhd, bi_points_trees$total_vol)
plot(bi_points_trees$bhd, bi_points_trees$merch_vol)
plot(bi_points_trees$bhd, bi_points_trees$agb)
par(mfrow = c(1,1))

rm(vol_agb)



# 04: aggregate volume and AGB per sample plot
#-------------------------------------------------------------------------------

# add column of leaf type
unique(bi_points_trees$bagr)
bi_points_trees$leaf_type <- ifelse(
  bi_points_trees$bagr %in% c('EI', 'ALN', 'BU', 'ALH'),
  'deciduous',
  'coniferous'
)

# group sums of volume, AGB,
# and other forest inventory attributes like
# tree density, basal area, and QMD
bi_points_trees <- bi_points_trees %>%
  dplyr::group_by(key, kspnr) %>%
  dplyr::mutate(
    total_vol_ha = sum(total_vol * nha),
    merch_vol_ha = sum(merch_vol * nha),
    agb_ha = sum(agb * nha) / 1000,
    tree_density = mean(nha),
    basal_area_tree = (pi / 4) * (bhd / 100)^2,
    basal_area_ha = sum(basal_area_tree * nha, na.rm = T),
    dg = sqrt(sum(bhd^2 * nha, na.rm = T) / sum(nha, na.rm = T)),
    # assign dominant leaf type to each plot
    total_deciduous = sum(dplyr::if_else(
      leaf_type == 'deciduous', nha, 0, missing = 0), na.rm = T),
    total_coniferous = sum(dplyr::if_else(
      leaf_type == 'coniferous', nha, 0, missing = 0), na.rm = T),
    dominant_leaf_type = dplyr::case_when(
      total_deciduous > total_coniferous ~ 'deciduous',
      total_coniferous > total_deciduous ~ 'coniferous',
      TRUE                               ~ 'tie'
    )) %>%
  dplyr::ungroup()

# extract unique forest inventory variables for all sample plots
inv_attr_plots <- unique(
  bi_points_trees[,c(
    'key', 'kspnr', 'abt', 'rw', 'hw', 'total_vol_ha', 'merch_vol_ha',
    'agb_ha', 'tree_density', 'basal_area_ha', 'dg', 'dominant_leaf_type'
    )]
  )

inv_attr_plots[is.na(inv_attr_plots)] <- 0

head(inv_attr_plots)
summary(inv_attr_plots)

par(mfrow = c(2,3))
boxplot(inv_attr_plots$total_vol_ha)
boxplot(inv_attr_plots$merch_vol_ha)
boxplot(inv_attr_plots$agb_ha)
boxplot(inv_attr_plots$tree_density)
boxplot(inv_attr_plots$basal_area_ha)
boxplot(inv_attr_plots$dg)
par(mfrow = c(1,1))

par(mfrow = c(2,3))
hist(inv_attr_plots$total_vol_ha)
hist(inv_attr_plots$merch_vol_ha)
hist(inv_attr_plots$agb_ha)
hist(inv_attr_plots$tree_density)
hist(inv_attr_plots$basal_area_ha)
hist(inv_attr_plots$dg)
par(mfrow = c(1,1))



# 05: include remeasured RTK-GNSS plots
#-------------------------------------------------------------------------------

# conversion to sf object (DHDN / 3-degree Gauss-Kruger zone 3)
inv_attr_plots_gk <- sf::st_as_sf(
  inv_attr_plots, coords = c('rw', 'hw'), crs = 31467
  )

# transformation to ETRS89 / UTM zone 32N
inv_attr_plots_utm <- sf::st_transform(inv_attr_plots_gk, crs = 25832)

# merge remeasured plots into inv_attr_plots_utm
inv_attr_plots_utm$remeasured <- 'no'

# identify matching plots based on kspnr column
matching_plots <- inv_attr_plots_utm$kspnr %in% bi_plots_rtk$kspnr

# mark remeasured plots
inv_attr_plots_utm$remeasured[matching_plots] <- 'yes'

# create sf object with non-RTK geometries (without RTK position)
remeasured_plots_non_rtk <- inv_attr_plots_utm[matching_plots, ]

# for plots that were remeasured,
# update their geometry with the more accurate RTK positions
if (any(matching_plots)) {
  
  # update geometry for remeasured plots
  for (i in which(matching_plots)) {
    kspnr_val <- inv_attr_plots_utm$kspnr[i]
    rtk_row <- which(bi_plots_rtk$kspnr == kspnr_val)
    if (length(rtk_row) > 0) {
      sf::st_geometry(inv_attr_plots_utm)[i] <- sf::st_geometry(bi_plots_rtk)[rtk_row[1]]
    }
  }
  
  # create sf object with RTK geometries (after RTK update)
  remeasured_plots_rtk <- inv_attr_plots_utm[matching_plots, ]
  
  cat('Updated', sum(matching_plots), 'plots with RTK-GNSS coordinates\n')
  cat('Created remeasured_plots_non_rtk:', nrow(remeasured_plots_non_rtk), 'plots with original geometries\n')
  cat('Created remeasured_plots_rtk:', nrow(remeasured_plots_rtk), 'plots with RTK geometries\n')
  
} else {
  
  cat('No matching plots found between inv_attr_plots_utm and bi_plots_rtk\n')
  remeasured_plots_non_rtk <- NULL
  remeasured_plots_rtk <- NULL
  
}



# 06: clip sample plots to the AOI
#-------------------------------------------------------------------------------

# quick plot
par(mfrow = c(1,2))
lidR::plot(pc_ctg_loff)
terra::plot(inv_attr_plots_utm$geometry, col = 'red', add = T)
lidR::plot(pc_ctg_lon)
terra::plot(inv_attr_plots_utm$geometry, col = 'red', add = T)
par(mfrow = c(1,1))

# clip BI plots to the area only covered by leaf-off point clouds
# leaf-off covers a slightly smaller area than leaf-on
aoi <- sf::st_as_sf(pc_ctg_loff)
inv_attr_plots_aoi <- sf::st_intersection(inv_attr_plots_utm, aoi)

# keep only the original columns from inv_attr_plots_utm 
# (remove LAScatalog columns)
original_cols <- names(inv_attr_plots_utm)
inv_attr_plots_aoi <- inv_attr_plots_aoi[, original_cols]

# clip remeasured plots (non-RTK and RTK) to AOI
original_cols_remeasured <- names(remeasured_plots_non_rtk)
remeasured_plots_non_rtk_aoi <- sf::st_intersection(remeasured_plots_non_rtk, aoi)
remeasured_plots_non_rtk_aoi <- remeasured_plots_non_rtk_aoi[, original_cols_remeasured]
remeasured_plots_rtk_aoi <- sf::st_intersection(remeasured_plots_rtk, aoi)
remeasured_plots_rtk_aoi <- remeasured_plots_rtk_aoi[, original_cols_remeasured]

cat('Remeasured plots (non-RTK) in AOI:', nrow(remeasured_plots_non_rtk_aoi), '\n')
cat('Remeasured plots (RTK) in AOI:', nrow(remeasured_plots_rtk_aoi), '\n')

summary(inv_attr_plots_aoi)
table(inv_attr_plots_aoi$dominant_leaf_type)
par(mfrow = c(2,3))
boxplot(inv_attr_plots_aoi$total_vol_ha)
boxplot(inv_attr_plots_aoi$merch_vol_ha)
boxplot(inv_attr_plots_aoi$agb_ha)
boxplot(inv_attr_plots_aoi$tree_density)
boxplot(inv_attr_plots_aoi$basal_area_ha)
boxplot(inv_attr_plots_aoi$dg)
par(mfrow = c(1,1))

# visualize locations of the BI plots
lidR::plot(pc_ctg_lon, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(inv_attr_plots_aoi, col.regions = 'black', cex = 5)

lidR::plot(pc_ctg_loff, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(inv_attr_plots_aoi, col.regions = 'black', cex = 5)



# 07: remove plots based on height differences between leaf-off and leaf-on
#-------------------------------------------------------------------------------

# define path for storing excluded plots
excluded_plots_file <- file.path(
  processed_data_dir, 'forest_inventory', 'excluded_plots_height_diff.csv'
)

# check if we already have the list of plots to exclude
if (file.exists(excluded_plots_file)) {
  
  cat('Loading previously identified plots to exclude...\n')
  plots_to_exclude_df <- read.csv(excluded_plots_file, stringsAsFactors = F)
  plots_to_exclude <- plots_to_exclude_df$kspnr
  cat('Loaded', length(plots_to_exclude), 'plots to exclude\n')
  
} else {
  
  cat('Running CHM analysis to identify plots to exclude...\n')
  
  extract_plot_heights_chm <- function(
    pc_path,
    plots,
    res = 0.5, 
    buffer_radius = 13)
    {
    
    # create buffered plots
    plots_buffered <- sf::st_buffer(plots, dist = buffer_radius)
    
    # get all LAZ files from the path
    laz_files <- list.files(pc_path, pattern = '\\.laz$', full.names = T)
    
    # process each file individually and collect CHMs
    chm_list <- list()
    
    for (i in seq_along(laz_files)) {
      
      cat(sprintf('  Processing file %d/%d\n', i, length(laz_files)))
      
      # create fresh pipeline for each file
      chm_stage <- lasR::rasterize(res = res, operators = max(HAG))
      na_fill <- lasR::focal(chm_stage, size = 3, fun = 'mean')
      pipeline <- chm_stage + na_fill
      
      # execute pipeline on single file
      ans <- lasR::exec(
        pipeline, 
        on = laz_files[i],
        with = list(ncores = lasR::half_cores(), progress = T)
      )
      
      # store the NA-filled CHM (second element)
      chm_list[[i]] <- ans[[2]]
    }
    
    # merge all CHMs into one
    if (length(chm_list) == 1) {
      chm <- chm_list[[1]]
    } else {
      chm <- do.call(terra::merge, chm_list)
    }
    
    # extract mean height within each plot
    plot_heights <- exactextractr::exact_extract(
      chm,
      plots_buffered,
      fun = 'mean'
    )
    
    # convert list to vector
    if (is.list(plot_heights)) {
      plot_heights <- unlist(plot_heights)
    }
    
    return(plot_heights)
  }
  
  # extract heights from CHMs of both point cloud catalogs (leaf-off/-on)
  cat('Extracting heights from leaf-off point cloud...\n')
  heights_loff <- extract_plot_heights_chm(pc_loff_path, inv_attr_plots_aoi)
  
  cat('Extracting heights from leaf-on point cloud...\n')
  heights_lon <- extract_plot_heights_chm(pc_lon_path, inv_attr_plots_aoi)
  
  # add heights to plots
  inv_attr_plots_aoi$height_loff <- heights_loff
  inv_attr_plots_aoi$height_lon <- heights_lon
  
  # calculate height difference (leaf-off - leaf-on)
  inv_attr_plots_aoi$height_diff <- 
    inv_attr_plots_aoi$height_loff - inv_attr_plots_aoi$height_lon
  
  # create filter based on height difference threshold
  height_diff_threshold <- -10
  valid_plots <- is.na(inv_attr_plots_aoi$height_diff) |
    inv_attr_plots_aoi$height_diff > height_diff_threshold
  table(valid_plots)
  
  # store the kspnr of plots to exclude
  plots_to_exclude <- inv_attr_plots_aoi$kspnr[!valid_plots]
  
  # save for future use
  plots_to_exclude_df <- data.frame(kspnr = plots_to_exclude)
  write.csv(plots_to_exclude_df, excluded_plots_file, row.names = F)
  cat('Saved', length(plots_to_exclude), 'plots to exclude to', excluded_plots_file, '\n')
  
  # save non-filtered dataset
  sf::st_write(
    inv_attr_plots_aoi,
    file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots_non_filtered.gpkg')
  )
  
}

# apply the exclusion filter
valid_plots <- !inv_attr_plots_aoi$kspnr %in% plots_to_exclude

# filter the main dataset
inv_attr_plots_aoi_filtered <- inv_attr_plots_aoi[valid_plots, ]

# filter remeasured plots
if (!is.null(remeasured_plots_non_rtk_aoi) && !is.null(remeasured_plots_rtk_aoi)) {
  remeasured_plots_non_rtk_aoi_filtered <- remeasured_plots_non_rtk_aoi[
    !remeasured_plots_non_rtk_aoi$kspnr %in% plots_to_exclude,
  ]
  remeasured_plots_rtk_aoi_filtered <- remeasured_plots_rtk_aoi[
    !remeasured_plots_rtk_aoi$kspnr %in% plots_to_exclude,
  ]
  
  cat('Filtered datasets:\n')
  cat('Main plots:', nrow(inv_attr_plots_aoi_filtered), '/', nrow(inv_attr_plots_aoi), '\n')
  cat('Remeasured (non-RTK):', nrow(remeasured_plots_non_rtk_aoi_filtered), '/', nrow(remeasured_plots_non_rtk_aoi), '\n')
  cat('Remeasured (RTK):', nrow(remeasured_plots_rtk_aoi_filtered), '/', nrow(remeasured_plots_rtk_aoi), '\n')
}

# summary statistics
summary(inv_attr_plots_aoi_filtered)
summary_df <- as.data.frame(
  do.call(cbind, lapply(inv_attr_plots_aoi_filtered, summary))
)
write.csv(
  summary_df, 
  file.path(
    processed_data_dir, 
    'forest_inventory', 
    'summary_stats_inv_attr_plots_aoi_filtered.csv')
)
table(inv_attr_plots_aoi_filtered$dominant_leaf_type)
table_df <- as.data.frame(table(inv_attr_plots_aoi_filtered$dominant_leaf_type))
write.csv(
  table_df,
  file.path(
    processed_data_dir, 
    'forest_inventory',
    'n_plots_dom_leaf_type.csv'), 
    row.names = F
  )

# boxplots
par(mfrow = c(2,3))
boxplot(inv_attr_plots_aoi_filtered$total_vol_ha)
boxplot(inv_attr_plots_aoi_filtered$merch_vol_ha)
boxplot(inv_attr_plots_aoi_filtered$agb_ha)
boxplot(inv_attr_plots_aoi_filtered$tree_density)
boxplot(inv_attr_plots_aoi_filtered$basal_area_ha)
boxplot(inv_attr_plots_aoi_filtered$dg)
par(mfrow = c(1,1))



# 08: save BI plots with the forest inventory attributes per sample plot
#-------------------------------------------------------------------------------

out_path <- file.path(processed_data_dir, 'forest_inventory')

# define file formats and corresponding save functions
file_formats <- list(
  rds = list(
    ext = '.RDS',
    save_func = saveRDS
  ),
  txt = list(
    ext = '.txt',
    save_func = function(data, file) 
      write.table(data, file, sep = '\t', row.names = F)
  ),
  gpkg = list(
    ext = '.gpkg',
    save_func = function(data, file) 
      sf::st_write(data, file, delete_dsn = T, quiet = T)
  )
)

# define datasets to save
datasets <- list(
  list(
    name = 'inv_attr_plots',
    data = inv_attr_plots_aoi_filtered,
    drop_geom = T
  ),
  list(
    name = 'inv_attr_plots_non_rtk',
    data = remeasured_plots_non_rtk_aoi_filtered,
    drop_geom = T
  ),
  list(
    name = 'inv_attr_plots_rtk',
    data = remeasured_plots_rtk_aoi_filtered,
    drop_geom = T
  )
)

# loop through each file format
for (format_name in names(file_formats)) {
  
  format_info <- file_formats[[format_name]]
  
  # check if all files exist for this format
  all_files_exist <- all(sapply(datasets, function(ds) {
    file.exists(file.path(out_path, paste0(ds$name, format_info$ext)))
  }))
  
  if (!all_files_exist) {
    
    cat('Saving files in', format_name, 'format...\n')
    
    # save each dataset
    for (ds in datasets) {
      
      file_path <- file.path(out_path, paste0(ds$name, format_info$ext))
      
      # prepare data based on format
      if (format_name %in% c('rds', 'txt') && ds$drop_geom) {
        data_to_save <- sf::st_drop_geometry(ds$data)
      } else {
        data_to_save <- ds$data
      }
      
      # save using appropriate function
      format_info$save_func(data_to_save, file_path)
      cat('  Saved:', basename(file_path), '\n')
      
    }
    
  } else {
    
    cat('All', format_name, 'files already exist. Skipping.\n')
    
  }
  
}
