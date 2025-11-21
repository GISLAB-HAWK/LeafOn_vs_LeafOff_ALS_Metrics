#-------------------------------------------------------------------------------
# Name:         vol_sample_plots.R
# Description:  Calculation of forest attributes in inventory plots
#               (Betriebsionventur (BI) Lower Saxony).
#               Inventory data is first pre-processed and then single tree volumes
#               are calculated. These tree volumes are aggregated per sample plot
#               to obtain the growing stock volume (GSV) [m³/ha]. Other attributes
#               which are calculated per sample plot include tree density [n/ha], 
#               basal area [m³/ha], and quadratic mean diameter (QMD) [cm].
#               The final plots with the calculated forest attributes are clipped
#               to the area of interest (AOI), which is the area covered by both
#               leaf-off and leaf-on airborne laser scanning (ALS) datasets.
#               Some plots were were remeasured with RTK-GNSS. Where available,
#               the corrected coordinates of these plots are used.
#               Plots with a defined vegetation height change between leaf-on 
#               and leaf-off that cannot be attributed to seasonal differences,
#               but rather to treefall (harvest, natural disturbance), are removed.
# Author:       Christoph Fischer, Georgia Reeves, Florian Franz
# Contact:      christoph.fischer@nw-fva.de
#-------------------------------------------------------------------------------


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
bi_plots_rtk <- sf::st_read(file.path(bi_path, 'bi_center_points_rtk.gpkg'))

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



# 03: calculate individual tree volume
#-------------------------------------------------------------------------------

vor <- bi_points_trees

# recode tree species before applying treegross
vor$ba1 <- ifelse(vor$ba == 561 | vor$ba == 526 , 511,
           ifelse(vor$ba == 21 , 211,
           ifelse(vor$ba == 716, 711, vor$ba)))

# add empty columns
vor$vol <- NA
vor$volC <- NA

# create empty table
cop <- vor[-(1:dim(vor)[1]),]

# calculation of volumes with tg_volume
# correction of negative volumes with GAM
# to avoid them
d = seq(7, 99, by = 2)
h = seq(1, 59,  by = 2)
data = expand.grid(d = d, h = h)

for(i in unique(vor$ba1)){
  
  print(i)
  
  vor2 <- vor[vor$ba1 == i,]
  vor2$vol <- TreeGrOSSinR::tg_volumen(ba=i, bhd=vor2$bhd, h=vor2$hoe_mod, info = F)
  
  data$vol <- TreeGrOSSinR::tg_volumen(ba=i, bhd=data$d, h=data$h, info = F)
  m <-mgcv::gam(vol ~ t2(d, h, k = 10), data = data, family = gaussian(link = 'log'))
  
  vor2$volC <- predict(m, newdata=data.frame(d=vor2$bhd, h=vor2$hoe_mod), type = 'response')
  
  cop <- rbind(cop, vor2)
  
}

cop$vol <- ifelse(cop$vol <= 0, cop$volC, cop$vol)

vor <- cop
rm(m, data, cop, i, vor2)
vor$ba1 <- NULL
vor$id2 <- NULL
vor$volC <- NULL 

bi_points_trees <- vor
summary(bi_points_trees$vol)
plot(bi_points_trees$bhd, bi_points_trees$vol)



# 04: calculate individual tree AGB
#-------------------------------------------------------------------------------

agb <- bi_points_trees

# recode tree species before applying rBDAT
agb$ba1 <- ifelse(agb$ba == 110 | agb$ba == 111 | agb$ba == 112, 17,
           ifelse(agb$ba == 113, 18,
           ifelse(agb$ba == 211, 15,
           ifelse(agb$ba == 221, 16,
           ifelse(agb$ba == 311, 21,
           ifelse(agb$ba == 320, 22,
           ifelse(agb$ba == 321, 23,
           ifelse(agb$ba == 322, 24,
           ifelse(agb$ba == 323, 25,
           ifelse(agb$ba == 330 | agb$ba == 331 | agb$ba == 332, 30,
           ifelse(agb$ba == 341 | agb$ba == 342, 27,
           ifelse(agb$ba == 351, 31,
           ifelse(agb$ba == 352 | agb$ba == 442, 33,
           ifelse(agb$ba == 353 | agb$ba == 355, 35,
           ifelse(agb$ba == 354 | agb$ba == 452, 29,
           ifelse(agb$ba == 357, 32,
           ifelse(agb$ba == 410 | agb$ba == 411 | agb$ba == 412 | agb$ba == 414, 26,
           ifelse(agb$ba == 420 | agb$ba == 421 | agb$ba == 422, 28,
           ifelse(agb$ba == 430 | agb$ba == 431, 19,
           ifelse(agb$ba == 441, 34,
           ifelse(agb$ba == 451, 36,
           ifelse(agb$ba == 511 | agb$ba == 513 | agb$ba == 551, 1,
           ifelse(agb$ba == 512, 2,
           ifelse(agb$ba == 521, 3,
           ifelse(agb$ba == 523, 4,
           ifelse(agb$ba == 541, 13,
           ifelse(agb$ba == 542, 12,
           ifelse(agb$ba == 611, 8,
           ifelse(agb$ba == 711, 5,
           ifelse(agb$ba == 712, 6,
           ifelse(agb$ba == 731, 7,
           ifelse(agb$ba == 810, 9,
           ifelse(agb$ba == 811, 10,
           ifelse(agb$ba == 812, 11, vor$ba)
           )))))))))))))))))))))))))))))))))

# add empty column
agb$agb <- NA

# calculation of AGB with rBDAT::getBiomass
for(i in unique(agb$ba1)){
  
  print(i)
  
  agb$agb <- rBDAT::getBiomass(
    agb,
    mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
  )
  
}

agb$ba1 <- NULL
bi_points_trees <- agb
summary(agb)
plot(agb$bhd, agb$agb)



# 05: aggregate volume and AGB per sample plot
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
    vol_ha = sum(vol * nha),
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
  bi_points_trees[,c("key", "kspnr", "abt", "rw", "hw", "vol_ha", "agb_ha",
                     "tree_density", "basal_area_ha", "dg", "dominant_leaf_type")])

#vol_stp <- merge(bi_points[,c("key", "kspnr", "abt", "rw", "hw")], 
#                 vol_stp, by=c("key", "kspnr", "abt", "rw", "hw"), 
#                 all.x=T)

inv_attr_plots[is.na(inv_attr_plots)] <- 0

head(inv_attr_plots)
summary(inv_attr_plots)

boxplot(inv_attr_plots$vol_ha)
boxplot(inv_attr_plots$agb_ha)
boxplot(inv_attr_plots$tree_density)
boxplot(inv_attr_plots$basal_area_ha)
boxplot(inv_attr_plots$dg)

hist(inv_attr_plots$vol_ha)
hist(inv_attr_plots$agb_ha)


# 06: clip sample plots to the AOI
#-------------------------------------------------------------------------------

# conversion to sf object (DHDN / 3-degree Gauss-Kruger zone 3)
inv_attr_plots_gk <- sf::st_as_sf(
  inv_attr_plots, coords = c('rw', 'hw'), crs = 31467
  )

# transformation to ETRS89 / UTM zone 32N
inv_attr_plots_utm <- sf::st_transform(inv_attr_plots_gk, crs = 25832)

# quick plot
par(mfrow = c(1,2))
lidR::plot(pc_ctg_loff)
terra::plot(inv_attr_plots_utm$geom, col = 'red', add = T)
lidR::plot(pc_ctg_lon)
terra::plot(inv_attr_plots_utm$geom, col = 'red', add = T)

# clip BI plots to the area only covered by leaf-off point clouds
# leaf-off covers a slightly smaller area than leaf-on
inv_attr_plots_aoi <- sf::st_intersection(
  inv_attr_plots_utm, sf::st_as_sf(pc_ctg_loff)
  )

# keep only the original columns from inv_attr_plots_utm 
# (remove LAScatalog columns)
original_cols <- names(inv_attr_plots_utm)
inv_attr_plots_aoi <- inv_attr_plots_aoi[, original_cols]

summary(inv_attr_plots_aoi)
table(inv_attr_plots_aoi$dominant_leaf_type)
par(mfrow = c(3,2))
boxplot(inv_attr_plots_aoi$vol_ha)
boxplot(inv_attr_plots_aoi$agb_ha)
boxplot(inv_attr_plots_aoi$tree_density)
boxplot(inv_attr_plots_aoi$basal_area_ha)
boxplot(inv_attr_plots_aoi$dg)

# visualize locations of the BI plots
lidR::plot(pc_ctg_lon, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(inv_attr_plots_aoi, col.regions = 'black', cex = 5)

lidR::plot(pc_ctg_loff, mapview = T, 
           map.type = 'OpenStreetMap',
           alpha.regions = 0) +
  
  mapview::mapview(inv_attr_plots_aoi, col.regions = 'black', cex = 5)



# 06: include remeasured RTK-GNSS plots
#-------------------------------------------------------------------------------

# merge remeasured plots into inv_attr_plots_aoi
inv_attr_plots_aoi$remeasured <- 'no'

# identify matching plots based on kspnr column
# note: inv_attr_plots_aoi has "kspnr", bi_plots_rtk has "KSPNR"
matching_plots <- inv_attr_plots_aoi$kspnr %in% bi_plots_rtk$KSPNR

# mark remeasured plots
inv_attr_plots_aoi$remeasured[matching_plots] <- 'yes'

# for plots that were remeasured,
# update their geometry with the more accurate RTK positions
if (any(matching_plots)) {
  
  # create a temporary data frame for merging
  rtk_temp <- bi_plots_rtk[, c('KSPNR')]
  rtk_temp$rtk_geometry <- sf::st_geometry(bi_plots_rtk)
  
  # merge RTK geometry data
  inv_attr_plots_aoi_temp <- merge(
    inv_attr_plots_aoi, 
    sf::st_drop_geometry(rtk_temp), 
    by.x = 'kspnr', 
    by.y = 'KSPNR', 
    all.x = T
  )
  
  # update geometry for remeasured plots
  for (i in which(matching_plots)) {
    kspnr_val <- inv_attr_plots_aoi$kspnr[i]
    rtk_row <- which(bi_plots_rtk$KSPNR == kspnr_val)
    if (length(rtk_row) > 0) {
      sf::st_geometry(inv_attr_plots_aoi)[i] <- sf::st_geometry(bi_plots_rtk)[rtk_row[1]]
    }
  }
  
  cat('Updated', sum(matching_plots), 'plots with RTK-GNSS coordinates\n')
  
} else {
  
  cat('No matching plots found between inv_attr_plots_aoi and bi_plots_rtk\n')
  
}



# 07: remove plots based on height differences between leaf-off and leaf-on
#-------------------------------------------------------------------------------

extract_plot_heights_chm <- function(catalog, plots, res = 0.5, buffer_radius = 13) {
  
  # create buffered plots
  plots_buffered <- sf::st_buffer(plots, dist = buffer_radius)
  
  # generate CHM for the area containing all plots
  chm_opt <- list(res = res, algorithm = lidR::p2r())
  chm <- lidR::rasterize_canopy(catalog, chm_opt$res, chm_opt$algorithm)
  
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

# extract heights from CHMs of both point cloud catalogs
cat('Extracting heights from leaf-off point cloud...\n')
heights_loff <- extract_plot_heights_chm(pc_ctg_loff, inv_attr_plots_aoi)

cat("Extracting heights from leaf-on point cloud...\n") 
heights_lon <- extract_plot_heights_chm(pc_ctg_lon, inv_attr_plots_aoi)

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

# remove corresponding plots
inv_attr_plots_aoi_filtered <- inv_attr_plots_aoi[valid_plots, ]



# 07: save BI plots with the forest inventory attributes per sample plot
#-------------------------------------------------------------------------------

# rds
out_path <- file.path(processed_data_dir, 'forest_inventory')

if (!file.exists(file.path(out_path, 'inv_attr_plots.RDS'))) {
  
  inv_attr_plots_aoi_no_geom <- sf::st_drop_geometry(inv_attr_plots_aoi_filtered)
  saveRDS(
    inv_attr_plots_aoi_no_geom, 
    file = file.path(out_path, 'inv_attr_plots.RDS')
    )
  
} else {
  
  print('File inv_attr_plots.RDS already exists.')
  
}

# txt
if (!file.exists(file.path(out_path, 'inv_attr_plots.txt'))) {
  
  inv_attr_plots_aoi_no_geom <- sf::st_drop_geometry(inv_attr_plots_aoi_filtered)
  write.table(
    inv_attr_plots_aoi_no_geom, 
    file = file.path(out_path, 'inv_attr_plots.txt'), 
    sep = ';',
    row.names = F
  )
  
} else {
  
  print('File inv_attr_plots.txt already exists.')
  
}

# gpkg
if (!file.exists(file.path(out_path, 'inv_attr_plots.gpkg'))) {
  
  sf::st_write(
    inv_attr_plots_aoi_filtered,
    dsn = file.path(out_path, 'inv_attr_plots.gpkg')
    )
  
} else {
  
  print('File inv_attr_plots.gpkg already exists.')
  
}


