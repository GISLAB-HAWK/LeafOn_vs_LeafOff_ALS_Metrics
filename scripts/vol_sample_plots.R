#-------------------------------------------------------------------------------
# Name:         vol_sample_plots.R
# Description:  Calculation of timber volume per ha in forest inventory plots
#               (Betriebsionventur (BI) Lower Saxony).
#               Inventory data is first pre-processed and then single tree volumes
#               are calculated. These tree volumes are aggregated per sample plot
#               to obtain the growing stock [m³/ha].
# Author:       Christoph Fischer, Georgia Reeves, Florian Franz
# Contact:      christoph.fischer@nw-fva.de
#-------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)



# 01: data reading
#-------------------------------------------------------------------------------

# input path to forest inventory data (BI)
bi_path <- file.path(raw_data_dir, 'forest_inventory')

# read BI data
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



# 02: data preparation
#-------------------------------------------------------------------------------

## source and apply function for data formatting ----
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
  bi_trees, bi_points[,c("key", "kspnr","hang_rad", "rw", "hw")],
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
  bi_points_trees,#
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



# 04: aggregate volume per sample plot
#-------------------------------------------------------------------------------

# group sums of volume
bi_points_trees <- bi_points_trees %>%
  dplyr::group_by(key, kspnr) %>%
  dplyr::mutate(vol_ha = sum(vol * nha)) %>%
  dplyr::ungroup()

# extract unique volumes for all sample points
vol_stp <- unique(bi_points_trees[,c("key", "kspnr", "rw", "hw", "vol_ha")])

vol_stp <- merge(bi_points[,c("key", "kspnr", "rw", "hw")], 
                 vol_stp, by=c("key", "kspnr", "rw", "hw"), 
                 all.x=T)

vol_stp[is.na(vol_stp)] <- 0

head(vol_stp)
summary(vol_stp)
boxplot(vol_stp$vol_ha)

# creation of an sf object and transformation to UTM

# conversion to sf object (DHDN / 3-degree Gauss-Kruger zone 3)
vol_stp_gk <- sf::st_as_sf(vol_stp, coords = c('rw', 'hw'), crs = 31467)

# transformation to ETRS89 / UTM zone 32N
vol_stp_utm <- sf::st_transform(vol_stp_gk, crs = 25832)



# 05: save data with the extracted volumes per sample plot
#-------------------------------------------------------------------------------

# rds
out_path <- file.path(processed_data_dir, 'forest_inventory')

if (!file.exists(file.path(out_path, 'vol_stp.RDS'))) {
  
  saveRDS(vol_stp, file = file.path(out_path, 'vol_stp.RDS'))
  
} else {
  
  print('File vol_stp.RDS already exists.')
  
}

# txt
if (!file.exists(file.path(out_path, 'vol_stp.txt'))) {
  
  write.table(
    vol_stp, 
    file = file.path(out_path, 'vol_stp.txt'), 
    sep = ';',
    row.names = F
    )
  
} else {
  
  print('File vol_stp.txt already exists.')
  
}

# gpkg
if (!file.exists(file.path(out_path, 'vol_stp.gpkg'))) {
  
  sf::st_write(vol_stp_utm, dsn = file.path(out_path, 'vol_stp.gpkg'))
  
} else {
  
  print('File vol_stp.gpkg already exists.')
  
}


