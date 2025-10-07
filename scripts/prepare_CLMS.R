#--------------------------------------------------------------------------
# Name:         prepare_CLMS.R
# Description:  merge, resample, and clip the CLMS dominant tree species
#               data to extent of Solling ALS dataset
# Author:       Svenja Dobelmann
#--------------------------------------------------------------------------


library(terra)

# Paths
t1 <- "C:/Users/sdobelma/Downloads/CLMS_HRLVLCC_DLT_S2020_R10m_E42N31_03035_V01_R00/CLMS_HRLVLCC_DLT_S2020_R10m_E42N31_03035_V01_R00.tif"
t2 <- "C:/Users/sdobelma/Downloads/CLMS_HRLVLCC_DLT_S2020_R10m_E43N31_03035_V01_R00/CLMS_HRLVLCC_DLT_S2020_R10m_E43N31_03035_V01_R00.tif"
t3 <- "R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/03_indices/Solling23_leafon_indices_harm.tiff"   # the file that provides the extent

out_clipped <- "R:/AG_Magdon/datensaetze/solling/dobelmann/tree_species_data/CLMS_DLT2020_clipped.tif"

# Read rasters
r1 <- rast(t1)
r2 <- rast(t2)
r3 <- rast(t3)[[1]]

# Merge the two tiles
mrg <- mosaic(r1, r2)

# Align CRSs: reproject everything to ESPG: 25832
if (!all(crs(r1) == crs(r3), crs(r2) == crs(r3))) {
  mrg <- project(mrg, "epsg:25832", method = "near")
  r3 <- project(r3, mrg) 
}

# mask to r3
mask <- is.na(r3)
masked <- mask(mrg, mask, maskvalues = TRUE)

# trim to remove NAs
trimmed <- trim(masked, padding=0, value=NA)

plot(trimmed)

# 6) Write outputs
writeRaster(trimmed, out_clipped, overwrite = TRUE)


