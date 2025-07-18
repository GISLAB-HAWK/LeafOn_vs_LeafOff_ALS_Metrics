# compare metrics

library(terra)
library(dplyr)
library(corrplot)
library(ggplot2)


# path to indices files 
loff_file <- "R:/AG_Magdon/projekte/foreslab/data/als/solling/leaf-on_leaf-off_data/solling24_leafoff_indices/data/solling24_leafoff_indices__0_0___.tiff"
lon_file <- "R:/AG_Magdon/projekte/foreslab/data/als/solling/leaf-on_leaf-off_data/solling23_leafon_indices/data/solling23_leafon_indices__0_0___.tiff"

# read raster files 
loff_r <- rast(loff_file)
lon_r <- rast(lon_file)

# crop to same extent 
lon_r <- crop(lon_r, loff_r)

#Create a mask: TRUE where either is NA or 0
mask <- is.na(lon_r) | is.na(loff_r) | (lon_r == 0) | (loff_r == 0)

# mask the two layer 
loff_r <- mask(loff_r, mask, maskvalue = TRUE)
lon_r <- mask(lon_r, mask, maskvalue = TRUE)

# rename bands for clarity 
names(loff_r) <-  paste0(names(loff_r), "_loff")
names(lon_r) <-  paste0(names(lon_r), "_lon")

# plot some bands 
par(mfrow = c(1,2))
plot(loff_r[[3]])
plot(lon_r[[3]])

# combine the two layer
r <- c(loff_r,lon_r)

# extract data
df <- as.data.frame(r, na.rm = TRUE)
lon_df <- as.data.frame(lon_r, na.rm = T)
loff_df <- as.data.frame(loff_r, na.rm = T)

# correlations
cor_matrix <- cor(loff_df,lon_df, use = "pairwise.complete.obs")
round(cor_matrix,2)
corrplot(cor_matrix,method = "color",        # color in upper
         type = "upper",
         tl.col = "black",        # text color
         tl.cex = 0.8,
         addCoef.col = "black",   # add numbers
         number.cex = 0.7,
         diag = TRUE)

## difference map 
diff_stack <- loff_r - lon_r
plot(diff_stack)


## correlation per layer
sapply(1:nlyr(loff_r), function(i)
  cor(values(loff_r[[i]]), values(lon_r[[i]]), use = "complete.obs")
)


## violon plot
df <- data.frame(
  val = c(values(loff_r[[5]]), values(lon_r[[5]])),
  stack = rep(c("loff", "lon"), each = ncell(loff_r[[1]]))
)
ggplot(df, aes(x = stack, y = val, fill = stack)) + geom_violin()

