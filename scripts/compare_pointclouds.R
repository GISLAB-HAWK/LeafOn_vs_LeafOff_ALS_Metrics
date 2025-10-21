#--------------------------------------------------------------------------
# Name:         compare_pointcloud.R
# Description:  comparing the harmonized laz files in the point cloud level
#               checking pulse count etc. 
# Author:       Svenja Dobelmann
#--------------------------------------------------------------------------

# TBC! 
# compare point clouds

library(lasR)
library(lidR)
library(ggplot2)
library(dplyr)
library(tidyverse)


folder = "R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/03_harmonized/loff24"
outdir = "C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/data/loff24_pc_description.rds"
ctg = readLAScatalog(folder)

read <- reader()
filter <- delete_points("Classification != 2") # keep only ground returns
summ <- summarise(metrics = "count")

pipeline <- read + summ


ans = exec(pipeline, on = ctg, progress = T)

saveRDS(ans, file = outdir)


###########################
loff24 <- readRDS(outdir)
lon23 <- readRDS("C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/data/lon23_pc_description.rds")

ggplot(ans$metrics, aes(x = count)) +
  geom_histogram(fill = "forestgreen", color = "white", bins = 20) +
  labs(
    title = "Distribution of Ground Points per Tile",
    x = "Ground Point Count",
    y = "Number of Tiles"
  ) +
  theme_minimal()

#### points per return
ppr <- bind_rows(
  data.frame(
    return = as.numeric(names(lon23$npoints_per_return)),
    npoints = as.numeric(lon23$npoints_per_return),
    dataset = "lon23"
  ),
  data.frame(
    return = as.numeric(names(loff24$npoints_per_return)),
    npoints = as.numeric(loff24$npoints_per_return),
    dataset = "loff24"
  )
)

p_ppr <- ggplot(df, aes(x = factor(return), y = npoints, fill = dataset)) +
  geom_col(position = "dodge") + 
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Comparison of number of points per return",
    x = "Return number",
    y = "npoints",
    fill = "Season"
  ) +
  theme_minimal(base_size = 14)

#### points per class 
ppc <- bind_rows(
  data.frame(class = names(lon23$npoints_per_class),
             npoints = as.numeric(lon23$npoints_per_class),
             dataset = "lon23"),
  data.frame(class = names(loff24$npoints_per_class),
             npoints = as.numeric(loff24$npoints_per_class),
             dataset = "loff24")) %>%
    mutate(
      class = recode(class,
                     `1` = "Vegetation",
                     `2` = "Ground"))

ggplot(ppc, aes(x = factor(class), y = npoints, fill = dataset)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "Points per class", x = "", y = "Number of points")




#### height distribution
# Combine both z_histograms into one data frame

z_lon23 <-   data.frame(
    z = as.numeric(names(lon23$z_histogram)),
    count = as.numeric(lon23$z_histogram),
    dataset = "lon23")
z_loff24 <-
  data.frame(
    z = as.numeric(names(loff24$z_histogram)),
    count = as.numeric(loff24$z_histogram),
    dataset = "loff24"
  )

z_all <- bind_rows(z_lon23, z_loff24)

# Plot
ggplot(z_all, aes(x = z, y = count, fill = dataset)) +
  geom_col(position = "identity", width = 2, alpha = 0.7) +
  xlim(c(150,500)) + 
  labs(
    title = "Z (Height) Distribution Comparison",
    subtitle = "Comparison of point counts per height bin",
    x = "Elevation (Z value)",
    y = "Number of points",
    color = "Season"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )


#### Canopy ground ratio
# Extract counts
veg_lon23 <- lon23$npoints_per_class["1"]
ground_lon23 <- lon23$npoints_per_class["2"]

veg_loff24 <- loff24$npoints_per_class["1"]
ground_loff24 <- loff24$npoints_per_class["2"]

# Compute ratios
ratio_lon23 <- veg_lon23 / ground_lon23
ratio_loff24 <- veg_loff24 / ground_loff24

# Compute fractions
frac_lon23 <- veg_lon23 / (veg_lon23 + ground_lon23)
frac_loff24 <- veg_loff24 / (veg_loff24 + ground_loff24)

# Combine results
ratios <- data.frame(
  dataset = c("lon23", "loff24"),
  vegetation_points = c(veg_lon23, veg_loff24),
  ground_points = c(ground_lon23, ground_loff24),
  canopy_to_ground_ratio = c(ratio_lon23, ratio_loff24),
  vegetation_fraction = c(frac_lon23, frac_loff24)
)

ratios

