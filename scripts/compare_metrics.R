#-------------------------------------------------------------------------------
# Name:         compare_metrics.R
# Description:  compare Metrics for the two rasters between leaf-on and leaf-off 
#               conditions.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


library(terra)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggfortify)

setwd("R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/04_indices/")

# path to indices files 
loff_file <- "solling24_leafoff_indices_harm.tiff"
lon_file <-  "solling23_leafon_indices_harm.tiff"


# read raster files  
loff_r <- rast(loff_file) 
lon_r <- rast(lon_file)

names(loff_r[[1]]) <- "species"
names(lon_r[[18]]) <- "species"

plot(lon_r$band1)
plot(loff_r)

# crop to same extent 
lon_r <- crop(lon_r, loff_r)

#Create a mask: TRUE where either is NA or 0
mask <- is.na(lon_r) | is.na(loff_r) | (lon_r == 0) | (loff_r == 0) 
mask_single <- app(mask, fun = function(x) any(x, na.rm = TRUE))

# mask the two layer 
loff_r <- mask(loff_r, mask_single, maskvalue = TRUE)
lon_r <- mask(lon_r, mask_single, maskvalue = TRUE)

# plot some bands 
par(mfrow = c(1,2))
plot(loff_r)
plot(lon_r[[2]])

diff <- loff_r$BE_H_MAX - lon_r$BE_H_MAX

mask <- diff > -10

# mask the two layer where height difference is more than 10 meter 
loff_r <- mask(loff_r, mask, maskvalue = FALSE)
lon_r <- mask(lon_r, mask, maskvalue = FALSE)

plot(lon_r[[2]])

############################

# Get all valid cell indices
valid_idx <- which(values(mask))

#  Randomly sample 5000 of those indices
set.seed(123)  # for reproducibility
sample_idx <- sample(valid_idx, size = 5000)

# Convert those cell indices to spatial coordinates
sample_pts <- xyFromCell(lon_r, sample_idx)

#  Extract corresponding values from both rasters
lon_vals  <- terra::extract(lon_r,  sample_pts)
loff_vals <- terra::extract(loff_r, sample_pts)

lon_vals  <- cbind(sample_pts, lon_vals)
loff_vals <- cbind(sample_pts, loff_vals)
############################

# extract data
lon_df <- as.data.frame(lon_vals,xy = TRUE, na.rm = T)
loff_df <- as.data.frame(loff_vals, xy = TRUE, na.rm = T)

df_long <- bind_rows(
  lon_df  %>% mutate(season = "leaf on"),
  loff_df %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(season, species, x, y), names_to = "variable", values_to = "value") %>%
  mutate(species = case_when(
    species == "1" ~ "decidious",
    species == "2" ~ "coniferous",
    TRUE ~ NA_character_
  ))

# summarize the data
df_summary <- df_long %>%
  group_by(season, variable) %>%
  summarise(mean = mean(value), sd = sd(value))

df_summary

# Pivot wide
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value)%>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`)

#### paired t.test ####

# shapiro test for normal distribution 
set.seed(42)  # for reproducibility

shapiro <- df_long %>%
  group_by(season, variable) %>%
  summarise(
    shapiro = list(shapiro.test(value)), # using a subsample n = 5000
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = shapiro$statistic,
    p.value   = shapiro$p.value,
    method    = shapiro$method
  ) %>%
  select(-shapiro)

shapiro

any(shapiro$p.value>0.05)
which(shapiro$p.value>0.05)


# data not normally distribute. using non-parametric wilcoxon-test instead

wilcox_test <- df_long %>%
  group_by(variable) %>%
  summarise(
    test = list(
      wilcox.test(value ~ season, data = cur_data())
    ),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = test$statistic,
    p.value   = test$p.value,
    method    = test$method
  ) %>%
  select(-test)

print(wilcox_test)

any(wilcox_test$p.value< 0.05)
which(wilcox_test$p.value< 0.05)

# Pivot data wider: one row per variable

df_wide <- df_long %>%
  group_by(variable, season) %>%
  mutate(pixel_id = row_number()) %>%   # per-variable pixel index
  ungroup() %>%
  pivot_wider(
    names_from  = season,
    values_from = value
  ) %>%
  arrange(variable, pixel_id) %>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`)

# change order of variables for later plotting
df_wide$variable <- factor(df_wide$variable, 
                           levels = c("BE_H_MAX","BE_H_P90","BE_H_P80","BE_H_P50","BE_H_P20",
                                      "BE_H_P10", "BE_PR_02","BE_PR_10", "BE_RD_02", "BE_RD_10", 
                                      "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
                                      "point_density", "pulse_returns_mean"))  
df_long$variable <- factor(df_long$variable, 
                           levels = c("BE_H_MAX","BE_H_P90","BE_H_P80","BE_H_P50","BE_H_P20",
                                      "BE_H_P10", "BE_PR_02","BE_PR_10", "BE_RD_02", "BE_RD_10", 
                                      "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
                                      "point_density", "pulse_returns_mean"))  

## scatterplot 
s <- ggplot(df_wide, aes(x = leaf_on, y = leaf_off, colour = species)) +
  geom_point(alpha = 0.8) +
  #geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~ variable, scales = "free") +
  ggpubr::stat_cor(
    aes(color = species, label = ..r.label..),        # <- compute a correlation per species
    method = "pearson",
    label.x.npc = "left",
    cor.coef.name = "r", 
    size = 6, 
    show.legend = FALSE,
    na.rm = T
  ) +
  theme_minimal() +
  scale_color_discrete(na.value = "lightgray") +
  labs(
    title = "Comparison of Metrices on the pixel level",
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  ) +
  theme(
    legend.text = element_text(size = 16),  
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 20),      
    strip.text = element_text(size = 17), 
    title = element_text(size = 20)
  )

s

ggsave("C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/plots/pixel_scatterplot.pdf",s,  dpi = 500, width = 15, height = 10)

df_summary <- df_long %>%
  group_by(season, variable) %>%
  summarise(mean = mean(value), sd = sd(value))

## violon plot

ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 5, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "none")

## density plot
ggplot(df_sub, aes(x = value, fill = season)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free", ncol = 4) +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "none")


# correlations
cor_matrix <- cor(loff_df,lon_df, use = "pairwise.complete.obs")
round(cor_matrix,2)
corrplot(cor_matrix,method = "color",        # color in upper
         type = "upper",
         tl.col = "black",        # text color
         tl.cex = 0.7,
         addCoef.col = "black",   # add numbers
         number.cex = 0.7,
         diag = TRUE)




## correlation per layer
sapply(1:nlyr(loff_r), function(i)
  cor(values(loff_r[[i]]), values(lon_r[[i]]), use = "complete.obs")
)




## PCA
combined <- bind_rows(loff_df %>% mutate(source = "df1"),
                      lon_df %>% mutate(source = "df2"))
pca <- prcomp(combined %>% select(-source), scale. = TRUE)
summary(pca)
autoplot(pca, data = combined, colour = "source", loadings = TRUE, loadings.label = TRUE)
