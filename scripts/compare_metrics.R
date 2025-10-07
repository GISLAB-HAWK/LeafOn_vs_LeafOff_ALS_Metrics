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

setwd("R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/03_indices/")

# path to indices files 
loff_file <- "solling24_leafoff_indices_harm.tiff"
lon_file <-  "solling23_leafon_indices_harm.tiff"


# read raster files 
loff_r <- rast(loff_file)[[-1]]
lon_r <- rast(lon_file)[[-1]]

plot(lon_r)
plot(loff_r)

# crop to same extent 
lon_r <- crop(lon_r, loff_r)

#Create a mask: TRUE where either is NA or 0
mask <- is.na(lon_r) | is.na(loff_r) | (lon_r == 0) | (loff_r == 0)

# mask the two layer 
loff_r <- mask(loff_r, mask, maskvalue = TRUE)
lon_r <- mask(lon_r, mask, maskvalue = TRUE)

# plot some bands 
par(mfrow = c(1,2))
plot(loff_r[[3]])
plot(lon_r[[3]])

diff <- loff_r$BE_H_MAX - lon_r$BE_H_MAX

mask <- diff > -10

# mask the two layer 
loff_r <- mask(loff_r, mask, maskvalue = FALSE)
lon_r <- mask(lon_r, mask, maskvalue = FALSE)

plot(lon_r[[2]])

# extract data
lon_df <- as.data.frame(lon_r, na.rm = T)
loff_df <- as.data.frame(loff_r, na.rm = T)

df_long <- bind_rows(
  lon_df  %>% mutate(season = "leaf on"),
  loff_df %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(season), names_to = "variable", values_to = "value")

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
# using a subsample of 5000
set.seed(42)
df_sub <- df_long %>% dplyr::sample_n(5000) 

# test for normal distribution 
set.seed(42)  # for reproducibility

shapiro <- df_sub %>%
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

wilcox_test <- df_sub %>%
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

df_wide <- df_sub %>%
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

## scatterplot 
ggplot(df_wide, aes(x = leaf_on, y = leaf_off)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal() +
  labs(
    title = "Leaf-On vs Leaf-Off Comparison per Variable",
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  )


df_summary <- df_sub %>%
  group_by(season, variable) %>%
  summarise(mean = mean(value), sd = sd(value))

## violon plot

ggplot(df_sub, aes(x = season, y = value, fill = season, alpha = 0.85)) +
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
