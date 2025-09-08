# compare metrics

library(terra)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggfortify)


# path to indices files 
loff_file <- "R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/03_indices/loff24__indices/data/solling24_leafoff_indices__0_0___.tiff"
lon_file <-  "R:/AG_Magdon/datensaetze/solling/dobelmann/leaf-on_leaf-off_data/03_indices/lon23_indices/data/solling23_leafon_indices__0_0___.tiff"


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

# plot some bands 
par(mfrow = c(1,2))
plot(loff_r[[3]])
plot(lon_r[[3]])


# extract data
lon_df <- as.data.frame(lon_r, na.rm = T)
loff_df <- as.data.frame(loff_r, na.rm = T)

## unpaired t-test
map_dfr(names(lon_df), function(var) {
  data.frame(
    variable = var,
    p_value = t.test(lon_df[[var]], loff_df[[var]])$p.value
  )
})


lon_df$season <- "leaf on"
loff_df$season<- "leaf off"


# Combine and pivot to long format
df_long <- bind_rows(lon_df, loff_df) %>%    
  pivot_longer(
    cols = -season,
    names_to = "variable",
    values_to = "value"
  )


# Pivot data wider: one row per variable
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value)%>%
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


df_summary <- df_long %>%
  group_by(season, variable) %>%
  summarise(mean = mean(value), sd = sd(value))

## violon plot

ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 3, ncol = 4, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "none")

## density plot
ggplot(df_long, aes(x = value, fill = season)) +
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
