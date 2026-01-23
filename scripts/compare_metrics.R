#-------------------------------------------------------------------------------
# Name:         compare_metrics.R
# Description:  compare Metrics for the two rasters between leaf-on and leaf-off 
#               conditions.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

# source setup script
source('src/setup.R', local = TRUE)

#### Data preparation ####

# read raster files  
loff_r <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling24_loff_ppm20_indices.tiff')
               ) 
lon_r <- rast(
  file.path(processed_data_dir, 'metrics','pix_level','solling23_lon_ppm20_indices.tiff')
  )
 

# crop to same extent 
lon_r <- crop(lon_r, loff_r)

# mask only valid cells
# 1. masking out NA
is_na <- is.na(lon_r) | is.na(loff_r)
mask_na <- app(is_na, fun = function(x) any(x, na.rm = TRUE))   

# 2. masking out height difference >10m (assuming clearcut)
diff <- loff_r$BE_H_MAX - lon_r$BE_H_MAX
mask_diff <- diff <= -10

# 3. masking out non forested areas 
mask_nontree <- lon_r$DLT == 0

# combine the three masks 
mask_combined <- mask_na | mask_diff | mask_nontree

# apply the mask 
loff_r <- mask(loff_r, mask_combined, maskvalue = TRUE)
lon_r <- mask(lon_r, mask_combined, maskvalue = TRUE)

# plot some bands 
par(mfrow = c(1,2))
plot(loff_r[[2]])
plot(lon_r[[2]])


# Get all valid cell indices
valid_idx <- which(!is.na(values(lon_r[[2]])))

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

# extract data
lon_df <- as.data.frame(lon_vals,xy = TRUE, na.rm = T) 
loff_df <- as.data.frame(loff_vals, xy = TRUE, na.rm = T) 

# combine datasets and bring into long format 
df_long <- bind_rows(
  lon_df  %>% mutate(season = "leaf on"),
  loff_df %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(season, DLT, x, y), names_to = "variable", values_to = "value") %>%
  mutate(species = case_when(
    DLT == "1" ~ "decidious",
    DLT == "2" ~ "coniferous",
    TRUE ~ NA_character_
  )) %>%
  select(-DLT)

# summarize the data
df_summary <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), sd = sd(value), min = min(value), max = max(value))

print(df_summary)

# brint data to wide format 
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value)%>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`)

# summarize the data
df_summary <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(mean_lon = mean(leaf_on), sd_lon = sd(leaf_on), mean_loff = mean(leaf_off), sd_loff = sd(leaf_off)) %>%
  mutate(diff = mean_lon-mean_loff)

print(df_summary)

#### wilcoxon-test ####

# first, testing for normal distribution 
shapiro <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(
    shapiro = list(shapiro.test(value)), 
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = shapiro$statistic,
    p.value   = shapiro$p.value,
    method    = shapiro$method
  ) %>%
  select(-shapiro, -method) %>%
  dplyr::mutate(signif = ifelse(p.value < 0.05, "***", ""))

print(shapiro, n = 34)
any(shapiro$p.value>0.05)

# data not normally distributed. using non-parametric wilcoxon-test instead of simple t-test
wilcox_test <- df_wide %>%
  tidyr::drop_na(leaf_on, leaf_off) %>%   # remove NAs
  group_by(variable) %>%
  dplyr::summarise(
    n_pairs  = n(),
    mean_diff = mean(leaf_on - leaf_off),  # mean difference (on - off)
    p.value  = wilcox.test(leaf_on, leaf_off, paired = TRUE)$p.value,
    .groups = "drop"
  )  %>%
  dplyr::mutate(signif = ifelse(p.value < 0.05, "***", ""))

print(wilcox_test)

any(wilcox_test$p.value > 0.05)

# significant differences for all metrics, wilcoxon test very sensitive especially with high sample sizes 

# change order of variables for later plotting
lvl = c("BE_H_MAX","BE_H_P90","BE_H_P80","BE_H_P50","BE_H_P20",
        "BE_H_P10", "BE_PR_02","BE_PR_10", "BE_RD_02", "BE_RD_10", 
        "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
        "point_density", "pulse_returns_mean")

df_wide$variable <- factor(df_wide$variable, 
                           levels = lvl )  
df_long$variable <- factor(df_long$variable, 
                           levels = lvl)  


#### scatterplot ####

panel_labels <- df_wide %>%
  distinct(variable) %>%
  arrange(variable) %>%              # ensures stable order
  mutate(
    label = paste0(letters[seq_along(variable)], ")")
  )

cols <- c(
  decidious  = "#1b9e77",
  coniferous = "#d95f02"
)

s <- ggplot(df_wide, aes(x = leaf_on, y = leaf_off, colour = species, linetype = species)) +
  geom_point(alpha = 0.3, size = 0.7,aes(shape = species)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "grey20") +
  facet_wrap(~ variable, labeller = as_labeller(letter_labeller), scales = "free", ncol = 3, nrow = 8) +
  geom_text(
    data = panel_labels,
    aes(label = label),
    x = -Inf, y = Inf,                # top-left corner
    hjust = -0.4, vjust = 1.4,
    size = 3,
    inherit.aes = FALSE
  ) + 
  ggpubr::stat_cor(
    aes(label = ..r.label..),
    method = "spearman",
    cor.coef.name = "rho",
    size = 2.5,
    show.legend = FALSE,
    na.rm = TRUE,
    geom = "text",          
    label.x.npc = 0.7,    
    label.y.npc = 0.2, 
    lineheight = 0.7 
  ) + 
  theme_classic(base_size = 12) +
  labs(
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  ) +
  theme(
    legend.text = element_text(size = 10),  
    axis.title = element_text(size = 10, color = "darkgrey"),
    axis.text = element_text(size = 8, color = "darkgrey"),
    axis.ticks = element_line(color = "darkgrey"),
    axis.line = element_line(color = "darkgrey"),
    legend.title = element_blank(), 
    legend.position = c(0.97, 0),
    legend.justification = c("right", "bottom"),
    strip.text = element_blank(), 
    panel.spacing = unit(0.6, "lines")
  )


print(s)
ggsave(paste0(output_dir,"/pix_scatterplot_ppm20.png"),s, units = "cm", dpi = 350, width = 14, height = 20)


#### violon plot ####

v <- ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 4, ncol = 5, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric comparison: area based ",
       subtitle = "n = 196, pulse density = 20ppm",
       x = "",
       y = "Value") +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c("leaf on" = "#1b9e77", 
               "leaf off" = "#d95f02"),
    na.value = "lightgray"
  )

print(v)
#ggsave(paste0(output_dir,"/pix_violinplot_ppm20.pdf"), v ,  dpi = 500, width = 15, height = 10)

#### density plot ####
d <- ggplot(df_long, aes(x = value, fill = season)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  theme_minimal() +
  labs(title = "Metric comparison for BI plots ",
       subtitle = "n = 196, pulse density = 20ppm",
       x = "",
       y = "Value") +
  theme(legend.position = "bottom") +
  scale_fill_manual(
    values = c("leaf on" = "#1b9e77", 
               "leaf off" = "#d95f02"),
    na.value = "lightgray"
  )

print(d)
#ggsave(paste0(output_dir,"/pix_densityplot_ppm20.pdf"), d ,  dpi = 500, width = 12, height = 10)


#### Principle Conponent Analysis ####
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



#### multiple linear regression ####

df_long <- df_wide %>%
  pivot_longer(
    c(leaf_on, leaf_off),
    names_to = "leaf_condition",
    values_to = "value"
  ) %>%
  mutate(leaf_condition = factor(leaf_condition,
                                 levels = c("leaf_off", "leaf_on")))

model <- lm(
  value ~ leaf_condition + species,
  data = df_long %>% filter(variable == "BE_H_MEAN")
)

summary(model)

#### Friedman test ####
friedman.test(cbind(read, write, math))

#### repeated measures logistic regression

