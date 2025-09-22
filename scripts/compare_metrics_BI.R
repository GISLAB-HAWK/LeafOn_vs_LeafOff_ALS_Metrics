#-------------------------------------------------------------------------------
# Name:         compare_emtrics_BI.R
# Description:  compare Metrics for the BI plots between leaf-on and leaf-off 
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
library(ggpubr)

setwd("C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/")

# files paths 
loff_df_file <- "data/BI_indices_2024_RTK_harm.csv"
lon_df_file <-  "data/BI_indices_2023_RTK_harm.csv"
tree_species_file <- "data/vol_stp_GR.csv"

# read the data 
loff_df <- read.csv(loff_df_file)#[,-1]
lon_df <- read.csv(lon_df_file)#[,-1]
ts <- read.csv(tree_species_file)#[,-1]

ts <- ts %>%
  mutate(dominant_species = case_when(
    dominant_species == "LB" ~ "decidious",
    dominant_species == "NB" ~ "coniferous",
    TRUE ~ NA_character_
  ))


keep <- intersect(lon_df$name, loff_df$name)

df_long <- bind_rows(
  lon_df  %>% filter(name %in% keep) %>% mutate(season = "leaf on"),
  loff_df %>% filter(name %in% keep) %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(name, season), names_to = "variable", values_to = "value") %>%
  left_join(ts %>% select(kspnr, dominant_species),
                                    by = c("name" = "kspnr"))  %>%
  rename(Tree_ID = name, species = dominant_species)


# summarize the data
df_summary <- df_long %>%
  group_by(season, variable) %>%
  summarise(mean = mean(value), sd = sd(value))

print(df_summary)

# Pivot wide
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value)%>%
  rename(leaf_on = `leaf on`,
         leaf_off = `leaf off`)

# summarize the data
df_summary <- df_wide %>%
  group_by(variable) %>%
  summarise(mean = mean(value), sd = sd(value))

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
#### paired t.test ####

# test for normal distribution 
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

print(shapiro, n = 34)

which(shapiro$p.value > 0.05) # lon skewness and SD normally distributed, rest non-normally

# data not normally distribute. using non-parametric wilcoxon-test instead
wilcox_test <- df_wide %>%
  tidyr::drop_na(leaf_on, leaf_off) %>%   # keep only complete pairs
  group_by(variable) %>%
  summarise(
    n_pairs  = n(),
    mean_diff = mean(leaf_on - leaf_off),  # mean difference (on - off)
    p_value  = wilcox.test(leaf_on, leaf_off, paired = TRUE)$p.value,
    .groups = "drop"
  )

print(wilcox_test)

which(wilcox_test$p_value > 0.05)

# significant differences for all pairs, exept BE_RD_10

#### scatterplot ####

df_wide_clean <- df_wide %>% 
  filter(!is.na(species))

s <- ggplot(df_wide_clean, aes(x = leaf_on, y = leaf_off, colour = species)) +
  geom_point() +
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
    title = "Comparison of Metrices on the plot level",
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


print(s)
ggsave("plots/BI_scatterplot.pdf",s,  dpi = 500, width = 15, height = 10)


#### violin plot ####
v <- ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 4, ncol = 5, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c("leaf on" = "#1b9e77", 
               "leaf off" = "#d95f02"),
    na.value = "lightgray"
  )

print(v)
ggsave("plots/BI_violinplot.png", v ,  dpi = 500, width = 12, height = 10)


#### density plot ####
d <- ggplot(df_long, aes(x = value, fill = season)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "bottom") +
  scale_fill_manual(
    values = c("leaf on" = "#1b9e77", 
               "leaf off" = "#d95f02"),
    na.value = "lightgray"
  )

print(d)

ggsave("plots/BI_densityplot.png", d ,  dpi = 500, width = 12, height = 10)

