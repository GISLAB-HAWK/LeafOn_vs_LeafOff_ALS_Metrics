# compare metrics

library(terra)
library(dplyr)
library(tidyverse)
library(corrplot)
library(ggplot2)
library(ggfortify)
library(ggpubr)

setwd("C:/Users/sdobelma/Documents/LeafOn_vs_LeafOff_ALS/")

# files paths 
loff_df_file <- "data/BI_indices_2024_RTK.csv"
lon_df_file <-  "data/BI_indices_2023_RTK.csv"


# read the data 
loff_df <- read.csv(loff_df_file)#[,-1]
lon_df <- read.csv(lon_df_file)#[,-1]

keep <- intersect(lon_df$name, loff_df$name)

df_long <- bind_rows(
  lon_df  %>% filter(name %in% keep) %>% mutate(season = "leaf on"),
  loff_df %>% filter(name %in% keep) %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(name, season), names_to = "variable", values_to = "value") %>%
  rename(Tree_ID = name)


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

# paired t.test
t_test <- df_wide %>%
  tidyr::drop_na(leaf_on, leaf_off) %>%   # keep only complete pairs
  group_by(variable) %>%
  summarise(
    n_pairs  = n(),
    mean_diff = mean(leaf_on - leaf_off),  # mean difference (on - off)
    p_value  = t.test(leaf_on, leaf_off, paired = TRUE)$p.value,
    .groups = "drop"
  )
print(t_test)

## scatterplot 
s <- ggplot(df_wide, aes(x = leaf_on, y = leaf_off)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free") +
  ggpubr::stat_cor(method = "pearson", 
           label.x.npc = "left", 
           label.y.npc = "top") +   # position of label
  theme_minimal() +
  labs(
    title = "Leaf-On vs Leaf-Off Comparison per Variable",
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  )


ggsave("plots/BI_scatterplot.png",s,  dpi = 500, width = 12, height = 10)


## violon plot
v <- ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 4, ncol = 5, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "none")

ggsave("plots/BI_violinplot.png", v ,  dpi = 500, width = 12, height = 10)


## density plot
d <- ggplot(df_long, aes(x = value, fill = season)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  theme_minimal() +
  labs(title = "Metric Comparison Leaf off vs. leaf on",
       x = "",
       y = "Value") +
  theme(legend.position = "bottom")

d

ggsave("plots/BI_densityplot.png", d ,  dpi = 500, width = 12, height = 10)

