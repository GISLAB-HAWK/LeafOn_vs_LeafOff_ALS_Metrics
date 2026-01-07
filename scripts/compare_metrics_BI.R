#-------------------------------------------------------------------------------
# Name:         compare_metrics_BI.R
# Description:  compare Metrics for the BI plots between leaf-on and leaf-off 
#               conditions.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

# source setup script
source('src/setup.R', local = TRUE)


#### Data preparation ####
# read the leaf-off data 
loff_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling24_loff_ppm20_indices_BI.csv')
  )

# read the leaf-on data 
lon_df <-read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling23_lon_ppm20_indices_BI.csv')
)

# get tree species information from plot file 
ts <-read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','vol_stp_GR.csv')
)


ts <- ts %>%
  mutate(dominant_species = case_when(
    dominant_species == "LB" ~ "decidious",
    dominant_species == "NB" ~ "coniferous",
    TRUE ~ NA_character_
  ))

# filter complete cases 
keep <- intersect(lon_df$name, loff_df$name)

# combine datasets in brint into long format 
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

# change order of variables for later plotting

lvl = c("BE_H_MAX","BE_H_P90","BE_H_P80","BE_H_P50","BE_H_P20",
        "BE_H_P10", "BE_PR_02","BE_PR_10", "BE_RD_02", "BE_RD_10", 
        "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
        "point_density", "pulse_returns_mean")

df_wide$variable <- factor(df_wide$variable, 
                           levels = lvl )  
df_long$variable <- factor(df_long$variable, 
                           levels = lvl)  


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

which(shapiro$p.value > 0.05) # lon skewness and SD normally distributed, rest non-normally

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

which(wilcox_test$p.value > 0.05)

# significant differences for all pairs, exept BE_RD_10

#### scatterplot ####

# keeping only samples where species is known 
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
    title = "Metric comparison for BI plots ",
    subtitle = "n = 196, pulse density = 20ppm",
    x = "Leaf-On Value",
    y = "Leaf-Off Value"
  ) +
  theme(
    legend.text = element_text(size = 14),  
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 18),      
    strip.text = element_text(size = 15), 
    title = element_text(size = 18)
  )


print(s)
#ggsave(paste0(output_dir,"/BI_scatterplot_ppm20.pdf"),s,  dpi = 500, width = 15, height = 10)


#### violin plot ####
v <- ggplot(df_long, aes(x = season, y = value, fill = season, alpha = 0.85)) +
  geom_violin(trim = FALSE) +
  facet_wrap(~ variable, nrow = 4, ncol = 5, scales = "free_y") +
  theme_minimal() +
  labs(title = "Metric comparison for BI plots ",
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
#ggsave(paste0(output_dir,"/BI_violinplot_ppm20.pdf"), v ,  dpi = 500, width = 12, height = 10)


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
#ggsave(paste0(output_dir,"/BI_densityplot_ppm20.pdf"), d ,  dpi = 500, width = 12, height = 10)




#################################### END OF SCRIPT #############################
