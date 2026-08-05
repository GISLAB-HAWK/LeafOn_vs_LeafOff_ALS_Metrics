#-------------------------------------------------------------------------------
# Name:         compare_metrics.R
# Description:  Compare ABA metrics between leaf-on and leaf-off conditions
#               at pixel level, using a balanced sample of validated pixels.
#               Produces descriptive summaries, paired Wilcoxon tests and
#               scatter / difference-density figures.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R')


# 01 - configuration
#-------------------

# Point density of the input rasters.
# Used for both input paths and output figure names.
PPM  <- "ppm20"

# Input rasters per season: subdirectory and file basename
inputs <- list(
  leaf_on = list(
    dir  = "metrics_leafon_2023",
    file = paste0("solling23_lon_", PPM, "_metrics.tiff")
  ),
  leaf_off = list(
    dir  = "metrics_leafoff_2024",
    file = paste0("solling24_loff_", PPM, "_metrics.tiff")
  )
)

# Sample of validated pixels
SAMPLE_FILE <- "sample_selection_n400_balanced.csv"


# Variable order used for tables and facets (24 metrics)
lvl <- c("BE_H_MAX", "BE_H_P95", "BE_H_P90", "BE_H_P80", "BE_H_P50", "BE_H_P20",
         "BE_H_P10", "BE_PR_02", "BE_PR_05", "BE_PR_10", "BE_PR_20", "BE_PR_30",
         "BE_RD_02", "BE_RD_05", "BE_RD_10", "BE_RD_20", "BE_RD_30",
         "BE_H_KURTOSIS", "BE_H_SKEW", "BE_H_VAR", "BE_H_SD", "BE_H_MEAN",
         "point_density", "pulse_returns_mean")


# 02 - data preparation
#----------------------

lon_r <- rast(
  file.path(processed_data_dir, inputs$leaf_on$dir, PPM, inputs$leaf_on$file)
)

loff_r <- rast(
  file.path(processed_data_dir, inputs$leaf_off$dir, PPM, inputs$leaf_off$file)
)

# Validated sample pixels
pix_valid <- read.csv(file.path(metadata_dir, SAMPLE_FILE)) %>%
  filter(status == "valid")


pts <- vect(pix_valid, geom = c("x", "y"), crs = as.character(crs(lon_r)))

# Extract corresponding values from both rasters
lon_vals  <- terra::extract(lon_r,  pts, xy = TRUE, cells = TRUE)
loff_vals <- terra::extract(loff_r, pts, xy = TRUE, cells = TRUE)

lon_df  <- as.data.frame(lon_vals,  xy = TRUE, na.rm = TRUE) %>% rename('ts' = class_name)
loff_df <- as.data.frame(loff_vals, xy = TRUE, na.rm = TRUE) %>% rename('ts' = class_name)

# Combine both seasons and reshape to long format
df_long <- bind_rows(
  lon_df  %>% mutate(season = "leaf on"),
  loff_df %>% mutate(season = "leaf off")
) %>%
  pivot_longer(cols = -c(season, ts, ID, x, y, cell),
               names_to = "variable", values_to = "value") %>%
  mutate(species = case_when(
    ts == "1" ~ "deciduous",
    ts == "2" ~ "coniferous",
    TRUE ~ NA_character_
  )) %>%
  dplyr::select(-ts) %>%
  mutate(variable = factor(variable, levels = lvl))

# Wide format with paired values and their difference
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value) %>%
  rename(leaf_on  = `leaf on`,
         leaf_off = `leaf off`) %>%
  mutate(diff = leaf_on - leaf_off)


# 03 - descriptive summaries
#---------------------------

# Per season and variable
df_summary_long <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value), sd = sd(value),
                   min = min(value), max = max(value),
                   .groups = "drop")

print(df_summary_long, n = Inf)

# Paired, per variable
df_summary_wide <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(mean_lon  = mean(leaf_on),  sd_lon  = sd(leaf_on),
                   mean_loff = mean(leaf_off), sd_loff = sd(leaf_off),
                   .groups = "drop") %>%
  mutate(diff = mean_lon - mean_loff)

print(df_summary_wide, n = Inf)


# 04 - statistical tests
#-----------------------

# Test for normality of the paired differences.
shapiro <- df_wide %>%
  #filter(species == "deciduous") %>%
  #filter(species == "coniferous") %>%
  group_by(variable) %>%
  dplyr::summarise(
    shapiro = list(shapiro.test(diff)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    statistic = shapiro$statistic,
    p.value   = shapiro$p.value,
    method    = shapiro$method
  ) %>%
  dplyr::select(-shapiro, -method) %>%
  mutate(
    signif = case_when(
      p.value < 0.01 ~ "***",
      p.value < 0.05 ~ "**",
      p.value < 0.1  ~ "*",
      TRUE ~ ""
    ))

print(shapiro, n = 34)
any(shapiro$p.value > 0.05)

# Differences are not normally distributed -> paired Wilcoxon test
# instead of a paired t-test.

#' Run paired Wilcoxon tests for all variables in one subset
#'
#' @param data Data frame in wide format with columns variable, leaf_on, leaf_off.
#'
#' @return Tibble with one row per variable: sample size, mean and median
#'   difference, test statistic, raw and BH-adjusted p-value, rank-biserial
#'   effect size and ratios.
run_wilcox <- function(data) {
  data %>%
    tidyr::drop_na(leaf_on, leaf_off) %>%
    group_by(variable) %>%
    dplyr::summarise(
      n_pairs     = n(),
      mean_diff   = mean(leaf_on - leaf_off),
      median_diff = median(leaf_on - leaf_off),
      med_on      = median(leaf_on),
      med_off     = median(leaf_off),
      test        = list(wilcox.test(leaf_on, leaf_off, paired = TRUE, exact = FALSE)),
      W           = test[[1]]$statistic,
      p.value     = test[[1]]$p.value,
      effsize     = list(effectsize::rank_biserial(
        x = leaf_on, y = leaf_off, paired = TRUE, ci = 0.95, verbose = FALSE
      )),
      r_rb        = effsize[[1]]$r_rank_biserial,
      .groups = "drop"
    ) %>%
    mutate(
      # Uniform ratio: > 1 = leaf-on higher, < 1 = leaf-off higher.
      # NA for signed distribution measures and denominators close to zero.
      ratio = if_else(
        str_detect(variable, "KURTOSIS|SKEW") | abs(med_off) < 1e-9,
        NA_real_,
        med_on / med_off
      ),
      # PR only: interception instead of transmission (for the text)
      int_ratio = if_else(
        str_detect(variable, "_PR_") & med_off < 1 - 1e-9,
        (1 - med_on) / (1 - med_off),
        NA_real_
      ),
      p_adj  = p.adjust(p.value, method = "BH"),
      signif = case_when(
        p_adj < 0.01 ~ "***",
        p_adj < 0.05 ~ "**",
        p_adj < 0.1  ~ "*",
        TRUE         ~ ""
      )
    ) %>%
    dplyr::select(-test, -effsize)
}

results <- list(
  all        = run_wilcox(df_wide),
  coniferous = run_wilcox(df_wide %>% filter(species == "coniferous")),
  deciduous  = run_wilcox(df_wide %>% filter(species == "deciduous"))  
)

invisible(lapply(results, function(df) print(df, n = Inf)))


#--- export test results ---------------------------------------------------
  

# Normality tests
write.csv(
  shapiro,
  file.path(output_dir, "stats", paste0("shapiro_", PPM, ".csv")),
  row.names = FALSE
)

# Wilcoxon results: one combined table with a subset column
results_tbl <- bind_rows(results, .id = "subset") 

write.csv(
  results_tbl,
  file.path(output_dir, "stats", paste0("wilcoxon_", PPM, ".csv")),
  row.names = FALSE
)

# Same content as RDS, keeps factor levels and full precision for reuse in R
saveRDS(
  results,
  file.path(output_dir, "stats", paste0("wilcoxon_", PPM, ".rds"))
)

saveRDS(
  df_wide,
  file.path(output_dir, "stats", paste0("df_wide_", PPM, ".rds"))
)

# Descriptive summaries
write.csv(
  df_summary_wide,
  file.path(output_dir, "stats", paste0("summary_paired_", PPM, ".csv")),
  row.names = FALSE
)


message("test results written to: ",  file.path(output_dir, "stats"))
