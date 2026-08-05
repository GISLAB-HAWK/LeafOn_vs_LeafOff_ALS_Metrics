#-------------------------------------------------------------------------------
# Name:         compare_structural_metrics.R
# Description:  Compare the structural metrics from calculate_structural_metrics.R
#               between leaf-on and leaf-off conditions at pixel level, using a
#               balanced sample of validated pixels. Produces descriptive
#               summaries and paired Wilcoxon tests.
#               Figures: compare_structural_metrics_viz.R
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

source('src/setup.R')
source('src/stats_func.R') # load the functions for the paired tests


# 01 - configuration
#-------------------

# Point density of the input rasters.
# Used for both input paths and output file names.
PPM <- "ppm20"

# Input rasters per season: subdirectory, file basename and season label
inputs <- list(
  leaf_on = list(
    dir   = "metrics_leafon_2023",
    file  = paste0("solling23_lon_", PPM, "_struct_metrics.tiff"),
    label = "leaf on"
  ),
  leaf_off = list(
    dir   = "metrics_leafoff_2024",
    file  = paste0("solling24_loff_", PPM, "_struct_metrics.tiff"),
    label = "leaf off"
  )
)

# Sample of validated pixels
SAMPLE_FILE <- "sample_selection_n400_balanced.csv"

# Variable order used for tables and facets
lvl <- c("ENL0D", "ENL1D", "ENL2D", "Box_dimension", "VCI",
         "CHM_mean", "Canopy_cover", "Rumple")

# Tree species classes of the template raster
species_lab <- c("1" = "deciduous", "2" = "coniferous")

stats_dir <- file.path(output_dir, "stats")


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

# Combine both seasons and reshape to long format.
# x, y and cell are dropped: they differ between the two rasters in the last
# decimals and would otherwise act as pivot keys, splitting the pairs.
df_long <- bind_rows(
  as.data.frame(lon_vals,  na.rm = TRUE) %>% mutate(season = inputs$leaf_on$label),
  as.data.frame(loff_vals, na.rm = TRUE) %>% mutate(season = inputs$leaf_off$label)
) %>%
  dplyr::select(-x, -y, -cell) %>%
  pivot_longer(cols = -c(season, class_name, ID),
               names_to = "variable", values_to = "value") %>%
  mutate(species = unname(species_lab[as.character(class_name)])) %>%
  dplyr::select(-class_name) %>%
  filter(variable %in% lvl) %>%
  mutate(variable = factor(variable, levels = lvl))

# Wide format with paired values and their difference
df_wide <- df_long %>%
  pivot_wider(names_from = season, values_from = value) %>%
  rename(leaf_on  = `leaf on`,
         leaf_off = `leaf off`) %>%
  mutate(diff = leaf_on - leaf_off)

# One row per sample pixel and metric, otherwise the pivot did not pair up
stopifnot(nrow(df_wide) == nrow(pix_valid) * length(lvl))


# 03 - descriptive summaries
#---------------------------

# Per season and variable
df_summary_long <- df_long %>%
  group_by(season, variable) %>%
  dplyr::summarise(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE),
                   min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE),
                   .groups = "drop")

print(df_summary_long, n = Inf)

# Paired, per variable
df_summary_wide <- df_wide %>%
  group_by(variable) %>%
  dplyr::summarise(mean_lon  = mean(leaf_on,  na.rm = TRUE),
                   sd_lon    = sd(leaf_on,    na.rm = TRUE),
                   mean_loff = mean(leaf_off, na.rm = TRUE),
                   sd_loff   = sd(leaf_off,   na.rm = TRUE),
                   .groups = "drop") %>%
  mutate(diff = mean_lon - mean_loff)

print(df_summary_wide, n = Inf)


# 04 - statistical tests
#-----------------------

# Test for normality of the paired differences
shapiro <- run_shapiro(df_wide)

print(shapiro, n = Inf)
any(shapiro$p.value > 0.05)

# Differences are not normally distributed -> paired Wilcoxon test
# instead of a paired t-test.
results <- run_wilcox_subsets(df_wide)

invisible(lapply(results, function(df) print(df, n = Inf)))


#--- export test results -------------------------------------------------------

# Normality tests
write.csv(
  shapiro,
  file.path(stats_dir, paste0("shapiro_structural_", PPM, ".csv")),
  row.names = FALSE
)

# Wilcoxon results: one combined table with a subset column
write.csv(
  bind_rows(results, .id = "subset"),
  file.path(stats_dir, paste0("wilcoxon_structural_", PPM, ".csv")),
  row.names = FALSE
)

# Descriptive summaries
write.csv(
  df_summary_wide,
  file.path(stats_dir, paste0("summary_paired_structural_", PPM, ".csv")),
  row.names = FALSE
)

# Same content as RDS, keeps factor levels and full precision for reuse in R
saveRDS(results, file.path(stats_dir, paste0("wilcoxon_structural_", PPM, ".rds")))
saveRDS(df_long, file.path(stats_dir, paste0("df_long_structural_",  PPM, ".rds")))
saveRDS(df_wide, file.path(stats_dir, paste0("df_wide_structural_",  PPM, ".rds")))

message("test results written to: ", stats_dir)
