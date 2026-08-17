#-------------------------------------------------------------------------------
# Name:         compare_point_density.R
# Description:  Compare the metrics across leaf-on / leaf-off conditions and
#               point densities at pixel level, on a balanced sample of
#               validated pixels. Fits one linear mixed model per metric, with
#               the sample pixel as random intercept to account for the repeated
#               measurements, and checks the model assumptions on the residuals.
#               Runs on either metric set, see METRIC_SET below.
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#
# Requires:     calculate_metrics.R and calculate_structural_metrics.R must have
#               been run for all season and density combinations.
#-------------------------------------------------------------------------------

source('src/setup.R', local = TRUE)


# 01 - configuration
#-------------------

# Metric set to analyse: "aba" for the 24 ABA metrics, "struct" for the eight
# structural complexity metrics. Run the script once per set.
METRIC_SET <- "aba"

metric_sets <- list(
  aba    = list(suffix = "metrics",        tag = "aba"),
  struct = list(suffix = "struct_metrics", tag = "structural")
)

stopifnot(METRIC_SET %in% names(metric_sets))
set_cfg <- metric_sets[[METRIC_SET]]

# seasons and point densities (points per square meter)
seasons <- list(
  solling23_lon  = list(dir = "metrics_leafon_2023"),
  solling24_loff = list(dir = "metrics_leafoff_2024")
)

ppms <- c("ppm20", "ppm10", "ppm4")

stats_dir <- file.path(output_dir, "stats")
fig_dir   <- file.path(output_dir, "figures")


# Sample of validated pixels
SAMPLE_FILE <- "sample_selection_n400_balanced.csv"

# Tree species classes of the template raster
species_lab <- c("1" = "deciduous", "2" = "coniferous")

# Columns that are design or diagnostic layers, not stand attributes
DESIGN_VARS <- c("point_density", "Zmax", "n_scales")

# Redundant metrics (e.g. SD vs VAR, percentile series) are dropped 
COR_CUTOFF <- 0.9


# 02 - data preparation
#----------------------

#' Path of the metric raster for one season / point density combination
metrics_file <- function(season, ppm) {
  file.path(processed_data_dir, seasons[[season]]$dir, ppm,
            paste0(season, "_", ppm, "_", set_cfg$suffix, ".tiff"))
}

# Validated sample pixels, read once and shared by all combinations.
# The CRS comes from the first metric raster; the tree species class is already
# a band in those rasters, so no separate template is needed.
pix_valid <- read.csv(file.path(metadata_dir, SAMPLE_FILE)) %>%
  filter(status == "valid")

ref_crs <- as.character(crs(rast(metrics_file(names(seasons)[1], ppms[1]))))

pts <- vect(pix_valid, geom = c("x", "y"), crs = ref_crs)

#' Extract the sample pixels from one season / point density combination
#'
#' The metric rasters already sit on the template grid, aligned in
#' calculate_metrics.R, so no extend or resample is needed here. Extracting
#' point values reads only the tiles the points fall into, which is far faster
#' than warping the whole raster for 400 pixels.
#'
#' @param season Character. Name of the season, must be one of names(seasons).
#' @param ppm Character. Point density, e.g. "ppm20".
#'
#' @return Data frame with one row per sample pixel and one column per metric,
#'   plus the grouping columns ID, leaf_con, dens and species.
read_metrics <- function(season, ppm) {
  
  r <- rast(metrics_file(season, ppm))
  
  terra::extract(r, pts) %>%
    as.data.frame(na.rm = TRUE) %>%
    mutate(
      leaf_con = season,
      dens     = ppm,
      species  = unname(species_lab[as.character(class_name)])
    ) %>%
    dplyr::select(-class_name)
}

# Run over all season / density combinations
results <- list()

for (season in names(seasons)) {
  for (ppm in ppms) {
    name <- paste(season, ppm, sep = "_")
    results[[name]] <- read_metrics(season = season, ppm = ppm)
  }
}

df_all <- bind_rows(results) %>%
  mutate(
    ID       = factor(ID),
    leaf_con = factor(leaf_con),
    dens     = factor(dens, levels = c("ppm4", "ppm10", "ppm20")), # ordered
    species  = factor(species),
    group    = interaction(leaf_con, dens)   # used for the variance check
  )

cat("metric set:", METRIC_SET, "\n")
cat("leaf_con levels:\n"); print(levels(df_all$leaf_con))  # exactly two
cat("dens levels:\n");     print(levels(df_all$dens))
cat("species levels:\n");  print(levels(df_all$species))

# Every pixel must carry all season / density combinations
stopifnot(all(table(df_all$ID) == length(seasons) * length(ppms)))


# 03 - metric selection
#----------------------

metric_names <- df_all %>%
  dplyr::select(where(is.numeric), -any_of(DESIGN_VARS)) %>%
  names()

cor_mat <- cor(df_all[metric_names], use = "pairwise.complete.obs")

# Report strongly correlated pairs before dropping any of them
cm <- cor_mat
diag(cm) <- NA
high <- which(abs(cm) > COR_CUTOFF, arr.ind = TRUE)
high <- high[high[, 1] < high[, 2], , drop = FALSE]

if (nrow(high) > 0) {
  cat("\ncorrelations >", COR_CUTOFF, ":\n")
  for (i in seq_len(nrow(high))) {
    cat("  ", rownames(cor_mat)[high[i, 1]], "<->",
        colnames(cor_mat)[high[i, 2]], ":",
        round(cor_mat[high[i, 1], high[i, 2]], 3), "\n")
  }
}

drop_vars  <- caret::findCorrelation(cor_mat, cutoff = COR_CUTOFF,
                                     names = TRUE, exact = TRUE)
metric_use <- setdiff(metric_names, drop_vars)

cat("\nmetrics dropped as redundant:\n"); print(drop_vars)
cat("\nmetrics used:\n");                 print(metric_use)


# 04 - linear mixed models
#-------------------------

# Sum contrasts for clean type III tests
options(contrasts = c("contr.sum", "contr.poly"))

#' Fit one mixed model per metric
#'
#' fixed:  leaf_con * dens * species
#' random: intercept per sample pixel, accounts for the repeated measurements
#'
#' @param v Metric name.
#'
#' @return merMod object.
fit_one <- function(v) {
  f <- as.formula(paste0("`", v, "` ~ leaf_con * dens * species + (1 + leaf_con|ID)"))
  lmerTest::lmer(f, data = df_all)
}

models <- setNames(lapply(metric_use, fit_one), metric_use)

#' Type III ANOVA table of one model
anova_one <- function(v) {
  a <- as.data.frame(car::Anova(models[[v]], type = 3, test.statistic = "F"))
  a$metric   <- v
  a$effect   <- rownames(a)
  rownames(a) <- NULL
  a
}

results_tab <- do.call(rbind, lapply(metric_use, anova_one))

# Correct for multiple testing across the metrics, separately per effect
p_col <- grep("^Pr", names(results_tab), value = TRUE)[1]

results_tab <- results_tab %>%
  filter(effect != "(Intercept)") %>%
  group_by(effect) %>%
  mutate(
    p_fdr  = p.adjust(.data[[p_col]], method = "BH"),
    signif = case_when(
      p_fdr < 0.01 ~ "***",
      p_fdr < 0.05 ~ "**",
      TRUE         ~ ""
    )
  ) %>%
  ungroup() %>%
  arrange(effect, p_fdr)

print(results_tab, n = Inf)

# 05 - model assumptions
#-----------------------
# Checked on the fitted models, not on the raw values:
#   - residuals approximately normal and homoscedastic
#   - random intercept variance not collapsed to zero (singular fit)
# Shapiro-Wilk is oversensitive at this sample size and will flag negligible
# deviations, so it is read together with the skewness and the plots below.

#' Assumption diagnostics for one model
#'
#' @param v Metric name.
#'
#' @return One-row data frame with the diagnostics.
diagnose_one <- function(v) {
  
  m <- models[[v]]
  r <- residuals(m)
  
  # Rows actually used by the model; lmer drops rows with NA in the response,
  # so df_all$group cannot be used directly
  used <- as.integer(rownames(model.frame(m)))
  grp  <- droplevels(df_all$group[used])
  
  # Random intercept share of the total variance
  vc  <- as.data.frame(VarCorr(m))
  v_id   <- vc$vcov[vc$grp == "ID"]
  v_res  <- vc$vcov[vc$grp == "Residual"]
  
  # Equal residual variance across the season / density groups
  lev <- car::leveneTest(r ~ grp)
  
  data.frame(
    metric        = v,
    n_obs         = length(r),
    shapiro_p     = shapiro.test(r)$p.value,
    skewness      = mean((r - mean(r))^3) / stats::sd(r)^3,
    levene_p      = lev$`Pr(>F)`[1],
    icc           = v_id / (v_id + v_res),
    singular      = isSingular(m),
    max_abs_scaled = max(abs(r / stats::sd(r)))
  )
}

diagnostics <- do.call(rbind, lapply(metric_use, diagnose_one)) %>%
  mutate(
    resid_normal = shapiro_p >= 0.05,
    var_equal    = levene_p  >= 0.05
  )

print(diagnostics)

cat("\nsingular fits:", sum(diagnostics$singular), "/", nrow(diagnostics), "\n")
cat("residuals normal (Shapiro p >= 0.05):",
    sum(diagnostics$resid_normal), "/", nrow(diagnostics), "\n")
cat("homoscedastic (Levene p >= 0.05):",
    sum(diagnostics$var_equal), "/", nrow(diagnostics), "\n")

# Residual plots for visual inspection: one page per metric
pdf(file.path(fig_dir, paste0("lmm_diagnostics_", set_cfg$tag, ".pdf")),
    width = 9, height = 4.5)

for (v in metric_use) {
  m <- models[[v]]
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(fitted(m), residuals(m), main = v,
       xlab = "fitted", ylab = "residuals", pch = 16, cex = 0.4,
       col = grDevices::adjustcolor("black", alpha.f = 0.3))
  abline(h = 0, col = "red", lty = 2)
  qqnorm(residuals(m), main = "", pch = 16, cex = 0.4,
         col = grDevices::adjustcolor("black", alpha.f = 0.3))
  qqline(residuals(m), col = "red", lty = 2)
  par(op)
}

dev.off()


#--- export results ------------------------------------------------------------

write.csv(
  results_tab,
  file.path(stats_dir, paste0("lmm_point_density_", set_cfg$tag, ".csv")),
  row.names = FALSE
)

write.csv(
  diagnostics,
  file.path(stats_dir, paste0("lmm_diagnostics_", set_cfg$tag, ".csv")),
  row.names = FALSE
)

# Keeps factor levels and full precision for reuse in R
saveRDS(results_tab,
        file.path(stats_dir, paste0("lmm_point_density_", set_cfg$tag, ".rds")))
saveRDS(df_all,
        file.path(stats_dir, paste0("df_all_point_density_", set_cfg$tag, ".rds")))

message("test results for ", set_cfg$tag, " written to: ", stats_dir)
message("diagnostic plots written to: ", fig_dir)
# 
