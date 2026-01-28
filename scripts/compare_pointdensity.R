#-------------------------------------------------------------------------------
# Name:         compare_pointdensity_BI.R
# Description:  compare Metrics for the BI plots between leaf-on and leaf-off 
#               conditions and between different point densities (20,10,4 pls/m^2)
#               using a MANOVA 
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------

# source setup script
source('src/setup.R', local = TRUE)


#### Data preparation ####
# read the leaf-off data 
loff_ppm4_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling24_loff_ppm4_indices_BI.csv')
)
loff_ppm10_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling24_loff_ppm10_indices_BI.csv')
)
loff_ppm20_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling24_loff_ppm20_indices_BI.csv')
)

# read the leaf-on data 
lon_ppm4_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling23_lon_ppm4_indices_BI.csv')
)
lon_ppm10_df <- read.csv(
  file.path(processed_data_dir, 'metrics','plt_level','solling23_lon_ppm10_indices_BI.csv')
)
lon_ppm20_df <- read.csv(
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
keep <- intersect(lon_ppm4_df$name, lon_ppm4_df$name)

## Combining the datasets
df <- bind_rows(
  lon_ppm4_df  %>% filter(name %in% keep) %>% mutate(season = "leaf on",  dens = "4ppm"),
  loff_ppm4_df %>% filter(name %in% keep) %>% mutate(season = "leaf off", dens = "4ppm"),
  lon_ppm10_df %>% filter(name %in% keep) %>% mutate(season = "leaf on",  dens = "10ppm"),
  loff_ppm10_df %>% filter(name %in% keep) %>% mutate(season = "leaf off", dens = "10ppm"),
  lon_ppm20_df %>% filter(name %in% keep) %>% mutate(season = "leaf on",  dens = "20ppm"),
  loff_ppm20_df %>% filter(name %in% keep) %>% mutate(season = "leaf off", dens = "20ppm")
) %>%
  # gather metrics
  pivot_longer(
    cols = -c(name, season, dens),
    names_to = "variable",
    values_to = "value"
  ) %>%
  # add species info
  left_join(
    ts %>% select(kspnr, dominant_species),
    by = c("name" = "kspnr")
  ) %>%
  rename(
    Tree_ID = name,
    species = dominant_species
  ) %>%
  # THIS is the key step for MANOVA
  pivot_wider(
    id_cols   = c(Tree_ID, season, dens, species),
    names_from  = variable,
    values_from = value
  ) %>%
  mutate(
    season = factor(season),
    dens   = factor(dens)
  )  %>%
  select(-Tree_ID) #ID not a variable


# Abhängige Variablen (Metriken)
dv_vars <- df %>% select(where(is.numeric))
# Unabhängige Variablen (Faktoren)
iv_vars <- df %>% select(season, dens)

cat("\nStichprobenumfang:\n")
cat("Gesamt:", nrow(df), "\n")
cat("Anzahl Metriken:", ncol(dv_vars), "\n\n")

cat("Gruppierungsvariablen:\n")
print(table(df$season, df$dens))

# Check requirements for MANOVA
#-------------------------------------------------------------------------------

# ============================================================================
#  Multivariate Normality
# ============================================================================

cat("Shapiro-Wilk Tests for single variables:\n")

# Testing for each group seperately
for(s in unique(df$season)) {
  for(d in unique(df$dens)) {
    subset_data <- df %>% 
      filter(season == s, dens == d) %>%
      select(where(is.numeric))
    
    if(nrow(subset_data) >= 3) {
      cat("Group:", s, "-", d, "(n =", nrow(subset_data), ")\n")
      
      normal_count <- 0
      for(var in names(subset_data)) {
        if(nrow(subset_data) >= 3) {
          sw_test <- shapiro.test(subset_data[[var]])
          if(sw_test$p.value >= 0.05) normal_count <- normal_count + 1
        }
      }
      cat("Normally distributed variables:", normal_count, "/", ncol(subset_data), "\n\n")
    }
  }
}

cat("NOTE: For large samples (n > 30 per group), MANOVA is robust to violations of the normal distribution assumption.\n")

# ============================================================================
# Homogeneity of covariance matrices (Box's M test)
# ============================================================================

cat("=" , rep("=", 70), "\n", sep="")
cat(" Homogeneity of covariance matrices (Box's M test)\n")
cat("=" , rep("=", 70), "\n", sep="")

# Box's M Test
# Achtung: Sehr sensitiv bei großen Stichproben!
tryCatch({
  # Kombiniere season und dens zu einer Gruppenvariable
  df$group <- interaction(df$season, df$dens)
  
  boxm_result <- boxM(dv_vars, df$group)
  print(boxm_result)
  
  cat("\nInterpretation:\n")
  cat("H0: covariance matrices are equal\n")
  if(boxm_result$p.value < 0.001) {
    cat("Box's M Test signifikant (p < 0.001)\n")
    cat("  → covariance matrices are not homogenous\n")
    cat("  → with large sample size: often not critical\n")
    cat("  → Pillai's Trace is more robust against violations\n\n")
  } else if(boxm_result$p.value < 0.05) {
    cat("Box's M Test signifikant (p < 0.05)\n")
  } else {
    cat("Kovarianzmatrizen are homogenous (p ≥ 0.05)\n\n")
  }
}, error = function(e) {
  cat("ERROR:\n")
  cat(e$message, "\n\n")
})

# ============================================================================
# No Multicollinearity of variables 
# ============================================================================

cat("=" , rep("=", 70), "\n", sep="")
cat("5. Multicollinearity\n")
cat("=" , rep("=", 70), "\n", sep="")

cor_matrix <- cor(dv_vars)
high_cor <- which(abs(cor_matrix) > 0.9 & cor_matrix != 1, arr.ind = TRUE)

cat("Correlations > 0.9:\n")
if(nrow(high_cor) > 0) {
  for(i in 1:nrow(high_cor)) {
    var1 <- rownames(cor_matrix)[high_cor[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor[i, 2]]
    cor_val <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
    if(high_cor[i, 1] < high_cor[i, 2]) {  # Nur einmal pro Paar ausgeben
      cat("  ", var1, "↔", var2, ":", round(cor_val, 3), "\n")
    }
  }
  cat("\n High correlations detected\n")
  cat("  → Consider removing redundant variables\n")
} else {
  cat("  No (all Correlations < 0.9)\n")
}

# select varibales to be removed 
remove <- c("BE_H_P20","BE_H_P50", "BE_H_MAX", "BE_H_P80","BE_H_VAR", "point_density")
dv_vars_noncorr <- dv_vars %>%
  dplyr::select(-all_of(remove))

# MANOVA 
#-------------------------------------------------------------------------------

cat("=" , rep("=", 70), "\n", sep="")
cat("MANOVA \n")
cat("=" , rep("=", 70), "\n", sep="")

# MANOVA
manova <- manova(as.matrix(dv_vars_noncorr) ~ season * dens, data = df)

# Pillai's Trace (recommended)
cat("Pillai's Trace:\n")
summary_pillai <- summary(manova, test = "Pillai")
print(summary_pillai)

cat("\n")

# Wilks' Lambda
cat("Wilks' Lambda:\n")
summary_wilks <- summary(manova, test = "Wilks")
print(summary_wilks)


# FOLLOW-UP: UNIVARIATE ANOVAS
#-------------------------------------------------------------------------------

response_vars <- colnames(dv_vars_noncorr)

for (var in response_vars) {
  
  cat("\n", rep("-", 60), "\n", sep = "")
  cat("Univariate ANOVA for:", var, "\n")
  cat(rep("-", 60), "\n", sep = "")
  
  formula_str <- paste(var, "~ season * dens")
  model <- aov(as.formula(formula_str), data = df)
  
  print(summary(model))
}

