#-------------------------------------------------------------------------------
# Name:         stats_func.R
# Description:  Shared paired tests for the leaf-on / leaf-off comparison,
#               used by compare_metrics.R (24 ABA metrics) and
#               compare_structural_metrics.R (8 structural complexity metrics).
# Author:       Svenja Dobelmann
# Contact:      svenja.dobelmann@hawk.de
#-------------------------------------------------------------------------------


#' Significance stars for a vector of p values
add_stars <- function(p) {
  case_when(
    p < 0.01 ~ "***",
    p < 0.05 ~ "**",
    p < 0.1  ~ "*",
    TRUE     ~ ""
  )
}


#' Shapiro-Wilk test on the paired differences, per variable
#'
#' @param data Wide format data frame with columns variable and diff.
#'
#' @return Tibble with statistic, p value and significance stars.
run_shapiro <- function(data) {
  data %>%
    tidyr::drop_na(diff) %>%
    group_by(variable) %>%
    dplyr::summarise(shapiro = list(shapiro.test(diff)), .groups = "drop") %>%
    rowwise() %>%
    mutate(statistic = shapiro$statistic,
           p.value   = shapiro$p.value) %>%
    ungroup() %>%
    dplyr::select(-shapiro) %>%
    mutate(signif = add_stars(p.value))
}


#' Run paired Wilcoxon tests for all variables in one subset
#'
#' The BH correction is applied within one call, i.e. over the metrics of the
#' subset passed in. Calling the function separately per species subset means
#' three correction families, which has to be stated in the methods.
#'
#' @param data Data frame in wide format with columns variable, leaf_on, leaf_off.
#' @param ratios Logical. Add the median ratio columns. Set FALSE for metric
#'   sets where a ratio of medians is not meaningful.
#'
#' @return Tibble with one row per variable: sample size, mean and median
#'   difference, medians per season, test statistic, raw and BH-adjusted
#'   p value, rank-biserial effect size and, optionally, ratios.
run_wilcox <- function(data, ratios = TRUE) {
  
  out <- data %>%
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
      # Rank-biserial correlation via effectsize; verbose = FALSE suppresses
      # one message per group
      effsize     = list(effectsize::rank_biserial(
        x = leaf_on, y = leaf_off, paired = TRUE, ci = 0.95, verbose = FALSE
      )),
      r_rb        = effsize[[1]]$r_rank_biserial,
      .groups = "drop"
    ) %>%
    dplyr::select(-test, -effsize)
  
  if (ratios) {
    out <- out %>%
      mutate(
        # Uniform ratio: > 1 = leaf-on higher, < 1 = leaf-off higher.
        # NA for signed distribution measures and denominators close to zero.
        ratio = if_else(
          str_detect(as.character(variable), "KURTOSIS|SKEW") | abs(med_off) < 1e-9,
          NA_real_,
          med_on / med_off
        ),
        # PR only: interception instead of transmission (for the text)
        int_ratio = if_else(
          str_detect(as.character(variable), "_PR_") & med_off < 1 - 1e-9,
          (1 - med_on) / (1 - med_off),
          NA_real_
        )
      )
  }
  
  out %>%
    mutate(
      # FDR correction (Benjamini-Hochberg) across all metrics
      p_adj  = p.adjust(p.value, method = "BH"),
      signif = add_stars(p_adj)
    )
}


#' Run the paired tests for all species subsets
#'
#' @param data Wide format data frame with a species column.
#' @param ... Passed on to run_wilcox().
#'
#' @return Named list with the subsets all, coniferous and deciduous.
run_wilcox_subsets <- function(data, ...) {
  list(
    all        = run_wilcox(data, ...),
    coniferous = run_wilcox(filter(data, species == "coniferous"), ...),
    deciduous  = run_wilcox(filter(data, species == "deciduous"), ...)
  )
}
