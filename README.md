# Leaf-on vs. Leaf-off ALS

Comparison of leaf-on vs. leaf-off airborne laser scanning (ALS) data for forestry applications.

## Abstract

Airborne laser scanning (ALS) campaigns of forests differ in canopy phenology, yet the metrics derived from them are increasingly compared and combined across acquisitions. We quantified how leaf condition (i.e. leaf-on vs leaf-off) affects area-based structural metrics, and whether pulse density modulates it. Two acquisitions of the same forest in Central Germany, six months apart and harmonised to a common pulse density and scan angle, were compared for 24 conventional and eight structural complexity metrics across 400 randomly selected raster cells stratified by forest type. Differences were assessed with paired Wilcoxon tests, rank-biserial effect sizes, pixel level rank correlations, and mixed-effects models at three pulse densities. Sensitivity increased with penetration depth into the stand, from a median difference of 0.038 m for maximum height to 1.66 m for the 10th height percentile and 4.75 m in deciduous stands, whereas effects in evergreen coniferous stands were an order of magnitude smaller. Effect size alone proved a poor criterion of robustness: mean canopy height differed strongly yet preserved the ranking of pixels, amounting to a systematic offset that can be corrected, whereas the fractal box dimension shifted only slightly but reordered the ranking. Seasonal differences were largely independent of pulse density. Maximum height was the only metric largely unaffected by leaf condition. Metrics describing the lower canopy, point densities and the box dimension should not be compared across leaf conditions in deciduous stands. Robustness thus follows from how a metric is calculated rather than from the structural property it describes.

## Citation

TODO: add citation once the manuscript is published (authors, journal, DOI).

<img src = "docs/pc_plot_transect.png" width = "800" height = "400">
*Cross sections of leaf-on and leaf-off point clouds at different points per squaremeter (ppm).*

## Data involved

- Leaf-on (September 2023) and leaf-off (March 2024) ALS data, acquired with a Riegl VQ-780 II-S / Riegl VQ-780i, downsampled three point densities (4, 10 and 20 pulses/m²).

## Data availability

The point clouds are **not** included in this repository or in the GRO deposit. TODO: state the reason here (e.g. file size, licensing of the flight campaign, availability on request) so users know whether to expect them elsewhere.

All derived data (metric rasters, metadata, and the validated pixel sample) are archived at GRO: TODO add GRO link.

After downloading, place the contents so that they match the structure `src/setup.R` creates automatically:

```
data/
├── raw_data/
│   ├── pc_leafon_2023/          # not included, see above
│   └── pc_leafoff_2024/         # not included, see above
├── processed_data/
│   ├── pc_leafon_2023/{ppm4,ppm10,ppm20}/  # not included, see above
│   ├── pc_leafoff_2024/{ppm4,ppm10,ppm20}/  # not included, see above
│   ├── metrics_leafon_2023/{ppm4,ppm10,ppm20}/
│   └── metrics_leafoff_2024/{ppm4,ppm10,ppm20}/
└── metadata/
    ├── tree_species/
    └── dsm/
```

You do not need to create these folders by hand. Running any script (they all source `src/setup.R` first) builds the full structure automatically; you only need to copy the downloaded files into it.

## Workflow

Run the scripts in `scripts/` in numerical order. Each script sources `src/setup.R` first, which builds the required folder structure and loads all dependencies.

Scripts `01-04` require the raw point clouds, which are not distributed (see Data availability above). If you only have the GRO deposit, start from `script 05`.

| # | Script | What it does |
|---|---|---|
| 01 | `pc_summary.R` | Clips the leaf-on and leaf-off point cloud catalogs to their common extent, height-normalizes them, and summarises both with a `lasR` pipeline. Results are stored as RDS per season. |
| 02 | `pc_summary_viz.R` | Visualizes the point-cloud-level comparison from script 01: a point count table, a return statistics table with conditional return probabilities, and mirrored height histograms for leaf-on and leaf-off. |
| 03 | `calculate_metrics.R` | Calculates the 24 ABA metrics from the ALS data using `lidR::pixel_metrics()`, for all season and point density combinations. |
| 04 | `calculate_structural_metrics.R` | Calculates the structural complexity metrics (CHM_mean, Canopy_cover, Rumple, ENL0D, ENL2D, Box_dimension, VCI) using `lidR::catalog_apply()`, for all season and point density combinations. |
| 05 | `pixel_selection.R` | Interactive script for selecting the validated pixel sample. Plots the extent of randomly selected 10 x 10 m pixels against a higher-resolution (0.5 m) surface model and lets the user mark each pixel as valid. Produces the balanced sample of 400 pixel locations used by scripts 06-10. |
| 06 | `compare_metrics.R` | Compares the ABA metrics between leaf-on and leaf-off at pixel level, on the validated sample. Produces descriptive summaries, Shapiro tests and paired Wilcoxon tests. |
| 07 | `compare_metrics_viz.R` | Figures for the leaf-on vs. leaf-off ABA metric comparison. Reads the pixel table written by script 06. |
| 08 | `compare_structural_metrics.R` | Same comparison as script 06, for the structural metrics from script 04. |
| 09 | `compare_structural_metrics_viz.R` | Figures for the structural metric comparison. One panel per metric, kernel density per season with a leaf-on vs. leaf-off scatter inset. |
| 10 | `compare_point_density.R` | Compares metrics across leaf-on / leaf-off conditions and point density (4, 10, 20 pulses/m²) on the validated sample. Fits one linear mixed model per metric with pixel as a random intercept and a random leaf-on / leaf-off slope, tests model assumptions on the residuals, and writes the type III ANOVA and diagnostic tables. Runs on either metric set (`METRIC_SET <- "aba"` or `"struct"`). |

Scripts 07, 09 and 10 depend on their corresponding upstream script having been run first for the same settings (see the `Requires:` note in each script header).

## Setup

1. Clone this repository and open `Repo_LeafOn_vs_LeafOff_ALS.Rproj` in RStudio.
2. Download the derived data from GRO and place it as described above.
3. Run any script in `scripts/`. The first call to `source('src/setup.R')` will install all required R packages automatically (from CRAN, GitHub and GitLab, see `src/setup.R` for the full list) and create the `data/` and `output/` subfolders.
4. Run the scripts in `scripts/` in numerical order (see Workflow above).

## Repository structure

| Path | Contents |
|---|---|
| `scripts/` | Numbered analysis pipeline (`01_` to `10_`), see Workflow above. |
| `src/` | Shared setup, helper functions and metric definitions sourced by the scripts. |
| `docs/` | Figures used in this README. |
| `data/`, `output/` | Not tracked in git, created automatically by `src/setup.R` on first run. |

## Outputs

Figures are written to `output/figures/`, statistical results (ANOVA tables, diagnostics) to `output/stats/`, both as CSV and RDS.

## License

GPL-3.0, see `LICENSE`.

## Contact

Svenja Dobelmann, svenja.dobelmann@hawk.de

