# An Adaptive $L_2$-type Test for High-dimensional White Noise

This repository contains the R code for the simulation studies and real-data analysis in the accompanying paper. The simulations provide the numerical results underlying Tables 1–2 in the main paper and Tables T1–T7 in the supplementary material. The real-data scripts implement the S&P 500 analysis in Section C of the supplementary material.

Run the code with the **repository root as the R working directory**, so that relative paths such as `./src/` and `./realdata/` resolve correctly. Results are collected in `result/`; some scripts currently write new outputs to the repository root instead, as documented below.

## Repository Structure

```text
.
├── github_version.Rproj
├── src/                       # Data-generating processes and test implementations
├── simulation/                # Simulation experiments and size plots
├── realdata/                  # S&P 500 analysis scripts and downloaded input data
└── result/                    # Saved numerical results and figures
    ├── plot_for_size/
    ├── modelx_sim
    ├── modelx_sim78
    └── realdata_sp500_gfc_2005_2010_rolling_pvalues/
```

### Helper code: `src/`

- `datamodel.R`: simulation models 1–6.
- `datamodel_add.R`: simulation models 7–8.
- `test_fun.R` implements the proposed and comparison tests. `test_fun_statistic.R` extends these implementations to return test statistics and p-values for the real-data analysis.
- `wnSVAR.R` and `HDWHtest.R` implement the comparison methods of [Wang et al. (2023)](https://github.com/YingcunXia/WhiteNoiseTest) and [Tsay (2020)](https://github.com/RTsay1/HDWNtest), respectively. `wnSVAR_statistic.R` extends `wnSVAR.R` to return test statistics and p-values for the real-data analysis.
- `tool.cpp`: C++ helper routines compiled through `Rcpp::sourceCpp()`.

### Scripts and paper outputs

| Script | Purpose | Paper output |
|---|---|---|
| `simulation/sim.R` | Main Monte Carlo experiments for models 1–6 | Numerical inputs for main-paper Tables 1–2 and supplementary Tables T2–T7 |
| `simulation/sim_model78.R` | Additional experiments for models 7–8 | Numerical inputs for supplementary Table T1 |
| `simulation/size_sim.R` | Empirical-size experiments | Data for supplementary Figure F1 (simulation study) |
| `simulation/size_plot.R` | Plot the empirical-size results | Supplementary Figure F1 (simulation study) |
| `realdata/realdata_sp500_rolling_pvalues.R` | Rolling-window tests and p-value plots | Figure F2 in supplementary Section C, Real data analysis |

The figure references distinguish the simulation study from Section C. The table scripts save numerical CSV results, rather than typeset tables. See the configuration notes below before rerunning the empirical analysis.

## Requirements

Install the required R packages:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "HDTSA", "ncvreg", "TauStar",
  "data.table", "ggplot2", "dplyr", "ggnewscale", "latex2exp"
))
```

The `parallel` package is included with R. A working C++ compilation toolchain is required for `src/tool.cpp`.

The results for the method of Chang et al. (2017) were obtained using **HDTSA version 1.0.6-2**.

The proposed test is also available through `HDTSA::WN_test(X_model, k, method = "L_2")`. Both implementations use the same test statistic and bootstrap procedure. The local function `newWn_test(..., resampling = TRUE)` returns rejection decisions in $res_resampling by comparing the observed statistic with the bootstrap order statistic at index `floor(B * (1 - alpha))`, where `B` is the number of bootstrap replications. In contrast, `HDTSA::WN_test()` returns bootstrap p-values. Applying `p.value < alpha` can give different decisions at boundary cases—for example, when `p.value = alpha`. To reproduce the rejection rates reported in the paper, use $res_resampling from the local implementation.

The simulation code uses `parallel::mclapply()`. Use Linux or macOS for multicore simulations; on Windows, set `cores = 1` in the simulation function calls. The standalone size simulation determines its worker count internally, so Windows users must also set its `num_cores` to `1` before running it. The rolling analysis contains a serial fallback on Windows.

The main simulations default to 5,000 replications and 200 workers. Adjust the worker count to available CPU and memory resources, especially for the largest dimensions. Smaller replication counts are useful for trial runs but do not reproduce the paper's Monte Carlo precision.

## Working Directory and Data

Open `github_version.Rproj` in RStudio, or set the working directory explicitly:

```r
setwd("/path/to/repository")
stopifnot(dir.exists("src"), dir.exists("simulation"), dir.exists("realdata"))
```

Keep this working directory when sourcing scripts; do not change into `simulation/` or `realdata/`.

The empirical analysis requires `kp9ncdlcmmbmoaby.csv`, downloaded from **Wharton Research Data Services (WRDS)** on **2026-08-06**.

**Google Drive download link: [`kp9ncdlcmmbmoaby.csv`](https://drive.google.com/file/d/1zNFlW9GtCAwDpt5e5ORVHd2rRqc4O1Qy/view?usp=sharing)**

Download the file and place it at:

```text
realdata/kp9ncdlcmmbmoaby.csv
```

Keep the original filename. Both empirical scripts read this location. The simulations do not require this dataset.

## Running the Simulations

### Main experiments: Tables 1–2 and T2–T7

Source the script to load the model generators and define `run_simulation()`, then call the function for the desired model and settings:

```r
source("simulation/sim.R")

run_simulation(
  data_model4, "model4",
  start_ind = 10, end_ind = 12,
  reptime = 5000, cores = 200
)
```

The example runs model 4 for settings 10–12. Use `data_model1` through `data_model6` with matching names (`"model1"` through `"model6"`). Set `start_ind = 1, end_ind = 12` to run the complete grid for a model.

| Setting indices | Sample sizes `n` | Dimensions `p` |
|---|---|---|
| 1–3 | 50, 100, 200 | 5, 5, 5 |
| 4–6 | 50, 100, 200 | 7, 10, 14 |
| 7–9 | 50, 100, 200 | 50, 100, 200 |
| 10–12 | 50, 100, 200 | 2,500, 10,000, 40,000 |

Results are written to `result/<model_name>_sim/`. Files follow the pattern `res_q1_<model_name>_n<n>p<p>.csv`; `q1`, `q2`, `q3`, and `q4` correspond to lag orders 2, 4, 6, and 8. Entries are rejection percentages at the 5% significance level.

Each CSV is a snapshot of the 12-row result matrix accumulated within that function call. Unrun settings and methods skipped at large dimensions remain `NA`; separate calls do not reload earlier matrices.

### Additional experiments: Table T1

```r
source("simulation/sim_model78.R")

run_simulation_model78_only(
  data_model7, "model7",
  start_ind = 1, end_ind = 12,
  reptime = 5000, cores = 200
)

run_simulation_model78_only(
  data_model8, "model8",
  start_ind = 1, end_ind = 12,
  reptime = 5000, cores = 200
)
```

The script saves results to `result/model7_sim78/` and `result/model8_sim78/`. It uses the same setting grid and lag-file convention as the main simulation. Only `Li2019`,  and `our_resampling` are computed; the remaining method columns are `NA`.

### Empirical-size plots: supplementary Figure F1

To regenerate the plots from the CSV already stored in `result/plot_for_size/`:

```r
source("simulation/size_plot.R")
```

This produces:

- `result/plot_for_size/plot_size_n_50400.pdf`
- `result/plot_for_size/plot_size_n2_50400.pdf`

To rerun the underlying experiment:

```r
source("simulation/size_sim.R")
```

Sourcing this script immediately starts 10,000 replications per setting, using `n = 50, 100, 200, 300, 400`, dimensions `p = n` and `p = n^2`, and t-distribution degrees of freedom 3, 5, 10, and infinity (Gaussian).

`size_sim.R` saves its results to `result/plot_for_size/figure_simulation_results_new.csv`, which is read by `size_plot.R`. The simulation appends to an existing output CSV; move or remove the existing CSV before a fresh full run to avoid duplicate rows.

## Running the Real-Data Analysis

### Rolling p-values: Section C, Figure F1

After downloading the input data, run:

```r
# Example worker limit; choose a suitable value for your machine.
Sys.setenv(ROLLING_WORKERS = "4") #parallel cores
source("realdata/realdata_sp500_rolling_pvalues.R")
```

Sourcing this script immediately runs the analysis through its final `run_rolling_application(cfg)` call. The analysis covers **January 1, 2005 through December 31, 2010**, using 252-observation rolling windows, a step of 5 observations, lag order 1, and a significance level of 0.05.

Results are saved to `result/realdata_sp500_gfc_2005_2010_rolling_pvalues/`. Outputs include:

- `method_results/rolling_<method>.csv`: per-method rolling results and checkpoints.
- `combined_rolling_results.csv` and `rolling_summary.csv`: combined results and summaries.
- `selected_stock_manifest.csv`, `rolling_analysis_design.csv`, and `rolling_window_index.csv`: stock selection and analysis design.
- `figures/`: individual and combined p-value plots, including `rolling_all_methods_K1.pdf`  with PNG versions.

Compatible checkpoints are reused on subsequent runs. The final `cfg$output_dir` assignment overrides `ROLLING_OUTPUT_DIR`, so setting that environment variable alone does not redirect the sourced script's final run. Use the `source()` entry point shown above; direct execution with `Rscript` also triggers an earlier default run before the final configured run.
