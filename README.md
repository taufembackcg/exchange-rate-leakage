# Data Leakage and the Apparent Predictability of Exchange Rates with Statistical Learning

Code and results for the paper *Data Leakage and the Apparent Predictability of Exchange Rates with Statistical Learning: Evidence from Developed and BRICS Economies*, by Cleiton Guollo Taufemback, Hudson da Silva Torrent, Wagner Scheid Mallmann and Bernardo Cainelli Gomes da Silva (Universidade Federal do Rio Grande do Sul).

The pipeline forecasts monthly exchange rates against the US dollar for six developed countries (Canada, Denmark, France, Japan, Sweden, United Kingdom) and six BRICS economies (Brazil, China, India, Indonesia, Russia, South Africa). It uses five economic models, first-difference and error correction specifications, twelve statistical learning methods and direct forecasts at horizons 1 to 12 months. Each forecast is evaluated against the driftless random walk with the relative RMSE, the Diebold–Mariano and Clark–West tests and the Model Confidence Set.

The exercise is run five times:

| Design | What changes with respect to the leak-free design |
|---|---|
| `main` | Leak-free design: every forecast uses only information available at its origin (the results of the paper). |
| `window` | The estimation window of the forecast made at *t* ends at *t*−1 instead of *t*−*h*. |
| `nolags` | Every series enters at its reference month, with no publication lag. |
| `filters` | The Chow–Lin disaggregation of GDP and the Hamilton filter are estimated once on the full sample. |
| `scaling` | Standardization and the spline basis are computed on the full sample. |

## Contents

| File | Step |
|---|---|
| `run_pipeline.R` | Runs every step below in order, for all designs and countries (see *Running*). |
| `run_pipeline.bat` | Windows shortcut for `run_pipeline.R`. |
| `01_functions.R` | Shared functions (installs the R packages it needs on first use). |
| `02_data_BRICS_plus_US_ref.r` | 02: downloads and builds the monthly panels. |
| `03_models.R` | 03: estimates every specification with rolling windows of 120 months. |
| `spline_glmnet.R` | Cubic spline regression (SP and SPI), used by 03. |
| `04_analysis_single.R` | 04: relative RMSE, Diebold–Mariano and Clark–West tests. |
| `forecast_alignment.R` | Common evaluation window of each country and horizon, used by 04 and 06. |
| `06_mcs_single.R` | 06: Model Confidence Set (fixed bootstrap seed for each country and horizon). |
| `05_export_results.r` | 05: tables and figures of each design. |
| `aggregate_stats.R` | 07: the aggregate numbers quoted in the paper. |
| `export_tuning.R` | Optional: the tuning parameters chosen in each yearly block. |
| `results/` | The results reported in the paper (see *Results*). |

Comments inside the R scripts are partly in Portuguese.

## Requirements

- R 4.3 or later (the results were produced with R 4.3.2 on Windows).
- The packages loaded in `01_functions.R`, `04_analysis_single.R`, `05_export_results.r`, `06_mcs_single.R` and `aggregate_stats.R`, among them `caret`, `glmnet`, `randomForest`, `kernlab`, `LiblineaR`, `earth`, `gam`, `mboost`, `pls`, `tempdisagg`, `neverhpfilter`, `fredr`, `rsdmx`, `MCS`, `tidyverse`, `readxl`, `writexl`, `kableExtra` and `tikzDevice`. `01_functions.R` installs missing packages from CRAN.
- Internet access for the data downloads (IMF, FRED and DBnomics).
- Disk space for the estimated models: about 10 to 15 GB per design (60 to 75 GB for the five).

## Data

Most series are downloaded when step 02 runs. Two inputs must be placed in this folder first:

1. **FRED API key.** Request a free key at <https://fred.stlouisfed.org/docs/api/api_key.html> and save it, alone on one line, in a file named `fred_key.txt`. This file is listed in `.gitignore`.
2. **OECD monetary aggregates (M1).** From the OECD Data Explorer, download the dataflow `OECD.SDD.STES,DSD_STES@DF_MONAGG,4.0` (Monthly monetary aggregates, all data) as CSV with labels, and keep the file name `OECD.SDD.STES,DSD_STES@DF_MONAGG,4.0+all.csv`. The file used for the paper was downloaded on 20 October 2026.

The series are revised over time, so a new download can produce results that differ slightly from those in `results/`.

## Running

From this folder:

```
Rscript run_pipeline.R
```

or, on Windows, double-click `run_pipeline.bat`. Options:

```
--models-dir=DIR    folder for the estimated models (default: models)
--cores=N           cores used in step 03 (default: all but one)
--designs=a,b       subset of main,window,nolags,filters,scaling
--countries=a,b     subset of ca,dk,fr,gb,jp,se,br,cn,id,in,ru,za
--tuning            also export the yearly tuning choices
--dry-run           list the jobs without running them
```

Each step runs in its own R process and writes its log to `logs/`. A finished step is recorded in `logs/done/` and skipped when the pipeline is started again, so an interrupted run resumes where it stopped.

The full run is long. With 15 cores, the leak-free design took about 66 hours, and the other designs took a similar time each, about two weeks in total. Step 03 accounts for almost all of this time. The results of each design are written to `outputs/` (main) or `outputs_<design>/`.

## Results

`results/<design>/` holds the results reported in the paper for each design:

- `Resultados_Todos_<country>.xlsx`: relative RMSE and test p-values of every forecast series (country, model, specification, horizon and method);
- `tables/`: tables by country and horizon, Model Confidence Set p-values (`mcs_rw_pvalue_<country>.csv`) and the LaTeX tables of the Supplementary Material;
- `aggregate_stats_for_paper.txt` and `.json`: the aggregate numbers quoted in the paper.

`results/tuning/` summarizes how often LASSO and Ridge Regression chose the largest penalty of the grid. It was exported before the leakage designs were extended from horizons 1, 3, 6 and 12 to all twelve horizons, so for those designs it covers only these four horizons.

## License

The code is released under the MIT License (see `LICENSE`). The data belong to their providers (IMF, FRED, OECD and DBnomics) and are not part of this repository.
