# Empirical likelihood weighting for counterfactual distributions

R code for ELW calibration, distributional effect estimation, application analyses,
and simulation studies. Supported tasks and entry points are listed below.
CSL data require authorized access and are not distributed with the code.

## Files

- `elw_core.R`: ELW weights, counterfactual CDFs and quantiles, DTE, QTE, and PTE.
- `example.R`: a small, entirely synthetic example. No data files are required or saved.
- `application.R`: self-contained application workflow, with no participant data.
- `simulation/`: separate drivers for estimation, prediction, testing, and sensitivity experiments.

Only base R is required for the core method and example. Run the example from this directory:

```r
source("example.R")
```

For your own data, source `elw_core.R` and use:

```r
fit <- elw(outcome = Y, treatment = T, propensity_score = ps)
elw_dte(fit, y = 2500)
elw_qte(fit, probs = c(0.10, 0.25, 0.50, 0.75))
elw_pte(fit, interval = c(-Inf, 2500))
```

`T` must be coded 0/1 and `ps` contains P(T = 1 | X). Alternatively,
`selection_probability` accepts two columns containing the selection
probabilities for potential outcomes 0 and 1, identified under the user's
assumptions.

DTE uses P(Y(1) <= y) - P(Y(0) <= y). PTE uses a lower-inclusive,
upper-exclusive interval, so `c(-Inf, 2500)` corresponds to birthweight
strictly below 2500. For strictly positive birthweights, `c(0, 2500)`
represents [0, 2500), equivalent to the clinical interval (0, 2500).
Quantiles use the inverse weighted empirical CDF.

## Simulations

Scripts are grouped by experiment:

- `simulation/point_estimation.R`: CE-L and CE-NL point-estimation comparisons.
- `simulation/ods_estimation.R`: outcome-dependent selection comparisons.
- `simulation/conformal/`: conformal coverage experiments.
- `simulation/bootstrap_testing.R`: bootstrap rejection-rate experiments.
- `simulation/gbm/`: application-calibrated GBM interval assessment.
- `simulation/known-bound/`: sensitivity assessment with a known odds-ratio bound.
- `simulation/shared/`: shared functions for GBM and sensitivity experiments.

Helpers next to the drivers contain experiment-specific implementations and are not
additional experiments. These implementations are preserved to retain their numerical settings.

The CSL-calibrated experiments require an authorized, locally held cleaned CSL file.
Set `CSL_SIM_DATA` to its path. No CSL records, fitted calibration objects, or results
are included. Packages used by these scripts include `gbm`, `MASS`, `rpart`,
`quantreg`, `ranger`, `kernlab`, `nnet`, `mgcv`, `ggplot2`, `patchwork`,
`dplyr`, and `tidyr`; requirements depend on the selected experiment.

For example, run a small point-estimation check from this directory:

```sh
CSL_SIM_DATA="/path/to/authorized_clean_data.csv" SIM_R=2 SIM_N=1000 SIM_CORES=1 SIM_OUT="/path/to/new_results" Rscript simulation/point_estimation.R
```

The conformal driver uses `EXACT_SIM_OUTPUT`, `EXACT_SIM_REPS`, `EXACT_SIM_N`,
and `EXACT_SIM_CORES`; its observed-range summary is generated separately with
`postprocess_observed_range.R` using the same `EXACT_SIM_OUTPUT`.
The GBM experiment first needs `gbm/prepare_calibration.R`; its `run.R` requires
an explicit `--run`. Use the environment-variable settings in each driver for
formal runs. Defaults can be computationally intensive; no experiment runs when
loading `elw_core.R`.


## Application: one code file, restricted data

`application.R` includes its helper functions, so no other R source files are
needed for its supported tasks. Supply an authorized cleaned CSL CSV yourself;
CSL data cannot be distributed with this repository.

```sh
Rscript application.R --list
Rscript application.R propensity /path/to/authorized_clean_data.csv /path/to/local_output --run
Rscript application.R distribution /path/to/authorized_clean_data.csv /path/to/local_output --run
```

Tasks: `propensity`, `descriptive`, `distribution`, `heterogeneity`,
`ps_robustness`, `weight_stress`, `msm`, `lmz`, and `lmz_gamma`.
Run `propensity` before `descriptive` or `ps_robustness`, and reuse the same
local output directory for related tasks. Without `--run`, the command prints
its plan only. Loading with `source("application.R")` does not run an analysis.
The input uses the cleaned schema: `row_id`, `treatment_smoking` (0/1),
`outcome_birthweight`, `Smoke` (1/2), `Birthweight`, `Momrace`, `Marital`,
`Momage`, `Parity`, `Hxpreterm`, `Height`, `PreBMI`, `Diabetes`, `Hxasthma`,
`HIV`, `Alcohol`, and `Sex` (1=female, 2=male).
The application additionally uses `writexl`; optional missing-data branches
are not selected by the supported entry point. Bootstrap tasks can take a long time.
The single-file application supports the tasks listed above. Conformal simulation
drivers are provided separately in `simulation/conformal/`.
