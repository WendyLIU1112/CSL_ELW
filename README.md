# Empirical likelihood weighting for counterfactual distributions

This repository provides the empirical likelihood weighting (ELW) method and the
CSL-calibrated data-generating process used in the accompanying paper. It contains two R
files and no data.

## ELW method

`elw_core.R` is a self-contained implementation with no package dependencies. Its main
function is `elw()`.

```r
source("elw_core.R")

fit <- elw(
  outcome = birthweight,
  treatment = smoking,
  propensity_score = fitted_propensity
)

elw_cdf(fit, y = 2500, level = 1)
elw_dte(fit, y = 2500)
elw_qte(fit, probs = c(0.10, 0.25, 0.50, 0.75))
elw_pte(fit, interval = c(1500, 2500))
```

The `selection_probability` argument may replace `propensity_score` when both potential-outcome
selection probabilities are available. See the function documentation for its accepted
two-column and list formats.

## Simulation DGP

`simulation_dgp.R` contains the exact CSL-calibrated data-generating functions used in the
paper:

```r
source("simulation_dgp.R")

calibration <- make_csl_simulation_base(private_csl_data)
simulated_data <- generate_csl_simulation_data(
  calibration,
  n = 10000,
  seed = 1
)
```

The DGP supports the conditional-exchangeability, nonlinear treatment-assignment,
limited-overlap, heterogeneous-effect, and LMZ outcome-dependent selection settings used in
the simulation study.

## Data

The de-identified Consortium on Safe Labor data are available through the
[NICHD Data and Specimen Hub (DASH)](https://dash.nichd.nih.gov/study/2331)
under controlled access. Because redistribution is not permitted under the DASH Data Use
Agreement, neither the individual-level data nor derived individual-level files are included
in this repository. Researchers may request access directly from NICHD DASH.

The CSL-calibrated DGP requires an authorized local copy of the data to create the calibration
object.
