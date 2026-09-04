# Empirical likelihood weighting for counterfactual distributions

R code for empirical likelihood weighting (ELW) and the CSL-calibrated simulation design.

## ELW method

`elw_core.R` implements ELW without package dependencies.

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

Use `selection_probability` instead of `propensity_score` when both selection probabilities
are available.

## Simulation DGP

`simulation_dgp.R` contains the CSL-calibrated data-generating functions used in the paper.

```r
source("simulation_dgp.R")

calibration <- make_csl_simulation_base(private_csl_data)
simulated_data <- generate_csl_simulation_data(
  calibration,
  n = 10000,
  seed = 1
)
```

The functions cover all simulation settings reported in the paper.

## Data

The de-identified CSL data are available by application through
[NICHD DASH](https://dash.nichd.nih.gov/study/2331). They are not redistributed here under
the DASH Data Use Agreement. An authorized local copy is required to calibrate the DGP.
