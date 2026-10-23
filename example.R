# Entirely synthetic illustration, unrelated to individual CSL records.
# Run from this directory with source("example.R") or Rscript example.R.
source("elw_core.R")
set.seed(20261023)
n <- 1000L
x <- rnorm(n)
ps <- plogis(-1 + 0.5 * x)
treatment <- rbinom(n, 1, ps)
outcome <- 3200 + 120 * x - 150 * treatment + rnorm(n, sd = 450)
fit <- elw(outcome, treatment, propensity_score = ps)
print(fit)
cat("DTE at 2500 g:\n")
print(elw_dte(fit, 2500))
cat("QTEs at 0.10, 0.25, 0.50, 0.75:\n")
print(elw_qte(fit, c(0.10, 0.25, 0.50, 0.75)))
cat("PTE below 2500 g:\n")
print(elw_pte(fit, c(-Inf, 2500)))
stopifnot(all(vapply(fit$arms, function(a)
  abs(sum(a$probability_mass) - 1) < 1e-10 &&
  a$calibration_error < 1e-7, logical(1))))
