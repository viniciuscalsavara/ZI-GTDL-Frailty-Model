# Zero-Inflated GTDL Frailty Model

## Overview

This repository contains the R code supporting the manuscript:

**A zero-inflated non-proportional hazards cure frailty model: An application to credit risk data**

**Authors:** Vinicius F. Calsavara, Amanda B. G. Milani, Heloiza O. Souza, Vera L. D. Tomazella, and Eder A. Milani.

The proposed zero-inflated generalized time-dependent logistic (ZI-GTDL) frailty model provides a likelihood-based framework for survival data characterized by:

- Time-zero events (including zero-adjusted observations).
- A potential long-term non-susceptible fraction (cure fraction).
- Non-proportional hazards.
- Unobserved heterogeneity modeled through gamma frailty.

The model allows covariates to affect the probability of a time-zero event, the event-time distribution, and the time-effect parameter.

## Repository Structure

The R scripts are organized into two main directories:

### `R-codes/Simulations/`

Contains the code used to evaluate the statistical properties of the proposed model through Monte Carlo simulation studies.

The simulation studies investigate:

1. Finite-sample performance of maximum likelihood estimators, including bias, RMSE, and 95% confidence interval coverage.
2. Estimation of the long-term non-susceptible fraction under different levels of unobserved heterogeneity.
3. Type I error and statistical power of likelihood ratio tests comparing survival functions between two groups.

The directory also includes supporting functions, parallel computation scripts, and code for generating simulation figures.

### `R-codes/Application/`

Contains the R scripts used for the credit-recovery application, including:

- Maximum likelihood estimation of the proposed model.
- Estimation and interpretation of model parameters.
- Comparison with alternative survival models.
- Supporting statistical analyses and visualization.

## Software

The analyses were developed using **R**.

Users should install the packages required by the relevant scripts before running the analyses.

## Reproducibility

The simulation scripts provide the computational procedures used to evaluate the proposed methodology.

The application scripts document the statistical modeling approach used in the credit-recovery analysis.

The original credit-recovery dataset is not publicly available.
