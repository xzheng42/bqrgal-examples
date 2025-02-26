## Generalized asymmetric Laplace distributions for Bayesian quantile regression

This repository contains the R package [**bqrgal**](https://github.com/xzheng42/bqrgal-examples-ba-2024/tree/main/bqrgal) (currently developer's version) 
and R scripts to reproduce the numerical results in

Yan, Y., Zheng, X., and Kottas, A. (2025). A new family of error distributions for Bayesian quantile regression. *Bayesian Analysis*. [DOI: 10.1214/25-BA1507](https://doi.org/10.1214/25-BA1507).

### Installing and using the **bqrgal** package

You can install the package with **devtools**
```
devtools::install_github("xzheng42/bqrgal-examples-ba-2024/", subdir = "bqrgal")
library(bqrgal)
```

Main functions of the package are `bgal` and `predict.bgal`:

- `bgal` fits a Bayesian linear quantile regression with generalized asymmetric Laplace (GAL) errors via Markov chain Monte Carlo (MCMC).
- `predict.bgal` (or simply `predict`) generates posterior predictive samples from a fitted GAL-based quantile regression model object with new covariates.

Notes: The current version was tested on Fedora Linux 38 under R version 4.3.3.

### Workflow to reproduce numerical results

R scripts to reproduce results in the paper are available in 
[*data-examples*](https://github.com/xzheng42/bqrgal-examples-ba-2024/tree/main/data-examples):

- Simulation study: `sim/run_all_rscripts.R` (Section 4.3 and Supplementary Material Section B; n = 100).
- Real data examples: `real_boston.R` (Section 5.1) and `real_immunoglobulinG.R` (Section 5.2).
