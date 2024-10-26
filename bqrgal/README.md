## Generalized asymmetric Laplace distributions for Bayesian quantile regression

This is the R package **bqrgal** (currently developer's version) for the paper:

Yan, Y., Zheng, X., and Kottas, A. (2024). A new family of error distributions for Bayesian quantile regression. UCSC Department of Statistics Technical Report UCSC-SOE-24-01 [[link](https://tr.soe.ucsc.edu/sites/default/files/technical-reports/UCSC-SOE-24-01.pdf)]

You can install the package with **devtools**
```
devtools::install_github("xzheng42/bqrgal-examples-ba-2024/", subdir = "bqrgal")
```

Main functions of the package are `bgal` and `predict.bgal`:

- `bgal` fits a Bayesian linear quantile regression with generalized asymmetric Laplace (GAL) errors via Markov chain Monte Carlo (MCMC).
- `predict.bgal` (or simply `predict`) generates posterior predictive samples from a fitted GAL-based quantile regression model object with new covariates.

Detailed guidelines for using the functions are referred to their help pages in R.
R scripts to reproduce results in the paper are available in [*data-examples/*](https://github.com/xzheng42/bqrgal-examples-ba-2024/tree/main/data-examples)
with instructions available in [*bqrgal-examples-ba-2024*](https://github.com/xzheng42/bqrgal-examples-ba-2024/).

Notes: The current version was tested on Fedora Linux 38 under R version 4.3.3.
