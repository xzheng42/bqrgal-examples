rm(list = ls())

set.seed(42)

inpath <- "path/to/output/folder/"
outpath <- "path/to/output/folder/"

library(bqrgal)

#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
load(paste0(inpath, "train_test_data.RData"))

nn <- true_params$nn
pp <- length(true_params$bb)
p0_vals <- true_params$p0_vals
train_kk <- true_params$train_kk

#-------------------------------------------------------------------------------
# Model settings
#-------------------------------------------------------------------------------
priors <- list("intercept_gaus_var" = 100,
               "sigma_invgamma" = c(2, 2),
               "eta_gamma" = c(0.1, 0.1))
starting <- list(be0 = 0, sigma = 1, vv = 1 * rexp(nn),
                 ss = TruncatedDistributions::rtnorm(nn, 0, 1, 0, Inf),
                 omega = rep(1, pp))
tuning <- list("step_size" = 0.01)
mcmc_settings <- list(n_iter = 150000, n_burn = 50000, n_thin = 20, n_report = 2000)

#-------------------------------------------------------------------------------
# Parallel computation settings
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Fit GAL models...set up multiple sessions for parallel computing\n")
cat("-------------------------------------------------------------------------\n")

suppressPackageStartupMessages({
  library(doFuture)
  library(flexiblas)
})

registerDoFuture()

ncores <- parallel::detectCores() - 1

num_of_workers <- min(ncores, train_kk)
plan(multisession, workers = num_of_workers)

# Switch the BLAS backend to NETLIB (single-thread BLAS) in each session
foreach(i = 1:num_of_workers) %dopar% {
  backend <- "NETLIB"
  idx <- flexiblas_load_backend(backend)
  flexiblas_switch(idx)
}

#-------------------------------------------------------------------------------
# Fit the model to normal data
#-------------------------------------------------------------------------------

cat("-------------------------------------------------------------------------\n")
cat("Fitting the model to normal data\n\n")

gal_normal_runtime <- system.time(
  gal_normal_out_list <- foreach(l = seq_along(p0_vals), .options.future = list(packages = "bqrgal", seed = TRUE)) %:%
    foreach(j = 1:train_kk) %dofuture% 
    {
      source("path/to/data-examples/run_gal_mcmc.R")
      run_gal_mcmc(normal_train_data_list[[l]][[j]], train_covar_list[[l]][[j]], p0_vals[l],
                   "laplace", priors, starting, tuning, mcmc_settings, "slice", FALSE)
    }
)

cat("done\n")
cat(paste0("Running time: ", round(gal_normal_runtime[3]/60), " minutes\n"))
cat("-------------------------------------------------------------------------\n")

#-------------------------------------------------------------------------------
# Fit the model to laplace data
#-------------------------------------------------------------------------------

cat("-------------------------------------------------------------------------\n")
cat("Fitting the model to laplace data\n\n")
gal_laplace_runtime <- system.time(
  gal_laplace_out_list <- foreach(l = seq_along(p0_vals), .options.future = list(packages = "bqrgal", seed = TRUE)) %:%
    foreach(j = 1:train_kk) %dofuture% 
    {
      source("~/bqr/code/data-examples/for-github/run_gal_mcmc.R")
      run_gal_mcmc(laplace_train_data_list[[l]][[j]], train_covar_list[[l]][[j]], p0_vals[l],
                   "laplace", priors, starting, tuning, mcmc_settings, "slice", FALSE)
    }
)

cat("done\n")
cat(paste0("Running time: ", round(gal_laplace_runtime[3]/60), " minutes\n"))
cat("-------------------------------------------------------------------------\n")

#-------------------------------------------------------------------------------
# Fit the model to lgpd data
#-------------------------------------------------------------------------------

cat("-------------------------------------------------------------------------\n")
cat("Fitting the model to lgpd data\n\n")

gal_lgpd_runtime <- system.time(
  gal_lgpd_out_list <- foreach(l = seq_along(p0_vals), .options.future = list(packages = "bqrgal", seed = TRUE)) %:%
    foreach(j = 1:train_kk) %dofuture% 
    {
      source("~/bqr/code/data-examples/for-github/run_gal_mcmc.R")
      run_gal_mcmc(lgpd_train_data_list[[l]][[j]], train_covar_list[[l]][[j]], p0_vals[l],
                   "laplace", priors, starting, tuning, mcmc_settings, "slice", FALSE)
    }
)

cat("done\n")
cat(paste0("Running time: ", round(gal_lgpd_runtime[3]/60), " minutes\n"))
cat("-------------------------------------------------------------------------\n")


#-------------------------------------------------------------------------------
# Fit the model to gaussian mixture data
#-------------------------------------------------------------------------------

cat("-------------------------------------------------------------------------\n")
cat("Fitting the model to gaussian mixture data\n\n")

gal_gausmix_runtime <- system.time(
  gal_gausmix_out_list <- foreach(l = seq_along(p0_vals), .options.future = list(packages = "bqrgal", seed = TRUE)) %:%
    foreach(j = 1:train_kk) %dofuture% 
    {
      source("~/bqr/code/data-examples/for-github/run_gal_mcmc.R")
      run_gal_mcmc(gausmix_train_data_list[[l]][[j]], train_covar_list[[l]][[j]], p0_vals[l],
                   "laplace", priors, starting, tuning, mcmc_settings, "slice", FALSE)
    }
)

cat("done\n")
cat(paste0("Running time: ", round(gal_gausmix_runtime[3]/60), " minutes\n"))
cat("-------------------------------------------------------------------------\n")

#-------------------------------------------------------------------------------
# Output
#-------------------------------------------------------------------------------
save(gal_normal_out_list, gal_normal_runtime,
     gal_laplace_out_list, gal_laplace_runtime,
     gal_lgpd_out_list, gal_lgpd_runtime,
     gal_gausmix_out_list, gal_gausmix_runtime,
     file = paste0(outpath, "gal_out.RData"))

cat("done fitting al models.\n")