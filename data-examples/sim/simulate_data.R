################################################################################
### Generating synthetic data
################################################################################
rm(list = ls())

set.seed(42)

outpath <- "path/to/output/folder/"

library(bqrgal)

#-------------------------------------------------------------------------------
# Set up parameters
#-------------------------------------------------------------------------------
nn <- 100
cov_mat <- 0.5^as.matrix(dist(1:8))
bb <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
true_ind <- c(1, 1, 0, 0, 1, 0, 0, 0)
p0_vals <- c(0.05, 0.25, 0.5)

# training and testing set replication sizes
train_kk <- 100
test_kk <- 100

# test set size
NN <- 100

#-------------------------------------------------------------------------------
# Source functions to simulate data
#-------------------------------------------------------------------------------
source("path/to/data-examples/simu_funcs.R")

#-------------------------------------------------------------------------------
# Covariate
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Simulating covariate \n")
cat("-------------------------------------------------------------------------\n")
train_covar_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  XX_train_list <- vector("list", length = train_kk)
  
  for (j in 1:train_kk) {
    XX <- mvtnorm::rmvnorm(nn, sigma = cov_mat)
    XX_train_list[[j]] <- XX
  }
  
  train_covar_list[[l]] <- XX_train_list
  
}

test_covar_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  XX_test_list <- vector("list", length = test_kk)
  
  for (j in 1:test_kk) {
    XX_test <- mvtnorm::rmvnorm(NN, sigma = cov_mat)
    XX_test_list[[j]] <- XX_test
  }
  
  test_covar_list[[l]] <- XX_test_list
  
}

#-------------------------------------------------------------------------------
# Normal
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Simulating data - normal\n")
cat("-------------------------------------------------------------------------\n")
normal_train_data_list <- vector("list", length = length(p0_vals))
sigma <- 3
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_train_list <- vector("list", length = train_kk)
  XX_train_list <- train_covar_list[[l]]
  
  for (j in 1:train_kk) {
    XX  <- XX_train_list[[j]]
    mu <- XX %*% bb
    yy_train_list[[j]] <- simGausQr(nn, p0, mu, sigma)
  }
  
  normal_train_data_list[[l]] <- yy_train_list
  
}

normal_test_data_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_test_list <- vector("list", length = test_kk)
  XX_test_list <- test_covar_list[[l]]
  
  for (j in 1:test_kk) {
    XX_test <- XX_test_list[[j]]
    mu_test <- XX_test %*% bb
    yy_test_list[[j]] <- simGausQr(NN, p0, mu_test, sigma)
  }
  
  normal_test_data_list[[l]] <- yy_test_list
  
}

#-------------------------------------------------------------------------------
# Laplace
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Simulating data - laplace\n")
cat("-------------------------------------------------------------------------\n")
laplace_train_data_list <- vector("list", length = length(p0_vals))
sigma <- 3
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_train_list <- vector("list", length = train_kk)
  XX_train_list <- train_covar_list[[l]]
  
  for (j in 1:train_kk) {
    XX  <- XX_train_list[[j]]
    mu <- XX %*% bb
    yy_train_list[[j]] <- simLaplaceQr(nn, p0, mu, sigma)
  }
  
  laplace_train_data_list[[l]] <- yy_train_list
  
}

laplace_test_data_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_test_list <- vector("list", length = test_kk)
  XX_test_list <- test_covar_list[[l]]
  
  for (j in 1:test_kk) {
    XX_test <- XX_test_list[[j]]
    mu_test <- XX_test %*% bb
    yy_test_list[[j]] <- simLaplaceQr(NN, p0, mu_test, sigma)
  }
  
  laplace_test_data_list[[l]] <- yy_test_list
  
}


#-------------------------------------------------------------------------------
# lgpd
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Simulating data - lgpd\n")
cat("-------------------------------------------------------------------------\n")
lgpd_train_data_list <- vector("list", length = length(p0_vals))
xi <- 3
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_train_list <- vector("list", length = train_kk)
  XX_train_list <- train_covar_list[[l]]
  
  for (j in 1:train_kk) {
    XX <- XX_train_list[[j]]
    mu <- XX %*% bb
    yy_train_list[[j]] <- simLogGPD(nn, p0, mu, xi)
  }
  
  lgpd_train_data_list[[l]] <- yy_train_list
  
}

lgpd_test_data_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_test_list <- vector("list", length = test_kk)
  XX_test_list <- test_covar_list[[l]]
  
  for (j in 1:test_kk) {
    XX_test <- XX_test_list[[j]]
    mu_test <- XX_test %*% bb
    yy_test_list[[j]] <- simLogGPD(NN, p0, mu_test, xi)
  }
  
  lgpd_test_data_list[[l]] <- yy_test_list
  
}

#-------------------------------------------------------------------------------
# Gaussian mixture
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Simulating data - gaussian mixture\n")
cat("-------------------------------------------------------------------------\n")
gausmix_train_data_list <- vector("list", length = length(p0_vals))
sigma <- c(sqrt(1), sqrt(5))
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_train_list <- vector("list", length = train_kk)
  XX_train_list <- train_covar_list[[l]]
  
  for (j in 1:train_kk) {
    XX <- XX_train_list[[j]] 
    mu <- XX %*% bb
    yy_train_list[[j]] <- simGausMixQr(nn, p0, mu, sigma)
  }
  
  gausmix_train_data_list[[l]] <- yy_train_list
  
}

gausmix_test_data_list <- vector("list", length = length(p0_vals))
for (l in seq_along(p0_vals)) {
  
  p0 <- p0_vals[l]
  yy_test_list <- vector("list", length = test_kk)
  XX_test_list <- test_covar_list[[l]]
  
  for (j in 1:test_kk) {
    XX_test <- XX_test_list[[j]]
    mu_test <- XX_test %*% bb
    yy_test_list[[j]] <- simGausMixQr(NN, p0, mu_test, sigma)
  }
  
  gausmix_test_data_list[[l]] <- yy_test_list
  
}


#-------------------------------------------------------------------------------
# Output
#-------------------------------------------------------------------------------
true_params <- list(nn = nn, NN = NN, train_kk = train_kk, test_kk = test_kk, 
                    cov_mat = cov_mat, bb = bb, true_ind = true_ind, p0_vals = p0_vals)

save(true_params, train_covar_list, test_covar_list,
     normal_train_data_list, normal_test_data_list,
     laplace_train_data_list, laplace_test_data_list,
     lgpd_train_data_list, lgpd_test_data_list,
     gausmix_train_data_list, gausmix_test_data_list,
     file = paste0(outpath, "train_test_data.RData"))

cat("done the simulation.\n")
