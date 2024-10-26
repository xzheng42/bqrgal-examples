rm(list = ls())

library(readr)

inpath <- "path/to/output/folder/"
outpath <- "path/to/output/folder/"

specRound <- function(x, k = 3) trimws(format(round(x, k), nsmall = k))

matCombine <- function(mat_a, mat_b) {
  
  nn_a <- nrow(mat_a)
  nn_b <- nrow(mat_b) 
  
  mat <- matrix(nrow = nn_a + nn_b, ncol = ncol(mat_a))
  ind_a <- 2 * seq_len(nn_a) - 1
  ind_b <- 2 * seq_len(nn_b)
  
  mat[ind_a, ] <- mat_a
  mat[ind_b, ] <- mat_b
  
  mat
  
}

matCombine2 <- function(mat_a, mat_b) {
  
  mat <- matrix(nrow = 18, ncol = ncol(mat_a))
  ind_a <- c(1:3, 7:9, 13:15)
  ind_b <- c(4:6, 10:12, 16:18)
  
  mat[ind_a, ] <- mat_a
  mat[ind_b, ] <- mat_b
  
  mat
  
}

matSplit <- function(mat, k) {

  pp <- nrow(mat) / k
  mats <- array(NA, dim = c(k, ncol(mat), pp))
  for (p in 1:pp) {
    mats[, , p] <- mat[((p-1)*k+1):(k*p), ]
  }
  
  mats
  
}

load(paste0(inpath, "all_al_out_metrics.RData"))
load(paste0(inpath, "all_gal_out_metrics.RData"))

cat("-----------------------------------------------------------------------\n")
cat("Producing tables for model evaluations and comparisons\n")
cat("-----------------------------------------------------------------------\n")
load(paste0(inpath, "train_test_data.RData"))
p0_vals <- true_params$p0_vals
np0 <- length(p0_vals)
nmod <- length(all_gal_cie_median)


#-------------------------------------------------------------------------------
# Correct inclusions and exclusions
#-------------------------------------------------------------------------------
gal_cie_median <- specRound(matrix(unlist(all_gal_cie_median), nrow = np0, byrow = FALSE))
al_cie_median <- specRound(matrix(unlist(all_al_cie_median), nrow = np0, byrow = FALSE))
gal_cie_sd <- specRound(matrix(unlist(all_gal_cie_sd), nrow = np0, byrow = FALSE))
al_cie_sd <- specRound(matrix(unlist(all_al_cie_sd), nrow = np0, byrow = FALSE))

cie_median <- matCombine(gal_cie_median, al_cie_median)
cie_sd <- matCombine(gal_cie_sd, al_cie_sd)
cie <- matrix(paste0(cie_median, " (", cie_sd, ")"), ncol = nmod)
cie <- cbind(rep(c("GAL", "AL"), np0), cie)
rownames(cie) <- c("0.05", " ", "0.25", " ", "0.50", " ")
colnames(cie) <- c("Model", "Normal", "Laplace", "Normal mixture", 
                   "Log-transformed generalized Pareto")
cie_latex <- knitr::kable(cie, format = "latex", align = "c")
write_file(cie_latex, paste0(outpath, "correct_inclusion_exclusion_", true_params$nn, ".txt"))


#-------------------------------------------------------------------------------
# Regression coefficients RMSE 
#-------------------------------------------------------------------------------
pp <- nrow(all_gal_rmse_mean[[1]])

# detailed version
gal_rmse_mean <- matSplit(sapply(1:nmod, function(x) all_gal_rmse_mean[[x]]), pp)
al_rmse_mean <- matSplit(sapply(1:nmod, function(x) all_al_rmse_mean[[x]]), pp)
rmse_mean <- sapply(1:np0, function(x) matCombine(gal_rmse_mean[,,x], al_rmse_mean[,,x]), simplify = FALSE)

gal_rmse_sd <- matSplit(sapply(1:nmod, function(x) all_gal_rmse_sd[[x]]), pp)
al_rmse_sd <- matSplit(sapply(1:nmod, function(x) all_al_rmse_sd[[x]]), pp)
rmse_sd <- sapply(1:np0, function(x) matCombine(gal_rmse_sd[,,x], al_rmse_sd[,,x]), simplify = FALSE)

for (j in 1:np0) {

  rmse_mean_j <- specRound(rmse_mean[[j]], 2)
  rmse_sd_j <- specRound(rmse_sd[[j]], 2)
  rmse_j <- matrix(paste0(rmse_mean_j, " (", rmse_sd_j, ")"), ncol = ncol(rmse_mean_j))
  rmse_j <- cbind(rep(c("GAL", "AL"), pp), rmse_j)
  colnames(rmse_j) <- c("Model", "Normal", "Laplace", "Normal mixture", "Log generalized Pareto")
  rn1 <- as.matrix(paste0("$\\beta_", 1:pp, "$"))
  rn2 <- as.matrix(rep("", pp))
  rownames(rmse_j) <- as.vector(matCombine(rn1, rn2))
  est_rmse_j_latex <- knitr::kable(rmse_j, format = "latex", align = "c", escape = FALSE)
  write_file(est_rmse_j_latex, paste0(outpath, "rmse_", p0_vals[j], "_", true_params$nn, ".txt"))

}

# aggregate 
gal_rmse_mean_agg <- sapply(all_gal_rmse_mean_2, function(x) apply(x, 2, median))
al_rmse_mean_agg <- sapply(all_al_rmse_mean_2, function(x) apply(x, 2, median))
gal_rmse_sd_agg <- sapply(all_gal_rmse_mean_2, function(x) apply(x, 2, sd))
al_rmse_sd_agg <- sapply(all_al_rmse_mean_2, function(x) apply(x, 2, sd))

rmse_mean <- specRound(matCombine(gal_rmse_mean_agg, al_rmse_mean_agg), 2)
rmse_sd <- specRound(matCombine(gal_rmse_sd_agg, al_rmse_sd_agg), 2)
rmse <- matrix(paste0(rmse_mean, " (", rmse_sd, ")"), ncol = nmod)
rmse <- cbind(rep(c("GAL", "AL"), np0), rmse)
rownames(rmse) <- c("0.05", " ", "0.25", " ", "0.50", " ")
colnames(rmse) <- c("Model", "Normal", "Laplace", "Normal mixture",
                   "Log-transformed generalized Pareto")
rmse_latex <- knitr::kable(rmse, format = "latex", align = "c")
write_file(rmse_latex, paste0(outpath, "rmse_", true_params$nn, ".txt"))


#-------------------------------------------------------------------------------
# Estimation coverage
#-------------------------------------------------------------------------------
pp <- nrow(all_gal_est_cover_mean[[1]])

# coverage
gal_est_cover <- matSplit(sapply(1:nmod, function(x) all_gal_est_cover_mean[[x]]), pp)
al_est_cover <- matSplit(sapply(1:nmod, function(x) all_al_est_cover_mean[[x]]), pp)
est_cover <- sapply(1:np0, function(x) matCombine(gal_est_cover[,,x], al_est_cover[,,x]), simplify = FALSE)
for (j in 1:np0) {
  est_cover_j <- cbind(rep(c("GAL", "AL"), pp), specRound(est_cover[[j]], 2))
  colnames(est_cover_j) <- c("Model", "Normal", "Laplace", "Normal mixture", "Log generalized Pareto")
  rn1 <- as.matrix(paste0("$\\beta_", 1:pp, "$"))
  rn2 <- as.matrix(rep("", pp))
  rownames(est_cover_j) <- as.vector(matCombine(rn1, rn2))
  est_cover_j_latex <- knitr::kable(est_cover_j, format = "latex", align = "c", escape = FALSE)
  write_file(est_cover_j_latex, paste0(outpath, "beta_est_cover_", p0_vals[j], "_", true_params$nn, ".txt"))
  
}

#-------------------------------------------------------------------------------
# Prediction interval score
#-------------------------------------------------------------------------------
# interval score
gal_pred_is_mean <- sapply(all_gal_pred_is_mean, function(x) apply(x, 2, median))
al_pred_is_mean <- sapply(all_al_pred_is_mean, function(x) apply(x, 2, median))
gal_pred_is_sd <- sapply(all_gal_pred_is_mean, function(x) apply(x, 2, sd))
al_pred_is_sd <- sapply(all_al_pred_is_mean, function(x) apply(x, 2, sd))

pred_is_mean <- specRound(matCombine(gal_pred_is_mean, al_pred_is_mean), 2)
pred_is_sd <- specRound(matCombine(gal_pred_is_sd, al_pred_is_sd), 2)
pred_is <- matrix(paste0(pred_is_mean, " (", pred_is_sd, ")"), ncol = nmod)
pred_is <- cbind(rep(c("GAL", "AL"), np0), pred_is)
rownames(pred_is) <- c("0.05", " ", "0.25", " ", "0.50", " ")
colnames(pred_is) <- c("Model", "Normal", "Laplace", "Normal mixture", 
                       "Log-transformed generalized Pareto")
pred_is_latex <- knitr::kable(pred_is, format = "latex", align = "c", escape = FALSE)
write_file(pred_is_latex, paste0(outpath, "pred_interval_score_", true_params$nn, ".txt"))


