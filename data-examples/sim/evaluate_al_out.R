rm(list = ls())

set.seed(42)

inpath <- "path/to/output/folder/"
outpath <- "path/to/output/folder/"

library(bqrgal)
library(scoringutils)

#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
load(paste0(inpath, "train_test_data.RData"))

nn <- true_params$nn
NN <- true_params$NN
pp <- length(true_params$bb)
p0_vals <- true_params$p0_vals
train_kk <- true_params$train_kk
test_kk <- true_params$test_kk
true_ind <- true_params$true_ind
bb <- true_params$bb

all_train_data <- list(normal = normal_train_data_list,
                       laplace = laplace_train_data_list,
                       gausmix = gausmix_train_data_list,
                       lgpd = lgpd_train_data_list)

all_test_data <- list(normal = normal_test_data_list,
                      laplace = laplace_test_data_list,
                      gausmix = gausmix_test_data_list,
                      lgpd = lgpd_test_data_list)

load(paste0(inpath, "al_out.RData"))
all_al_out <- list(normal = al_normal_out_list,
                   laplace = al_laplace_out_list,
                   gausmix = al_gausmix_out_list,
                   lgpd = al_lgpd_out_list)

#-------------------------------------------------------------------------------
# Compute correct inclusion and exclusions (CIE)
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Compute CIE\n")
cat("-------------------------------------------------------------------------\n")

all_al_cie_median <- vector("list", length = length(all_al_out))
all_al_cie_sd <- vector("list", length = length(all_al_out))

al_cie <- vector("list", length = length(p0_vals))

for (m in seq_along(all_al_out)) {
  
  al_out_list <- all_al_out[[m]]
  train_data_list <- all_train_data[[m]]
  
  for (l in seq_along(p0_vals)) {
    
    cie_l <- rep(NA, train_kk)
    
    for (j in 1:train_kk) {
      
      yy <- train_data_list[[l]][[j]]
      XX <- train_covar_list[[l]][[j]]
      sd_xx <- apply(XX, 2, sd)
      sd_yy <- sd(yy)
      be_sams <- al_out_list[[l]][[j]]$post_sams$be
      be_star_sams <- apply(be_sams, 2, function(x) x * sd_xx / sd_yy)
      be_ind <- abs(be_star_sams) > 0.1
      cie_l[j] <- mean(apply(be_ind, 2, function(x) sum(x == true_ind)) / pp)
      
    }
    
    al_cie[[l]] <- cie_l
    
  }
  
  all_al_cie_median[[m]] <- sapply(al_cie, function(x) median(x))
  all_al_cie_sd[[m]] <- sapply(al_cie, function(x) sd(x))
  
}


#-------------------------------------------------------------------------------
# Compute root mean squared error (RMSE)
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Compute RMSE\n")
cat("-------------------------------------------------------------------------\n")

all_al_rmse_mean <- vector("list", length = length(all_al_out))
all_al_rmse_sd <- vector("list", length = length(all_al_out))
all_al_rmse_mean_2 <- vector("list", length = length(all_al_out))

al_rmse <- vector("list", length = length(p0_vals))

for (m in seq_along(all_al_out)) {
  
  al_out_list <- all_al_out[[m]]
  train_data_list <- all_train_data[[m]]
  
  for (l in seq_along(p0_vals)) {
    
    rmse_mat <- array(NA, dim = c(train_kk, pp))
    
    for (j in 1:train_kk) {
      
      yy <- train_data_list[[l]][[j]]
      XX <- train_covar_list[[l]][[j]]
      be_sams <- al_out_list[[l]][[j]]$post_sams$be
      be_err <- apply(be_sams, 2, function(x) x - bb)
      rmse_mat[j, ] <- apply(be_err, 1, function(x) sqrt(mean(x^2)))
      
    }
    
    al_rmse[[l]] <- rmse_mat
    
  }
  
  all_al_rmse_mean[[m]] <- sapply(al_rmse, function(x) colMeans(x))
  all_al_rmse_sd[[m]] <- sapply(al_rmse, function(x) apply(x, 2, sd))
  all_al_rmse_mean_2[[m]] <- sapply(al_rmse, function(x) rowMeans(x))
  
}


#-------------------------------------------------------------------------------
# Coverage probability (CP) for beta and prediction interval score (IS)
#-------------------------------------------------------------------------------
cat("-------------------------------------------------------------------------\n")
cat("Compute CP and IS\n")
cat("-------------------------------------------------------------------------\n")

all_al_est_cover_mean <- vector("list", length = length(all_al_out))
all_al_pred_is_mean <- vector("list", length = length(all_al_out))

al_est_cover <- vector("list", length = length(p0_vals))
al_pred_is <- vector("list", length = length(p0_vals))

for (m in seq_along(all_al_out)) {
  
  al_out_list <- all_al_out[[m]]
  test_data_list <- all_test_data[[m]]
  
  for (l in seq_along(p0_vals)) {
    
    est_cover_mat <- array(NA, dim = c(train_kk, length(bb)))
    pred_is_mat <- array(NA, dim = c(test_kk, NN))
    
    for (j in 1:train_kk) {
      
      bqral_out <- al_out_list[[l]][[j]]
      
      # estimation coverage
      be_sams <- bqral_out$post_sams$be
      be_qq <- apply(be_sams, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
      est_cover_bool <- sapply(1:(ncol(be_qq)), function(x) {
        bb[x] >= be_qq[1,x] & bb[x] <= be_qq[2,x]
      })
      est_cover_mat[j, ] <- est_cover_bool
      
      # prediction interval score
      yy_test <- test_data_list[[l]][[j]]
      XX_test <- test_covar_list[[l]][[j]]
      pred_sams <- predict(bqral_out, XX_test, probs = c(0.025, 0.975), FALSE)
      pred_is_mat[j, ] <- interval_score(yy_test, pred_sams$yy_qq[1, ], pred_sams$yy_qq[2, ], 95, weigh = FALSE)
      
    }
    
    al_est_cover[[l]] <- est_cover_mat
    al_pred_is[[l]] <- pred_is_mat
    
  }
  
  all_al_est_cover_mean[[m]] <- sapply(al_est_cover, function(x) colMeans(x))
  all_al_pred_is_mean[[m]] <- sapply(al_pred_is, function(x) rowMeans(x))
  
}


#-------------------------------------------------------------------------------
# Output
#-------------------------------------------------------------------------------
save(all_al_cie_median, all_al_cie_sd,
     all_al_rmse_mean, all_al_rmse_sd, all_al_rmse_mean_2, 
     all_al_est_cover_mean, all_al_pred_is_mean, 
     file = paste0(outpath, "all_al_out_metrics.RData"))

cat("done the evaluation.\n\n")
