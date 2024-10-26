################################################################################
# Section 5.2: Boston housing data
################################################################################

rm(list = ls())

set.seed(42)

outpath <- "path/to/save/results/"

library(scales)
library(readr)
library(spData)
library(MCMCglmm)

library(bqrgal)

#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
data("boston")

yy <- log(boston.c$CMEDV)
XX <- cbind(boston.c$LON, boston.c$LAT, boston.c$CRIM, boston.c$ZN, boston.c$INDUS,
            boston.c$CHAS, boston.c$NOX, boston.c$RM, boston.c$AGE, boston.c$DIS,
            boston.c$RAD, boston.c$TAX, boston.c$PTRATIO, boston.c$B, boston.c$LSTAT)
XX <- apply(XX, 2, scale)

nn <- length(yy)

#-------------------------------------------------------------------------------
# Fit the GAL and AL models
#-------------------------------------------------------------------------------
p0_vals <- c(0.1, 0.9)
pp <- ncol(XX)

priors <- list("intercept_gaus_var" = 100,
               "sigma_invgamma" = c(2, 2),
               "eta_gamma" = c(0.1, 0.1))
starting <- list(be0 = 0, sigma = 1, vv = 1 * rexp(nn),
                 ss = TruncatedDistributions::rtnorm(nn, 0, 1, 0, Inf),
                 omega = rep(1, pp))
tuning <- list("step_size" = 0.01)
mcmc_settings <- list(n_iter = 150000, n_burn = 50000, n_thin = 20, n_report = 5000)


gal_out_list <- vector("list", length = length(p0_vals))
al_out_list <- vector("list", length = length(p0_vals))

for (l in seq_along(p0_vals)) {

  p0 <- p0_vals[l]

  u_ga <- find_ga_lb(p0, c(-100, 0))
  v_ga <- find_ga_ub(p0, c(0, 100))

  priors[["ga_uniform"]] <- c(u_ga, v_ga)
  starting$ga <- (u_ga + v_ga) / 2

  bqrgal_out <- bgal(yy, XX, p0, "laplace", priors, starting, tuning,
                     mcmc_settings, "slice", TRUE)
  bqral_out <- bal(yy, XX, p0, "laplace", priors, starting,
                   mcmc_settings, TRUE)
  cat(paste0("done fitting models for the ", l, "-th p0 value.\n"))

  gal_out_list[[l]] <- bqrgal_out
  al_out_list[[l]] <- bqral_out

}

#-------------------------------------------------------------------------------
# Log-likelihood and BIC
#-------------------------------------------------------------------------------
tbl <- array(NA, dim = c(length(p0_vals) * 2, 5))
colnames(tbl) <- c("Quantile", "Model", "Mean (95\\% CrI) for $\\gamma$", "Log-likelihood", "BIC")

for (l in seq_along(p0_vals)) {

  p0 <- p0_vals[l]
  bqrgal_out <- gal_out_list[[l]]
  bqral_out <- al_out_list[[l]]

  gal_sams <- bqrgal_out$post_sams
  al_sams <- bqral_out$post_sams

  gal_mu <- mean(gal_sams$be0) + XX %*% rowMeans(gal_sams$be)
  al_mu <- mean(al_sams$be0) + XX %*% rowMeans(al_sams$be)
  gal_loglik <- sum(sapply(1:nn, function(x) dgal(yy[x], p0, gal_mu[x], mean(gal_sams$sigma), mean(gal_sams$ga), TRUE)))
  al_loglik <- sum(sapply(1:nn, function(x) dgal(yy[x], p0, al_mu[x], mean(al_sams$sigma), 0, TRUE)))

  tbl[1 + 2*(l-1), 4] <- al_loglik
  tbl[2*l, 4] <- gal_loglik

  gal_bic <- -2 * gal_loglik + log(nn) * (ncol(XX) + 3)
  al_bic <- -2 * al_loglik + log(nn) * (ncol(XX) + 2)

  tbl[1 + 2*(l-1), 5] <- al_bic
  tbl[2*l, 5] <- gal_bic

}

tbl <- round(tbl)

for (l in seq_along(p0_vals)) {

  bqrgal_out <- gal_out_list[[l]]
  gal_sams <- bqrgal_out$post_sams

  tbl[1 + 2*(l-1), 3] <- ""
  tbl[2*l, 3] <- summarize(gal_sams$ga)

}

tbl[, 1] <- ""
tbl[c(1,3), 1] <- paste0("$p_0 = ", p0_vals, "$")
tbl[, 2] <- rep(c("$M_0$", "$M_1$"), length(p0_vals))
tbl_latex <- knitr::kable(tbl, format = "latex", align = "c", escape = FALSE)
readr::write_file(tbl_latex, paste0(outpath, "loglik_bic.txt"))


#-------------------------------------------------------------------------------
# Coefficients
#-------------------------------------------------------------------------------
covar_names <- c("LON","LAT","CRIM","ZN","INDUS","CHAS","NOX",
                 "RM","AGE","DIS","RAD","TAX","PTRATIO","B",
                 "LSTAT")
id <- 1:ncol(XX)

col_al <- rgb(80, 160, 200,maxColorValue = 255)
col_gal <- rgb(4, 90, 141,maxColorValue = 255)

for (l in seq_along(p0_vals)) {

  p0 <- p0_vals[l]
  bqrgal_out <- gal_out_list[[l]]
  bqral_out <- al_out_list[[l]]

  gal_sams <- bqrgal_out$post_sams
  al_sams <- bqral_out$post_sams

  gal_be_mode <- MCMCglmm::posterior.mode(mcmc(t(gal_sams$be)))
  al_be_mode <- MCMCglmm::posterior.mode(mcmc(t(al_sams$be)))
  gal_be_qq <- apply(gal_sams$be, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  al_be_qq <- apply(al_sams$be, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

  png(paste0(outpath, "boston_coef_p0_", p0, ".png"), width=1200, height=800, pointsize=20)
  par(mar = c(4.1, 4.1, 2.1, 2.1))
  plot(id - .15, gal_be_mode,
       xli = c(0.5, max(id) + .5), ylim = c(-.2, .2 + .025 * (p0 ==.9)),
       lwd = 2, pch = 24, xlab = "ID", ylab = "Effect", col = col_gal)
  arrows(id - .15, gal_be_qq[1,], id - .15, gal_be_qq[2,], length = 0.05, angle = 90,
         code = 3, col = col_gal, lwd = 2)
  points(id + .15, al_be_mode, pch = 19, col = col_al, lwd = 2)
  arrows(id + .15, al_be_qq[1,], id + .15, al_be_qq[2,], length = 0.05, angle = 90,
         code = 3, col = col_al, lwd = 2)
  abline(h = 0, lty = 2)
  text(id, apply(cbind(gal_be_qq[2,], al_be_qq[2,]), 1, max) + .012,
       covar_names, cex = .8, col = "black")
  legend("topleft", c("GAL","AL"), lty = c(1,1), pch = c(24,19), lwd = 2, cex = 1.1,
         col=c(col_gal, col_al), bty = "n")
  dev.off()
  
}

