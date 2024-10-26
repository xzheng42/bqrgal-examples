################################################################################
# Section 5.1: Immunoglobulin-G data
################################################################################

rm(list = ls())

set.seed(42)

outpath <- "path/to/save/results/"

library(scales)
library(Brq)
library(readr)

library(bqrgal)

#-------------------------------------------------------------------------------
# Input data
#-------------------------------------------------------------------------------
data("ImmunogG")

yy <- ImmunogG$IgG
xx <- ImmunogG$Age
XX <- cbind(1, xx, xx^2)

nn <- length(yy)

#-------------------------------------------------------------------------------
# Fit the GAL and AL models
#-------------------------------------------------------------------------------
p0_vals <- c(0.05, 0.25, 0.5, 0.75, 0.85, 0.95)
pp <- ncol(XX)

priors <- list("beta_gaus" = list(mean_vec = rep(0, pp), var_mat = 100 * diag(pp)),
                "sigma_invgamma" = c(2, 2))
starting <- list(sigma = 1, vv = 1 * rexp(nn),
                 ss = TruncatedDistributions::rtnorm(nn, 0, 1, 0, Inf))
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

  bqrgal_out <- bgal(yy, XX, p0, "gaussian", priors, starting, tuning,
                     mcmc_settings, "slice", TRUE)
  bqral_out <- bal(yy, XX, p0, "gaussian", priors, starting,
                   mcmc_settings, TRUE)
  cat(paste0("done fitting models for the ", l, "-th p0 value.\n"))

  gal_out_list[[l]] <- bqrgal_out
  al_out_list[[l]] <- bqral_out

}


#-------------------------------------------------------------------------------
# Posterior error densities and regression functionals
#-------------------------------------------------------------------------------
eps_mat <- rbind(seq(-15, 20, length = 100),
                 seq(-10, 15, length = 100),
                 seq(-10, 10, length = 100),
                 seq(-15, 10, length = 100),
                 seq(-20, 10, length = 100),
                 seq(-20, 10, length = 100))

for (l in seq_along(p0_vals)) {

  p0 <- p0_vals[l]
  bqrgal_out <- gal_out_list[[l]]
  bqral_out <- al_out_list[[l]]

  eps <- eps_mat[l, ]
  gal_sams <- bqrgal_out$post_sams
  al_sams <- bqral_out$post_sams
  nsams <- length(gal_sams$ga)

  gal_err_dens <- sapply(1:nsams, function(x) dgal(eps, p0, 0, gal_sams$sigma[x], gal_sams$ga[x]))
  al_err_dens <- sapply(1:nsams, function(x) dgal(eps, p0, 0, al_sams$sigma[x], 0))
  gal_err_dens_mean <- rowMeans(gal_err_dens)
  al_err_dens_mean <- rowMeans(al_err_dens)
  gal_err_dens_qq <- apply(gal_err_dens, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  al_err_dens_qq <- apply(al_err_dens, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

  png(paste0(outpath, "igg_err_dens_p0_", p0, ".png"), width=600, height=500, pointsize=20)
  par(mar = c(4.1, 4.1, 2.1, 2.1), lwd = 3)
  plot(NULL, xlim = range(eps),
       ylim = c(0, max(max(gal_err_dens_qq), max(al_err_dens_qq)) + 0.01),
       xlab = expression(epsilon), ylab = "Density",
       main = substitute(paste(p[0], " = ", k), list(k = p0)))
  polygon(c(eps, rev(eps)), c(gal_err_dens_qq[1,], rev(gal_err_dens_qq[2,])),
          border = NA, col = alpha('blue', 0.4), lwd = 2)
  lines(eps, gal_err_dens_mean, col = "blue", lty = 2)
  polygon(c(eps, rev(eps)), c(al_err_dens_qq[1,], rev(al_err_dens_qq[2,])),
          border = NA, col = alpha('red', 0.4), lwd = 2)
  lines(eps, al_err_dens_mean, col = "red", lty = 4)
  dev.off()

  xx_grid <- seq(min(xx), max(xx), length = 1000)
  XX_grid <- cbind(1, xx_grid, xx_grid^2)
  gal_reg <- apply(gal_sams$be, 2, function(bb) XX_grid %*% bb)
  al_reg <- apply(al_sams$be, 2, function(bb) XX_grid %*% bb)
  gal_reg_mean <- rowMeans(gal_reg)
  al_reg_mean <- rowMeans(al_reg)
  gal_reg_qq <- apply(gal_reg, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  al_reg_qq <- apply(al_reg, 1, function(x) quantile(x, probs = c(0.025, 0.975)))

  png(paste0(outpath, "igg_reg_p0_", p0, ".png"), width=600, height=500, pointsize=20)
  par(mar = c(4.1, 4.1, 2.1, 2.1), lwd = 3)
  plot(xx, yy, xlab = "Age (Years)", ylab = "IgG (g/l)", pch = 19, col = 'black', cex = 0.2,
       main = substitute(paste(p[0], " = ", k), list(k = p0)))
  # plot(NULL, xlim = range(xx), ylim = range(yy), xlab = "Age (Years)", ylab = "IgG (g/l)")
  polygon(c(xx_grid, rev(xx_grid)), c(gal_reg_qq[1,], rev(gal_reg_qq[2,])),
          border = NA, col = alpha('blue', 0.4), lwd = 2)
  lines(xx_grid, gal_reg_mean, col = "blue", lty = 2)
  polygon(c(xx_grid, rev(xx_grid)), c(al_reg_qq[1,], rev(al_reg_qq[2,])),
          border = NA, col = alpha('red', 0.4), lwd = 2)
  lines(xx_grid, al_reg_mean, col = "red", lty = 4)
  dev.off()


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

  gal_mu <- XX %*% rowMeans(gal_sams$be)
  al_mu <- XX %*% rowMeans(al_sams$be)
  gal_loglik <- sum(sapply(1:nn, function(x) dgal(yy[x], p0, gal_mu[x], mean(gal_sams$sigma), mean(gal_sams$ga), TRUE)))
  al_loglik <- sum(sapply(1:nn, function(x) dgal(yy[x], p0, al_mu[x], mean(al_sams$sigma), 0, TRUE)))

  tbl[1 + 2*(l-1), 4] <- al_loglik
  tbl[2*l, 4] <- gal_loglik

  gal_bic <- -2 * gal_loglik + log(nn) * (ncol(XX) + 2)
  al_bic <- -2 * al_loglik + log(nn) * (ncol(XX) + 1)

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
tbl[c(1,3,5,7,9,11), 1] <- paste0("$p_0 = ", p0_vals, "$")
tbl[, 2] <- rep(c("$M_0$", "$M_1$"), length(p0_vals))
tbl_latex <- knitr::kable(tbl, format = "latex", align = "c", escape = FALSE)
readr::write_file(tbl_latex, paste0(outpath, "loglik_bic.txt"))
