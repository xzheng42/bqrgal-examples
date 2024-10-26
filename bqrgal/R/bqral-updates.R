updateBetaAL <- function(be_mu0, be_Sigma0, yy, XX, sigma, vv, AA, BB) {

  be_Sigma0_inv <- chol2inv(chol(be_Sigma0))
  XX_vv <- apply(XX, 2, function(x) x / sqrt(BB * sigma * vv))
  be_Sigma1 <- chol2inv(chol(be_Sigma0_inv + t(XX_vv) %*% XX_vv))
  weights <- (yy - AA * vv) / (BB * sigma * vv)
  resids <- apply(XX, 2, function(x) x * weights)
  be_mu1 <- be_Sigma1 %*% (be_Sigma0_inv %*% be_mu0 + colSums(resids))

  t(mvtnorm::rmvnorm(1, be_mu1, be_Sigma1))

}


updateBeta0AL <- function(be0_sigmasq0, yy, XX, be, sigma, vv, AA, BB) {

  be0_sigmasq1 <- 1 / (1 / be0_sigmasq0 + sum(1 / vv) / (BB * sigma))
  resid_sum <- sum((yy - XX %*% be - AA * vv) / vv)
  be0_mu1 <- be0_sigmasq1 * resid_sum / (BB * sigma)

  rnorm(1, be0_mu1, be0_sigmasq1)

}


updateBetaLassoAL <- function(Omega, yy, XX, be0, sigma, vv, AA, BB) {

  XX_vv <- apply(XX, 2, function(x) x / sqrt(BB * sigma * vv))
  be_Sigma1 <- chol2inv(chol(Omega + t(XX_vv) %*% XX_vv))
  weights <- (yy - be0 - AA * vv) / (BB * sigma * vv)
  resids <- apply(XX, 2, function(x) x * weights)
  be_mu1 <- be_Sigma1 %*% (colSums(resids))

  t(mvtnorm::rmvnorm(1, be_mu1, be_Sigma1))

}


updateSigmaAL <- function(u_sigma, v_sigma, yy, mu, vv, AA, BB, nn) {

  nu <- u_sigma + 1.5 * nn
  cc <- v_sigma + sum(vv) + 0.5 * sum((yy - mu - AA * vv)^2 / vv) / BB
  1 / rgamma(1, nu, cc)

}


updateVVAL <- function(yy, mu, sigma, AA, BB, nn) {

  a_i <- (yy - mu)^2 / (BB * sigma)
  b_i <- 2 / sigma + AA^2 / (BB * sigma)
  vv <- rep(NA, nn)
  for (i in 1:nn) {
    vv[i] <- GIGrvg::rgig(1, 0.5, a_i[i], b_i)
  }

  vv

}


updateOmegaAL <- function(kk, be, etasq) {

  omega <- rep(NA, kk)
  besq <- be^2

  for (j in 1:kk) {
    omega[j] <-  GIGrvg::rgig(1, 0.5, besq[j], etasq)
  }

  omega

}
