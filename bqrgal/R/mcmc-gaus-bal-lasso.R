mcmc_gaus_bal_lasso <- function(yy, XX, p0,
                                 be0_sigmasq0,
                                 u_eta, v_eta,
                                 u_sigma, v_sigma,
                                 be0, sigma, vv, omega,
                                 mcmc_settings,
                                 verbose) {

  #--------------------------
  # MCMC settings
  #--------------------------
  if (verbose) {

    if (is.null(mcmc_settings$n_report)) {
      nreport <- 1000
    } else {
      nreport <- mcmc_settings$n_report
    }

    cat("--------------------------------------------\n")
    cat("\t  Running MCMC\n");
    cat("--------------------------------------------\n")

  }
  niter <- mcmc_settings$n_iter
  nburn <- mcmc_settings$n_burn
  nthin <- mcmc_settings$n_thin

  #--------------------------
  # Initialization
  #--------------------------
  nn <- length(yy)
  kk <- ncol(XX)

  pp <- p0
  AA <- A_func(pp)
  BB <- B_func(pp)
  Omega <- diag(1 / omega)

  # empty stuff to save samples
  be0_save <- array(NA, dim = niter - nburn)
  be_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  vv_save <- array(NA, dim = c(nn, niter - nburn))
  sigma_save <- array(NA, dim = niter - nburn)
  etasq_save <- array(NA, dim = niter - nburn)
  omega_save <- array(NA, dim = c(length(omega), niter - nburn))

  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0

  for (iter in 1:niter) {

    start_time <- Sys.time()

    #--------------------------
    # update beta
    #--------------------------
    Omega <- diag(1 / as.vector(omega))
    be <- updateBetaLassoAL(Omega, yy, XX, be0, sigma, vv, AA, BB)
    mu <- be0 + XX %*% be

    #--------------------------
    # update intercept
    #--------------------------
    be0 <- updateBeta0AL(be0_sigmasq0, yy, XX, be, sigma, vv, AA, BB)

    #--------------------------
    # update eta
    #--------------------------
    etasq <- rgamma(1, u_eta + kk, v_eta + sum(omega) / 2)

    #--------------------------
    # update omega
    #--------------------------
    omega <- updateOmegaAL(kk, be, etasq)

    #--------------------------
    # update vv
    #--------------------------
    vv <- updateVVAL(yy, mu, sigma, AA, BB, nn)

    #--------------------------
    # update sigma
    #--------------------------
    sigma <- updateSigmaAL(u_sigma, v_sigma, yy, mu, vv, AA, BB, nn)

    #-------------------------------
    # calculate block running time
    #-------------------------------
    end_time <- Sys.time()
    block_runtime <- block_runtime + as.numeric(end_time - start_time)

    #--------------------------
    # Print MCMC progress
    #--------------------------
    if (verbose) {

      if (iter %% nreport == 0) {

        cat(paste0("Iterations: ", iter, "/", niter,
                   "  Percentage: ", specRound(iter / niter * 100, 2), "%\n"))

        ert <- (block_runtime / nreport) * (niter - iter)
        cat(paste0("Estimated remaining time: ",
                   specRound(ert / 60, 2), " minutes \n"))

        block_runtime <- 0

      }

    }

    #--------------------------
    # Save samples
    #--------------------------
    if (iter > nburn) {
      be0_save[iter - nburn] <- be0
      be_save[, iter - nburn] <- be
      vv_save[, iter - nburn] <- vv
      sigma_save[iter - nburn] <- sigma
      etasq_save[iter - nburn] <- etasq
      omega_save[, iter - nburn] <- omega
    }

  }

  #--------------------------
  # Thinning
  #--------------------------
  selc_idx <- seq(1, niter - nburn, by = nthin)

  post_sams <- list(be0 = be0_save[selc_idx],
                    be = be_save[, selc_idx],
                    vv = vv_save[, selc_idx],
                    sigma = sigma_save[selc_idx],
                    etasq = etasq_save[selc_idx],
                    omega = omega_save[, selc_idx])

  list(post_sams = post_sams)


}
