mcmc_gaus_bal <- function(yy, XX, p0,
                           be_mu0, be_Sigma0,
                           u_sigma, v_sigma,
                           sigma, vv,
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
  pp <- p0
  AA <- A_func(pp)
  BB <- B_func(pp)

  # empty stuff to save samples
  be_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  vv_save <- array(NA, dim = c(nn, niter - nburn))
  sigma_save <- array(NA, dim = niter - nburn)

  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0

  for (iter in 1:niter) {

    start_time <- Sys.time()

    #--------------------------
    # update beta
    #--------------------------
    be <- updateBetaAL(be_mu0, be_Sigma0, yy, XX, sigma, vv, AA, BB)
    mu <- XX %*% be

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
      be_save[, iter - nburn] <- be
      vv_save[, iter - nburn] <- vv
      sigma_save[iter - nburn] <- sigma
    }

  }

  #--------------------------
  # Thinning
  #--------------------------
  selc_idx <- seq(1, niter - nburn, by = nthin)

  post_sams <- list(be = be_save[, selc_idx],
                    vv = vv_save[, selc_idx],
                    sigma = sigma_save[selc_idx])

  list(post_sams = post_sams)


}


