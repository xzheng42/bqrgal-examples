mcmc_gaus_bgal <- function(yy, XX, p0,
                           be_mu0, be_Sigma0,
                           u_sigma, v_sigma,
                           u_ga, v_ga, se_ga,
                           ga, sigma, vv, ss,
                           mcmc_settings, step_size,
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
  pp <- p_func(ga, p0)
  AA <- A_func(pp)
  BB <- B_func(pp)
  CC <- C_func(ga, pp)

  # empty stuff to save samples
  be_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  vv_save <- array(NA, dim = c(nn, niter - nburn))
  ss_save <- array(NA, dim = c(nn, niter - nburn))
  sigma_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = niter - nburn)

  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0

  for (iter in 1:niter) {

    start_time <- Sys.time()

    #--------------------------
    # update beta
    #--------------------------
    be <- updateBeta(be_mu0, be_Sigma0, yy, XX, sigma, ga, ss, vv, AA, BB, CC)
    mu <- XX %*% be

    #--------------------------
    # update vv
    #--------------------------
    vv <- updateVV(yy, mu, sigma, ga, ss, AA, BB, CC, nn)

    #--------------------------
    # update ss
    #--------------------------
    ss <- updateSS(yy, mu, sigma, ga, vv, AA, BB, CC, nn)

    #--------------------------
    # update sigma
    #--------------------------
    sigma <- updateSigma(u_sigma, v_sigma, yy, mu, ga, ss, vv, AA, BB, CC, nn)

    #--------------------------
    # update ga
    #--------------------------
    ga <- uniSlice(x0 = ga, g = logcondpga, w = step_size, m = Inf,
                   lower = u_ga, upper = v_ga, yy, p0, mu, sigma, ss, vv)

    pp <- p_func(ga, p0)
    AA <- A_func(pp)
    BB <- B_func(pp)
    CC <- C_func(ga, pp)

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
      ss_save[, iter - nburn] <- ss
      sigma_save[iter - nburn] <- sigma
      ga_save[iter - nburn] <- ga
    }

  }

  #--------------------------
  # Thinning
  #--------------------------
  selc_idx <- seq(1, niter - nburn, by = nthin)

  post_sams <- list(be = be_save[, selc_idx],
                    vv = vv_save[, selc_idx],
                    ss = ss_save[, selc_idx],
                    sigma = sigma_save[selc_idx],
                    ga = ga_save[selc_idx])

  list(post_sams = post_sams)


}


mcmc_gaus_bgal_mh <- function(yy, XX, p0,
                              be_mu0, be_Sigma0,
                              u_sigma, v_sigma,
                              u_ga, v_ga, se_ga,
                              ga, sigma, vv, ss,
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
  pp <- p_func(ga, p0)
  AA <- A_func(pp)
  BB <- B_func(pp)
  CC <- C_func(ga, pp)

  # empty stuff to save samples
  be_save <- array(NA, dim = c(ncol(XX), niter - nburn))
  vv_save <- array(NA, dim = c(nn, niter - nburn))
  ss_save <- array(NA, dim = c(nn, niter - nburn))
  sigma_save <- array(NA, dim = niter - nburn)
  ga_save <- array(NA, dim = niter - nburn)
  ga_acct_save <- array(0, dim = niter)

  #--------------------------
  # MCMC
  #--------------------------
  block_runtime <- 0

  for (iter in 1:niter) {

    start_time <- Sys.time()

    #--------------------------
    # update beta
    #--------------------------
    be <- updateBeta(be_mu0, be_Sigma0, yy, XX, sigma, ga, ss, vv, AA, BB, CC)
    mu <- XX %*% be

    #--------------------------
    # update vv
    #--------------------------
    vv <- updateVV(yy, mu, sigma, ga, ss, AA, BB, CC, nn)

    #--------------------------
    # update ss
    #--------------------------
    ss <- updateSS(yy, mu, sigma, ga, vv, AA, BB, CC, nn)

    #--------------------------
    # update sigma
    #--------------------------
    sigma <- updateSigma(u_sigma, v_sigma, yy, mu, ga, ss, vv, AA, BB, CC, nn)

    #--------------------------
    # update ga
    #--------------------------
    ga_res <- updateGa(ga, se_ga, u_ga, v_ga, yy, mu, sigma, p0, ss, vv, AA, BB, CC)
    ga_acct_save[iter] <- ga_res$accept
    ga <- ga_res$ga
    AA <- ga_res$AA
    BB <- ga_res$BB
    CC <- ga_res$CC

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

        cat(paste0("Metropolis sampler acceptance rates: ",
                   specRound(sum(ga_acct_save[1:iter]) / iter * 100, 2), "%\n\n"))

      }

    }

    #--------------------------
    # Save samples
    #--------------------------
    if (iter > nburn) {
      be_save[, iter - nburn] <- be
      vv_save[, iter - nburn] <- vv
      ss_save[, iter - nburn] <- ss
      sigma_save[iter - nburn] <- sigma
      ga_save[iter - nburn] <- ga
    }

  }

  #--------------------------
  # Thinning
  #--------------------------
  selc_idx <- seq(1, niter - nburn, by = nthin)

  post_sams <- list(be = be_save[, selc_idx],
                    vv = vv_save[, selc_idx],
                    ss = ss_save[, selc_idx],
                    sigma = sigma_save[selc_idx],
                    ga = ga_save[selc_idx])

  list(post_sams = post_sams, mh_acct_save = ga_acct_save)


}