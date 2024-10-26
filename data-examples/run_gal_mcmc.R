run_gal_mcmc <- function(yy, XX, p0, beta_prior, priors, 
                         starting, tuning, mcmc_settings, 
                         ga_sampler, verbose) {
  
  u_ga <- find_ga_lb(p0, c(-100, 0))
  v_ga <- find_ga_ub(p0, c(0, 100))
  
  priors[["ga_uniform"]] <- c(u_ga, v_ga)
  starting$ga <- (u_ga + v_ga) / 2
  
  bqrgal_out <- bgal(yy, 
                     XX, 
                     p0, 
                     beta_prior, 
                     priors, 
                     starting, 
                     tuning,
                     mcmc_settings, 
                     ga_sampler,
                     verbose)
  
  bqrgal_out
  
}
