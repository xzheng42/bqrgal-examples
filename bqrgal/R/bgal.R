#' @title Function for fitting a Bayesian linear quantile regression with generalized asymmetric Laplace (GAL) errors
#'        via Markov chain Monte Carlo (MCMC).
#'
#' @param resp a vector of response data of length \eqn{n}.
#' @param covars an \eqn{n\times p} matrix of \eqn{p} covariates.
#' @param prob a quantile probability with value in (0,1). For example, if prob is 0.5, the model corresponds to a median quantile regression.
#' @param beta_prior a quoted keyword that specifies the prior for regression parameter \eqn{\beta}.
#'                   Supported key words are \code{gaussian} and \code{laplace}.
#' @param priors a list of priors. Each element of the list corresponds to a combination of an unknown parameter
#'               and its prior distribution; for example, "\code{phi_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param tuning a list of tuning parameter for sampling \eqn{\gamma} (the skewness parameter of GAL) in MCMC. 
#'               If \code{ga_sampler = metropolis}, use \code{se_ga}: the standard deviation of the random
#'               walk Metropolis sampler proposal distribution. If \code{ga_sampler} = slice, use \code{step_size}: the step size or initial interval width 
#'               of the stepping out procedure of the slice sampler.
#' @param mcmc_settings a list of MCMC simulation parameters.
#'        \code{n_iter}: number of iterations; \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If   \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress and tropolis acceptance rate (if \code{ga_sampler = metropolis}).
#' @param ga_sampler a quoted keyword that specifies the MCMC sampler for the skewness parameter \eqn{\gamma}.
#'                   Supported key words are \code{metropolis} and \code{slice}. The default is \code{slice}.
#' @param verbose logical; if true, model specification, MCMC progress, and Metropolis acceptance rate (if \code{ga_sampler = metropolis})
#'                are printed to the screen.
#' 
#' @return 
#' An object of class \code{bgal}. The return object is a list that comprises:
#' 
#' \tabular{ll}{
#' 
#' \code{post_sams} \tab a list of posterior samples. Each element of the list corresponds to an unknown parameter. \cr
#' \tab \cr
#' \code{response} \tab the response vector. \cr
#' \tab \cr
#' \code{covariates} \tab the covariate matrix. \cr
#' \tab \cr
#' \code{prob} \tab the quantile probability. \cr
#' \tab \cr
#' \code{beta_prior} \tab the prior for regression parameter \eqn{\beta}. \cr
#' \tab \cr
#' \code{starting} \tab a list of starting values. \cr
#' \tab \cr
#' \code{tuning} \tab a list of tuning parameters. \cr
#' \tab \cr
#' \code{mcmc_settings}: \tab a list of MCMC simulation parameters:
#'        \code{n_iter}: number of iterations; 
#'        \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If   \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rate (if \code{ga_sampler = metropolis}).
#'  \cr
#' \tab \cr
#' \code{mh_data}: \tab Null if \code{ga_sampler = slice}; 
#'                      a vector of Metropolis sampler sampler acceptance and rejection indicators if \code{ga_sampler = metropolis}. \cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}. \cr
#' 
#' }
#'
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Yan, Y., Zheng, X., and Kottas, A. (2024). 
#' "A new family of error distributions for Bayesian quantile regression".
#' \href{https://tr.soe.ucsc.edu/research/technical-reports/UCSC-SOE-24-01}{https://tr.soe.ucsc.edu/research/technical-reports/UCSC-SOE-24-01}.
#' 
#' @export
#'
bgal <- function(resp,
                 covars,
                 prob,
                 beta_prior,
                 priors,
                 starting,
                 tuning,
                 mcmc_settings,
                 ga_sampler = "slice",
                 verbose = TRUE) {

  yy <- as.numeric(resp)
  XX <- as.matrix(covars)
  p0 <- prob

  #----------------------------------------------------
  # check prob
  #----------------------------------------------------
  if (prob <= 0 | prob >= 1) {
    stop("error: the argument prob must be strictly between 0 and 1")
  }

  #----------------------------------------------------
  # verbose
  #----------------------------------------------------
  if (verbose) {

    cat("--------------------------------------------------------------\n")
    cat(" << Fitting a Bayesian linear quantile regression with GAL errors >>\n")
    cat("--------------------------------------------------------------\n")
    cat(paste0("Number of observations: ", length(resp), "\n"))
    cat(paste0("The model corresponds to a quantile-fixed GAL distribution with quantile probability ", p0, "\n"))

    cat("--------------------------------------------\n")
    cat("\t  MCMC settings\n");
    cat("--------------------------------------------\n")
    cat(paste0("Number of iterations: ", mcmc_settings$n_iter, "\n"))
    cat(paste0("Number of burn-in samples: ", mcmc_settings$n_burn, "\n"))
    cat(paste0("Thinning degree: ", mcmc_settings$n_thin, "\n"))

  }

  #----------------------------------------------------
  # MCMC
  #----------------------------------------------------
  mcmc_out <- bqrgal(yy, XX, p0, beta_prior, priors, starting, tuning,
                     mcmc_settings, ga_sampler, verbose)


  #----------------------------------------------------
  # output
  #----------------------------------------------------
  post_sams <- mcmc_out$post_sams
  runtime <- mcmc_out$runtim
  mh_data <- mcmc_out$mh_acct_save

  bgal_out <- list(
    post_sams = post_sams,
    response = yy,
    covariates = XX,
    prob = p0,
    beta_prior = beta_prior,
    starting = starting,
    tuning = tuning,
    mcmc_settings = mcmc_settings,
    mh_data = mh_data,
    runtime = runtime
  )

  class(bgal_out) <- "bgal"

  bgal_out

}


bqrgal <- function(yy, XX, p0, beta_prior, priors, starting, tuning,
                   mcmc_settings, ga_sampler, verbose) {


  #----------------------------------------------------
  # priors
  #----------------------------------------------------
  u_sigma <- priors$sigma_invgamma[1]
  v_sigma <- priors$sigma_invgamma[2]

  u_ga <- priors$ga_uniform[1]
  v_ga <- priors$ga_uniform[2]

  #----------------------------------------------------
  # Starting values
  #----------------------------------------------------
  ga <- starting$ga
  sigma <- starting$sigma
  vv <- starting$vv
  ss <- starting$ss
  omega <- starting$omega

  #----------------------------------------------------
  # fit models
  #----------------------------------------------------
  # GAL model with Gaussian priors
  if (beta_prior == "gaussian") {

    be_mu0 <- priors$beta_gaus$mean_vec
    be_Sigma0 <- priors$beta_gaus$var_mat

    # GAL model using Metropolis to sample gamma
    if (ga_sampler == "metropolis") {
      se_ga <- tuning$se_ga

      runtime <- system.time(
      mcmc_out <- mcmc_gaus_bgal_mh(yy, XX, p0,
                                    be_mu0, be_Sigma0,
                                    u_sigma, v_sigma,
                                    u_ga, v_ga, se_ga,
                                    ga, sigma, vv, ss,
                                    mcmc_settings,
                                    verbose)
      )
    }
    # GAL model using slice sampling to sample gamma
    else if (ga_sampler == "slice") {

      # Check step size
      if (is.null(tuning$step_size)) {
        step_size <- 0.1
      } else{
        step_size <- tuning$step_size
      }

      runtime <- system.time(
        mcmc_out <- mcmc_gaus_bgal(yy, XX, p0,
                                   be_mu0, be_Sigma0,
                                   u_sigma, v_sigma,
                                   u_ga, v_ga, se_ga,
                                   ga, sigma, vv, ss,
                                   mcmc_settings, step_size,
                                   verbose)
      )
    } else {

    stop("error: the specified sampler for gamma is not supported.
           Available options are metropolis and slice.")
    }
  }
  # GAL model with laplace priors
  else if (beta_prior == "laplace") {

    be0 <- starting$be0
    be0_sigmasq0 <- priors$intercept_gaus_var
    u_eta <- priors$eta_gamma[1]
    v_eta <- priors$eta_gamma[2]

    # GAL model using Metropolis to sample gamma
    if (ga_sampler == "metropolis") {

      se_ga <- tuning$se_ga

      runtime <- system.time(
        mcmc_out <- mcmc_gaus_bgal_lasso_mh(yy, XX, p0,
                                            be0_sigmasq0,
                                            u_eta, v_eta,
                                            u_sigma, v_sigma,
                                            u_ga, v_ga, se_ga,
                                            be0, ga, sigma, vv, ss, omega,
                                            mcmc_settings,
                                            verbose)
      )
    }
    # GAL model using slice sampling to sample gamma
    else if (ga_sampler == "slice") {

      if (is.null(tuning$step_size)) {
        step_size <- 0.1
      } else{
        step_size <- tuning$step_size
      }

      runtime <- system.time(
        mcmc_out <- mcmc_gaus_bgal_lasso(yy, XX, p0,
                                         be0_sigmasq0,
                                         u_eta, v_eta,
                                         u_sigma, v_sigma,
                                         u_ga, v_ga, se_ga,
                                         be0, ga, sigma, vv, ss, omega,
                                         mcmc_settings, step_size,
                                         verbose)
      )
    } else {

      stop("error: the specified sampler for gamma is not supported.
           Available options are metropolis and slice.")

    }

  }
  else {

    stop("error: this family of prior for beta is not supported.")

  }

  mcmc_out$runtime <- runtime

  mcmc_out

}