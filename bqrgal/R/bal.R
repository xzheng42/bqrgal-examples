#' @title Function for fitting a Bayesian linear quantile regression with asymmetric Laplace (AL) errors 
#'        via Markov chain Monte Carlo (MCMC).
#'
#' @param resp a vector of response data of length \eqn{n}
#' @param covars an \eqn{n\times p} matrix of \eqn{p} covariates.
#' @param prob a quantile probability with value in (0,1). For example, if prob is 0.5, the model corresponds to a median quantile regression.
#' @param beta_prior a quoted keyword that specifies the prior for regression parameter \eqn{\beta}.
#'                   Supported key words are \code{gaussian} and \code{laplace}.
#' @param priors a list of priors. Each element of the list corresponds to a combination of an unknown parameter
#'               and its prior distribution; for example, "\code{phi_invgamma}".
#' @param starting a list of starting values. Each element of the list corresponds to an unknown parameter.
#' @param mcmc_settings a list of MCMC simulation parameters.
#'        \code{n_iter}: number of iterations; \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If   \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress.
#' @param verbose logical; if true, model specification and MCMC progress are printed to the screen.
#'
#' 
#' @return 
#' An object of class \code{bal}. The return object is a list that comprises:
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
#' \code{mcmc_settings}: \tab a list of MCMC simulation parameters:
#'        \code{n_iter}: number of iterations; 
#'        \code{n_burn}: number of burn-in samples;
#'        \code{n_thin}: thinning degree. If   \code{verbose = TRUE},
#'        \code{n_report} defines the interval to report MCMC progress and Metropolis acceptance rate (if \code{ga_sampler = metropolis}).
#'  \cr
#' \tab \cr
#' \code{runtime}: \tab running time for MCMC simulation calculated using \code{system.time}. \cr
#' 
#' }
#'
#' @author Xiaotian Zheng \email{xiaotianzheng42@gmail.com}
#' 
#' @references 
#' Kozumi, H. and Kobayashi, G. (2011), 
#' "Gibbs sampling methods for Bayesian quantile regression." 
#' Journal of Statistical Computation and Simulation, 81, 1565–1578.
#' DOI: \href{https://doi.org/10.1080/00949655.2010.496117}{https://doi.org/10.1080/00949655.2010.496117}.
#' 
#'@export
#'
bal <- function(resp,
                covars,
                prob,
                beta_prior,
                priors,
                starting,
                mcmc_settings,
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
    cat(" << Fitting a Bayesian quantile regression with AL errors >>\n")
    cat("--------------------------------------------------------------\n")
    cat(paste0("Number of observations: ", length(resp), "\n"))
    cat(paste0("The model corresponds to an AL distribution with probability ", p0, "\n"))

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
  mcmc_out <- bqral(yy, XX, p0, beta_prior, priors, starting, mcmc_settings, verbose)


  #----------------------------------------------------
  # output
  #----------------------------------------------------
  post_sams <- mcmc_out$post_sams
  runtime <- mcmc_out$runtim
  mh_data <- mcmc_out$mh_acct_save

  bal_out <- list(
    post_sams = post_sams,
    response = yy,
    covariates = XX,
    prob = p0,
    beta_prior = beta_prior,
    starting = starting,
    mcmc_settings = mcmc_settings,
    runtime = runtime
  )

  class(bal_out) <- "bal"

  bal_out

}


bqral <- function(yy, XX, p0, beta_prior, priors, starting, mcmc_settings, verbose) {


  #----------------------------------------------------
  # priors
  #----------------------------------------------------
  u_sigma <- priors$sigma_invgamma[1]
  v_sigma <- priors$sigma_invgamma[2]

  #----------------------------------------------------
  # Starting values
  #----------------------------------------------------
  sigma <- starting$sigma
  vv <- starting$vv

  #----------------------------------------------------
  # fit models
  #----------------------------------------------------
  # GAL model with Gaussian priors
  if (beta_prior == "gaussian") {

    be_mu0 <- priors$beta_gaus$mean_vec
    be_Sigma0 <- priors$beta_gaus$var_mat

    runtime <- system.time(
      mcmc_out <- mcmc_gaus_bal(yy, XX, p0,
                                be_mu0, be_Sigma0,
                                u_sigma, v_sigma,
                                sigma, vv,
                                mcmc_settings,
                                verbose)
    )
  }
  # GAL model with laplace priors
  else if (beta_prior == "laplace") {

    be0 <- starting$be0
    omega <- starting$omega

    be0_sigmasq0 <- priors$intercept_gaus_var
    u_eta <- priors$eta_gamma[1]
    v_eta <- priors$eta_gamma[2]

      runtime <- system.time(
        mcmc_out <- mcmc_gaus_bal_lasso(yy, XX, p0,
                                        be0_sigmasq0,
                                        u_eta, v_eta,
                                        u_sigma, v_sigma,
                                        be0, sigma, vv, omega,
                                        mcmc_settings,
                                        verbose)
      )

  }

  mcmc_out$runtime <- runtime

  mcmc_out

}
