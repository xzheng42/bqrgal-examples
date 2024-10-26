#' @title The generalized asymmetric distribution
#'
#' @description Density for the four-parameter generalized asymmetric distribution with quantile probability \code{p0}, 
#'              location parameter \code{mu}, scale parameter \code{sigma}, and skewness parameter \code{ga}.
#' 
#' @param y vector of quantiles.
#' @param p0 numeric value between 0 and 1.
#' @param mu real-valued location parameter.
#' @param sigma positive-valued scale parameter.
#' @param ga numeric value bounded between L and U, where L is the negative root of 
#'           g(ga) = 1-p0 and U is the positive root of g(ga) = p0.
#' @param logd logical; if true, then densities d are given as log(d).
#' 
#' @export
dgal <- function(y, p0, mu, sigma, ga, logd = FALSE) {

  dd <- dgal_cpp(y, p0, mu, sigma, ga, logd)
  dd <- as.vector(dd)
 
  if (any(is.na(dd))) {
    warning("warning: the return value is NaN. 
             Make sure the argument ga is specified between L and U. 
             Use functions find_ga_lb and find_ga_ub to find L and U.")
  }

  dd

}


#' @title Find the negative root of \eqn{g(\gamma) = 1 - p_0}
#'
#' @description This function finds the negative root of \eqn{g(\gamma) = 1 - p_0} using R's uniroot function.
#' 
#' @param p0 numeric value between 0 and 1.
#' @param interval vector containing the end-points of the interval to be searched for the root.
#' @param ... additional arguments to be passed to uniroot
#' 
#' @export
find_ga_lb <- function(p0, interval, ...) {

  find_L <- function(ga, p0) ga_func(ga) - (1 - p0)

  uniroot(f = find_L, interval = interval, p0, ...)$root

}


#' @title Find the positive root of \eqn{g(\gamma) = p_0}
#'
#' @description This function finds the positive root of \eqn{g(\gamma) = p_0} using R's uniroot function.
#' 
#' @param p0 numeric value between 0 and 1.
#' @param interval vector containing the end-points of the interval to be searched for the root.
#' @param ... additional arguments to be passed to uniroot
#' 
#' @export
find_ga_ub <- function(p0, interval, ...) {

  find_U <- function(ga, p0) ga_func(ga) - p0

  uniroot(f = find_U, interval = interval, p0, ...)$root

}


ga_func <- function(ga, log = FALSE) {
  log_g <- log(2) + pnorm(-abs(ga), log.p = TRUE) + .5 * ga^2
  if (log) {
    log_g
  } else {
    exp(log_g)
  }
}


p_func <- function(ga, p0) {
  ind <- as.numeric(ga < 0)
  p <- ind + (p0 - ind) / ga_func(ga)
  p
}

A_func <- function(p) {
  (1 - 2 * p) / (p * (1 - p))
}

B_func <- function(p) {
  2 / (p * (1 - p))
}

C_func <- function(ga, p) {
  ind <- as.numeric(ga > 0)
  1 / (ind - p)
}


logcondpga <- function(ga, yy, p0, mu, sigma, ss, vv) {

  pp <- p_func(ga, p0)
  AA <- A_func(pp)
  BB <- B_func(pp)
  CC <- C_func(ga, pp)

  yy_mean <- mu + sigma * CC * abs(ga) * ss + AA * vv
  yy_sigma <- sqrt(sigma * BB * vv)
  sum(dnorm(yy, yy_mean, yy_sigma, log = TRUE))

}