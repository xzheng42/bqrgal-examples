#' @title Prediction methods for GAL-model fits
#'
#' @description This function obtains out-of-sample predictions from a fitted GAL-based quantile regression model object
#' 
#' @param object a fitted object of class inheriting from "bgal".
#' @param covars an \eqn{m\times p} matrix of \eqn{p} covariates with which to predict. 
#'               If \code{err_dens} is true, this argument is ignored. 
#' @param probs a vector of quantile probabilities with values in \eqn{[0,1]} to generate the corresponding quantiles.
#' @param err_dens logical; if true, the \code{covars} argument is ignored, and samples or quantiles of the 
#'                 posterior predictive error density are returned.
#' @param predict_sam logical; if true, posterior predictive samples are returned, 
#'                    otherwise only quantiles (specified according to the probs argument) 
#'                    of the predictive distribution are returned.
#' @param ... additional arguments. (No additional arguments are supported currently)
#' 
#' @export
predict.bgal <- function(object, covars, probs = c(0.025, 0.5, 0.975),
                         err_dens = FALSE, predict_sam = FALSE, ...) {


  if (missing(object)) {
    stop("error: need bgal object for prediction\n")
  }
  if (!inherits(object, "bgal")) {
    stop("error: require bgal object for prediction\n")
  }

  post_sams <- object$post_sams

  if (!err_dens) {

    if (nrow(post_sams$be) != ncol(covars)) {
      stop("error: the dimension of the covariate space does not match the regression coefficients\n")
    }
    XX <- covars

  } else {

    pp <- nrow(post_sams$be)
    XX <- t(as.matrix(rep(0, pp)))

  }
  
  p0 <- object$prob

  pp <- p_func(post_sams$ga, p0)
  AA <- A_func(pp)
  BB <- B_func(pp)
  CC <- C_func(post_sams$ga, pp)

  nsams <- length(post_sams$ga)
  if (is.null(post_sams$be0)) {
    mu_sams <- sapply(1:nsams, function(x) XX %*% post_sams$be[, x])
  } else {
    mu_sams <- sapply(1:nsams, function(x) XX %*% post_sams$be[, x] + post_sams$be0[x])
  }

  nn <- nrow(XX)
  mm <- length(post_sams$sigma)
  ss <- matrix(truncnorm::rtruncnorm(nn * mm, 0, Inf), nrow = nn, ncol = mm)
  zz <- matrix(rexp(nn * mm), nrow = nn, ncol = mm)
  vv <- t(apply(zz, 1, function(x) x * post_sams$sigma))
  pred_out <- predBQRGAL_cpp(mu_sams, post_sams$ga, post_sams$sigma,
                             ss, vv, AA, BB, CC, probs, err_dens, predict_sam)

  pred_out

}

#' @title Prediction methods for AL-model fits
#'
#' @description This function obtains out-of-sample predictions from a fitted AL-based quantile regression model object
#' 
#' @param object a fitted object of class inheriting from "bgal".
#' @param covars an \eqn{m\times p} matrix of \eqn{p} covariates with which to predict. 
#'               If \code{err_dens} is true, this argument is ignored. 
#' @param probs a vector of quantile probabilities with values in \eqn{[0,1]} to generate the corresponding quantiles.
#' @param err_dens logical; if true, the \code{covars} argument is ignored, and samples or quantiles of the 
#'                 posterior predictive error density are returned.
#' @param predict_sam logical; if true, posterior predictive samples are returned, 
#'                    otherwise only quantiles (specified according to the probs argument) 
#'                    of the predictive distribution are returned.
#' @param ... additional arguments. (No additional arguments are supported currently)
#' 
#' @export
predict.bal <- function(object, covars, probs = c(0.025, 0.5, 0.975),
                        err_dens = FALSE, predict_sam = FALSE, ...) {


  if (missing(object)) {
    stop("error: need bal object for prediction\n")
  }
  if (!inherits(object, "bal")) {
    stop("error: require bal object for prediction\n")
  }

  post_sams <- object$post_sams

  if (!err_dens) {

    if (nrow(post_sams$be) != ncol(covars)) {
      stop("error: the dimension of the covariate space does not match the regression coefficients\n")
    }
    XX <- covars

  } else {

    pp <- nrow(post_sams$be)
    XX <- t(as.matrix(rep(0, pp)))

  }

  pp <- object$prob
  AA <- A_func(pp)
  BB <- B_func(pp)

  nn <- nrow(XX)
  mm <- length(post_sams$sigma)
  if (is.null(post_sams$be0)) {
    mu_sams <- sapply(1:mm, function(x) XX %*% post_sams$be[, x])
  } else {
    mu_sams <- sapply(1:mm, function(x) XX %*% post_sams$be[, x] + post_sams$be0[x])
  }

  zz <- matrix(rexp(nn * mm), nrow = nn, ncol = mm)
  vv <- t(apply(zz, 1, function(x) x * post_sams$sigma))
  pred_out <- predBQRAL_cpp(mu_sams, post_sams$sigma,
                            vv, AA, BB, probs, err_dens, predict_sam)

  pred_out

}
