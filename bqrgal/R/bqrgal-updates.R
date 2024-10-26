updateBeta <- function(be_mu0, be_Sigma0, yy, XX, sigma, ga, ss, vv,
                       AA, BB, CC) {

  be_Sigma0_inv <- chol2inv(chol(be_Sigma0))
  XX_vv <- apply(XX, 2, function(x) x / sqrt(BB * sigma * vv))
  be_Sigma1 <- chol2inv(chol(be_Sigma0_inv + t(XX_vv) %*% XX_vv))
  weights <- (yy - sigma * CC * abs(ga) * ss - AA * vv) / (BB * sigma * vv)
  resids <- apply(XX, 2, function(x) x * weights)
  be_mu1 <- be_Sigma1 %*% (be_Sigma0_inv %*% be_mu0 + colSums(resids))

  t(mvtnorm::rmvnorm(1, be_mu1, be_Sigma1))

}


updateBeta0 <- function(be0_sigmasq0, yy, XX, be, sigma, ga, ss, vv, AA, BB, CC) {

  be0_sigmasq1 <- 1 / (1 / be0_sigmasq0 + sum(1 / vv) / (BB * sigma))
  resid_sum <- sum((yy - XX %*% be - sigma * CC * abs(ga) * ss - AA * vv) / vv)
  be0_mu1 <- be0_sigmasq1 * resid_sum / (BB * sigma)

  rnorm(1, be0_mu1, be0_sigmasq1)

}


updateBetaLasso <- function(Omega, yy, XX, be0, sigma, ga, ss, vv, AA, BB, CC) {

  XX_vv <- apply(XX, 2, function(x) x / sqrt(BB * sigma * vv))
  be_Sigma1 <- chol2inv(chol(Omega + t(XX_vv) %*% XX_vv))
  weights <- (yy - be0 - sigma * CC * abs(ga) * ss - AA * vv) / (BB * sigma * vv)
  resids <- apply(XX, 2, function(x) x * weights)
  be_mu1 <- be_Sigma1 %*% (colSums(resids))

  t(mvtnorm::rmvnorm(1, be_mu1, be_Sigma1))

}


updateSigma <- function(u_sigma, v_sigma, yy, mu, ga, ss, vv,
                        AA, BB, CC, nn) {

  nu <- -(u_sigma + 1.5 * nn)
  cc <- 2 * v_sigma + 2 * sum(vv) + sum((yy - mu - AA * vv)^2 / vv) / BB
  dd <- sum((CC * ga * ss)^2 / vv) / BB
  GIGrvg::rgig(1, nu, cc, dd)

}


updateSS <- function(yy, mu, sigma, ga, vv, AA, BB, CC, nn) {

  sigma2_ss <- 1 / (1 + (CC * ga)^2 * sigma / (BB * vv))
  mu_ss <- sigma2_ss * CC * abs(ga) * (yy - mu - AA * vv) / (BB * vv)

  truncnorm::rtruncnorm(nn, rep(0, nn), rep(Inf, nn), mu_ss, sqrt(sigma2_ss))

}


updateVV <- function(yy, mu, sigma, ga, ss, AA, BB, CC, nn) {

  a_i <- (yy - mu - sigma * CC * abs(ga) * ss)^2 / (BB * sigma)
  b_i <- 2 / sigma + AA^2 / (BB * sigma)
  vv <- rep(NA, nn)
  for (i in 1:nn) {
    vv[i] <- GIGrvg::rgig(1, 0.5, a_i[i], b_i)
  }
  vv
}


updateGa <- function(ga, se_ga, u_ga, v_ga,
                     yy, mu, sigma, p0, ss, vv, AA, BB, CC) {

  logit_ga <- glogit(ga, u_ga, v_ga)
  prop_logit_ga <- rnorm(1, logit_ga, se_ga)
  prop_ga <- gexpit(prop_logit_ga, u_ga, v_ga)

  prop_pp <- p_func(prop_ga, p0)
  prop_AA <- A_func(prop_pp)
  prop_BB <- B_func(prop_pp)
  prop_CC <- C_func(prop_ga, prop_pp)

  prop_yy_mean <- mu + sigma * prop_CC * abs(prop_ga) * ss + prop_AA * vv
  prop_yy_sigma <- sqrt(sigma * prop_BB * vv)
  cur_yy_mean <- mu + sigma * CC * abs(ga) * ss + AA * vv
  cur_yy_sigma <- sqrt(sigma * BB * vv)

  prop_loglik <- sum(dnorm(yy, prop_yy_mean, prop_yy_sigma, log = TRUE)) +
    log(prop_ga - u_ga) + log(v_ga - prop_ga)
  cur_loglik <- sum(dnorm(yy, cur_yy_mean, cur_yy_sigma, log = TRUE)) +
    log(ga - u_ga) + log(v_ga - ga)
  diff <- prop_loglik - cur_loglik

  if (diff > log(runif(1))) {
    return(list(ga = prop_ga, AA = prop_AA, BB = prop_BB, CC = prop_CC, accept = TRUE))
  } else {
    return(list(ga = ga, AA = AA, BB = BB, CC = CC, accept = FALSE))
  }

}


updateOmega <- function(kk, be, etasq) {

  omega <- rep(NA, kk)
  besq <- be^2

  for (j in 1:kk) {
    omega[j] <-  GIGrvg::rgig(1, 0.5, besq[j], etasq)
  }

  omega

}


#############################################################################
### slice sampler
### the code is obtained from the https://www.cs.toronto.edu/~radford/ftp/slice-R-prog.
### the code was developed by Radford Neal, and is available for free use;
### see https://www.cs.toronto.edu/~radford/software-online.html.
#############################################################################
uniSlice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, ...)
{

  uni.slice.calls <- 0	# Number of calls of the slice sampling function
  uni.slice.evals <- 0	# Number of density evaluations done in these calls

  # Check the validity of the arguments.

  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g)
      || !is.numeric(w) || length(w)!=1 || w<=0
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower
      # || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1)
  )
  {
    stop ("Invalid slice sampling argument")
  }

  # Keep track of the number of calls made to this function.

  uni.slice.calls <<- uni.slice.calls + 1

  # Find the log density at the initial point, if not already known.

  # if (is.null(gx0))
  # { uni.slice.evals <<- uni.slice.evals + 1
  # gx0 <- g(x0)
  # }
  gx0 <- g(x0, ...)

  # Determine the slice level, in log terms.

  logy <- gx0 - rexp(1)

  # Find the initial interval to sample from.

  u <- runif(1, 0, w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff

  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.

  if (is.infinite(m))  # no limit on number of steps
  {
    repeat
    { if (L<=lower) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(L, ...) <= logy) break
      L <- L - w
    }

    repeat
    { if (R>=upper) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(R, ...) <= logy) break
      R <- R + w
    }
  }

  else if (m>1)  # limit on steps, bigger than one
  {
    J <- floor(runif(1,0,m))
    K <- (m-1) - J

    while (J>0)
    { if (L<=lower) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(L, ...)<=logy) break
      L <- L - w
      J <- J - 1
    }

    while (K>0)
    { if (R>=upper) break
      uni.slice.evals <<- uni.slice.evals + 1
      if (g(R, ...)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }


  # Shrink interval to lower and upper bounds.

  if (L<lower)
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }

  # Sample from the interval, shrinking it on each rejection.


  repeat
  {
    x1 <- runif(1, L, R)

    uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1, ...)

    if (gx1 >= logy) break

    if (x1 > x0)
    { R <- x1
    }
    else
    { L <- x1
    }
  }

  # Return the point sampled, with its log density attached as an attribute.
  #attr(x1,"log.density") <- gx1
  return (x1)
}
