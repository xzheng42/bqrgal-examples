simGausQr <- function(n, p0, mu, sigma) {

  mu0 <- -sigma * qnorm(p0)
  mu + rnorm(n, mu0, sigma)

}

simLaplaceQr <- function(n, p0, mu, sigma) {

  if (p0 <= 0.5) {
    mu0 <- -sigma * log(2 * p0)
  } else {
    mu0 <- sigma * log(2 * (1 - p0))
  }

  mu + nimble::rdexp(n, mu0, sigma)

}

simGausMixQr <- function(n, p0, mu, sigma) {
  
  findGausMixMu <- function(x0, p0, sigma) {
    0.1 * pnorm(0, x0, sigma[1]) + 0.9 * pnorm(0, x0 + 1, sigma[2]) - p0
  }
  
  sol <- uniroot(findGausMixMu, c(-100, 100), p0, sigma)
  
  lab <- sample(1:2, n, replace = TRUE, prob = c(0.1, 0.9))
  
  mu0 <- c(sol$root, sol$root + 1)
  loc_vec <- mu0[lab]
  scal_vec <- sigma[lab]
  
  mu + rnorm(n, loc_vec, scal_vec)
  
}

simLogGPD <- function(n, p0, mu, xi) {

  sigma0 <- xi * (1-p0)^xi / (1 - (1 - p0)^xi)

  mu + log(SpatialExtremes::rgpd(n, loc = 0, scale = sigma0, shape = xi))

}
