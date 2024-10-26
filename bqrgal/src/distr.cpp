#include <iostream>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double ind(const double& a,
           const double& b) {

    double z;
    if (a < b) {
        z = 1.0;
    } else {
        z = 0.0;
    }

    return z;

}

// [[Rcpp::export]]
double gfunc_cpp(const double& ga) {

    double log_g = log(2.0) + R::pnorm(-abs(ga), 0, 1, true, true) + .5 * ga * ga;

    return exp(log_g);
}


// [[Rcpp::export]]
arma::colvec dgal_cpp(const arma::colvec& y,
                      const double& p0,
                      const double& mu,
                      const double& sigma,
                      const double& ga,
                      const bool& logd = false) {

  double p, p_ga1, p_ga0, log_scal, y_star, p_ratio;
  double term1, log_term2, log_term3, log_term4;
  double t1q, t3q;
  int nn = y.n_rows;
  arma::colvec log_dens(nn);


  if (ga == 0) {

    p = p0;
    log_scal = log(p) + log(1.0 - p) - log(sigma);

    for (int i = 0; i < nn; ++i) {

      log_dens(i) = log_scal - (y(i) - mu) * (p - ind(y(i), mu)) / sigma;

    }

  } else {

    p = ind(ga, 0) + (p0 - ind(ga, 0)) / gfunc_cpp(ga);
    log_scal = log(2.0) + log(p) + log(1.0 - p) - log(sigma);
    p_ga1 = p - ind(0, ga);
    p_ga0 = p - ind(ga, 0);

    for (int i = 0; i < nn; ++i) {

      y_star = (y(i) - mu) / sigma;
      log_term4 = -p_ga1 * y_star + .5 * ga * ga;

      if (y_star / ga > 0) {

        p_ratio = p_ga0 / p_ga1;

        t1q = -p_ga1 * y_star / abs(ga) + p_ratio * abs(ga);
        term1 = R::pnorm(t1q, 0, 1, true, false) -
          R::pnorm(p_ratio * abs(ga), 0, 1, true, false);
        log_term2 = -p_ga0 * y_star + .5 * ga * ga * p_ratio * p_ratio;
        t3q = -abs(ga) + p_ga1 * y_star / abs(ga);
        log_term3 = R::pnorm(t3q, 0, 1, true, true);

        log_dens(i) = log_scal + log(exp(log(term1) + log_term2) + exp(log_term3 + log_term4));

      } else {

        log_term3 = R::pnorm(-abs(ga), 0, 1, true, true);
        log_dens(i) = log_scal + log_term3 + log_term4;

      }

    }

  }

  if (logd) {
    return log_dens;
  } else {
    return arma::exp(log_dens);
  }

}


