#include <iostream>
#include <RcppArmadillo.h>
#include <truncnorm.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

//////////////////////////////////////////////////////////////////////////
// Posterior prediction
/////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List predBQRGAL_cpp(const arma::mat& mu,
                    const arma::colvec& ga,
                    const arma::colvec& sigma,
                    const arma::mat& ss,
                    const arma::mat& vv,
                    const arma::colvec& AA,
                    const arma::colvec& BB,
                    const arma::colvec& CC,
                    const arma::colvec& probs,
                    const bool err = false,
                    const bool sam = false) {

  int data_size = mu.n_rows;
  int sample_size = ga.n_rows;
  int probs_size = probs.n_rows;

  double ijmu, tausq;
  arma::mat pred_sam(data_size, sample_size);
  arma::mat pred_yy_qq(probs_size, data_size);
  arma::colvec ipred_yy(sample_size), pred_yy_mean(data_size);
  arma::rowvec iXX;

  for (int i = 0; i < data_size; ++i) {

    for (int j = 0; j < sample_size; ++j) {

      ijmu = arma::as_scalar(sigma(j) * CC(j) * abs(ga(j)) * ss(i,j) + vv(i,j) * AA(j));
      tausq = sigma(j) * vv(i,j) * BB(j);
      if (!err) {
        ijmu += mu(i, j);
      }
      ipred_yy(j) = R::rnorm(ijmu, sqrt(tausq));

      if (sam) {
        pred_sam(i, j) = ipred_yy(j);
      }

    }

    pred_yy_mean(i) = mean(ipred_yy);
    pred_yy_qq.col(i) = arma::quantile(ipred_yy, probs);

  }

  if (sam) {
    return(List::create(Named("yy_mu") = pred_yy_mean,
                        Named("yy_qq") = pred_yy_qq,
                        Named("yy_sam") = pred_sam));
  } else {
    return(List::create(Named("yy_mu") = pred_yy_mean,
                        Named("yy_qq") = pred_yy_qq));
  }

}


// [[Rcpp::export]]
List predBQRAL_cpp(const arma::mat& mu,
                   const arma::colvec& sigma,
                   const arma::mat& vv,
                   const double& AA,
                   const double& BB,
                   const arma::colvec& probs,
                   const bool err = false,
                   const bool sam = false) {

  int data_size = mu.n_rows;
  int sample_size = sigma.n_rows;
  int probs_size = probs.n_rows;

  double ijmu, tausq;
  arma::mat pred_sam(data_size, sample_size);
  arma::mat pred_yy_qq(probs_size, data_size);
  arma::colvec ipred_yy(sample_size), pred_yy_mean(data_size);
  arma::rowvec iXX;

  for (int i = 0; i < data_size; ++i) {

    for (int j = 0; j < sample_size; ++j) {

      ijmu = arma::as_scalar(vv(i,j) * AA);
      tausq = sigma(j) * vv(i,j) * BB;
      if (!err) {
        ijmu += mu(i, j);
      }
      ipred_yy(j) = R::rnorm(ijmu, sqrt(tausq));

      if (sam) {
        pred_sam(i, j) = ipred_yy(j);
      }

    }

    pred_yy_mean(i) = mean(ipred_yy);
    pred_yy_qq.col(i) = arma::quantile(ipred_yy, probs);

  }

  if (sam) {
    return(List::create(Named("yy_mu") = pred_yy_mean,
                        Named("yy_qq") = pred_yy_qq,
                        Named("yy_sam") = pred_sam));
  } else {
    return(List::create(Named("yy_mu") = pred_yy_mean,
                        Named("yy_qq") = pred_yy_qq));
  }

}
