#include <RcppArmadillo.h>
#include "ot_common.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Note: log_sum_exp functions moved to ot_common.h

// Note: solve_sinkhorn_stabilized_cpp moved to sinkhorn_unified.cpp

// Note: normalize_doubly_stochastic_cpp moved to ot_common.h as normalize_doubly_stochastic

// Keep the export for backward compatibility
// [[Rcpp::export]]
arma::mat normalize_doubly_stochastic_cpp(arma::mat S, const arma::vec& mu, 
                                         const arma::vec& nu, double tol = 1e-9, 
                                         int max_iter = 20) {
  return normalize_doubly_stochastic(S, mu, nu, tol, max_iter);
}