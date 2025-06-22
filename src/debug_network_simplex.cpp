#include <RcppArmadillo.h>
#include "ot_common.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Debug version without LEMON
// [[Rcpp::export]] 
arma::mat debug_network_simplex(arma::mat cost, arma::vec a, arma::vec b) {
  Rcpp::Rcout << "Debug network simplex\n";
  Rcpp::Rcout << "Cost dimensions: " << cost.n_rows << " x " << cost.n_cols << "\n";
  Rcpp::Rcout << "Supply a size: " << a.n_elem << ", sum: " << arma::sum(a) << "\n";
  Rcpp::Rcout << "Demand b size: " << b.n_elem << ", sum: " << arma::sum(b) << "\n";
  
  // Check marginals
  try {
    check_marginals(a, b);
    Rcpp::Rcout << "Marginals check passed\n";
  } catch (std::exception& e) {
    Rcpp::Rcout << "Marginals check failed: " << e.what() << "\n";
  }
  
  // Simple northwest corner for now
  int n = a.n_elem;
  int m = b.n_elem;
  arma::mat P(n, m, arma::fill::zeros);
  
  arma::vec a_rem = a;
  arma::vec b_rem = b;
  
  int i = 0, j = 0;
  while (i < n && j < m) {
    double mass = std::min(a_rem(i), b_rem(j));
    P(i, j) = mass;
    a_rem(i) -= mass;
    b_rem(j) -= mass;
    
    if (a_rem(i) < 1e-10) i++;
    if (b_rem(j) < 1e-10) j++;
  }
  
  Rcpp::Rcout << "Northwest corner solution sum: " << arma::sum(arma::sum(P)) << "\n";
  
  return P;
}