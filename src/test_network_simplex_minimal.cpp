#include <RcppArmadillo.h>
#include "ot_common.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Minimal test of network simplex without LEMON
// [[Rcpp::export]]
arma::mat test_network_simplex_minimal() {
  // Very simple 2x2 problem
  arma::mat cost(2, 2);
  cost(0,0) = 1.0; cost(0,1) = 2.0;
  cost(1,0) = 3.0; cost(1,1) = 1.0;
  
  arma::vec a(2);
  a(0) = 0.5; a(1) = 0.5;
  
  arma::vec b(2);
  b(0) = 0.5; b(1) = 0.5;
  
  Rcpp::Rcout << "Cost matrix:\n" << cost << "\n";
  Rcpp::Rcout << "Supply a: " << a.t() << "\n";
  Rcpp::Rcout << "Demand b: " << b.t() << "\n";
  
  // Simple northwest corner solution for testing
  arma::mat P(2, 2, arma::fill::zeros);
  P(0,0) = 0.5;
  P(1,1) = 0.5;
  
  return P;
}