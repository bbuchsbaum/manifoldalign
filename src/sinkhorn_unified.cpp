#include <RcppArmadillo.h>
#include <algorithm>
#include "ot_common.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Unified Sinkhorn algorithm with optional log-domain stabilization
//' 
//' @param cost Cost matrix (n x m)
//' @param a Source marginal (n x 1)
//' @param b Target marginal (m x 1)
//' @param epsilon Entropic regularization parameter
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param stabilized Use log-domain stabilization for numerical stability
//' @return Transport plan matrix
// [[Rcpp::export]]
arma::mat sinkhorn_unified(const arma::mat& cost,
                           const arma::vec& a,
                           const arma::vec& b,
                           double epsilon,
                           int max_iter = 1000,
                           double tol = 1e-9,
                           bool stabilized = true) {
  
  // Check marginals
  check_marginals(a, b);
  
  int n = a.n_elem;
  int m = b.n_elem;
  
  if (stabilized && epsilon < 0.1) {
    // Use log-domain stabilization for small epsilon
    
    // Initialize log-potentials
    arma::vec log_u(n, arma::fill::zeros);
    arma::vec log_v(m, arma::fill::zeros);
    
    // Pre-compute log-kernel
    arma::mat log_K = -cost / epsilon;
    double max_log_K = log_K.max();
    log_K -= max_log_K;  // Improve numerical stability
    
    // Sinkhorn iterations in log domain
    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec log_u_prev = log_u;
      
      // Update log_v
      arma::vec log_sum_u_K = log_sum_exp_cols(log_K, log_u);
      log_v = log(b) - log_sum_u_K;
      
      // Update log_u
      arma::vec log_sum_v_KT = log_sum_exp_rows(log_K, log_v);
      log_u = log(a) - log_sum_v_KT;
      
      // Check convergence
      if (arma::max(arma::abs(log_u - log_u_prev)) < tol) {
        break;
      }
    }
    
    // Reconstruct transport plan
    arma::mat log_P = arma::repmat(log_u, 1, m) + 
                      arma::repmat(log_v.t(), n, 1) + log_K;
    arma::mat P = arma::exp(log_P);
    
    // Ensure exact marginal constraints
    P = normalize_doubly_stochastic(P, a, b, tol/10, 20);
    
    return P;
    
  } else {
    // Use standard Sinkhorn for large epsilon or when not stabilized
    
    // Gibbs kernel
    arma::mat K = arma::exp(-cost / epsilon);
    
    // Initialize scaling factors
    arma::vec u = arma::ones(n);
    arma::vec v = arma::ones(m);
    
    // Sinkhorn iterations
    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec u_old = u;
      
      // Update scaling factors
      u = a / (K * v + 1e-20);
      v = b / (K.t() * u + 1e-20);
      
      // Check convergence
      double err = arma::norm(u - u_old, "inf");
      if (err < tol) {
        break;
      }
    }
    
    // Compute transport plan
    arma::mat P = arma::diagmat(u) * K * arma::diagmat(v);
    
    return P;
  }
}

// Maintain backward compatibility - redirect old functions to unified version

// [[Rcpp::export]]
arma::mat solve_sinkhorn_stabilized_cpp(const arma::mat& C_in, double tau, 
                                        int max_iter, double tol) {
  // For backward compatibility with parrot code
  int n = C_in.n_rows;
  int m = C_in.n_cols;
  arma::vec mu(n, arma::fill::value(1.0 / std::max(1, n)));
  arma::vec nu(m, arma::fill::value(1.0 / std::max(1, m)));
  
  return sinkhorn_unified(C_in, mu, nu, tau, max_iter, tol, true);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sinkhorn_ot_cpp(Rcpp::NumericMatrix cost_r,
                                    Rcpp::NumericVector a_r,
                                    Rcpp::NumericVector b_r,
                                    double epsilon,
                                    int max_iter = 1000,
                                    double tol = 1e-9) {
  // For backward compatibility with FPGW code
  arma::mat cost = Rcpp::as<arma::mat>(cost_r);
  arma::vec a = Rcpp::as<arma::vec>(a_r);
  arma::vec b = Rcpp::as<arma::vec>(b_r);
  
  arma::mat P = sinkhorn_unified(cost, a, b, epsilon, max_iter, tol, epsilon < 0.1);
  
  return Rcpp::wrap(P);
}