#ifndef OT_COMMON_H
#define OT_COMMON_H

#include <RcppArmadillo.h>

// Common utility functions for optimal transport

// Check if marginals are valid for OT
inline void check_marginals(const arma::vec& a, const arma::vec& b, double tol = 1e-9) {
  double sum_a = arma::sum(a);
  double sum_b = arma::sum(b);
  
  if (std::abs(sum_a - sum_b) > tol) {
    Rcpp::stop("Marginals must have equal mass");
  }
  
  if (arma::any(a < 0) || arma::any(b < 0)) {
    Rcpp::stop("Marginals must be non-negative");
  }
}

// Normalize matrix to be doubly stochastic
inline arma::mat normalize_doubly_stochastic(arma::mat S, 
                                            const arma::vec& mu, 
                                            const arma::vec& nu, 
                                            double tol = 1e-9, 
                                            int max_iter = 20) {
  for (int iter = 0; iter < max_iter; ++iter) {
    // Row normalization
    arma::vec row_sums = arma::sum(S, 1);
    row_sums.replace(0, 1e-16);  // Avoid division by zero
    S.each_col() %= (mu / row_sums);
    
    // Column normalization
    arma::rowvec col_sums = arma::sum(S, 0);
    col_sums.replace(0, 1e-16);
    S.each_row() %= (nu.t() / col_sums);
    
    // Check convergence
    double row_error = arma::max(arma::abs(arma::sum(S, 1) - mu));
    double col_error = arma::max(arma::abs(arma::sum(S, 0).t() - nu));
    
    if (row_error < tol && col_error < tol) {
      break;
    }
  }
  
  return S;
}

// Log-sum-exp for numerical stability
inline double log_sum_exp(const arma::vec& x) {
  double max_val = x.max();
  return max_val + std::log(arma::sum(arma::exp(x - max_val)));
}

// Log-sum-exp for columns
inline arma::vec log_sum_exp_cols(const arma::mat& log_K, const arma::vec& log_u) {
  int n_cols = log_K.n_cols;
  arma::vec result(n_cols);
  
  for (int j = 0; j < n_cols; ++j) {
    arma::vec temp = log_u + log_K.col(j);
    result(j) = log_sum_exp(temp);
  }
  
  return result;
}

// Log-sum-exp for rows
inline arma::vec log_sum_exp_rows(const arma::mat& log_K, const arma::vec& log_v) {
  int n_rows = log_K.n_rows;
  arma::vec result(n_rows);
  
  for (int i = 0; i < n_rows; ++i) {
    arma::vec temp = log_v + log_K.row(i).t();
    result(i) = log_sum_exp(temp);
  }
  
  return result;
}

#endif // OT_COMMON_H