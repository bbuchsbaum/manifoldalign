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
    const arma::vec log_a = arma::log(a);
    const arma::vec log_b = arma::log(b);
    
    // Pre-compute log-kernel
    arma::mat log_K = -cost / epsilon;
    double max_log_K = log_K.max();
    log_K -= max_log_K;  // Improve numerical stability
    
    arma::mat tmp(n, m, arma::fill::none);

    // Sinkhorn iterations in log domain
    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec log_u_prev = log_u;
      
      // Update log_v
      tmp = log_K;
      tmp.each_col() += log_u;
      arma::rowvec maxv = arma::max(tmp, 0);
      tmp.each_row() -= maxv;
      tmp = arma::exp(tmp);
      arma::rowvec sumexp_v = arma::sum(tmp, 0);
      sumexp_v += 1e-300;
      arma::vec log_sum_u_K = (maxv + arma::log(sumexp_v)).t();
      log_v = log_b - log_sum_u_K;
      
      // Update log_u
      tmp = log_K;
      tmp.each_row() += log_v.t();
      arma::colvec maxu = arma::max(tmp, 1);
      tmp.each_col() -= maxu;
      tmp = arma::exp(tmp);
      arma::colvec sumexp_u = arma::sum(tmp, 1);
      sumexp_u += 1e-300;
      arma::vec log_sum_v_KT = maxu + arma::log(sumexp_u);
      log_u = log_a - log_sum_v_KT;
      
      // Check convergence
      if (arma::max(arma::abs(log_u - log_u_prev)) < tol) {
        break;
      }
    }
    
    // Reconstruct transport plan
    tmp = log_K;
    tmp.each_col() += log_u;
    tmp.each_row() += log_v.t();
    arma::mat P = arma::exp(tmp);
    
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

//' Unified Sinkhorn algorithm returning dual potentials for warm-starting
//'
//' This variant returns both the transport plan and the (log) scaling vectors
//' that can be reused to warm-start subsequent Sinkhorn solves.
//'
//' @param cost Cost matrix (n x m)
//' @param a Source marginal (n x 1)
//' @param b Target marginal (m x 1)
//' @param epsilon Entropic regularization parameter
//' @param log_u0 Optional warm-start log scaling for rows (length n)
//' @param log_v0 Optional warm-start log scaling for columns (length m)
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param stabilized Use log-domain stabilization for small epsilon
//' @return A list with elements pi (plan), log_u, log_v
// [[Rcpp::export]]
Rcpp::List sinkhorn_unified_potentials(const arma::mat& cost,
                                      const arma::vec& a,
                                      const arma::vec& b,
                                      double epsilon,
                                      Rcpp::Nullable<Rcpp::NumericVector> log_u0 = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericVector> log_v0 = R_NilValue,
                                      int max_iter = 1000,
                                      double tol = 1e-9,
                                      bool stabilized = true) {

  check_marginals(a, b);

  const int n = a.n_elem;
  const int m = b.n_elem;

  arma::vec log_u(n, arma::fill::zeros);
  arma::vec log_v(m, arma::fill::zeros);

  if (log_u0.isNotNull()) {
    Rcpp::NumericVector u0(log_u0);
    if (u0.size() != n) Rcpp::stop("log_u0 must have length nrow(cost).");
    log_u = Rcpp::as<arma::vec>(u0);
  }
  if (log_v0.isNotNull()) {
    Rcpp::NumericVector v0(log_v0);
    if (v0.size() != m) Rcpp::stop("log_v0 must have length ncol(cost).");
    log_v = Rcpp::as<arma::vec>(v0);
  }

  arma::mat P;

  if (stabilized && epsilon < 0.1) {
    // Log-domain stabilization
    const arma::vec log_a = arma::log(a);
    const arma::vec log_b = arma::log(b);

    arma::mat log_K = -cost / epsilon;
    double max_log_K = log_K.max();
    log_K -= max_log_K;

    arma::mat tmp(n, m, arma::fill::none);

    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec log_u_prev = log_u;

      // Update log_v
      tmp = log_K;
      tmp.each_col() += log_u;
      arma::rowvec maxv = arma::max(tmp, 0);
      tmp.each_row() -= maxv;
      tmp = arma::exp(tmp);
      arma::rowvec sumexp_v = arma::sum(tmp, 0);
      sumexp_v += 1e-300;
      arma::vec log_sum_u_K = (maxv + arma::log(sumexp_v)).t();
      log_v = log_b - log_sum_u_K;

      // Update log_u
      tmp = log_K;
      tmp.each_row() += log_v.t();
      arma::colvec maxu = arma::max(tmp, 1);
      tmp.each_col() -= maxu;
      tmp = arma::exp(tmp);
      arma::colvec sumexp_u = arma::sum(tmp, 1);
      sumexp_u += 1e-300;
      arma::vec log_sum_v_KT = maxu + arma::log(sumexp_u);
      log_u = log_a - log_sum_v_KT;

      if (arma::max(arma::abs(log_u - log_u_prev)) < tol) {
        break;
      }
    }

    // Reconstruct plan
    tmp = log_K;
    tmp.each_col() += log_u;
    tmp.each_row() += log_v.t();
    P = arma::exp(tmp);

    // Enforce marginals (important for downstream usage)
    P = normalize_doubly_stochastic(P, a, b, tol / 10, 20);

  } else {
    // Standard Sinkhorn scaling (optionally warm-started via log_u/log_v)
    arma::mat K = arma::exp(-cost / epsilon);

    arma::vec u = arma::ones(n);
    arma::vec v = arma::ones(m);

    if (log_u0.isNotNull()) {
      double s = log_u.max();
      u = arma::exp(log_u - s);
    }
    if (log_v0.isNotNull()) {
      double s = log_v.max();
      v = arma::exp(log_v - s);
    }

    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec u_old = u;
      u = a / (K * v + 1e-20);
      v = b / (K.t() * u + 1e-20);
      if (arma::norm(u - u_old, "inf") < tol) {
        break;
      }
    }

    P = arma::diagmat(u) * K * arma::diagmat(v);

    log_u = arma::log(u + 1e-300);
    log_v = arma::log(v + 1e-300);
  }

  return Rcpp::List::create(
    Rcpp::_["pi"] = P,
    Rcpp::_["log_u"] = log_u,
    Rcpp::_["log_v"] = log_v
  );
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
