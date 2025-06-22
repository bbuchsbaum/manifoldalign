#include "ot_solvers.h"
#include "ot_common.h"
#include "network_simplex_adapter.h"
#include <queue>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// Network simplex for optimal transport
// Now uses LEMON-based implementation when available
arma::mat network_simplex_ot(const arma::mat& cost, 
                              const arma::vec& a, 
                              const arma::vec& b,
                              double eps) {
  return NetworkSimplexAdapter::solve(cost, a, b, eps);
}

// Note: Sinkhorn algorithm moved to sinkhorn_unified.cpp for consolidation

// Mass-constrained partial OT
// Proper implementation using network simplex with dummy sink
arma::mat partial_ot_mass_cpp(const arma::mat& cost,
                               const arma::vec& a,
                               const arma::vec& b,
                               double mass,
                               double eps) {
  
  int n = a.n_elem;
  int m = b.n_elem;
  
  // If requested mass is >= total mass, solve standard OT
  double total_mass_a = arma::sum(a);
  double total_mass_b = arma::sum(b);
  if (mass >= std::min(total_mass_a, total_mass_b) - eps) {
    return network_simplex_ot(cost, a, b, eps);
  }
  
  // Create extended problem with dummy nodes
  // Add one dummy source and one dummy sink
  arma::vec a_ext(n + 1);
  arma::vec b_ext(m + 1);
  
  // Original marginals with slack
  a_ext.head(n) = a;
  b_ext.head(m) = b;
  
  // Dummy source can supply the deficit (total_b - mass)
  a_ext(n) = std::max(0.0, total_mass_b - mass);
  
  // Dummy sink can absorb the deficit (total_a - mass)
  b_ext(m) = std::max(0.0, total_mass_a - mass);
  
  // Extended cost matrix
  arma::mat cost_ext(n + 1, m + 1);
  
  // Original costs
  cost_ext.submat(0, 0, n-1, m-1) = cost;
  
  // Zero cost for dummy edges (allows free disposal of excess mass)
  cost_ext.col(m).zeros();  // Costs to dummy sink
  cost_ext.row(n).zeros();  // Costs from dummy source
  
  // Solve extended problem
  arma::mat P_ext = network_simplex_ot(cost_ext, a_ext, b_ext, eps);
  
  // Extract the original transport plan (excluding dummy nodes)
  arma::mat P = P_ext.submat(0, 0, n-1, m-1);
  
  // Verify mass constraint
  double actual_mass = arma::accu(P);
  if (std::abs(actual_mass - mass) > eps) {
    // If mass constraint not satisfied exactly, rescale
    P *= (mass / actual_mass);
  }
  
  return P;
}

// Export functions to R
// [[Rcpp::export]]
Rcpp::NumericMatrix network_simplex_ot_cpp(Rcpp::NumericMatrix cost_r,
                                           Rcpp::NumericVector a_r,
                                           Rcpp::NumericVector b_r,
                                           double eps = 1e-9) {
  // Convert to Armadillo
  arma::mat cost = Rcpp::as<arma::mat>(cost_r);
  arma::vec a = Rcpp::as<arma::vec>(a_r);
  arma::vec b = Rcpp::as<arma::vec>(b_r);
  
  // Check inputs
  ::check_marginals(a, b);
  
  // Solve
  arma::mat P = network_simplex_ot(cost, a, b, eps);
  
  // Convert back to R
  return Rcpp::wrap(P);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix fallback_network_simplex_ot_cpp(Rcpp::NumericMatrix cost_r,
                                                    Rcpp::NumericVector a_r,
                                                    Rcpp::NumericVector b_r) {
  arma::mat cost = Rcpp::as<arma::mat>(cost_r);
  arma::vec a = Rcpp::as<arma::vec>(a_r);
  arma::vec b = Rcpp::as<arma::vec>(b_r);

  arma::mat P = fallback::network_simplex_ot(cost, a, b);

  return Rcpp::wrap(P);
}

// Note: sinkhorn_ot_cpp moved to sinkhorn_unified.cpp for consolidation

// [[Rcpp::export]]
Rcpp::NumericMatrix partial_ot_mass_rcpp(Rcpp::NumericMatrix cost_r,
                                         Rcpp::NumericVector a_r,
                                         Rcpp::NumericVector b_r,
                                         double mass,
                                         double eps = 1e-9) {
  // Convert to Armadillo
  arma::mat cost = Rcpp::as<arma::mat>(cost_r);
  arma::vec a = Rcpp::as<arma::vec>(a_r);
  arma::vec b = Rcpp::as<arma::vec>(b_r);
  
  // Solve
  arma::mat P = partial_ot_mass_cpp(cost, a, b, mass, eps);
  
  // Convert back to R
  return Rcpp::wrap(P);
}