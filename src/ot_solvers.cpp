#include "ot_solvers.h"
#include "ot_common.h"
#include "network_simplex_adapter.h"
#include <queue>
#include <limits>
#include <cmath>

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

  if (mass < 0.0) {
    Rcpp::stop("partial_ot_mass_cpp: mass must be non-negative");
  }

  double total_mass_a = arma::sum(a);
  double total_mass_b = arma::sum(b);

  double max_mass = std::min(total_mass_a, total_mass_b);
  if (mass > max_mass + eps) {
    Rcpp::stop("partial_ot_mass_cpp: mass exceeds available marginals");
  }

  if (std::abs(mass - max_mass) <= eps) {
    // Degenerates to classical OT
    return network_simplex_ot(cost, a, b, eps);
  }

  double row_slack = std::max(0.0, total_mass_a - mass);
  double col_slack = std::max(0.0, total_mass_b - mass);

  arma::vec a_ext(n + 1);
  arma::vec b_ext(m + 1);
  a_ext.head(n) = a;
  b_ext.head(m) = b;
  a_ext(n) = col_slack;
  b_ext(m) = row_slack;

  arma::mat cost_ext(n + 1, m + 1, arma::fill::zeros);
  cost_ext.submat(0, 0, n - 1, m - 1) = cost;
  cost_ext.submat(0, m, n - 1, m).zeros();
  cost_ext.submat(n, 0, n, m - 1).zeros();

  // Prevent dummy row/column from cancelling each other out
  double max_cost = cost.is_empty() ? 1.0 : cost.max();
  if (!std::isfinite(max_cost)) {
    max_cost = 1.0;
  }
  double penalty = std::max(1.0, std::abs(max_cost)) * 1e6;
  cost_ext(n, m) = penalty;

  arma::mat P_ext = network_simplex_ot(cost_ext, a_ext, b_ext, eps);
  arma::mat P = P_ext.submat(0, 0, n - 1, m - 1);

  double actual_mass = arma::accu(P);
  double mass_error = actual_mass - mass;
  if (std::abs(mass_error) > 10 * eps) {
    // Small numerical drift - correct by trimming the dummy column contribution.
    arma::vec col_sums = arma::sum(P, 0).t();
    arma::vec col_cap = b;
    arma::vec excess = arma::clamp(col_sums - col_cap, 0.0, arma::datum::inf);
    if (arma::accu(excess) > 0) {
      // Remove excess proportionally from affected columns
      for (int j = 0; j < m; ++j) {
        if (excess(j) <= 0) continue;
        double scale = std::max(0.0, (col_sums(j) - excess(j)) / std::max(1e-12, col_sums(j)));
        P.col(j) *= scale;
      }
    } else if (mass_error != 0.0) {
      double scale = mass / std::max(1e-12, actual_mass);
      P *= scale;
    }
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
