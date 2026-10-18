#include "ot_solvers.h"
#include "ot_common.h"
#include "network_simplex_adapter.h"
#include <limits>
#include <cmath>
#include <algorithm>

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

// Mass-constrained partial OT:
//   min <cost, P>  s.t.  P >= 0, P 1 <= a, P' 1 <= b, sum(P) = mass.
// Solved exactly as a balanced problem on an extended graph: a dummy row with
// mass sum(b) - mass and a dummy column with mass sum(a) - mass absorb the
// untransported mass at zero cost. The dummy-to-dummy arc gets a cost M with
// M + c_ij > 0 for every real arc, so any flow through it can be strictly
// improved and every optimum transports exactly `mass` between real nodes.
arma::mat partial_ot_mass_cpp(const arma::mat& cost,
                               const arma::vec& a,
                               const arma::vec& b,
                               double mass,
                               double eps) {

  const int n = a.n_elem;
  const int m = b.n_elem;

  if (cost.n_rows != a.n_elem || cost.n_cols != b.n_elem) {
    Rcpp::stop("partial_ot_mass_cpp: cost must be %d x %d", n, m);
  }
  if (!cost.is_finite() || !a.is_finite() || !b.is_finite()) {
    Rcpp::stop("partial_ot_mass_cpp: cost and marginals must be finite");
  }
  if (arma::any(a < 0) || arma::any(b < 0)) {
    Rcpp::stop("partial_ot_mass_cpp: marginals must be non-negative");
  }
  if (!std::isfinite(mass) || mass < 0.0) {
    Rcpp::stop("partial_ot_mass_cpp: mass must be finite and non-negative");
  }

  const double total_mass_a = arma::sum(a);
  const double total_mass_b = arma::sum(b);
  const double max_mass = std::min(total_mass_a, total_mass_b);
  if (mass > max_mass + eps) {
    Rcpp::stop("partial_ot_mass_cpp: mass (%g) exceeds min(sum(a), sum(b)) = %g",
               mass, max_mass);
  }
  mass = std::min(mass, max_mass);
  if (n == 0 || m == 0 || mass == 0.0) {
    return arma::mat(n, m, arma::fill::zeros);
  }

  arma::vec a_ext(n + 1);
  arma::vec b_ext(m + 1);
  a_ext.head(n) = a;
  b_ext.head(m) = b;
  a_ext(n) = std::max(0.0, total_mass_b - mass);  // dummy source
  b_ext(m) = std::max(0.0, total_mass_a - mass);  // dummy sink

  const double cmin = cost.min();
  const double cmax = cost.max();
  arma::mat cost_ext(n + 1, m + 1, arma::fill::zeros);
  cost_ext.submat(0, 0, n - 1, m - 1) = cost;
  cost_ext(n, m) = std::abs(cmin) + std::abs(cmax) + 1.0;

  arma::mat P_ext = network_simplex_ot(cost_ext, a_ext, b_ext, eps);
  arma::mat P = P_ext.submat(0, 0, n - 1, m - 1);

  // Validate against the original constraints. The plan is never rescaled:
  // rescaling would leave the LP vertex and can break the row/column bounds.
  const double total_ext = arma::sum(a_ext);
  const double tol = 1e-10 * std::max(1.0, total_ext);
  const double mass_err = std::abs(arma::accu(P) - mass);
  const double row_excess = arma::max(arma::sum(P, 1) - a);
  const double col_excess = arma::max(arma::sum(P, 0).t() - b);
  if (mass_err > tol || row_excess > tol || col_excess > tol || P.min() < 0.0) {
    Rcpp::stop("partial_ot_mass_cpp: solver returned an infeasible plan "
               "(mass error %g, row excess %g, column excess %g)",
               mass_err, row_excess, col_excess);
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
