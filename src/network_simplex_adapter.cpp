#include "network_simplex_adapter.h"
#include "ot_common.h"
#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>

// Check if we have the LEMON headers available
#ifdef __has_include
  #if __has_include("network_simplex_simple.h") && __has_include("full_bipartitegraph.h")
    #define LEMON_AVAILABLE 1
    #include "network_simplex_simple.h"
    #include "full_bipartitegraph.h"
  #else
    #define LEMON_AVAILABLE 0
  #endif
#else
  // Try to include anyway for older compilers
  #define LEMON_AVAILABLE 1
  #include "network_simplex_simple.h"
  #include "full_bipartitegraph.h"
#endif

std::vector<int64_t> NetworkSimplexAdapter::quantize_masses(const arma::vec& w,
                                                            double scale,
                                                            int64_t total) {
  const arma::uword n = w.n_elem;
  std::vector<int64_t> out(n, 0);
  std::vector<double> frac(n, -1.0);
  int64_t assigned = 0;
  arma::uword n_pos = 0;
  for (arma::uword i = 0; i < n; ++i) {
    if (w(i) <= 0.0) continue;  // zero masses stay exactly zero
    ++n_pos;
    double x = w(i) * scale;
    double fl = std::floor(x);
    out[i] = static_cast<int64_t>(fl);
    frac[i] = x - fl;
    assigned += out[i];
  }
  int64_t rem = total - assigned;
  if (rem == 0 || n_pos == 0) {
    return out;
  }

  // Largest remainder: hand out the missing units to the largest fractional
  // parts (or take surplus units from the smallest ones). Only positive masses
  // take part, so a zero mass never receives flow.
  std::vector<arma::uword> order;
  order.reserve(n_pos);
  for (arma::uword i = 0; i < n; ++i) {
    if (frac[i] >= 0.0) order.push_back(i);
  }
  if (rem > 0) {
    std::stable_sort(order.begin(), order.end(), [&](arma::uword x, arma::uword y) {
      return frac[x] > frac[y];
    });
    while (rem > 0) {
      for (std::size_t k = 0; k < order.size() && rem > 0; ++k) {
        out[order[k]] += 1;
        --rem;
      }
    }
  } else {
    std::stable_sort(order.begin(), order.end(), [&](arma::uword x, arma::uword y) {
      return frac[x] < frac[y];
    });
    bool progressed = true;
    while (rem < 0 && progressed) {
      progressed = false;
      for (std::size_t k = 0; k < order.size() && rem < 0; ++k) {
        if (out[order[k]] > 0) {
          out[order[k]] -= 1;
          ++rem;
          progressed = true;
        }
      }
    }
    if (rem != 0) {
      Rcpp::stop("Network simplex: internal error while quantising masses.");
    }
  }
  return out;
}

namespace {

void validate_inputs(const arma::mat& cost, const arma::vec& a, const arma::vec& b) {
  if (cost.n_rows != a.n_elem || cost.n_cols != b.n_elem) {
    Rcpp::stop("Network simplex: cost must be %d x %d to match the marginals.",
               static_cast<int>(a.n_elem), static_cast<int>(b.n_elem));
  }
  if (!a.is_finite() || !b.is_finite()) {
    Rcpp::stop("Network simplex: marginals must be finite.");
  }
  if (!cost.is_finite()) {
    Rcpp::stop("Network simplex: cost matrix must be finite.");
  }
}

// Check the returned plan against the original (unquantised) marginals.
void validate_plan(const arma::mat& P, const arma::vec& a, const arma::vec& b,
                   double tol) {
  arma::vec rs = arma::sum(P, 1);
  arma::vec cs = arma::sum(P, 0).t();
  double row_err = a.n_elem ? arma::max(arma::abs(rs - a)) : 0.0;
  double col_err = b.n_elem ? arma::max(arma::abs(cs - b)) : 0.0;
  if (!P.is_finite() || P.min() < 0.0 || row_err > tol || col_err > tol) {
    Rcpp::stop("Network simplex returned a plan that violates the marginals "
               "(row error %g, column error %g, tolerance %g).",
               row_err, col_err, tol);
  }
}

}  // namespace

// Implementation of the adapter
arma::mat NetworkSimplexAdapter::solve(const arma::mat& cost,
                                       const arma::vec& a,
                                       const arma::vec& b,
                                       double eps) {
  validate_inputs(cost, a, b);
  check_marginals(a, b, eps);

  const int n = static_cast<int>(a.n_elem);
  const int m = static_cast<int>(b.n_elem);
  if (n == 0 || m == 0) {
    return arma::mat(n, m, arma::fill::zeros);
  }

  const double sum_a = arma::sum(a);
  const double sum_b = arma::sum(b);
  const double total = 0.5 * (sum_a + sum_b);
  if (!(total > 0.0)) {
    return arma::mat(n, m, arma::fill::zeros);
  }

#if LEMON_AVAILABLE
  typedef lemon::FullBipartiteDigraph Digraph;
  typedef int64_t FlowType;   // exact integer flows
  typedef double CostType;    // costs are not quantised
  typedef lemon::NetworkSimplexSimple<Digraph, FlowType, CostType, int64_t> NetworkSimplex;

  // --- Masses: integer units, scaled independently of the costs ------------
  const double units = mass_units();
  const int64_t total_units = static_cast<int64_t>(units);
  std::vector<int64_t> sup_a = quantize_masses(a, units / sum_a, total_units);
  std::vector<int64_t> sup_b = quantize_masses(b, units / sum_b, total_units);
  std::vector<FlowType> supplies_src(n);
  std::vector<FlowType> supplies_dst(m);
  for (int i = 0; i < n; ++i) supplies_src[i] = sup_a[i];
  for (int j = 0; j < m; ++j) supplies_dst[j] = -sup_b[j];

  // --- Costs: shift/scale to [0, 1] (does not change the optimal plans) ----
  const double cmin = cost.min();
  const double cmax = cost.max();
  const double crange = cmax - cmin;
  const double cscale = (crange > 0.0) ? 1.0 / crange : 0.0;

  Digraph graph(n, m);
  const int64_t n_arcs = static_cast<int64_t>(n) * static_cast<int64_t>(m);
  NetworkSimplex ns(graph, true, n + m, n_arcs);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      int64_t arc_id = static_cast<int64_t>(i) * m + j;
      ns.setCost(graph.arcFromId(arc_id), (cost(i, j) - cmin) * cscale);
    }
  }
  ns.supplyMap(&supplies_src[0], n, &supplies_dst[0], m);

  int result = ns.run();
  if (result != NetworkSimplex::OPTIMAL) {
    if (result == NetworkSimplex::INFEASIBLE) {
      Rcpp::stop("Network simplex: problem is infeasible (%d x %d).", n, m);
    } else if (result == NetworkSimplex::UNBOUNDED) {
      Rcpp::stop("Network simplex: problem is unbounded (%d x %d).", n, m);
    }
    Rcpp::stop("Network simplex: error code %d for problem %d x %d.", result, n, m);
  }

  // Integer vertex -> plan. Rows are divided by the row scale, which keeps
  // each row sum within one unit of the original marginal.
  const double inv_scale = total / units;
  arma::mat P(n, m, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      int64_t arc_id = static_cast<int64_t>(i) * m + j;
      FlowType flow = ns.flow(graph.arcFromId(arc_id));
      if (flow != 0) {
        P(i, j) = static_cast<double>(flow) * inv_scale;
      }
    }
  }

  // Quantisation bound: one unit per node plus the input mass imbalance and
  // floating-point summation error.
  const double tol = std::abs(sum_a - sum_b) +
    total * (4.0 / units + 4.0 * (n + m) * std::numeric_limits<double>::epsilon());
  validate_plan(P, a, b, tol);
  return P;

#else
  Rcpp::warning("LEMON network simplex not available, using lpSolve fallback");
  arma::mat P = fallback::network_simplex_ot(cost, a, b, eps);
  validate_plan(P, a, b, std::abs(sum_a - sum_b) + 1e-9 * total);
  return P;
#endif
}

bool NetworkSimplexAdapter::isAvailable() {
#if LEMON_AVAILABLE
  return true;
#else
  return false;
#endif
}

// Fallback implementation using lpSolve via the R helper
namespace fallback {

arma::mat network_simplex_ot(const arma::mat& cost,
                              const arma::vec& a,
                              const arma::vec& b,
                              double /*eps*/) {
  Rcpp::Environment pkg = Rcpp::Environment::namespace_env("manifoldalign");
  Rcpp::Function lp_fun = pkg["classical_ot_lp_r"];

  Rcpp::NumericMatrix cost_r = Rcpp::wrap(cost);
  Rcpp::NumericVector a_r = Rcpp::wrap(a);
  Rcpp::NumericVector b_r = Rcpp::wrap(b);

  Rcpp::NumericMatrix sol = lp_fun(cost_r, a_r, b_r);

  return Rcpp::as<arma::mat>(sol);
}

} // namespace fallback
