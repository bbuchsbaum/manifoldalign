#ifndef NETWORK_SIMPLEX_ADAPTER_H
#define NETWORK_SIMPLEX_ADAPTER_H

#include <RcppArmadillo.h>
#include <vector>
#include <cstdint>

// Adapter around the vendored LEMON / Bonneel network simplex
// (network_simplex_simple.h, full_bipartitegraph.h) for balanced discrete
// optimal transport.
//
// Numerical contract
//  * Costs are passed to the solver as doubles (no integer quantisation); they
//    are shifted and scaled to [0, 1], which leaves the set of optimal plans of
//    a balanced problem unchanged.
//  * Masses are quantised to 64-bit integers with a scale chosen separately
//    from the costs: the common total mass is mapped to `mass_units()` integer
//    units and each marginal is rounded with the largest-remainder method so
//    that both integer marginals sum exactly to that total. No residual is
//    dumped on a single node and no masses are rounded to zero unless they are
//    below one unit (about 1e-15 of the total).
//  * The returned plan is the solver's integer vertex divided by the scale; it
//    is not rebalanced afterwards. It is validated against the original
//    (unquantised) marginals and an error is raised if the deviation exceeds the
//    quantisation bound.
class NetworkSimplexAdapter {
public:
  // Solve min <cost, P> s.t. P 1 = a, P' 1 = b, P >= 0.
  // `eps` is the tolerance used when checking that sum(a) == sum(b).
  static arma::mat solve(const arma::mat& cost,
                         const arma::vec& a,
                         const arma::vec& b,
                         double eps = 1e-9);

  // Check if LEMON library is available at compile time
  static bool isAvailable();

  // Number of integer units the total mass is mapped to.
  static double mass_units() { return 1125899906842624.0; }  // 2^50

  // Largest-remainder rounding of `w * scale` to non-negative integers that
  // sum exactly to `total`. Exposed for testing.
  static std::vector<int64_t> quantize_masses(const arma::vec& w,
                                              double scale,
                                              int64_t total);
};

// Fallback implementation using lpSolve via the R helper classical_ot_lp_r()
// (used if LEMON is not available).
namespace fallback {
  arma::mat network_simplex_ot(const arma::mat& cost,
                               const arma::vec& a,
                               const arma::vec& b,
                               double eps = 1e-9);
}

#endif // NETWORK_SIMPLEX_ADAPTER_H
