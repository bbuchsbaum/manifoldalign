#ifndef OT_SOLVERS_H
#define OT_SOLVERS_H

#include <RcppArmadillo.h>

// Network simplex solver for optimal transport
arma::mat network_simplex_ot(const arma::mat& cost, 
                              const arma::vec& a, 
                              const arma::vec& b,
                              double eps = 1e-9);

// Mass-constrained optimal transport solver
arma::mat partial_ot_mass_cpp(const arma::mat& cost,
                               const arma::vec& a,
                               const arma::vec& b,
                               double mass,
                               double eps = 1e-9);

// Note: Sinkhorn algorithm moved to sinkhorn_unified.cpp

// Note: check_marginals moved to ot_common.h

#endif // OT_SOLVERS_H