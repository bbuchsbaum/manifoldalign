#ifndef TV_PROXIMAL_H
#define TV_PROXIMAL_H

#include <RcppArmadillo.h>

// Total Variation proximal operator using Pool Adjacent Violators Algorithm (PAVA)
// Solves: argmin_x { 0.5 * ||x - y||^2 + lambda * TV(x) }
// where TV(x) = sum_i |x[i+1] - x[i]|
class TVProximal {
public:
  // Apply 1D TV proximal operator to a vector
  // Uses isotonic regression (PAVA) algorithm
  static arma::vec prox_tv_1d(const arma::vec& y, double lambda);
  
  // Apply TV proximal operator to rows of a matrix
  static arma::mat prox_tv_rows(const arma::mat& Y, double lambda);
  
  // Apply TV proximal operator to columns of a matrix
  static arma::mat prox_tv_cols(const arma::mat& Y, double lambda);
  
  // Apply 2D TV proximal operator (alternating rows and cols)
  static arma::mat prox_tv_2d(const arma::mat& Y, double lambda, int max_iter = 10);
  
private:
  // Helper: Pool Adjacent Violators Algorithm
  static void pava_increasing(arma::vec& y, const arma::vec& w);
  
  // Helper: Compute subgradient for TV
  static arma::vec tv_subgradient(const arma::vec& x);
};

#endif // TV_PROXIMAL_H