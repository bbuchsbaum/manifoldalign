#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Compute squared Euclidean distance matrix
//' 
//' Efficient computation of pairwise squared distances
//' 
//' @param X1 First matrix (n1 x d)
//' @param X2 Second matrix (n2 x d)
//' @return Distance matrix (n1 x n2)
// [[Rcpp::export]]
arma::mat compute_squared_distances_cpp(const arma::mat& X1, const arma::mat& X2) {
  // Compute squared norms
  arma::vec X1_sq = arma::sum(X1 % X1, 1);
  arma::vec X2_sq = arma::sum(X2 % X2, 1);
  
  // Compute distance matrix using broadcasting
  arma::mat D = arma::repmat(X1_sq, 1, X2.n_rows) + 
                arma::repmat(X2_sq.t(), X1.n_rows, 1) - 
                2 * X1 * X2.t();
  
  // Ensure non-negative (numerical issues can cause small negative values)
  D.transform([](double val) { return std::max(0.0, val); });
  
  return D;
}

//' Solve Sylvester equation for RWR cost
//' 
//' Iterative solution of: C_rwr = (1+β)C_node + (1-β)γ * W1 * C_rwr * W2^T
//' 
//' @param W1 Transition matrix of network 1
//' @param W2T Transpose of transition matrix of network 2
//' @param C_node Node-level cost matrix
//' @param beta RWR restart probability
//' @param gamma Cross-graph discount factor
//' @param tol Convergence tolerance
//' @param max_iter Maximum iterations
//' @return Solution matrix C_rwr
// [[Rcpp::export]]
arma::mat solve_sylvester_rwr_cpp(const arma::mat& W1,
                                  const arma::mat& W2T,
                                  const arma::mat& C_node,
                                  double beta = 0.15,
                                  double gamma = 0.1,
                                  double tol = 1e-6,
                                  int max_iter = 50) {
  arma::mat X = C_node;
  arma::mat X_old;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    X_old = X;
    X = (1 + beta) * C_node + (1 - beta) * gamma * (W1 * X * W2T);
    
    // Check convergence
    double max_diff = arma::abs(X - X_old).max();
    if (max_diff < tol) {
      break;
    }
  }
  
  return X;
}

//' Vectorized RWR computation
//' 
//' Compute Random Walk with Restart for multiple restart vectors simultaneously
//' 
//' @param W_transpose Transpose of transition matrix
//' @param E Matrix of restart vectors (columns are different restart distributions)
//' @param sigma Restart probability
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @return RWR result matrix
// [[Rcpp::export]]
arma::mat compute_rwr_vectorized_cpp(const arma::mat& W_transpose,
                                     const arma::mat& E,
                                     double sigma,
                                     int max_iter,
                                     double tol) {
  int n = W_transpose.n_cols;
  int n_anchors = E.n_cols;
  
  // Initialize with uniform distribution
  arma::mat R = E / n;
  arma::mat R_old;
  
  for (int iter = 0; iter < max_iter; ++iter) {
    R_old = R;
    R = (1 - sigma) * (W_transpose * R) + sigma * E;
    
    // Check convergence
    double max_diff = arma::abs(R - R_old).max();
    if (max_diff < tol) {
      break;
    }
  }
  
  return R;
}

//' Compute PARROT cost matrix
//' 
//' Complete cost matrix computation including Sylvester equation
//' 
//' @param X1 Features of network 1
//' @param X2 Features of network 2
//' @param R1 RWR descriptors of network 1
//' @param R2 RWR descriptors of network 2
//' @param W1 Transition matrix of network 1
//' @param W2 Transition matrix of network 2
//' @param alpha Weight for attribute vs RWR cost
//' @param sigma RWR restart probability
//' @param gamma Cross-graph discount factor
//' @return Position-aware cost matrix
// [[Rcpp::export]]
arma::mat compute_parrot_cost_cpp(const arma::mat& X1,
                                  const arma::mat& X2,
                                  const arma::mat& R1,
                                  const arma::mat& R2,
                                  const arma::mat& W1,
                                  const arma::mat& W2,
                                  double alpha = 0.5,
                                  double sigma = 0.15,
                                  double gamma = 0.1) {
  auto row_normalize = [](const arma::mat& M) {
    arma::mat out = M;
    arma::vec norms = arma::sqrt(arma::sum(arma::square(M), 1));
    for (arma::uword i = 0; i < M.n_rows; ++i) {
      double denom = norms(i);
      if (denom < 1e-12) denom = 1e-12;
      out.row(i) /= denom;
    }
    return out;
  };

  arma::mat R1n = row_normalize(R1);
  arma::mat R2n = row_normalize(R2);
  arma::mat X1n = row_normalize(X1);
  arma::mat X2n = row_normalize(X2);

  auto clamp_cosine = [](arma::mat& M) {
    M.transform([](double val) {
      if (val > 1.0) return 1.0;
      if (val < -1.0) return -1.0;
      return val;
    });
  };

  arma::mat sim_rwr = R1n * R2n.t();
  arma::mat sim_attr = X1n * X2n.t();
  clamp_cosine(sim_rwr);
  clamp_cosine(sim_attr);

  double kappa = 1.0;
  arma::mat cost_rwr = arma::exp(-kappa * sim_rwr);
  arma::mat cost_attr = arma::exp(-kappa * sim_attr);

  arma::mat C_node = (1 - alpha) * cost_rwr + alpha * cost_attr;

  arma::mat W2T = W2.t();
  arma::mat C_rwr = solve_sylvester_rwr_cpp(W1, W2T, C_node, sigma, gamma);

  double min_c_rwr = C_rwr.min();
  if (min_c_rwr <= 0) {
    C_rwr = C_rwr - min_c_rwr + 1e-6;
  }

  return C_rwr;
}