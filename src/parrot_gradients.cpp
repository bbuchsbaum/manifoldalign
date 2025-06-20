#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Compute pairwise squared distances for connected nodes
//' 
//' @param X Feature matrix (n x d)
//' @param A Adjacency matrix (sparse, n x n)
//' @return Matrix of squared distances for edges
// [[Rcpp::export]]
arma::mat compute_edge_distances_cpp(const arma::mat& X, const arma::sp_mat& A) {
  int n = X.n_rows;
  arma::mat D(n, n, arma::fill::zeros);
  
  // Iterate over non-zero elements of adjacency matrix
  arma::sp_mat::const_iterator it = A.begin();
  arma::sp_mat::const_iterator it_end = A.end();
  
  for (; it != it_end; ++it) {
    int i = it.row();
    int j = it.col();
    
    // Compute squared Euclidean distance
    arma::rowvec diff = X.row(i) - X.row(j);
    D(i, j) = arma::dot(diff, diff);
  }
  
  return D;
}

//' Compute edge consistency gradient
//' 
//' Efficient computation of the gradient of edge consistency regularizer
//' 
//' @param S Current transport plan (n1 x n2)
//' @param A1 Adjacency matrix of network 1 (sparse)
//' @param A2 Adjacency matrix of network 2 (sparse)
//' @param X1 Features of network 1 (n1 x d)
//' @param X2 Features of network 2 (n2 x d)
//' @return Gradient matrix (n1 x n2)
// [[Rcpp::export]]
arma::mat compute_edge_gradient_cpp(const arma::mat& S,
                                   const arma::sp_mat& A1,
                                   const arma::sp_mat& A2,
                                   const arma::mat& X1,
                                   const arma::mat& X2) {
  int n1 = S.n_rows;
  int n2 = S.n_cols;
  arma::mat grad(n1, n2, arma::fill::zeros);

  arma::sp_mat::const_iterator it1_end = A1.end();
  arma::sp_mat::const_iterator it2_end = A2.end();

  for (arma::sp_mat::const_iterator it1 = A1.begin(); it1 != it1_end; ++it1) {
    int i = it1.row();
    int j = it1.col();
    
    // To match R, only process lower-triangular part to avoid double counting
    if (i >= j) continue;

    arma::rowvec diff1 = X1.row(i) - X1.row(j);
    double d1_sq = arma::dot(diff1, diff1);

    for (arma::sp_mat::const_iterator it2 = A2.begin(); it2 != it2_end; ++it2) {
      int a = it2.row();
      int b = it2.col();
      
      if (a >= b) continue;

      arma::rowvec diff2 = X2.row(a) - X2.row(b);
      double d2_sq = arma::dot(diff2, diff2);
      
      double term = 2 * (d1_sq - d2_sq);
      
      grad(i, a) += S(j, b) * term;
      grad(j, b) += S(i, a) * term;
    }
  }

  return grad;
}

//' Compute neighborhood consistency gradient
//' 
//' Gradient of KL divergence term for neighborhood consistency
//' 
//' @param S_prev Previous transport plan
//' @param W1 Transition matrix of network 1
//' @param W2 Transition matrix of network 2
//' @param eps Small epsilon for numerical stability
//' @return Gradient matrix
// [[Rcpp::export]]
arma::mat compute_neighborhood_gradient_cpp(const arma::mat& S_prev,
                                           const arma::mat& W1,
                                           const arma::mat& W2,
                                           double eps = 1e-16) {
  // Compute S_hat = W1^T * S * W2
  arma::mat S_hat = W1.t() * S_prev * W2;
  
  // Compute ratio with numerical stability
  arma::mat S_safe = S_prev + eps;
  arma::mat S_hat_safe = S_hat + eps;
  arma::mat ratio = S_hat_safe / S_safe;
  
  // Gradient: -W1 * (ratio - 1) * W2^T
  arma::mat grad = -W1 * (ratio - arma::ones(ratio.n_rows, ratio.n_cols)) * W2.t();
  
  return grad;
}

//' Compute anchor prior gradient
//' 
//' Gradient of KL divergence from anchor prior
//' 
//' @param S_prev Previous transport plan
//' @param anchor_idx1 Indices of anchors in network 1
//' @param anchor_idx2 Indices of anchors in network 2
//' @param anchor_vals1 Anchor values for network 1
//' @param anchor_vals2 Anchor values for network 2
//' @param eps Small epsilon for numerical stability
//' @return Gradient matrix
// [[Rcpp::export]]
arma::mat compute_anchor_gradient_cpp(const arma::mat& S_prev,
                                     const arma::uvec& anchor_idx1,
                                     const arma::uvec& anchor_idx2,
                                     const arma::vec& anchor_vals1,
                                     const arma::vec& anchor_vals2,
                                     double eps = 1e-16) {
  int n1 = S_prev.n_rows;
  int n2 = S_prev.n_cols;
  
  // Initialize H matrix with small values
  arma::mat H(n1, n2, arma::fill::value(eps));
  
  // Set high probability for anchor correspondences
  for (size_t i = 0; i < anchor_idx1.n_elem; ++i) {
    int idx1 = anchor_idx1(i) - 1; // Convert to 0-based indexing
    double anchor_val = anchor_vals1(i);
    
    for (size_t j = 0; j < anchor_idx2.n_elem; ++j) {
      int idx2 = anchor_idx2(j) - 1;
      if (std::abs(anchor_vals2(j) - anchor_val) < eps) {
        H(idx1, idx2) = 1.0;
      }
    }
  }
  
  // Normalize H to be a valid probability distribution
  H = H / arma::accu(H);
  
  // Gradient of KL(S || H) is log(S/H)
  arma::mat S_safe = S_prev + eps;
  arma::mat grad = arma::log(S_safe / H);
  
  return grad;
}