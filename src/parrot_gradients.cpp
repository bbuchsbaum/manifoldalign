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

namespace {

// Strictly upper-triangular edge-distance matrix D (squared Euclidean feature
// distance on each undirected edge i < j) and its 0/1 pattern U.
void parrot_edge_mats(const arma::sp_mat& A, const arma::mat& X,
                      arma::sp_mat& D, arma::sp_mat& U) {
  std::vector<arma::uword> rows, cols;
  std::vector<double> dist;
  for (arma::sp_mat::const_iterator it = A.begin(); it != A.end(); ++it) {
    const arma::uword i = it.row();
    const arma::uword j = it.col();
    if (i >= j || (*it) == 0.0) continue;
    double d = 0.0;
    for (arma::uword k = 0; k < X.n_cols; ++k) {
      const double diff = X(i, k) - X(j, k);
      d += diff * diff;
    }
    rows.push_back(i);
    cols.push_back(j);
    dist.push_back(d);
  }
  const arma::uword nnz = rows.size();
  arma::umat loc(2, nnz);
  arma::vec dv(nnz);
  for (arma::uword e = 0; e < nnz; ++e) {
    loc(0, e) = rows[e];
    loc(1, e) = cols[e];
    dv(e) = dist[e];
  }
  D = arma::sp_mat(loc, dv, A.n_rows, A.n_cols);
  U = arma::sp_mat(loc, arma::vec(nnz, arma::fill::ones), A.n_rows, A.n_cols);
}

} // namespace

//' Compute edge consistency gradient
//'
//' Edge-length consistency term over undirected edges (each counted once,
//' \eqn{i < j}, \eqn{a < b}) with squared Euclidean feature distances
//' \eqn{d}: for every edge pair, \eqn{G_{ia} += 2 (d^1_{ij} - d^2_{ab}) S_{jb}}
//' and \eqn{G_{jb} += 2 (d^1_{ij} - d^2_{ab}) S_{ia}}. Instead of looping over
//' all \eqn{|E_1| |E_2|} edge pairs, this is evaluated with sparse products:
//' with \eqn{D_k} the strictly upper-triangular edge-distance matrix and
//' \eqn{U_k} its 0/1 pattern,
//' \eqn{G = 2 (D_1 S U_2^T - U_1 S D_2^T + D_1^T S U_2 - U_1^T S D_2)},
//' at cost \eqn{O((|E_1| + |E_2|) \max(n_1, n_2))}.
//'
//' @param S Current transport plan (n1 x n2)
//' @param A1 Adjacency matrix of network 1 (sparse; nonzero pattern is used)
//' @param A2 Adjacency matrix of network 2 (sparse; nonzero pattern is used)
//' @param X1 Features of network 1 (n1 x d)
//' @param X2 Features of network 2 (n2 x d)
//' @return Gradient matrix (n1 x n2)
// [[Rcpp::export]]
arma::mat compute_edge_gradient_cpp(const arma::mat& S,
                                   const arma::sp_mat& A1,
                                   const arma::sp_mat& A2,
                                   const arma::mat& X1,
                                   const arma::mat& X2) {
  const arma::uword n1 = S.n_rows;
  const arma::uword n2 = S.n_cols;
  if (A1.n_rows != n1 || A1.n_cols != n1 || X1.n_rows != n1) {
    Rcpp::stop("compute_edge_gradient_cpp: A1/X1 do not match nrow(S).");
  }
  if (A2.n_rows != n2 || A2.n_cols != n2 || X2.n_rows != n2) {
    Rcpp::stop("compute_edge_gradient_cpp: A2/X2 do not match ncol(S).");
  }

  arma::sp_mat D1, U1, D2, U2;
  parrot_edge_mats(A1, X1, D1, U1);
  parrot_edge_mats(A2, X2, D2, U2);

  // (sparse * dense) products; S * sp_mat^T computed as (sp_mat * S^T)^T.
  const arma::mat St = S.t();
  const arma::mat SU2t = arma::mat(U2 * St).t();   // S U2^T
  const arma::mat SD2t = arma::mat(D2 * St).t();   // S D2^T
  const arma::mat SU2 = arma::mat(U2.t() * St).t(); // S U2
  const arma::mat SD2 = arma::mat(D2.t() * St).t(); // S D2

  arma::mat grad = arma::mat(D1 * SU2t);
  grad -= arma::mat(U1 * SD2t);
  grad += arma::mat(D1.t() * SU2);
  grad -= arma::mat(U1.t() * SD2);
  return 2.0 * grad;
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