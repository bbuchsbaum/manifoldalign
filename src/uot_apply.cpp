#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

inline double neg_inf() {
  return -std::numeric_limits<double>::infinity();
}

inline bool is_neg_inf(double x) {
  return x == neg_inf();
}

inline arma::vec safe_log_weights(const arma::vec& w) {
  arma::vec out(w.n_elem);
  double ninf = neg_inf();
  for (arma::uword i = 0; i < w.n_elem; ++i) {
    out(i) = (w(i) > 0.0) ? std::log(w(i)) : ninf;
  }
  return out;
}

} // namespace

//' Apply UOT coupling as a barycentric map (dense, vector signal)
//'
//' Computes \eqn{\hat s_j = (\sum_i \pi_{ij} s_i) / (\pi_{2,j} + \delta)} without
//' materializing \eqn{\pi}, using translation-invariant potentials.
//'
//' @param cost Dense cost matrix (n x m).
//' @param alpha Source masses (length n).
//' @param beta Target masses (length m).
//' @param fbar Translation-invariant source potential (length n).
//' @param gbar Translation-invariant target potential (length m).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param signal Source signal vector (length n).
//' @param delta Stabilizer added to the denominator.
//' @return A numeric vector (length m) in template space.
//' @export
// [[Rcpp::export]]
arma::vec uot_apply_map_dense_vec_cpp(const arma::mat& cost,
                                     const arma::vec& alpha,
                                     const arma::vec& beta,
                                     const arma::vec& fbar,
                                     const arma::vec& gbar,
                                     double epsilon,
                                     const arma::vec& signal,
                                     double delta = 1e-8) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (delta < 0.0) stop("delta must be >= 0");
  const int n = static_cast<int>(cost.n_rows);
  const int m = static_cast<int>(cost.n_cols);
  if (alpha.n_elem != static_cast<arma::uword>(n)) stop("alpha length must match nrow(cost)");
  if (beta.n_elem != static_cast<arma::uword>(m)) stop("beta length must match ncol(cost)");
  if (fbar.n_elem != static_cast<arma::uword>(n)) stop("fbar length must match nrow(cost)");
  if (gbar.n_elem != static_cast<arma::uword>(m)) stop("gbar length must match ncol(cost)");
  if (signal.n_elem != static_cast<arma::uword>(n)) stop("signal length must match nrow(cost)");

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);
  const double inv_eps = 1.0 / epsilon;

  arma::vec out(m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
    double lb = logbeta(j);
    if (is_neg_inf(lb)) continue;
    double max_logw = neg_inf();
    double sum_scaled = 0.0;
    double num_scaled = 0.0;
    for (int i = 0; i < n; ++i) {
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost(i, j)) * inv_eps;
      if (!std::isfinite(logw)) continue;
      if (logw > max_logw) {
        double scale = std::isfinite(max_logw) ? std::exp(max_logw - logw) : 0.0;
        sum_scaled = sum_scaled * scale + 1.0;
        num_scaled = num_scaled * scale + signal(i);
        max_logw = logw;
      } else {
        double e = std::exp(logw - max_logw);
        sum_scaled += e;
        num_scaled += e * signal(i);
      }
    }
    if (sum_scaled <= 0.0) continue;
    double denom = sum_scaled;
    if (delta > 0.0 && std::isfinite(max_logw)) {
      denom += delta * std::exp(-max_logw);
    }
    out(j) = (denom > 0.0) ? (num_scaled / denom) : 0.0;
  }

  return out;
}

//' Apply UOT coupling as a barycentric map (dense, matrix signal)
//'
//' Matrix version for time series or multi-contrast signals.
//'
//' @param cost Dense cost matrix (n x m).
//' @param alpha Source masses (length n).
//' @param beta Target masses (length m).
//' @param fbar Translation-invariant source potential (length n).
//' @param gbar Translation-invariant target potential (length m).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param signal Source signal matrix (T x n).
//' @param delta Stabilizer added to the denominator.
//' @return A numeric matrix (T x m) in template space.
//' @export
// [[Rcpp::export]]
arma::mat uot_apply_map_dense_mat_cpp(const arma::mat& cost,
                                     const arma::vec& alpha,
                                     const arma::vec& beta,
                                     const arma::vec& fbar,
                                     const arma::vec& gbar,
                                     double epsilon,
                                     const arma::mat& signal,
                                     double delta = 1e-8) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (delta < 0.0) stop("delta must be >= 0");
  const int n = static_cast<int>(cost.n_rows);
  const int m = static_cast<int>(cost.n_cols);
  const int T = static_cast<int>(signal.n_rows);
  if (alpha.n_elem != static_cast<arma::uword>(n)) stop("alpha length must match nrow(cost)");
  if (beta.n_elem != static_cast<arma::uword>(m)) stop("beta length must match ncol(cost)");
  if (fbar.n_elem != static_cast<arma::uword>(n)) stop("fbar length must match nrow(cost)");
  if (gbar.n_elem != static_cast<arma::uword>(m)) stop("gbar length must match ncol(cost)");
  if (signal.n_cols != static_cast<arma::uword>(n)) stop("signal must have n columns");

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);
  const double inv_eps = 1.0 / epsilon;

  arma::mat out(T, m, arma::fill::zeros);

  for (int j = 0; j < m; ++j) {
    double lb = logbeta(j);
    if (is_neg_inf(lb)) continue;
    double max_logw = neg_inf();
    double sum_scaled = 0.0;
    arma::vec num_scaled(T, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost(i, j)) * inv_eps;
      if (!std::isfinite(logw)) continue;
      if (logw > max_logw) {
        double scale = std::isfinite(max_logw) ? std::exp(max_logw - logw) : 0.0;
        sum_scaled = sum_scaled * scale + 1.0;
        num_scaled = num_scaled * scale + signal.col(i);
        max_logw = logw;
      } else {
        double e = std::exp(logw - max_logw);
        sum_scaled += e;
        num_scaled += e * signal.col(i);
      }
    }
    if (sum_scaled <= 0.0) continue;
    double denom = sum_scaled;
    if (delta > 0.0 && std::isfinite(max_logw)) {
      denom += delta * std::exp(-max_logw);
    }
    if (denom > 0.0) {
      out.col(j) = num_scaled / denom;
    }
  }

  return out;
}

//' Apply UOT coupling as a barycentric map (sparse, vector signal)
//'
//' Sparse version using a CSR-like row pointer representation.
//'
//' @param row_ptr Integer vector of length n+1 with CSR row offsets (0-based or
//'   1-based accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).
//' @param col_idx Integer vector of length nnz giving 1-based column indices.
//' @param cost Numeric vector of length nnz aligned with `col_idx`.
//' @param n_rows Number of source nodes.
//' @param n_cols Number of target nodes.
//' @param alpha Source masses (length n_rows).
//' @param beta Target masses (length n_cols).
//' @param fbar Translation-invariant source potential (length n_rows).
//' @param gbar Translation-invariant target potential (length n_cols).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param signal Source signal vector (length n_rows).
//' @param delta Stabilizer added to the denominator.
//' @return A numeric vector (length n_cols) in template space.
//' @export
// [[Rcpp::export]]
arma::vec uot_apply_map_sparse_vec_cpp(const Rcpp::IntegerVector& row_ptr,
                                      const Rcpp::IntegerVector& col_idx,
                                      const Rcpp::NumericVector& cost,
                                      int n_rows,
                                      int n_cols,
                                      const arma::vec& alpha,
                                      const arma::vec& beta,
                                      const arma::vec& fbar,
                                      const arma::vec& gbar,
                                      double epsilon,
                                      const arma::vec& signal,
                                      double delta = 1e-8) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (delta < 0.0) stop("delta must be >= 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (fbar.n_elem != static_cast<arma::uword>(n_rows)) stop("fbar length must match n_rows");
  if (gbar.n_elem != static_cast<arma::uword>(n_cols)) stop("gbar length must match n_cols");
  if (signal.n_elem != static_cast<arma::uword>(n_rows)) stop("signal length must match n_rows");
  if (row_ptr.size() != n_rows + 1) stop("row_ptr must have length n_rows + 1");
  if (col_idx.size() != cost.size()) stop("col_idx and cost must have the same length");

  const int nnz = col_idx.size();

  // Determine whether row_ptr is 0-based or 1-based.
  int rp0 = row_ptr[0];
  int offset = (rp0 == 0) ? 0 : 1;
  if (!(rp0 == 0 || rp0 == 1)) stop("row_ptr[0] must be 0 or 1");

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);
  const double inv_eps = 1.0 / epsilon;

  arma::vec max_logw(n_cols);
  arma::vec sum_scaled(n_cols);
  arma::vec num_scaled(n_cols);
  max_logw.fill(neg_inf());
  sum_scaled.zeros();
  num_scaled.zeros();

  for (int i = 0; i < n_rows; ++i) {
    double la = logalpha(i);
    if (is_neg_inf(la)) continue;
    int start = row_ptr[i] - offset;
    int end   = row_ptr[i + 1] - offset;
    if (start < 0 || end < start || end > nnz) stop("Invalid row_ptr offsets");
    for (int idx = start; idx < end; ++idx) {
      int j = col_idx[idx] - 1;
      if (j < 0 || j >= n_cols) stop("col_idx out of bounds");
      double lb = logbeta(j);
      if (is_neg_inf(lb)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost[idx]) * inv_eps;
      if (!std::isfinite(logw)) continue;

      if (logw > max_logw(j)) {
        double scale = std::isfinite(max_logw(j)) ? std::exp(max_logw(j) - logw) : 0.0;
        sum_scaled(j) = sum_scaled(j) * scale + 1.0;
        num_scaled(j) = num_scaled(j) * scale + signal(i);
        max_logw(j) = logw;
      } else {
        double e = std::exp(logw - max_logw(j));
        sum_scaled(j) += e;
        num_scaled(j) += e * signal(i);
      }
    }
  }

  arma::vec out(n_cols, arma::fill::zeros);
  for (int j = 0; j < n_cols; ++j) {
    if (sum_scaled(j) <= 0.0) continue;
    double denom = sum_scaled(j);
    if (delta > 0.0 && std::isfinite(max_logw(j))) {
      denom += delta * std::exp(-max_logw(j));
    }
    out(j) = (denom > 0.0) ? (num_scaled(j) / denom) : 0.0;
  }

  return out;
}

//' Apply UOT coupling as a barycentric map (sparse, matrix signal via CSC)
//'
//' Sparse matrix version for time series or multi-contrast signals. Requires a
//' CSC representation of the cost graph (column pointers + row indices).
//'
//' This computes \eqn{\hat S = S W / (\pi_2 + \delta)} without materializing
//' the full \eqn{\pi}. The implementation uses per-column log-stabilization to
//' avoid overflow in the exponential weights.
//'
//' @param col_ptr Integer vector of length m+1 with CSC column offsets (0-based
//'   or 1-based accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).
//' @param row_idx Integer vector of length nnz giving 1-based row indices.
//' @param cost_csc Numeric vector of length nnz aligned with `row_idx`.
//' @param n_rows Number of source nodes.
//' @param n_cols Number of target nodes.
//' @param alpha Source masses (length n_rows).
//' @param beta Target masses (length n_cols).
//' @param fbar Translation-invariant source potential (length n_rows).
//' @param gbar Translation-invariant target potential (length n_cols).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param signal Source signal matrix (T x n_rows).
//' @param delta Stabilizer added to the denominator.
//' @return A numeric matrix (T x n_cols) in template space.
//' @export
// [[Rcpp::export]]
arma::mat uot_apply_map_sparse_mat_cpp(const Rcpp::IntegerVector& col_ptr,
                                      const Rcpp::IntegerVector& row_idx,
                                      const Rcpp::NumericVector& cost_csc,
                                      int n_rows,
                                      int n_cols,
                                      const arma::vec& alpha,
                                      const arma::vec& beta,
                                      const arma::vec& fbar,
                                      const arma::vec& gbar,
                                      double epsilon,
                                      const arma::mat& signal,
                                      double delta = 1e-8) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (delta < 0.0) stop("delta must be >= 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (fbar.n_elem != static_cast<arma::uword>(n_rows)) stop("fbar length must match n_rows");
  if (gbar.n_elem != static_cast<arma::uword>(n_cols)) stop("gbar length must match n_cols");
  if (signal.n_cols != static_cast<arma::uword>(n_rows)) stop("signal must have n_rows columns");
  if (col_ptr.size() != n_cols + 1) stop("col_ptr must have length n_cols + 1");
  if (row_idx.size() != cost_csc.size()) stop("row_idx and cost_csc must have the same length");

  const int nnz = row_idx.size();
  const int T = static_cast<int>(signal.n_rows);
  const double inv_eps = 1.0 / epsilon;

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);

  // Determine whether col_ptr is 0-based or 1-based.
  int cp0 = col_ptr[0];
  int offset = (cp0 == 0) ? 0 : 1;
  if (!(cp0 == 0 || cp0 == 1)) stop("col_ptr[0] must be 0 or 1");

  arma::vec max_logw(n_cols);
  arma::vec sum_scaled(n_cols);
  max_logw.fill(neg_inf());
  sum_scaled.zeros();

  // Pass 1: compute per-column (max_logw, sum_scaled) for stable normalization.
  for (int j = 0; j < n_cols; ++j) {
    int start = col_ptr[j] - offset;
    int end   = col_ptr[j + 1] - offset;
    if (start < 0 || end < start || end > nnz) stop("Invalid col_ptr offsets");
    double lb = logbeta(j);
    if (is_neg_inf(lb)) continue;

    double mx = neg_inf();
    double ss = 0.0;
    for (int idx = start; idx < end; ++idx) {
      int i = row_idx[idx] - 1;
      if (i < 0 || i >= n_rows) stop("row_idx out of bounds");
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost_csc[idx]) * inv_eps;
      if (!std::isfinite(logw)) continue;
      if (logw > mx) {
        double scale = std::isfinite(mx) ? std::exp(mx - logw) : 0.0;
        ss = ss * scale + 1.0;
        mx = logw;
      } else {
        ss += std::exp(logw - mx);
      }
    }
    max_logw(j) = mx;
    sum_scaled(j) = ss;
  }

  arma::mat out(T, n_cols, arma::fill::zeros);

  // Pass 2: accumulate numerators in the scaled domain and divide by denom_scaled.
  for (int j = 0; j < n_cols; ++j) {
    double mx = max_logw(j);
    double ss = sum_scaled(j);
    if (!(ss > 0.0) || !std::isfinite(mx)) continue;

    double denom_scaled = ss;
    if (delta > 0.0) {
      denom_scaled += delta * std::exp(-mx);
    }
    if (!(denom_scaled > 0.0)) continue;

    int start = col_ptr[j] - offset;
    int end   = col_ptr[j + 1] - offset;

    double* out_col = out.colptr(j);
    for (int idx = start; idx < end; ++idx) {
      int i = row_idx[idx] - 1;
      if (i < 0 || i >= n_rows) continue;
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      double lb = logbeta(j);
      if (is_neg_inf(lb)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost_csc[idx]) * inv_eps;
      if (!std::isfinite(logw)) continue;
      double e = std::exp(logw - mx);
      if (!(e > 0.0)) continue;

      const double* sig_col = signal.colptr(i);
      for (int t = 0; t < T; ++t) {
        out_col[t] += e * sig_col[t];
      }
    }

    double inv_denom = 1.0 / denom_scaled;
    for (int t = 0; t < T; ++t) {
      out_col[t] *= inv_denom;
    }
  }

  return out;
}
