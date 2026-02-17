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

//' Compute log second marginal and transported feature means (dense)
//'
//' Given translation-invariant dual potentials (fbar,gbar), compute the log of
//' the second marginal (column sums) of the entropic UOT coupling and the
//' feature barycentric projection means per target node.
//'
//' @param cost Dense cost matrix (n x m).
//' @param alpha Source masses (length n).
//' @param beta Target masses (length m).
//' @param fbar Translation-invariant source potential (length n).
//' @param gbar Translation-invariant target potential (length m).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param F Optional source features (n x D). If provided, returns mean features.
//' @return A list with `log_pi2` (length m) and optionally `mean_F` (m x D).
//' @export
// [[Rcpp::export]]
Rcpp::List uot_kl_logpi2_meanF_dense_cpp(const arma::mat& cost,
                                        const arma::vec& alpha,
                                        const arma::vec& beta,
                                        const arma::vec& fbar,
                                        const arma::vec& gbar,
                                        double epsilon,
                                        Rcpp::Nullable<arma::mat> F = R_NilValue) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  const int n = static_cast<int>(cost.n_rows);
  const int m = static_cast<int>(cost.n_cols);
  if (alpha.n_elem != static_cast<arma::uword>(n)) stop("alpha length must match nrow(cost)");
  if (beta.n_elem != static_cast<arma::uword>(m)) stop("beta length must match ncol(cost)");
  if (fbar.n_elem != static_cast<arma::uword>(n)) stop("fbar length must match nrow(cost)");
  if (gbar.n_elem != static_cast<arma::uword>(m)) stop("gbar length must match ncol(cost)");

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);
  const double inv_eps = 1.0 / epsilon;

  arma::vec log_pi2(m);
  log_pi2.fill(neg_inf());

  bool have_F = F.isNotNull();
  arma::mat mean_F;
  int D = 0;
  if (have_F) {
    arma::mat Fmat = Rcpp::as<arma::mat>(F);
    if (Fmat.n_rows != static_cast<arma::uword>(n)) stop("F must have n rows");
    D = static_cast<int>(Fmat.n_cols);
    mean_F.set_size(m, D);
    mean_F.zeros();

    for (int j = 0; j < m; ++j) {
      double max_logw = neg_inf();
      double sum_scaled = 0.0;
      arma::rowvec num_scaled(D, arma::fill::zeros);
      double lb = logbeta(j);
      if (is_neg_inf(lb)) {
        log_pi2(j) = neg_inf();
        continue;
      }
      for (int i = 0; i < n; ++i) {
        double la = logalpha(i);
        if (is_neg_inf(la)) continue;
        double logw = la + lb + (fbar(i) + gbar(j) - cost(i, j)) * inv_eps;
        if (!std::isfinite(logw)) continue;
        if (logw > max_logw) {
          double scale = std::isfinite(max_logw) ? std::exp(max_logw - logw) : 0.0;
          sum_scaled = sum_scaled * scale + 1.0;
          num_scaled = num_scaled * scale + Fmat.row(i);
          max_logw = logw;
        } else {
          double e = std::exp(logw - max_logw);
          sum_scaled += e;
          num_scaled += e * Fmat.row(i);
        }
      }
      if (sum_scaled > 0.0) {
        log_pi2(j) = max_logw + std::log(sum_scaled);
        mean_F.row(j) = num_scaled / sum_scaled;
      } else {
        log_pi2(j) = neg_inf();
      }
    }

    return Rcpp::List::create(
      _["log_pi2"] = log_pi2,
      _["mean_F"] = mean_F
    );
  }

  // No features: only log_pi2
  for (int j = 0; j < m; ++j) {
    double max_logw = neg_inf();
    double sum_scaled = 0.0;
    double lb = logbeta(j);
    if (is_neg_inf(lb)) {
      log_pi2(j) = neg_inf();
      continue;
    }
    for (int i = 0; i < n; ++i) {
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      double logw = la + lb + (fbar(i) + gbar(j) - cost(i, j)) * inv_eps;
      if (!std::isfinite(logw)) continue;
      if (logw > max_logw) {
        double scale = std::isfinite(max_logw) ? std::exp(max_logw - logw) : 0.0;
        sum_scaled = sum_scaled * scale + 1.0;
        max_logw = logw;
      } else {
        sum_scaled += std::exp(logw - max_logw);
      }
    }
    if (sum_scaled > 0.0) {
      log_pi2(j) = max_logw + std::log(sum_scaled);
    } else {
      log_pi2(j) = neg_inf();
    }
  }

  return Rcpp::List::create(_["log_pi2"] = log_pi2);
}

//' Compute log second marginal and transported feature means (sparse CSR)
//'
//' Same as `uot_kl_logpi2_meanF_dense_cpp` but for a sparse cost graph.
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
//' @param F Optional source features (n_rows x D).
//' @return A list with `log_pi2` and optionally `mean_F`.
//' @export
// [[Rcpp::export]]
Rcpp::List uot_kl_logpi2_meanF_sparse_cpp(const Rcpp::IntegerVector& row_ptr,
                                         const Rcpp::IntegerVector& col_idx,
                                         const Rcpp::NumericVector& cost,
                                         int n_rows,
                                         int n_cols,
                                         const arma::vec& alpha,
                                         const arma::vec& beta,
                                         const arma::vec& fbar,
                                         const arma::vec& gbar,
                                         double epsilon,
                                         Rcpp::Nullable<arma::mat> F = R_NilValue) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (fbar.n_elem != static_cast<arma::uword>(n_rows)) stop("fbar length must match n_rows");
  if (gbar.n_elem != static_cast<arma::uword>(n_cols)) stop("gbar length must match n_cols");
  if (row_ptr.size() != n_rows + 1) stop("row_ptr must have length n_rows + 1");
  if (col_idx.size() != cost.size()) stop("col_idx and cost must have the same length");

  const int nnz = col_idx.size();
  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);
  const double inv_eps = 1.0 / epsilon;

  // Determine whether row_ptr is 0-based or 1-based.
  int rp0 = row_ptr[0];
  int offset = (rp0 == 0) ? 0 : 1;
  if (!(rp0 == 0 || rp0 == 1)) stop("row_ptr[0] must be 0 or 1");

  arma::vec max_logw(n_cols);
  arma::vec sum_scaled(n_cols);
  max_logw.fill(neg_inf());
  sum_scaled.zeros();

  bool have_F = F.isNotNull();
  arma::mat mean_F;
  arma::mat num_scaled;
  int D = 0;
  arma::mat Fmat;
  if (have_F) {
    Fmat = Rcpp::as<arma::mat>(F);
    if (Fmat.n_rows != static_cast<arma::uword>(n_rows)) stop("F must have n_rows rows");
    D = static_cast<int>(Fmat.n_cols);
    num_scaled.set_size(n_cols, D);
    num_scaled.zeros();
  }

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
        if (have_F) {
          num_scaled.row(j) = num_scaled.row(j) * scale + Fmat.row(i);
        }
        max_logw(j) = logw;
      } else {
        double e = std::exp(logw - max_logw(j));
        sum_scaled(j) += e;
        if (have_F) {
          num_scaled.row(j) += e * Fmat.row(i);
        }
      }
    }
  }

  arma::vec log_pi2(n_cols);
  log_pi2.fill(neg_inf());
  if (have_F) {
    mean_F.set_size(n_cols, D);
    mean_F.zeros();
  }

  for (int j = 0; j < n_cols; ++j) {
    if (sum_scaled(j) > 0.0) {
      log_pi2(j) = max_logw(j) + std::log(sum_scaled(j));
      if (have_F) {
        mean_F.row(j) = num_scaled.row(j) / sum_scaled(j);
      }
    }
  }

  if (have_F) {
    return Rcpp::List::create(
      _["log_pi2"] = log_pi2,
      _["mean_F"] = mean_F
    );
  }

  return Rcpp::List::create(_["log_pi2"] = log_pi2);
}

