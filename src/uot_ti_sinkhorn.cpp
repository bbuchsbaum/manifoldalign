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

inline double logsumexp_affine(const arma::vec& logw,
                               const arma::vec& x,
                               double a) {
  // returns log sum_i exp(logw[i] + a * x[i]), treating logw=-inf as zero weight.
  double max_val = neg_inf();
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double lw = logw(i);
    if (is_neg_inf(lw)) continue;
    double v = lw + a * x(i);
    if (v > max_val) max_val = v;
  }
  if (is_neg_inf(max_val)) return neg_inf();
  double sum = 0.0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    double lw = logw(i);
    if (is_neg_inf(lw)) continue;
    double v = lw + a * x(i);
    sum += std::exp(v - max_val);
  }
  return max_val + std::log(sum);
}

inline double softmin_scalar(const arma::vec& logw,
                             const arma::vec& v,
                             double temp) {
  // -temp * log sum_i w_i * exp(-v_i/temp)
  double inv = -1.0 / temp;
  double lse = logsumexp_affine(logw, v, inv);
  if (is_neg_inf(lse)) return std::numeric_limits<double>::infinity();
  return -temp * lse;
}

inline double max_abs_diff(const arma::vec& a, const arma::vec& b) {
  if (a.n_elem != b.n_elem) {
    stop("Internal error: vector size mismatch");
  }
  double mx = 0.0;
  for (arma::uword i = 0; i < a.n_elem; ++i) {
    double d = std::abs(a(i) - b(i));
    if (d > mx) mx = d;
  }
  return mx;
}

} // namespace

//' KL translation-invariant Sinkhorn (dense cost)
//'
//' CPU implementation of the KL translation-invariant Sinkhorn updates for
//' entropic unbalanced OT with KL marginal penalties.
//'
//' @param cost Dense cost matrix (n x m).
//' @param alpha Source masses (length n, nonnegative, not all zero).
//' @param beta Target masses (length m, nonnegative, not all zero).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param rho1 KL penalty on the first marginal (> 0).
//' @param rho2 KL penalty on the second marginal (> 0).
//' @param max_iter Maximum number of iterations.
//' @param tol Stopping tolerance on the infinity-norm iterate difference.
//' @return A list with translation-invariant potentials `fbar`, `gbar`,
//'   translated dual potentials `f`, `g`, translation `lambda`, iteration
//'   count, convergence flag, and last residual.
//' @export
// [[Rcpp::export]]
Rcpp::List uot_ti_sinkhorn_kl_dense_cpp(const arma::mat& cost,
                                       const arma::vec& alpha,
                                       const arma::vec& beta,
                                       double epsilon,
                                       double rho1,
                                       double rho2,
                                       int max_iter = 2000,
                                       double tol = 1e-6) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (rho1 <= 0.0) stop("rho1 must be > 0");
  if (rho2 <= 0.0) stop("rho2 must be > 0");
  if (alpha.n_elem != cost.n_rows) stop("alpha length must match nrow(cost)");
  if (beta.n_elem != cost.n_cols) stop("beta length must match ncol(cost)");

  const int n = static_cast<int>(cost.n_rows);
  const int m = static_cast<int>(cost.n_cols);

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);

  arma::vec fbar(n, arma::fill::zeros);
  arma::vec gbar(m, arma::fill::zeros);
  arma::vec ftmp(n);
  arma::vec gtmp(m);
  arma::vec fbar_new(n);
  arma::vec gbar_new(m);

  const double denom = epsilon + rho1 + rho2;
  const double xi12 = epsilon * rho2 / (rho1 * denom);
  const double xi21 = epsilon * rho1 / (rho2 * denom);
  const double k11 = (epsilon / (epsilon + rho1)) * (rho1 / (rho1 + rho2));
  const double k22 = (epsilon / (epsilon + rho2)) * (rho2 / (rho1 + rho2));
  const double a1 = rho1 / (rho1 + epsilon);
  const double a2 = rho2 / (rho2 + epsilon);
  const double inv_eps = 1.0 / epsilon;

  bool converged = false;
  double residual = NA_REAL;
  int it = 0;

  for (it = 0; it < max_iter; ++it) {
    // scalar softmin_beta_rho2(gbar)
    double smin_beta_rho2_gbar = softmin_scalar(logbeta, gbar, rho2);

    // row-wise softmin_beta_eps(C - gbar)
    for (int i = 0; i < n; ++i) {
      double max_z = neg_inf();
      double sum_exp = 0.0;
      for (int j = 0; j < m; ++j) {
        double lb = logbeta(j);
        if (is_neg_inf(lb)) continue;
        double z = lb + (gbar(j) - cost(i, j)) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z) {
          double scale = std::isfinite(max_z) ? std::exp(max_z - z) : 0.0;
          sum_exp = sum_exp * scale + 1.0;
          max_z = z;
        } else {
          sum_exp += std::exp(z - max_z);
        }
      }
      double logsum = (sum_exp > 0.0) ? (max_z + std::log(sum_exp)) : neg_inf();
      double smin_eps_row = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      ftmp(i) = a1 * smin_eps_row - k11 * smin_beta_rho2_gbar;
    }

    double smin_alpha_rho1_ftmp = softmin_scalar(logalpha, ftmp, rho1);
    fbar_new = ftmp + xi12 * smin_alpha_rho1_ftmp;

    // scalar softmin_alpha_rho1(fbar_new)
    double smin_alpha_rho1_fbar = softmin_scalar(logalpha, fbar_new, rho1);

    // col-wise softmin_alpha_eps(C^T - fbar_new)
    for (int j = 0; j < m; ++j) {
      double max_z = neg_inf();
      double sum_exp = 0.0;
      for (int i = 0; i < n; ++i) {
        double la = logalpha(i);
        if (is_neg_inf(la)) continue;
        double z = la + (fbar_new(i) - cost(i, j)) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z) {
          double scale = std::isfinite(max_z) ? std::exp(max_z - z) : 0.0;
          sum_exp = sum_exp * scale + 1.0;
          max_z = z;
        } else {
          sum_exp += std::exp(z - max_z);
        }
      }
      double logsum = (sum_exp > 0.0) ? (max_z + std::log(sum_exp)) : neg_inf();
      double smin_eps_col = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      gtmp(j) = a2 * smin_eps_col - k22 * smin_alpha_rho1_fbar;
    }

    double smin_beta_rho2_gtmp = softmin_scalar(logbeta, gtmp, rho2);
    gbar_new = gtmp + xi21 * smin_beta_rho2_gtmp;

    residual = std::max(max_abs_diff(fbar_new, fbar), max_abs_diff(gbar_new, gbar));
    fbar.swap(fbar_new);
    gbar.swap(gbar_new);

    if (residual <= tol) {
      converged = true;
      break;
    }
  }

  // translation parameter lambda*
  double La = logsumexp_affine(logalpha, fbar, -1.0 / rho1);
  double Lb = logsumexp_affine(logbeta, gbar, -1.0 / rho2);
  double lambda = (rho1 * rho2) / (rho1 + rho2) * (La - Lb);

  arma::vec f = fbar + lambda;
  arma::vec g = gbar - lambda;

  return Rcpp::List::create(
    _["fbar"] = fbar,
    _["gbar"] = gbar,
    _["lambda"] = lambda,
    _["f"] = f,
    _["g"] = g,
    _["iterations"] = it + 1,
    _["converged"] = converged,
    _["residual"] = residual
  );
}

//' KL translation-invariant Sinkhorn (sparse kNN cost)
//'
//' Sparse CPU implementation using a CSR-like row pointer representation.
//'
//' @param row_ptr Integer vector of length n+1 with CSR row offsets (0-based or
//'   1-based accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).
//' @param col_idx Integer vector of length nnz giving 1-based column indices.
//' @param cost Numeric vector of length nnz giving cost values aligned with
//'   `col_idx`.
//' @param n_rows Number of source nodes (n).
//' @param n_cols Number of target nodes (m).
//' @param alpha Source masses (length n).
//' @param beta Target masses (length m).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param rho1 KL penalty on the first marginal (> 0).
//' @param rho2 KL penalty on the second marginal (> 0).
//' @param max_iter Maximum number of iterations.
//' @param tol Stopping tolerance on the infinity-norm iterate difference.
//' @return A list with translation-invariant potentials `fbar`, `gbar`,
//'   translated dual potentials `f`, `g`, translation `lambda`, iteration
//'   count, convergence flag, and last residual.
//' @export
// [[Rcpp::export]]
Rcpp::List uot_ti_sinkhorn_kl_sparse_cpp(const Rcpp::IntegerVector& row_ptr,
                                        const Rcpp::IntegerVector& col_idx,
                                        const Rcpp::NumericVector& cost,
                                        int n_rows,
                                        int n_cols,
                                        const arma::vec& alpha,
                                        const arma::vec& beta,
                                        double epsilon,
                                        double rho1,
                                        double rho2,
                                        int max_iter = 2000,
                                        double tol = 1e-6) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (rho1 <= 0.0) stop("rho1 must be > 0");
  if (rho2 <= 0.0) stop("rho2 must be > 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (row_ptr.size() != n_rows + 1) stop("row_ptr must have length n_rows + 1");
  if (col_idx.size() != cost.size()) stop("col_idx and cost must have the same length");

  const int nnz = col_idx.size();

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);

  arma::vec fbar(n_rows, arma::fill::zeros);
  arma::vec gbar(n_cols, arma::fill::zeros);
  arma::vec ftmp(n_rows);
  arma::vec gtmp(n_cols);
  arma::vec fbar_new(n_rows);
  arma::vec gbar_new(n_cols);

  const double denom = epsilon + rho1 + rho2;
  const double xi12 = epsilon * rho2 / (rho1 * denom);
  const double xi21 = epsilon * rho1 / (rho2 * denom);
  const double k11 = (epsilon / (epsilon + rho1)) * (rho1 / (rho1 + rho2));
  const double k22 = (epsilon / (epsilon + rho2)) * (rho2 / (rho1 + rho2));
  const double a1 = rho1 / (rho1 + epsilon);
  const double a2 = rho2 / (rho2 + epsilon);
  const double inv_eps = 1.0 / epsilon;

  bool converged = false;
  double residual = NA_REAL;
  int it = 0;

  // Determine whether row_ptr is 0-based or 1-based.
  int rp0 = row_ptr[0];
  int offset = (rp0 == 0) ? 0 : 1;
  if (!(rp0 == 0 || rp0 == 1)) stop("row_ptr[0] must be 0 or 1");

  for (it = 0; it < max_iter; ++it) {
    double smin_beta_rho2_gbar = softmin_scalar(logbeta, gbar, rho2);

    // fbar update
    for (int i = 0; i < n_rows; ++i) {
      int start = row_ptr[i] - offset;
      int end   = row_ptr[i + 1] - offset;
      if (start < 0 || end < start || end > nnz) stop("Invalid row_ptr offsets");
      double max_z = neg_inf();
      double sum_exp = 0.0;
      for (int idx = start; idx < end; ++idx) {
        int j = col_idx[idx] - 1;
        if (j < 0 || j >= n_cols) stop("col_idx out of bounds");
        double lb = logbeta(j);
        if (is_neg_inf(lb)) continue;
        double z = lb + (gbar(j) - cost[idx]) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z) {
          double scale = std::isfinite(max_z) ? std::exp(max_z - z) : 0.0;
          sum_exp = sum_exp * scale + 1.0;
          max_z = z;
        } else {
          sum_exp += std::exp(z - max_z);
        }
      }
      double logsum = (sum_exp > 0.0) ? (max_z + std::log(sum_exp)) : neg_inf();
      double smin_eps_row = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      ftmp(i) = a1 * smin_eps_row - k11 * smin_beta_rho2_gbar;
    }

    double smin_alpha_rho1_ftmp = softmin_scalar(logalpha, ftmp, rho1);
    fbar_new = ftmp + xi12 * smin_alpha_rho1_ftmp;

    double smin_alpha_rho1_fbar = softmin_scalar(logalpha, fbar_new, rho1);

    // gbar update: column-wise softmin via a single scan over edges
    arma::vec max_z(n_cols);
    arma::vec sum_exp(n_cols);
    max_z.fill(neg_inf());
    sum_exp.zeros();

    for (int i = 0; i < n_rows; ++i) {
      double la = logalpha(i);
      if (is_neg_inf(la)) continue;
      int start = row_ptr[i] - offset;
      int end   = row_ptr[i + 1] - offset;
      for (int idx = start; idx < end; ++idx) {
        int j = col_idx[idx] - 1;
        if (j < 0 || j >= n_cols) stop("col_idx out of bounds");
        double z = la + (fbar_new(i) - cost[idx]) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z(j)) {
          double scale = std::isfinite(max_z(j)) ? std::exp(max_z(j) - z) : 0.0;
          sum_exp(j) = sum_exp(j) * scale + 1.0;
          max_z(j) = z;
        } else {
          sum_exp(j) += std::exp(z - max_z(j));
        }
      }
    }

    for (int j = 0; j < n_cols; ++j) {
      double logsum = (sum_exp(j) > 0.0) ? (max_z(j) + std::log(sum_exp(j))) : neg_inf();
      double smin_eps_col = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      gtmp(j) = a2 * smin_eps_col - k22 * smin_alpha_rho1_fbar;
    }

    double smin_beta_rho2_gtmp = softmin_scalar(logbeta, gtmp, rho2);
    gbar_new = gtmp + xi21 * smin_beta_rho2_gtmp;

    residual = std::max(max_abs_diff(fbar_new, fbar), max_abs_diff(gbar_new, gbar));
    fbar.swap(fbar_new);
    gbar.swap(gbar_new);

    if (residual <= tol) {
      converged = true;
      break;
    }
  }

  double La = logsumexp_affine(logalpha, fbar, -1.0 / rho1);
  double Lb = logsumexp_affine(logbeta, gbar, -1.0 / rho2);
  double lambda = (rho1 * rho2) / (rho1 + rho2) * (La - Lb);

  arma::vec f = fbar + lambda;
  arma::vec g = gbar - lambda;

  return Rcpp::List::create(
    _["fbar"] = fbar,
    _["gbar"] = gbar,
    _["lambda"] = lambda,
    _["f"] = f,
    _["g"] = g,
    _["iterations"] = it + 1,
    _["converged"] = converged,
    _["residual"] = residual,
    _["backend"] = "sparse_csr_scan"
  );
}

//' KL translation-invariant Sinkhorn (sparse CSR+CSC cost)
//'
//' Sparse CPU implementation using both CSR and CSC representations.
//' This improves cache locality and enables per-column parallelism.
//'
//' @param row_ptr Integer vector of length n+1 with CSR row offsets (0-based or
//'   1-based accepted; if 1-based, it must start at 0 or 1 and be nondecreasing).
//' @param col_idx Integer vector of length nnz giving 1-based column indices.
//' @param cost Numeric vector of length nnz giving cost values aligned with
//'   `col_idx`.
//' @param col_ptr Integer vector of length m+1 with CSC column offsets (0-based
//'   or 1-based accepted).
//' @param row_idx Integer vector of length nnz giving 1-based row indices.
//' @param cost_csc Numeric vector of length nnz aligned with `row_idx`.
//' @param n_rows Number of source nodes (n).
//' @param n_cols Number of target nodes (m).
//' @param alpha Source masses (length n).
//' @param beta Target masses (length m).
//' @param epsilon Entropic regularization parameter (> 0).
//' @param rho1 KL penalty on the first marginal (> 0).
//' @param rho2 KL penalty on the second marginal (> 0).
//' @param max_iter Maximum number of iterations.
//' @param tol Stopping tolerance on the infinity-norm iterate difference.
//' @return A list with translation-invariant potentials `fbar`, `gbar`,
//'   translated dual potentials `f`, `g`, translation `lambda`, iteration
//'   count, convergence flag, last residual, and backend tag.
//' @export
// [[Rcpp::export]]
Rcpp::List uot_ti_sinkhorn_kl_sparse_csr_csc_cpp(const Rcpp::IntegerVector& row_ptr,
                                                const Rcpp::IntegerVector& col_idx,
                                                const Rcpp::NumericVector& cost,
                                                const Rcpp::IntegerVector& col_ptr,
                                                const Rcpp::IntegerVector& row_idx,
                                                const Rcpp::NumericVector& cost_csc,
                                                int n_rows,
                                                int n_cols,
                                                const arma::vec& alpha,
                                                const arma::vec& beta,
                                                double epsilon,
                                                double rho1,
                                                double rho2,
                                                int max_iter = 2000,
                                                double tol = 1e-6) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (rho1 <= 0.0) stop("rho1 must be > 0");
  if (rho2 <= 0.0) stop("rho2 must be > 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (row_ptr.size() != n_rows + 1) stop("row_ptr must have length n_rows + 1");
  if (col_ptr.size() != n_cols + 1) stop("col_ptr must have length n_cols + 1");
  if (col_idx.size() != cost.size()) stop("col_idx and cost must have the same length");
  if (row_idx.size() != cost_csc.size()) stop("row_idx and cost_csc must have the same length");
  if (col_idx.size() != row_idx.size()) stop("CSR and CSC must have same nnz");

  const int nnz = col_idx.size();

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);

  arma::vec fbar(n_rows, arma::fill::zeros);
  arma::vec gbar(n_cols, arma::fill::zeros);
  arma::vec ftmp(n_rows);
  arma::vec gtmp(n_cols);
  arma::vec fbar_new(n_rows);
  arma::vec gbar_new(n_cols);

  const double denom = epsilon + rho1 + rho2;
  const double xi12 = epsilon * rho2 / (rho1 * denom);
  const double xi21 = epsilon * rho1 / (rho2 * denom);
  const double k11 = (epsilon / (epsilon + rho1)) * (rho1 / (rho1 + rho2));
  const double k22 = (epsilon / (epsilon + rho2)) * (rho2 / (rho1 + rho2));
  const double a1 = rho1 / (rho1 + epsilon);
  const double a2 = rho2 / (rho2 + epsilon);
  const double inv_eps = 1.0 / epsilon;

  bool converged = false;
  double residual = NA_REAL;
  int it = 0;

  // Determine whether row_ptr / col_ptr are 0-based or 1-based.
  int rp0 = row_ptr[0];
  int row_offset = (rp0 == 0) ? 0 : 1;
  if (!(rp0 == 0 || rp0 == 1)) stop("row_ptr[0] must be 0 or 1");
  int cp0 = col_ptr[0];
  int col_offset = (cp0 == 0) ? 0 : 1;
  if (!(cp0 == 0 || cp0 == 1)) stop("col_ptr[0] must be 0 or 1");

  for (it = 0; it < max_iter; ++it) {
    double smin_beta_rho2_gbar = softmin_scalar(logbeta, gbar, rho2);

    // fbar update (CSR)
    for (int i = 0; i < n_rows; ++i) {
      int start = row_ptr[i] - row_offset;
      int end   = row_ptr[i + 1] - row_offset;
      if (start < 0 || end < start || end > nnz) stop("Invalid row_ptr offsets");
      double max_z = neg_inf();
      double sum_exp = 0.0;
      for (int idx = start; idx < end; ++idx) {
        int j = col_idx[idx] - 1;
        if (j < 0 || j >= n_cols) stop("col_idx out of bounds");
        double lb = logbeta(j);
        if (is_neg_inf(lb)) continue;
        double z = lb + (gbar(j) - cost[idx]) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z) {
          double scale = std::isfinite(max_z) ? std::exp(max_z - z) : 0.0;
          sum_exp = sum_exp * scale + 1.0;
          max_z = z;
        } else {
          sum_exp += std::exp(z - max_z);
        }
      }
      double logsum = (sum_exp > 0.0) ? (max_z + std::log(sum_exp)) : neg_inf();
      double smin_eps_row = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      ftmp(i) = a1 * smin_eps_row - k11 * smin_beta_rho2_gbar;
    }

    double smin_alpha_rho1_ftmp = softmin_scalar(logalpha, ftmp, rho1);
    fbar_new = ftmp + xi12 * smin_alpha_rho1_ftmp;

    double smin_alpha_rho1_fbar = softmin_scalar(logalpha, fbar_new, rho1);

    // gbar update (CSC)
    for (int j = 0; j < n_cols; ++j) {
      int start = col_ptr[j] - col_offset;
      int end   = col_ptr[j + 1] - col_offset;
      if (start < 0 || end < start || end > nnz) stop("Invalid col_ptr offsets");
      double max_z = neg_inf();
      double sum_exp = 0.0;
      for (int idx = start; idx < end; ++idx) {
        int i = row_idx[idx] - 1;
        if (i < 0 || i >= n_rows) stop("row_idx out of bounds");
        double la = logalpha(i);
        if (is_neg_inf(la)) continue;
        double z = la + (fbar_new(i) - cost_csc[idx]) * inv_eps;
        if (!std::isfinite(z)) continue;
        if (z > max_z) {
          double scale = std::isfinite(max_z) ? std::exp(max_z - z) : 0.0;
          sum_exp = sum_exp * scale + 1.0;
          max_z = z;
        } else {
          sum_exp += std::exp(z - max_z);
        }
      }
      double logsum = (sum_exp > 0.0) ? (max_z + std::log(sum_exp)) : neg_inf();
      double smin_eps_col = is_neg_inf(logsum) ? std::numeric_limits<double>::infinity()
                                               : (-epsilon * logsum);
      gtmp(j) = a2 * smin_eps_col - k22 * smin_alpha_rho1_fbar;
    }

    double smin_beta_rho2_gtmp = softmin_scalar(logbeta, gtmp, rho2);
    gbar_new = gtmp + xi21 * smin_beta_rho2_gtmp;

    residual = std::max(max_abs_diff(fbar_new, fbar), max_abs_diff(gbar_new, gbar));
    fbar.swap(fbar_new);
    gbar.swap(gbar_new);

    if (residual <= tol) {
      converged = true;
      break;
    }
  }

  double La = logsumexp_affine(logalpha, fbar, -1.0 / rho1);
  double Lb = logsumexp_affine(logbeta, gbar, -1.0 / rho2);
  double lambda = (rho1 * rho2) / (rho1 + rho2) * (La - Lb);

  arma::vec f = fbar + lambda;
  arma::vec g = gbar - lambda;

  return Rcpp::List::create(
    _["fbar"] = fbar,
    _["gbar"] = gbar,
    _["lambda"] = lambda,
    _["f"] = f,
    _["g"] = g,
    _["iterations"] = it + 1,
    _["converged"] = converged,
    _["residual"] = residual,
    _["backend"] = "sparse_csr_csc"
  );
}
