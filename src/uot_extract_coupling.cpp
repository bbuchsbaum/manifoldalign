#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

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

inline double log_double_max() {
  static const double v = std::log(std::numeric_limits<double>::max());
  return v;
}

} // namespace

//' Extract a sparse coupling / map operator from UOT dual potentials (CSC)
//'
//' Given translation-invariant dual potentials (fbar,gbar), compute edge weights
//' on a sparse cost graph in CSC form and return a `dgCMatrix`.
//'
//' Supported weight types:
//' - weight_type = 0: map operator weights, w_ij / (pi2_j + delta)
//' - weight_type = 1: conditional weights, w_ij / pi2_j (column-stochastic)
//' - weight_type = 2: raw coupling weights, w_ij
//'
//' Pruning is column-wise:
//' - keep only weights >= threshold (after weight_type transform)
//' - optionally keep only the topk_col largest weights per column
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
//' @param weight_type See above (0,1,2).
//' @param delta Stabilizer for weight_type=0.
//' @param topk_col Optional top-k pruning per column (0 keeps all).
//' @param threshold Optional threshold pruning (<= 0 keeps all).
//' @return A list with `coupling` (dgCMatrix) and `log_pi2` (length n_cols).
//' @export
// [[Rcpp::export]]
Rcpp::List uot_extract_coupling_sparse_csc_cpp(const Rcpp::IntegerVector& col_ptr,
                                              const Rcpp::IntegerVector& row_idx,
                                              const Rcpp::NumericVector& cost_csc,
                                              int n_rows,
                                              int n_cols,
                                              const arma::vec& alpha,
                                              const arma::vec& beta,
                                              const arma::vec& fbar,
                                              const arma::vec& gbar,
                                              double epsilon,
                                              int weight_type = 0,
                                              double delta = 1e-8,
                                              int topk_col = 0,
                                              double threshold = 0.0) {

  if (epsilon <= 0.0) stop("epsilon must be > 0");
  if (delta < 0.0) stop("delta must be >= 0");
  if (n_rows <= 0 || n_cols <= 0) stop("n_rows and n_cols must be positive");
  if (alpha.n_elem != static_cast<arma::uword>(n_rows)) stop("alpha length must match n_rows");
  if (beta.n_elem != static_cast<arma::uword>(n_cols)) stop("beta length must match n_cols");
  if (fbar.n_elem != static_cast<arma::uword>(n_rows)) stop("fbar length must match n_rows");
  if (gbar.n_elem != static_cast<arma::uword>(n_cols)) stop("gbar length must match n_cols");
  if (col_ptr.size() != n_cols + 1) stop("col_ptr must have length n_cols + 1");
  if (row_idx.size() != cost_csc.size()) stop("row_idx and cost_csc must have the same length");
  if (!(weight_type == 0 || weight_type == 1 || weight_type == 2)) {
    stop("weight_type must be 0 (map), 1 (conditional), or 2 (coupling)");
  }
  if (topk_col < 0) stop("topk_col must be >= 0");

  const int nnz = row_idx.size();
  const double inv_eps = 1.0 / epsilon;
  const double log_max = log_double_max();

  arma::vec logalpha = safe_log_weights(alpha);
  arma::vec logbeta  = safe_log_weights(beta);

  // Determine whether col_ptr is 0-based or 1-based.
  int cp0 = col_ptr[0];
  int offset = (cp0 == 0) ? 0 : 1;
  if (!(cp0 == 0 || cp0 == 1)) stop("col_ptr[0] must be 0 or 1");

  arma::vec max_logw(n_cols);
  arma::vec sum_scaled(n_cols);
  arma::vec log_pi2(n_cols);
  max_logw.fill(neg_inf());
  sum_scaled.zeros();
  log_pi2.fill(neg_inf());

  // Pass 1: compute per-column log_pi2 using stable max+sum scaling.
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
    if (ss > 0.0 && std::isfinite(mx)) {
      log_pi2(j) = mx + std::log(ss);
    }
  }

  std::vector<int> out_p(n_cols + 1);
  out_p[0] = 0;
  std::vector<int> out_i;
  std::vector<double> out_x;
  out_i.reserve(static_cast<size_t>(nnz));
  out_x.reserve(static_cast<size_t>(nnz));

  int pos = 0;
  for (int j = 0; j < n_cols; ++j) {
    int start = col_ptr[j] - offset;
    int end   = col_ptr[j + 1] - offset;
    if (start < 0 || end < start || end > nnz) stop("Invalid col_ptr offsets");
    if (start == end) {
      out_p[j + 1] = pos;
      continue;
    }

    double lb = logbeta(j);
    if (is_neg_inf(lb)) {
      out_p[j + 1] = pos;
      continue;
    }

    double lp2 = log_pi2(j);
    if (!std::isfinite(lp2)) {
      out_p[j + 1] = pos;
      continue;
    }

    double inv_denom_factor = 1.0;
    if (weight_type == 0 && delta > 0.0) {
      double z = -lp2;
      // If exp(z) overflows, denom_factor is huge and inv factor ~ 0.
      if (z > 700.0) {
        inv_denom_factor = 0.0;
      } else {
        inv_denom_factor = 1.0 / (1.0 + delta * std::exp(z));
      }
    }

    std::vector<std::pair<int, double>> items;
    items.reserve(static_cast<size_t>(end - start));

    for (int idx = start; idx < end; ++idx) {
      int i0 = row_idx[idx] - 1;
      if (i0 < 0 || i0 >= n_rows) continue;
      double la = logalpha(i0);
      if (is_neg_inf(la)) continue;
      double logw = la + lb + (fbar(i0) + gbar(j) - cost_csc[idx]) * inv_eps;
      if (!std::isfinite(logw)) continue;

      double val = 0.0;
      if (weight_type == 0) {
        val = std::exp(logw - lp2) * inv_denom_factor;
      } else if (weight_type == 1) {
        val = std::exp(logw - lp2);
      } else { // raw coupling
        if (logw >= log_max) {
          val = std::numeric_limits<double>::max();
        } else {
          val = std::exp(logw);
        }
      }

      if (threshold > 0.0 && val < threshold) continue;
      if (!(val > 0.0) || !std::isfinite(val)) continue;
      items.emplace_back(i0, val);
    }

    if (items.empty()) {
      out_p[j + 1] = pos;
      continue;
    }

    if (topk_col > 0 && static_cast<int>(items.size()) > topk_col) {
      auto mid = items.begin() + topk_col;
      std::nth_element(items.begin(), mid, items.end(),
                       [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                         return a.second > b.second;
                       });
      items.resize(static_cast<size_t>(topk_col));
    }

    std::sort(items.begin(), items.end(),
              [](const std::pair<int, double>& a, const std::pair<int, double>& b) {
                return a.first < b.first;
              });

    // Deduplicate within column (in case of repeated edges).
    std::vector<std::pair<int, double>> dedup;
    dedup.reserve(items.size());
    for (const auto& it : items) {
      if (!dedup.empty() && dedup.back().first == it.first) {
        dedup.back().second += it.second;
      } else {
        dedup.push_back(it);
      }
    }

    for (const auto& it : dedup) {
      out_i.push_back(it.first);  // 0-based for dgCMatrix
      out_x.push_back(it.second);
      pos += 1;
    }

    out_p[j + 1] = pos;
  }

  Rcpp::S4 mat("dgCMatrix");
  mat.slot("i") = Rcpp::wrap(out_i);
  mat.slot("p") = Rcpp::wrap(out_p);
  mat.slot("x") = Rcpp::wrap(out_x);
  mat.slot("Dim") = Rcpp::IntegerVector::create(n_rows, n_cols);
  mat.slot("Dimnames") = Rcpp::List::create(R_NilValue, R_NilValue);

  return Rcpp::List::create(
    _["coupling"] = mat,
    _["log_pi2"] = log_pi2
  );
}

