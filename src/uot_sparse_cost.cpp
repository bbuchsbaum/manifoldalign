#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Build a sparse cost graph from neighbour indices/distances (CSR+CSC)
//'
//' Helper for CPU pipelines: given neighbor indices and distances (e.g. from
//' `RANN::nn2`), compute a sparse bipartite cost graph between source nodes and
//' a fixed template.
//'
//' @param nn_idx Integer matrix (n x k) of 1-based neighbor indices into the
//'   template support.
//' @param nn_dists Numeric matrix (n x k) of Euclidean distances aligned with
//'   `nn_idx`.
//' @param n_cols Number of template nodes (m).
//' @param lambda_anat Weight for anatomical distance-squared term.
//' @param F Optional source features (n x D).
//' @param G Optional template features (m x D).
//' @param lambda_feat Weight for feature distance-squared term.
//' @param max_dist Optional maximum distance threshold; edges with
//'   `nn_dists > max_dist` are dropped.
//' @param extra_i Optional 1-based row indices for additional edges.
//' @param extra_j Optional 1-based column indices for additional edges.
//' @param extra_dists Optional distances for additional edges (aligned to
//'   `extra_i`, `extra_j`).
//' @param Y_template Optional template coordinates/features (m x d). When
//'   provided together with `prior_map` and `lambda_prior > 0`, adds a soft
//'   diagonal/identity bias term based on squared distances in template space.
//' @param prior_map Optional integer vector (length n) giving the preferred
//'   template index (1..m) for each source node. Use NA for “no prior” on a row.
//' @param lambda_prior Nonnegative weight for the soft prior term.
//' @param prior_sigma Optional positive scale (in the units of `Y_template`).
//'   When finite, the prior penalty uses `||y_j - y_pref||^2 / prior_sigma^2`.
//' @return A list with CSR fields `row_ptr`, `col_idx`, `cost` and CSC fields
//'   `col_ptr`, `row_idx`, `cost_csc`, plus dimensions.
//' @export
// [[Rcpp::export]]
Rcpp::List uot_build_sparse_cost_knn_cpp(const Rcpp::IntegerMatrix& nn_idx,
                                        const Rcpp::NumericMatrix& nn_dists,
                                        int n_cols,
                                        double lambda_anat = 1.0,
                                        Rcpp::Nullable<arma::mat> F = R_NilValue,
                                        Rcpp::Nullable<arma::mat> G = R_NilValue,
                                        double lambda_feat = 0.0,
                                        double max_dist = std::numeric_limits<double>::infinity(),
                                        Rcpp::Nullable<Rcpp::IntegerVector> extra_i = R_NilValue,
                                        Rcpp::Nullable<Rcpp::IntegerVector> extra_j = R_NilValue,
                                        Rcpp::Nullable<Rcpp::NumericVector> extra_dists = R_NilValue,
                                        Rcpp::Nullable<arma::mat> Y_template = R_NilValue,
                                        Rcpp::Nullable<Rcpp::IntegerVector> prior_map = R_NilValue,
                                        double lambda_prior = 0.0,
                                        double prior_sigma = NA_REAL) {

  const int n_rows = nn_idx.nrow();
  const int k = nn_idx.ncol();
  if (nn_dists.nrow() != n_rows || nn_dists.ncol() != k) {
    stop("nn_dists must have the same dimensions as nn_idx");
  }
  if (n_cols <= 0) stop("n_cols must be positive");
  if (!std::isfinite(lambda_anat) || lambda_anat < 0.0) stop("lambda_anat must be >= 0");
  if (!std::isfinite(lambda_feat) || lambda_feat < 0.0) stop("lambda_feat must be >= 0");
  if (!std::isfinite(lambda_prior) || lambda_prior < 0.0) stop("lambda_prior must be >= 0");
  if (!(std::isfinite(max_dist) || std::isinf(max_dist)) || max_dist <= 0.0) {
    stop("max_dist must be a positive finite value or Inf");
  }

  bool have_feat = (F.isNotNull() && G.isNotNull() && lambda_feat > 0.0);
  arma::mat Fmat;
  arma::mat Gmat;
  int D = 0;
  if (have_feat) {
    Fmat = Rcpp::as<arma::mat>(F);
    Gmat = Rcpp::as<arma::mat>(G);
    if (Fmat.n_rows != static_cast<arma::uword>(n_rows)) stop("F must have n rows");
    if (Gmat.n_rows != static_cast<arma::uword>(n_cols)) stop("G must have n_cols rows");
    if (Fmat.n_cols != Gmat.n_cols) stop("F and G must have the same number of columns");
    D = static_cast<int>(Fmat.n_cols);
    if (D == 0) have_feat = false;
  }

  bool have_prior = (prior_map.isNotNull() && lambda_prior > 0.0);
  arma::mat Ymat;
  Rcpp::IntegerVector prior_vec;
  int DY = 0;
  double prior_scale = 1.0;
  if (have_prior) {
    if (Y_template.isNull()) stop("Y_template must be provided when lambda_prior > 0");
    Ymat = Rcpp::as<arma::mat>(Y_template);
    if (Ymat.n_rows != static_cast<arma::uword>(n_cols)) {
      stop("Y_template must have n_cols rows");
    }
    DY = static_cast<int>(Ymat.n_cols);
    if (DY == 0) stop("Y_template must have at least one column");

    prior_vec = Rcpp::as<Rcpp::IntegerVector>(prior_map);
    if (prior_vec.size() != n_rows) stop("prior_map must have length n_rows");

    if (std::isfinite(prior_sigma)) {
      if (!(prior_sigma > 0.0)) stop("prior_sigma must be > 0");
      prior_scale = 1.0 / (prior_sigma * prior_sigma);
    }
  }

  // Validate extra edges, if provided.
  Rcpp::IntegerVector extra_i_vec;
  Rcpp::IntegerVector extra_j_vec;
  Rcpp::NumericVector extra_d_vec;
  bool have_extra = extra_i.isNotNull() || extra_j.isNotNull() || extra_dists.isNotNull();
  if (have_extra) {
    if (extra_i.isNull() || extra_j.isNull() || extra_dists.isNull()) {
      stop("extra_i, extra_j, and extra_dists must be provided together");
    }
    extra_i_vec = Rcpp::as<Rcpp::IntegerVector>(extra_i);
    extra_j_vec = Rcpp::as<Rcpp::IntegerVector>(extra_j);
    extra_d_vec = Rcpp::as<Rcpp::NumericVector>(extra_dists);
    if (extra_i_vec.size() != extra_j_vec.size() || extra_i_vec.size() != extra_d_vec.size()) {
      stop("extra_i, extra_j, and extra_dists must have the same length");
    }
  }

  // First pass: count nnz per row/column (after filtering).
  std::vector<int> row_counts(n_rows, 0);
  std::vector<int> col_counts(n_cols, 0);

  for (int i = 0; i < n_rows; ++i) {
    for (int t = 0; t < k; ++t) {
      int j = nn_idx(i, t);
      if (IntegerVector::is_na(j) || j <= 0) continue;
      if (j > n_cols) stop("nn_idx contains index > n_cols");
      double d = nn_dists(i, t);
      if (!std::isfinite(d)) continue;
      if (d > max_dist) continue;
      row_counts[i] += 1;
      col_counts[j - 1] += 1;
    }
  }

  std::vector<int> extra_row_counts(n_rows, 0);
  std::vector<int> extra_i0;
  std::vector<int> extra_j0;
  std::vector<double> extra_d0;
  if (have_extra) {
    extra_i0.reserve(extra_i_vec.size());
    extra_j0.reserve(extra_i_vec.size());
    extra_d0.reserve(extra_i_vec.size());
    for (R_xlen_t e = 0; e < extra_i_vec.size(); ++e) {
      int ii = extra_i_vec[e];
      int jj = extra_j_vec[e];
      if (IntegerVector::is_na(ii) || IntegerVector::is_na(jj)) continue;
      if (ii <= 0 || ii > n_rows) stop("extra_i out of bounds");
      if (jj <= 0 || jj > n_cols) stop("extra_j out of bounds");
      double d = extra_d_vec[e];
      if (!std::isfinite(d)) continue;
      if (d > max_dist) continue;
      int i0 = ii - 1;
      int j0 = jj - 1;
      extra_row_counts[i0] += 1;
      col_counts[j0] += 1;
      extra_i0.push_back(i0);
      extra_j0.push_back(j0);
      extra_d0.push_back(d);
    }
  }

  Rcpp::IntegerVector row_ptr(n_rows + 1);
  row_ptr[0] = 0;
  for (int i = 0; i < n_rows; ++i) {
    row_ptr[i + 1] = row_ptr[i] + row_counts[i] + extra_row_counts[i];
  }
  const int nnz_total = row_ptr[n_rows];

  Rcpp::IntegerVector col_ptr(n_cols + 1);
  col_ptr[0] = 0;
  for (int j = 0; j < n_cols; ++j) {
    col_ptr[j + 1] = col_ptr[j] + col_counts[j];
  }

  Rcpp::IntegerVector col_idx(nnz_total);
  Rcpp::NumericVector cost(nnz_total);

  std::vector<int> row_next(n_rows);
  for (int i = 0; i < n_rows; ++i) {
    row_next[i] = row_ptr[i];
  }

  auto compute_cost = [&](int i0, int j0, double d) -> double {
    double c = 0.0;
    if (lambda_anat > 0.0) {
      c += lambda_anat * d * d;
    }
    if (have_feat) {
      double s = 0.0;
      for (int dd = 0; dd < D; ++dd) {
        double diff = Fmat(i0, dd) - Gmat(j0, dd);
        s += diff * diff;
      }
      c += lambda_feat * s;
    }
    if (have_prior) {
      int pj = prior_vec[i0];
      if (!Rcpp::IntegerVector::is_na(pj) && pj > 0 && pj <= n_cols) {
        int p0 = pj - 1;
        double s = 0.0;
        for (int dd = 0; dd < DY; ++dd) {
          double diff = Ymat(j0, dd) - Ymat(p0, dd);
          s += diff * diff;
        }
        c += lambda_prior * prior_scale * s;
      }
    }
    return c;
  };

  // Second pass: fill CSR from nn matrices.
  for (int i = 0; i < n_rows; ++i) {
    for (int t = 0; t < k; ++t) {
      int j = nn_idx(i, t);
      if (IntegerVector::is_na(j) || j <= 0) continue;
      double d = nn_dists(i, t);
      if (!std::isfinite(d)) continue;
      if (d > max_dist) continue;
      int pos = row_next[i]++;
      int j0 = j - 1;
      col_idx[pos] = j;
      cost[pos] = compute_cost(i, j0, d);
    }
  }

  // Fill CSR from extra edges.
  if (!extra_i0.empty()) {
    for (size_t e = 0; e < extra_i0.size(); ++e) {
      int i0 = extra_i0[e];
      int j0 = extra_j0[e];
      double d = extra_d0[e];
      int pos = row_next[i0]++;
      col_idx[pos] = j0 + 1;
      cost[pos] = compute_cost(i0, j0, d);
    }
  }

  // Build CSC by transposing CSR.
  Rcpp::IntegerVector row_idx(nnz_total);
  Rcpp::NumericVector cost_csc(nnz_total);
  std::vector<int> col_next(n_cols);
  for (int j = 0; j < n_cols; ++j) {
    col_next[j] = col_ptr[j];
  }
  for (int i = 0; i < n_rows; ++i) {
    int start = row_ptr[i];
    int end = row_ptr[i + 1];
    for (int p = start; p < end; ++p) {
      int j0 = col_idx[p] - 1;
      int q = col_next[j0]++;
      row_idx[q] = i + 1; // 1-based
      cost_csc[q] = cost[p];
    }
  }

  return Rcpp::List::create(
    _["row_ptr"] = row_ptr,
    _["col_idx"] = col_idx,
    _["cost"] = cost,
    _["col_ptr"] = col_ptr,
    _["row_idx"] = row_idx,
    _["cost_csc"] = cost_csc,
    _["n_rows"] = n_rows,
    _["n_cols"] = n_cols,
    _["nnz"] = nnz_total
  );
}
