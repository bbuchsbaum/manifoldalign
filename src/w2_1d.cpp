// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

struct SortedVec {
  std::vector<double> x;
  std::vector<double> w;
};

// Sort the values of one column ascending, carrying the (validated) weights.
static inline SortedVec sort_with_weights(const NumericVector &vals, const NumericVector &weights) {
  int n = vals.size();
  std::vector<int> idx(n);
  for (int i = 0; i < n; ++i) idx[i] = i;
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b){ return vals[a] < vals[b]; });
  SortedVec out; out.x.resize(n); out.w.resize(n);
  for (int k = 0; k < n; ++k) {
    int i = idx[k];
    out.x[k] = vals[i];
    out.w[k] = weights[i];
  }
  return out;
}

// Compute 1-D weighted W2^2 between two sorted empirical measures (x,w)
static inline double w2_pair_sorted(const std::vector<double> &xa, const std::vector<double> &wa,
                                    const std::vector<double> &xb, const std::vector<double> &wb) {
  int na = (int)xa.size();
  int nb = (int)xb.size();
  int i = 0, j = 0;
  double wi = (na > 0 ? wa[0] : 0.0);
  double wj = (nb > 0 ? wb[0] : 0.0);
  double xi = (na > 0 ? xa[0] : 0.0);
  double xj = (nb > 0 ? xb[0] : 0.0);
  double res = 0.0;
  const double eps = 1e-15;
  while (i < na && j < nb) {
    double m = wi < wj ? wi : wj;
    double dx = xi - xj;
    res += (dx * dx) * m;
    wi -= m; wj -= m;
    if (wi <= eps) { ++i; if (i < na) { xi = xa[i]; wi = wa[i]; } }
    if (wj <= eps) { ++j; if (j < nb) { xj = xb[j]; wj = wb[j]; } }
  }
  return res;
}

// Validate a weight vector (finite, non-negative, positive total) and
// normalise it to sum to one.
static NumericVector normalised_weights(const NumericVector &w, const char *name) {
  double s = 0.0;
  for (R_xlen_t i = 0; i < w.size(); ++i) {
    double v = w[i];
    if (!R_finite(v)) stop("%s contains non-finite values", name);
    if (v < 0) stop("%s must be non-negative", name);
    s += v;
  }
  if (!(s > 0)) stop("%s must have a positive sum", name);
  NumericVector out = clone(w);
  for (R_xlen_t i = 0; i < out.size(); ++i) out[i] /= s;
  return out;
}

static void check_finite_matrix(const NumericMatrix &M, const char *name) {
  for (R_xlen_t i = 0; i < M.size(); ++i) {
    if (!R_finite(M[i])) stop("%s contains non-finite values (NA, NaN or Inf)", name);
  }
}

// [[Rcpp::export]]
Rcpp::List w2_cost_matrix_1d(Rcpp::NumericMatrix Uref,
                             Rcpp::NumericMatrix Ug,
                             Rcpp::NumericVector w_ref,
                             Rcpp::NumericVector w_g) {
  int n1 = Uref.nrow();
  int n2 = Ug.nrow();
  int K = Uref.ncol();
  if (Ug.ncol() != K) stop("Ug must have same number of columns as Uref");
  if (w_ref.size() != n1) stop("w_ref length mismatch");
  if (w_g.size() != n2) stop("w_g length mismatch");
  check_finite_matrix(Uref, "Uref");
  check_finite_matrix(Ug, "Ug");

  // Normalize weights to sum to 1
  NumericVector wr = normalised_weights(w_ref, "w_ref");
  NumericVector wg = normalised_weights(w_g, "w_g");

  // Pre-sort each column
  std::vector<SortedVec> Rsorted(K);
  std::vector<SortedVec> Gsorted(K);
  std::vector<SortedVec> Gsorted_neg(K);
  for (int k = 0; k < K; ++k) {
    Rsorted[k] = sort_with_weights(Uref(_, k), wr);
    Gsorted[k] = sort_with_weights(Ug(_, k), wg);
    // Negated sorted: ascending(-Ug) equals reversed ascending(Ug) with negated values
    int m = Gsorted[k].x.size();
    Gsorted_neg[k].x.resize(m);
    Gsorted_neg[k].w.resize(m);
    for (int i = 0; i < m; ++i) {
      int j = m - 1 - i;
      Gsorted_neg[k].x[i] = -Gsorted[k].x[j];
      Gsorted_neg[k].w[i] =  Gsorted[k].w[j];
    }
  }

  NumericMatrix C(K, K);
  IntegerMatrix sign(K, K);
  for (int k = 0; k < K; ++k) {
    for (int l = 0; l < K; ++l) {
      double cpos = w2_pair_sorted(Rsorted[k].x, Rsorted[k].w,
                                   Gsorted[l].x, Gsorted[l].w);
      double cneg = w2_pair_sorted(Rsorted[k].x, Rsorted[k].w,
                                   Gsorted_neg[l].x, Gsorted_neg[l].w);
      if (cpos <= cneg) { C(k,l) = cpos; sign(k,l) = +1; }
      else              { C(k,l) = cneg; sign(k,l) = -1; }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("C") = C,
    Rcpp::Named("sign") = sign
  );
}
