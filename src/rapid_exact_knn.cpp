// Exact bounded self-kNN for RAPID-MA relation construction.
//
// Every query row receives its `k` nearest other rows under Euclidean
// distance. Candidates are ordered by (distance, tie_rank), where `tie_rank`
// is a caller-supplied permutation (RAPID-MA passes the lexicographic rank of
// each row), so the result is unique and independent of row order even when
// distances tie.
//
// Rows are grouped into spatially coherent leaves by a median-split tree.
// For each query leaf, reference leaves are visited in increasing order of a
// conservative ball lower bound; a BLAS tile of inner products screens
// candidates, and every candidate that can enter a query's top-k is
// re-evaluated with the direct difference formula sum_d (q_d - r_d)^2, the
// same arithmetic as the ANN kd-tree used by RANN::nn2. Screening and pruning
// use bounds inflated well beyond floating-point error, so they never discard
// a row that could belong to the exact answer. Memory is O(n * k) for the
// answer plus one leaf-by-leaf tile; no n-by-n object is formed.

#include <RcppArmadillo.h>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct KnnLeaf {
  int start;
  int end;
};

// Recursive median split on the dimension of largest spread. The partition
// only affects speed: the returned neighbours do not depend on it.
void rapid_knn_build_leaves(const double* X, int n, int d,
                            std::vector<int>& order, int lo, int hi,
                            int leaf_size, std::vector<KnnLeaf>& leaves) {
  std::vector<std::pair<int, int>> stack;
  stack.push_back(std::make_pair(lo, hi));
  while (!stack.empty()) {
    std::pair<int, int> range = stack.back();
    stack.pop_back();
    int a = range.first;
    int b = range.second;
    if (b - a <= leaf_size) {
      leaves.push_back(KnnLeaf{a, b});
      continue;
    }
    int best_dim = 0;
    double best_spread = -1.0;
    for (int j = 0; j < d; ++j) {
      double mn = std::numeric_limits<double>::infinity();
      double mx = -std::numeric_limits<double>::infinity();
      const double* col = X + static_cast<std::size_t>(j) * n;
      for (int t = a; t < b; ++t) {
        double v = col[order[t]];
        if (v < mn) mn = v;
        if (v > mx) mx = v;
      }
      if (mx - mn > best_spread) {
        best_spread = mx - mn;
        best_dim = j;
      }
    }
    if (!(best_spread > 0.0)) {
      // All rows in this range coincide; keep them together in chunks.
      for (int s = a; s < b; s += leaf_size) {
        leaves.push_back(KnnLeaf{s, std::min(b, s + leaf_size)});
      }
      continue;
    }
    const double* col = X + static_cast<std::size_t>(best_dim) * n;
    int mid = a + (b - a) / 2;
    std::nth_element(
      order.begin() + a, order.begin() + mid, order.begin() + b,
      [col](int u, int v) {
        return col[u] < col[v] || (col[u] == col[v] && u < v);
      }
    );
    stack.push_back(std::make_pair(mid, b));
    stack.push_back(std::make_pair(a, mid));
  }
}

struct KnnHeap {
  // Max-heap per query on (dist, rank); slots [q * k, q * k + size[q]).
  int k;
  std::vector<double> dist;
  std::vector<double> d2;
  std::vector<int> rank;
  std::vector<int> index;
  std::vector<int> size;

  KnnHeap(int n, int k_) : k(k_), dist(static_cast<std::size_t>(n) * k_),
    d2(static_cast<std::size_t>(n) * k_), rank(static_cast<std::size_t>(n) * k_),
    index(static_cast<std::size_t>(n) * k_), size(n, 0) {}

  inline bool worse(std::size_t a, std::size_t b) const {
    return dist[a] > dist[b] || (dist[a] == dist[b] && rank[a] > rank[b]);
  }

  inline void swap_slots(std::size_t a, std::size_t b) {
    std::swap(dist[a], dist[b]);
    std::swap(d2[a], d2[b]);
    std::swap(rank[a], rank[b]);
    std::swap(index[a], index[b]);
  }

  inline bool full(int q) const { return size[q] == k; }

  inline const double& top_d2(int q) const {
    return d2[static_cast<std::size_t>(q) * k];
  }

  // Candidate (dist, rank) beats the current worst kept entry?
  inline bool improves(int q, double cand_dist, int cand_rank) const {
    std::size_t top = static_cast<std::size_t>(q) * k;
    return cand_dist < dist[top] ||
      (cand_dist == dist[top] && cand_rank < rank[top]);
  }

  void push(int q, double cand_d2, double cand_dist, int cand_rank, int cand_index) {
    std::size_t base = static_cast<std::size_t>(q) * k;
    if (size[q] < k) {
      std::size_t pos = base + size[q];
      dist[pos] = cand_dist;
      d2[pos] = cand_d2;
      rank[pos] = cand_rank;
      index[pos] = cand_index;
      ++size[q];
      // Sift up.
      std::size_t child = pos - base;
      while (child > 0) {
        std::size_t parent = (child - 1) / 2;
        if (worse(base + child, base + parent)) {
          swap_slots(base + child, base + parent);
          child = parent;
        } else {
          break;
        }
      }
      return;
    }
    dist[base] = cand_dist;
    d2[base] = cand_d2;
    rank[base] = cand_rank;
    index[base] = cand_index;
    std::size_t node = 0;
    const std::size_t m = static_cast<std::size_t>(k);
    while (true) {
      std::size_t left = 2 * node + 1;
      if (left >= m) break;
      std::size_t right = left + 1;
      std::size_t largest = left;
      if (right < m && worse(base + right, base + left)) largest = right;
      if (worse(base + largest, base + node)) {
        swap_slots(base + largest, base + node);
        node = largest;
      } else {
        break;
      }
    }
  }
};

}  // namespace

// Exact k-nearest other rows for every row of `X` (rows are points).
//
// Returns 1-based neighbour indices and Euclidean distances as query-major
// vectors of length nrow(X) * k: entries (i - 1) * k + 1, ..., i * k belong
// to row i, sorted by increasing (distance, tie_rank).
// [[Rcpp::export]]
Rcpp::List rapid_exact_knn_cpp(const Rcpp::NumericMatrix& X,
                               const Rcpp::IntegerVector& tie_rank,
                               int k,
                               int leaf_size = 256) {
  const int n = X.nrow();
  const int d = X.ncol();
  if (n < 2) Rcpp::stop("At least two rows are required.");
  if (d < 1) Rcpp::stop("At least one column is required.");
  if (k < 1 || k > n - 1) Rcpp::stop("`k` must lie in 1..nrow(X) - 1.");
  if (tie_rank.size() != n) Rcpp::stop("`tie_rank` must have one entry per row.");
  if (leaf_size < 8) leaf_size = 8;
  const double* Xp = X.begin();
  for (R_xlen_t t = 0; t < X.size(); ++t) {
    if (!std::isfinite(Xp[t])) Rcpp::stop("`X` must contain only finite values.");
  }

  std::vector<int> order(n);
  for (int i = 0; i < n; ++i) order[i] = i;
  std::vector<KnnLeaf> leaves;
  rapid_knn_build_leaves(Xp, n, d, order, 0, n, leaf_size, leaves);
  const int n_leaves = static_cast<int>(leaves.size());

  // Leaf-ordered, point-contiguous copy (d x n, column = point).
  arma::mat P(d, n);
  for (int t = 0; t < n; ++t) {
    int row = order[t];
    for (int j = 0; j < d; ++j) {
      P(j, t) = Xp[static_cast<std::size_t>(j) * n + row];
    }
  }
  std::vector<double> norm2(n);
  for (int t = 0; t < n; ++t) {
    const double* p = P.colptr(t);
    double s = 0.0;
    for (int j = 0; j < d; ++j) s += p[j] * p[j];
    norm2[t] = s;
  }

  // Leaf centres and radii; the slack keeps ball bounds conservative.
  const double slack = 1e-9;
  arma::mat centre(d, n_leaves, arma::fill::zeros);
  std::vector<double> radius(n_leaves, 0.0);
  for (int l = 0; l < n_leaves; ++l) {
    const int a = leaves[l].start;
    const int b = leaves[l].end;
    double* c = centre.colptr(l);
    for (int t = a; t < b; ++t) {
      const double* p = P.colptr(t);
      for (int j = 0; j < d; ++j) c[j] += p[j];
    }
    for (int j = 0; j < d; ++j) c[j] /= static_cast<double>(b - a);
    double rmax = 0.0;
    for (int t = a; t < b; ++t) {
      const double* p = P.colptr(t);
      double s = 0.0;
      for (int j = 0; j < d; ++j) {
        double diff = p[j] - c[j];
        s += diff * diff;
      }
      if (s > rmax) rmax = s;
    }
    radius[l] = std::sqrt(rmax) * (1.0 + slack) + slack;
  }

  KnnHeap heap(n, k);
  std::vector<int> rank_by_position(n);
  for (int t = 0; t < n; ++t) rank_by_position[t] = tie_rank[order[t]];

  // Screening margin covers the rounding of both the Gram-based distance and
  // the direct difference formula (each below (d + 2) * eps * (|q|^2 + |r|^2)).
  const double margin_scale = 16.0 * (d + 4) * DBL_EPSILON;
  const double inflate = 1.0 + 8.0 * DBL_EPSILON;
  const double inf = std::numeric_limits<double>::infinity();

  std::vector<double> leaf_bound(n_leaves);
  std::vector<int> leaf_order(n_leaves);
  std::vector<double> query_threshold;
  std::vector<double> query_centre_bound;
  arma::mat G;
  const double keep_scale = 1.0 - margin_scale;
  std::vector<double> ref_scaled;
  std::vector<int> candidates;
  double tiles = 0.0;

  for (int ql = 0; ql < n_leaves; ++ql) {
    Rcpp::checkUserInterrupt();
    const int qa = leaves[ql].start;
    const int qb = leaves[ql].end;
    const int nq = qb - qa;
    const double* cq = centre.colptr(ql);
    for (int rl = 0; rl < n_leaves; ++rl) {
      const double* cr = centre.colptr(rl);
      double s = 0.0;
      for (int j = 0; j < d; ++j) {
        double diff = cq[j] - cr[j];
        s += diff * diff;
      }
      double gap = std::sqrt(s) * (1.0 - slack) - radius[ql] - radius[rl];
      leaf_bound[rl] = gap > 0.0 ? gap * gap * (1.0 - slack) : 0.0;
      leaf_order[rl] = rl;
    }
    // Visit the query's own leaf first: it fills every heap with close rows.
    leaf_bound[ql] = -1.0;
    std::sort(leaf_order.begin(), leaf_order.end(), [&](int u, int v) {
      return leaf_bound[u] < leaf_bound[v] ||
        (leaf_bound[u] == leaf_bound[v] && u < v);
    });
    query_threshold.assign(nq, inf);
    query_centre_bound.assign(nq, 0.0);
    arma::mat Aq(P.colptr(qa), d, nq, false, true);

    for (int pos = 0; pos < n_leaves; ++pos) {
      const int rl = leaf_order[pos];
      double worst = 0.0;
      for (int u = 0; u < nq; ++u) {
        if (query_threshold[u] > worst) worst = query_threshold[u];
      }
      if (leaf_bound[rl] > worst) break;
      const int ra = leaves[rl].start;
      const int rb = leaves[rl].end;
      const int nr = rb - ra;
      const double* cr = centre.colptr(rl);
      bool any_active = false;
      for (int u = 0; u < nq; ++u) {
        const double* q = P.colptr(qa + u);
        double s = 0.0;
        for (int j = 0; j < d; ++j) {
          double diff = q[j] - cr[j];
          s += diff * diff;
        }
        double gap = std::sqrt(s) * (1.0 - slack) - radius[rl];
        double bound = gap > 0.0 ? gap * gap * (1.0 - slack) : 0.0;
        query_centre_bound[u] = bound;
        if (bound <= query_threshold[u]) any_active = true;
      }
      if (!any_active) continue;

      arma::mat Ar(P.colptr(ra), d, nr, false, true);
      // Column u of G holds query u's inner products with every reference.
      G = Ar.t() * Aq;
      ++tiles;
      if (static_cast<int>(ref_scaled.size()) < nr) {
        ref_scaled.resize(nr);
        candidates.resize(nr);
      }
      for (int v = 0; v < nr; ++v) ref_scaled[v] = keep_scale * norm2[ra + v];
      for (int u = 0; u < nq; ++u) {
        if (query_centre_bound[u] > query_threshold[u]) continue;
        const int qpos = qa + u;
        const int qrow = order[qpos];
        const double qn = norm2[qpos];
        const double* q = P.colptr(qpos);
        double threshold = query_threshold[u];
        const double* g = G.colptr(u);
        // Screen: (1 - m)(|q|^2 + |r|^2) - 2 q.r <= threshold. A candidate
        // failing it has an exact squared distance above the threshold.
        double cutoff = threshold - keep_scale * qn;
        int n_candidates = 0;
        for (int v = 0; v < nr; ++v) {
          candidates[n_candidates] = v;
          n_candidates += (ref_scaled[v] - 2.0 * g[v] <= cutoff);
        }
        for (int t = 0; t < n_candidates; ++t) {
          const int v = candidates[t];
          const int rpos = ra + v;
          if (rpos == qpos) continue;
          if (ref_scaled[v] - 2.0 * g[v] > cutoff) continue;
          const double* r = P.colptr(rpos);
          double s = 0.0;
          for (int j = 0; j < d; ++j) {
            double diff = q[j] - r[j];
            s += diff * diff;
          }
          const double dist = std::sqrt(s);
          const int cand_rank = rank_by_position[rpos];
          if (heap.full(qrow) && !heap.improves(qrow, dist, cand_rank)) continue;
          heap.push(qrow, s, dist, cand_rank, order[rpos]);
          if (heap.full(qrow)) {
            threshold = heap.top_d2(qrow) * inflate;
            cutoff = threshold - keep_scale * qn;
          }
        }
        query_threshold[u] = threshold;
      }
    }
  }

  // Sort each query's heap into ascending (distance, rank) order.
  Rcpp::IntegerVector out_index(static_cast<R_xlen_t>(n) * k);
  Rcpp::NumericVector out_dist(static_cast<R_xlen_t>(n) * k);
  std::vector<int> slot(k);
  for (int q = 0; q < n; ++q) {
    if (heap.size[q] != k) Rcpp::stop("Internal error: incomplete neighbour list.");
    const std::size_t base = static_cast<std::size_t>(q) * k;
    for (int t = 0; t < k; ++t) slot[t] = t;
    std::sort(slot.begin(), slot.end(), [&](int u, int v) {
      std::size_t a = base + u;
      std::size_t b = base + v;
      return heap.dist[a] < heap.dist[b] ||
        (heap.dist[a] == heap.dist[b] && heap.rank[a] < heap.rank[b]);
    });
    for (int t = 0; t < k; ++t) {
      out_index[base + t] = heap.index[base + slot[t]] + 1;
      out_dist[base + t] = heap.dist[base + slot[t]];
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("index") = out_index,
    Rcpp::Named("distance") = out_dist,
    Rcpp::Named("leaves") = n_leaves,
    Rcpp::Named("tiles") = tiles
  );
}
