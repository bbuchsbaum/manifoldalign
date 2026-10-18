#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include "ot_common.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

namespace {

struct SinkhornResult {
  arma::mat P;
  arma::vec log_u;
  arma::vec log_v;
  int iterations;
  bool converged;
  double marginal_error;
};

void validate_sinkhorn_inputs(const arma::mat& cost,
                              const arma::vec& a,
                              const arma::vec& b,
                              double epsilon,
                              int max_iter,
                              double tol) {
  check_marginals(a, b);
  if (cost.n_rows != a.n_elem || cost.n_cols != b.n_elem) {
    Rcpp::stop("Sinkhorn: cost must be length(a) x length(b).");
  }
  if (!std::isfinite(epsilon) || epsilon <= 0.0) {
    Rcpp::stop("Sinkhorn: epsilon must be a finite positive number.");
  }
  if (max_iter < 1) {
    Rcpp::stop("Sinkhorn: max_iter must be >= 1.");
  }
  if (!std::isfinite(tol) || tol <= 0.0) {
    Rcpp::stop("Sinkhorn: tol must be a finite positive number.");
  }
  if (!cost.is_finite()) {
    Rcpp::stop("Sinkhorn: cost matrix contains non-finite values.");
  }
}

// Column-wise log-sum-exp of (log_K + log_u 1^T); returns a length-m vector.
arma::vec lse_cols(const arma::mat& log_K, const arma::vec& log_u, arma::mat& tmp) {
  tmp = log_K;
  tmp.each_col() += log_u;
  arma::rowvec mx = arma::max(tmp, 0);
  mx.replace(-arma::datum::inf, 0.0);
  tmp.each_row() -= mx;
  arma::rowvec s = arma::sum(arma::exp(tmp), 0);
  return (mx + arma::log(s)).t();
}

// Row-wise log-sum-exp of (log_K + 1 log_v^T); returns a length-n vector.
arma::vec lse_rows(const arma::mat& log_K, const arma::vec& log_v, arma::mat& tmp) {
  tmp = log_K;
  tmp.each_row() += log_v.t();
  arma::vec mx = arma::max(tmp, 1);
  mx.replace(-arma::datum::inf, 0.0);
  tmp.each_col() -= mx;
  arma::vec s = arma::sum(arma::exp(tmp), 1);
  return mx + arma::log(s);
}

double marginal_residual(const arma::mat& P, const arma::vec& a, const arma::vec& b) {
  double row_err = arma::abs(arma::sum(P, 1) - a).max();
  double col_err = arma::abs(arma::sum(P, 0).t() - b).max();
  double err = std::max(row_err, col_err);
  if (!std::isfinite(err)) err = std::numeric_limits<double>::infinity();
  return err;
}

// Balanced entropic OT:  min <C, P> + epsilon * sum P (log P - 1)  s.t. P1 = a, P^T1 = b.
//
// Both paths share one convention: the Gibbs kernel is exp(-cost / epsilon) on
// the raw (unscaled) cost, P = diag(exp(log_u)) K diag(exp(log_v)), and the
// iteration stops when the marginal residual max(|P1 - a|, |P^T 1 - b|) is
// below `tol`.  Each sweep updates v (columns exact) and then u (rows exact),
// so the residual is the column error, which is available for free from the
// next sweep's column log-sum-exp.  No post-hoc rebalancing is applied: the
// returned plan is exactly the plan implied by the returned potentials, and
// `converged` reports whether the tolerance was reached within max_iter.
SinkhornResult sinkhorn_core(const arma::mat& cost,
                             const arma::vec& a,
                             const arma::vec& b,
                             double epsilon,
                             int max_iter,
                             double tol,
                             bool stabilized,
                             const arma::vec* log_u0,
                             const arma::vec* log_v0) {
  validate_sinkhorn_inputs(cost, a, b, epsilon, max_iter, tol);

  const arma::uword n = a.n_elem;
  const arma::uword m = b.n_elem;

  SinkhornResult out;
  out.converged = false;
  out.iterations = 0;
  out.log_u = log_u0 ? *log_u0 : arma::vec(n, arma::fill::zeros);
  out.log_v = log_v0 ? *log_v0 : arma::vec(m, arma::fill::zeros);

  if (stabilized) {
    const arma::vec log_a = arma::log(a);
    const arma::vec log_b = arma::log(b);
    const arma::mat log_K = -cost / epsilon;
    arma::mat tmp(n, m, arma::fill::none);

    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec lc = lse_cols(log_K, out.log_u, tmp);
      if (iter > 0) {
        // Rows are exact after the previous u-update; check the columns.
        double col_err = arma::abs(arma::exp(out.log_v + lc) - b).max();
        if (std::isfinite(col_err) && col_err < tol) {
          out.converged = true;
          break;
        }
      }
      out.log_v = log_b - lc;
      arma::vec lr = lse_rows(log_K, out.log_v, tmp);
      out.log_u = log_a - lr;
      out.iterations = iter + 1;
    }

    tmp = log_K;
    tmp.each_col() += out.log_u;
    tmp.each_row() += out.log_v.t();
    out.P = arma::exp(tmp);
  } else {
    const arma::mat K = arma::exp(-cost / epsilon);
    arma::vec u = arma::ones(n);
    arma::vec v = arma::ones(m);
    if (log_u0) u = arma::exp(out.log_u - out.log_u.max());
    if (log_v0) v = arma::exp(out.log_v - out.log_v.max());

    for (int iter = 0; iter < max_iter; ++iter) {
      arma::vec KTu = K.t() * u;
      if (iter > 0) {
        double col_err = arma::abs(v % KTu - b).max();
        if (std::isfinite(col_err) && col_err < tol) {
          out.converged = true;
          break;
        }
      }
      v = b / (KTu + 1e-300);
      u = a / (K * v + 1e-300);
      out.iterations = iter + 1;
    }

    out.P = arma::diagmat(u) * K * arma::diagmat(v);
    out.log_u = arma::log(u);
    out.log_v = arma::log(v);
  }

  // Report convergence from the plan actually returned.
  out.marginal_error = marginal_residual(out.P, a, b);
  out.converged = out.marginal_error < tol;
  return out;
}

arma::vec nullable_to_vec(const Rcpp::Nullable<Rcpp::NumericVector>& x,
                          arma::uword expected,
                          const char* what) {
  Rcpp::NumericVector v(x);
  if (static_cast<arma::uword>(v.size()) != expected) {
    Rcpp::stop("%s has the wrong length.", what);
  }
  arma::vec out = Rcpp::as<arma::vec>(v);
  if (!out.is_finite()) {
    Rcpp::stop("%s must be finite.", what);
  }
  return out;
}

} // namespace

//' Unified Sinkhorn algorithm with optional log-domain stabilization
//'
//' Solves balanced entropic OT with Gibbs kernel \code{exp(-cost / epsilon)}
//' on the raw cost. With \code{stabilized = TRUE} the iterations run in the
//' log domain for every \code{epsilon}; with \code{stabilized = FALSE} they
//' use plain kernel scaling, which underflows when \code{cost / epsilon} is
//' large. Iterations stop when the marginal residual is below \code{tol}. The
//' plan is returned as computed (no post-hoc rebalancing); use
//' \code{sinkhorn_unified_potentials()} to obtain the convergence flag and
//' marginal residual.
//'
//' @param cost Cost matrix (n x m)
//' @param a Source marginal (n x 1)
//' @param b Target marginal (m x 1)
//' @param epsilon Entropic regularization parameter
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance on the marginal residual
//' @param stabilized Use log-domain iterations
//' @return Transport plan matrix
// [[Rcpp::export]]
arma::mat sinkhorn_unified(const arma::mat& cost,
                           const arma::vec& a,
                           const arma::vec& b,
                           double epsilon,
                           int max_iter = 1000,
                           double tol = 1e-9,
                           bool stabilized = true) {
  SinkhornResult res = sinkhorn_core(cost, a, b, epsilon, max_iter, tol,
                                     stabilized, nullptr, nullptr);
  return res.P;
}

//' Unified Sinkhorn algorithm returning dual potentials for warm-starting
//'
//' Same solver as \code{sinkhorn_unified()}, additionally returning the log
//' scaling vectors (\code{pi = exp(-cost/epsilon + log_u + log_v^T)} exactly)
//' and convergence diagnostics.
//'
//' @param cost Cost matrix (n x m)
//' @param a Source marginal (n x 1)
//' @param b Target marginal (m x 1)
//' @param epsilon Entropic regularization parameter
//' @param log_u0 Optional warm-start log scaling for rows (length n)
//' @param log_v0 Optional warm-start log scaling for columns (length m)
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance on the marginal residual
//' @param stabilized Use log-domain iterations
//' @return A list with elements \code{pi} (plan), \code{log_u}, \code{log_v},
//'   \code{converged}, \code{iterations} and \code{marginal_error}
//'   (max absolute row/column marginal violation of \code{pi}).
// [[Rcpp::export]]
Rcpp::List sinkhorn_unified_potentials(const arma::mat& cost,
                                      const arma::vec& a,
                                      const arma::vec& b,
                                      double epsilon,
                                      Rcpp::Nullable<Rcpp::NumericVector> log_u0 = R_NilValue,
                                      Rcpp::Nullable<Rcpp::NumericVector> log_v0 = R_NilValue,
                                      int max_iter = 1000,
                                      double tol = 1e-9,
                                      bool stabilized = true) {
  arma::vec u0, v0;
  const arma::vec* pu = nullptr;
  const arma::vec* pv = nullptr;
  if (log_u0.isNotNull()) {
    u0 = nullable_to_vec(log_u0, cost.n_rows, "log_u0");
    pu = &u0;
  }
  if (log_v0.isNotNull()) {
    v0 = nullable_to_vec(log_v0, cost.n_cols, "log_v0");
    pv = &v0;
  }

  SinkhornResult res = sinkhorn_core(cost, a, b, epsilon, max_iter, tol,
                                     stabilized, pu, pv);

  return Rcpp::List::create(
    Rcpp::_["pi"] = res.P,
    Rcpp::_["log_u"] = res.log_u,
    Rcpp::_["log_v"] = res.log_v,
    Rcpp::_["converged"] = res.converged,
    Rcpp::_["iterations"] = res.iterations,
    Rcpp::_["marginal_error"] = res.marginal_error
  );
}

// Maintain backward compatibility - redirect old functions to unified version

// [[Rcpp::export]]
arma::mat solve_sinkhorn_stabilized_cpp(const arma::mat& C_in, double tau,
                                        int max_iter, double tol) {
  // Uniform marginals, log-domain iterations on the raw cost.
  int n = C_in.n_rows;
  int m = C_in.n_cols;
  arma::vec mu(n, arma::fill::value(1.0 / std::max(1, n)));
  arma::vec nu(m, arma::fill::value(1.0 / std::max(1, m)));

  return sinkhorn_unified(C_in, mu, nu, tau, max_iter, tol, true);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sinkhorn_ot_cpp(Rcpp::NumericMatrix cost_r,
                                    Rcpp::NumericVector a_r,
                                    Rcpp::NumericVector b_r,
                                    double epsilon,
                                    int max_iter = 1000,
                                    double tol = 1e-9) {
  arma::mat cost = Rcpp::as<arma::mat>(cost_r);
  arma::vec a = Rcpp::as<arma::vec>(a_r);
  arma::vec b = Rcpp::as<arma::vec>(b_r);

  arma::mat P = sinkhorn_unified(cost, a, b, epsilon, max_iter, tol, true);

  return Rcpp::wrap(P);
}
