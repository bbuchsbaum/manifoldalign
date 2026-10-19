#include <RcppArmadillo.h>
#include <roptim.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

class LinearSimEmbed : public roptim::Functor {
public:
  LinearSimEmbed(const arma::mat& X, const arma::mat& T, const arma::mat& M, 
                 const double sigma_P, const int ncomp, const double alpha_p)
    : X_(X), T_(T), M_(M), sigma_P_(sigma_P), ncomp_(ncomp), alpha_p_(alpha_p) {
    
    if (!std::isfinite(sigma_P_) || sigma_P_ <= 0) {
      Rcpp::stop("sigma_P must be finite and positive");
    }
    
    // Pre-calculate sum of Mask for normalization
    sum_M_ = arma::accu(M_);
    if (sum_M_ < 1e-12) {
      Rcpp::stop("Mask matrix M sums to zero. Cannot optimize.");
    }
    
    // Cache dimensions
    n_ = X_.n_rows;
    d_ = X_.n_cols;
    
    // Initialize cache variables
    cached_w_.set_size(d_ * ncomp_);
    cached_w_.fill(arma::datum::nan);  // Invalid cache initially
  }
  
  // Efficiently calculate similarity matrix P = exp(-D²/sigma) using vectorized operations
  arma::mat calculate_P(const arma::mat& Y) {
    // Fast exit for single point case to prevent broadcasting bug
    if (n_ == 1) {
      return arma::ones<arma::mat>(1, 1);
    }
    
    // Vectorized computation of squared distances following R pairwise_sqdist pattern
    arma::vec Y_rss = arma::sum(arma::square(Y), 1);  // Row sum of squares
    arma::mat D2 = arma::repmat(Y_rss, 1, n_) + arma::repmat(Y_rss.t(), n_, 1) - 2.0 * Y * Y.t();
    
    // Ensure numerical stability (handle floating point errors)
    D2.elem(arma::find(D2 < 0)).zeros();
    
    // Compute similarities and enforce symmetry
    arma::mat P = arma::exp(-D2 / sigma_P_);
    P = 0.5 * (P + P.t());  // Enforce exact symmetry as required by paper
    
    return P;
  }
  
  // CORRECTED Objective Function: (1-alpha_p)*Js + alpha_p*Jp
  double operator()(const arma::vec& w) override {
    // Check for user interrupt periodically
    if (++iter_count_ % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // Cache Y and P if w has changed (tighter tolerance for consistency)
    if (!arma::approx_equal(w, cached_w_, "absdiff", 1e-15)) {
      cached_w_ = w;
      arma::mat W = arma::reshape(w, d_, ncomp_);
      cached_Y_ = X_ * W;
      cached_P_ = calculate_P(cached_Y_);
    }
    
    // Similarity loss Js (using consistent normalization)
    double Js = arma::accu(M_ % arma::square(cached_P_ - T_)) / (2.0 * sum_M_);
    
    // Orthogonality loss Jp (following paper's formula: ||W'W - I||²_F / (2m²))
    arma::mat W = arma::reshape(w, d_, ncomp_);
    arma::mat WtW = W.t() * W;
    // FIXED: Don't mutate WtW in-place
    arma::mat WWmI = WtW;
    WWmI.diag() -= 1.0;
    double Jp = arma::accu(arma::square(WWmI)) / (2.0 * ncomp_ * ncomp_);
    
    // CORRECTED: Use convex combination as per paper Eq. 5
    double objective = (1.0 - alpha_p_) * Js + alpha_p_ * Jp;
    
    return objective;
  }
  
  // CORRECTED and VECTORIZED Gradient Function
  void Gradient(const arma::vec& w, arma::vec& grad) override {
    arma::mat W = arma::reshape(w, d_, ncomp_);
    
    // Use cached Y and P if available, otherwise compute
    arma::mat Y, P;
    if (arma::approx_equal(w, cached_w_, "absdiff", 1e-15)) {
      Y = cached_Y_;
      P = cached_P_;
    } else {
      Y = X_ * W;
      P = calculate_P(Y);
    }
    
    // For Js = sum(M * (P - T)^2) / (2 * sum(M)), differentiating
    // the symmetric pairwise-distance matrix contributes both (i, j) and
    // (j, i). Keep this expression tied directly to that objective.
    arma::mat A = (P - T_) % M_ % P;
    
    // Efficiently compute gradient components using matrix operations
    arma::mat XtAy = X_.t() * (A * Y);                    // d × m
    arma::vec A_rowsums = arma::sum(A, 1);                // n × 1
    arma::mat XtDX_W = (X_.t() * (X_.each_col() % A_rowsums)) * W; // d × m
    
    arma::mat grad_Js = (4.0 / (sigma_P_ * sum_M_)) * (XtAy - XtDX_W);
    
    // d ||W'W - I||_F^2 / dW = 4 W(W'W - I).
    arma::mat WtW = W.t() * W;
    arma::mat WWmI = WtW;
    WWmI.diag() -= 1.0;
    arma::mat grad_Jp = (2.0 / (ncomp_ * ncomp_)) * W * WWmI;
    
    // CORRECTED: Use exact same convex combination as objective
    grad = arma::vectorise((1.0 - alpha_p_) * grad_Js + alpha_p_ * grad_Jp);
  }
  
private:
  const arma::mat& X_;
  const arma::mat& T_;
  const arma::mat& M_;
  int ncomp_;
  double sigma_P_;
  double alpha_p_;
  double sum_M_;
  size_t n_, d_;
  mutable int iter_count_ = 0;  // For progress reporting and interrupt checking
  
  // Cache variables for efficiency
  mutable arma::vec cached_w_;
  mutable arma::mat cached_Y_;
  mutable arma::mat cached_P_;
};

// [[Rcpp::export]]
Rcpp::List linear_sim_embed_cpp(const arma::mat& X, const arma::mat& T, const arma::mat& M,
                                double sigma_P, int ncomp, double alpha_p, 
                                int maxit = 500, double tol = 1e-6) {
  
  // Input validation
  if (X.n_rows != T.n_rows || X.n_rows != M.n_rows || T.n_rows != T.n_cols || M.n_rows != M.n_cols) {
    Rcpp::stop("Dimension mismatch: X, T, M must have consistent dimensions");
  }
  
  if (sigma_P <= 0) {
    Rcpp::stop("sigma_P must be positive");
  }
  
  if (ncomp <= 0 || ncomp > X.n_cols) {
    Rcpp::stop("ncomp must be positive and <= ncol(X)");
  }
  
  if (alpha_p < 0 || alpha_p > 1) {
    Rcpp::stop("alpha_p must be in [0, 1]");
  }

  if (maxit < 1) {
    Rcpp::stop("maxit must be >= 1");
  }

  if (!std::isfinite(tol) || tol <= 0) {
    Rcpp::stop("tol must be finite and positive");
  }

  if (!X.is_finite() || !T.is_finite() || !M.is_finite()) {
    Rcpp::stop("X, T, and M must contain only finite values");
  }

  if (arma::min(arma::vectorise(M)) < 0) {
    Rcpp::stop("M must contain non-negative weights");
  }

  const double symmetry_tol = 1e-10;
  const double T_scale = std::max(1.0, arma::norm(T, "inf"));
  const double M_scale = std::max(1.0, arma::norm(M, "inf"));
  if (arma::norm(T - T.t(), "inf") > symmetry_tol * T_scale) {
    Rcpp::stop("T must be symmetric for Linear Similarity Embedding");
  }
  if (arma::norm(M - M.t(), "inf") > symmetry_tol * M_scale) {
    Rcpp::stop("M must be symmetric for Linear Similarity Embedding");
  }
  
  // FIXED: Proper SVD initialization for n < d case
  arma::mat initial_W;
  if (ncomp <= std::min(X.n_rows, X.n_cols)) {
    // Use SVD for initialization with correct orientation
    arma::mat U, V;
    arma::vec s;
    bool svd_ok = arma::svd_econ(U, s, V, X);  // X not X.t() - V becomes d×n
    if (svd_ok && V.n_cols >= static_cast<arma::uword>(ncomp)) {
      initial_W = V.cols(0, ncomp - 1);  // Now shapes match: d × ncomp
    } else {
      // Deterministic orthonormal fallback; do not mutate either R or
      // Armadillo RNG state merely to initialize an optimizer.
      initial_W = arma::zeros<arma::mat>(X.n_cols, ncomp);
      initial_W.submat(0, 0, ncomp - 1, ncomp - 1).eye();
    }
  } else {
    Rcpp::stop("ncomp is too large for proper initialization");
  }
  
  // Create optimization instance
  LinearSimEmbed task(X, T, M, sigma_P, ncomp, alpha_p);
  
  // Create L-BFGS-B optimizer with proper settings
  roptim::Roptim<LinearSimEmbed> opt("L-BFGS-B");
  
  // Set optimization parameters
  opt.control.maxit = maxit;
  opt.control.trace = 0;  // We handle progress reporting in the objective function
  // L-BFGS-B ignores reltol. pgtol is its projected-gradient stopping
  // tolerance and is the quantity documented by the R interface.
  opt.control.pgtol = tol;
  opt.control.fnscale = 1.0;  // Ensure proper sign convention
  
  // FIXED: Disable Hessian for L-BFGS
  opt.set_hessian(false);
  
  // Perform optimization
  arma::vec initial_w = arma::vectorise(initial_W);
  const double initial_value = task(initial_w);
  arma::vec initial_gradient_vec;
  task.Gradient(initial_w, initial_gradient_vec);
  
  try {
    opt.minimize(task, initial_w);
  } catch (const std::exception& e) {
    Rcpp::stop("L-BFGS-B optimization failed: " + std::string(e.what()));
  }
  
  // Extract results
  arma::mat optimized_W = arma::reshape(opt.par(), X.n_cols, ncomp);
  arma::vec final_gradient_vec;
  task.Gradient(opt.par(), final_gradient_vec);
  
  return Rcpp::List::create(
    Rcpp::Named("W") = optimized_W,
    Rcpp::Named("convergence") = opt.convergence(),
    Rcpp::Named("message") = opt.message(),
    Rcpp::Named("value") = opt.value(),
    Rcpp::Named("iterations") = opt.grcount(),
    Rcpp::Named("function_evaluations") = opt.fncount(),
    Rcpp::Named("gradient_evaluations") = opt.grcount(),
    Rcpp::Named("gradient_norm") = arma::norm(final_gradient_vec, 2),
    Rcpp::Named("initial_W") = initial_W,
    Rcpp::Named("initial_value") = initial_value,
    Rcpp::Named("initial_gradient") = arma::reshape(initial_gradient_vec, X.n_cols, ncomp)
  );
}
