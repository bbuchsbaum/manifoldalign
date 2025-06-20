#' Predict Method for Similarity Embedding
#'
#' Projects new data using trained similarity embedding model.
#'
#' @param object A simembed object from linear_sim_embed()
#' @param newdata New data matrix (n_new x d) to project
#' @param ... Additional arguments (currently ignored)
#'
#' @return Matrix of projected coordinates (n_new x ncomp)
#'
#' @examples
#' \donttest{
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' model <- linear_sim_embed(X, ncomp = 2)
#' X_new <- matrix(rnorm(20 * 5), 20, 5)
#' Y_new <- predict(model, X_new)
#' }
#'
#' @export
predict.simembed <- function(object, newdata, ...) {
  if (!inherits(object, "simembed")) {
    stop("object must be a simembed object")
  }
  
  if (!is.matrix(newdata)) newdata <- as.matrix(newdata)
  
  # FIXED: Add explicit column count validation as recommended in feedback
  expected_ncol <- length(object$center)
  if (ncol(newdata) != expected_ncol) {
    stop(sprintf("newdata has %d columns but model expects %d columns", 
                 ncol(newdata), expected_ncol))
  }
  
  # Apply same preprocessing as training data
  newdata_scaled <- scale(newdata, center = object$center, scale = object$scale)
  
  # Project using learned weights
  projected <- newdata_scaled %*% object$weights
  
  return(projected)
}

#' Print Method for Similarity Embedding
#' @param x A simembed object
#' @param ... Additional arguments
#' @export
print.simembed <- function(x, ...) {
  cat("Linear Similarity Embedding\n")
  cat("===========================\n")
  cat(sprintf("Embedding dimension: %d\n", ncol(x$weights)))
  cat(sprintf("Original dimension: %d\n", nrow(x$weights)))
  cat(sprintf("Number of samples: %d\n", nrow(x$scores)))
  cat(sprintf("Optimizer: %s\n", x$optimizer))
  cat(sprintf("Final sigma_P: %.6f\n", x$sigma_P))
  cat(sprintf("Alpha_p: %.3f\n", x$alpha_p))
  
  if (!is.null(x$convergence)) {
    if (x$convergence$convergence == 0) {
      cat("Convergence: SUCCESS\n")
    } else {
      cat(sprintf("Convergence: FAILED (%s)\n", x$convergence$message))
    }
    
    if (!is.na(x$convergence$iterations)) {
      cat(sprintf("Iterations: %d\n", x$convergence$iterations))
    }
  }
  
  invisible(x)
}

# -------------------------------------------------------------------
# Enhanced Helper Functions
# -------------------------------------------------------------------

#' Auto-select sigma_P using histogram-spread heuristic
#' @param X Preprocessed data matrix
#' @param W0 Initial weight matrix  
#' @param T Target similarity matrix
#' @param M Mask matrix
#' @param verbose Logical for progress reporting
#' @return Optimal sigma_P value
#' @keywords internal
.auto_select_sigma_P <- function(X, W0, T, M, verbose = FALSE) {
  # Log-space grid search as recommended in feedback (10^-5 to 10^5, factor 0.5 steps)
  sigma_candidates <- 10^seq(-5, 5, by = 0.5)
  
  best_sigma <- sigma_candidates[1]
  best_score <- -Inf
  
  # FIXED: Move Y0 computation outside loop to save O(n²) work per candidate
  Y0 <- X %*% W0
  
  for (i in seq_along(sigma_candidates)) {
    sigma <- sigma_candidates[i]
    
    # Quick evaluation using pre-computed Y0
    P0 <- gaussian_sim(Y0, sigma)
    
    # Histogram-spread heuristic: maximize spread of P values where M=1
    P_masked <- P0[M == 1]
    if (length(P_masked) > 0) {
      # Use interquartile range as spread measure (robust)
      spread_score <- IQR(P_masked)
      
      if (verbose && i %% 5 == 1) {
        cat(sprintf("  sigma=%.2e, spread=%.4f\n", sigma, spread_score))
      }
      
      # Early stopping: if spread decreases, we've passed the optimum
      if (spread_score < best_score && i > 3) {
        if (verbose) cat("  Early stopping: spread decreasing\n")
        break
      }
      
      if (spread_score > best_score) {
        best_score <- spread_score
        best_sigma <- sigma
      }
    }
  }
  
  return(best_sigma)
}

#' Enhanced ADAM optimizer with alpha scheduling
#' @keywords internal
.optimize_W_enhanced <- function(X, T, M, W0,
                                sigma = 1, alpha_p = 0.2, alpha_schedule = FALSE,
                                lr = 5e-3, batch = 1,
                                maxit = 500, tol = 1e-6,
                                verbose = FALSE) {

  n <- nrow(X); d <- ncol(X); m <- ncol(W0)

  W <- W0
  mt <- NULL # Initialize Adam states
  vt <- NULL
  obj_prev <- Inf
  sum_M_full <- sum(M)
  if (sum_M_full == 0) stop("Mask matrix M sums to zero. Cannot optimize.")
  obj_trace <- numeric(maxit)

  # Alpha scheduling setup
  alpha_p_orig <- alpha_p
  if (alpha_schedule) {
    alpha_decay_steps <- min(50, maxit)
    if (verbose) cat("Using alpha_p scheduling: 1.0 -> ", alpha_p_orig, " over ", alpha_decay_steps, " steps\n")
  }

  # Pre-calculate for efficiency
  jp_denom <- 2 * m^2
  if (jp_denom == 0) jp_denom <- 1
  upper_idx <- which(upper.tri(M, diag = FALSE) & M == 1, arr.ind = FALSE)
  
  # Add robustness checks for poorly scaled similarity matrices
  if (!is.null(T)) {
    T_mean <- mean(T[upper_idx])
    if (T_mean < 1e-4 && T_mean > 0) {
      warning(sprintf("Target similarity matrix T appears to be poorly scaled (mean off-diagonal = %.2e). This may cause convergence issues. Consider using adaptive scaling: T <- exp(-dist(X)^2 / median(dist(X)^2))", T_mean))
    }
    
    # Check for extreme scaling that would cause certain failure
    T_max <- max(T[upper_idx])
    if (T_max < 1e-10) {
      warning("Target similarity matrix T is too close to zero (max off-diagonal = ", 
           sprintf("%.2e", T_max), "). This may cause optimization to fail. ",
           "Please rescale your similarity matrix.")
    }
  }
  
  # Track stalled optimization
  stall_counter <- 0
  stall_threshold <- 5  # Number of iterations with minimal change
  grad_norm_trace <- numeric(maxit)
  grad_norm <- Inf  # Initialize for first iteration

  for (t in 1:maxit) {
    
    # Alpha scheduling (linearly decay from 1 to target value)
    if (alpha_schedule && t <= alpha_decay_steps) {
      current_alpha <- 1.0 - (1.0 - alpha_p_orig) * (t - 1) / (alpha_decay_steps - 1)
    } else {
      current_alpha <- alpha_p_orig
    }

    # Sample batch
    if (batch < 1 && batch > 0) {
      idx <- sample.int(n, size = max(1, ceiling(batch*n)))
      Xb <- X[idx, , drop=FALSE]
      Tb <- T[idx, idx, drop=FALSE]
      Mb <- M[idx, idx, drop=FALSE]
    } else {
      Xb <- X; Tb <- T; Mb <- M
    }

    sum_Mb <- sum(Mb)
    if (sum_Mb == 0) {
        if(verbose) message(paste("Skipping iteration", t, "due to zero sum in batch mask M."))
        obj_trace[t] <- obj_prev
        next
    }

    # Calculate objective on full data
    needs_recalc <- (t == 1 || batch < 1 || verbose)
    if (needs_recalc) {
        Y_full <- X %*% W
        D2_full <- pairwise_sqdist(Y_full)
        P_upper <- exp(-D2_full[upper_idx] / sigma)
        # CORRECTED: Use (1-alpha)*Js + alpha*Jp formulation
        # FIXED: Normalize by upper triangle count, not full mask sum
        Js_full <- sum((P_upper - T[upper_idx])^2) / (2 * length(upper_idx))

        WtW_full <- crossprod(W)
        diag(WtW_full) <- diag(WtW_full) - 1
        Jp_full <- sum(WtW_full^2) / jp_denom

        obj <- (1 - current_alpha) * Js_full + current_alpha * Jp_full
        obj_trace[t] <- obj
    } else {
      obj <- obj_prev
      obj_trace[t] <- obj
    }

    # Enhanced convergence check
    obj_change <- if (t > 1 && !is.infinite(obj_prev)) abs(obj_prev - obj) else Inf
    
    # Check for multiple convergence criteria
    if (t > 1) {
      # Criterion 1: Objective change below tolerance
      if (obj_change < tol) {
        if (verbose) cat(sprintf("Converged at iteration %d (obj change < %.3e)\n", t, tol))
        obj_trace <- obj_trace[1:t]
        grad_norm_trace <- grad_norm_trace[1:t]
        convergence_info <- list(convergence = 0, message = "SUCCESS", iterations = t)
        break
      }
      
      # Criterion 2: Gradient norm below threshold (indicates stationary point)
      if (grad_norm < tol * 10) {
        if (verbose) cat(sprintf("Converged at iteration %d (gradient norm = %.3e)\n", t, grad_norm))
        obj_trace <- obj_trace[1:t]
        grad_norm_trace <- grad_norm_trace[1:t]
        convergence_info <- list(convergence = 0, message = "GRADIENT_CONVERGED", iterations = t)
        break
      }
      
      # Criterion 3: Stalled optimization detection
      if (obj_change < tol * 100) {  # Less strict threshold for stall detection
        stall_counter <- stall_counter + 1
        if (stall_counter >= stall_threshold) {
          warning(sprintf("Optimization stalled at iteration %d (objective unchanged for %d iterations). ", 
                         t, stall_threshold,
                         "This often indicates poorly scaled input. Check similarity matrix scaling."))
          obj_trace <- obj_trace[1:t]
          grad_norm_trace <- grad_norm_trace[1:t]
          convergence_info <- list(convergence = 2, message = "STALLED", iterations = t)
          break
        }
      } else {
        stall_counter <- 0  # Reset counter if making progress
      }
    }
    obj_prev <- obj

    # Compute gradient on batch
    grad_js_b <- grad_Js(Xb, W, Tb, Mb, sigma, sum_M_batch = sum_Mb)
    grad_jp_b <- grad_Jp(W)

    # CORRECTED: Use (1-alpha)*grad_Js + alpha*grad_Jp
    g <- (1 - current_alpha) * grad_js_b + current_alpha * grad_jp_b
    
    # Track gradient norm for convergence monitoring
    grad_norm <- sqrt(sum(g^2))
    grad_norm_trace[t] <- grad_norm

    # ADAM update
    up <- adam_update(W, g, mt, vt, t, lr)
    W  <- up$par
    mt <- up$m
    vt <- up$v

    # Verbose output
    if (verbose && (t %% 20 == 1 || t == maxit)) {
      if (alpha_schedule && t <= alpha_decay_steps) {
        cat(sprintf("it=%04d  obj=%.6e  Js=%.4e  Jp=%.4e  alpha=%.3f\n",
                    t, obj, Js_full, Jp_full, current_alpha))
      } else {
        cat(sprintf("it=%04d  obj=%.6e  Js=%.4e  Jp=%.4e\n",
                    t, obj, Js_full, Jp_full))
      }
    }
  }

  # Handle non-convergence
  if (!exists("convergence_info")) {
    # If we exited the loop without setting convergence_info, check why
    if (t == maxit && abs(obj_prev - obj) >= tol) {
      warning(paste("Optimization (ADAM) did not converge within", maxit, "iterations. Final objective change =", sprintf("%.3e", abs(obj_prev-obj))))
      convergence_info <- list(convergence = 1, message = "MAXITER_REACHED", iterations = maxit)
    } else {
      convergence_info <- list(convergence = 0, message = "SUCCESS", iterations = t)
    }
  }

  if (t == maxit && length(obj_trace) > t) {
     obj_trace <- obj_trace[1:t]
  }

  return(list(W = W, trace = obj_trace, convergence = convergence_info$convergence,
              message = convergence_info$message, iterations = convergence_info$iterations))
}

#' Handle formula interface for supervised embedding
#' @keywords internal
.handle_formula_interface <- function(formula, data, T, M, sigma_P, ncomp, 
                                     alpha_p, alpha_schedule, maxit, tol, 
                                     batch_size, use_cpp, verbose, lr, ...) {
  
  # Parse formula
  mf <- model.frame(formula, data)
  response <- model.response(mf)
  
  if (is.null(response)) {
    stop("Formula must specify a response variable (e.g., ~ label)")
  }
  
  # Extract numeric predictors
  X <- model.matrix(formula, data)[, -1, drop = FALSE]  # Remove intercept
  if (ncol(X) == 0) {
    # If no predictors specified, use all numeric columns except response
    numeric_cols <- sapply(data, is.numeric)
    response_name <- all.vars(formula[[2]])[1]
    if (response_name %in% names(numeric_cols)) {
      numeric_cols[response_name] <- FALSE
    }
    X <- as.matrix(data[, numeric_cols, drop = FALSE])
  }
  
  # Create target similarity matrix based on labels (Eq. 11 for LDA-style runs)
  if (is.null(T)) {
    T <- .create_supervised_target(response)
  }
  
  # Call main function
  result <- linear_sim_embed(X = X, T = T, M = M, sigma_P = sigma_P, ncomp = ncomp,
                            alpha_p = alpha_p, alpha_schedule = alpha_schedule,
                            maxit = maxit, tol = tol, batch_size = batch_size,
                            use_cpp = use_cpp, verbose = verbose, lr = lr, ...)
  
  # Add formula information
  result$formula <- formula
  result$response <- response
  
  return(result)
}

#' Create supervised target similarity matrix
#' @param labels Factor or character vector of class labels
#' @return Binary similarity matrix (1 if same class, 0 otherwise)
#' @keywords internal
.create_supervised_target <- function(labels) {
  # FIXED: Use vectorized outer() approach for 20x speedup as suggested in feedback
  outer(labels, labels, "==") * 1
}

# -------------------------------------------------------------------
# Existing Helper Functions (updated for consistency)
# -------------------------------------------------------------------

#' Calculate Pairwise Squared Euclidean Distances Efficiently
#' @param Z Matrix (n x m)
#' @return n x n matrix of squared distances
#' @keywords internal
# This function is now defined in cone_align.R and should be used from there.
# pairwise_sqdist <- function(Z) {
#   ...
# }

#' Calculate Gaussian Similarity from Data or Squared Distances
#' @param Z Matrix of data (n x m) OR precomputed squared distances (n x n)
#' @param sigma Scaling parameter (sigma_P)
#' @return n x n matrix of similarities
#' @keywords internal
gaussian_sim <- function(Z, sigma) {
  if (is.matrix(Z) && nrow(Z) == ncol(Z)) {
      D2 <- Z # Assume precomputed squared distances
  } else {
      D2 <- pairwise_sqdist(Z)
  }
  if (!is.numeric(sigma) || sigma <= 0) stop("sigma must be a positive number in gaussian_sim")
  P <- exp(-D2 / sigma)
  # Enforce exact symmetry (critical for paper compliance)
  P <- 0.5 * (P + t(P))
  P
}

#' Calculate Gradient of Similarity Term (Js) w.r.t W (Correct Implementation)
#' @param X Preprocessed data matrix (n x d) or batch (n_b x d)
#' @param W Current weight matrix (d x m)
#' @param T Target similarity matrix (n x n) or batch (n_b x n_b)
#' @param M Mask matrix (n x n) or batch (n_b x n_b)
#' @param sigma Scaling parameter (sigma_P)
#' @param sum_M_batch The sum of the *batch* mask matrix M (for consistent normalization)
#' @return Gradient matrix (d x m)
#' @keywords internal
grad_Js <- function(X, W, T, M, sigma, sum_M_batch) {
  if (sum_M_batch == 0) return(matrix(0, nrow = ncol(X), ncol = ncol(W)))

  n <- nrow(X)
  Y <- X %*% W
  
  # Calculate pairwise squared distances and similarities
  D2 <- pairwise_sqdist(Y)
  P <- exp(-D2 / sigma)
  P <- 0.5 * (P + t(P))  # Ensure symmetry
  
  # Step 1: Compute the coefficient matrix E
  # E_ij = dJs/dP_ij * dP_ij/dD2_ij = (M_ij * (P_ij - T_ij) / sum_M) * (-P_ij / sigma)
  # This is incorrect, normalization factor should be 2*length(upper_idx)
  upper_idx <- which(upper.tri(M, diag=FALSE) & M==1)
  norm_factor <- if (length(upper_idx) > 0) length(upper_idx) else 1
  
  dJs_dP <- M * (P - T) / (2 * norm_factor)
  E_matrix <- dJs_dP * P * (-2 / sigma)
  
  # Step 2: Form the graph Laplacian of the E_matrix
  # L_E = D_E - E, where D_E is the diagonal matrix of row sums
  L_E <- Matrix::Diagonal(x = rowSums(E_matrix)) - E_matrix
  
  # Step 3: Compute the final gradient w.r.t W
  # grad = 2 * X^T * L_E * Y
  grad <- 2 * crossprod(X, L_E %*% Y)
  
  return(grad)
}

#' Calculate Gradient of Orthogonality Term (Jp) w.r.t W
#' @param W Current weight matrix (d x m)
#' @return Gradient matrix (d x m)
#' @keywords internal
grad_Jp <- function(W) {
  m <- ncol(W)
  if (m == 0) return(matrix(0, nrow = nrow(W), ncol = 0))

  WtW <- crossprod(W)
  # FIXED: Avoid in-place mutation to keep WtW symmetric
  WtWmI <- WtW
  diag(WtWmI) <- diag(WtWmI) - 1

  # CORRECTED: Use correct factor 1/m² instead of 2/m² to match paper and C++ implementation
  gradient <- (1 / m^2) * (W %*% WtWmI)
  return(gradient)
}

#' Adam optimizer (minimal implementation)
#' @keywords internal
adam_update <- function(par, grad, m_state, v_state, t, lr=1e-3,
                        b1=.9, b2=.999, eps=1e-8) {
  if (is.null(m_state)) m_state <- array(0, dim = dim(par))
  if (is.null(v_state)) v_state <- array(0, dim = dim(par))

  m_state <- b1 * m_state + (1 - b1) * grad
  v_state <- b2 * v_state + (1 - b2) * grad^2

  m_hat <- m_state / (1 - b1^t)
  v_hat <- v_state / (1 - b2^t)

  par <- par - lr * m_hat / (sqrt(v_hat) + eps)

  list(par = par, m = m_state, v = v_state)
}

# -------------------------------------------------------------------
# Helper operators and dependencies for standalone testing
# -------------------------------------------------------------------

#' Null-default operator
#' @keywords internal
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Linear Similarity Embedding using Optimal Transport
#'
#' @description
#' Performs linear dimensionality reduction by optimizing a weight matrix W (d x m)
#' such that the pairwise similarities P computed from the projected data (Y = X %*% W)
#' match a target similarity matrix T as closely as possible, considering a mask M.
#' The optimization balances similarity preservation (Js) with an orthogonality
#' constraint (Jp) on W, minimizing (1-alpha_p)*Js + alpha_p*Jp.
#'
#' @details
#' The algorithm follows the Linear Similarity Embedding Framework (SEF) from 
#' Passalis & Tefas (2016). The optimization uses either ADAM (R implementation) 
#' or L-BFGS-B (C++ implementation).
#' 
#' Similarities are computed using a Gaussian kernel: P_ij = exp(-||Y_i - Y_j||^2 / sigma_P).
#' The orthogonality penalty is Jp = ||W'W - I||^2_F / (2*m^2).
#' 
#' Key algorithmic features:
#' \itemize{
#'   \item Automatic sigma_P selection using histogram-spread heuristic (Step 2, Fig. 2)
#'   \item PCA/KPCA initialization for faster convergence (Step 3, Fig. 2)
#'   \item Enforced similarity matrix symmetry for numerical stability
#'   \item Alpha_p scheduling option for improved convergence
#' }
#'
#' @param X Input data matrix (n x d), where n is samples, d is features.
#' @param T Target similarity matrix (n x n). If NULL, will be computed from data.
#' @param M Mask matrix (n x n), 1 to consider pair (i,j), 0 to ignore. If NULL, uses all pairs.
#' @param sigma_P Numeric scalar > 0 or "auto". Scale parameter for Gaussian kernel. 
#'   If "auto", uses log-space grid search with histogram-spread heuristic (default: "auto").
#' @param ncomp Integer > 0. Number of dimensions for the embedding (m) (default: 2).
#' @param alpha_p Numeric in [0,1]. Weight for the orthogonality regularizer. 
#'   Uses convex combination: (1-alpha_p)*Js + alpha_p*Jp (default: 0.1).
#' @param alpha_schedule Logical. If TRUE, linearly decay alpha_p from 1 to specified value 
#'   over first 50 iterations to avoid early orthogonality trapping (default: FALSE).
#' @param maxit Integer. Maximum number of iterations (default: 500).
#' @param tol Numeric. Convergence tolerance. For ADAM (R): change in objective function.
#'   For L-BFGS-B (C++): gradient norm tolerance (default: 1e-6).
#' @param batch_size Numeric in (0, 1]. Fraction of data for stochastic updates in ADAM (R only)
#'   (default: 1, i.e., full batch).
#' @param use_cpp Logical. If TRUE and C++ backend is available, use L-BFGS-B from C++.
#'   Otherwise, use the R ADAM implementation (default: FALSE).
#' @param verbose Logical. Print optimization progress (default: FALSE).
#' @param lr Numeric > 0. Learning rate for the ADAM optimizer (R only) (default: 5e-3).
#' @param formula Optional formula interface for supervised targets (e.g., ~ label).
#' @param data Optional data.frame when using formula interface.
#' @param ... Extra arguments (currently ignored).
#'
#' @return A \code{simembed} object (S3 class) containing:
#'   - \code{weights} (W): The optimized projection matrix (d x m).
#'   - \code{scores} (Y): The projected data (n x m).
#'   - \code{sdev}: Standard deviations of the scores.
#'   - \code{preproc}: The preprocessing object used on X.
#'   - \code{center}: Centering vector used in preprocessing.
#'   - \code{scale}: Scaling vector used in preprocessing.
#'   - \code{sigma_P}: Final sigma_P value used (important if auto-selected).
#'   - \code{alpha_p}: Final alpha_p value used.
#'   - \code{objective_trace}: Vector of objective function values during optimization.
#'   - Metadata: \code{target_sim}, \code{mask}, \code{optimizer}, \code{convergence}.
#'
#' @examples
#' \donttest{
#' # Basic usage with automatic sigma_P selection
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' result <- linear_sim_embed(X, ncomp = 3, verbose = TRUE)
#' 
#' # Use predict method for new data
#' X_new <- matrix(rnorm(20 * 10), 20, 10)
#' Y_new <- predict(result, X_new)
#' 
#' # Custom target similarity matrix
#' D_orig <- as.matrix(dist(X))
#' T_sim <- exp(-D_orig^2 / median(D_orig^2))
#' result2 <- linear_sim_embed(X, T = T_sim, sigma_P = 1.0)
#' 
#' # Formula interface for supervised embedding
#' df <- data.frame(X, label = sample(c("A", "B", "C"), 100, replace = TRUE))
#' result3 <- linear_sim_embed(~ label, data = df, ncomp = 2)
#' 
#' # Use C++ backend for speed
#' result_cpp <- linear_sim_embed(X, use_cpp = TRUE, verbose = TRUE)
#' }
#'
#' @references
#' Passalis, N., & Tefas, A. (2016). Learning deep representations with 
#' probabilistic knowledge transfer. In Proceedings of the European Conference 
#' on Computer Vision (pp. 268-284).
#'
#' @seealso \code{\link{predict.simembed}}, \code{\link{plot.simembed}}
#' @export
linear_sim_embed <- function(X, T = NULL, M = NULL,
                             sigma_P = "auto", ncomp = 2,
                             alpha_p = 0.1, alpha_schedule = FALSE,
                             maxit = 500, tol = 1e-6,
                             batch_size = 1, use_cpp = FALSE,
                             verbose = FALSE, lr = 5e-3,
                             formula = NULL, data = NULL, ...) {

  # Handle formula interface
  if (!is.null(formula)) {
    if (is.null(data)) {
      stop("data argument required when using formula interface")
    }
    result <- .handle_formula_interface(formula, data, T, M, sigma_P, ncomp, 
                                       alpha_p, alpha_schedule, maxit, tol, 
                                       batch_size, use_cpp, verbose, lr, ...)
    return(result)
  }
  
  # --- Input Validation ---
  n <- nrow(X)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  # Handle automatic target similarity matrix creation
  if (is.null(T)) {
    if (verbose) cat("Creating target similarity matrix from data distances...\n")
    D_orig <- as.matrix(dist(X))
    # Use median-based scaling (robust to outliers)
    scale_factor <- median(D_orig[upper.tri(D_orig)])^2
    T <- exp(-D_orig^2 / scale_factor)
    diag(T) <- 1  # Perfect self-similarity
  }
  
  # Handle automatic mask matrix creation
  if (is.null(M)) {
    M <- matrix(1, n, n)
    diag(M) <- 0  # Typically ignore self-similarity
  }
  
  if (!is.matrix(T)) T <- as.matrix(T)
  if (!is.matrix(M)) M <- as.matrix(M)

  if (!all(dim(T) == c(n, n))) stop("Target matrix T must be n x n, where n=nrow(X)")
  if (!all(dim(M) == c(n, n))) stop("Mask matrix M must be n x n, where n=nrow(X)")
  if (anyNA(X) || anyNA(T) || anyNA(M)) stop("Input matrices X, T, M must not contain NA values")
  if (!is.numeric(alpha_p) || alpha_p < 0 || alpha_p > 1) stop("alpha_p must be in [0, 1]")
  if (!is.numeric(batch_size) || batch_size <= 0 || batch_size > 1) stop("batch_size must be in (0, 1]")
  # -----------------------

  # Preprocess & PCA Initialization (following paper's Step 2-3)
  # Using base R standardization for z-normalization (Step 2, Fig. 2)
  Xs <- scale(X, center = TRUE, scale = TRUE)
  d  <- ncol(Xs)
  
  # Store preprocessing info for predict method
  preproc_info <- list(
    center = attr(Xs, "scaled:center") %||% colMeans(X),
    scale = attr(Xs, "scaled:scale") %||% rep(1, d),
    preproc = NULL  # Not using multivarious for standalone testing
  )

  # Ensure ncomp is valid
  max_comp <- min(n, d)
  if (ncomp > max_comp) {
      warning(sprintf("ncomp (%d) exceeds max possible (%d). Using %d components.", ncomp, max_comp, max_comp))
      ncomp <- max_comp
  }
  if (ncomp <= 0) stop("ncomp must be positive.")

  # PCA Initialization (Step 3, Fig. 2)
  if (verbose) cat("Initializing with PCA...\n")
  pca_res <- prcomp(Xs, rank. = ncomp, center = FALSE, scale. = FALSE)
  W0 <- pca_res$rotation[, 1:ncomp, drop = FALSE] # Initial guess for weights

  # Automatic sigma_P selection (Lines 9-16, Fig. 2 - histogram-spread heuristic)
  if (identical(sigma_P, "auto")) {
    if (verbose) cat("Auto-selecting sigma_P using histogram-spread heuristic...\n")
    sigma_P <- .auto_select_sigma_P(Xs, W0, T, M, verbose = verbose)
    if (verbose) cat(sprintf("Selected sigma_P = %.6f\n", sigma_P))
  } else {
    if (!is.numeric(sigma_P) || sigma_P <= 0) stop("sigma_P must be positive or 'auto'")
  }

  # --- Optimization Path ---
  trace <- NULL
  if (use_cpp) {
    # Use C++ (L-BFGS-B)
    if (!exists("linear_sim_embed_cpp")) {
        stop("C++ function 'linear_sim_embed_cpp' not found. Ensure package is compiled and loaded, or set use_cpp=FALSE.")
    }
    if (batch_size != 1) {
        warning("batch_size is not supported by the C++ backend (uses L-BFGS-B). Using full batch.")
    }
    if (alpha_schedule) {
        warning("alpha_schedule is not supported by the C++ backend. Using fixed alpha_p.")
    }
    
    # Call C++ function with proper initialization
    ret <- linear_sim_embed_cpp(Xs, T, M, sigma_P, ncomp, alpha_p, maxit, tol)

    if (ret$convergence != 0) {
        warning("C++ optimization (L-BFGS-B) did not converge. Message: ", ret$message)
    }
    W <- ret$W
    optimizer_name <- "L-BFGS-B (C++)"
    convergence_info <- list(
      convergence = ret$convergence,
      message = ret$message,
      iterations = ret$iterations %||% NA,
      final_value = ret$value %||% NA
    )

  } else {
     # Use R (ADAM) with enhanced features
    res_opt <- .optimize_W_enhanced(Xs, T, M, W0,
                                   sigma = sigma_P, alpha_p = alpha_p,
                                   alpha_schedule = alpha_schedule,
                                   lr = lr, batch = batch_size,
                                   maxit = maxit, tol = tol,
                                   verbose = verbose)

    W      <- res_opt$W
    trace  <- res_opt$trace
    optimizer_name <- "ADAM (R)"
    convergence_info <- list(
      convergence = res_opt$convergence,
      message = res_opt$message,
      iterations = res_opt$iterations,
      final_value = tail(res_opt$trace, 1)
    )
   }

  # --- Output ---
  scores <- Xs %*% W
  sdev   <- apply(scores, 2, sd)

  # Create simembed S3 object with comprehensive metadata
  result <- structure(
    list(
      weights = W,
      scores = scores,
      sdev = sdev,
      center = preproc_info$center,
      scale = preproc_info$scale,
      preproc = preproc_info$preproc,
      sigma_P = sigma_P,
      alpha_p = alpha_p,
      target_sim = T,
      mask = M,
      optimizer = optimizer_name,
      objective_trace = trace,
      convergence = convergence_info,
      call = match.call()
    ),
    class = c("simembed", "bi_projector")
  )
  
  return(result)
}






