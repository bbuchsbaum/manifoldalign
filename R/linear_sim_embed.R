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
#'
#' @param x A simembed object
#' @param ... Additional arguments
#'
#' @return The simembed object x, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(100), 20, 5)
#' model <- linear_sim_embed(X, ncomp = 2)
#' print(model)
#' }
#'
#' @method print simembed
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
#' @return Numeric scalar giving the optimal sigma_P value
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
#' @return List with components W (optimized weight matrix), trace (objective values),
#'   convergence (status code), message (convergence message), and iterations (count)
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
  upper_idx <- which(upper.tri(M, diag = FALSE) & M > 0, arr.ind = FALSE)
  
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
        P_full <- exp(-D2_full / sigma)
        P_full <- 0.5 * (P_full + t(P_full))
        Js_full <- sum(M * (P_full - T)^2) / (2 * sum_M_full)

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
          warning(sprintf(
            "Optimization stalled at iteration %d (objective unchanged for %d iterations). This often indicates poorly scaled input. Check similarity matrix scaling.",
            t, stall_threshold
          ))
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
    # Exhausting the loop is never convergence: the last update has not yet
    # been checked against either the objective or gradient criterion.
    warning(paste("Optimization (ADAM) did not converge within", maxit, "iterations."))
    convergence_info <- list(convergence = 1, message = "MAXITER_REACHED", iterations = maxit)
  }

  if (t == maxit && length(obj_trace) > t) {
     obj_trace <- obj_trace[1:t]
  }

  return(list(W = W, trace = obj_trace, convergence = convergence_info$convergence,
              message = convergence_info$message, iterations = convergence_info$iterations))
}

#' Handle formula interface for supervised embedding
#' @return A simembed object with formula and response components added
#' @keywords internal
.handle_formula_interface <- function(formula, data, T, M, sigma_P, ncomp, 
                                     alpha_p, alpha_schedule, maxit, tol, 
                                     batch_size, use_cpp, verbose, lr, ...) {
  
  if (!inherits(formula, "formula")) {
    stop("formula must be a formula object", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("data must be a data.frame", call. = FALSE)
  }

  # Support simple response-only form: ~ label
  response <- NULL
  if (length(formula) == 2L) {
    response_name <- all.vars(formula[[2]])
    if (length(response_name) != 1L) {
      stop("Formula must be of the form ~ label or label ~ predictors", call. = FALSE)
    }
    response_name <- response_name[[1]]
    if (!response_name %in% names(data)) {
      stop("Response variable '", response_name, "' not found in data", call. = FALSE)
    }
    response <- data[[response_name]]
    numeric_cols <- vapply(data, is.numeric, logical(1))
    numeric_cols[response_name] <- FALSE
    X <- as.matrix(data[, numeric_cols, drop = FALSE])
  } else {
    # Standard response ~ predictors
    mf <- model.frame(formula, data)
    response <- model.response(mf)
    if (is.null(response)) {
      stop("Formula must specify a response variable", call. = FALSE)
    }
    X <- model.matrix(formula, data)
    if (ncol(X) > 0 && identical(colnames(X)[1], "(Intercept)")) {
      X <- X[, -1, drop = FALSE]
    }
  }

  if (!is.matrix(X) || ncol(X) == 0) {
    stop("No numeric predictors found for embedding", call. = FALSE)
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
  
  # With a symmetric mask, both ordered entries (i, j) and (j, i)
  # contribute to dD^2/dW. This factor of four is part of the derivative,
  # not an optimizer tuning constant.
  A <- M * (P - T) * P
  L_A <- Matrix::Diagonal(x = rowSums(A)) - A
  grad <- as.matrix((-4 / (sigma * sum_M_batch)) * crossprod(X, L_A %*% Y))
  
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

  gradient <- (2 / m^2) * (W %*% WtWmI)
  return(gradient)
}

.lse_objective_gradient <- function(W, X, T, M, sigma_P, alpha_p) {
  Y <- X %*% W
  D2 <- pairwise_sqdist(Y)
  P <- exp(-D2 / sigma_P)
  P <- 0.5 * (P + t(P))

  Js <- sum(M * (P - T)^2) / (2 * sum(M))
  WtWmI <- crossprod(W)
  diag(WtWmI) <- diag(WtWmI) - 1
  Jp <- sum(WtWmI^2) / (2 * ncol(W)^2)

  gradient <- (1 - alpha_p) * grad_Js(
    X, W, T, M, sigma_P, sum_M_batch = sum(M)
  ) + alpha_p * grad_Jp(W)

  list(
    value = (1 - alpha_p) * Js + alpha_p * Jp,
    gradient = as.matrix(gradient),
    similarity_loss = Js,
    orthogonality_loss = Jp
  )
}

.optimize_W_reference <- function(X, T, M, W0, sigma_P, alpha_p, maxit, tol) {
  dims <- dim(W0)
  evaluate <- function(parameters) {
    W <- matrix(parameters, nrow = dims[[1L]], ncol = dims[[2L]])
    .lse_objective_gradient(W, X, T, M, sigma_P, alpha_p)
  }
  fit <- stats::optim(
    par = as.numeric(W0),
    fn = function(parameters) evaluate(parameters)$value,
    gr = function(parameters) as.numeric(evaluate(parameters)$gradient),
    method = "L-BFGS-B",
    control = list(maxit = maxit, pgtol = tol)
  )

  list(
    W = matrix(fit$par, nrow = dims[[1L]], ncol = dims[[2L]]),
    trace = fit$value,
    convergence = as.integer(fit$convergence),
    message = fit$message %||% .lse_cpp_status(fit$convergence),
    iterations = unname(fit$counts[["gradient"]]),
    final_value = fit$value,
    function_evaluations = unname(fit$counts[["function"]]),
    gradient_evaluations = unname(fit$counts[["gradient"]])
  )
}

#' Adam optimizer update step
#'
#' @param par Current parameter matrix
#' @param grad Gradient matrix
#' @param m_state First moment state (NULL on first call)
#' @param v_state Second moment state (NULL on first call)
#' @param t Iteration number (1-indexed)
#' @param lr Learning rate
#' @param b1 Exponential decay rate for first moment
#' @param b2 Exponential decay rate for second moment
#' @param eps Small constant for numerical stability
#' @return A list with updated parameters and optimizer state
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
#' @param x First value
#' @param y Default value if x is NULL
#' @return The value of x if not NULL, otherwise y
#' @name grapes-or-or-grapes
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

.lse_cpp_status <- function(code) {
  switch(
    as.character(as.integer(code)),
    `0` = "converged",
    `1` = "iteration_limit",
    `51` = "optimizer_warning",
    `52` = "optimizer_error",
    "unknown_failure"
  )
}

#' Linear Similarity Embedding using Optimal Transport
#'
#' @description
#' Performs linear dimensionality reduction by optimizing a weight matrix W (d x m)
#' such that the pairwise similarities P computed from the projected data
#' (`Y = X %*% W`)
#' match a target similarity matrix T as closely as possible, considering a mask M.
#' The optimization balances similarity preservation (Js) with an orthogonality
#' constraint (Jp) on W, minimizing (1-alpha_p)*Js + alpha_p*Jp.
#'
#' @details
#' The algorithm follows the Linear Similarity Embedding Framework (SEF) from 
#' Passalis & Tefas (2016). The optimization uses either ADAM (R implementation)
#' or L-BFGS-B (C++ implementation).
#'
#' Compiled L-BFGS-B status codes are reported as `converged` (0),
#' `iteration_limit` (1), `optimizer_warning` (51), `optimizer_error` (52), or
#' `unknown_failure`. A nonzero compiled status triggers an exact-objective
#' base-R L-BFGS-B retry. The function returns a fit only if one backend
#' converges; otherwise it stops instead of returning unchecked weights.
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
#' @param T Symmetric target similarity matrix (n x n). If NULL, will be
#'   computed from data.
#' @param M Symmetric non-negative mask matrix (n x n), with positive entries
#'   selecting or weighting pairs. If NULL, uses all off-diagonal pairs.
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
#'   A failed compiled fit is retried with exact-objective base-R L-BFGS-B.
#'   If FALSE, use the R ADAM implementation (default: FALSE).
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
#'   - Metadata: \code{target_sim}, \code{mask}, \code{optimizer}, and
#'     \code{convergence}; the latter records normalized status, backend,
#'     objective, evaluation counts when available, and any fallback source.
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
#' @seealso \code{\link{predict.simembed}}
#' @export
linear_sim_embed <- function(X, T = NULL, M = NULL,
                             sigma_P = "auto", ncomp = 2,
                             alpha_p = 0.1, alpha_schedule = FALSE,
                             maxit = 500, tol = 1e-6,
                             batch_size = 1, use_cpp = FALSE,
                             verbose = FALSE, lr = 5e-3,
                             formula = NULL, data = NULL, ...) {

  # Allow formula as the first argument (typical R formula interface)
  if (inherits(X, "formula") && is.null(formula)) {
    formula <- X
    X <- NULL
  }

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
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  
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
  if (!is.numeric(X) || !is.numeric(T) || !(is.numeric(M) || is.logical(M))) {
    stop("X, T, and M must be numeric matrices")
  }
  storage.mode(M) <- "double"
  if (any(!is.finite(X)) || any(!is.finite(T)) || any(!is.finite(M))) {
    stop("Input matrices X, T, M must contain only finite values")
  }
  if (!isTRUE(isSymmetric(T, tol = 1e-10))) {
    stop("Target matrix T must be symmetric")
  }
  if (!isTRUE(isSymmetric(M, tol = 1e-10))) {
    stop("Mask matrix M must be symmetric")
  }
  if (any(M < 0)) stop("Mask matrix M must contain non-negative weights")
  M_offdiag <- M
  diag(M_offdiag) <- 0
  if (sum(M_offdiag) <= 0) stop("Mask matrix M sums to zero off the diagonal. Cannot optimize.")
  if (!is.numeric(alpha_p) || length(alpha_p) != 1L || !is.finite(alpha_p) || alpha_p < 0 || alpha_p > 1) {
    stop("alpha_p must be in [0, 1]")
  }
  if (!is.numeric(batch_size) || length(batch_size) != 1L || !is.finite(batch_size) || batch_size <= 0 || batch_size > 1) {
    stop("batch_size must be in (0, 1]")
  }
  if (!is.numeric(maxit) || length(maxit) != 1L || !is.finite(maxit) || maxit < 1 || maxit != as.integer(maxit)) {
    stop("maxit must be a positive integer")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("tol must be finite and positive")
  }
  if (!is.numeric(lr) || length(lr) != 1L || !is.finite(lr) || lr <= 0) {
    stop("lr must be finite and positive")
  }
  if (!is.logical(use_cpp) || length(use_cpp) != 1L || is.na(use_cpp)) {
    stop("use_cpp must be TRUE or FALSE")
  }
  if (!is.logical(alpha_schedule) || length(alpha_schedule) != 1L || is.na(alpha_schedule)) {
    stop("alpha_schedule must be TRUE or FALSE")
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE")
  }
  # -----------------------

  # Preprocess & PCA Initialization (following paper's Step 2-3)
  # Using base R standardization for z-normalization (Step 2, Fig. 2)
  Xs <- scale(X, center = TRUE, scale = TRUE)
  if (any(!is.finite(Xs))) {
    stop("X contains a constant column; remove constant features before fitting.")
  }
  d  <- ncol(Xs)
  
  # Store preprocessing info for predict method
  preproc_info <- list(
    center = attr(Xs, "scaled:center") %||% colMeans(X),
    scale = attr(Xs, "scaled:scale") %||% rep(1, d),
    preproc = NULL  # Not using multivarious for standalone testing
  )

  # Ensure ncomp is valid
  if (!is.numeric(ncomp) || length(ncomp) != 1L || !is.finite(ncomp) ||
      ncomp < 1 || ncomp != as.integer(ncomp)) {
    stop("ncomp must be a positive integer.")
  }
  ncomp <- as.integer(ncomp)
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
    if (!is.numeric(sigma_P) || length(sigma_P) != 1L || !is.finite(sigma_P) || sigma_P <= 0) {
      stop("sigma_P must be finite and positive or 'auto'")
    }
  }

  # --- Optimization Path ---
  trace <- NULL
  cpp_failure <- NULL
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
    
    # Call C++ function with proper initialization, fallback to R path on failure
    ret <- tryCatch(
      linear_sim_embed_cpp(Xs, T, M, sigma_P, ncomp, alpha_p, maxit, tol),
      error = function(e) {
        cpp_failure <<- list(
          convergence = NA_integer_,
          status = "backend_error",
          message = conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(ret) && identical(as.integer(ret$convergence), 0L)) {
      W <- ret$W
      optimizer_name <- "L-BFGS-B (C++)"
      convergence_info <- list(
        convergence = 0L,
        status = .lse_cpp_status(ret$convergence),
        message = ret$message,
        iterations = ret$iterations %||% NA,
        final_value = ret$value %||% NA,
        gradient_norm = ret$gradient_norm %||% NA,
        function_evaluations = ret$function_evaluations %||% NA,
        gradient_evaluations = ret$gradient_evaluations %||% NA,
        backend = "cpp"
      )
      use_cpp <- TRUE
    } else {
      if (!is.null(ret)) {
        cpp_failure <- list(
          convergence = as.integer(ret$convergence),
          status = .lse_cpp_status(ret$convergence),
          message = ret$message
        )
      }
      warning(
        "C++ optimization (L-BFGS-B) did not converge [",
        cpp_failure$status,
        if (!is.na(cpp_failure$convergence)) paste0("; code ", cpp_failure$convergence) else "",
        "]: ", cpp_failure$message,
        ". Retrying with the R optimizer using the same objective."
      )
      use_cpp <- FALSE
    }

  }

  if (!use_cpp) {
    if (is.null(cpp_failure)) {
      # User-selected R path: retain the stochastic/full-batch ADAM backend.
      res_opt <- .optimize_W_enhanced(
        Xs, T, M, W0,
        sigma = sigma_P, alpha_p = alpha_p,
        alpha_schedule = alpha_schedule,
        lr = lr, batch = batch_size,
        maxit = maxit, tol = tol,
        verbose = verbose
      )
    } else {
      # A failed compiled fit is retried only with an optimizer that evaluates
      # the exact same objective and analytic gradient.
      res_opt <- .optimize_W_reference(
        Xs, T, M, W0,
        sigma_P = sigma_P,
        alpha_p = alpha_p,
        maxit = maxit,
        tol = tol
      )
    }

    if (!identical(as.integer(res_opt$convergence), 0L)) {
      fallback_context <- if (is.null(cpp_failure)) {
        ""
      } else {
        paste0(
          " after C++ ", cpp_failure$status,
          if (!is.na(cpp_failure$convergence)) paste0(" (code ", cpp_failure$convergence, ")") else ""
        )
      }
      stop(
        "Linear Similarity Embedding failed to converge", fallback_context,
        ": ", res_opt$message,
        call. = FALSE
      )
    }

    W      <- res_opt$W
    trace  <- res_opt$trace
    optimizer_name <- if (is.null(cpp_failure)) "ADAM (R)" else "L-BFGS-B (R fallback)"
    convergence_info <- list(
      convergence = 0L,
      status = if (is.null(cpp_failure)) "converged" else "fallback_converged",
      message = res_opt$message,
      iterations = res_opt$iterations,
      final_value = tail(res_opt$trace, 1),
      backend = "r",
      fallback_from = cpp_failure
    )
   }

  if (any(!is.finite(W))) {
    stop("Linear Similarity Embedding produced non-finite weights.", call. = FALSE)
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
