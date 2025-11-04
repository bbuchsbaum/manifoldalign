#' Coupled Diagonalization for Hyperdesign Objects
#'
#' Performs Coupled Diagonalization alignment on hyperdesign data
#' structures containing multiple modalities. Finds coupled eigenvectors
#' that jointly diagonalize graph
#' Laplacians while respecting cross-modal correspondences.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param correspondence Optional list of correspondence matrices
#'   between modalities. If NULL (default), assumes full correspondence
#'   based on sample order.
#'   Each matrix should be n_i x n_j sparse matrix indicating correspondences 
#'   between modality i and j.
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of coupled eigenvectors to extract (default: 10)
#' @param ncomp_per_domain Number of eigenvectors per domain to compute before 
#'   coupling (default: 20). Must be >= ncomp.
#' @param mu_coupling Coupling weight controlling alignment vs diagonalization 
#'   trade-off (default: 1). Higher values emphasize cross-modal alignment.
#' @param alpha_match Optional weight for matching the diagonal of 
#'   \code{A_i^T Lambda_i A_i} to target eigenvalues (default: 0).
#' @param lambda_target Optional list of length \code{m} providing target eigenvalue
#'   vectors for each modality (only the first \code{ncomp} entries are used).
#' @param knn Number of nearest neighbors for graph construction (default: 10)
#' @param sigma Gaussian kernel bandwidth for graph weights. If NULL (default), 
#'   uses median of nearest neighbor distances.
#' @param max_iter Maximum iterations for Stiefel manifold optimization
#'   (default: 200)
#' @param step_size Base step size for gradient descent (default: 0.3)
#' @param tol Convergence tolerance for relative cost change (default: 1e-6)
#' @param verbose Whether to print convergence information (default: FALSE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing:
#' \itemize{
#'   \item \code{s}: Coupled eigenvectors (scores) concatenated across domains
#'   \item \code{v}: Projection matrices for out-of-sample extension  
#'   \item \code{coupled_bases}: List of coupled bases V_i for each modality
#'   \item \code{sdev}: Standard deviations of coupled components
#'   \item \code{eigenvalues}: List of eigenvalues for each domain
#'   \item \code{diagnostics}: List with per-iteration cost diagnostics
#'   \item \code{converged}: Logical indicating convergence
#'   \item \code{iterations}: Number of iterations performed
#'   \item \code{final_cost}: Final objective value
#' }
#'
#' @details
#' The algorithm minimizes the objective
#' \preformatted{
#' F = sum_i ||A_i^T Lambda_i A_i - diag(A_i^T Lambda_i A_i)||_F^2 +
#'     mu_c sum_{i<j} ||F_i^T U_i A_i - F_j^T U_j A_j||_F^2
#' }
#' 
#' It uses projected gradient descent on the Stiefel manifold to maintain
#' orthogonality of the coupling matrices A_i.
#' When \code{alpha_match > 0}, the diagonal term also nudges the joint modes
#' toward user-specified target eigenvalues as described in Eynard et al. (2015, Eq. 9).
#'
#' @examples
#' \dontrun{
#' library(multidesign)
#' # Create synthetic multi-modal data
#' X1 <- matrix(rnorm(150), 50, 3)
#' X2 <- matrix(rnorm(200), 50, 4)
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(
#'   domain1 = multidesign(X1, data.frame(id = 1:50)),
#'   domain2 = multidesign(X2, data.frame(id = 1:50))
#' ))
#' 
#' # Run coupled diagonalization
#' result <- coupled_diagonalization(hd, ncomp = 3)
#' }
#'
#' @references
#' Eynard et al. (2015). Multimodal manifold analysis by simultaneous 
#' diagonalization of Laplacians. IEEE TPAMI, 37(12), 2505-2517.
#'
#' @rdname coupled_diagonalization
#' @method coupled_diagonalization hyperdesign
#' @export
#' @importFrom Matrix sparseMatrix Diagonal crossprod colSums rowSums bdiag
#' @importFrom neighborweights graph_weights
#' @importFrom multivarious center init_transform multiblock_biprojector
#' @importFrom chk chk_number chk_true
coupled_diagonalization.hyperdesign <- function(data,
                                              correspondence = NULL,
                                              preproc = center(),
                                              ncomp = 10,
                                              ncomp_per_domain = 20,
                                              mu_coupling = 1,
                                              alpha_match = 0,
                                              lambda_target = NULL,
                                              knn = 10,
                                              sigma = NULL,
                                              max_iter = 200,
                                              step_size = 0.3,
                                              tol = 1e-6,
                                              verbose = FALSE,
                                              ...) {
  
  # Input validation
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(ncomp_per_domain)
  chk::chk_true(ncomp_per_domain >= ncomp)
  chk::chk_number(mu_coupling)
  chk::chk_true(mu_coupling >= 0)
  chk::chk_number(alpha_match)
  chk::chk_true(alpha_match >= 0)
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(step_size)
  chk::chk_true(step_size > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  # Extract number of domains
  m <- length(data)
  if (m < 2) {
    stop("Coupled diagonalization requires at least 2 domains", call. = FALSE)
  }
  if (!is.null(lambda_target)) {
    if (!is.list(lambda_target) || length(lambda_target) != m) {
      stop("lambda_target must be a list with one element per domain", call. = FALSE)
    }
    lambda_target <- lapply(lambda_target, function(x) {
      if (is.matrix(x)) {
        x <- diag(x)
      }
      as.numeric(x)
    })
  }
  
  # Apply preprocessing to each domain
  pdata <- multivarious::init_transform(data, preproc)
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  # Get sample block indices using the exported function
  sample_block_idx <- block_indices(pdata)
  
  sample_blocks_list <- lapply(seq_len(nrow(sample_block_idx)), function(i) {
    out <- as.integer(sample_block_idx[i, ])
    names(out) <- colnames(sample_block_idx)
    out
  })
  names(sample_blocks_list) <- names(pdata)
  
  # This call should now succeed.
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)
  
  feature_blocks <- feature_block_indices(pdata)
  
  # Extract data matrices
  X_list <- lapply(pdata, function(x) x$x)
  
  # Build correspondence matrices if not provided
  if (is.null(correspondence)) {
    # Assume full correspondence based on sample order
    n_samples <- sapply(X_list, nrow)
    
    if (length(unique(n_samples)) > 1) {
      warning("Domains have different numbers of samples. Using identity ",
              "correspondence for matching samples only.", call. = FALSE)
    }
    
    # Create identity correspondence matrices
    n_min <- min(n_samples)
    correspondence <- lapply(seq_len(m), function(i) {
      n_i <- n_samples[i]
      # Identity matrix truncated to minimum size
      F_i <- Matrix::sparseMatrix(i = seq_len(n_min), 
                                  j = seq_len(n_min), 
                                  x = 1, 
                                  dims = c(n_i, n_min))
      F_i
    })
  }
  
  # Validate correspondence matrices
  if (!is.list(correspondence) || length(correspondence) != m) {
    stop("correspondence must be a list of length ", m, " (one per domain)", 
         call. = FALSE)
  }
  corr_cols <- vapply(correspondence, ncol, integer(1))
  if (length(unique(corr_cols)) != 1) {
    stop("All correspondence matrices must share the same number of columns", call. = FALSE)
  }
  for (i in seq_len(m)) {
    if (nrow(correspondence[[i]]) != nrow(X_list[[i]])) {
      stop("correspondence[[", i, "]] has ", nrow(correspondence[[i]]),
           " rows but domain ", i, " has ", nrow(X_list[[i]]), " samples", call. = FALSE)
    }
  }
  
  # Build graph Laplacians using neighborweights
  if (verbose) message("Building graph Laplacians...")
  
  Laplacians <- vector("list", m)
  eigendecomps <- vector("list", m)
  
  for (i in seq_len(m)) {
    domain_name <- if (!is.null(names(data)) && length(names(data)) >= i) {
      names(data)[i]
    } else {
      paste0("domain_", i)
    }
    # Use neighborweights for consistency with package
    W_graph <- neighborweights::graph_weights(X_list[[i]], 
                                            k = knn,
                                            weight_mode = "heat",
                                            sigma = sigma)
    W <- neighborweights::adjacency(W_graph)
    # Symmetrize weights in case the builder returned a directed graph
    W <- 0.5 * (W + Matrix::t(W))
    Matrix::diag(W) <- 0
    
    # Compute normalized Laplacian
    d <- Matrix::rowSums(W)
    d[d == 0] <- 1
    Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
    L <- Matrix::forceSymmetric(Matrix::Diagonal(nrow(W)) - Dm12 %*% W %*% Dm12)
    
    # Store Laplacian
    Laplacians[[i]] <- L
    
    # Compute eigendecomposition
    # Request one additional eigenpair so we can drop the trivial λ≈0 mode
    n_eigs_total <- min(ncomp_per_domain + 1L, nrow(L))
    eig <- tryCatch({
      if (requireNamespace("RSpectra", quietly = TRUE) && n_eigs_total < nrow(L)) {
        RSpectra::eigs_sym(L, n_eigs_total, which = "SM")
      } else if (n_eigs_total < nrow(L) / 2 && nrow(L) > 50) {
        PRIMME::eigs_sym(L, NEig = n_eigs_total, which = "SA",
                        method = "PRIMME_DEFAULT_MIN_MATVECS")
      } else {
        eig_full <- eigen(as.matrix(L), symmetric = TRUE)
        idx <- order(eig_full$values)[seq_len(n_eigs_total)]
        list(
          vectors = eig_full$vectors[, idx, drop = FALSE],
          values = eig_full$values[idx]
        )
      }
    }, error = function(e) {
      if (verbose) message("Sparse eigensolver failed, using base eigen: ", e$message)
      eig_full <- eigen(as.matrix(L), symmetric = TRUE)
      idx <- order(eig_full$values)[seq_len(n_eigs_total)]
      list(
        vectors = eig_full$vectors[, idx, drop = FALSE],
        values = eig_full$values[idx]
      )
    })
    
    ord <- order(eig$values)
    vals <- eig$values[ord]
    vecs <- eig$vectors[, ord, drop = FALSE]
    non_trivial <- vals > 1e-12
    if (!any(non_trivial)) {
      stop("No non-trivial Laplacian eigenvalues found for domain ", domain_name,
           call. = FALSE)
    }
    keep <- head(which(non_trivial), ncomp_per_domain)
    if (length(keep) < ncomp_per_domain) {
      warning("Domain ", domain_name, ": only ", length(keep),
              " non-trivial eigenvectors available; reducing ncomp_per_domain.",
              call. = FALSE)
    }
    keep <- keep[keep <= ncol(vecs)]
    U_i <- vecs[, keep, drop = FALSE]
    lambda_i <- vals[keep]
    # Normalize spectra so diagonalization terms are comparable across domains
    lambda_ref_idx <- max(1L, min(length(lambda_i), ncomp))
    lambda_ref <- lambda_i[lambda_ref_idx]
    scale_factor <- max(lambda_ref, 1e-8)
    lambda_i <- lambda_i / scale_factor
    attr(lambda_i, "scale_factor") <- scale_factor
    
    eigendecomps[[i]] <- list(
      vectors = U_i,
      values = lambda_i
    )
  }
  
  # Extract eigenvectors and eigenvalues
  Ubar <- lapply(eigendecomps, function(x) x$vectors)
  Lambda <- lapply(eigendecomps, function(x) x$values)
  if (alpha_match > 0 && is.null(lambda_target)) {
    lambda_target <- lapply(seq_len(m), function(i) {
      vals <- Lambda[[i]]
      if (length(vals) == 0) {
        numeric(0)
      } else {
        vals[seq_len(min(ncomp, length(vals)))]
      }
    })
  }
  
  # Precompute F_i^T U_i for efficiency
  FiUbar <- lapply(seq_len(m), function(i) {
    crossprod(correspondence[[i]], Ubar[[i]])
  })
  coupling_sizes <- vapply(FiUbar, nrow, integer(1))
  coupling_scale <- sqrt(pmax(1, coupling_sizes))
  FiUbar <- Map(function(mat, s) {
    if (is.null(mat) || length(mat) == 0L || s <= 0) {
      return(mat)
    }
    mat / s
  }, FiUbar, coupling_scale)

  sum_matrix_list <- function(lst) {
    lst <- Filter(function(x) !is.null(x), lst)
    if (length(lst) == 0) return(NULL)
    acc <- lst[[1]]
    if (length(lst) > 1) {
      for (k in 2:length(lst)) {
        acc <- acc + lst[[k]]
      }
    }
    acc
  }
  
  # Initialize coupling matrices A_i on Stiefel manifold
  # Check that we have enough eigenvectors for the requested components
  min_eigenvecs <- min(sapply(eigendecomps, function(x) ncol(x$vectors)))
  if (min_eigenvecs < ncomp) {
    warning("Requested ", ncomp, " components but only ", min_eigenvecs, 
            " eigenvectors available. Reducing ncomp.", call. = FALSE)
    ncomp <- min_eigenvecs
  }
  
  A <- lapply(seq_len(m), function(i) {
    n_vecs <- ncol(eigendecomps[[i]]$vectors)
    M <- matrix(0, n_vecs, ncomp)
    M[seq_len(ncomp), ] <- diag(ncomp)
    M
  })
  
  # Optimization on Stiefel manifold with backtracking line search
  converged <- FALSE
  beta_backtrack <- 0.5
  armijo_c <- 1e-4
  max_backtracks <- 12L
  current_step_size <- step_size
  best_cost <- Inf
  best_A <- A
  cost_history <- numeric(0)
  component_history <- vector("list", max_iter)

  cost_val <- compute_cost_cd(A, Lambda, FiUbar, mu_coupling, lambda_target, alpha_match)
  cost_current <- as.numeric(cost_val)
  component_history[[1]] <- attr(cost_val, "components")
  cost_history[1] <- cost_current
  
  for (iter in seq_len(max_iter)) {
    cost_before_iter <- cost_current
    A_prev <- A
    B <- lapply(seq_len(m), function(j) {
      if (!is.null(FiUbar[[j]]) && length(FiUbar[[j]]) > 0L) {
        FiUbar[[j]] %*% A[[j]]
      } else {
        NULL
      }
    })
    sumB <- sum_matrix_list(B)
    any_update <- FALSE

    for (i in seq_len(m)) {
      grad_i <- compute_gradient_cd(i, A, Lambda, FiUbar, mu_coupling,
                                    B, sumB, lambda_target, alpha_match)
      if (any(!is.finite(grad_i))) {
        warning("Invalid gradient at iteration ", iter, " for domain ", i,
                ". Skipping update.", call. = FALSE)
        next
      }
      gnorm <- sqrt(sum(grad_i * grad_i))
      if (gnorm < 1e-12) {
        next
      }
      direction <- grad_i / gnorm
      step_try <- current_step_size
      accepted <- FALSE
      old_Ai <- A[[i]]
      old_Bi <- B[[i]]

      for (ls in seq_len(max_backtracks)) {
        candidate <- old_Ai - step_try * direction
        candidate <- qr.Q(qr(candidate))[, seq_len(ncomp), drop = FALSE]
        overlap <- crossprod(A_prev[[i]], candidate)
        align <- sign(diag(overlap))
        align[align == 0] <- 1
        candidate <- sweep(candidate, 2L, align, `*`)

        B_candidate <- if (!is.null(FiUbar[[i]]) && length(FiUbar[[i]]) > 0L) {
          FiUbar[[i]] %*% candidate
        } else {
          NULL
        }

        A_trial <- A
        A_trial[[i]] <- candidate
        cost_trial_val <- compute_cost_cd(A_trial, Lambda, FiUbar, mu_coupling,
                                          lambda_target, alpha_match)
        cost_trial <- as.numeric(cost_trial_val)

        if (is.finite(cost_trial) &&
            cost_trial <= cost_current - armijo_c * step_try * gnorm * gnorm) {
          A[[i]] <- candidate
          B[[i]] <- B_candidate
          sumB <- sum_matrix_list(B)
          cost_current <- cost_trial
          component_history[[iter]] <- attr(cost_trial_val, "components")
          accepted <- TRUE
          any_update <- TRUE
          break
        }
        step_try <- step_try * beta_backtrack
      }

      if (!accepted) {
        # restore on rejection
        A[[i]] <- old_Ai
        B[[i]] <- old_Bi
      }
    }

    if (!any_update) {
      if (verbose) message("No successful updates at iteration ", iter)
      if (iter > 1) {
        converged <- TRUE
        break
      }
    }

    if (cost_current < best_cost) {
      best_cost <- cost_current
      best_A <- A
    }

    cost_history[iter + 1] <- cost_current
    if (is.null(component_history[[iter]])) {
      component_history[[iter]] <- attr(
        compute_cost_cd(A, Lambda, FiUbar, mu_coupling, lambda_target, alpha_match),
        "components"
      )
    }

    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %3d   cost = %.6e", iter, cost_current))
    }

    rel_change <- abs(cost_before_iter - cost_current) / (abs(cost_before_iter) + 1e-9)
    if (rel_change < tol) {
      converged <- TRUE
      if (verbose) message("Converged at iteration ", iter)
      break
    }

    if (cost_current > cost_before_iter) {
      current_step_size <- max(current_step_size * beta_backtrack, 1e-6)
    } else {
      current_step_size <- min(current_step_size * 1.05, step_size * 5)
    }
  }

  if (!converged && verbose) {
    warning("Did not converge within ", max_iter, " iterations", call. = FALSE)
  }

  if (best_cost < cost_current) {
    A <- best_A
    cost_current <- best_cost
  }

  final_cost <- cost_current
  iterations_run <- if (exists("iter")) iter else max_iter
  cost_history_valid <- cost_history[!is.na(cost_history)]
  component_history_valid <- component_history[seq_len(min(iterations_run, length(component_history)))]
  component_history_valid <- Filter(function(x) !is.null(x), component_history_valid)
  if (length(cost_history_valid) > 0) {
    names(cost_history_valid) <- paste0("iter_", seq_along(cost_history_valid) - 1L)
  }
  if (length(component_history_valid) > 0) {
    names(component_history_valid) <- paste0("iter_", seq_along(component_history_valid))
  }

  # Compute final coupled bases V_i = U_i A_i
  V <- lapply(seq_len(m), function(i) Ubar[[i]] %*% A[[i]])
  names(V) <- names(data)
  
  # Concatenate scores for multiblock structure
  s <- do.call(rbind, V)
  
  # Prepare v matrix for projections using a ridge-regularized blockwise fit
  ridge_penalty <- 1e-3
  v_blocks <- vector("list", m)
  for (i in seq_len(m)) {
    Xi <- X_list[[i]]
    Vi <- V[[i]]
    if (ncol(Xi) == 0) {
      v_blocks[[i]] <- matrix(0, 0, ncol(Vi))
      next
    }
    gram <- as.matrix(crossprod(Xi))
    rhs <- as.matrix(crossprod(Xi, Vi))
    penalty <- ridge_penalty * diag(ncol(Xi))
    v_blocks[[i]] <- solve(gram + penalty, rhs)
  }
  v <- do.call(rbind, v_blocks)
  
  diagnostics <- list(
    cost_history = cost_history_valid,
    component_history = component_history_valid,
    coupling_scale = coupling_scale,
    stiefel_factors = A,
    subspace_bases = Ubar
  )

  new_alignment_result(
    scores = s,
    loadings = v,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "coupled_diagonalization",
    extras = list(
      coupled_bases = V,
      eigenvalues = Lambda,
      converged = converged,
      iterations = iterations_run,
      final_cost = final_cost,
      diagnostics = diagnostics
    )
  )
}

# Helper function: compute gradient for domain i
compute_gradient_cd <- function(idx,
                                A,
                                Lambda,
                                FiUbar,
                                mu_coupling,
                                FiUbarA,
                                sumB,
                                lambda_target = NULL,
                                alpha_match = 0) {
  m <- length(A)
  Ai <- A[[idx]]
  lambda_i <- pmax(Lambda[[idx]], 1e-8)
  LiAi <- sweep(Ai, 1L, lambda_i, `*`)
  Mi <- crossprod(Ai, LiAi)
  S <- Mi - diag(diag(Mi))
  
  # Gradient of off-diagonalization term
  grad_diag <- 4 * (LiAi %*% S)
  
  grad_match <- 0
  if (alpha_match > 0 && !is.null(lambda_target) && !is.null(lambda_target[[idx]])) {
    target_vec <- lambda_target[[idx]]
    if (is.matrix(target_vec)) {
      target_vec <- diag(target_vec)
    }
    if (length(target_vec) < ncol(Ai)) {
      target_vec <- c(target_vec, rep(tail(target_vec, 1L), ncol(Ai) - length(target_vec)))
    }
    target_vec <- target_vec[seq_len(ncol(Ai))]
    match_residual <- Mi
    diag(match_residual) <- diag(match_residual) - target_vec
    grad_match <- 4 * (LiAi %*% match_residual)
  }
  
  Bi <- FiUbarA[[idx]]
  grad_coupling <- 0
  if (mu_coupling > 0 && m > 1 && !is.null(Bi) && !is.null(sumB) && !is.null(FiUbar[[idx]])) {
    grad_coupling <- 2 * mu_coupling * (
      m * crossprod(FiUbar[[idx]], Bi) - crossprod(FiUbar[[idx]], sumB)
    )
  }
  
  grad_diag + grad_coupling + alpha_match * grad_match
}

# Helper function: compute objective cost
compute_cost_cd <- function(A,
                            Lambda,
                            FiUbar,
                            mu_coupling,
                            lambda_target = NULL,
                            alpha_match = 0) {
  m <- length(A)
  cost_d <- 0
  cost_match <- 0
  B <- vector("list", m)
  norms_B <- numeric(m)
  
  for (i in seq_len(m)) {
    Ai <- A[[i]]
    lambda_i <- pmax(Lambda[[i]], 1e-8)
    LiAi <- sweep(Ai, 1L, lambda_i, `*`)
    Mi <- crossprod(Ai, LiAi)
    off_diag <- Mi - diag(diag(Mi))
    cost_d <- cost_d + sum(off_diag * off_diag)
    if (alpha_match > 0 && !is.null(lambda_target) && !is.null(lambda_target[[i]])) {
      target_vec <- lambda_target[[i]]
      if (is.matrix(target_vec)) {
        target_vec <- diag(target_vec)
      }
      if (length(target_vec) < ncol(Ai)) {
        target_vec <- c(target_vec, rep(tail(target_vec, 1L), ncol(Ai) - length(target_vec)))
      }
      target_vec <- target_vec[seq_len(ncol(Ai))]
      match_residual <- Mi
      diag(match_residual) <- diag(match_residual) - target_vec
      cost_match <- cost_match + sum(match_residual * match_residual)
    }
    if (!is.null(FiUbar[[i]]) && length(FiUbar[[i]]) > 0L) {
      B[[i]] <- FiUbar[[i]] %*% Ai
    } else {
      B[[i]] <- NULL
    }
    norms_B[i] <- if (!is.null(B[[i]])) sum(B[[i]] * B[[i]]) else 0
  }
  
  cost_c <- 0
  if (mu_coupling > 0 && m > 1) {
    non_null_B <- Filter(function(x) !is.null(x), B)
    if (length(non_null_B) > 0) {
      sumB <- non_null_B[[1]]
      if (length(non_null_B) > 1) {
        for (idx in 2:length(non_null_B)) {
          sumB <- sumB + non_null_B[[idx]]
        }
      }
      cost_c <- m * sum(norms_B) - sum(sumB * sumB)
    }
  }

  total <- cost_d + mu_coupling * cost_c + alpha_match * cost_match
  attr(total, "components") <- list(
    diagonal = cost_d,
    coupling = cost_c,
    match = cost_match
  )
  total
}

#' @rdname coupled_diagonalization
#' @method coupled_diagonalization default
#' @export
coupled_diagonalization.default <- function(data, ...) {
  stop("No applicable method for coupled_diagonalization. ",
       "data must be a hyperdesign object.", call. = FALSE)
}

# For backward compatibility - keep the original function as internal
#' @keywords internal
coupled_diagonalisation_internal <- function(
    X_list,                # list of data matrices (modality i: n_i × d_i)
    F_list,                # list of correspondence matrices (n_i × q)
    k        = 10,         # # joint eigenvectors wanted
    kprime   = 20,         # # low-frequency eigenvectors per modality (≥k)
    mu_c     = 1,          # coupling weight
    knn      = 10,         # k-NN for graph
    sigma    = NULL,       # Gaussian σ (NULL => median NN dist)
    maxIter  = 200,        # outer iterations
    step0    = 0.3,        # base step size
    tol      = 1e-6) {     # relative cost tolerance
  
  # This is the original implementation preserved for reference
  # Users should use the new coupled_diagonalization() function instead
  
  m <- length(X_list)
  stopifnot(length(F_list) == m, kprime >= k)
  
  # Build Laplacians and get first k′ eigenpairs
  Ubar    <- vector("list", m)
  Lambda  <- vector("list", m)
  FiUbar  <- vector("list", m)
  
  for (i in seq_len(m)) {
    # Use FNN for backward compatibility in internal function
    nn <- FNN::get.knn(X_list[[i]], k = knn)
    n  <- nrow(X_list[[i]])
    i_idx  <- rep(seq_len(n), each = knn)
    j_idx  <- as.vector(nn$nn.index)
    d2 <- as.vector(nn$nn.dist)^2
    
    sigma_i <- sigma
    if (is.null(sigma_i)) {
      sigma_i <- median(sqrt(d2))
    }
    
    w  <- exp(-d2 / (2 * sigma_i^2))
    W  <- Matrix::sparseMatrix(i = c(i_idx, j_idx), j = c(j_idx, i_idx), 
                              x = c(w, w), dims = c(n, n))
    diag(W) <- 0
    
    # Laplacian
    d <- rowSums(W)
    d[d == 0] <- 1
    Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
    L <- Matrix::Diagonal(n = nrow(W)) - Dm12 %*% W %*% Dm12
    
    eig <- RSpectra::eigs_sym(L, kprime, which = "SM")
    Ubar[[i]]   <- eig$vectors
    Lambda[[i]] <- eig$values
    FiUbar[[i]] <- crossprod(F_list[[i]], Ubar[[i]])
  }
  
  # Initialize A_i
  A <- lapply(seq_len(m), function(i) {
    M <- matrix(0, kprime, k)
    M[seq_len(k), ] <- diag(k)
    M
  })
  
  # Cost function
  cost_function <- function(A) {
    cost_d <- 0
    cost_c <- 0
    for (i in seq_len(m)) {
      Ai <- A[[i]]
      Li_diag <- diag(Lambda[[i]])
      Mi <- crossprod(Ai, Li_diag %*% Ai)
      cost_d <- cost_d + sum((Mi - diag(diag(Mi)))^2)
    }
    
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        diff <- FiUbar[[i]] %*% A[[i]] - FiUbar[[j]] %*% A[[j]]
        cost_c <- cost_c + sum(diff * diff)
      }
    }
    cost_d + mu_c * cost_c
  }
  
  # Gradient
  gradient_block <- function(idx, A) {
    Ai      <- A[[idx]]
    Li_diag <- diag(Lambda[[idx]])
    Mi      <- crossprod(Ai, Li_diag %*% Ai)
    S       <- Mi - diag(diag(Mi))
    
    grad_diag <- 4 * (Li_diag %*% Ai %*% S)
    
    accum <- matrix(0, nrow = kprime, ncol = k)
    for (j in seq_len(m)) {
      if (j != idx) {
        diff   <- FiUbar[[idx]] %*% Ai - FiUbar[[j]] %*% A[[j]]
        accum  <- accum + crossprod(FiUbar[[idx]], diff)
      }
    }
    grad_coupling <- 2 * mu_c * accum
    
    grad_diag + grad_coupling
  }
  
  # Optimization
  cost_old <- Inf
  for (iter in seq_len(maxIter)) {
    for (i in seq_len(m)) {
      grad_i <- gradient_block(i, A)
      gnorm  <- sqrt(sum(grad_i * grad_i))
      eta    <- step0 / max(gnorm, 1e-9)
      A_new  <- A[[i]] - eta * grad_i
      A[[i]] <- qr.Q(qr(A_new))[, seq_len(k)]
    }
    
    cost_new <- cost_function(A)
    
    if ((cost_old - cost_new) / (abs(cost_old) + 1e-9) < tol) {
      message("Convergence criterion met.")
      break
    }
    cost_old <- cost_new
  }
  
  # Final coupled bases
  V <- lapply(seq_len(m), function(i) Ubar[[i]] %*% A[[i]])
  names(V) <- paste0("modality_", seq_len(m))
  invisible(V)
}
