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
#'   \item \code{converged}: Logical indicating convergence
#'   \item \code{iterations}: Number of iterations performed
#'   \item \code{final_cost}: Final objective value
#' }
#'
#' @details
#' The algorithm minimizes:
#' \deqn{F = \sum_i ||A_i^T \Lambda_i A_i - diag(A_i^T \Lambda_i A_i)||_F^2 + 
#'       \mu_c \sum_{i<j} ||F_i^T U_i A_i - F_j^T U_j A_j||_F^2}
#'
#' Using projected gradient descent on the Stiefel manifold to maintain 
#' orthogonality of the coupling matrices A_i.
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
  
  # Apply preprocessing to each domain
  pdata <- multivarious::init_transform(data, preproc)
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  # Get sample block indices using the exported function
  sample_block_idx <- block_indices(pdata)
  
  # CRITICAL: Use the 'split' method to create the list of block indices.
  # This preserves the 'start' and 'end' names within each list element,
  # which concat_pre_processors appears to require.
  sample_blocks_list <- split(sample_block_idx, row(sample_block_idx))
  
  # The names of the list elements should match the data block names
  names(sample_blocks_list) <- names(pdata)
  
  # This call should now succeed.
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)
  
  # Compute feature block indices for multiblock_biprojector
  feat_per_block <- vapply(pdata, function(b) ncol(b$x), integer(1))
  end_idx   <- cumsum(feat_per_block)
  start_idx <- c(1L, head(end_idx, -1) + 1L)
  
  feature_block_idx <- lapply(seq_along(feat_per_block), function(i) {
    start_idx[i]:end_idx[i]
  })
  names(feature_block_idx) <- names(pdata)
  
  # Use feature block indices for multiblock_biprojector
  block_indices <- feature_block_idx
  
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
  
  # Build graph Laplacians using neighborweights
  if (verbose) message("Building graph Laplacians...")
  
  Laplacians <- vector("list", m)
  eigendecomps <- vector("list", m)
  
  for (i in seq_len(m)) {
    # Use neighborweights for consistency with package
    W_graph <- neighborweights::graph_weights(X_list[[i]], 
                                            k = knn,
                                            weight_mode = "heat",
                                            sigma = sigma)
    W <- neighborweights::adjacency(W_graph)
    
    # Compute normalized Laplacian
    d <- Matrix::rowSums(W)
    d[d == 0] <- 1
    Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
    L <- Matrix::Diagonal(nrow(W)) - Dm12 %*% W %*% Dm12
    
    # Ensure Laplacian is exactly symmetric (numerical precision issues)
    L <- Matrix::forceSymmetric(L)
    
    # Store Laplacian
    Laplacians[[i]] <- L
    
    # Compute eigendecomposition
    # Don't ask for more eigenvectors than the matrix dimension
    n_eigs <- min(ncomp_per_domain, nrow(L) - 1)
    
    # Use PRIMME for eigendecomposition - it handles sparse matrices well
    eig <- tryCatch({
      if (n_eigs < nrow(L) / 2 && nrow(L) > 50) {
        # Use PRIMME for larger sparse matrices
        PRIMME::eigs_sym(L, NEig = n_eigs, which = "SA", 
                        method = 'PRIMME_DEFAULT_MIN_MATVECS')
      } else {
        # Use base eigen for small matrices
        eig_full <- eigen(as.matrix(L), symmetric = TRUE)
        idx <- order(abs(eig_full$values))[1:n_eigs]  # Smallest magnitude
        list(
          vectors = eig_full$vectors[, idx, drop = FALSE],
          values = eig_full$values[idx]
        )
      }
    }, error = function(e) {
      # Fallback to base eigen if PRIMME fails
      if (verbose) message("PRIMME failed, using base eigen: ", e$message)
      eig_full <- eigen(as.matrix(L), symmetric = TRUE)
      idx <- order(abs(eig_full$values))[1:n_eigs]
      list(
        vectors = eig_full$vectors[, idx, drop = FALSE],
        values = eig_full$values[idx]
      )
    })
    
    # Clean up very small eigenvalues to avoid numerical issues
    # Eigenvalues of normalized Laplacian should be in [0, 2]
    eig$values[abs(eig$values) < 1e-10] <- 1e-10
    
    eigendecomps[[i]] <- list(
      vectors = eig$vectors,
      values = eig$values
    )
  }
  
  # Extract eigenvectors and eigenvalues
  Ubar <- lapply(eigendecomps, function(x) x$vectors)
  Lambda <- lapply(eigendecomps, function(x) x$values)
  
  # Precompute F_i^T U_i for efficiency
  FiUbar <- lapply(seq_len(m), function(i) {
    crossprod(correspondence[[i]], Ubar[[i]])
  })
  
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
  
  # Optimization on Stiefel manifold with momentum
  cost_old <- Inf
  cost_new <- Inf  # Initialize in case loop breaks early
  converged <- FALSE
  
  # Initialize momentum terms
  momentum <- lapply(A, function(a) matrix(0, nrow(a), ncol(a)))
  momentum_rate <- 0.9
  
  # Track best solution
  best_cost <- Inf
  best_A <- A
  
  # Adaptive step size
  current_step_size <- step_size
  
  for (iter in seq_len(max_iter)) {
    # Store previous A for potential rollback
    A_prev <- A
    
    # Update each domain's coupling matrix
    for (i in seq_len(m)) {
      # Compute gradient
      grad_i <- compute_gradient_cd(i, A, Lambda, FiUbar, mu_coupling)
      
      # Check for valid gradient
      if (any(is.na(grad_i)) || any(!is.finite(grad_i))) {
        warning("Invalid gradient at iteration ", iter, " for domain ", i, 
                ". Using smaller step size.", call. = FALSE)
        grad_i[!is.finite(grad_i)] <- 0
      }
      
      # Apply momentum
      momentum[[i]] <- momentum_rate * momentum[[i]] + (1 - momentum_rate) * grad_i
      
      # Adaptive step size with bounds
      gnorm <- sqrt(sum(momentum[[i]] * momentum[[i]]))
      if (gnorm < 1e-12) {
        # Gradient is essentially zero, skip update
        next
      }
      eta <- current_step_size / gnorm
      eta <- min(eta, 1.0)  # Cap maximum step size
      
      # Gradient step with momentum
      A_new <- A[[i]] - eta * momentum[[i]]
      
      # Project back to Stiefel manifold via QR
      A[[i]] <- qr.Q(qr(A_new))[, seq_len(ncomp), drop = FALSE]
    }
    
    # Compute cost
    cost_new <- compute_cost_cd(A, Lambda, FiUbar, mu_coupling)
    
    # If cost increased significantly, reduce step size and rollback
    if (cost_new > cost_old * 1.1 && iter > 1) {
      current_step_size <- current_step_size * 0.5
      # Don't let step size get too small
      if (current_step_size < 1e-6) {
        if (verbose) message("Step size too small, stopping")
        break
      }
      A <- A_prev
      momentum <- lapply(momentum, function(m) m * 0.5)  # Reduce momentum too
      cost_new <- cost_old
      if (verbose && iter %% 10 == 0) message("  Reducing step size to ", current_step_size)
      next
    }
    
    # Increase step size if cost is decreasing steadily
    if (iter > 5 && cost_new < cost_old * 0.99) {
      current_step_size <- min(current_step_size * 1.1, step_size)
    }
    
    # Track best solution
    if (cost_new < best_cost) {
      best_cost <- cost_new
      best_A <- A
    }
    
    # Check for valid cost
    if (is.na(cost_new) || !is.finite(cost_new)) {
      warning("Invalid cost computed at iteration ", iter, ". Stopping.", call. = FALSE)
      break
    }
    
    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %3d   cost = %.6e", iter, cost_new))
    }
    
    # Check convergence (skip on first iteration)
    if (iter > 1) {
      rel_change <- abs(cost_old - cost_new) / (abs(cost_old) + 1e-9)
      if (rel_change < tol) {
        converged <- TRUE
        if (verbose) message("Converged at iteration ", iter)
        break
      }
      
      # Also check if we've been stuck at the same cost for multiple iterations
      if (iter > 10 && abs(cost_new - cost_old) < 1e-12) {
        converged <- TRUE
        if (verbose) message("Cost stabilized at iteration ", iter)
        break
      }
    }
    cost_old <- cost_new
  }
  
  if (!converged && verbose) {
    warning("Did not converge within ", max_iter, " iterations", call. = FALSE)
  }
  
  # Use best solution found
  if (best_cost < cost_new) {
    A <- best_A
    cost_new <- best_cost
    if (verbose) message("Using best solution with cost ", best_cost)
  }
  
  # Compute final coupled bases V_i = U_i A_i
  V <- lapply(seq_len(m), function(i) Ubar[[i]] %*% A[[i]])
  names(V) <- names(data)
  
  # Concatenate scores for multiblock structure
  s <- do.call(rbind, V)
  
  # Compute standard deviations
  sdev <- apply(s, 2, sd)
  
  # Prepare v matrix for projections
  # For coupled diagonalization, we need to compute projection coefficients
  # This is done by projecting the data onto the coupled bases
  X_concat <- Matrix::bdiag(X_list)
  
  # Simple projection approach - can be refined later
  v <- as.matrix(crossprod(X_concat, s))
  
  # Return multiblock_biprojector
  multivarious::multiblock_biprojector(
    v = v,
    s = s,
    sdev = sdev,
    preproc = proc,
    block_indices = block_indices,
    coupled_bases = V,
    eigenvalues = Lambda,
    converged = converged,
    iterations = if(converged) iter else max_iter,
    final_cost = cost_new,
    classes = "coupled_diagonalization"
  )
}

# Helper function: compute gradient for domain i
compute_gradient_cd <- function(idx, A, Lambda, FiUbar, mu_coupling) {
  m <- length(A)
  Ai <- A[[idx]]
  # Ensure eigenvalues are valid (small positive values)
  Li_diag <- diag(pmax(Lambda[[idx]], 1e-8))
  Mi <- crossprod(Ai, Li_diag %*% Ai)
  S <- Mi - diag(diag(Mi))
  
  # Gradient of off-diagonalization term
  grad_diag <- 4 * (Li_diag %*% Ai %*% S)
  
  # Gradient of coupling term
  accum <- matrix(0, nrow = nrow(Ai), ncol = ncol(Ai))
  for (j in seq_len(m)) {
    if (j != idx) {
      diff <- FiUbar[[idx]] %*% Ai - FiUbar[[j]] %*% A[[j]]
      accum <- accum + crossprod(FiUbar[[idx]], diff)
    }
  }
  grad_coupling <- 2 * mu_coupling * accum
  
  grad_diag + grad_coupling
}

# Helper function: compute objective cost
compute_cost_cd <- function(A, Lambda, FiUbar, mu_coupling) {
  m <- length(A)
  cost_d <- 0
  cost_c <- 0
  
  # Diagonalization cost
  for (i in seq_len(m)) {
    # Ensure eigenvalues are valid (small positive values)
    Li_diag <- diag(pmax(Lambda[[i]], 1e-8))
    Mi <- crossprod(A[[i]], Li_diag %*% A[[i]])
    cost_d <- cost_d + sum((Mi - diag(diag(Mi)))^2)
  }
  
  # Coupling cost
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      diff <- FiUbar[[i]] %*% A[[i]] - FiUbar[[j]] %*% A[[j]]
      cost_c <- cost_c + sum(diff * diff)
    }
  }
  
  cost_d + mu_coupling * cost_c
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
    
    if (is.null(sigma))
      sigma <- median(sqrt(d2))
    
    w  <- exp(-d2 / (2 * sigma^2))
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