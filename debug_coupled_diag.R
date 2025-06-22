# Debug script for coupled_diagonalization numerical issues
devtools::load_all(".")
library(multidesign)
library(multivarious)
library(Matrix)
library(tibble)  # Required for multidesign

# Create minimal test case
set.seed(123)
n <- 30
X1 <- matrix(rnorm(n * 3), n, 3)
X2 <- matrix(rnorm(n * 4), n, 4)

design1 <- data.frame(sample_id = 1:n)
design2 <- data.frame(sample_id = 1:n)

md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

# Add debugging to the functions
debug_compute_cost <- function(A, Lambda, FiUbar, mu_coupling) {
  m <- length(A)
  cost_d <- 0
  cost_c <- 0
  
  cat("\n=== Computing cost ===\n")
  
  # Diagonalization cost
  for (i in seq_len(m)) {
    cat("Domain", i, ":\n")
    cat("  Lambda range:", range(Lambda[[i]]), "\n")
    cat("  A dimensions:", dim(A[[i]]), "\n")
    
    # Ensure eigenvalues are valid
    Li_diag <- diag(pmax(Lambda[[i]], 1e-8))
    cat("  Li_diag range:", range(diag(Li_diag)), "\n")
    
    Mi <- crossprod(A[[i]], Li_diag %*% A[[i]])
    cat("  Mi range:", range(Mi), "\n")
    
    diag_cost_i <- sum((Mi - diag(diag(Mi)))^2)
    cat("  Diag cost:", diag_cost_i, "\n")
    
    if (is.na(diag_cost_i) || !is.finite(diag_cost_i)) {
      cat("  WARNING: Invalid diag cost!\n")
      cat("  Mi:\n")
      print(Mi)
    }
    
    cost_d <- cost_d + diag_cost_i
  }
  
  # Coupling cost
  for (i in 1:(m - 1)) {
    for (j in (i + 1):m) {
      diff <- FiUbar[[i]] %*% A[[i]] - FiUbar[[j]] %*% A[[j]]
      coupling_cost_ij <- sum(diff * diff)
      cat("  Coupling", i, "-", j, ":", coupling_cost_ij, "\n")
      
      if (is.na(coupling_cost_ij) || !is.finite(coupling_cost_ij)) {
        cat("  WARNING: Invalid coupling cost!\n")
        cat("  diff range:", range(diff), "\n")
      }
      
      cost_c <- cost_c + coupling_cost_ij
    }
  }
  
  total_cost <- cost_d + mu_coupling * cost_c
  cat("Total cost:", total_cost, "(diag:", cost_d, ", coupling:", mu_coupling * cost_c, ")\n")
  
  total_cost
}

# Run with debugging
cat("Running coupled_diagonalization with debugging...\n")

# First, let's check the eigendecomposition step
result <- tryCatch({
  # Extract data and preprocess
  pdata <- multivarious::init_transform(hd, multivarious::center())
  X_list <- lapply(pdata, function(x) x$x)
  
  cat("\nData dimensions:\n")
  for (i in seq_along(X_list)) {
    cat("  Domain", i, ":", dim(X_list[[i]]), "\n")
  }
  
  # Build graphs and compute eigendecompositions
  cat("\nBuilding graphs and computing eigendecompositions...\n")
  
  Laplacians <- vector("list", 2)
  eigendecomps <- vector("list", 2)
  
  for (i in 1:2) {
    cat("\nDomain", i, ":\n")
    
    # Build graph
    W_graph <- neighborweights::graph_weights(X_list[[i]], 
                                            k = 10,
                                            weight_mode = "heat",
                                            sigma = NULL)
    W <- neighborweights::adjacency(W_graph)
    cat("  W dimensions:", dim(W), ", nnz:", Matrix::nnzero(W), "\n")
    
    # Compute normalized Laplacian
    d <- Matrix::rowSums(W)
    cat("  Degree range:", range(d), "\n")
    
    d[d == 0] <- 1
    Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
    L <- Matrix::Diagonal(nrow(W)) - Dm12 %*% W %*% Dm12
    
    cat("  Laplacian dimensions:", dim(L), "\n")
    cat("  Laplacian type:", class(L), "\n")
    
    # Check if Laplacian is valid
    L_dense <- as.matrix(L)
    cat("  Laplacian range:", range(L_dense), "\n")
    cat("  Is symmetric:", isSymmetric(L_dense), "\n")
    
    # Compute eigendecomposition
    n_eigs <- min(20, nrow(L) - 1)
    
    eig <- eigen(L_dense, symmetric = TRUE)
    idx <- order(abs(eig$values))[1:n_eigs]
    
    eigendecomps[[i]] <- list(
      vectors = eig$vectors[, idx, drop = FALSE],
      values = eig$values[idx]
    )
    
    cat("  Eigenvalues:", head(eigendecomps[[i]]$values, 10), "\n")
    cat("  Min eigenvalue:", min(eigendecomps[[i]]$values), "\n")
    cat("  Eigenvector dimensions:", dim(eigendecomps[[i]]$vectors), "\n")
    
    Laplacians[[i]] <- L
  }
  
  # Extract eigenvectors and eigenvalues
  Ubar <- lapply(eigendecomps, function(x) x$vectors)
  Lambda <- lapply(eigendecomps, function(x) x$values)
  
  # Build correspondence matrices
  F1 <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1, dims = c(n, n))
  F2 <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1, dims = c(n, n))
  correspondence <- list(F1, F2)
  
  # Precompute F_i^T U_i
  FiUbar <- lapply(seq_len(2), function(i) {
    crossprod(correspondence[[i]], Ubar[[i]])
  })
  
  cat("\nFiUbar dimensions:\n")
  for (i in 1:2) {
    cat("  Domain", i, ":", dim(FiUbar[[i]]), "\n")
  }
  
  # Initialize coupling matrices
  ncomp <- 3
  A <- lapply(seq_len(2), function(i) {
    n_vecs <- ncol(eigendecomps[[i]]$vectors)
    M <- matrix(0, n_vecs, ncomp)
    M[seq_len(ncomp), ] <- diag(ncomp)
    M
  })
  
  cat("\nInitial A dimensions:\n")
  for (i in 1:2) {
    cat("  Domain", i, ":", dim(A[[i]]), "\n")
  }
  
  # Test cost computation
  cat("\nTesting initial cost computation:\n")
  initial_cost <- debug_compute_cost(A, Lambda, FiUbar, mu_coupling = 1)
  
  if (is.na(initial_cost) || !is.finite(initial_cost)) {
    stop("Initial cost is invalid!")
  }
  
  # Test gradient computation
  cat("\nTesting gradient computation:\n")
  
  # Simplified gradient function with debugging
  debug_gradient <- function(idx, A, Lambda, FiUbar, mu_coupling) {
    cat("\nComputing gradient for domain", idx, "\n")
    
    Ai <- A[[idx]]
    Li_diag <- diag(pmax(Lambda[[idx]], 1e-8))
    Mi <- crossprod(Ai, Li_diag %*% Ai)
    S <- Mi - diag(diag(Mi))
    
    cat("  Ai range:", range(Ai), "\n")
    cat("  Li_diag range:", range(diag(Li_diag)), "\n")
    cat("  Mi range:", range(Mi), "\n")
    cat("  S range:", range(S), "\n")
    
    # Gradient of off-diagonalization term
    grad_diag <- 4 * (Li_diag %*% Ai %*% S)
    cat("  grad_diag range:", range(grad_diag), "\n")
    
    # Gradient of coupling term
    accum <- matrix(0, nrow = nrow(Ai), ncol = ncol(Ai))
    for (j in seq_len(2)) {
      if (j != idx) {
        diff <- FiUbar[[idx]] %*% Ai - FiUbar[[j]] %*% A[[j]]
        accum <- accum + crossprod(FiUbar[[idx]], diff)
      }
    }
    grad_coupling <- 2 * mu_coupling * accum
    cat("  grad_coupling range:", range(grad_coupling), "\n")
    
    total_grad <- grad_diag + grad_coupling
    cat("  total gradient range:", range(total_grad), "\n")
    
    total_grad
  }
  
  # Test gradient for each domain
  for (i in 1:2) {
    grad_i <- debug_gradient(i, A, Lambda, FiUbar, mu_coupling = 1)
    if (any(is.na(grad_i)) || any(!is.finite(grad_i))) {
      cat("WARNING: Invalid gradient for domain", i, "\n")
    }
  }
  
  # Test one optimization step
  cat("\n\nTesting one optimization step:\n")
  step_size <- 0.3
  
  for (i in 1:2) {
    grad_i <- debug_gradient(i, A, Lambda, FiUbar, mu_coupling = 1)
    gnorm <- sqrt(sum(grad_i * grad_i))
    eta <- step_size / max(gnorm, 1e-9)
    
    cat("\nDomain", i, "step:\n")
    cat("  Gradient norm:", gnorm, "\n")
    cat("  Step size eta:", eta, "\n")
    
    # Gradient step
    A_new <- A[[i]] - eta * grad_i
    cat("  A_new range:", range(A_new), "\n")
    
    # Project back to Stiefel manifold
    qr_result <- qr(A_new)
    A[[i]] <- qr.Q(qr_result)[, seq_len(ncomp), drop = FALSE]
    cat("  A after projection range:", range(A[[i]]), "\n")
    
    # Check orthogonality
    ATA <- crossprod(A[[i]])
    cat("  A^T A diagonal:", diag(ATA), "\n")
    cat("  A^T A off-diagonal max:", max(abs(ATA - diag(diag(ATA)))), "\n")
  }
  
  # Compute new cost
  cat("\nComputing cost after one step:\n")
  new_cost <- debug_compute_cost(A, Lambda, FiUbar, mu_coupling = 1)
  
  cat("\nCost change:", initial_cost, "->", new_cost, "\n")
  
  "Debug completed"
  
}, error = function(e) {
  cat("\nERROR:", e$message, "\n")
  traceback()
  NULL
})