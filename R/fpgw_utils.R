#' FPGW utility functions
#' @keywords internal

#' Project onto doubly stochastic matrices
#' @keywords internal
project_onto_simplex_doubly_stochastic <- function(P, p, q, max_iter = 100, tol = 1e-9) {
  # Use Sinkhorn-Knopp algorithm to project onto doubly stochastic matrices
  # with specified marginals p and q
  
  n <- length(p)
  m <- length(q)
  
  # Ensure P is non-negative
  P[P < 0] <- 0
  
  # Add small epsilon to avoid division by zero
  P <- P + 1e-20
  
  # Sinkhorn-Knopp iterations
  for (iter in seq_len(max_iter)) {
    # Row normalization
    row_sums <- rowSums(P)
    P <- sweep(P, 1, p / row_sums, "*")
    
    # Column normalization
    col_sums <- colSums(P)
    P <- sweep(P, 2, q / col_sums, "*")
    
    # Check convergence
    if (max(abs(rowSums(P) - p)) < tol && max(abs(colSums(P) - q)) < tol) {
      break
    }
  }
  
  return(P)
}

#' Prepare OT data from hyperdesign
#' @keywords internal
prepare_ot_data <- function(hyperdesign) {
  resolved <- resolve_hyperdesign(hyperdesign)

  if (resolved$n_domains < 2) {
    stop("FPGW requires at least 2 domains", call. = FALSE)
  }

  X_list <- lapply(resolved$domains, function(domain) domain$x)
  n_samples <- vapply(X_list, nrow, integer(1))

  list(
    X_list = X_list,
    n_samples = n_samples,
    n_domains = resolved$n_domains,
    domain_names = resolved$domain_names
  )
}

#' Compute distance matrix
#' @keywords internal
compute_distance_matrix <- function(X, metric = "euclidean") {
  if (metric == "euclidean") {
    as.matrix(dist(X, method = "euclidean"))
  } else if (metric == "manhattan") {
    as.matrix(dist(X, method = "manhattan"))
  } else if (metric == "cosine") {
    # Cosine distance = 1 - cosine similarity
    X_norm <- X / sqrt(rowSums(X^2))
    1 - tcrossprod(X_norm)
  } else {
    stop("Unsupported metric: ", metric, call. = FALSE)
  }
}

#' Compute feature cost between domains
#' @keywords internal
compute_feature_cost <- function(X1, X2, metric = "euclidean") {
  n1 <- nrow(X1)
  n2 <- nrow(X2)

  if (ncol(X1) != ncol(X2)) {
    d_max <- max(ncol(X1), ncol(X2))
    X1_pad <- matrix(0, n1, d_max)
    X2_pad <- matrix(0, n2, d_max)
    X1_pad[, seq_len(ncol(X1))] <- X1
    X2_pad[, seq_len(ncol(X2))] <- X2
    X1 <- X1_pad
    X2 <- X2_pad
  }
  
  if (metric == "euclidean") {
    # Efficient squared Euclidean distance computation
    # ||x - y||^2 = ||x||^2 + ||y||^2 - 2<x,y>
    X1_sqnorms <- rowSums(X1^2)
    X2_sqnorms <- rowSums(X2^2)
    XY <- tcrossprod(X1, X2)
    
    C <- outer(X1_sqnorms, rep(1, n2)) + 
         outer(rep(1, n1), X2_sqnorms) - 
         2 * XY
    
    # Ensure non-negative (numerical issues)
    C[C < 0] <- 0
    sqrt(C)
  } else if (metric == "manhattan") {
    # Manhattan distance
    C <- matrix(0, n1, n2)
    for (i in seq_len(n1)) {
      diffs <- abs(matrix(X1[i, ], n2, ncol(X1), byrow = TRUE) - X2)
      C[i, ] <- rowSums(diffs)
    }
    C
  } else if (metric == "cosine") {
    # Cosine distance
    X1_norm <- X1 / sqrt(rowSums(X1^2) + 1e-10)
    X2_norm <- X2 / sqrt(rowSums(X2^2) + 1e-10)
    1 - tcrossprod(X1_norm, X2_norm)
  } else {
    stop("Unsupported metric: ", metric, call. = FALSE)
  }
}
