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
  if (!inherits(hyperdesign, "hyperdesign")) {
    stop("Input must be a hyperdesign object", call. = FALSE)
  }
  
  # Extract data blocks
  blocks <- hyperdesign$blocks
  n_domains <- length(blocks)
  
  if (n_domains < 2) {
    stop("FPGW requires at least 2 domains", call. = FALSE)
  }
  
  # Extract data matrices and ensure they're in the right format
  X_list <- lapply(blocks, function(block) {
    if (inherits(block, "multidesign")) {
      as.matrix(block$X)
    } else if (is.matrix(block) || is.data.frame(block)) {
      as.matrix(block)
    } else {
      stop("Each block must be a multidesign object, matrix, or data frame", 
           call. = FALSE)
    }
  })
  
  # Get sample sizes
  n_samples <- sapply(X_list, nrow)
  
  # Get domain names
  domain_names <- names(blocks)
  if (is.null(domain_names)) {
    domain_names <- paste0("Domain", seq_len(n_domains))
  }
  
  list(
    X_list = X_list,
    n_samples = n_samples,
    n_domains = n_domains,
    domain_names = domain_names
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
    for (i in 1:n1) {
      for (j in 1:n2) {
        C[i, j] <- sum(abs(X1[i, ] - X2[j, ]))
      }
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