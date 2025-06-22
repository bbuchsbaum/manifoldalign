#' Shared Utilities for Optimal Transport Methods
#' 
#' Internal functions shared between gromov_wasserstein and fpgw implementations.
#' @keywords internal
#' @name utils-ot
NULL

#' Compute Distance Matrix
#' 
#' Memory-efficient distance computation that avoids creating full distance matrix
#' @param X Data matrix (samples x features)
#' @param metric Distance metric to use (default: "euclidean")
#' @param chunk_size Chunk size for memory-efficient computation
#' @keywords internal
compute_distance_matrix <- function(X, metric = "euclidean", chunk_size = 1000) {
  n <- nrow(X)
  
  # For small matrices, use standard approach
  if (n <= chunk_size || n * n * 8 < 1e8) { # < 100MB
    D <- as.matrix(dist(X, method = metric))
    return(D)
  }
  
  # For large matrices, compute in chunks
  D <- matrix(0, n, n)
  
  # Process in chunks to control memory usage
  for (i in seq(1, n, by = chunk_size)) {
    i_end <- min(i + chunk_size - 1, n)
    for (j in seq(i, n, by = chunk_size)) {
      j_end <- min(j + chunk_size - 1, n)
      
      # Compute distance block
      if (metric == "euclidean") {
        # Optimized Euclidean distance
        Xi <- X[i:i_end, , drop = FALSE]
        Xj <- X[j:j_end, , drop = FALSE]
        
        # ||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 * x_i' * x_j
        Xi_norm2 <- rowSums(Xi^2)
        Xj_norm2 <- rowSums(Xj^2)
        XiXj <- tcrossprod(Xi, Xj)
        
        D_block <- outer(Xi_norm2, rep(1, length(Xj_norm2))) +
                   outer(rep(1, length(Xi_norm2)), Xj_norm2) -
                   2 * XiXj
        D_block <- sqrt(pmax(D_block, 0))
      } else {
        # Safe computation for other metrics
        if (i == j) {
          # Diagonal block - use standard dist
          D_block <- as.matrix(dist(X[i:i_end, , drop = FALSE], method = metric))
        } else {
          # Off-diagonal block - compute pairwise distances correctly
          Xi <- X[i:i_end, , drop = FALSE]
          Xj <- X[j:j_end, , drop = FALSE]
          ni <- nrow(Xi)
          nj <- nrow(Xj)
          
          # Compute distances between all pairs
          D_block <- matrix(0, ni, nj)
          
          # If proxy package is available, use it for efficiency
          if (requireNamespace("proxy", quietly = TRUE)) {
            D_block <- as.matrix(proxy::dist(Xi, Xj, method = metric))
          } else {
            # Manual computation as fallback
            for (ii in seq_len(ni)) {
              for (jj in seq_len(nj)) {
                # Create temporary matrix and compute distance
                pair <- rbind(Xi[ii, , drop = FALSE], Xj[jj, , drop = FALSE])
                D_block[ii, jj] <- as.matrix(dist(pair, method = metric))[1, 2]
              }
            }
          }
        }
      }
      
      # Fill in the matrix
      D[i:i_end, j:j_end] <- D_block
      if (i != j) {
        D[j:j_end, i:i_end] <- t(D_block)
      }
    }
  }
  
  D
}

#' Validate marginal distributions
#' @keywords internal
validate_marginals <- function(marginals, n_samples, n_domains) {
  if (is.null(marginals)) {
    # Default to uniform marginals
    marginals <- lapply(n_samples, function(n) rep(1/n, n))
  } else {
    # Validate provided marginals
    if (!is.list(marginals) || length(marginals) != n_domains) {
      stop("marginals must be a list of length ", n_domains, call. = FALSE)
    }
    
    # Check and normalize marginals
    for (i in seq_len(n_domains)) {
      if (length(marginals[[i]]) != n_samples[i]) {
        stop("marginals[[", i, "]] must have length ", n_samples[i], call. = FALSE)
      }
      
      # Ensure non-negative
      if (any(marginals[[i]] < 0)) {
        stop("marginals must be non-negative", call. = FALSE)
      }
      
      # Check for zero masses
      if (any(marginals[[i]] == 0)) {
        stop("marginals must be strictly positive (no zero masses)", call. = FALSE)
      }
      
      # Normalize to sum to 1
      marg_sum <- sum(marginals[[i]])
      if (abs(marg_sum - 1) > 1e-10) {
        marginals[[i]] <- marginals[[i]] / marg_sum
      }
    }
  }
  
  marginals
}

#' Extract data from hyperdesign and apply preprocessing
#' @keywords internal
prepare_ot_data <- function(data, preproc = NULL) {
  if (!inherits(data, "hyperdesign")) {
    stop("Input must be a hyperdesign object", call. = FALSE)
  }
  
  n_domains <- length(data)
  if (n_domains < 2) {
    stop("Optimal transport requires at least 2 domains", call. = FALSE)
  }
  
  # Apply preprocessing if provided
  if (!is.null(preproc)) {
    data <- multivarious::init_transform(data, preproc)
  }
  
  # Extract data matrices
  X_list <- lapply(data, function(d) d$x)
  n_samples <- sapply(X_list, nrow)
  
  # Get domain names with fallback
  domain_names <- names(data)
  if (is.null(domain_names)) {
    domain_names <- paste0("domain", seq_len(n_domains))
  }
  
  list(
    X_list = X_list,
    n_samples = n_samples,
    n_domains = n_domains,
    domain_names = domain_names,
    data = data
  )
}

#' Compute feature cost matrix between two domains
#' @keywords internal
compute_feature_cost <- function(X1, X2, metric = "euclidean") {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  if (metric == "euclidean") {
    # Efficient computation for Euclidean distance
    # Handle case where X1 and X2 have different number of columns
    if (ncol(X1) != ncol(X2)) {
      # Compute pairwise distances manually
      D <- matrix(0, n1, n2)
      for (i in seq_len(n1)) {
        for (j in seq_len(n2)) {
          # This will give NA if dimensions don't match, which is what we want
          D[i,j] <- sqrt(sum((X1[i,] - X2[j,])^2, na.rm = TRUE))
        }
      }
      return(D)
    }
    
    # Same dimensions - use efficient computation
    as.matrix(dist(rbind(X1, X2)))[seq_len(n1), (n1 + 1):(n1 + n2)]
  } else {
    # Use compute_distance_matrix for other metrics
    if (ncol(X1) != ncol(X2)) {
      stop("Different number of features requires euclidean metric", call. = FALSE)
    }
    D_full <- compute_distance_matrix(rbind(X1, X2), metric = metric)
    D_full[seq_len(n1), (n1 + 1):(n1 + n2)]
  }
}