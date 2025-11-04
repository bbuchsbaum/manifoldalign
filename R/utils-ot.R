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
    X1_sq <- rowSums(X1^2)
    X2_sq <- rowSums(X2^2)
    XY <- tcrossprod(X1, X2)
    C <- outer(X1_sq, rep(1, n2)) + outer(rep(1, n1), X2_sq) - 2 * XY
    C[C < 0] <- 0
    sqrt(C)
  } else if (metric == "manhattan") {
    C <- matrix(0, n1, n2)
    for (i in seq_len(n1)) {
      diffs <- abs(matrix(X1[i, ], n2, ncol(X1), byrow = TRUE) - X2)
      C[i, ] <- rowSums(diffs)
    }
    C
  } else if (metric == "cosine") {
    X1_norm <- X1 / sqrt(rowSums(X1^2) + 1e-10)
    X2_norm <- X2 / sqrt(rowSums(X2^2) + 1e-10)
    1 - tcrossprod(X1_norm, X2_norm)
  } else {
    stop("Unsupported metric: ", metric, call. = FALSE)
  }
}

#' Normalize transport plan rows to probability distributions
#' @keywords internal
normalize_plan_rows <- function(plan) {
  if (!is.matrix(plan)) {
    plan <- as.matrix(plan)
  }
  if (nrow(plan) == 0 || ncol(plan) == 0) {
    return(plan)
  }
  row_sums <- rowSums(plan)
  plan_norm <- plan
  positive <- row_sums > 0
  if (any(positive)) {
    plan_norm[positive, ] <- plan_norm[positive, , drop = FALSE] / row_sums[positive]
  }
  if (any(!positive)) {
    plan_norm[!positive, ] <- 1 / ncol(plan_norm)
  }
  plan_norm
}

#' Compute K-nearest neighbours under supported metrics
#' @keywords internal
find_knn_indices <- function(query, data, k = 5, metric = "euclidean") {
  if (!is.matrix(query)) query <- as.matrix(query)
  if (!is.matrix(data)) data <- as.matrix(data)
  if (nrow(data) == 0) {
    stop("Data matrix must contain at least one sample", call. = FALSE)
  }
  k <- max(1L, min(as.integer(k), nrow(data)))
  n_query <- nrow(query)
  idx <- matrix(NA_integer_, n_query, k)
  dists <- matrix(NA_real_, n_query, k)

  data_norm <- NULL
  if (metric == "cosine") {
    data_norm <- sqrt(rowSums(data^2))
    data_norm[data_norm < 1e-12] <- 1e-12
  }

  for (i in seq_len(n_query)) {
    d_vec <- compute_point_distance(query[i, , drop = FALSE], data, metric, data_norm)
    ord <- order(d_vec, na.last = NA)[seq_len(k)]
    idx[i, ] <- ord
    dists[i, ] <- d_vec[ord]
  }

  list(idx = idx, dists = dists)
}

#' Compute distances between a point and all rows of a matrix
#' @keywords internal
compute_point_distance <- function(point, data, metric, data_norm = NULL) {
  if (metric == "euclidean") {
    diffs <- sweep(data, 2, point, FUN = "-")
    sqrt(rowSums(diffs^2))
  } else if (metric == "manhattan") {
    diffs <- sweep(data, 2, point, FUN = "-")
    rowSums(abs(diffs))
  } else if (metric == "cosine") {
    point_norm <- sqrt(sum(point^2))
    if (point_norm < 1e-12) {
      point_norm <- 1e-12
    }
    sims <- (data %*% as.numeric(point)) / (data_norm * point_norm)
    sims <- pmin(pmax(sims, -1), 1)
    as.numeric(1 - sims)
  } else {
    stop("Unsupported metric for neighbour search: ", metric, call. = FALSE)
  }
}

#' Combine neighbour plans into barycentric weights
#' @keywords internal
aggregate_barycentric_weights <- function(plan_rows, neighbor_idx, neighbor_dists) {
  n_query <- nrow(neighbor_idx)
  n_target <- ncol(plan_rows)
  weights <- matrix(0, n_query, n_target)

  for (i in seq_len(n_query)) {
    idx <- neighbor_idx[i, ]
    dists <- neighbor_dists[i, ]
    if (any(is.na(idx))) {
      next
    }
    if (any(dists < 1e-12)) {
      mask <- dists < 1e-12
      local_w <- mask / sum(mask)
    } else {
      local_w <- 1 / (dists + 1e-8)
      local_w <- local_w / sum(local_w)
    }
    combined <- local_w %*% plan_rows[idx, , drop = FALSE]
    total <- sum(combined)
    if (total > 0) {
      combined <- combined / total
    } else {
      combined <- rep(1 / n_target, n_target)
    }
    weights[i, ] <- combined
  }
  weights
}

#' Resolve domain input (index or name)
#' @keywords internal
resolve_domain_index <- function(x, domain_names) {
  if (is.character(x)) {
    idx <- match(x, domain_names)
    if (is.na(idx)) {
      stop("Domain '", x, "' not found", call. = FALSE)
    }
    idx
  } else {
    idx <- as.integer(x)
    if (length(idx) != 1 || is.na(idx) || idx < 1 || idx > length(domain_names)) {
      stop("Domain index must be between 1 and ", length(domain_names), call. = FALSE)
    }
    idx
  }
}

#' Compute linear index for upper-triangular pair
#' @keywords internal
get_plan_index <- function(n_domains, from_idx, to_idx) {
  if (from_idx >= to_idx) {
    stop("from_idx must be strictly less than to_idx", call. = FALSE)
  }
  idx <- (from_idx - 1) * n_domains - (from_idx * (from_idx - 1)) / 2 + (to_idx - from_idx)
  as.integer(idx)
}

#' Extract transport plan for an ordered domain pair
#' @keywords internal
extract_transport_plan <- function(plans, from_idx, to_idx, n_domains) {
  if (from_idx == to_idx) {
    stop("Source and target domains must be different", call. = FALSE)
  }
  if (from_idx < to_idx) {
    idx <- get_plan_index(n_domains, from_idx, to_idx)
    plan <- plans[[idx]]
    if (is.null(plan)) {
      stop("Transport plan not found for specified domain pair", call. = FALSE)
    }
    plan
  } else {
    t(extract_transport_plan(plans, to_idx, from_idx, n_domains))
  }
}
