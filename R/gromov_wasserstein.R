# Generic is defined in all_generic.R
# Note: compute_distance_matrix has been moved to utils-ot.R

#' @rdname gromov_wasserstein
#' @method gromov_wasserstein hyperdesign
#' @export
gromov_wasserstein.hyperdesign <- function(data,
                                         epsilon = 0.1,
                                         metric = c("euclidean", "maximum", "manhattan", 
                                                   "canberra", "binary", "minkowski"),
                                         init = "uniform",
                                         max_iter = 100,
                                         tol = 1e-9,
                                         verbose = FALSE,
                                         marginals = NULL,
                                         scale = c("per_domain", "global", "none"),
                                         inner_max_iter = 30,
                                         inner_tol = 1e-9,
                                         dist_chunk_size = 1000,
                                         ...) {
  
  # Input validation
  if (!inherits(data, "hyperdesign")) {
    stop("gromov_wasserstein requires a hyperdesign object", call. = FALSE)
  }
  
  metric <- match.arg(metric)
  scale <- match.arg(scale)
  
  if (epsilon <= 0) {
    stop("epsilon must be positive", call. = FALSE)
  }
  
  # Extract data matrices
  n_domains <- length(data)
  if (n_domains < 2) {
    stop("Gromov-Wasserstein requires at least 2 domains", call. = FALSE)
  }
  
  X_list <- lapply(data, function(d) d$x)
  n_samples <- sapply(X_list, nrow)
  
  # Get domain names with fallback
  domain_names <- names(data)
  if (is.null(domain_names)) {
    domain_names <- paste0("domain", seq_len(n_domains))
  }
  
  # Handle marginals
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
  
  # Compute distance matrices for each domain
  if (verbose) message("Computing distance matrices...")
  
  C_list <- lapply(X_list, function(X) {
    compute_distance_matrix(X, metric = metric, chunk_size = dist_chunk_size)
  })
  
  # Apply scaling based on user preference
  if (scale == "per_domain") {
    # Scale each domain independently
    C_list <- lapply(C_list, function(D) {
      max_val <- max(D)
      if (max_val > 0) D / max_val else D
    })
  } else if (scale == "global") {
    # Find global maximum across all domains
    global_max <- max(sapply(C_list, max))
    if (global_max > 0) {
      C_list <- lapply(C_list, function(D) D / global_max)
    }
  }
  # else scale == "none", keep original distances
  
  # Initialize results
  n_pairs <- n_domains * (n_domains - 1) / 2
  transport_plans <- vector("list", n_pairs)
  distances <- matrix(0, n_domains, n_domains)
  losses <- matrix(NA, n_domains, n_domains)
  converged <- matrix(TRUE, n_domains, n_domains)
  
  # Compute pairwise Gromov-Wasserstein distances
  pair_idx <- 1
  for (i in 1:(n_domains - 1)) {
    for (j in (i + 1):n_domains) {
      if (verbose) {
        message(sprintf("\nComputing GW distance between domain %d and %d", i, j))
      }
      
      # Get distance matrices
      C1 <- C_list[[i]]
      C2 <- C_list[[j]]
      n1 <- n_samples[i]
      n2 <- n_samples[j]
      
      # Get marginals for this pair
      p <- marginals[[i]]
      q <- marginals[[j]]
      
      # Initialize transport plan
      if (is.character(init)) {
        if (init == "uniform") {
          P <- matrix(1 / (n1 * n2), n1, n2)
        } else if (init == "random") {
          P <- matrix(runif(n1 * n2), n1, n2)
          P <- P / sum(P)
        } else {
          stop("Unknown init method: ", init, call. = FALSE)
        }
      } else if (is.list(init)) {
        if (pair_idx > length(init)) {
          stop("init list has ", length(init), " elements but pair ", 
               pair_idx, " was requested", call. = FALSE)
        }
        P <- init[[pair_idx]]
        if (!is.matrix(P) || nrow(P) != n1 || ncol(P) != n2) {
          stop("init[[", pair_idx, "]] must be a ", n1, "x", n2, " matrix", 
               call. = FALSE)
        }
      } else {
        stop("init must be 'uniform', 'random', or a list of matrices", call. = FALSE)
      }
      
      # Solve entropic Gromov-Wasserstein
      result <- solve_entropic_gw(
        C1, C2, p, q, P,
        epsilon = epsilon,
        max_iter = max_iter,
        tol = tol,
        inner_max_iter = inner_max_iter,
        inner_tol = inner_tol,
        verbose = verbose
      )
      
      # Store results
      transport_plans[[pair_idx]] <- result$P
      distances[i, j] <- distances[j, i] <- result$distance
      losses[i, j] <- losses[j, i] <- result$loss
      converged[i, j] <- converged[j, i] <- result$converged
      
      pair_idx <- pair_idx + 1
    }
  }
  
  # Create return object
  structure(
    list(
      transport_plans = transport_plans,
      distances = distances,
      cost_matrices = C_list,
      losses = losses,
      converged = converged,
      epsilon = epsilon,
      loss_function = "square_loss",
      n_samples = n_samples,
      domain_names = domain_names
    ),
    class = c("gromov_wasserstein", "multiblock_biprojector")
  )
}

# Simplified and stable entropic GW solver
solve_entropic_gw <- function(C1, C2, p, q, P_init,
                            epsilon = 0.1,
                            max_iter = 100,
                            tol = 1e-9,
                            inner_max_iter = 30,
                            inner_tol = 1e-9,
                            verbose = FALSE) {
  
  n1 <- nrow(C1)
  n2 <- nrow(C2)
  
  # Initialize
  P <- P_init
  
  # Precompute constants that don't change during iterations
  C1_sqr <- C1^2
  C2_sqr <- C2^2
  
  # These terms only depend on marginals, not on P
  # Pre-compute them once outside the loop
  constC1 <- as.vector(C1_sqr %*% p)  # n1 x 1 vector
  constC2 <- as.vector(C2_sqr %*% q)  # n2 x 1 vector
  
  # Pre-compute the constant part of the linear term
  const_term <- outer(constC1, rep(1, n2)) + outer(rep(1, n1), constC2)
  
  loss_old <- Inf
  consecutive_inf_count <- 0  # Track consecutive Inf losses
  
  for (iter in 1:max_iter) {
    # Compute the linear term for kernel computation
    # Only the P-dependent part changes: -2 * C1 * P * C2^T
    linear_term <- const_term - 2 * C1 %*% P %*% t(C2)
    
    # Ensure non-negative
    linear_term <- pmax(linear_term, 0)
    
    # Numerical stabilization: shift before exp to prevent overflow/underflow
    max_val <- max(linear_term)
    stable_linear_term <- linear_term - max_val
    
    # Sinkhorn iterations with stabilized kernel
    K <- exp(-stable_linear_term / epsilon)
    
    # Update P using Sinkhorn scaling
    u <- p
    v <- q
    
    # Inner Sinkhorn iterations with configurable limit
    u_old <- u
    v_old <- v
    
    for (sinkhorn_iter in 1:inner_max_iter) {
      u <- p / (K %*% v + 1e-20)
      v <- q / (t(K) %*% u + 1e-20)
      
      # Check for numerical issues and break if detected
      if (any(!is.finite(u)) || any(!is.finite(v))) {
        if (verbose) {
          message(sprintf("  Warning: Numerical instability at iter %d, sinkhorn_iter %d", 
                         iter, sinkhorn_iter))
        }
        # Reset to safe values
        u <- p
        v <- q
        break
      }
      
      # Check inner loop convergence
      if (inner_tol > 0) {
        u_change <- max(abs(u - u_old) / (abs(u_old) + 1e-15))
        v_change <- max(abs(v - v_old) / (abs(v_old) + 1e-15))
        if (max(u_change, v_change) < inner_tol) {
          break
        }
        u_old <- u
        v_old <- v
      }
    }
    
    # Update transport plan
    P <- sweep(sweep(K, 1, u, "*"), 2, v, "*")
    
    # Check for numerical underflow
    P_sum <- sum(P)
    if (!is.finite(P_sum) || P_sum < 1e-300) {
      # Numerical underflow or NA - reset to product of marginals
      if (verbose) {
        message(sprintf("  Warning: Transport plan %s at iter %d, resetting", 
                       ifelse(is.na(P_sum), "NA", "underflow"), iter))
      }
      P <- outer(p, q)  # Reset to product of marginals
    }
    
    # CRITICAL FIX: Recompute linear term with updated P for loss calculation
    # This ensures the loss reflects the current transport plan
    linear_term_updated <- const_term - 2 * C1 %*% P %*% t(C2)
    linear_term_updated <- pmax(linear_term_updated, 0)
    
    # Compute loss with the updated linear term
    loss_new <- sum(P * linear_term_updated) - epsilon * sum(P * log(P + 1e-20))
    
    # Enhanced handling of Inf loss cases
    if (!is.finite(loss_new)) {
      consecutive_inf_count <- consecutive_inf_count + 1
      
      if (verbose) {
        message(sprintf("  Warning: Non-finite loss at iter %d (count: %d)", 
                       iter, consecutive_inf_count))
      }
      
      # If we get too many consecutive Inf losses, terminate early
      if (consecutive_inf_count >= 3) {
        if (verbose) {
          message("  Stopping due to repeated non-finite losses")
        }
        return(list(
          P = P,
          distance = Inf,  # Use Inf instead of NA for symmetry
          loss = Inf,
          converged = FALSE,
          iterations = iter
        ))
      }
      
      loss_new <- Inf
    } else {
      consecutive_inf_count <- 0  # Reset counter on valid loss
    }
    
    if (verbose && iter %% 10 == 0) {
      message(sprintf("  Iteration %3d: loss = %.6e", iter, loss_new))
    }
    
    # Check convergence - only if both losses are finite
    if (is.finite(loss_new) && is.finite(loss_old)) {
      rel_change <- abs(loss_old - loss_new) / (abs(loss_old) + 1e-15)
      if (is.finite(rel_change) && rel_change < tol) {
        if (verbose) message(sprintf("  Converged at iteration %d", iter))
        return(list(
          P = P,
          distance = sqrt(max(0, sum(P * linear_term_updated))),
          loss = loss_new,
          converged = TRUE,
          iterations = iter
        ))
      }
    } else if (!is.finite(loss_new) && !is.finite(loss_old)) {
      # Both losses are Inf - check if we should stop
      if (iter > 5) {
        if (verbose) message("  Stopping: unable to achieve finite loss")
        return(list(
          P = P,
          distance = Inf,  # Use Inf instead of NA for symmetry
          loss = Inf,
          converged = FALSE,
          iterations = iter
        ))
      }
    }
    
    loss_old <- loss_new
  }
  
  if (verbose) message("  Did not converge within max_iter iterations")
  
  # Final computation of linear term for accurate distance
  linear_term_final <- const_term - 2 * C1 %*% P %*% t(C2)
  linear_term_final <- pmax(linear_term_final, 0)
  
  list(
    P = P,
    distance = sqrt(max(0, sum(P * linear_term_final))),
    loss = loss_new,
    converged = FALSE,
    iterations = max_iter
  )
}

print.gromov_wasserstein <- function(x, ...) {
  cat("Gromov-Wasserstein Alignment\n")
  cat("============================\n")
  cat("Number of domains:", length(x$domain_names), "\n")
  cat("Domain names:", paste(x$domain_names, collapse = ", "), "\n")
  cat("Samples per domain:", paste(x$n_samples, collapse = ", "), "\n")
  cat("\nParameters:\n")
  cat("  Loss function:", x$loss_function, "\n")
  cat("  Epsilon (regularization):", x$epsilon, "\n")
  cat("\nPairwise GW distances:\n")
  print(round(x$distances, 4))
  
  if (!all(x$converged)) {
    cat("\nWarning: Some optimizations did not converge\n")
  }
  
  invisible(x)
}

#' Predict method for Gromov-Wasserstein
#' 
#' Transport new samples using learned GW alignment
#' 
#' @param object A gromov_wasserstein object
#' @param newdata New data from one domain to transport to another
#' @param from Source domain index or name
#' @param to Target domain index or name
#' @param ... Additional arguments
#' 
#' @return Transported samples in the target domain space
#' @export
predict.gromov_wasserstein <- function(object, newdata, from, to, ...) {
  # This would implement out-of-sample extension
  # For now, return a placeholder
  stop("Out-of-sample prediction not yet implemented for Gromov-Wasserstein", 
       call. = FALSE)
}

#' @rdname gromov_wasserstein
#' @method gromov_wasserstein default
#' @export
gromov_wasserstein.default <- function(data, ...) {
  stop("gromov_wasserstein requires a hyperdesign object", call. = FALSE)
}