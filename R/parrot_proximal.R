#' Compute PARROT objective
#'
#' Computes the full objective:
#' \eqn{<C_rwr, S> + lambda_e * L_e(S) + lambda_n * L_n(S) + lambda_p * L_p(S)}.
#'
#' @param S Transport plan matrix
#' @param C_rwr Position-aware cost matrix
#' @param networks List of network structures  
#' @param anchor_info Anchor information
#' @param lambda_e Edge consistency weight
#' @param lambda_n Neighborhood consistency weight
#' @param lambda_p Anchor prior weight
#' @return Objective function value
#' @keywords internal
compute_parrot_objective <- function(S, C_rwr, networks, anchor_info, 
                                   lambda_e, lambda_n, lambda_p) {
  
  # Base transport cost
  obj <- sum(C_rwr * S)

  # Edge consistency regularizer L_e(S) from Eq. 8
  if (lambda_e > 0) {
    row_normalize <- function(M) {
      norms <- sqrt(rowSums(M^2))
      norms[norms < 1e-12] <- 1e-12
      M / norms
    }

    A1 <- networks[[1]]$adjacency
    A2 <- networks[[2]]$adjacency
    if (inherits(A1, "Matrix")) A1 <- as.matrix(A1)
    if (inherits(A2, "Matrix")) A2 <- as.matrix(A2)

    X1n <- row_normalize(networks[[1]]$features)
    X2n <- row_normalize(networks[[2]]$features)

    C1 <- (1 - X1n %*% t(X1n)) * (A1 != 0)
    C2 <- (1 - X2n %*% t(X2n)) * (A2 != 0)

    C1_sq <- rowSums(C1^2)
    C2_sq <- rowSums(C2^2)

    L_base <- outer(C1_sq, rep(1, length(C2_sq))) +
      outer(rep(1, length(C1_sq)), C2_sq) -
      2 * (C1 %*% S %*% t(C2))

    obj <- obj + lambda_e * sum(L_base * S)
  }

  
  # Neighborhood consistency regularizer L_n(S) from Eq. 9
  if (lambda_n > 0) {
    W1 <- networks[[1]]$transition
    W2 <- networks[[2]]$transition

    S_hat <- t(W1) %*% S %*% W2

    eps <- 1e-12
    S_safe <- S + eps
    S_hat_safe <- S_hat + eps

    kl_div <- sum(S_safe * (log(S_safe) - log(S_hat_safe)))
    obj <- obj + lambda_n * kl_div
  }

  # Anchor prior regularizer L_p(S)
  if (lambda_p > 0) {
    n1 <- nrow(S)
    n2 <- ncol(S)
    eps <- 1e-12
    H <- matrix(eps, n1, n2)

    if (length(anchor_info$idx1) > 0 && length(anchor_info$idx2) > 0) {
      for (i in seq_along(anchor_info$idx1)) {
        idx1 <- anchor_info$idx1[i]
        anchor_val <- anchor_info$vec1[idx1]

        for (j in seq_along(anchor_info$idx2)) {
          idx2 <- anchor_info$idx2[j]
          if (anchor_info$vec2[idx2] == anchor_val) {
            H[idx1, idx2] <- 1.0
          }
        }
      }
    } else {
      H[] <- 1.0
    }

    H <- H / sum(H)
    S_safe <- S + eps
    obj <- obj + lambda_p * sum(S_safe * (log(S_safe) - log(H)))
  }
  
  obj
}

#' Compute Gradient of Edge Consistency Regularizer (R implementation)
#' 
#' @param S Current transport plan
#' @param networks Network structures
#' @return Gradient matrix
#' @keywords internal
compute_edge_gradient_r <- function(S, networks) {
  A1 <- networks[[1]]$adjacency
  A2 <- networks[[2]]$adjacency
  X1 <- networks[[1]]$features
  X2 <- networks[[2]]$features
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  # Initialize gradient
  grad <- matrix(0, n1, n2)
  
  # Compute pairwise edge distances
  edges1 <- which(A1 != 0, arr.ind = TRUE)
  edges1 <- edges1[edges1[,1] < edges1[,2], ]
  
  if (nrow(edges1) > 0) {
    # Pre-compute all edge distances in G1
    d1_edges <- numeric(nrow(edges1))
    for (e in 1:nrow(edges1)) {
      i <- edges1[e, 1]
      j <- edges1[e, 2]
      d1_edges[e] <- sum((X1[i,] - X1[j,])^2)
    }
    
    # Pre-compute all edge distances in G2
    edges2 <- which(A2 != 0, arr.ind = TRUE)
    edges2 <- edges2[edges2[,1] < edges2[,2], ]
    
    if (nrow(edges2) > 0) {
      d2_edges <- numeric(nrow(edges2))
      for (e in 1:nrow(edges2)) {
        a <- edges2[e, 1]
        b <- edges2[e, 2]
        d2_edges[e] <- sum((X2[a,] - X2[b,])^2)
      }
      
      # Compute gradient contributions
      for (e1 in 1:nrow(edges1)) {
        i <- edges1[e1, 1]
        j <- edges1[e1, 2]
        d1 <- d1_edges[e1]
        
        for (e2 in 1:nrow(edges2)) {
          a <- edges2[e2, 1]
          b <- edges2[e2, 2]
          d2 <- d2_edges[e2]
          
          diff <- 2 * (d1 - d2)
          grad[i, a] <- grad[i, a] + S[j, b] * diff
          grad[j, b] <- grad[j, b] + S[i, a] * diff
        }
      }
    }
  }
  
  grad
}

#' Compute Gradient of Edge Consistency Regularizer (Dispatcher)
#' 
#' @param S Current transport plan
#' @param networks Network structures
#' @param use_cpp Whether to use C++ implementation
#' @return Gradient matrix
#' @keywords internal
compute_edge_gradient <- function(S, networks, use_cpp = FALSE) {
  
  if (isTRUE(use_cpp)) {
    # Extract components for C++ function
    A1 <- networks[[1]]$adjacency
    A2 <- networks[[2]]$adjacency
    X1 <- networks[[1]]$features
    X2 <- networks[[2]]$features
    
    # Ensure matrices are in correct format
    if (!inherits(A1, "dgCMatrix")) A1 <- as(A1, "dgCMatrix")
    if (!inherits(A2, "dgCMatrix")) A2 <- as(A2, "dgCMatrix")
    if (inherits(X1, "Matrix")) X1 <- as.matrix(X1)
    if (inherits(X2, "Matrix")) X2 <- as.matrix(X2)
    if (inherits(S, "Matrix")) S <- as.matrix(S)
    
    compute_edge_gradient_cpp_fallback(S, A1, A2, X1, X2)
  } else {
    # Use existing R implementation
    compute_edge_gradient_r(S, networks)
  }
}

#' Solve PARROT with Proper Proximal Point Method
#' 
#' Implements Algorithm 2 from the paper with monotonicity checks
#' 
#' @inheritParams solve_parrot_transport
#' @return List with transport plan and convergence info
#' @keywords internal
solve_parrot_proximal <- function(networks, rwr_features, anchor_info,
                                 lambda_e, lambda_n, lambda_p, 
                                 tau, alpha, sigma, gamma, solver,
                                 max_iter, max_outer = 10, tol, use_cpp = FALSE) {

  n1 <- nrow(networks[[1]]$features)
  n2 <- nrow(networks[[2]]$features)

  build_anchor_prior <- function() {
    eps <- 1e-12
    H <- matrix(eps, n1, n2)

    if (length(anchor_info$idx1) > 0 && length(anchor_info$idx2) > 0) {
      for (i in seq_along(anchor_info$idx1)) {
        idx1 <- anchor_info$idx1[i]
        anchor_val <- anchor_info$vec1[idx1]

        for (j in seq_along(anchor_info$idx2)) {
          idx2 <- anchor_info$idx2[j]
          if (anchor_info$vec2[idx2] == anchor_val) {
            H[idx1, idx2] <- 1.0
          }
        }
      }
    } else {
      H[] <- 1.0
    }

    H / sum(H)
  }

  eps <- 1e-12
  W1 <- networks[[1]]$transition
  W2 <- networks[[2]]$transition
  H_prior <- build_anchor_prior()

  C_rwr <- compute_parrot_cost(networks, rwr_features, anchor_info, alpha, sigma, gamma, use_cpp)

  epsilon <- tau + lambda_n + lambda_p
  if (epsilon <= 0) {
    epsilon <- max(tau, 1e-3)
  }

  if (solver != "sinkhorn") {
    stop("Proximal solver currently supports sinkhorn solver only")
  }

  S <- solve_sinkhorn_stabilized(C_rwr, epsilon, max_iter, tol, use_cpp)
  C_eff_prev <- C_rwr
  phi_prev <- sum(C_eff_prev * S)
  obj_prev <- compute_parrot_objective(S, C_rwr, networks, anchor_info,
                                      lambda_e, lambda_n, lambda_p)

  obj_history <- numeric(max_outer)

  for (t in seq_len(max_outer)) {
    S_prev <- S

    C_candidate <- C_rwr

    if (lambda_e > 0) {
      L_edge <- compute_edge_consistency(networks, S_prev, 1)
      C_candidate <- C_candidate + lambda_e * L_edge
    }

    if (lambda_n > 0) {
      Shat_prev <- t(W1) %*% S_prev %*% W2
      C_candidate <- C_candidate - lambda_n * log(Shat_prev + eps)
    }

    if (lambda_p > 0) {
      C_candidate <- C_candidate - lambda_p * log(H_prior + eps)
    }

    if (tau > 0) {
      C_candidate <- C_candidate - tau * log(S_prev + eps)
    }

    phi_candidate <- sum(C_candidate * S_prev)

    if (phi_candidate <= phi_prev) {
      C_eff <- C_candidate
      C_eff_prev <- C_candidate
      phi_prev <- phi_candidate
    } else {
      C_eff <- C_eff_prev
    }

    min_cost <- min(C_eff)
    if (min_cost <= 0) {
      C_eff <- C_eff - min_cost + eps
      C_eff_prev <- C_eff
    }

    S <- solve_sinkhorn_stabilized(C_eff, epsilon, max_iter, tol, use_cpp)
    phi_prev <- sum(C_eff_prev * S)

    obj_prev <- compute_parrot_objective(S, C_rwr, networks, anchor_info,
                                        lambda_e, lambda_n, lambda_p)
    obj_history[t] <- obj_prev

    if (max(abs(S - S_prev)) < tol) {
      obj_history <- obj_history[seq_len(t)]
      break
    }

    if (t == max_outer) {
      obj_history <- obj_history[seq_len(t)]
    }
  }

  list(
    transport_plan = S,
    cost_matrix = C_eff_prev,
    objective_history = obj_history,
    outer_iterations = length(obj_history),
    final_objective = tail(obj_history, 1)
  )
}


#' Log-domain Stabilized Sinkhorn Algorithm (R implementation)
#' 
#' Implements numerically stable Sinkhorn in log domain.
#' For a square matrix (n1=n2), it produces a doubly-stochastic matrix (row/col sums = 1).
#' For a rectangular matrix, it produces a transport plan with uniform probability marginals.
#' 
#' @param C Cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Transport plan matrix
#' @keywords internal
solve_sinkhorn_stabilized_r <- function(C, tau, max_iter, tol) {
  n1 <- nrow(C)
  n2 <- ncol(C)
  
  # Normalize cost matrix to [0, 1] range for numerical stability
  C_min <- min(C)
  C_max <- max(C)
  if (C_max - C_min > 1e-9) {
    C <- (C - C_min) / (C_max - C_min)
  } else {
    # If all costs are the same, return uniform transport
    return(matrix(1/(n1*n2), n1, n2))
  }
  
  # Uniform probability marginals (Eq. 17--18)
  mu <- rep(1 / n1, n1)
  nu <- rep(1 / n2, n2)
  
  # Initialize log-potentials
  log_u <- rep(0, n1)
  log_v <- rep(0, n2)
  
  # Pre-compute log-kernel
  log_K <- -C / tau
  
  # Stabilize using the log-sum-exp trick
  max_log_K <- max(log_K)
  log_K <- log_K - max_log_K
  
  # More iterations needed for small tau values
  actual_max_iter <- max(max_iter, min(1000, 100 / tau))
  
  for (iter in 1:actual_max_iter) {
    log_u_prev <- log_u
    
    # Update log_v: log(v) = log(nu) - log(K^T * exp(log_u))
    # This is the log-sum-exp operation across rows of log_K
    log_sum_u_K <- apply(log_K, 2, function(col) {
      max_val <- max(log_u + col)
      max_val + log(sum(exp(log_u + col - max_val)))
    })
    log_v <- log(nu) - log_sum_u_K
    
    # Update log_u: log(u) = log(mu) - log(K * exp(log_v))
    # This is the log-sum-exp operation across columns of log_K
    log_sum_v_KT <- apply(log_K, 1, function(row) {
      max_val <- max(log_v + row)
      max_val + log(sum(exp(log_v + row - max_val)))
    })
    log_u <- log(mu) - log_sum_v_KT
    
    # Check for convergence in the potentials
    if (max(abs(log_u - log_u_prev)) < tol) {
      break
    }
  }
  
  # Reconstruct the transport plan: S = u * K * v
  # which in log-space is log(S) = log(u) + log(K) + log(v)
  log_S <- outer(log_u, log_v, "+") + log_K
  
  # Convert back from log domain. This is the correct OT solution.
  S <- exp(log_S)
  
  # --- VERIFICATION: Check doubly stochastic constraints ---
  # Verify that the Sinkhorn algorithm produced correct marginals
  row_sums <- rowSums(S)
  col_sums <- colSums(S)
  
  # Check if constraints are satisfied within reasonable tolerance
  row_error <- max(abs(row_sums - mu))
  col_error <- max(abs(col_sums - nu))
  
  # Always apply correction to ensure constraints are met
  # This is necessary for numerical stability with small tau values
  if (row_error > 1e-6 || col_error > 1e-6) {
    # Apply alternating row/column normalization (Sinkhorn-style balancing)
    for (iter in 1:20) {
      # Normalize rows to match mu
      row_sums_curr <- rowSums(S)
      row_sums_curr[row_sums_curr < 1e-16] <- 1e-16  # Avoid division by zero
      S <- S / row_sums_curr * mu
      
      # Normalize columns to match nu
      col_sums_curr <- colSums(S)
      col_sums_curr[col_sums_curr < 1e-16] <- 1e-16
      S <- t(t(S) / col_sums_curr * nu)
      
      # Check convergence
      row_error_new <- max(abs(rowSums(S) - mu))
      col_error_new <- max(abs(colSums(S) - nu))
      
      if (row_error_new < 1e-9 && col_error_new < 1e-9) {
        break
      }
    }
  }
  
  return(S)
}

#' Log-domain Stabilized Sinkhorn Algorithm (Dispatcher)
#' 
#' Implements numerically stable Sinkhorn in log domain.
#' Routes to either R or C++ implementation based on use_cpp flag.
#' 
#' @param C Cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param use_cpp Whether to use C++ implementation
#' @return Transport plan matrix
#' @keywords internal
solve_sinkhorn_stabilized <- function(C, tau, max_iter, tol, use_cpp = FALSE) {
  if (inherits(C, "Matrix")) {
    C <- as.matrix(C)
  }

  if (isTRUE(use_cpp)) {
    solve_sinkhorn_stabilized_cpp_fallback(C, tau, max_iter, tol)
  } else {
    solve_sinkhorn_stabilized_r(C, tau, max_iter, tol)
  }
}
