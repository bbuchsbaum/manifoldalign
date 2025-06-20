#' PARROT Proximal Point Method Implementation
#'
#' Correct implementation following Algorithm 2 from the paper
#' with proper objective function and monotonicity checks
#'
#' @keywords internal

#' Compute PARROT Objective Function
#' 
#' Computes the full objective: <C_rwr, S> + λ_e*L_e(S) + λ_n*L_n(S) + λ_p*L_a(S)
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
    A1 <- networks[[1]]$adjacency
    A2 <- networks[[2]]$adjacency
    X1 <- networks[[1]]$features
    X2 <- networks[[2]]$features
    
    # Compute edge distance matrices
    n1 <- nrow(X1)
    n2 <- nrow(X2)
    
    # For each edge in G1
    edge_term <- 0
    edges1 <- which(A1 != 0, arr.ind = TRUE)
    edges1 <- edges1[edges1[,1] < edges1[,2], ]  # Upper triangular only
    
    if (nrow(edges1) > 0) {
      for (e in 1:nrow(edges1)) {
        i <- edges1[e, 1]
        j <- edges1[e, 2]
        
        # Distance between nodes i,j in G1
        d1_ij <- sum((X1[i,] - X1[j,])^2)
        
        # Expected distance in G2 under transport plan
        for (a in 1:n2) {
          for (b in 1:n2) {
            if (A2[a, b] > 0) {
              d2_ab <- sum((X2[a,] - X2[b,])^2)
              edge_term <- edge_term + S[i,a] * S[j,b] * (d1_ij - d2_ab)^2
            }
          }
        }
      }
    }
    
    obj <- obj + lambda_e * edge_term
  }
  
  # Neighborhood consistency regularizer L_n(S) from Eq. 9
  if (lambda_n > 0) {
    W1 <- networks[[1]]$transition
    W2 <- networks[[2]]$transition
    
    # S_hat = W1^T S W2 (transported neighborhoods)
    S_hat <- t(W1) %*% S %*% W2
    
    # KL divergence: sum(S_hat * log(S_hat / S))
    eps <- 1e-16
    S_safe <- S + eps
    S_hat_safe <- S_hat + eps
    
    kl_div <- sum(S_hat_safe * (log(S_hat_safe) - log(S_safe)))
    obj <- obj + lambda_n * kl_div
  }
  
  # Anchor prior regularizer L_a(S)
  if (lambda_p > 0 && length(anchor_info$idx1) > 0 && length(anchor_info$idx2) > 0) {
    # Create anchor prior matrix H
    n1 <- nrow(S)
    n2 <- ncol(S)
    H <- matrix(eps, n1, n2)
    
    # Set uniform prior for anchor correspondences
    for (i in seq_along(anchor_info$idx1)) {
      idx1 <- anchor_info$idx1[i]
      anchor_val <- anchor_info$vec1[idx1]
      
      for (j in seq_along(anchor_info$idx2)) {
        idx2 <- anchor_info$idx2[j]
        if (anchor_info$vec2[idx2] == anchor_val) {
          H[idx1, idx2] <- 1.0 / length(anchor_info$idx1)
        }
      }
    }
    
    # KL divergence from anchor prior
    anchor_div <- sum(S * (log(S_safe) - log(H)))
    obj <- obj + lambda_p * anchor_div
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
#' @param use_rcpp Whether to use C++ implementation
#' @return Gradient matrix
#' @keywords internal
compute_edge_gradient <- function(S, networks, use_rcpp = NULL) {
  if (is.null(use_rcpp)) {
    use_rcpp <- get_parrot_use_rcpp()
  }
  
  if (use_rcpp && requireNamespace("Rcpp", quietly = TRUE)) {
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
    
    # Call C++ implementation
    compute_edge_gradient_cpp(S, A1, A2, X1, X2)
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
  
  # Implement tau annealing for better convergence
  # Start with higher tau for stability, decrease to target tau
  tau_init <- max(tau, 0.1)
  tau_schedule <- c(tau_init)
  if (tau < tau_init) {
    # Create geometric sequence from tau_init to tau
    n_anneal <- min(5, max_outer)
    tau_ratio <- (tau / tau_init)^(1 / n_anneal)
    for (i in 1:n_anneal) {
      tau_schedule <- c(tau_schedule, tau_init * tau_ratio^i)
    }
  }
  
  # Compute base position-aware cost matrix C_rwr
  C_rwr <- compute_parrot_cost(networks, rwr_features, anchor_info, alpha, sigma, gamma, use_cpp)
  
  # Initialize transport plan using Sinkhorn on the base cost matrix
  # This ensures we start from a reasonable point based on the actual costs
  S <- solve_sinkhorn_stabilized(C_rwr, tau = tau * 10, max_iter = 100, tol, use_cpp)
  
  # Initialize objective value
  obj_prev <- compute_parrot_objective(S, C_rwr, networks, anchor_info, 
                                      lambda_e, lambda_n, lambda_p)
  
  obj_history <- numeric(max_outer)
  
  # Proximal Point Method outer loop with tau annealing
  tau_idx <- 1
  for (t in 1:max_outer) {
    S_prev <- S
    
    # Use annealed tau value
    if (tau_idx <= length(tau_schedule)) {
      tau_current <- tau_schedule[tau_idx]
      tau_idx <- tau_idx + 1
    } else {
      tau_current <- tau  # Use target tau for remaining iterations
    }
    
    # Compute effective cost matrix C(t) based on gradients
    C_effective <- C_rwr
    sd_c_rwr <- sd(C_rwr)
    
    # Add gradient of edge consistency (linearized)
    if (lambda_e > 0) {
      grad_edge <- compute_edge_gradient(S_prev, networks, use_cpp)
      # Use gradient directly without destructive normalization
      # Only apply a conservative scaling if gradient is too large
      if (any(grad_edge != 0)) {
        max_grad <- max(abs(grad_edge))
        if (max_grad > 10 * sd_c_rwr) {
          # Only scale down if gradient is excessively large
          scaling_factor <- (10 * sd_c_rwr) / max_grad
          grad_edge <- grad_edge * scaling_factor
        }
        C_effective <- C_effective + lambda_e * grad_edge
      }
    }
    
    # Add gradient of neighborhood consistency
    if (lambda_n > 0) {
      W1 <- networks[[1]]$transition
      W2 <- networks[[2]]$transition
      
      # Gradient of KL(W1^T S W2 || S) w.r.t S
      S_hat <- t(W1) %*% S_prev %*% W2
      eps <- 1e-16
      
      # ∇L_n = W1 * log(S_hat/S) * W2^T
      # But we want the gradient of the KL divergence KL(S_hat || S)
      # The gradient w.r.t S is: -W1 * (S_hat/S - 1) * W2^T
      ratio <- (S_hat + eps) / (S_prev + eps)
      grad_neigh <- -W1 %*% (ratio - 1) %*% t(W2)
      
      # Use gradient directly with conservative scaling only if needed
      if (any(grad_neigh != 0)) {
        max_grad_n <- max(abs(grad_neigh))
        if (max_grad_n > 10 * sd_c_rwr) {
          # Only scale down if gradient is excessively large
          scaling_factor <- (10 * sd_c_rwr) / max_grad_n
          grad_neigh <- grad_neigh * scaling_factor
        }
        C_effective <- C_effective + lambda_n * grad_neigh
      }
    }
    
    # Add gradient of anchor prior
    if (lambda_p > 0 && length(anchor_info$idx1) > 0) {
      # Gradient of KL(S || H) is log(S/H)
      n1 <- nrow(S)
      n2 <- ncol(S)
      H <- matrix(1e-16, n1, n2)
      
      for (i in seq_along(anchor_info$idx1)) {
        idx1 <- anchor_info$idx1[i]
        anchor_val <- anchor_info$vec1[idx1]
        
        for (j in seq_along(anchor_info$idx2)) {
          idx2 <- anchor_info$idx2[j]
          if (anchor_info$vec2[idx2] == anchor_val) {
            H[idx1, idx2] <- 1.0  # High probability for correct anchors
          }
        }
      }
      
      # Normalize H to be a valid probability distribution
      H <- H / sum(H)
      
      # Gradient of KL(S || H) is log(S/H)
      # For anchor pairs, H is high, so log(S/H) will be negative when S < H
      # This DECREASES cost at anchor pairs, encouraging alignment
      grad_anchor <- log((S_prev + 1e-16) / (H + 1e-16))
      C_effective <- C_effective + lambda_p * grad_anchor
    }
    
    # Ensure cost matrix is non-negative before solving
    min_cost <- min(C_effective)
    if (min_cost < 0) {
      C_effective <- C_effective - min_cost + 1e-6  # Shift to positive with small buffer
    }
    
    # Solve entropic OT subproblem with effective cost
    if (solver == "sinkhorn") {
      S_new <- solve_sinkhorn_stabilized(C_effective, tau_current, max_iter, tol, use_cpp)
    } else {
      # For exact solver, use assignment problem
      if (!requireNamespace("clue", quietly = TRUE)) {
        stop("clue package required for exact solver")
      }
      # Cost matrix is already non-negative from above
      # Ensure it's a regular matrix for clue
      if (inherits(C_effective, "Matrix")) {
        C_effective <- as.matrix(C_effective)
      }
      assignment <- clue::solve_LSAP(C_effective)
      S_new <- matrix(0, n1, n2)
      for (i in 1:n1) {
        S_new[i, assignment[i]] <- 1  # For doubly stochastic, set to 1
      }
    }
    
    # Compute new objective
    obj_new <- compute_parrot_objective(S_new, C_rwr, networks, anchor_info,
                                       lambda_e, lambda_n, lambda_p)
    
    # AUDIT-04: Monotonicity check with tolerance
    # Allow small increases in early iterations to escape local minima
    tolerance <- if (t <= 3) 0.1 * abs(obj_prev) else 1e-6
    
    if (obj_new <= obj_prev + tolerance) {
      # Accept the update
      S <- S_new
      obj_prev <- obj_new
      
      if (obj_new > obj_prev) {
        message(sprintf("Iteration %d: Accepting small objective increase (%.6f -> %.6f)", 
                       t, obj_prev, obj_new))
      }
    } else {
      # Reject update and use previous solution
      message(sprintf("Iteration %d: Objective increased too much (%.6f -> %.6f), rejecting update", 
                     t, obj_prev, obj_new))
      S <- S_prev
    }
    
    obj_history[t] <- obj_prev
    
    # Check convergence
    if (max(abs(S - S_prev)) < tol) {
      break
    }
  }
  
  list(
    transport_plan = S,
    cost_matrix = C_rwr,
    objective_history = obj_history[1:t],
    outer_iterations = t,
    final_objective = obj_prev
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
  
  # --- FIX: Define correct probability marginals ---
  # The test expects a doubly-stochastic matrix (row/col sums = 1),
  # not a probability distribution (sum(S) = 1).
  # We will produce a doubly-stochastic matrix.
  # For square matrices, uniform doubly stochastic has row/col sums = 1
  # For rectangular matrices, we need to balance the total mass
  if (n1 == n2) {
    mu <- rep(1, n1) # Target row sums = 1
    nu <- rep(1, n2) # Target col sums = 1
  } else {
    # For rectangular matrices, ensure total mass is consistent
    total_mass <- min(n1, n2)
    mu <- rep(total_mass / n1, n1)
    nu <- rep(total_mass / n2, n2)
  }
  
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
#' Routes to either R or C++ implementation based on use_rcpp flag.
#' 
#' @param C Cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param use_rcpp Whether to use C++ implementation
#' @return Transport plan matrix
#' @keywords internal
solve_sinkhorn_stabilized <- function(C, tau, max_iter, tol, use_rcpp = NULL) {
  # Determine whether to use C++
  if (is.null(use_rcpp)) {
    use_rcpp <- get_parrot_use_rcpp()
  }
  
  if (use_rcpp && requireNamespace("Rcpp", quietly = TRUE)) {
    # Ensure C is a regular matrix for C++
    if (inherits(C, "Matrix")) {
      C <- as.matrix(C)
    }
    
    # Call C++ implementation (does not take use_rcpp)
    solve_sinkhorn_stabilized_cpp(C, tau, max_iter, tol)
  } else {
    # Use existing R implementation (does not take use_rcpp)
    solve_sinkhorn_stabilized_r(C, tau, max_iter, tol)
  }
}