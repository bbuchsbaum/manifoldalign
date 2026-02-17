call_cpp_or <- function(call_cpp, fallback) {
  res <- tryCatch(
    call_cpp(),
    error = function(e) structure(list(error = e), class = "cpp_error")
  )
  if (!inherits(res, "cpp_error")) {
    return(res)
  }
  fallback()
}

compute_edge_gradient_cpp_fallback <- function(S, A1, A2, X1, X2) {
  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_compute_edge_gradient_cpp`, S, A1, A2, X1, X2)
    },
    fallback = function() {
      compute_edge_gradient_r(S, list(
        list(adjacency = A1, features = X1),
        list(adjacency = A2, features = X2)
      ))
    }
  )
}

compute_squared_distances_cpp_fallback <- function(X1, X2) {
  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_compute_squared_distances_cpp`, X1, X2)
    },
    fallback = function() {
      X1_sq <- rowSums(X1^2)
      X2_sq <- rowSums(X2^2)
      outer(X1_sq, X2_sq, "+") - 2 * X1 %*% t(X2)
    }
  )
}

solve_sinkhorn_stabilized_cpp_fallback <- function(C_in, tau, max_iter, tol) {
  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_solve_sinkhorn_stabilized_cpp`, C_in, tau, max_iter, tol)
    },
    fallback = function() {
      solve_sinkhorn_stabilized_r(C_in, tau, max_iter, tol)
    }
  )
}

compute_rwr_vectorized_r_single <- function(W_transpose, E, sigma, max_iter, tol) {
  R <- as.matrix(E) / nrow(W_transpose)
  for (iter in seq_len(max_iter)) {
    R_prev <- R
    R <- (1 - sigma) * (W_transpose %*% R) + sigma * E
    if (max(abs(R - R_prev)) < tol) {
      break
    }
  }
  as.matrix(R)
}

compute_rwr_vectorized_cpp_fallback <- function(W_transpose, E, sigma, max_iter, tol) {
  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_compute_rwr_vectorized_cpp`, W_transpose, E, sigma, max_iter, tol)
    },
    fallback = function() {
      compute_rwr_vectorized_r_single(W_transpose, E, sigma, max_iter, tol)
    }
  )
}

solve_sylvester_rwr_cpp_fallback <- function(W1, W2T, C_node, beta, gamma, tol, max_iter) {
  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_solve_sylvester_rwr_cpp`, W1, W2T, C_node, beta, gamma, tol, max_iter)
    },
    fallback = function() {
      solve_sylvester_rwr_r(W1, W2T, C_node, beta, gamma, tol, max_iter)
    }
  )
}

compute_parrot_cost_cpp_fallback <- function(networks, rwr_features, anchor_info,
                                             alpha, sigma, gamma) {
  X1 <- networks[[1]]$features
  X2 <- networks[[2]]$features
  R1 <- rwr_features[[1]]
  R2 <- rwr_features[[2]]
  W1 <- networks[[1]]$transition
  W2 <- networks[[2]]$transition

  to_dense <- function(mat) {
    if (inherits(mat, "Matrix")) {
      return(as.matrix(mat))
    }
    mat
  }

  X1 <- to_dense(X1)
  X2 <- to_dense(X2)
  R1 <- to_dense(R1)
  R2 <- to_dense(R2)
  W1 <- to_dense(W1)
  W2 <- to_dense(W2)

  call_cpp_or(
    call_cpp = function() {
      .Call(`_manifoldalign_compute_parrot_cost_cpp`, X1, X2, R1, R2, W1, W2, alpha, sigma, gamma)
    },
    fallback = function() {
      compute_parrot_cost_r(networks, rwr_features, anchor_info, alpha, sigma, gamma)
    }
  )
}
