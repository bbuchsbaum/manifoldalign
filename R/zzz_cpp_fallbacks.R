# Compiled-kernel dispatch with R fallbacks
#
# `call_cpp_or()` runs the compiled routine whenever its native symbol is
# registered and loaded. Errors raised by the compiled code (bad input,
# numerical failure, ...) are re-signalled unchanged: they are never turned
# into a silent rerun of the R code. The R fallback runs only when the symbol
# is genuinely unavailable (e.g. the shared library failed to load), or when
# the internal option `manifoldalign.force_r_fallback = TRUE` is set (used by
# the equivalence tests). Every fallback below is required to compute the same
# quantity as its compiled counterpart. The backend that served the most
# recent call of each kernel is recorded and can be read back with
# `.last_cpp_backend()`.

.cpp_backend_state <- new.env(parent = emptyenv())

# The dense kernels below take base matrices; Matrix inputs are converted at
# this boundary so that both backends see identical values.
.dense_if_matrix <- function(x) {
  if (inherits(x, "Matrix")) as.matrix(x) else x
}

.cpp_symbol_available <- function(symbol) {
  ns <- asNamespace("manifoldalign")
  if (!exists(symbol, envir = ns, inherits = FALSE)) {
    return(FALSE)
  }
  isTRUE(tryCatch(
    is.loaded(symbol, PACKAGE = "manifoldalign", type = "Call"),
    error = function(e) FALSE
  ))
}

.last_cpp_backend <- function(kernel) {
  if (!exists(kernel, envir = .cpp_backend_state, inherits = FALSE)) {
    return(NA_character_)
  }
  get(kernel, envir = .cpp_backend_state, inherits = FALSE)
}

call_cpp_or <- function(symbol, call_cpp, fallback, kernel = symbol) {
  use_cpp <- !isTRUE(getOption("manifoldalign.force_r_fallback", FALSE)) &&
    .cpp_symbol_available(symbol)
  if (use_cpp) {
    assign(kernel, "cpp", envir = .cpp_backend_state)
    return(call_cpp())
  }
  assign(kernel, "r", envir = .cpp_backend_state)
  fallback()
}

compute_edge_gradient_cpp_fallback <- function(S, A1, A2, X1, X2) {
  call_cpp_or(
    "_manifoldalign_compute_edge_gradient_cpp",
    call_cpp = function() {
      .Call(`_manifoldalign_compute_edge_gradient_cpp`, S, A1, A2, X1, X2)
    },
    fallback = function() {
      compute_edge_gradient_r(S, list(
        list(adjacency = A1, features = X1),
        list(adjacency = A2, features = X2)
      ))
    },
    kernel = "compute_edge_gradient"
  )
}

compute_squared_distances_cpp_fallback <- function(X1, X2) {
  # The distance matrix is dense either way; accept Matrix inputs explicitly.
  X1 <- .dense_if_matrix(X1)
  X2 <- .dense_if_matrix(X2)
  call_cpp_or(
    "_manifoldalign_compute_squared_distances_cpp",
    call_cpp = function() {
      .Call(`_manifoldalign_compute_squared_distances_cpp`, X1, X2)
    },
    fallback = function() {
      compute_squared_distances_r(X1, X2)
    },
    kernel = "compute_squared_distances"
  )
}

# Uniform-marginal log-domain Sinkhorn on the raw cost. Both backends return
# list(plan, log_u, log_v, converged, iterations, marginal_error).
solve_sinkhorn_stabilized_cpp_fallback <- function(C_in, tau, max_iter, tol,
                                                   log_u0 = NULL, log_v0 = NULL) {
  n <- nrow(C_in)
  m <- ncol(C_in)
  call_cpp_or(
    "_manifoldalign_sinkhorn_unified_potentials",
    call_cpp = function() {
      res <- .Call(`_manifoldalign_sinkhorn_unified_potentials`, C_in,
                   rep(1 / n, n), rep(1 / m, m), tau, log_u0, log_v0,
                   as.integer(max_iter), tol, TRUE)
      list(
        plan = res$pi,
        log_u = as.numeric(res$log_u),
        log_v = as.numeric(res$log_v),
        converged = isTRUE(res$converged),
        iterations = as.integer(res$iterations),
        marginal_error = as.numeric(res$marginal_error)
      )
    },
    fallback = function() {
      solve_sinkhorn_stabilized_r(C_in, tau, max_iter, tol, log_u0, log_v0)
    },
    kernel = "solve_sinkhorn_stabilized"
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
  W_transpose <- .dense_if_matrix(W_transpose)
  E <- .dense_if_matrix(E)
  call_cpp_or(
    "_manifoldalign_compute_rwr_vectorized_cpp",
    call_cpp = function() {
      .Call(`_manifoldalign_compute_rwr_vectorized_cpp`, W_transpose, E, sigma, max_iter, tol)
    },
    fallback = function() {
      compute_rwr_vectorized_r_single(W_transpose, E, sigma, max_iter, tol)
    },
    kernel = "compute_rwr_vectorized"
  )
}

solve_sylvester_rwr_cpp_fallback <- function(W1, W2T, C_node, beta, gamma, tol, max_iter) {
  W1 <- .dense_if_matrix(W1)
  W2T <- .dense_if_matrix(W2T)
  C_node <- .dense_if_matrix(C_node)
  call_cpp_or(
    "_manifoldalign_solve_sylvester_rwr_cpp",
    call_cpp = function() {
      .Call(`_manifoldalign_solve_sylvester_rwr_cpp`, W1, W2T, C_node, beta, gamma, tol, max_iter)
    },
    fallback = function() {
      solve_sylvester_rwr_r(W1, W2T, C_node, beta, gamma, tol, max_iter)
    },
    kernel = "solve_sylvester_rwr"
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
    "_manifoldalign_compute_parrot_cost_cpp",
    call_cpp = function() {
      .Call(`_manifoldalign_compute_parrot_cost_cpp`, X1, X2, R1, R2, W1, W2, alpha, sigma, gamma)
    },
    fallback = function() {
      compute_parrot_cost_r(networks, rwr_features, anchor_info, alpha, sigma, gamma)
    },
    kernel = "compute_parrot_cost"
  )
}
