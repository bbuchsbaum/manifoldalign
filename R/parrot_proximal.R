# PARROT objective, gradient and proximal solver
#
# Estimand: a transport plan S in Pi(mu, nu) (uniform marginals) minimising
#
#   F(S) = <C_rwr, S> + lambda_e L_e(S)
#          + lambda_n KL(S || Shat_delta(S)) + lambda_p KL(S || H)
#
# where C_rwr is the position-aware cost, L_e the edge-consistency term
#   L_e(S) = sum_ij (a_i + b_j) S_ij - 2 <C1 S C2^T, S>
# (C1, C2 are cosine feature distances restricted to graph edges,
# a = rowSums(C1^2), b = rowSums(C2^2)), Shat(S) = W1^T S W2 is the plan
# propagated through the row-stochastic transition matrices, and H is the
# normalised anchor prior (1 on anchor pairs, 1e-6 elsewhere, before
# normalisation). Shat is mixed with the product plan,
# Shat_delta = (1 - delta) Shat + delta / (n1 n2), so the KL term is finite
# for every feasible S.
#
# Solver: Bregman (KL) proximal gradient. Writing F = f + h with
# h(S) = (lambda_n + lambda_p) sum S log S kept exact and f smooth, each outer
# step solves
#   S+ = argmin_{S in Pi} <grad f(S_t), S> + h(S) + tau_t KL(S || S_t),
# an entropic OT problem with cost grad f(S_t) - tau_t log S_t and
# regularisation tau_t + lambda_n + lambda_p. A step is accepted only if it
# does not increase F; otherwise tau_t is doubled (a stronger proximal pull
# towards S_t) and the step retried.

.parrot_kl_delta <- 1e-8

# Coerce dense or sparse (possibly symmetric/pattern) input to dgCMatrix.
.parrot_as_dgc <- function(A) {
  if (inherits(A, "dgCMatrix")) return(A)
  A <- methods::as(A, "dMatrix")
  A <- methods::as(A, "generalMatrix")
  methods::as(A, "CsparseMatrix")
}

.parrot_row_normalize <- function(M) {
  norms <- sqrt(rowSums(M^2))
  norms[norms < 1e-12] <- 1e-12
  M / norms
}

# Cosine feature distance on graph edges, as a sparse matrix with the
# adjacency pattern: C[i, k] = 1 - cos(x_i, x_k) for A[i, k] != 0.
.parrot_edge_cost <- function(A, X) {
  X <- as.matrix(X)
  Xn <- .parrot_row_normalize(X)
  A <- .parrot_as_dgc(A)
  trip <- Matrix::summary(A)
  keep <- trip$x != 0
  i <- trip$i[keep]
  k <- trip$j[keep]
  x <- 1 - rowSums(Xn[i, , drop = FALSE] * Xn[k, , drop = FALSE])
  Matrix::sparseMatrix(i = i, j = k, x = x, dims = dim(A))
}

# Anchor prior H: weight 1 on declared anchor pairs and `eps` elsewhere,
# normalised to sum 1 (uniform when there are no anchors). The floor 1e-6
# gives anchor pairs prior odds of 1e6 while keeping -log H bounded; a much
# smaller floor (the former 1e-12) makes the inner Sinkhorn problems nearly
# decomposable and their convergence sublinear without changing the matches.
.parrot_anchor_prior <- function(anchor_info, n1, n2, eps = 1e-6) {
  H <- matrix(eps, n1, n2)
  idx1 <- anchor_info$idx1
  idx2 <- anchor_info$idx2
  if (length(idx1) > 0 && length(idx2) > 0) {
    vals1 <- anchor_info$vec1[idx1]
    vals2 <- anchor_info$vec2[idx2]
    match_mat <- outer(as.character(vals1), as.character(vals2), "==")
    hit <- which(match_mat, arr.ind = TRUE)
    if (nrow(hit) > 0) {
      H[cbind(idx1[hit[, 1]], idx2[hit[, 2]])] <- 1
    }
  } else {
    H[] <- 1
  }
  H / sum(H)
}

.parrot_xlogx <- function(S) {
  S * log(pmax(S, .Machine$double.xmin))
}

# Precompute everything the objective and gradient need.
.parrot_problem <- function(networks, C_rwr, anchor_info,
                            lambda_e, lambda_n, lambda_p,
                            kl_delta = .parrot_kl_delta) {
  C_rwr <- as.matrix(C_rwr)
  n1 <- nrow(C_rwr)
  n2 <- ncol(C_rwr)

  prob <- list(
    C_rwr = C_rwr, n1 = n1, n2 = n2,
    lambda_e = lambda_e, lambda_n = lambda_n, lambda_p = lambda_p,
    kl_delta = kl_delta
  )

  if (lambda_e > 0) {
    C1 <- .parrot_edge_cost(networks[[1]]$adjacency, networks[[1]]$features)
    C2 <- .parrot_edge_cost(networks[[2]]$adjacency, networks[[2]]$features)
    prob$C1 <- C1
    prob$C2 <- C2
    prob$C1t <- Matrix::t(C1)
    prob$C2t <- Matrix::t(C2)
    prob$a_edge <- Matrix::rowSums(C1^2)
    prob$b_edge <- Matrix::rowSums(C2^2)
  }
  if (lambda_n > 0) {
    prob$W1 <- networks[[1]]$transition
    prob$W2 <- networks[[2]]$transition
    prob$W1t <- Matrix::t(prob$W1)
    prob$W2t <- Matrix::t(prob$W2)
  }
  if (lambda_p > 0) {
    prob$log_H <- log(.parrot_anchor_prior(anchor_info, n1, n2))
  }
  prob
}

.parrot_edge_quadratic <- function(prob, S) {
  as.matrix(prob$C1 %*% S %*% prob$C2t)
}

.parrot_shat <- function(prob, S) {
  Shat <- as.matrix(prob$W1t %*% S %*% prob$W2)
  (1 - prob$kl_delta) * Shat + prob$kl_delta / (prob$n1 * prob$n2)
}

# Full objective F(S).
.parrot_objective <- function(prob, S) {
  obj <- sum(prob$C_rwr * S)
  if (prob$lambda_e > 0) {
    lin <- sum(prob$a_edge * rowSums(S)) + sum(prob$b_edge * colSums(S))
    obj <- obj + prob$lambda_e * (lin - 2 * sum(S * .parrot_edge_quadratic(prob, S)))
  }
  if (prob$lambda_n > 0) {
    Shat <- .parrot_shat(prob, S)
    obj <- obj + prob$lambda_n * (sum(.parrot_xlogx(S)) - sum(S * log(Shat)))
  }
  if (prob$lambda_p > 0) {
    obj <- obj + prob$lambda_p * (sum(.parrot_xlogx(S)) - sum(S * prob$log_H))
  }
  obj
}

# Gradient of the smooth part f = F - (lambda_n + lambda_p) sum S log S.
.parrot_smooth_gradient <- function(prob, S) {
  G <- prob$C_rwr
  if (prob$lambda_e > 0) {
    quad <- as.matrix(prob$C1 %*% S %*% prob$C2t) +
      as.matrix(prob$C1t %*% S %*% prob$C2)
    G <- G + prob$lambda_e * (outer(prob$a_edge, prob$b_edge, "+") - 2 * quad)
  }
  if (prob$lambda_n > 0) {
    Shat <- .parrot_shat(prob, S)
    back <- as.matrix(prob$W1 %*% (S / Shat) %*% prob$W2t)
    G <- G - prob$lambda_n * (log(Shat) + (1 - prob$kl_delta) * back)
  }
  if (prob$lambda_p > 0) {
    G <- G - prob$lambda_p * prob$log_H
  }
  G
}

#' Compute PARROT objective
#'
#' Evaluates the PARROT objective
#' \deqn{F(S) = \langle C_{rwr}, S\rangle + \lambda_e L_e(S) +
#'   \lambda_n KL(S \| \hat S_\delta(S)) + \lambda_p KL(S \| H),}
#' with \eqn{L_e(S) = \sum_{ij} (a_i + b_j) S_{ij} - 2 \langle C_1 S C_2^T, S\rangle}
#' (cosine feature distances on graph edges), \eqn{\hat S(S) = W_1^T S W_2}
#' mixed with the product plan, \eqn{\hat S_\delta = (1-\delta)\hat S +
#' \delta/(n_1 n_2)} with \eqn{\delta = 10^{-8}}, and \eqn{H} the anchor
#' prior (1 on declared anchor pairs, \eqn{10^{-6}} elsewhere, normalised to
#' sum 1; uniform when there are no anchors).
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
  prob <- .parrot_problem(networks, C_rwr, anchor_info,
                          lambda_e, lambda_n, lambda_p)
  .parrot_objective(prob, as.matrix(S))
}

#' Compute Gradient of Edge Consistency Regularizer (R implementation)
#'
#' Sparse-product evaluation of the edge-length consistency term used by
#' \code{compute_edge_gradient_cpp()}: for undirected edges \eqn{i<j} of
#' network 1 and \eqn{a<b} of network 2 with squared feature distances
#' \eqn{d^1_{ij}}, \eqn{d^2_{ab}}, each edge pair adds
#' \eqn{2 (d^1_{ij} - d^2_{ab}) S_{jb}} to \eqn{G_{ia}} and
#' \eqn{2 (d^1_{ij} - d^2_{ab}) S_{ia}} to \eqn{G_{jb}}.
#'
#' @param S Current transport plan
#' @param networks Network structures
#' @return Gradient matrix
#' @keywords internal
compute_edge_gradient_r <- function(S, networks) {
  S <- as.matrix(S)
  edge_mats <- function(A, X) {
    X <- as.matrix(X)
    A <- .parrot_as_dgc(A)
    trip <- Matrix::summary(A)
    keep <- trip$x != 0 & trip$i < trip$j
    i <- trip$i[keep]
    j <- trip$j[keep]
    d <- rowSums((X[i, , drop = FALSE] - X[j, , drop = FALSE])^2)
    list(
      D = Matrix::sparseMatrix(i = i, j = j, x = d, dims = dim(A)),
      U = Matrix::sparseMatrix(i = i, j = j, x = rep(1, length(i)), dims = dim(A))
    )
  }
  e1 <- edge_mats(networks[[1]]$adjacency, networks[[1]]$features)
  e2 <- edge_mats(networks[[2]]$adjacency, networks[[2]]$features)

  G <- e1$D %*% S %*% Matrix::t(e2$U) -
    e1$U %*% S %*% Matrix::t(e2$D) +
    Matrix::t(e1$D) %*% S %*% e2$U -
    Matrix::t(e1$U) %*% S %*% e2$D
  2 * as.matrix(G)
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
    A1 <- .parrot_as_dgc(A1)
    A2 <- .parrot_as_dgc(A2)
    if (inherits(X1, "Matrix")) X1 <- as.matrix(X1)
    if (inherits(X2, "Matrix")) X2 <- as.matrix(X2)
    if (inherits(S, "Matrix")) S <- as.matrix(S)

    compute_edge_gradient_cpp_fallback(S, A1, A2, X1, X2)
  } else {
    compute_edge_gradient_r(S, networks)
  }
}

#' Solve PARROT with a Bregman Proximal Gradient Method
#'
#' Minimises the PARROT objective (see \code{compute_parrot_objective()}) over
#' transport plans with uniform marginals. Each outer step linearises the
#' smooth part of the objective at the current plan \eqn{S_t}, keeps the
#' entropic parts of the two KL terms exact, and adds the proximal term
#' \eqn{\tau_t KL(S \| S_t)}; the step is an entropic OT problem with cost
#' \eqn{\nabla f(S_t) - \tau_t \log S_t} and regularisation
#' \eqn{\tau_t + \lambda_n + \lambda_p}, solved by log-domain Sinkhorn. A
#' step is accepted only if the objective does not increase; otherwise
#' \eqn{\tau_t} is doubled (up to 30 times) and the step retried. The loop
#' starts from the product plan and stops when the largest entrywise plan
#' change falls below \code{tol}, when \code{max_outer} steps were taken, or
#' when no descent step can be found.
#'
#' @inheritParams solve_parrot_transport
#' @return List with \code{transport_plan}, \code{cost_matrix} (the
#'   position-aware cost), \code{objective_history} (objective at the initial
#'   product plan followed by the objective after each accepted step),
#'   \code{outer_iterations}, \code{final_objective}, \code{converged},
#'   \code{status} (\code{"converged"}, \code{"max_outer"} or
#'   \code{"no_descent"}), \code{backtracks}, and \code{sinkhorn} (convergence
#'   diagnostics of the inner solve that produced the returned plan).
#' @keywords internal
solve_parrot_proximal <- function(networks, rwr_features, anchor_info,
                                  lambda_e, lambda_n, lambda_p,
                                  tau, alpha, sigma, gamma, solver,
                                  max_iter, max_outer = 10, tol, use_cpp = FALSE) {

  if (!identical(solver, "sinkhorn")) {
    stop("PARROT supports solver = \"sinkhorn\" only.", call. = FALSE)
  }
  max_outer <- as.integer(max_outer)
  if (length(max_outer) != 1L || is.na(max_outer) || max_outer < 1L) {
    stop("max_outer must be a positive integer.", call. = FALSE)
  }

  C_rwr <- compute_parrot_cost(networks, rwr_features, anchor_info, alpha, sigma, gamma, use_cpp)
  prob <- .parrot_problem(networks, C_rwr, anchor_info, lambda_e, lambda_n, lambda_p)

  n1 <- prob$n1
  n2 <- prob$n2
  ent <- lambda_n + lambda_p
  max_backtrack <- 30L

  S <- matrix(1 / (n1 * n2), n1, n2)
  F_cur <- .parrot_objective(prob, S)
  obj_history <- F_cur
  step <- tau
  backtracks <- 0L
  status <- "max_outer"
  converged <- FALSE
  outer_iterations <- 0L
  sink_info <- list(converged = NA, iterations = 0L, marginal_error = 0)
  warm_u <- NULL
  warm_v <- NULL

  for (t in seq_len(max_outer)) {
    G <- .parrot_smooth_gradient(prob, S)
    log_S <- log(pmax(S, .Machine$double.xmin))

    accepted <- FALSE
    for (bt in 0:max_backtrack) {
      sk <- solve_sinkhorn_stabilized(G - step * log_S, step + ent,
                                      max_iter, tol, use_cpp,
                                      log_u0 = warm_u, log_v0 = warm_v)
      F_new <- .parrot_objective(prob, sk$plan)
      if (is.finite(F_new) && F_new <= F_cur + 1e-12 * max(1, abs(F_cur))) {
        accepted <- TRUE
        break
      }
      step <- 2 * step
      backtracks <- backtracks + 1L
    }

    if (!accepted) {
      status <- "no_descent"
      break
    }

    delta <- max(abs(sk$plan - S))
    S <- sk$plan
    F_cur <- F_new
    obj_history <- c(obj_history, F_cur)
    outer_iterations <- t
    sink_info <- sk[c("converged", "iterations", "marginal_error")]
    # Warm-start the next inner solve from the accepted potentials.
    warm_u <- sk$log_u
    warm_v <- sk$log_v

    if (delta < tol) {
      converged <- TRUE
      status <- "converged"
      break
    }
    step <- max(tau, step / 2)
  }

  list(
    transport_plan = S,
    cost_matrix = C_rwr,
    objective_history = obj_history,
    outer_iterations = outer_iterations,
    final_objective = F_cur,
    converged = converged,
    status = status,
    backtracks = backtracks,
    sinkhorn = sink_info
  )
}


.lse_cols <- function(M) {
  mx <- M[cbind(max.col(t(M), ties.method = "first"), seq_len(ncol(M)))]
  mx[!is.finite(mx)] <- 0
  mx + log(colSums(exp(M - rep(mx, each = nrow(M)))))
}

.lse_rows <- function(M) {
  mx <- M[cbind(seq_len(nrow(M)), max.col(M, ties.method = "first"))]
  mx[!is.finite(mx)] <- 0
  mx + log(rowSums(exp(M - mx)))
}

#' Log-domain Sinkhorn with uniform marginals (R implementation)
#'
#' Solves \eqn{\min_{S \in \Pi(\mu,\nu)} \langle C, S\rangle + \tau \sum S(\log S - 1)}
#' with uniform marginals on the raw (unscaled) cost, using the same
#' iteration and stopping rule as the compiled \code{sinkhorn_unified()}:
#' each sweep updates the column then the row potentials, and the loop stops
#' when the marginal residual is below \code{tol} or after \code{max_iter}
#' sweeps. The plan is returned as computed (no rebalancing).
#'
#' @param C Cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance on the marginal residual
#' @param log_u0,log_v0 Optional warm-start log scalings (lengths
#'   \code{nrow(C)} and \code{ncol(C)})
#' @return List with \code{plan}, \code{log_u}, \code{log_v} (the plan is
#'   exactly \code{exp(-C/tau + log_u + t(log_v))}), \code{converged},
#'   \code{iterations} and \code{marginal_error}.
#' @keywords internal
solve_sinkhorn_stabilized_r <- function(C, tau, max_iter, tol,
                                        log_u0 = NULL, log_v0 = NULL) {
  C <- as.matrix(C)
  if (any(!is.finite(C))) {
    stop("Sinkhorn: cost matrix contains non-finite values.", call. = FALSE)
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau <= 0) {
    stop("Sinkhorn: epsilon must be a finite positive number.", call. = FALSE)
  }
  n1 <- nrow(C)
  n2 <- ncol(C)
  mu <- rep(1 / n1, n1)
  nu <- rep(1 / n2, n2)
  log_mu <- log(mu)
  log_nu <- log(nu)

  log_K <- -C / tau
  log_u <- if (is.null(log_u0)) numeric(n1) else as.numeric(log_u0)
  log_v <- if (is.null(log_v0)) numeric(n2) else as.numeric(log_v0)
  if (length(log_u) != n1 || length(log_v) != n2 ||
      any(!is.finite(log_u)) || any(!is.finite(log_v))) {
    stop("Sinkhorn: warm-start potentials must be finite with lengths nrow(C), ncol(C).",
         call. = FALSE)
  }
iterations <- 0L

  for (iter in seq_len(max_iter)) {
    lc <- .lse_cols(log_K + log_u)
    if (iter > 1L) {
      col_err <- max(abs(exp(log_v + lc) - nu))
      if (is.finite(col_err) && col_err < tol) break
    }
    log_v <- log_nu - lc
    lr <- .lse_rows(log_K + rep(log_v, each = n1))
    log_u <- log_mu - lr
    iterations <- iter
  }

  S <- exp(log_K + log_u + rep(log_v, each = n1))
  marginal_error <- max(abs(rowSums(S) - mu), abs(colSums(S) - nu))
  if (!is.finite(marginal_error)) marginal_error <- Inf

  list(
    plan = S,
    log_u = log_u,
    log_v = log_v,
    converged = marginal_error < tol,
    iterations = as.integer(iterations),
    marginal_error = marginal_error
  )
}

#' Log-domain Sinkhorn with uniform marginals (Dispatcher)
#'
#' Routes to the R or the compiled implementation; both solve the same problem
#' on the raw cost and return the same structure.
#'
#' @param C Cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance on the marginal residual
#' @param use_cpp Whether to use C++ implementation
#' @param log_u0,log_v0 Optional warm-start log scalings
#' @return List with \code{plan}, \code{log_u}, \code{log_v},
#'   \code{converged}, \code{iterations} and \code{marginal_error}.
#' @keywords internal
solve_sinkhorn_stabilized <- function(C, tau, max_iter, tol, use_cpp = FALSE,
                                      log_u0 = NULL, log_v0 = NULL) {
  if (inherits(C, "Matrix")) {
    C <- as.matrix(C)
  }

  if (isTRUE(use_cpp)) {
    solve_sinkhorn_stabilized_cpp_fallback(C, tau, max_iter, tol, log_u0, log_v0)
  } else {
    solve_sinkhorn_stabilized_r(C, tau, max_iter, tol, log_u0, log_v0)
  }
}
