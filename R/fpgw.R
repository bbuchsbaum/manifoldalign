#' Fused-Partial Gromov-Wasserstein Distance
#'
#' Computes the Fused-Partial Gromov-Wasserstein (FPGW) distance between domains
#' in a hyperdesign object. This method extends Gromov-Wasserstein by combining
#' feature and structural information. Partial transport is supported through
#' an explicit transported-mass budget, \code{rho}.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param ... Additional arguments passed to methods:
#'   \itemize{
#'     \item \code{omega1}: Weight of the feature term (0 ≤ omega1 ≤ 1).
#'       The structural term weight omega2 = 1 - omega1 is computed
#'       automatically (default: 0.001)
#'     \item \code{lambda}: Retained for API compatibility and required to be
#'       zero. The former \code{lambda > 0} route is withdrawn because it used
#'       a balanced transport oracle and therefore could not produce a partial
#'       coupling. Use \code{rho} instead (default: 0).
#'     \item \code{rho}: Mass budget in (0, min(|μ|,|ν|)] for the
#'       mass-constrained variant. If \code{rho} is supplied, the solver
#'       transports exactly that mass; otherwise it uses classical balanced
#'       FGW.
#'     \item \code{epsilon}: Entropic regularization of the warm start
#'       (default: 0.01). Set to 0 to disable the warm start
#'     \item \code{metric}: Distance metric, one of \code{"euclidean"},
#'       \code{"manhattan"} or \code{"cosine"} (default: "euclidean")
#'     \item \code{max_iter}: Maximum iterations for Frank-Wolfe optimization
#'       (default: 200)
#'     \item \code{tol}: Convergence tolerance for the Frank-Wolfe gap
#'       (default: 1e-6)
#'     \item \code{inner_max_iter}: Number of Sinkhorn iterations of the
#'       entropic warm start (default: 50)
#'     \item \code{verbose}: Print convergence information (default: FALSE)
#'   }
#'
#' @return An object of class \code{fpgw} containing:
#' \itemize{
#'   \item \code{transport_plans}: List of optimal transport plans between domains
#'     (upper-triangular pair order)
#'   \item \code{distances}: Matrix of pairwise FPGW distances between domains
#'     (square root of the minimised objective)
#'   \item \code{converged}: Convergence status for each optimization (Frank-Wolfe
#'     gap below \code{tol})
#'   \item \code{n_samples}: Number of samples per domain
#'   \item \code{omega1}: Feature weight used
#'   \item \code{lambda}: Compatibility field; always zero for new fits
#'   \item \code{rho}: Mass budget (if mass-constrained variant)
#'   \item \code{domain_names}: Names of the domains
#'   \item \code{training_data}: Original domain matrices for prediction
#'   \item \code{metric}: Distance metric used
#' }
#' The object is not a projector: it has no scores, loadings or
#' preprocessing pipeline, and it does not inherit from
#' \code{multiblock_biprojector}.
#'
#' @details
#' With uniform marginals \eqn{\mu, \nu} the verified fixed-mass FPGW
#' objective is
#' \deqn{L(\gamma) = \omega_1 \langle C, \gamma \rangle +
#'       \omega_2 \sum_{i,j,i',j'} |C_{X,ii'} - C_{Y,jj'}|^2 \gamma_{ij} \gamma_{i'j'}}
#'
#' where C is the feature cost matrix and C_X, C_Y are within-domain distance
#' matrices (scaled to a maximum of one).
#'
#' The implementation supports two verified modes:
#' \itemize{
#'   \item Classical FGW (\code{rho = NULL}): \eqn{\gamma} has marginals
#'     exactly \eqn{\mu} and \eqn{\nu}.
#'   \item Mass-constrained FPGW (\code{rho}): \eqn{\gamma \ge 0} with row sums
#'     at most \eqn{\mu}, column sums at most \eqn{\nu} and total mass exactly
#'     \code{rho}.
#' }
#'
#' The former penalized mode (\code{lambda > 0}) selected balanced Frank-Wolfe
#' vertices, fixing the transported mass at one and making \code{lambda}
#' ineffective as a partial transport control. It is rejected rather than
#' silently returning a balanced result.
#'
#' Each pair is solved by Frank-Wolfe started from a feasible plan (the
#' rounded entropic warm start, scaled to mass \code{rho} for the
#' mass-constrained variant). The linear minimisation oracle is an exact
#' network-simplex LP over the variant's feasible set and the step size is the
#' exact minimiser of the quadratic objective along the search direction. The
#' returned plans are checked against the marginal bounds and total mass of
#' their variant; a violation is an error.
#'
#' @references
#' Bai et al. (2025). Fused Partial Gromov-Wasserstein for Structured Objects.
#' arXiv:2502.09934.
#'
#' @examples
#' \donttest{
#' library(multidesign)
#'
#' # Example 1: Basic FPGW between two domains with different dimensions
#' set.seed(123)
#' n <- 30
#'
#' # Domain 1: 3D data with two clusters
#' X1 <- matrix(rnorm(n * 3), n, 3)
#' X1[1:15, ] <- X1[1:15, ] + 2
#'
#' # Domain 2: 5D data with similar structure
#' X2 <- matrix(rnorm(n * 5), n, 5)
#' X2[1:15, ] <- X2[1:15, ] + 2
#'
#' # Create hyperdesign
#' design <- data.frame(id = 1:n, cluster = rep(1:2, each = 15))
#' hd <- hyperdesign(list(
#'   visual = multidesign(X1, design),
#'   semantic = multidesign(X2, design)
#' ))
#'
#' # Classical FGW with balanced feature/structure weight
#' result <- fpgw(hd, omega1 = 0.5)
#' print(result)
#'
#' # Example 2: Mass-constrained for noisy data
#' # Add outliers to domain 2
#' X2_noisy <- rbind(X2, matrix(rnorm(10 * 5, sd = 5), 10, 5))
#' design_noisy <- data.frame(id = 1:(n + 10))
#'
#' hd_noisy <- hyperdesign(list(
#'   clean = multidesign(X1, design[1:n,]),
#'   noisy = multidesign(X2_noisy, design_noisy)
#' ))
#'
#' # Transport only 75% of mass to avoid outliers
#' result_partial <- fpgw(hd_noisy, omega1 = 0.3, rho = 0.75)
#'
#' # Check transported mass
#' P <- result_partial$transport_plans[[1]]
#' cat("Transported mass:", sum(P), "\n")
#'
#' # Example 3: Multi-domain alignment
#' X3 <- matrix(rnorm(n * 4), n, 4)
#' X3[16:30, ] <- X3[16:30, ] + 1.5
#'
#' hd_multi <- hyperdesign(list(
#'   modality1 = multidesign(X1, design),
#'   modality2 = multidesign(X2, design),
#'   modality3 = multidesign(X3, design)
#' ))
#'
#' # Compute all pairwise alignments
#' result_multi <- fpgw(hd_multi, omega1 = 0.2, max_iter = 50)
#'
#' # Examine distance matrix
#' print(result_multi$distances)
#' }
#'
#' @export
fpgw <- function(data, ...) {
  UseMethod("fpgw")
}

validate_fpgw_lambda <- function(lambda) {
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      is.na(lambda) || !is.finite(lambda) || lambda < 0) {
    stop("`lambda` must be one finite non-negative number", call. = FALSE)
  }
  if (lambda > 0) {
    stop(
      "`lambda > 0` is withdrawn: the former solver used a balanced ",
      "transport oracle and could not produce a partial coupling. Use `rho` ",
      "for a fixed transported-mass budget.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @rdname fpgw
#' @method fpgw hyperdesign
#' @param omega1 Weight of the feature term (0 ≤ omega1 ≤ 1). Default: 0.001
#' @param lambda Compatibility parameter that must be zero. Use \code{rho} for
#'   partial transport. Default: 0
#' @param rho Mass budget in (0, min(|μ|,|ν|)] for mass-constrained variant
#' @param epsilon Entropic regularization of the warm start (>= 0; 0 disables
#'   it). Default: 0.01
#' @param metric Distance metric: "euclidean", "manhattan" or "cosine".
#'   Default: "euclidean"
#' @param max_iter Maximum iterations for Frank-Wolfe (>= 1). Default: 200
#' @param tol Convergence tolerance for Frank-Wolfe gap (> 0). Default: 1e-6
#' @param inner_max_iter Number of Sinkhorn iterations of the entropic warm
#'   start (>= 1). Default: 50
#' @param verbose Print convergence information. Default: FALSE
#' @section Memory:
#' Every pair of domains stores a dense \eqn{n_i \times n_j} plan and the
#' solver keeps dense within-domain distance matrices. The estimated size is
#' checked against \code{getOption("manifoldalign.ot_max_dense_bytes")}
#' (default 4e9 bytes) before any work starts; larger problems are refused with
#' an error that reports the estimate.
#' @export
fpgw.hyperdesign <- function(data,
                            omega1 = 0.001,
                            lambda = 0,
                            rho = NULL,
                            epsilon = 1e-2,
                            metric = "euclidean",
                            max_iter = 200,
                            tol = 1e-6,
                            inner_max_iter = 50,
                            verbose = FALSE,
                            ...) {

  if (!inherits(data, "hyperdesign")) {
    stop("fpgw requires a hyperdesign object", call. = FALSE)
  }
  .ot_reject_dots(list(...), "fpgw")
  .fpgw_check_args(omega1 = omega1, lambda = lambda, rho = rho,
                   epsilon = epsilon, metric = metric, max_iter = max_iter,
                   tol = tol, inner_max_iter = inner_max_iter,
                   verbose = verbose)
  if (!is.null(rho) && rho > 1) {
    stop("`rho` must be one finite number in (0, 1]", call. = FALSE)
  }

  # Prepare data using shared utilities
  prep <- prepare_ot_data(data)
  X_list <- prep$X_list
  n_samples <- prep$n_samples
  n_domains <- prep$n_domains
  domain_names <- prep$domain_names

  .ot_check_dense_memory(n_samples, "fpgw")

  # With uniform marginals both sides carry mass one.
  if (!is.null(rho) && rho > 1 + 1e-12) {
    stop("rho (", rho, ") exceeds min(sum(p), sum(q)) = 1", call. = FALSE)
  }

  # Compute cost matrices (feature & geometry)
  C_list <- lapply(X_list, function(X) {
    D <- compute_distance_matrix(X, metric)
    D / (max(D) + 1e-10)
  })

  # Pairwise loop -> call FW solver
  n_pairs <- n_domains * (n_domains - 1) / 2
  transport_list <- vector("list", n_pairs)
  dist_mat <- matrix(0, n_domains, n_domains)
  conv <- matrix(TRUE, n_domains, n_domains)

  pair <- 1
  for (i in 1:(n_domains - 1)) {
    for (j in (i + 1):n_domains) {
      if (verbose) cat(sprintf("FPGW: Computing distance between domain %d and %d\n", i, j))

      # Uniform marginals (can be customized later)
      p <- rep(1 / n_samples[i], n_samples[i])
      q <- rep(1 / n_samples[j], n_samples[j])

      # Short-circuit perfectly identical domains in the classical setting
      if (lambda == 0 && is.null(rho) &&
          n_samples[i] == n_samples[j] &&
          isTRUE(all.equal(X_list[[i]], X_list[[j]], tolerance = 1e-10))) {
        ident_plan <- diag(1 / n_samples[i], n_samples[i])
        transport_list[[pair]] <- ident_plan
        dist_mat[i, j] <- dist_mat[j, i] <- 0
        conv[i, j] <- conv[j, i] <- TRUE
        pair <- pair + 1
        next
      }

      C_feat <- compute_feature_cost(X_list[[i]], X_list[[j]], metric)
      Cx <- C_list[[i]]
      Cy <- C_list[[j]]

      ## Feasible starting plan: the (rounded) entropic warm start lies in
      ## Pi(p, q); the mass-constrained variant scales it to mass rho.
      P0 <- outer(p, q)
      if (epsilon > 0) {
        P0 <- sinkhorn_fgw_warmstart(C_feat, Cx, Cy, p, q, omega1, epsilon,
                                     inner_max_iter)
      }
      if (!is.null(rho)) {
        P0 <- P0 * (rho / sum(p))
      }

      ## ---- SOLVER ---- Frank-Wolfe
      sol <- fw_fpgw(C_feat, Cx, Cy, p, q,
                     omega1 = omega1,
                     lambda = lambda,
                     rho = rho,
                     P0 = P0,
                     max_iter = max_iter,
                     tol = tol,
                     verbose = verbose)

      transport_list[[pair]] <- sol$P
      dist_mat[i, j] <- dist_mat[j, i] <- sol$distance
      conv[i, j] <- conv[j, i] <- sol$converged
      pair <- pair + 1
    }
  }

  structure(
    list(
      transport_plans = transport_list,
      distances = dist_mat,
      converged = conv,
      n_samples = n_samples,
      omega1 = omega1,
      lambda = lambda,
      rho = rho,
      domain_names = domain_names,
      training_data = X_list,
      metric = metric
    ),
    class = "fpgw"
  )
}

# Argument validation for fpgw.hyperdesign().
.fpgw_check_args <- function(omega1, lambda, rho, epsilon, metric, max_iter,
                             tol, inner_max_iter, verbose) {
  is_num1 <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!is_num1(omega1) || omega1 < 0 || omega1 > 1) {
    stop("omega1 must be a single number in [0, 1]", call. = FALSE)
  }
  validate_fpgw_lambda(lambda)
  if (!is.null(rho) && (!is_num1(rho) || rho <= 0)) {
    stop("rho must be NULL or a single positive number", call. = FALSE)
  }
  if (!is_num1(epsilon) || epsilon < 0) {
    stop("epsilon must be a single non-negative number (0 disables the warm start)",
         call. = FALSE)
  }
  if (!is.character(metric) || length(metric) != 1L ||
      !metric %in% c("euclidean", "manhattan", "cosine")) {
    stop("metric must be one of 'euclidean', 'manhattan' or 'cosine'",
         call. = FALSE)
  }
  if (!is_num1(max_iter) || max_iter < 1 || max_iter != round(max_iter)) {
    stop("max_iter must be a positive integer", call. = FALSE)
  }
  if (!is_num1(tol) || tol <= 0) {
    stop("tol must be a single positive number", call. = FALSE)
  }
  if (!is_num1(inner_max_iter) || inner_max_iter < 1 ||
      inner_max_iter != round(inner_max_iter)) {
    stop("inner_max_iter must be a positive integer", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE", call. = FALSE)
  }
  invisible(TRUE)
}

# Check that P is feasible for the FPGW variant; error otherwise.
.fpgw_check_plan <- function(P, p, q, variant, rho = NULL, tol = 1e-8,
                             what = "plan") {
  if (!is.matrix(P) || nrow(P) != length(p) || ncol(P) != length(q)) {
    stop("FPGW ", what, " has the wrong dimensions", call. = FALSE)
  }
  if (any(!is.finite(P))) {
    stop("FPGW ", what, " contains non-finite entries", call. = FALSE)
  }
  scale <- max(1, sum(p), sum(q))
  tol <- tol * scale
  rs <- rowSums(P)
  cs <- colSums(P)
  problems <- character(0)
  if (min(P) < -tol) problems <- c(problems, sprintf("negative entry %.3g", min(P)))
  if (variant == "balanced") {
    err <- max(abs(rs - p), abs(cs - q))
    if (err > tol) problems <- c(problems, sprintf("marginal error %.3g", err))
  } else {
    excess <- max(rs - p, cs - q)
    if (excess > tol) problems <- c(problems, sprintf("marginal excess %.3g", excess))
    if (variant == "mass") {
      merr <- abs(sum(P) - rho)
      if (merr > tol) {
        problems <- c(problems, sprintf("total mass %.12g differs from rho = %.12g",
                                        sum(P), rho))
      }
    }
  }
  if (length(problems)) {
    stop("FPGW ", what, " is infeasible for the ", variant, " variant: ",
         paste(problems, collapse = "; "), call. = FALSE)
  }
  invisible(TRUE)
}

#' Frank-Wolfe solver for FPGW
#'
#' @return A list with transport plan P, distance, objective, convergence
#'   status, iterations, final Frank-Wolfe gap and final gradient. The
#'   objective trace is attached as attribute \code{"obj_trace"}.
#'
#' @details
#' Minimises
#' \eqn{f(\gamma) = \omega_1\langle C,\gamma\rangle + \omega_2\langle M(\gamma),\gamma\rangle}
#' where \eqn{M(\gamma) = C_X^2 \gamma 1 \oplus C_Y^2 \gamma^T 1 - 2 C_X \gamma C_Y^T}.
#' The gradient is \eqn{\omega_1 C + 2\omega_2 M(\gamma)}.
#' The feasible set and linear oracle depend on the variant: \eqn{\Pi(p, q)}
#' with \code{classical_ot_lp()} (\code{rho = NULL}) and
#' \eqn{\{\gamma \le (p, q), |\gamma| = \rho\}} with \code{partial_ot_mass()}
#' (\code{rho}). \code{lambda} is retained in this internal signature for
#' compatibility and must be zero. \code{P0} must lie in the feasible set. The
#' step size minimises the quadratic \eqn{f(P + \alpha d)} exactly over
#' \eqn{[0, 1]}, so the objective never increases. Convergence means a
#' Frank-Wolfe gap \eqn{-\langle\nabla f(P), d\rangle} of at most \code{tol}.
#'
#' @keywords internal
fw_fpgw <- function(C, Cx, Cy, p, q,
                    omega1, lambda, rho,
                    P0,
                    max_iter = 200,
                    tol = 1e-6,
                    verbose = FALSE,
                    ...) {
  .ot_reject_dots(list(...), "fw_fpgw")
  validate_fpgw_lambda(lambda)

  variant <- if (!is.null(rho)) "mass" else "balanced"
  n1 <- length(p)
  n2 <- length(q)
  P <- as.matrix(P0)
  .fpgw_check_plan(P, p, q, variant, rho, what = "starting plan")

  Cx2 <- Cx^2
  Cy2 <- Cy^2
  feature_weight <- omega1
  struct_weight <- 1 - omega1

  gw_M <- function(P) {
    outer(as.vector(Cx2 %*% rowSums(P)), as.vector(Cy2 %*% colSums(P)), "+") -
      2 * (Cx %*% P %*% t(Cy))
  }
  objective <- function(P, M) {
    sum(P * (feature_weight * C + struct_weight * M))
  }
  oracle <- function(G) {
    switch(variant,
      balanced = classical_ot_lp(G - min(G), p, q),
      mass = partial_ot_mass(G - min(G), p, q, rho)
    )
  }

  M_P <- gw_M(P)
  obj <- objective(P, M_P)
  obj_trace <- obj
  converged <- FALSE
  gap <- NA_real_
  grad <- NULL
  iter <- 0L

  while (TRUE) {
    grad <- feature_weight * C + struct_weight * 2 * M_P

    gamma_star <- oracle(grad)
    dP <- gamma_star - P
    b <- sum(grad * dP)
    gap <- -b

    if (verbose && iter %% 10 == 0) {
      cat(sprintf("  Iteration %3d: objective = %.6e, FW gap = %.3e\n", iter, obj, gap))
    }
    if (gap <= tol) {
      converged <- TRUE
      break
    }
    if (iter >= max_iter) break
    iter <- iter + 1L

    # f(P + alpha dP) = f(P) + b alpha + a alpha^2
    M_d <- gw_M(dP)
    a <- struct_weight * sum(M_d * dP)

    candidates <- c(0, 1)
    if (a > 0) candidates <- c(candidates, min(1, max(0, -b / (2 * a))))
    phi <- b * candidates + a * candidates^2
    alpha <- candidates[which.min(phi)]
    if (alpha <= 0) {
      # Cannot happen for a positive gap (phi(1) or the interior minimiser is
      # negative); guard against floating-point ties anyway.
      break
    }

    P <- P + alpha * dP
    M_P <- M_P + alpha * M_d
    obj <- objective(P, M_P)
    obj_trace <- c(obj_trace, obj)
  }

  if (!converged && verbose) message("  Did not converge within max_iter iterations")

  .fpgw_check_plan(P, p, q, variant, rho, what = "returned plan")

  # Recompute from scratch for the reported value (avoids incremental drift).
  obj <- objective(P, gw_M(P))
  distance <- sqrt(max(0, obj))
  if (distance < 1e-8) {
    distance <- 0
  }

  result <- list(P = P, distance = distance, objective = obj,
                 converged = converged, iterations = iter, gap = gap,
                 grad = grad, variant = variant)
  attr(result, "obj_trace") <- obj_trace
  result
}

#' Mass-constrained partial OT oracle
#'
#' Exact solution of \eqn{\min \langle C, \gamma\rangle} over
#' \eqn{\gamma \ge 0}, \eqn{\gamma 1 \le p}, \eqn{\gamma^T 1 \le q},
#' \eqn{|\gamma| = \rho} (network simplex on an extended graph).
#' @return A numeric matrix representing the optimal transport plan.
#' @keywords internal
partial_ot_mass <- function(cost, p, q, rho) {
  partial_ot_mass_rcpp(as.matrix(cost), as.numeric(p), as.numeric(q), rho)
}

#' Classical OT linear program
#'
#' Exact network-simplex solution of \eqn{\min \langle C, \gamma\rangle} over
#' \eqn{\Pi(p, q)}. Solver failures are errors (no silent fallback).
#' @return A numeric matrix representing the optimal transport plan.
#' @keywords internal
classical_ot_lp <- function(cost, p, q) {
  network_simplex_ot_cpp(as.matrix(cost), as.numeric(p), as.numeric(q))
}

# Dense LP size guard shared by the lpSolve reference implementations.
.ot_lp_size_guard <- function(n, m) {
  nvar <- as.numeric(n) * as.numeric(m)
  limit <- getOption("manifoldalign.lp_max_variables", 250000)
  if (nvar > limit) {
    stop(sprintf(paste0(
      "The lpSolve transport LP would have %.0f variables (%d x %d), above ",
      "getOption('manifoldalign.lp_max_variables') = %.0f. Use the network ",
      "simplex solver (classical_ot_lp()/partial_ot_mass()) or raise the option."),
      nvar, as.integer(n), as.integer(m), limit), call. = FALSE)
  }
  invisible(nvar)
}

# Sparse (row, column, value) constraint triplets for the n x m transport LP
# with column-major variables x[i + (j - 1) n]: rows 1..n are row sums,
# rows n+1..n+m are column sums.
.ot_lp_marginal_triplets <- function(n, m) {
  var <- seq_len(n * m)
  i <- (var - 1L) %% n + 1L
  j <- (var - 1L) %/% n + 1L
  rbind(cbind(i, var, 1), cbind(n + j, var, 1))
}

.ot_lp_solve <- function(f.obj, dir, rhs, dense, what) {
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Install 'lpSolve' for the ", what, call. = FALSE)
  }
  lp <- lpSolve::lp("min", f.obj, const.dir = dir, const.rhs = rhs,
                    dense.const = dense)
  if (lp$status != 0) {
    stop("lpSolve failed to solve the ", what, " (status ", lp$status, ")",
         call. = FALSE)
  }
  sol <- lp$solution
  sol[abs(sol) < 1e-12] <- 0
  sol
}

#' Mass-constrained partial OT via lpSolve (reference implementation)
#' @return A numeric matrix representing the optimal transport plan.
#' @keywords internal
partial_ot_mass_lp_r <- function(cost, p, q, rho) {
  n <- length(p)
  m <- length(q)
  .ot_lp_size_guard(n, m)
  dense <- rbind(.ot_lp_marginal_triplets(n, m),
                 cbind(n + m + 1, seq_len(n * m), 1))
  sol <- .ot_lp_solve(as.vector(cost), c(rep("<=", n + m), "="), c(p, q, rho),
                      dense, "partial-OT LP")
  matrix(sol, n, m)
}

#' Classical OT linear program (lpSolve reference implementation)
#' @return A numeric matrix representing the optimal transport plan.
#' @keywords internal
classical_ot_lp_r <- function(cost, p, q) {
  n <- length(p)
  m <- length(q)
  .ot_lp_size_guard(n, m)
  sol <- .ot_lp_solve(as.vector(cost), rep("=", n + m), c(p, q),
                      .ot_lp_marginal_triplets(n, m), "OT LP")
  matrix(sol, n, m)
}

# Round a nonnegative matrix onto Pi(p, q) (Altschuler, Weed & Rigollet 2017):
# scale down rows and columns that exceed their marginals, then add the rank-one
# correction outer(err_r, err_c) / sum(err_r). Requires sum(p) == sum(q).
.round_to_transport_polytope <- function(P, p, q) {
  r <- rowSums(P)
  x <- ifelse(r > p, p / r, 1)
  P <- P * x
  cc <- colSums(P)
  y <- ifelse(cc > q, q / cc, 1)
  P <- sweep(P, 2, y, "*")
  err_r <- pmax(p - rowSums(P), 0)
  err_c <- pmax(q - colSums(P), 0)
  s <- sum(err_r)
  if (s > 0) {
    P <- P + outer(err_r, err_c) / s
  }
  P
}

#' Withdrawn TV-penalized transport oracle
#' @return This function always errors; use \code{partial_ot_mass()}.
#' @keywords internal
partial_ot_tv <- function(cost, p, q, lambda) {
  stop(
    "`partial_ot_tv()` is withdrawn because it did not implement the FPGW ",
    "penalized subcoupling objective. Use `partial_ot_mass()` with an explicit ",
    "`rho` budget.",
    call. = FALSE
  )
}

#' Warm-start with entropic OT on the feature cost
#'
#' Log-domain Sinkhorn on \code{omega1 * Cfeat} with regularisation \code{eps},
#' followed by rounding onto \eqn{\Pi(p, q)} so that the result is a feasible
#' Frank-Wolfe starting point.
#' @return A numeric matrix in \eqn{\Pi(p, q)}.
#' @keywords internal
sinkhorn_fgw_warmstart <- function(Cfeat, Cx, Cy, p, q, omega1, eps, iters) {
  K_log <- -(omega1 * Cfeat) / eps
  log_p <- log(p)
  log_q <- log(q)
  lse_rows <- function(A) {
    mx <- apply(A, 1L, max)
    mx + log(rowSums(exp(A - mx)))
  }
  f <- numeric(length(p))
  g <- numeric(length(q))
  for (t in seq_len(iters)) {
    f <- log_p - lse_rows(sweep(K_log, 2, g, "+"))
    g <- log_q - lse_rows(t(K_log + f))
  }
  P <- exp(sweep(K_log + f, 2, g, "+"))
  .round_to_transport_polytope(P, p, q)
}

#' Print Method for FPGW Objects
#'
#' @param x An fpgw object
#' @param ... Additional arguments (currently unused)
#' @return The fpgw object x, invisibly.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' print(res)
#' }
#' @method print fpgw
#' @export
print.fpgw <- function(x, ...) {
  cat("Fused-Partial Gromov-Wasserstein\n")
  cat("================================\n")
  cat("Number of domains:", length(x$domain_names), "\n")
  cat("Domain names:", paste(x$domain_names, collapse = ", "), "\n")
  cat("Feature weight (omega1):", x$omega1, "\n")

  if (!is.null(x$rho)) {
    cat("Mode: Mass-constrained (rho =", x$rho, ")\n")
  } else if (x$lambda > 0) {
    cat("Mode: Withdrawn legacy TV-penalized result (lambda =", x$lambda, ")\n")
  } else {
    cat("Mode: Classical Fused GW\n")
  }

  cat("\nPairwise distances:\n")
  print(round(x$distances, 4))

  if (!all(x$converged)) {
    cat("\nWarning: Some optimizations did not converge\n")
    cat("Non-converged pairs:\n")
    for (i in 1:(nrow(x$converged)-1)) {
      for (j in (i+1):ncol(x$converged)) {
        if (!x$converged[i,j]) {
          cat(sprintf("  (%s, %s)\n", x$domain_names[i], x$domain_names[j]))
        }
      }
    }
  }

  invisible(x)
}

#' Predict method for FPGW
#'
#' Transport new samples using learned FPGW alignment
#'
#' @importFrom RANN nn2
#'
#' @param object An fpgw object
#' @param newdata New data from one domain to transport
#' @param from Source domain index or name
#' @param to Target domain index or name
#' @param type Type of prediction: "transport" (default), "weights", or "embedding"
#' @param k Number of nearest neighbors used for barycentric interpolation
#' @param ... Additional arguments
#'
#' @return Transported samples, barycentric weights, or embeddings (when implemented)
#'
#' @details
#' This method provides out-of-sample extension for FPGW by:
#' \itemize{
#'   \item Finding nearest neighbors in the source domain
#'   \item Using the transport plan to map to target domain
#'   \item Interpolating based on transport weights
#' }
#' The object must contain the training data used for fitting,
#' otherwise an error is thrown.
#'
#' @examples
#' \dontrun{
#' # After computing FPGW alignment
#' result <- fpgw(hd, omega1 = 0.1)
#'
#' # Transport new samples from domain 1 to domain 2
#' new_samples <- matrix(rnorm(10 * 5), 10, 5)
#' transported <- predict(result, new_samples, from = 1, to = 2)
#' head(transported)
#' }
#'
#' @export
predict.fpgw <- function(object, newdata, from, to,
                         type = c("transport", "weights", "embedding"),
                         k = NULL, ...) {
  # Validate inputs
  if (!inherits(object, "fpgw")) {
    stop("object must be an fpgw result", call. = FALSE)
  }

  type <- match.arg(type)

  # Convert domain names to indices if needed
  if (is.character(from)) {
    from_idx <- which(object$domain_names == from)
    if (length(from_idx) == 0) {
      stop("Domain '", from, "' not found", call. = FALSE)
    }
    from <- from_idx
  }

  if (is.character(to)) {
    to_idx <- which(object$domain_names == to)
    if (length(to_idx) == 0) {
      stop("Domain '", to, "' not found", call. = FALSE)
    }
    to <- to_idx
  }

  # Get the appropriate transport plan
  n_domains <- length(object$domain_names)
  if (from < to) {
    pair_idx <- (from - 1) * n_domains -
                from * (from - 1) / 2 + (to - from)
  } else {
    pair_idx <- (to - 1) * n_domains -
                to * (to - 1) / 2 + (from - to)
  }

  if (pair_idx > length(object$transport_plans)) {
    stop("Transport plan not found for specified domain pair", call. = FALSE)
  }

  P <- object$transport_plans[[pair_idx]]

  # If from > to, we need to transpose the plan
  if (from > to) {
    P <- t(P)
  }

  if (type %in% c("transport", "weights")) {
    if (is.null(object$training_data)) {
      stop("Training data not stored in object; cannot perform barycentric mapping",
           call. = FALSE)
    }

    X_source <- object$training_data[[from]]
    if (ncol(newdata) != ncol(X_source)) {
      stop("New data has different number of features than source domain",
           call. = FALSE)
    }

    if (is.null(k)) {
      k <- min(5L, nrow(X_source))
    } else {
      k <- as.integer(k)
      if (is.na(k) || k <= 0) {
        stop("k must be a positive integer", call. = FALSE)
      }
      k <- min(k, nrow(X_source))
    }
    if (k < 1) {
      stop("k must be at least 1", call. = FALSE)
    }

    nn <- RANN::nn2(X_source, newdata, k = k)
    idx <- as.matrix(nn$nn.idx)
    dists <- as.matrix(nn$nn.dists)
    w <- 1 / (dists + 1e-8)
    w_den <- rowSums(w)
    w_den[w_den == 0] <- 1
    w <- w / w_den

    n_new <- nrow(newdata)
    plan_mass <- rowSums(P)
    safe_mass <- plan_mass
    safe_mass[safe_mass == 0] <- 1
    normalized_plan <- sweep(P, 1, safe_mass, "/")
    normalized_plan[plan_mass == 0, ] <- 0

    bary_weights <- matrix(0, n_new, ncol(P))
    for (i in seq_len(n_new)) {
      plan_rows <- normalized_plan[idx[i, ], , drop = FALSE]
      bary_weights[i, ] <- w[i, ] %*% plan_rows
      weights_sum <- sum(bary_weights[i, ])
      if (weights_sum > 0) {
        bary_weights[i, ] <- bary_weights[i, ] / weights_sum
      } else {
        bary_weights[i, ] <- rep(1 / ncol(P), ncol(P))
      }
    }

    if (type == "weights") {
      return(bary_weights)
    }

    X_target <- object$training_data[[to]]
    transported <- bary_weights %*% X_target
    return(transported)
  }

  stop("Embedding prediction not yet implemented for FPGW", call. = FALSE)
}
