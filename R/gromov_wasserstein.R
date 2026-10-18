# Generic is defined in all_generic.R
# Note: compute_distance_matrix has been moved to utils-ot.R

#' @rdname gromov_wasserstein
#' @method gromov_wasserstein hyperdesign
#' @param epsilon Entropic regularization parameter. Default: 0.1
#' @param metric Distance metric for within-domain distances. Default: "euclidean"
#' @param init Initialization method: "uniform" or "random". Default: "uniform"
#' @param max_iter Maximum iterations for outer loop. Default: 100
#' @param tol Convergence tolerance for outer loop. Default: 1e-9
#' @param verbose Print convergence information. Default: FALSE
#' @param marginals Optional list of marginal distributions per domain
#' @param scale Scaling strategy: "per_domain", "global", or "none". Default: "per_domain"
#' @param inner_max_iter Maximum Sinkhorn iterations. Default: 30
#' @param inner_tol Sinkhorn convergence tolerance. Default: 1e-9
#' @param dist_chunk_size Chunk size for distance computation. Default: 1000
#' @details
#' For each pair of domains the solver minimises the entropic GW objective
#' \deqn{E(P) = \sum_{i,j,k,l} (C^1_{ik} - C^2_{jl})^2 P_{ij} P_{kl} +
#'   \epsilon \sum_{ij} P_{ij} \log P_{ij}}
#' over couplings with marginals \code{p} and \code{q} by the projected
#' mirror-descent scheme of Peyre, Cuturi & Solomon (2016): each outer step
#' solves an entropic OT problem whose cost is the gradient \eqn{2M(P)} of the
#' quadratic term, with regularisation \code{epsilon}. Its fixed points are
#' stationary points of \eqn{E} with the same \code{epsilon}, and the
#' reported \code{losses} are \eqn{E(P)}. The inner entropic OT problems use a
#' log-stabilised Sinkhorn iteration; if it cannot produce a finite coupling
#' the function stops with an error rather than substituting another plan.
#' \code{distances} are \eqn{\sqrt{\langle M(P), P\rangle}} (without the
#' entropy term).
#'
#' The returned object is not a projector and does not inherit from
#' \code{multiblock_biprojector}. Dense plans for all pairs are held in memory;
#' the estimated size is checked against
#' \code{getOption("manifoldalign.ot_max_dense_bytes")} (default 4e9 bytes)
#' before any work starts.
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
  
  .ot_reject_dots(list(...), "gromov_wasserstein")

  X_list <- lapply(data, function(d) d$x)
  n_samples <- sapply(X_list, nrow)
  .ot_check_dense_memory(n_samples, "gromov_wasserstein")

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
      domain_names = domain_names,
      domain_data = X_list,
      metric = metric
    ),
    class = "gromov_wasserstein"
  )
}

# Row-wise log-sum-exp of a matrix.
.lse_rows <- function(A) {
  mx <- apply(A, 1L, max)
  mx + log(rowSums(exp(A - mx)))
}

# Entropic OT  min <cost, P> + eps * sum(P log P)  over Pi(p, q), solved by a
# stabilised Sinkhorn iteration. Dual potentials (f, g) are kept in the log
# domain and the kernel exp((f_i + g_j - cost_ij) / eps) is rebuilt (scalings
# absorbed into the potentials) whenever the scalings leave [e^-50, e^50] or
# become non-finite; the exact log-domain updates used on rebuild cannot
# underflow. `g_init` warm-starts the column potential. Returns the coupling,
# the potentials, the final column-marginal error, and the iteration count.
.sinkhorn_stabilized <- function(cost, p, q, eps, g_init = NULL,
                                 max_iter = 1000, tol = 1e-9) {
  log_p <- log(p)
  log_q <- log(q)
  g <- if (is.null(g_init)) numeric(length(q)) else g_init

  f_update <- function(g) eps * (log_p - .lse_rows(sweep(-cost, 2, g, "+") / eps))
  g_update <- function(f) eps * (log_q - .lse_rows(t((f - cost) / eps)))
  kernel <- function(f, g) exp((sweep(-cost, 2, g, "+") + f) / eps)
  ok <- function(x) all(is.finite(x)) && all(x > 0) && max(abs(log(x))) < 50

  f <- f_update(g)
  K <- kernel(f, g)  # row sums are exactly p, so no entry exceeds 1
  u <- rep(1, length(p))
  v <- rep(1, length(q))
  err <- Inf
  it <- 0L
  while (it < max_iter) {
    it <- it + 1L
    v_new <- q / as.vector(crossprod(K, u))
    u_new <- p / as.vector(K %*% v_new)
    if (!ok(u_new) || !ok(v_new)) {
      # absorb the last good scalings, then exact log-domain updates
      f <- f + eps * log(u)
      g <- g + eps * log(v)
      g <- g_update(f)
      f <- f_update(g)
      K <- kernel(f, g)
      u <- rep(1, length(p))
      v <- rep(1, length(q))
    } else {
      u <- u_new
      v <- v_new
    }
    # rows are exact after the u-update; measure the column error
    err <- sum(abs(v * as.vector(crossprod(K, u)) - q))
    if (err <= tol) break
  }

  P <- K * outer(u, v)
  if (any(!is.finite(P))) {
    stop("Entropic GW: the inner Sinkhorn iteration produced a non-finite ",
         "coupling (epsilon = ", eps, ")", call. = FALSE)
  }
  list(P = P, f = f + eps * log(u), g = g + eps * log(v), err = err,
       iterations = it)
}

# Entropic Gromov-Wasserstein (square loss) by projected mirror descent
# (Peyre, Cuturi & Solomon 2016): P <- argmin_{Pi(p,q)} <2 M(P_old), P> +
# epsilon * sum(P log P), where 2 M(P) is the gradient of <M(P), P>. Fixed
# points are stationary for the reported loss <M(P), P> + epsilon * sum(P log P).
solve_entropic_gw <- function(C1, C2, p, q, P_init,
                              epsilon = 0.1,
                              max_iter = 100,
                              tol = 1e-9,
                              inner_max_iter = 30,
                              inner_tol = 1e-9,
                              verbose = FALSE) {

  n1 <- nrow(C1)
  n2 <- nrow(C2)

  P <- P_init

  C1_sqr <- C1^2
  C2_sqr <- C2^2

  constC1 <- as.vector(C1_sqr %*% p)
  constC2 <- as.vector(C2_sqr %*% q)
  const_term <- outer(constC1, rep(1, n2)) + outer(rep(1, n1), constC2)
  gw_M <- function(P) const_term - 2 * (C1 %*% P %*% t(C2))
  neg_entropy <- function(P) sum(P[P > 0] * log(P[P > 0]))

  M <- gw_M(P)
  loss_old <- Inf
  g_warm <- NULL
  converged <- FALSE
  iter <- 0L

  for (iter in seq_len(max_iter)) {
    sk <- .sinkhorn_stabilized(2 * M, p, q, epsilon, g_init = g_warm,
                               max_iter = inner_max_iter, tol = inner_tol)
    P <- sk$P
    g_warm <- sk$g

    M <- gw_M(P)
    loss_new <- sum(P * M) + epsilon * neg_entropy(P)

    if (verbose && iter %% 10 == 0) {
      message(sprintf("  Iteration %3d: loss = %.6e", iter, loss_new))
    }

    rel_change <- abs(loss_old - loss_new) / (abs(loss_old) + 1e-15)
    loss_old <- loss_new
    if (is.finite(rel_change) && rel_change < tol) {
      converged <- TRUE
      if (verbose) {
        message(sprintf("  Converged at iteration %d", iter))
      }
      break
    }
  }

  if (!converged && verbose) {
    message("  Did not converge within max_iter iterations")
  }

  list(
    P = P,
    distance = sqrt(max(0, sum(P * M))),
    loss = loss_old,
    converged = converged,
    iterations = iter
  )
}


#' @method print gromov_wasserstein
#' @export
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
#' Transport new samples using the learned GW alignment.
#' 
#' @param object A gromov_wasserstein object
#' @param newdata New data from the source domain (samples x features)
#' @param from Source domain index or name
#' @param to Target domain index or name
#' @param type Prediction type. `"weights"` (default) returns barycentric
#'   weights over target-domain samples: the inverse-distance weighted
#'   average of the row-normalised coupling rows of the `k` nearest
#'   source-domain training samples. `"transport"` returns those weights
#'   multiplied by the target-domain features (the barycentric projection), so
#'   both types always agree; there is no implicit fallback to another map.
#' @param k Number of nearest neighbours used to build the barycentric
#'   combination (default: 5).
#' @param ... Reserved for future use.
#' 
#' @return Matrix of barycentric weights (`type = "weights"`) or transported
#'   samples (`type = "transport"`).
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' fit <- gromov_wasserstein(hd, epsilon = 0.5, max_iter = 10)
#' newX <- matrix(rnorm(6), 2, 3)
#' pred <- predict(fit, newX, from = 1, to = 2)
#' }
#' @export
predict.gromov_wasserstein <- function(object, newdata, from, to,
                                       type = c("weights", "transport"),
                                       k = 5, ...) {
  if (!inherits(object, "gromov_wasserstein")) {
    stop("object must be a gromov_wasserstein result", call. = FALSE)
  }

  type <- match.arg(type)

  if (missing(from) || missing(to)) {
    stop("Both 'from' and 'to' domains must be specified", call. = FALSE)
  }

  from_idx <- resolve_domain_index(from, object$domain_names)
  to_idx <- resolve_domain_index(to, object$domain_names)

  source_data <- object$domain_data[[from_idx]]
  target_data <- object$domain_data[[to_idx]]

  if (is.null(source_data)) {
    stop("Source domain data not available for out-of-sample prediction", call. = FALSE)
  }

  if (!is.matrix(newdata)) {
    newdata <- as.matrix(newdata)
  }

  if (ncol(newdata) != ncol(source_data)) {
    stop("newdata has ", ncol(newdata), " columns but source domain has ",
         ncol(source_data), call. = FALSE)
  }

  plan <- extract_transport_plan(object$transport_plans, from_idx, to_idx,
                                 length(object$domain_names))
  plan_rows <- normalize_plan_rows(plan)

  metric <- if (!is.null(object$metric)) object$metric else "euclidean"
  nn <- find_knn_indices(newdata, source_data, k = k, metric = metric)
  bary_weights <- aggregate_barycentric_weights(plan_rows, nn$idx, nn$dists)

  if (type == "weights" || is.null(target_data)) {
    if (!is.null(target_data)) {
      colnames(bary_weights) <- rownames(target_data)
    }
    return(bary_weights)
  }

  # The transport map is always the barycentric projection through the
  # learned coupling, so type = "transport" equals the "weights" output
  # multiplied by the target-domain features.
  transported <- bary_weights %*% as.matrix(target_data)
  colnames(transported) <- colnames(target_data)
  transported
}

#' @rdname gromov_wasserstein
#' @method gromov_wasserstein default
#' @export
gromov_wasserstein.default <- function(data, ...) {
  stop("gromov_wasserstein requires a hyperdesign object", call. = FALSE)
}
