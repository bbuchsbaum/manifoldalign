#' OT + Procrustes aligner adapter for align_many()
#'
#' Provides a pairwise interface implementing the orthogonal-invariant OT
#' objective of Alvarez-Melis et al. (AISTATS 2019), specialized to learning an
#' orthogonal map between feature spaces.

#' Construct an OT-Procrustes aligner descriptor
#'
#' This aligner estimates an orthogonal transform between two domains by
#' alternating between an entropic OT coupling (Sinkhorn) and a Procrustes
#' update (SVD) under the current coupling. It is intended as a lightweight,
#' initialization-robust alternative to adversarial alignment for feature
#' spaces that are identifiable only up to global orthogonal transforms.
#'
#' @return An object of class `c("ot_procrustes_aligner", "aligner")`.
#' @examples
#' algo <- ot_procrustes_aligner()
#' aligner_capabilities(algo)
#' @export
ot_procrustes_aligner <- function() {
  new_aligner("ot_procrustes", group = "O", supports_multi = FALSE)
}

.ot_procrustes_random_orthogonal <- function(d, force_SO = FALSE) {
  if (d <= 0) stop("d must be positive", call. = FALSE)
  if (d == 1L) {
    return(matrix(if (isTRUE(force_SO)) 1 else sample(c(-1, 1), 1), 1, 1))
  }
  A <- matrix(stats::rnorm(d * d), d, d)
  sv <- tryCatch(svd(A), error = function(e) NULL)
  if (is.null(sv)) {
    return(diag(d))
  }
  R <- sv$u %*% t(sv$v)
  if (isTRUE(force_SO) && det(R) < 0) {
    sv$u[, d] <- -sv$u[, d]
    R <- sv$u %*% t(sv$v)
  }
  R
}

.ot_procrustes_normalize_marginal <- function(w, n, name) {
  if (is.null(w)) {
    return(rep(1 / n, n))
  }
  w <- as.numeric(w)
  if (length(w) != n) {
    stop("`", name, "` must have length ", n, ".", call. = FALSE)
  }
  if (any(!is.finite(w))) stop("`", name, "` must be finite.", call. = FALSE)
  if (any(w < 0)) stop("`", name, "` must be nonnegative.", call. = FALSE)
  s <- sum(w)
  if (!is.finite(s) || s <= 0) {
    stop("`", name, "` must have positive total mass.", call. = FALSE)
  }
  w / s
}

.ot_procrustes_distance_signature <- function(X, k) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2) {
    stop("distance signature requires at least 2 points.", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("k must be a positive integer.", call. = FALSE)
  }
  if (k > (n - 1L)) k <- n - 1L

  D <- compute_squared_distances_cpp_fallback(X, X)
  sig <- matrix(0, n, k)
  for (ii in seq_len(n)) {
    r <- D[ii, ]
    # First entry is (usually) the self-distance; drop it.
    sig[ii, ] <- sort.int(r)[2:(k + 1L)]
  }
  sig
}

.ot_procrustes_signature_init_rotation <- function(X_i, X_j, a, b,
                                                  epsilon0, epsilon_min,
                                                  sinkhorn_max_iter, sinkhorn_tol,
                                                  stabilized, cost_scale,
                                                  force_SO) {
  n1 <- nrow(X_i)
  n2 <- nrow(X_j)
  d <- ncol(X_i)

  k_sig <- min(10L, n1 - 1L, n2 - 1L)
  if (!is.finite(k_sig) || k_sig < 1L) return(diag(d))

  eps_sig <- max(epsilon_min, 0.01 * epsilon0)
  if (!is.finite(eps_sig) || eps_sig <= 0) eps_sig <- epsilon_min

  sig_i <- .ot_procrustes_distance_signature(X_i, k_sig)
  sig_j <- .ot_procrustes_distance_signature(X_j, k_sig)
  C_sig <- compute_squared_distances_cpp_fallback(sig_i, sig_j)

  sink <- .ot_procrustes_solve_sinkhorn(
    C_sig,
    a,
    b,
    epsilon = eps_sig,
    max_iter = max(as.integer(sinkhorn_max_iter), 5000L),
    tol = sinkhorn_tol,
    stabilized = stabilized,
    cost_scale = cost_scale
  )
  pi0 <- sink$pi

  if (!all(is.finite(pi0))) {
    return(diag(d))
  }

  M <- t(X_i) %*% pi0 %*% X_j
  sv <- tryCatch(svd(M, nu = d, nv = d), error = function(e) NULL)
  if (is.null(sv) || any(!is.finite(sv$u)) || any(!is.finite(sv$v))) {
    return(diag(d))
  }
  R0 <- sv$u %*% t(sv$v)
  if (isTRUE(force_SO) && det(R0) < 0) {
    sv$u[, d] <- -sv$u[, d]
    R0 <- sv$u %*% t(sv$v)
  }
  R0
}

.ot_procrustes_run <- function(X_i, X_j, a, b,
                               R_init,
                               epsilon0,
                               epsilon_min,
                               decay,
                               max_iter,
                               tol,
                               sinkhorn_max_iter,
                               sinkhorn_tol,
                               stabilized,
                               cost_scale,
                               force_SO,
                               verbose) {
  d <- ncol(X_i)
  eps <- epsilon0
  R <- R_init
  converged <- FALSE
  max_iter <- as.integer(max_iter)
  loss_trace <- numeric(max_iter)
  eps_trace <- numeric(max_iter)
  trace_len <- 0L
  last_pi <- NULL
  log_u <- NULL
  log_v <- NULL

  Xt <- t(X_i)

  for (iter in seq_len(max_iter)) {
    trace_len <- trace_len + 1L
    eps_trace[trace_len] <- eps

    Xr <- X_i %*% R
    C_raw <- compute_squared_distances_cpp_fallback(Xr, X_j)

    sink <- .ot_procrustes_solve_sinkhorn(C_raw, a, b, epsilon = eps,
                                          max_iter = sinkhorn_max_iter,
                                          tol = sinkhorn_tol,
                                          stabilized = stabilized,
                                          cost_scale = cost_scale,
                                          log_u0 = log_u,
                                          log_v0 = log_v)
    pi <- sink$pi
    log_u <- sink$log_u
    log_v <- sink$log_v

    if (!all(is.finite(pi))) {
      warning("ot_procrustes: Sinkhorn produced non-finite plan; using product plan.",
              call. = FALSE)
      pi <- outer(a, b)
    }

    # Procrustes update: maximize trace(R^T X^T pi Y)
    M <- Xt %*% (pi %*% X_j)
    sv <- tryCatch(svd(M, nu = d, nv = d), error = function(e) NULL)
    if (is.null(sv) || any(!is.finite(sv$u)) || any(!is.finite(sv$v))) {
      R_new <- R
    } else {
      R_new <- sv$u %*% t(sv$v)
      if (isTRUE(force_SO) && det(R_new) < 0) {
        sv$u[, d] <- -sv$u[, d]
        R_new <- sv$u %*% t(sv$v)
      }
    }

    dR <- max(abs(R_new - R))
    R <- R_new

    # Loss: expected squared distance under the current plan (no entropy term).
    loss <- sum(pi * (C_raw))
    loss_trace[trace_len] <- loss

    if (isTRUE(verbose)) {
      msg <- sprintf("iter=%d eps=%.4g loss=%.6g dR=%.3g", iter, eps, loss, dR)
      message(msg)
    }

    last_pi <- pi

    if (eps <= epsilon_min + 1e-15) {
      if (is.finite(dR) && dR < tol) {
        converged <- TRUE
        break
      }
    }

    eps <- max(epsilon_min, eps * decay)
  }

  list(
    rotation = R,
    loss = if (trace_len) loss_trace[trace_len] else NA_real_,
    loss_trace = loss_trace[seq_len(trace_len)],
    epsilon_trace = eps_trace[seq_len(trace_len)],
    converged = converged,
    iterations = trace_len,
    transport = last_pi
  )
}

.ot_procrustes_solve_sinkhorn <- function(cost, a, b, epsilon,
                                         max_iter, tol, stabilized,
                                         cost_scale = c("median", "max", "none"),
                                         log_u0 = NULL,
                                         log_v0 = NULL) {
  cost_scale <- match.arg(cost_scale)
  if (!is.matrix(cost)) cost <- as.matrix(cost)
  if (any(!is.finite(cost))) {
    stop("Sinkhorn cost matrix contains non-finite values.", call. = FALSE)
  }
  # Stabilize exponentiation by shifting min(cost) to 0 (keeps OT solution).
  c0 <- min(cost)
  if (is.finite(c0) && c0 != 0) {
    cost <- cost - c0
  }

  scale_factor <- 1
  if (cost_scale == "median") {
    pos <- cost[cost > 0]
    s <- stats::median(pos, na.rm = TRUE)
    if (is.finite(s) && s > 0) scale_factor <- s
  } else if (cost_scale == "max") {
    s <- max(cost)
    if (is.finite(s) && s > 0) scale_factor <- s
  }
  if (!is.finite(scale_factor) || scale_factor <= 0) scale_factor <- 1
  if (scale_factor != 1) {
    cost <- cost / scale_factor
  }

  if (isTRUE(stabilized)) {
    res <- sinkhorn_unified_potentials(
      cost,
      a,
      b,
      epsilon = epsilon,
      log_u0 = log_u0,
      log_v0 = log_v0,
      max_iter = as.integer(max_iter),
      tol = tol,
      stabilized = TRUE
    )

    return(list(
      pi = res$pi,
      log_u = res$log_u,
      log_v = res$log_v,
      cost_shift = c0,
      cost_scale = scale_factor
    ))
  }

  pi <- sinkhorn_unified(
    cost,
    a,
    b,
    epsilon = epsilon,
    max_iter = as.integer(max_iter),
    tol = tol,
    stabilized = FALSE
  )

  list(
    pi = pi,
    log_u = NULL,
    log_v = NULL,
    cost_shift = c0,
    cost_scale = scale_factor
  )
}

#' Fit OT-Procrustes on a pair of domains
#'
#' Estimates an orthogonal map between two domains without known correspondences
#' by alternating between an entropic OT coupling and a Procrustes update.
#'
#' @param algo An ot_procrustes_aligner object.
#' @param X_i First domain data matrix (samples x features).
#' @param X_j Second domain data matrix (samples x features).
#' @param links Optional correspondence links (currently unused).
#' @param i Internal domain index used by [align_many()] (unused).
#' @param j Internal domain index used by [align_many()] (unused).
#' @param a Optional source marginal (length nrow(X_i)); defaults to uniform.
#' @param b Optional target marginal (length nrow(X_j)); defaults to uniform.
#' @param epsilon0 Initial entropic regularization (> 0).
#' @param epsilon_min Minimum entropic regularization (> 0).
#' @param decay Multiplicative decay factor in (0,1] applied each iteration.
#' @param n_init Number of random restarts. When `init="identity"` and
#'   `n_init>1`, the method tries the identity, then a deterministic
#'   distance-signature initialization, and uses remaining runs for random
#'   orthogonal restarts. The best run (lowest loss) is returned.
#' @param seed Optional seed for reproducible initialization. When called via
#'   [align_many()], the internal `i` and `j` indices are used to derive a
#'   per-edge seed from this value.
#' @param max_iter Maximum number of outer (alternating) iterations.
#' @param tol Convergence tolerance on max entrywise change in the rotation.
#' @param sinkhorn_max_iter Maximum number of Sinkhorn iterations.
#' @param sinkhorn_tol Sinkhorn convergence tolerance.
#' @param stabilized Logical; use log-domain stabilization in Sinkhorn (slower
#'   but helpful when `epsilon_min` is very small).
#' @param cost_scale Cost normalization used before Sinkhorn: `"median"` divides
#'   by the median positive cost, `"max"` divides by the maximum cost, and
#'   `"none"` disables scaling.
#' @param init Initialization method: `"identity"`, `"signature"`, or
#'   `"random_orthogonal"`.
#' @param force_SO Logical; if TRUE enforce det(R)=+1 for the learned rotation.
#' @param store_transport Logical; if TRUE store the final transport plan.
#' @param verbose Logical; print basic progress information.
#' @param ... Additional arguments (unused).
#' @return An object of class `ot_procrustes_pair_fit` containing the learned
#'   rotation (i->j), loss, convergence diagnostics, and optionally the final
#'   transport plan.
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(150), 50, 3)
#' algo <- ot_procrustes_aligner()
#' fit <- fit_pair(algo, X, X)
#' fit$converged
#' }
#' @export
fit_pair.ot_procrustes_aligner <- function(algo, X_i, X_j, links = NULL,
                                           i = NULL,
                                           j = NULL,
                                           a = NULL,
                                           b = NULL,
                                           epsilon0 = 1,
                                           epsilon_min = 0.05,
                                           decay = 0.95,
                                           n_init = 1L,
                                           seed = NULL,
                                           max_iter = 50,
                                           tol = 1e-6,
                                           sinkhorn_max_iter = 200,
                                           sinkhorn_tol = 1e-6,
                                           stabilized = (epsilon_min < 0.01),
                                           cost_scale = c("median", "max", "none"),
                                           init = c("identity", "signature", "random_orthogonal"),
                                           force_SO = FALSE,
                                           store_transport = FALSE,
                                           verbose = FALSE,
                                           ...) {

  cost_scale <- match.arg(cost_scale)
  init <- match.arg(init)

  chk::chk_number(epsilon0)
  chk::chk_true(epsilon0 > 0)
  chk::chk_number(epsilon_min)
  chk::chk_true(epsilon_min > 0)
  chk::chk_number(decay)
  chk::chk_true(decay > 0 && decay <= 1)
  chk::chk_number(n_init)
  chk::chk_true(n_init >= 1)
  if (!is.null(seed)) {
    chk::chk_number(seed)
    chk::chk_true(is.finite(seed))
  }
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  chk::chk_number(sinkhorn_max_iter)
  chk::chk_true(sinkhorn_max_iter > 0)
  chk::chk_number(sinkhorn_tol)
  chk::chk_true(sinkhorn_tol > 0)
  chk::chk_logical(stabilized)
  chk::chk_true(is.character(cost_scale) && length(cost_scale) == 1L)
  chk::chk_logical(force_SO)
  chk::chk_logical(store_transport)
  chk::chk_logical(verbose)

  X_i <- as.matrix(X_i)
  X_j <- as.matrix(X_j)
  if (!is.numeric(X_i) || !is.numeric(X_j)) {
    stop("X_i and X_j must be numeric matrices.", call. = FALSE)
  }
  if (nrow(X_i) < 1 || nrow(X_j) < 1) {
    stop("X_i and X_j must have at least one row.", call. = FALSE)
  }
  if (ncol(X_i) != ncol(X_j)) {
    stop("ot_procrustes requires matching feature dimensions: ncol(X_i) == ncol(X_j).",
         call. = FALSE)
  }
  if (any(!is.finite(X_i)) || any(!is.finite(X_j))) {
    stop("X_i and X_j must contain only finite values.", call. = FALSE)
  }

  n1 <- nrow(X_i)
  n2 <- nrow(X_j)
  d <- ncol(X_i)

  a <- .ot_procrustes_normalize_marginal(a, n1, "a")
  b <- .ot_procrustes_normalize_marginal(b, n2, "b")

  if (!is.null(seed)) {
    ii <- suppressWarnings(as.integer(i %||% 0L))
    jj <- suppressWarnings(as.integer(j %||% 0L))
    if (!is.finite(ii)) ii <- 0L
    if (!is.finite(jj)) jj <- 0L
    set.seed(as.integer(seed) + 1000L * ii + jj)
  }

  best <- NULL
  best_loss <- Inf
  best_idx <- NA_integer_

  n_init <- as.integer(n_init)
  R_sig <- NULL
  R0_list <- vector("list", n_init)
  for (k in seq_len(n_init)) {
    R0_list[[k]] <- if (init == "identity") {
      if (k == 1L) {
        diag(d)
      } else if (k == 2L) {
        if (is.null(R_sig)) {
          R_sig <- tryCatch(
            .ot_procrustes_signature_init_rotation(
              X_i = X_i,
              X_j = X_j,
              a = a,
              b = b,
              epsilon0 = epsilon0,
              epsilon_min = epsilon_min,
              sinkhorn_max_iter = sinkhorn_max_iter,
              sinkhorn_tol = sinkhorn_tol,
              stabilized = stabilized,
              cost_scale = cost_scale,
              force_SO = force_SO
            ),
            error = function(e) NULL
          )
        }
        if (is.null(R_sig)) diag(d) else R_sig
      } else {
        .ot_procrustes_random_orthogonal(d, force_SO = force_SO)
      }
    } else if (init == "signature") {
      if (k == 1L) {
        if (is.null(R_sig)) {
          R_sig <- tryCatch(
            .ot_procrustes_signature_init_rotation(
              X_i = X_i,
              X_j = X_j,
              a = a,
              b = b,
              epsilon0 = epsilon0,
              epsilon_min = epsilon_min,
              sinkhorn_max_iter = sinkhorn_max_iter,
              sinkhorn_tol = sinkhorn_tol,
              stabilized = stabilized,
              cost_scale = cost_scale,
              force_SO = force_SO
            ),
            error = function(e) NULL
          )
        }
        if (is.null(R_sig)) diag(d) else R_sig
      } else {
        .ot_procrustes_random_orthogonal(d, force_SO = force_SO)
      }
    } else if (init == "random_orthogonal") {
      .ot_procrustes_random_orthogonal(d, force_SO = force_SO)
    } else {
      diag(d)
    }
  }

  # When using the common two-start pattern (identity + signature), do a short
  # warm-up run for each and only continue the better one when clearly better.
  run_indices <- seq_len(n_init)
  warm_runs <- vector("list", n_init)
  cont_eps0 <- NULL
  cont_max_iter <- 0L
  if (n_init == 2L && init == "identity" && !isTRUE(stabilized)) {
    warm_iter <- min(5L, as.integer(max_iter))
    if (warm_iter >= 1L) {
      warm_runs[[1L]] <- .ot_procrustes_run(
        X_i = X_i,
        X_j = X_j,
        a = a,
        b = b,
        R_init = R0_list[[1L]],
        epsilon0 = epsilon0,
        epsilon_min = epsilon_min,
        decay = decay,
        max_iter = warm_iter,
        tol = tol,
        sinkhorn_max_iter = sinkhorn_max_iter,
        sinkhorn_tol = sinkhorn_tol,
        stabilized = stabilized,
        cost_scale = cost_scale,
        force_SO = force_SO,
        verbose = FALSE
      )
      warm_runs[[2L]] <- .ot_procrustes_run(
        X_i = X_i,
        X_j = X_j,
        a = a,
        b = b,
        R_init = R0_list[[2L]],
        epsilon0 = epsilon0,
        epsilon_min = epsilon_min,
        decay = decay,
        max_iter = warm_iter,
        tol = tol,
        sinkhorn_max_iter = sinkhorn_max_iter,
        sinkhorn_tol = sinkhorn_tol,
        stabilized = stabilized,
        cost_scale = cost_scale,
        force_SO = force_SO,
        verbose = FALSE
      )

      s1 <- warm_runs[[1L]]$loss
      s2 <- warm_runs[[2L]]$loss
      score_margin <- 0.02

      if (is.finite(s1) && is.finite(s2)) {
        if (s1 <= s2 * (1 - score_margin)) {
          run_indices <- 1L
        } else if (s2 <= s1 * (1 - score_margin)) {
          run_indices <- 2L
        }
      } else if (is.finite(s1) && !is.finite(s2)) {
        run_indices <- 1L
      } else if (!is.finite(s1) && is.finite(s2)) {
        run_indices <- 2L
      }

      cont_eps0 <- max(epsilon_min, epsilon0 * decay^warm_iter)
      cont_max_iter <- max(0L, as.integer(max_iter) - warm_iter)
    }
  }

  for (k in run_indices) {
    run0 <- warm_runs[[k]]
    if (!is.null(run0) && cont_max_iter > 0L) {
      cont <- .ot_procrustes_run(
        X_i = X_i,
        X_j = X_j,
        a = a,
        b = b,
        R_init = run0$rotation,
        epsilon0 = cont_eps0,
        epsilon_min = epsilon_min,
        decay = decay,
        max_iter = cont_max_iter,
        tol = tol,
        sinkhorn_max_iter = sinkhorn_max_iter,
        sinkhorn_tol = sinkhorn_tol,
        stabilized = stabilized,
        cost_scale = cost_scale,
        force_SO = force_SO,
        verbose = verbose
      )
      run <- cont
      run$loss_trace <- c(run0$loss_trace, run$loss_trace)
      run$epsilon_trace <- c(run0$epsilon_trace, run$epsilon_trace)
      run$iterations <- run0$iterations + run$iterations
    } else if (!is.null(run0)) {
      run <- run0
    } else {
      run <- .ot_procrustes_run(
        X_i = X_i,
        X_j = X_j,
        a = a,
        b = b,
        R_init = R0_list[[k]],
        epsilon0 = epsilon0,
        epsilon_min = epsilon_min,
        decay = decay,
        max_iter = max_iter,
        tol = tol,
        sinkhorn_max_iter = sinkhorn_max_iter,
        sinkhorn_tol = sinkhorn_tol,
        stabilized = stabilized,
        cost_scale = cost_scale,
        force_SO = force_SO,
        verbose = verbose
      )
    }

    if (is.null(best)) {
      best <- run
      best_loss <- run$loss
      best_idx <- k
    } else if (is.finite(run$loss) && run$loss < best_loss) {
      best <- run
      best_loss <- run$loss
      best_idx <- k
    }
    # Early exit when we effectively hit the optimum.
    if (is.finite(best_loss) && best_loss < 1e-12) break
  }

  structure(
    list(
      rotation = best$rotation,
      loss = best$loss,
      loss_trace = best$loss_trace,
      epsilon_trace = best$epsilon_trace,
      converged = best$converged,
      iterations = best$iterations,
      n1 = n1,
      n2 = n2,
      d = d,
      transport = if (isTRUE(store_transport)) best$transport else NULL,
      n_init = n_init,
      best_init = best_idx
    ),
    class = "ot_procrustes_pair_fit"
  )
}

#' @rdname relative_transform
#' @export
relative_transform.ot_procrustes_pair_fit <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) stop("from and to must differ", call. = FALSE)
  op <- if (from == "i" && to == "j") fit$rotation else t(fit$rotation)
  new_align_transform("O", op, from = from, to = to, k = fit$d, dim = c(fit$n1, fit$n2))
}

#' @rdname pair_loss
#' @export
pair_loss.ot_procrustes_pair_fit <- function(fit, X_i = NULL, X_j = NULL, ...) {
  fit$loss %||% NA_real_
}

#' @rdname latent_dim
#' @export
latent_dim.ot_procrustes_pair_fit <- function(fit, ...) fit$d
