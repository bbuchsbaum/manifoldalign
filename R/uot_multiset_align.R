#' Translation-invariant unbalanced OT (KL) via TI-Sinkhorn
#'
#' @description
#' Solve entropic unbalanced OT with KL marginal penalties using the
#' translation-invariant (TI) Sinkhorn updates from Séjourné, Vialard, Peyré
#' (AISTATS 2022). Supports dense costs and sparse neighbourhood graphs.
#'
#' @param cost A dense numeric matrix (n x m) or a sparse cost list. Sparse
#'   costs must include CSR fields `row_ptr`, `col_idx`, `cost`, `n_rows`,
#'   `n_cols`. If CSC fields `col_ptr`, `row_idx`, `cost_csc` are present they
#'   are used for faster updates.
#' @param alpha Source masses (length n, nonnegative).
#' @param beta Target masses (length m, nonnegative).
#' @param epsilon Entropic regularization parameter (> 0).
#' @param rho1 KL penalty on the first marginal (> 0).
#' @param rho2 KL penalty on the second marginal (> 0).
#' @param max_iter Maximum number of TI-Sinkhorn iterations.
#' @param tol Stopping tolerance on iterate change.
#'
#' @return
#' A list with translation-invariant potentials `fbar`, `gbar`, translated dual
#' potentials `f`, `g`, translation `lambda`, `iterations`, `converged`, and
#' `residual`.
#'
#' @examples
#' set.seed(1)
#' C <- matrix(runif(30), 5, 6)
#' a <- rep(1, 5)
#' b <- rep(1, 6)
#' fit <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.5, rho1 = 1, rho2 = 1,
#'                           max_iter = 200, tol = 1e-8)
#' str(fit)
#'
#' @export
uot_ti_sinkhorn_kl <- function(cost,
                               alpha,
                               beta,
                               epsilon,
                               rho1,
                               rho2,
                               max_iter = 2000,
                               tol = 1e-6) {

  chk::chk_number(epsilon)
  chk::chk_true(epsilon > 0)
  chk::chk_number(rho1)
  chk::chk_true(rho1 > 0)
  chk::chk_number(rho2)
  chk::chk_true(rho2 > 0)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)

  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  chk::chk_true(all(alpha >= 0))
  chk::chk_true(all(beta >= 0))
  if (sum(alpha) <= 0) stop("alpha must have positive total mass", call. = FALSE)
  if (sum(beta) <= 0) stop("beta must have positive total mass", call. = FALSE)

  if (is.matrix(cost)) {
    return(uot_ti_sinkhorn_kl_dense_cpp(
      cost = cost,
      alpha = alpha,
      beta = beta,
      epsilon = epsilon,
      rho1 = rho1,
      rho2 = rho2,
      max_iter = as.integer(max_iter),
      tol = tol
    ))
  }

  if (!is.list(cost)) {
    stop("cost must be a matrix or a sparse cost list", call. = FALSE)
  }
  needed <- c("row_ptr", "col_idx", "cost", "n_rows", "n_cols")
  miss <- setdiff(needed, names(cost))
  if (length(miss)) stop("Missing sparse cost fields: ", paste(miss, collapse = ", "),
                         call. = FALSE)

  has_csc <- all(c("col_ptr", "row_idx", "cost_csc") %in% names(cost))
  if (isTRUE(has_csc)) {
    return(uot_ti_sinkhorn_kl_sparse_csr_csc_cpp(
      row_ptr = cost$row_ptr,
      col_idx = cost$col_idx,
      cost = cost$cost,
      col_ptr = cost$col_ptr,
      row_idx = cost$row_idx,
      cost_csc = cost$cost_csc,
      n_rows = as.integer(cost$n_rows),
      n_cols = as.integer(cost$n_cols),
      alpha = alpha,
      beta = beta,
      epsilon = epsilon,
      rho1 = rho1,
      rho2 = rho2,
      max_iter = as.integer(max_iter),
      tol = tol
    ))
  }

  uot_ti_sinkhorn_kl_sparse_cpp(
    row_ptr = cost$row_ptr,
    col_idx = cost$col_idx,
    cost = cost$cost,
    n_rows = as.integer(cost$n_rows),
    n_cols = as.integer(cost$n_cols),
    alpha = alpha,
    beta = beta,
    epsilon = epsilon,
    rho1 = rho1,
    rho2 = rho2,
    max_iter = as.integer(max_iter),
    tol = tol
  )
}

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) {
    return(m)
  }
  m + log(sum(exp(x - m)))
}

.normalize_from_log <- function(logw, mass = 1) {
  logZ <- .logsumexp(logw)
  if (!is.finite(logZ) && logZ < 0) {
    return(rep(0, length(logw)))
  }
  out <- exp(logw - logZ) * mass
  out[!is.finite(out)] <- 0
  as.numeric(out)
}

.uot_logpi2_meanF <- function(cost, alpha, beta, fbar, gbar, epsilon, F = NULL) {
  if (is.matrix(cost)) {
    return(uot_kl_logpi2_meanF_dense_cpp(
      cost = cost,
      alpha = alpha,
      beta = beta,
      fbar = fbar,
      gbar = gbar,
      epsilon = epsilon,
      F = F
    ))
  }
  uot_kl_logpi2_meanF_sparse_cpp(
    row_ptr = cost$row_ptr,
    col_idx = cost$col_idx,
    cost = cost$cost,
    n_rows = as.integer(cost$n_rows),
    n_cols = as.integer(cost$n_cols),
    alpha = alpha,
    beta = beta,
    fbar = fbar,
    gbar = gbar,
    epsilon = epsilon,
    F = F
  )
}

#' Build dense or sparse costs for multi-set UOT alignment
#'
#' @description
#' Cost builder for the multi-set alignment loop. Supports dense costs for
#' small parcellations and sparse neighborhood costs (kNN, radius, or hybrid)
#' for large voxel/vertex problems.
#'
#' @param X Source coordinates (n x 3) or any numeric matrix.
#' @param Y Template coordinates (m x 3) or any numeric matrix.
#' @param F Optional source features (n x D).
#' @param G Optional template features (m x D).
#' @param lambda_anat Weight for anatomical squared distance.
#' @param lambda_feat Weight for feature squared distance.
#' @param prior_map Optional integer vector (length n) giving a preferred target
#'   index (1..m) for each source node. When provided with `lambda_prior > 0`,
#'   adds a soft identity/anatomical bias term based on squared distances in
#'   template space between `Y[j,]` and `Y[prior_map[i],]`.
#' @param lambda_prior Nonnegative weight for the soft prior term.
#' @param prior_sigma Optional positive scale (in the units of `Y`). When set,
#'   the prior penalty uses `||Y[j]-Y[pref]||^2 / prior_sigma^2`.
#' @param neighbor_mode Neighborhood mode: `"auto"`, `"dense"`, `"knn"`,
#'   `"radius"`, or `"hybrid"`.
#' @param k_neighbors Integer k for kNN neighborhoods (used by `"knn"` and
#'   `"hybrid"`).
#' @param radius Positive scalar radius for `"radius"` / `"hybrid"` modes.
#' @param maxk Maximum number of candidate neighbors requested in `"radius"` /
#'   `"hybrid"` mode (prevents quadratic blowups).
#' @param min_neighbors Minimum number of within-radius neighbors required
#'   per row in `"hybrid"` mode before falling back to kNN.
#' @param group_X Optional vector (length n) of hard constraint group labels for
#'   source nodes. When provided with `group_Y`, neighborhoods are computed
#'   within matching groups only.
#' @param group_Y Optional vector (length m) of hard constraint group labels for
#'   template nodes. Must share the same label alphabet as `group_X`.
#' @param dense_max_bytes Maximum dense cost size (in bytes) allowed when
#'   `neighbor_mode="auto"`. Above this threshold, `"auto"` selects a sparse
#'   neighborhood.
#' @param ensure_cols Logical; if TRUE (default for sparse kNN/hybrid modes),
#'   adds a minimal set of reverse 1NN edges to ensure every template column has
#'   at least one incoming edge in the sparse graph.
#'
#' @return A dense numeric matrix (n x m) when `neighbor_mode="dense"` (or when
#'   `neighbor_mode="auto"` selects dense); otherwise a sparse cost list with
#'   CSR+CSC fields.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(30), 10, 3)
#' Y <- matrix(rnorm(36), 12, 3)
#' C <- uot_build_cost(X, Y, neighbor_mode = "dense")
#' str(C)
#'
#' @export
uot_build_cost <- function(X, Y,
                           F = NULL,
                           G = NULL,
                           lambda_anat = 1,
                           lambda_feat = 0,
                           prior_map = NULL,
                           lambda_prior = 0,
                           prior_sigma = NULL,
                           neighbor_mode = c("auto", "dense", "knn", "radius", "hybrid"),
                           k_neighbors = NULL,
                           radius = NULL,
                           maxk = 128L,
                           min_neighbors = 1L,
                           group_X = NULL,
                           group_Y = NULL,
                           dense_max_bytes = 256e6,
                           ensure_cols = NULL) {

  X <- as.matrix(X)
  Y <- as.matrix(Y)
  chk::chk_true(nrow(X) > 0)
  chk::chk_true(nrow(Y) > 0)
  if (!all(is.finite(X)) || !all(is.finite(Y))) {
    stop("X and Y must contain only finite values", call. = FALSE)
  }
  chk::chk_number(lambda_anat)
  chk::chk_true(lambda_anat >= 0)
  chk::chk_number(lambda_feat)
  chk::chk_true(lambda_feat >= 0)
  chk::chk_number(lambda_prior)
  chk::chk_true(lambda_prior >= 0)
  if (!is.null(prior_sigma)) {
    chk::chk_number(prior_sigma)
    chk::chk_true(is.finite(prior_sigma))
    chk::chk_true(prior_sigma > 0)
  }
  if (!is.null(prior_map)) {
    prior_map <- as.integer(prior_map)
    if (length(prior_map) != nrow(X)) {
      stop("prior_map length must match nrow(X)", call. = FALSE)
    }
    bad <- !is.na(prior_map) & !(prior_map > 0 & prior_map <= nrow(Y))
    if (any(bad)) {
      stop("prior_map must contain only NA or indices in 1..nrow(Y)", call. = FALSE)
    }
  } else if (lambda_prior > 0) {
    stop("lambda_prior > 0 requires prior_map", call. = FALSE)
  }

  neighbor_mode <- match.arg(neighbor_mode)

  if (!is.null(group_X) || !is.null(group_Y)) {
    if (is.null(group_X) || is.null(group_Y)) {
      stop("group_X and group_Y must be provided together", call. = FALSE)
    }
    if (length(group_X) != nrow(X)) stop("group_X length must match nrow(X)",
                                         call. = FALSE)
    if (length(group_Y) != nrow(Y)) stop("group_Y length must match nrow(Y)",
                                         call. = FALSE)
  }

  if (neighbor_mode == "auto") {
    dense_bytes <- as.double(nrow(X)) * as.double(nrow(Y)) * 8
    if (!is.null(radius) && is.finite(radius)) {
      neighbor_mode <- if (is.null(k_neighbors)) "radius" else "hybrid"
    } else if (!is.null(k_neighbors)) {
      neighbor_mode <- "knn"
    } else if (dense_bytes <= dense_max_bytes) {
      neighbor_mode <- "dense"
    } else {
      neighbor_mode <- "knn"
      k_neighbors <- 64L
    }
  }

  if (neighbor_mode == "dense") {
    if (!is.null(ensure_cols) && isTRUE(ensure_cols)) {
      warning("ensure_cols ignored for dense costs", call. = FALSE)
    }
    C <- matrix(0, nrow(X), nrow(Y))
    if (lambda_anat > 0) {
      C <- C + lambda_anat * compute_squared_distances_cpp_fallback(X, Y)
    }
    if (lambda_feat > 0 && !is.null(F) && !is.null(G)) {
      F <- as.matrix(F)
      G <- as.matrix(G)
      if (!all(is.finite(F)) || !all(is.finite(G))) {
        stop("F and G must contain only finite values", call. = FALSE)
      }
      C <- C + lambda_feat * compute_squared_distances_cpp_fallback(F, G)
    }
    if (lambda_prior > 0) {
      pref <- prior_map
      valid <- !is.na(pref) & pref > 0L & pref <= nrow(Y)
      if (any(valid)) {
        pref2 <- pref
        pref2[!valid] <- 1L
        Ypref <- Y[pref2, , drop = FALSE]
        P <- compute_squared_distances_cpp_fallback(Ypref, Y)
        if (!is.null(prior_sigma)) {
          P <- P / (as.numeric(prior_sigma)^2)
        }
        P[!valid, ] <- 0
        C <- C + lambda_prior * P
      }
    }
    if (!is.null(group_X)) {
      gx <- as.character(group_X)
      gy <- as.character(group_Y)
      ok <- outer(gx, gy, FUN = `==`)
      if (any(!ok)) {
        C[!ok] <- Inf
      }
      if (any(rowSums(is.finite(C)) == 0)) {
        stop("Dense constraints left at least one row with no finite costs",
             call. = FALSE)
      }
    }
    return(C)
  }

  if (is.null(ensure_cols)) {
    ensure_cols <- neighbor_mode %in% c("knn", "hybrid")
  }

  k_neighbors <- if (!is.null(k_neighbors)) as.integer(k_neighbors) else NULL
  maxk <- as.integer(maxk)
  min_neighbors <- as.integer(min_neighbors)
  chk::chk_true(is.null(k_neighbors) || k_neighbors > 0)
  chk::chk_true(maxk > 0)
  chk::chk_true(min_neighbors >= 0)

  if (neighbor_mode %in% c("radius", "hybrid")) {
    chk::chk_number(radius)
    chk::chk_true(radius > 0)
    if (!requireNamespace("Rnanoflann", quietly = TRUE)) {
      stop("Rnanoflann is required for radius neighborhoods", call. = FALSE)
    }
  }

  # -- Neighbour search helpers --------------------------------------------
  nn_knn <- function(Xq, Yd, k) {
    safe_compute(
      RANN::nn2(data = Yd, query = Xq, k = min(k, nrow(Yd))),
      "kNN neighbour search failed."
    )
  }

  nn_radius <- function(Xq, Yd, k, r) {
    safe_compute(
      Rnanoflann::nn(data = Yd, points = Xq, k = min(k, nrow(Yd)), radius = r),
      "Radius neighbour search failed."
    )
  }

  run_neighbors <- function(Xq, Yd, mode) {
    if (mode == "knn") {
      if (is.null(k_neighbors)) stop("k_neighbors must be set for knn mode",
                                     call. = FALSE)
      out <- nn_knn(Xq, Yd, k_neighbors)
      list(idx = out$nn.idx, dists = out$nn.dists, max_dist = Inf)
    } else if (mode == "radius") {
      out <- nn_radius(Xq, Yd, maxk, radius)
      idx <- out$indices
      d <- out$distances
      bad <- !is.finite(d) | d > radius
      idx[bad] <- NA_integer_
      d[bad] <- NA_real_
      list(idx = idx, dists = d, max_dist = radius)
    } else { # hybrid
      if (is.null(k_neighbors)) stop("k_neighbors must be set for hybrid mode",
                                     call. = FALSE)
      k_r <- max(maxk, k_neighbors)
      out_r <- nn_radius(Xq, Yd, k_r, radius)
      idx <- out_r$indices
      d <- out_r$distances
      bad <- !is.finite(d) | d > radius
      idx[bad] <- NA_integer_
      d[bad] <- NA_real_

      n_ok <- rowSums(!is.na(idx) & idx > 0)
      need <- which(n_ok < min_neighbors)
      if (length(need)) {
        out_k <- nn_knn(Xq[need, , drop = FALSE], Yd, k_neighbors)
        idx[need, seq_len(ncol(out_k$nn.idx))] <- out_k$nn.idx
        d[need, seq_len(ncol(out_k$nn.dists))] <- out_k$nn.dists
      }
      list(idx = idx, dists = d, max_dist = Inf)
    }
  }

  build_sparse <- function(Xq, Yd, group_x = NULL, group_y = NULL) {
    check_isolates <- function(out) {
      if (neighbor_mode != "radius") return(out)
      nnz_row <- diff(out$row_ptr)
      n0 <- sum(nnz_row == 0L)
      if (n0 > 0L) {
        stop(
          "Radius neighborhood produced ", n0, " source rows with zero edges. ",
          "Increase `radius`, increase `maxk`, or use `neighbor_mode='hybrid'`.",
          call. = FALSE
        )
      }
      out
    }

    if (is.null(group_x)) {
      nn <- run_neighbors(Xq, Yd, neighbor_mode)
      idx <- nn$idx
      d <- nn$dists
      max_dist <- nn$max_dist

      # Optional: ensure every column has an incoming edge (kNN/hybrid only).
      extra_i <- integer(0)
      extra_j <- integer(0)
      extra_d <- numeric(0)
      if (isTRUE(ensure_cols) && neighbor_mode %in% c("knn", "hybrid")) {
        covered <- tabulate(as.integer(idx[!is.na(idx) & idx > 0]), nbins = nrow(Yd))
        miss <- which(covered == 0)
        if (length(miss)) {
          out_rev <- nn_knn(Yd[miss, , drop = FALSE], Xq, k = 1L)
          extra_i <- as.integer(out_rev$nn.idx[, 1])
          extra_j <- as.integer(miss)
          extra_d <- as.numeric(out_rev$nn.dists[, 1])
        }
      }

      out <- uot_build_sparse_cost_knn_cpp(
        nn_idx = idx,
        nn_dists = d,
        n_cols = nrow(Yd),
        lambda_anat = lambda_anat,
        F = if (!is.null(F)) as.matrix(F) else NULL,
        G = if (!is.null(G)) as.matrix(G) else NULL,
        lambda_feat = lambda_feat,
        max_dist = max_dist,
        extra_i = extra_i,
        extra_j = extra_j,
        extra_dists = extra_d,
        Y_template = if (lambda_prior > 0) as.matrix(Yd) else NULL,
        prior_map = if (lambda_prior > 0) prior_map else NULL,
        lambda_prior = lambda_prior,
        prior_sigma = if (is.null(prior_sigma)) NA_real_ else as.numeric(prior_sigma)
      )
      return(check_isolates(out))
    }

    # Grouped neighborhoods: compute within matching groups only.
    gx <- as.character(group_x)
    gy <- as.character(group_y)
    groups <- sort(unique(gx))

    k_cols <- if (neighbor_mode == "knn") k_neighbors else maxk
    if (neighbor_mode == "hybrid") k_cols <- max(maxk, k_neighbors)

    idx_all <- matrix(NA_integer_, nrow(Xq), k_cols)
    d_all <- matrix(NA_real_, nrow(Xq), k_cols)
    extra_i <- integer(0)
    extra_j <- integer(0)
    extra_d <- numeric(0)

    for (g in groups) {
      rows <- which(gx == g)
      cols <- which(gy == g)
      if (!length(rows)) next
      if (!length(cols)) {
        stop("No template nodes for group '", g, "'", call. = FALSE)
      }

      nn <- run_neighbors(Xq[rows, , drop = FALSE], Yd[cols, , drop = FALSE],
                          neighbor_mode)
      idx_sub <- nn$idx
      d_sub <- nn$dists
      max_dist <- nn$max_dist

      # Map subset indices back to global template indices.
      idx_sub[!is.na(idx_sub) & idx_sub > 0] <- cols[idx_sub[!is.na(idx_sub) & idx_sub > 0]]

      idx_all[rows, seq_len(ncol(idx_sub))] <- idx_sub
      d_all[rows, seq_len(ncol(d_sub))] <- d_sub

      if (isTRUE(ensure_cols) && neighbor_mode %in% c("knn", "hybrid")) {
        covered <- tabulate(as.integer(idx_sub[!is.na(idx_sub) & idx_sub > 0]),
                            nbins = nrow(Yd))
        miss <- cols[covered[cols] == 0]
        if (length(miss)) {
          out_rev <- nn_knn(Yd[miss, , drop = FALSE], Xq[rows, , drop = FALSE],
                            k = 1L)
          extra_i <- c(extra_i, as.integer(rows[out_rev$nn.idx[, 1]]))
          extra_j <- c(extra_j, as.integer(miss))
          extra_d <- c(extra_d, as.numeric(out_rev$nn.dists[, 1]))
        }
      }
    }

    out <- uot_build_sparse_cost_knn_cpp(
      nn_idx = idx_all,
      nn_dists = d_all,
      n_cols = nrow(Yd),
      lambda_anat = lambda_anat,
      F = if (!is.null(F)) as.matrix(F) else NULL,
      G = if (!is.null(G)) as.matrix(G) else NULL,
      lambda_feat = lambda_feat,
      max_dist = if (neighbor_mode == "radius") radius else Inf,
      extra_i = extra_i,
      extra_j = extra_j,
      extra_dists = extra_d,
      Y_template = if (lambda_prior > 0) as.matrix(Yd) else NULL,
      prior_map = if (lambda_prior > 0) prior_map else NULL,
      lambda_prior = lambda_prior,
      prior_sigma = if (is.null(prior_sigma)) NA_real_ else as.numeric(prior_sigma)
    )
    check_isolates(out)
  }

  build_sparse(X, Y, group_x = group_X, group_y = group_Y)
}

#' Apply a UOT coupling to map signals into template space
#'
#' @description
#' Applies the entropic UOT coupling implied by translation-invariant dual
#' potentials, without materializing the full transport plan.
#'
#' @param cost A dense numeric matrix (n x m) or a sparse cost list with CSR
#'   fields from \code{\link{uot_build_cost}}.
#' @param alpha Source masses (length n).
#' @param beta Target masses (length m).
#' @param fbar Translation-invariant source potential (length n).
#' @param gbar Translation-invariant target potential (length m).
#' @param epsilon Entropic regularization parameter (> 0).
#' @param signal Source signal: numeric vector (length n) or matrix (T x n).
#' @param delta Stabilizer added to the denominator.
#'
#' @return A numeric vector (length m) or matrix (T x m).
#'
#' @examples
#' set.seed(1)
#' C <- matrix(runif(30), 5, 6)
#' a <- rep(1, 5)
#' b <- rep(1, 6)
#' fit <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.5, rho1 = 1, rho2 = 1,
#'                           max_iter = 200, tol = 1e-8)
#' s <- rnorm(5)
#' y <- uot_apply_map(C, a, b, fit$fbar, fit$gbar, epsilon = 0.5, signal = s)
#' str(y)
#'
#' @export
uot_apply_map <- function(cost,
                          alpha,
                          beta,
                          fbar,
                          gbar,
                          epsilon,
                          signal,
                          delta = 1e-8) {

  chk::chk_number(epsilon)
  chk::chk_true(epsilon > 0)
  chk::chk_number(delta)
  chk::chk_true(delta >= 0)

  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  fbar <- as.numeric(fbar)
  gbar <- as.numeric(gbar)

  if (is.matrix(cost)) {
    if (is.matrix(signal)) {
      return(uot_apply_map_dense_mat_cpp(
        cost = cost,
        alpha = alpha,
        beta = beta,
        fbar = fbar,
        gbar = gbar,
        epsilon = epsilon,
        signal = as.matrix(signal),
        delta = delta
      ))
    }
    return(as.numeric(uot_apply_map_dense_vec_cpp(
      cost = cost,
      alpha = alpha,
      beta = beta,
      fbar = fbar,
      gbar = gbar,
      epsilon = epsilon,
      signal = as.numeric(signal),
      delta = delta
    )))
  }

  if (!is.list(cost)) stop("cost must be a matrix or a sparse cost list",
                           call. = FALSE)
  if (is.matrix(signal)) {
    has_csc <- all(c("col_ptr", "row_idx", "cost_csc") %in% names(cost))
    if (!isTRUE(has_csc)) {
      stop(
        "Sparse mapping for matrix signals requires CSC fields ",
        "(`col_ptr`, `row_idx`, `cost_csc`). Rebuild the cost with `uot_build_cost()`.",
        call. = FALSE
      )
    }
    return(uot_apply_map_sparse_mat_cpp(
      col_ptr = cost$col_ptr,
      row_idx = cost$row_idx,
      cost_csc = cost$cost_csc,
      n_rows = as.integer(cost$n_rows),
      n_cols = as.integer(cost$n_cols),
      alpha = alpha,
      beta = beta,
      fbar = fbar,
      gbar = gbar,
      epsilon = epsilon,
      signal = as.matrix(signal),
      delta = delta
    ))
  }

  as.numeric(uot_apply_map_sparse_vec_cpp(
    row_ptr = cost$row_ptr,
    col_idx = cost$col_idx,
    cost = cost$cost,
    n_rows = as.integer(cost$n_rows),
    n_cols = as.integer(cost$n_cols),
    alpha = alpha,
    beta = beta,
    fbar = fbar,
    gbar = gbar,
    epsilon = epsilon,
    signal = as.numeric(signal),
    delta = delta
  ))
}

#' Extract a sparse coupling or mapping operator (dgCMatrix)
#'
#' @description
#' Given translation-invariant dual potentials from [uot_ti_sinkhorn_kl()], this
#' function computes edge weights on a sparse cost graph and returns a
#' `Matrix::dgCMatrix`. This is useful for downstream pipelines that want to:
#' (1) prune the coupling and/or (2) reuse the same operator across multiple
#' signals without recomputing exponentials.
#'
#' @param cost A sparse cost list as returned by [uot_build_cost()] containing
#'   CSC fields `col_ptr`, `row_idx`, `cost_csc`, plus dimensions `n_rows`,
#'   `n_cols`.
#' @param alpha Source masses (length `n_rows`, nonnegative, not all zero).
#' @param beta Target masses (length `n_cols`, nonnegative, not all zero).
#' @param fbar Translation-invariant source potential (length `n_rows`).
#' @param gbar Translation-invariant target potential (length `n_cols`).
#' @param epsilon Entropic regularization parameter (> 0).
#' @param weights Which weights to compute:
#'   - `"map"`: barycentric map weights \eqn{w_{ij} / (\pi_{2,j} + \delta)}
#'   - `"cond"`: conditional weights \eqn{w_{ij} / \pi_{2,j}} (column-stochastic)
#'   - `"coupling"`: raw coupling weights \eqn{w_{ij}} (may overflow for extreme
#'     parameter regimes).
#' @param delta Stabilizer for `"map"` weights.
#' @param prune_topk_col Optional integer top-k pruning per template column.
#'   `NULL` keeps all edges.
#' @param prune_threshold Optional nonnegative threshold pruning applied after
#'   the weight transform. `NULL` keeps all edges.
#'
#' @return A list with:
#' - `coupling`: a `dgCMatrix` of dimension `n_rows x n_cols`
#' - `log_pi2`: numeric vector (length `n_cols`) with `log(pi2)` on the sparse graph
#'
#' @export
uot_extract_coupling <- function(cost,
                                 alpha,
                                 beta,
                                 fbar,
                                 gbar,
                                 epsilon,
                                 weights = c("map", "cond", "coupling"),
                                 delta = 1e-8,
                                 prune_topk_col = NULL,
                                 prune_threshold = NULL) {

  chk::chk_number(epsilon)
  chk::chk_true(epsilon > 0)
  chk::chk_number(delta)
  chk::chk_true(delta >= 0)

  if (!is.list(cost)) stop("cost must be a sparse cost list", call. = FALSE)
  needed <- c("col_ptr", "row_idx", "cost_csc", "n_rows", "n_cols")
  miss <- setdiff(needed, names(cost))
  if (length(miss)) {
    stop("Missing sparse cost fields: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)
  fbar <- as.numeric(fbar)
  gbar <- as.numeric(gbar)
  chk::chk_true(all(alpha >= 0))
  chk::chk_true(all(beta >= 0))
  if (sum(alpha) <= 0) stop("alpha must have positive total mass", call. = FALSE)
  if (sum(beta) <= 0) stop("beta must have positive total mass", call. = FALSE)

  weights <- match.arg(weights)
  weight_type <- switch(weights,
    map = 0L,
    cond = 1L,
    coupling = 2L
  )

  topk_col <- if (is.null(prune_topk_col)) 0L else as.integer(prune_topk_col)
  thresh <- if (is.null(prune_threshold)) 0 else as.numeric(prune_threshold)
  chk::chk_true(topk_col >= 0)
  chk::chk_true(thresh >= 0)

  out <- uot_extract_coupling_sparse_csc_cpp(
    col_ptr = cost$col_ptr,
    row_idx = cost$row_idx,
    cost_csc = cost$cost_csc,
    n_rows = as.integer(cost$n_rows),
    n_cols = as.integer(cost$n_cols),
    alpha = alpha,
    beta = beta,
    fbar = fbar,
    gbar = gbar,
    epsilon = epsilon,
    weight_type = as.integer(weight_type),
    delta = delta,
    topk_col = as.integer(topk_col),
    threshold = thresh
  )

  out$weights <- weights
  out
}

#' Fit a pairwise UOT alignment to a fixed template support
#'
#' @description
#' Convenience wrapper that builds a dense/sparse cost graph via [uot_build_cost()]
#' and then solves entropic KL-UOT via [uot_ti_sinkhorn_kl()].
#'
#' This is the recommended entry point for **pairwise** workflows; it bundles the
#' pieces you need for mapping and coupling extraction into a single object.
#'
#' @param X Source coordinates (n x 3) or any numeric matrix.
#' @param Y Template coordinates (m x 3) or any numeric matrix.
#' @param alpha Source masses (length n, nonnegative, not all zero).
#' @param beta Target masses (length m, nonnegative, not all zero).
#' @param F Optional source features (n x D).
#' @param G Optional template features (m x D).
#' @param lambda_anat Weight for anatomical squared distance.
#' @param lambda_feat Weight for feature squared distance.
#' @param prior_map Optional integer vector (length n) giving a preferred target
#'   index (1..m) for each source node. Passed to [uot_build_cost()].
#' @param lambda_prior Nonnegative weight for the soft prior term.
#' @param prior_sigma Optional positive scale for the prior term.
#' @param neighbor_mode Neighborhood mode passed to [uot_build_cost()].
#' @param k_neighbors Integer k for kNN neighborhoods (used by `"knn"` and
#'   `"hybrid"`).
#' @param radius Positive scalar radius for `"radius"` / `"hybrid"` modes.
#' @param maxk Maximum number of candidate neighbors requested in `"radius"` /
#'   `"hybrid"` mode (prevents quadratic blowups).
#' @param min_neighbors Minimum number of within-radius neighbors required per
#'   row in `"hybrid"` mode before falling back to kNN.
#' @param group_X Optional vector (length n) of hard constraint group labels for
#'   source nodes. When provided with `group_Y`, neighborhoods are computed
#'   within matching groups only.
#' @param group_Y Optional vector (length m) of hard constraint group labels for
#'   template nodes.
#' @param dense_max_bytes Maximum dense cost size (in bytes) allowed when
#'   `neighbor_mode="auto"`.
#' @param ensure_cols Logical; if TRUE, adds reverse 1NN edges to ensure every
#'   template column has at least one incoming edge (kNN/hybrid modes only).
#' @param epsilon Entropic regularization parameter (> 0).
#' @param rho1 KL penalty on the first marginal (> 0).
#' @param rho2 KL penalty on the second marginal (> 0).
#' @param max_iter Maximum number of TI-Sinkhorn iterations.
#' @param tol TI-Sinkhorn tolerance.
#'
#' @return An object of class `uot_pair_fit` with fields `cost`, `sol`, `alpha`,
#'   `beta`, `epsilon`, `rho1`, `rho2`.
#'
#' @export
uot_fit_pair <- function(X, Y,
                         alpha,
                         beta,
                         F = NULL,
                         G = NULL,
                         lambda_anat = 1,
                         lambda_feat = 0,
                         prior_map = NULL,
                         lambda_prior = 0,
                         prior_sigma = NULL,
                         neighbor_mode = c("auto", "dense", "knn", "radius", "hybrid"),
                         k_neighbors = NULL,
                         radius = NULL,
                         maxk = 128L,
                         min_neighbors = 1L,
                         group_X = NULL,
                         group_Y = NULL,
                         dense_max_bytes = 256e6,
                         ensure_cols = NULL,
                         epsilon,
                         rho1,
                         rho2,
                         max_iter = 2000,
                         tol = 1e-6) {

  neighbor_mode <- match.arg(neighbor_mode)
  alpha <- as.numeric(alpha)
  beta <- as.numeric(beta)

  cost <- uot_build_cost(
    X = X, Y = Y,
    F = F, G = G,
    lambda_anat = lambda_anat,
    lambda_feat = lambda_feat,
    prior_map = prior_map,
    lambda_prior = lambda_prior,
    prior_sigma = prior_sigma,
    neighbor_mode = neighbor_mode,
    k_neighbors = k_neighbors,
    radius = radius,
    maxk = maxk,
    min_neighbors = min_neighbors,
    group_X = group_X,
    group_Y = group_Y,
    dense_max_bytes = dense_max_bytes,
    ensure_cols = ensure_cols
  )

  sol <- uot_ti_sinkhorn_kl(
    cost = cost,
    alpha = alpha,
    beta = beta,
    epsilon = epsilon,
    rho1 = rho1,
    rho2 = rho2,
    max_iter = max_iter,
    tol = tol
  )

  structure(
    list(
      cost = cost,
      sol = sol,
      alpha = alpha,
      beta = beta,
      epsilon = epsilon,
      rho1 = rho1,
      rho2 = rho2
    ),
    class = "uot_pair_fit"
  )
}

#' Map signals into template space using a `uot_pair_fit`
#'
#' @param fit An object returned by [uot_fit_pair()].
#' @param signal Numeric vector (length n) or matrix (T x n).
#' @param delta Stabilizer for the barycentric denominator.
#'
#' @return A numeric vector (length m) or matrix (T x m).
#' @export
uot_map <- function(fit, signal, delta = 1e-8) {
  if (!inherits(fit, "uot_pair_fit")) {
    stop("fit must be a uot_pair_fit object (from uot_fit_pair())", call. = FALSE)
  }
  uot_apply_map(
    cost = fit$cost,
    alpha = fit$alpha,
    beta = fit$beta,
    fbar = fit$sol$fbar,
    gbar = fit$sol$gbar,
    epsilon = fit$epsilon,
    signal = signal,
    delta = delta
  )
}

#' Extract a reusable map/coupling operator from a `uot_pair_fit`
#'
#' @param fit An object returned by [uot_fit_pair()].
#' @param weights Passed to [uot_extract_coupling()].
#' @param delta Stabilizer for `"map"` weights.
#' @param prune_topk_col Optional top-k pruning per template column.
#' @param prune_threshold Optional threshold pruning.
#'
#' @return A `dgCMatrix` of dimension `n x m`.
#' @export
uot_operator <- function(fit,
                         weights = c("map", "cond", "coupling"),
                         delta = 1e-8,
                         prune_topk_col = NULL,
                         prune_threshold = NULL) {
  if (!inherits(fit, "uot_pair_fit")) {
    stop("fit must be a uot_pair_fit object (from uot_fit_pair())", call. = FALSE)
  }
  weights <- match.arg(weights)

  if (is.matrix(fit$cost)) {
    if (identical(weights, "coupling")) {
      stop("weights='coupling' is not supported for dense costs; use uot_apply_map() or sparse costs",
           call. = FALSE)
    }

    C <- fit$cost
    n <- nrow(C)
    m <- ncol(C)
    alpha <- as.numeric(fit$alpha)
    beta <- as.numeric(fit$beta)
    fbar <- as.numeric(fit$sol$fbar)
    gbar <- as.numeric(fit$sol$gbar)
    eps <- fit$epsilon

    loga <- ifelse(alpha > 0, log(alpha), -Inf)
    logb <- ifelse(beta > 0, log(beta), -Inf)

    Q <- matrix(0, n, m)
    inv_eps <- 1 / eps
    for (j in seq_len(m)) {
      if (!is.finite(logb[j])) next
      logw <- loga + logb[j] + (fbar + gbar[j] - C[, j]) * inv_eps
      mx <- max(logw)
      if (!is.finite(mx)) next
      e <- exp(logw - mx)
      e[!is.finite(e)] <- 0
      s <- sum(e)
      if (!(s > 0)) next
      denom <- if (identical(weights, "map")) {
        s + delta * exp(-mx)
      } else {
        s
      }
      if (!(denom > 0)) next
      Q[, j] <- e / denom
    }

    if (is.null(prune_topk_col) && is.null(prune_threshold)) {
      return(Q)
    }

    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Matrix package is required for pruning operators", call. = FALSE)
    }
    topk <- if (is.null(prune_topk_col)) Inf else as.integer(prune_topk_col)
    thresh <- if (is.null(prune_threshold)) 0 else as.numeric(prune_threshold)

    ii <- integer(0)
    jj <- integer(0)
    xx <- numeric(0)
    for (j in seq_len(m)) {
      v <- Q[, j]
      keep <- which(v >= thresh & v > 0)
      if (!length(keep)) next
      if (is.finite(topk) && length(keep) > topk) {
        ord <- order(v[keep], decreasing = TRUE)[seq_len(topk)]
        keep <- keep[ord]
      }
      ii <- c(ii, keep)
      jj <- c(jj, rep.int(j, length(keep)))
      xx <- c(xx, v[keep])
    }

    return(Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, m)))
  }

  out <- uot_extract_coupling(
    cost = fit$cost,
    alpha = fit$alpha,
    beta = fit$beta,
    fbar = fit$sol$fbar,
    gbar = fit$sol$gbar,
    epsilon = fit$epsilon,
    weights = weights,
    delta = delta,
    prune_topk_col = prune_topk_col,
    prune_threshold = prune_threshold
  )
  out$coupling
}

#' Map multiple subjects into template space from a multi-set fit
#'
#' @description
#' Convenience helper for [multiset_uot_align()]. It applies each subject's
#' fitted UOT map to the corresponding signal, using the **final** template
#' weights and the pairwise potentials stored in `fit$fits`.
#'
#' @param fit A result from [multiset_uot_align()].
#' @param datasets The datasets list originally passed to [multiset_uot_align()].
#'   Each dataset must contain `alpha`.
#' @param signals A list of length K of signals, each a numeric vector (length N_k)
#'   or matrix (T x N_k).
#' @param epsilon Optional override of the entropic parameter. If NULL, uses
#'   `fit$params$epsilon`.
#' @param delta Stabilizer for the barycentric denominator.
#'
#' @return A list of mapped signals (length K), each in template space.
#' @export
multiset_uot_map <- function(fit,
                             datasets,
                             signals,
                             epsilon = NULL,
                             delta = 1e-8) {

  chk::chk_true(is.list(fit) && !is.null(fit$fits) && !is.null(fit$template))
  chk::chk_true(is.list(datasets) && length(datasets) == length(fit$fits))
  chk::chk_true(is.list(signals) && length(signals) == length(fit$fits))

  if (is.null(epsilon)) {
    epsilon <- fit$params$epsilon
  }
  chk::chk_number(epsilon)
  chk::chk_true(epsilon > 0)
  chk::chk_number(delta)
  chk::chk_true(delta >= 0)

  beta <- as.numeric(fit$template$beta)

  lapply(seq_along(fit$fits), function(k) {
    fk <- fit$fits[[k]]
    alpha <- as.numeric(datasets[[k]]$alpha)
    uot_apply_map(
      cost = fk$cost,
      alpha = alpha,
      beta = beta,
      fbar = fk$sol$fbar,
      gbar = fk$sol$gbar,
      epsilon = epsilon,
      signal = signals[[k]],
      delta = delta
    )
  })
}

#' Multi-subject UOT alignment to a shared template
#'
#' @description
#' Iteratively aligns multiple subject measures to a fixed-support template by
#' solving subject-to-template UOT problems and updating template weights (and
#' optionally template features) from transported marginals.
#'
#' @param datasets A list of datasets. Each element must be a list with fields
#'   `X` (n x 3 coordinates), `alpha` (length n masses), and optional `F`
#'   (n x D features).
#' @param template A list with fields `Y` (m x 3 coordinates), `beta` (length m
#'   masses), and optional `G` (m x D features).
#' @param omega Optional nonnegative weights over datasets (length K). If NULL,
#'   uses uniform weights.
#' @param epsilon Entropic regularization parameter (> 0).
#' @param rho1 KL penalty on the first marginal (> 0).
#' @param rho2 KL penalty on the second marginal (> 0).
#' @param lambda_anat Weight for anatomical squared distance.
#' @param lambda_feat Weight for feature squared distance.
#' @param neighbor_mode Neighborhood mode passed to [uot_build_cost()].
#' @param k_neighbors Integer k for kNN neighborhoods (used by `"knn"` and
#'   `"hybrid"` modes).
#' @param radius Positive scalar radius for `"radius"` / `"hybrid"` modes.
#' @param maxk Maximum number of candidate neighbors requested in `"radius"` /
#'   `"hybrid"` mode (prevents quadratic blowups).
#' @param min_neighbors Minimum number of within-radius neighbors required per
#'   row in `"hybrid"` mode before falling back to kNN.
#' @param constraint_fields Optional character vector of field names present on
#'   each dataset and the template that define **hard transport constraints**.
#'   Edges are only allowed when all fields match (implemented via grouped
#'   neighbor search). Typical examples: `c("hemi")` or `c("hemi","network")`.
#' @param prior_fields Optional character vector of field names present on each
#'   dataset and the template that define a **soft identity prior**. When
#'   provided together with `lambda_prior > 0`, each dataset builds a
#'   `prior_map` by matching these fields to the template and adds a soft bias
#'   term via [uot_build_cost()].
#' @param lambda_prior Nonnegative weight for the soft prior term.
#' @param prior_sigma Optional positive scale for the prior term (in the units
#'   of `template$Y`).
#' @param dense_max_bytes Maximum dense cost size (in bytes) allowed when
#'   `neighbor_mode="auto"`.
#' @param ensure_cols Logical; if TRUE, adds reverse 1NN edges to ensure every
#'   template column has at least one incoming edge (kNN/hybrid modes only).
#' @param max_outer Maximum number of template update iterations.
#' @param max_inner Maximum number of TI-Sinkhorn iterations per subject.
#' @param tol_inner TI-Sinkhorn tolerance.
#' @param tol_outer Template weight tolerance.
#' @param rescale_pi2 How to scale each subject's transported mass before
#'   aggregating (`"alpha"` or `"unit"`).
#' @param target_mass Target total mass for template weights. If `"mean_alpha"`,
#'   uses mean subject mass; if `"none"`, no rescaling; if numeric, uses that.
#' @param learn_template_features Logical; update template features `G`.
#' @param delta Small stabilizer for feature update denominators.
#' @param parallel Logical; if TRUE, uses `parallel::mclapply` on Unix.
#' @param ncores Number of cores for parallel subject solves.
#' @param verbose Logical; print outer-loop progress.
#'
#' @return A list with updated `template` and per-subject `fits`.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' K <- 3
#' datasets <- lapply(seq_len(K), function(k) {
#'   X <- matrix(rnorm(30), 10, 3)
#'   F <- matrix(rnorm(40), 10, 4)
#'   list(X = X, alpha = rep(1, 10), F = F)
#' })
#' template <- list(Y = matrix(rnorm(36), 12, 3),
#'                  beta = rep(1, 12),
#'                  G = matrix(0, 12, 4))
#' fit <- multiset_uot_align(datasets, template,
#'   epsilon = 0.5, rho1 = 1, rho2 = 1,
#'   lambda_anat = 1, lambda_feat = 0.1,
#'   k_neighbors = 6, max_outer = 3, max_inner = 200
#' )
#' str(fit$template$beta)
#' }
#'
#' @export
multiset_uot_align <- function(datasets,
                               template,
                               omega = NULL,
                               epsilon,
                               rho1,
                               rho2,
                               lambda_anat = 1,
                               lambda_feat = 0,
                               neighbor_mode = c("auto", "dense", "knn", "radius", "hybrid"),
                               k_neighbors = NULL,
                               radius = NULL,
                               maxk = 128L,
                               min_neighbors = 1L,
                               constraint_fields = NULL,
                               prior_fields = NULL,
                               lambda_prior = 0,
                               prior_sigma = NULL,
                               dense_max_bytes = 256e6,
                               ensure_cols = NULL,
                               max_outer = 10,
                               max_inner = 2000,
                               tol_inner = 1e-6,
                               tol_outer = 1e-4,
                               rescale_pi2 = c("alpha", "unit"),
                               target_mass = "mean_alpha",
                               learn_template_features = FALSE,
                               delta = 1e-8,
                               parallel = FALSE,
                               ncores = 1,
                               verbose = FALSE) {

  rescale_pi2 <- match.arg(rescale_pi2)
  neighbor_mode <- match.arg(neighbor_mode)
  if (is.character(target_mass)) {
    target_mass <- match.arg(target_mass, c("mean_alpha", "none"))
  } else {
    chk::chk_number(target_mass)
    chk::chk_true(target_mass > 0)
  }

  chk::chk_true(is.list(datasets) && length(datasets) > 0)
  chk::chk_true(is.list(template))
  chk::chk_number(epsilon)
  chk::chk_true(epsilon > 0)
  chk::chk_number(rho1)
  chk::chk_true(rho1 > 0)
  chk::chk_number(rho2)
  chk::chk_true(rho2 > 0)
  chk::chk_number(max_outer)
  chk::chk_true(max_outer > 0)
  chk::chk_number(max_inner)
  chk::chk_true(max_inner > 0)
  chk::chk_number(tol_outer)
  chk::chk_true(tol_outer > 0)
  chk::chk_number(tol_inner)
  chk::chk_true(tol_inner > 0)
  chk::chk_logical(learn_template_features)
  chk::chk_number(delta)
  chk::chk_true(delta >= 0)
  chk::chk_logical(parallel)
  chk::chk_number(ncores)
  chk::chk_true(ncores >= 1)
  chk::chk_number(lambda_prior)
  chk::chk_true(lambda_prior >= 0)
  if (!is.null(prior_sigma)) {
    chk::chk_number(prior_sigma)
    chk::chk_true(is.finite(prior_sigma))
    chk::chk_true(prior_sigma > 0)
  }

  if (!is.null(constraint_fields)) {
    chk::chk_character(constraint_fields)
    if (length(constraint_fields) == 0L || anyNA(constraint_fields) || any(constraint_fields == "")) {
      stop("constraint_fields must be a non-empty character vector", call. = FALSE)
    }
  }
  if (!is.null(prior_fields)) {
    chk::chk_character(prior_fields)
    if (length(prior_fields) == 0L || anyNA(prior_fields) || any(prior_fields == "")) {
      stop("prior_fields must be a non-empty character vector", call. = FALSE)
    }
  } else if (lambda_prior > 0) {
    # Still allow per-dataset `prior_map`, but enforce that something is provided.
    if (!all(vapply(datasets, function(d) !is.null(d$prior_map), logical(1)))) {
      stop("lambda_prior > 0 requires prior_fields or per-dataset prior_map", call. = FALSE)
    }
  }

  if (is.null(template$Y) || is.null(template$beta)) {
    stop("template must contain Y and beta", call. = FALSE)
  }
  template$Y <- as.matrix(template$Y)
  template$beta <- as.numeric(template$beta)
  M <- nrow(template$Y)
  if (length(template$beta) != M) stop("template$beta length must match nrow(Y)",
                                       call. = FALSE)

  if (isTRUE(learn_template_features)) {
    if (any(vapply(datasets, function(d) is.null(d$F), logical(1)))) {
      stop("learn_template_features=TRUE requires each dataset to provide F",
           call. = FALSE)
    }
    D <- ncol(as.matrix(datasets[[1]]$F))
    for (k in seq_along(datasets)) {
      Fk <- as.matrix(datasets[[k]]$F)
      if (ncol(Fk) != D) stop("All dataset F matrices must have same ncol",
                              call. = FALSE)
    }
    if (is.null(template$G)) {
      template$G <- matrix(0, M, D)
    } else {
      template$G <- as.matrix(template$G)
      if (nrow(template$G) != M || ncol(template$G) != D) {
        stop("template$G must have dimensions nrow(Y) x D", call. = FALSE)
      }
    }
  }

  if (is.null(omega)) {
    omega <- rep(1 / length(datasets), length(datasets))
  } else {
    omega <- as.numeric(omega)
    if (length(omega) != length(datasets)) stop("omega length must match datasets",
                                                call. = FALSE)
    if (any(omega < 0)) stop("omega must be nonnegative", call. = FALSE)
    s <- sum(omega)
    if (s <= 0) stop("omega must have positive sum", call. = FALSE)
    omega <- omega / s
  }

  subj_mass <- vapply(datasets, function(d) sum(as.numeric(d$alpha)), numeric(1))
  m0 <- if (is.character(target_mass)) {
    if (target_mass == "mean_alpha") mean(subj_mass) else NULL
  } else {
    as.numeric(target_mass)
  }

  fits <- vector("list", length(datasets))

  group_Y <- NULL
  if (!is.null(constraint_fields)) {
    missing_template <- setdiff(constraint_fields, names(template))
    if (length(missing_template)) {
      stop("template is missing constraint field(s): ",
           paste(missing_template, collapse = ", "),
           call. = FALSE)
    }
    group_terms <- lapply(constraint_fields, function(f) {
      v <- template[[f]]
      if (length(v) != M) stop("template$", f, " must have length nrow(template$Y)",
                               call. = FALSE)
      v
    })
    group_Y <- do.call(interaction, c(group_terms, list(drop = TRUE, lex.order = TRUE)))
  }

  prior_key_Y <- NULL
  if (!is.null(prior_fields)) {
    missing_template <- setdiff(prior_fields, names(template))
    if (length(missing_template)) {
      stop("template is missing prior field(s): ",
           paste(missing_template, collapse = ", "),
           call. = FALSE)
    }
    prior_terms <- lapply(prior_fields, function(f) {
      v <- template[[f]]
      if (length(v) != M) stop("template$", f, " must have length nrow(template$Y)",
                               call. = FALSE)
      v
    })
    prior_key_Y <- do.call(interaction, c(prior_terms, list(drop = TRUE, lex.order = TRUE)))
  }

  subject_solve_one <- function(d, beta, G) {
    X <- as.matrix(d$X)
    alpha <- as.numeric(d$alpha)
    if (sum(alpha) <= 0) stop("dataset alpha must have positive mass", call. = FALSE)
    F <- if (!is.null(d$F)) as.matrix(d$F) else NULL
    group_X <- NULL
    if (!is.null(constraint_fields)) {
      missing_dataset <- setdiff(constraint_fields, names(d))
      if (length(missing_dataset)) {
        stop("dataset is missing constraint field(s): ",
             paste(missing_dataset, collapse = ", "),
             call. = FALSE)
      }
      group_terms <- lapply(constraint_fields, function(f) {
        v <- d[[f]]
        if (length(v) != nrow(X)) stop("dataset$", f, " must have length nrow(dataset$X)",
                                       call. = FALSE)
        v
      })
      group_X <- do.call(interaction, c(group_terms, list(drop = TRUE, lex.order = TRUE)))
    }

    prior_map_k <- NULL
    if (lambda_prior > 0) {
      if (!is.null(prior_fields)) {
        missing_dataset <- setdiff(prior_fields, names(d))
        if (length(missing_dataset)) {
          stop("dataset is missing prior field(s): ",
               paste(missing_dataset, collapse = ", "),
               call. = FALSE)
        }
        prior_terms <- lapply(prior_fields, function(f) {
          v <- d[[f]]
          if (length(v) != nrow(X)) stop("dataset$", f, " must have length nrow(dataset$X)",
                                         call. = FALSE)
          v
        })
        key_X <- do.call(interaction, c(prior_terms, list(drop = TRUE, lex.order = TRUE)))
        prior_map_k <- match(key_X, prior_key_Y)
      } else {
        prior_map_k <- d$prior_map
      }
    }

    cost <- uot_build_cost(
      X = X, Y = template$Y,
      F = F, G = G,
      lambda_anat = lambda_anat,
      lambda_feat = lambda_feat,
      prior_map = prior_map_k,
      lambda_prior = lambda_prior,
      prior_sigma = prior_sigma,
      neighbor_mode = neighbor_mode,
      k_neighbors = k_neighbors,
      radius = radius,
      maxk = maxk,
      min_neighbors = min_neighbors,
      group_X = group_X,
      group_Y = group_Y,
      dense_max_bytes = dense_max_bytes,
      ensure_cols = ensure_cols
    )

    sol <- uot_ti_sinkhorn_kl(
      cost = cost,
      alpha = alpha,
      beta = beta,
      epsilon = epsilon,
      rho1 = rho1,
      rho2 = rho2,
      max_iter = max_inner,
      tol = tol_inner
    )

    mom <- .uot_logpi2_meanF(
      cost = cost,
      alpha = alpha,
      beta = beta,
      fbar = sol$fbar,
      gbar = sol$gbar,
      epsilon = epsilon,
      F = if (learn_template_features) F else NULL
    )

    # Rescale pi2 for aggregation
    mass_k <- if (rescale_pi2 == "alpha") sum(alpha) else 1.0
    pi2 <- .normalize_from_log(mom$log_pi2, mass = mass_k)

    out <- list(
      cost = cost,
      sol = sol,
      pi2 = pi2
    )
    if (learn_template_features) {
      out$mean_F <- mom$mean_F
    }
    out
  }

  beta <- template$beta
  G <- template$G

  for (outer in seq_len(max_outer)) {
    if (isTRUE(verbose)) {
      message("multiset_uot_align: outer iter ", outer, "/", max_outer)
    }

    do_parallel <- isTRUE(parallel) && (.Platform$OS.type != "windows") &&
      requireNamespace("parallel", quietly = TRUE) && ncores > 1

    results <- if (do_parallel) {
      parallel::mclapply(
        datasets,
        subject_solve_one,
        beta = beta,
        G = G,
        mc.cores = as.integer(ncores)
      )
    } else {
      lapply(datasets, subject_solve_one, beta = beta, G = G)
    }

    pi2_list <- lapply(results, `[[`, "pi2")
    beta_new <- Reduce(`+`, Map(`*`, omega, pi2_list))

    if (!is.null(m0)) {
      beta_new <- beta_new * (m0 / (sum(beta_new) + 1e-16))
    }

    if (learn_template_features) {
      # Compute numerator = sum_k omega_k * pi2_k * mean_F_k
      D <- ncol(results[[1]]$mean_F)
      num <- matrix(0, M, D)
      den <- beta_new
      for (k in seq_along(results)) {
        mean_F <- results[[k]]$mean_F
        pi2_k <- as.numeric(pi2_list[[k]])
        num <- num + omega[k] * sweep(mean_F, 1, pi2_k, `*`)
      }
      G_new <- num / (den + delta)
      G <- G_new
    }

    d_beta <- max(abs(beta_new - beta))
    beta <- beta_new
    fits <- results

    if (d_beta <= tol_outer) break
  }

  template$beta <- beta
  if (learn_template_features) template$G <- G

  list(
    template = template,
    fits = fits,
    params = list(
      epsilon = epsilon,
      rho1 = rho1,
      rho2 = rho2,
      lambda_anat = lambda_anat,
      lambda_feat = lambda_feat,
      neighbor_mode = neighbor_mode,
      k_neighbors = k_neighbors,
      radius = radius,
      maxk = maxk,
      min_neighbors = min_neighbors,
      constraint_fields = constraint_fields,
      dense_max_bytes = dense_max_bytes,
      ensure_cols = ensure_cols,
      rescale_pi2 = rescale_pi2,
      target_mass = target_mass,
      learn_template_features = learn_template_features,
      delta = delta
    )
  )
}
