#' Multiscale manifold alignment (Wang-Mahadevan style)
#'
#' Experimental multiset manifold alignment with a QR-free spectral hierarchy.
#' The implementation builds a joint graph Laplacian, compresses each domain
#' with truncated SVD, forms a reduced operator, and extracts nested multiscale
#' bases from its spectrum.
#'
#' @param data A `hyperdesign` object or list of matrices.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `multiblock_biprojector`-compatible object with subclass
#'   `multiscale_manifold_align`.
#' @export
multiscale_manifold_align <- function(data, ...) {
  UseMethod("multiscale_manifold_align")
}

.msa_preprocess_domains <- function(domains, preproc) {
  n_domains <- length(domains)
  out <- domains
  fitted_preproc <- vector("list", n_domains)

  if (is.null(preproc)) {
    for (i in seq_along(domains)) {
      out[[i]]$x <- as.matrix(domains[[i]]$x)
      fitted_preproc[i] <- list(NULL)
    }
    names(fitted_preproc) <- names(domains)
    return(list(domains = out, preproc = fitted_preproc))
  }

  is_stateful <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")
  if (is_stateful) {
    preproc_list <- replicate(
      n_domains,
      unserialize(serialize(preproc, connection = NULL)),
      simplify = FALSE
    )
  } else if (is.function(preproc)) {
    preproc_list <- rep(list(preproc), n_domains)
  } else if (is.list(preproc) && length(preproc) == n_domains) {
    preproc_list <- lapply(preproc, function(pre_i) {
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        unserialize(serialize(pre_i, connection = NULL))
      } else {
        pre_i
      }
    })
  } else {
    stop("`preproc` must be NULL, a function/preprocessor, or a list matching domain count.", call. = FALSE)
  }

  for (i in seq_along(domains)) {
    Xi <- as.matrix(domains[[i]]$x)
    pre_i <- preproc_list[[i]]
    if (is.null(pre_i)) {
      out[[i]]$x <- Xi
      fitted_preproc[i] <- list(NULL)
      next
    }
    if (is.function(pre_i)) {
      out[[i]]$x <- pre_i(Xi)
      fitted_preproc[i] <- list(pre_i)
      next
    }
    if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
      ft <- multivarious::fit_transform(pre_i, Xi)
      out[[i]]$x <- ft$transformed
      fitted_preproc[i] <- list(ft$preproc)
    } else {
      templ <- multivarious::prep(pre_i)
      Xi_tf <- multivarious::init_transform(templ, Xi)
      out[[i]]$x <- Xi_tf
      fitted_preproc[i] <- list(attr(Xi_tf, "preproc"))
      if (is.null(fitted_preproc[[i]])) {
        fitted_preproc[i] <- list(templ)
      }
    }
  }
  names(fitted_preproc) <- names(domains)
  list(domains = out, preproc = fitted_preproc)
}

.msa_build_within_graph <- function(X, knn, sigma = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n <= 1L) {
    A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    return(list(A = A, degree = numeric(n)))
  }

  k_use <- min(as.integer(knn), n - 1L)
  sigma_use <- if (is.null(sigma)) choose_sigma(X) else sigma

  gw <- safe_compute(
    neighborweights::graph_weights(
      X,
      k = k_use,
      weight_mode = "heat",
      sigma = sigma_use,
      neighbor_mode = "knn"
    ),
    "Within-domain graph construction failed in multiscale_manifold_align()."
  )

  A <- neighborweights::adjacency(gw)
  A <- (A + Matrix::t(A)) / 2
  Matrix::diag(A) <- 0
  A <- Matrix::drop0(methods::as(A, "dgCMatrix"))
  list(A = A, degree = as.numeric(Matrix::rowSums(A)))
}

.msa_extract_labels <- function(domains, y_quo) {
  if (rlang::quo_is_null(y_quo)) return(NULL)
  y_expr <- rlang::get_expr(y_quo)
  y_name <- if (is.character(y_expr) && length(y_expr) == 1L) {
    y_expr
  } else {
    rlang::as_name(y_expr)
  }
  lapply(domains, function(dom) {
    if (is.null(dom$design) || !is.data.frame(dom$design)) {
      stop("`y` provided but at least one domain has no design data.frame.", call. = FALSE)
    }
    if (!y_name %in% names(dom$design)) {
      stop("Label column `", y_name, "` not found in a domain design frame.", call. = FALSE)
    }
    dom$design[[y_name]]
  })
}

.msa_build_label_correspondences <- function(labels_list, max_pairs_per_label, weight, seed) {
  set.seed(seed)
  m <- length(labels_list)
  rows <- vector("list", m * (m - 1L) / 2L)
  ptr <- 0L
  for (i in seq_len(m - 1L)) {
    li <- labels_list[[i]]
    for (j in (i + 1L):m) {
      lj <- labels_list[[j]]
      shared <- intersect(unique(li[!is.na(li)]), unique(lj[!is.na(lj)]))
      if (!length(shared)) next
      for (lab in shared) {
        idx_i <- which(li == lab)
        idx_j <- which(lj == lab)
        n_take <- min(length(idx_i), length(idx_j), as.integer(max_pairs_per_label))
        if (n_take < 1L) next
        ii <- sample(idx_i, size = n_take, replace = FALSE)
        jj <- sample(idx_j, size = n_take, replace = FALSE)
        ptr <- ptr + 1L
        rows[[ptr]] <- data.frame(
          domain_i = i,
          index_i = ii,
          domain_j = j,
          index_j = jj,
          weight = rep(weight, n_take),
          source = rep("label", n_take),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (ptr == 0L) {
    return(data.frame(
      domain_i = integer(0), index_i = integer(0), domain_j = integer(0), index_j = integer(0),
      weight = numeric(0), source = character(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows[seq_len(ptr)])
}

.msa_build_mnn_correspondences <- function(domains, mnn_k, max_pairs_per_pair, weight) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    return(data.frame(
      domain_i = integer(0), index_i = integer(0), domain_j = integer(0), index_j = integer(0),
      weight = numeric(0), source = character(0), stringsAsFactors = FALSE
    ))
  }

  m <- length(domains)
  rows <- vector("list", m * (m - 1L) / 2L)
  ptr <- 0L
  for (i in seq_len(m - 1L)) {
    Xi <- scale(as.matrix(domains[[i]]$x), center = TRUE, scale = FALSE)
    for (j in (i + 1L):m) {
      Xj <- scale(as.matrix(domains[[j]]$x), center = TRUE, scale = FALSE)
      if (!nrow(Xi) || !nrow(Xj)) next
      k_ij <- min(as.integer(mnn_k), nrow(Xj))
      k_ji <- min(as.integer(mnn_k), nrow(Xi))
      nn_ij <- RANN::nn2(Xj, Xi, k = k_ij)$nn.idx
      nn_ji <- RANN::nn2(Xi, Xj, k = k_ji)$nn.idx

      pairs <- vector("list", nrow(Xi))
      n_pairs <- 0L
      for (a in seq_len(nrow(Xi))) {
        b <- nn_ij[a, 1L]
        if (a %in% nn_ji[b, ]) {
          n_pairs <- n_pairs + 1L
          pairs[[n_pairs]] <- c(a, b)
        }
      }
      if (n_pairs < 1L) next
      mat <- do.call(rbind, pairs[seq_len(n_pairs)])
      if (nrow(mat) > max_pairs_per_pair) {
        keep <- sample.int(nrow(mat), size = as.integer(max_pairs_per_pair), replace = FALSE)
        mat <- mat[keep, , drop = FALSE]
      }
      ptr <- ptr + 1L
      rows[[ptr]] <- data.frame(
        domain_i = i,
        index_i = mat[, 1L],
        domain_j = j,
        index_j = mat[, 2L],
        weight = rep(weight, nrow(mat)),
        source = rep("mnn", nrow(mat)),
        stringsAsFactors = FALSE
      )
    }
  }

  if (ptr == 0L) {
    return(data.frame(
      domain_i = integer(0), index_i = integer(0), domain_j = integer(0), index_j = integer(0),
      weight = numeric(0), source = character(0), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, rows[seq_len(ptr)])
}

.msa_validate_correspondences <- function(corr, domain_sizes, default_weight = 1) {
  if (is.null(corr)) {
    return(data.frame(
      domain_i = integer(0), index_i = integer(0), domain_j = integer(0), index_j = integer(0),
      weight = numeric(0), source = character(0), stringsAsFactors = FALSE
    ))
  }
  req <- c("domain_i", "index_i", "domain_j", "index_j")
  if (!is.data.frame(corr) || !all(req %in% names(corr))) {
    stop("`correspondences` must be a data.frame with columns domain_i,index_i,domain_j,index_j.", call. = FALSE)
  }
  out <- corr
  out$domain_i <- as.integer(out$domain_i)
  out$domain_j <- as.integer(out$domain_j)
  out$index_i <- as.integer(out$index_i)
  out$index_j <- as.integer(out$index_j)
  out$weight <- if ("weight" %in% names(out)) as.numeric(out$weight) else rep(default_weight, nrow(out))
  if (!("source" %in% names(out))) out$source <- "provided"

  m <- length(domain_sizes)
  keep <- is.finite(out$weight) & out$weight > 0 &
    out$domain_i >= 1L & out$domain_i <= m &
    out$domain_j >= 1L & out$domain_j <= m &
    out$domain_i != out$domain_j
  out <- out[keep, , drop = FALSE]
  if (!nrow(out)) return(out)

  valid_idx <- vapply(seq_len(nrow(out)), function(r) {
    out$index_i[r] >= 1L && out$index_i[r] <= domain_sizes[out$domain_i[r]] &&
      out$index_j[r] >= 1L && out$index_j[r] <= domain_sizes[out$domain_j[r]]
  }, logical(1))
  out[valid_idx, , drop = FALSE]
}

.msa_build_joint_laplacian <- function(within_graphs, correspondences, domain_sizes) {
  n_total <- sum(domain_sizes)
  L_within <- Matrix::bdiag(lapply(within_graphs, function(g) {
    Matrix::Diagonal(x = g$degree) - g$A
  }))

  if (is.null(correspondences) || !nrow(correspondences)) {
    return(methods::as(L_within, "dgCMatrix"))
  }

  offsets <- c(0L, cumsum(domain_sizes))
  gi <- offsets[correspondences$domain_i] + correspondences$index_i
  gj <- offsets[correspondences$domain_j] + correspondences$index_j
  w <- correspondences$weight

  i <- c(gi, gj, gi, gj)
  j <- c(gi, gj, gj, gi)
  x <- c(w, w, -w, -w)
  L_cross <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n_total, n_total))
  methods::as(L_within + L_cross, "dgCMatrix")
}

.msa_domain_svd <- function(X, degree, rank_per_domain, svd_tol) {
  X <- as.matrix(X)
  dhalf <- sqrt(pmax(as.numeric(degree), 1e-12))
  A <- Matrix::t(X) %*% Matrix::Diagonal(x = dhalf)
  k <- min(as.integer(rank_per_domain), nrow(A), ncol(A))
  if (k < 1L) {
    return(list(V = matrix(0, nrow(X), 1L), d = 0))
  }

  sv <- if (requireNamespace("irlba", quietly = TRUE) && k < min(dim(A))) {
    irlba::irlba(A, nv = k, nu = 0)
  } else {
    s <- svd(as.matrix(A), nu = 0, nv = k)
    list(v = s$v, d = s$d[seq_len(k)])
  }

  keep <- which(is.finite(sv$d) & sv$d >= svd_tol * max(sv$d))
  if (!length(keep)) keep <- 1L
  list(V = as.matrix(sv$v[, keep, drop = FALSE]), d = as.numeric(sv$d[keep]))
}

.msa_basis_rank <- function(ctrl, r, ncomp) {
  max_rank <- max(1L, r - 1L)
  base_rank <- min(max_rank, max(as.integer(ncomp), 2L))
  if (identical(ctrl$backend, "randomized_dwt")) {
    target <- max(base_rank, as.integer(ctrl$randomized_sketch_size))
  } else {
    target <- max(base_rank, min(as.integer(ctrl$eigen_k_max), max_rank))
  }
  min(target, as.integer(ctrl$eigen_k_max), max_rank)
}

.msa_smallest_eig_exact <- function(A, k) {
  A <- as.matrix((A + Matrix::t(A)) / 2)
  n <- nrow(A)
  if (k < 1L) {
    stop("Need at least one eigenpair.", call. = FALSE)
  }

  if (requireNamespace("RSpectra", quietly = TRUE) && k < (n - 1L)) {
    eig <- safe_compute(
      RSpectra::eigs_sym(A, k = k, which = "SA"),
      "Spectral decomposition failed in multiscale_manifold_align()."
    )
    vals <- as.numeric(eig$values)
    vecs <- as.matrix(eig$vectors)
  } else {
    eig <- eigen(as.matrix(A), symmetric = TRUE)
    vals <- as.numeric(eig$values[seq_len(k)])
    vecs <- as.matrix(eig$vectors[, seq_len(k), drop = FALSE])
  }

  ord <- order(vals, decreasing = FALSE)
  list(values = vals[ord], vectors = vecs[, ord, drop = FALSE])
}

.msa_smallest_eig_randomized <- function(A, k, ctrl) {
  A <- as.matrix((A + Matrix::t(A)) / 2)
  n <- nrow(A)
  if (k < 1L) {
    stop("Need at least one eigenpair.", call. = FALSE)
  }

  ell <- min(n, k + as.integer(ctrl$randomized_oversample))
  set.seed(as.integer(ctrl$seed))
  Omega <- matrix(stats::rnorm(n * ell), n, ell)

  shifted <- A + diag(as.numeric(ctrl$shift_invert_delta), n)
  R <- chol(shifted)
  solve_shift <- function(B) {
    backsolve(R, forwardsolve(t(R), B))
  }

  Y <- solve_shift(Omega)
  if (as.integer(ctrl$randomized_power_iters) > 0L) {
    for (iter in seq_len(as.integer(ctrl$randomized_power_iters))) {
      Y <- solve_shift(qr.Q(qr(Y)))
    }
  }

  Q <- qr.Q(qr(Y))
  B <- crossprod(Q, A %*% Q)
  B <- (B + t(B)) / 2
  eig <- eigen(B, symmetric = TRUE)
  ord <- order(eig$values, decreasing = FALSE)
  take <- seq_len(min(k, length(ord)))
  vals <- as.numeric(eig$values[ord][take])
  vecs <- Q %*% eig$vectors[, ord, drop = FALSE][, take, drop = FALSE]
  list(values = vals, vectors = vecs)
}

.msa_build_filter_bank <- function(node_eig, eigvals, ctrl) {
  n_nodes <- nrow(node_eig)
  scaling_levels <- list()
  wavelet_levels <- list()
  scaling_weights <- list()
  wavelet_weights <- list()
  level_dims <- integer(0)
  times <- numeric(0)

  make_level <- function(weights) {
    weights <- as.numeric(weights)
    if (!length(weights) || !any(is.finite(weights) & weights > 0)) {
      return(matrix(0, n_nodes, 0L))
    }
    wmax <- max(weights, na.rm = TRUE)
    active <- which(is.finite(weights) & weights > ctrl$eps_rank * wmax & weights > 0)
    if (!length(active)) {
      return(matrix(0, n_nodes, 0L))
    }
    sweep(node_eig[, active, drop = FALSE], 2, sqrt(weights[active]), `*`)
  }

  for (j in seq_len(as.integer(ctrl$max_levels))) {
    t_j <- 2^(j - 1L)
    h_j <- 1 / (1 + t_j * eigvals)
    h_next <- 1 / (1 + (2 * t_j) * eigvals)
    g_j <- pmax(h_j - h_next, 0)

    scaling_levels[[j]] <- make_level(h_j)
    wavelet_levels[[j]] <- make_level(g_j)
    scaling_weights[[j]] <- h_j
    wavelet_weights[[j]] <- g_j
    level_dims[j] <- ncol(scaling_levels[[j]])
    times[j] <- t_j

    if (ncol(scaling_levels[[j]]) <= 1L) {
      break
    }
  }

  if (!length(scaling_levels)) {
    scaling_levels[[1L]] <- node_eig[, seq_len(min(1L, ncol(node_eig))), drop = FALSE]
    wavelet_levels[[1L]] <- matrix(0, n_nodes, 0L)
    scaling_weights[[1L]] <- rep(1, ncol(node_eig))
    wavelet_weights[[1L]] <- rep(0, ncol(node_eig))
    level_dims[1L] <- ncol(scaling_levels[[1L]])
    times[1L] <- 1
  }

  scaling_levels <- scaling_levels[seq_along(level_dims)]
  wavelet_levels <- wavelet_levels[seq_along(level_dims)]
  scaling_weights <- scaling_weights[seq_along(level_dims)]
  wavelet_weights <- wavelet_weights[seq_along(level_dims)]
  times <- as.numeric(times)
  level_dims <- as.integer(level_dims)

  parts <- c(
    Filter(function(x) ncol(x) > 0L, wavelet_levels),
    list(scaling_levels[[length(scaling_levels)]])
  )
  features <- do.call(cbind, Filter(function(x) ncol(x) > 0L, parts))
  if (is.null(features) || !ncol(features)) {
    features <- node_eig[, seq_len(min(ncol(node_eig), 1L)), drop = FALSE]
  }

  list(
    levels = scaling_levels,
    wavelet_levels = wavelet_levels,
    scaling_weights = scaling_weights,
    wavelet_weights = wavelet_weights,
    level_dims = level_dims,
    times = times,
    features = as.matrix(features)
  )
}

.msa_compress_multiscale_features <- function(features, ncomp) {
  features <- as.matrix(features)
  k <- min(as.integer(ncomp), ncol(features), nrow(features))
  if (k < 1L) {
    stop("No non-empty multiscale features were constructed.", call. = FALSE)
  }
  if (ncol(features) == k) {
    return(features[, seq_len(k), drop = FALSE])
  }
  if (requireNamespace("irlba", quietly = TRUE) && k < min(dim(features))) {
    sv <- irlba::irlba(features, nv = k, nu = 0)
    return(as.matrix(features %*% sv$v[, seq_len(k), drop = FALSE]))
  }
  sv <- svd(features, nu = 0, nv = k)
  as.matrix(features %*% sv$v[, seq_len(k), drop = FALSE])
}

.msa_fit_loadings <- function(X, scores) {
  X <- as.matrix(X)
  scores <- as.matrix(scores)
  p <- ncol(X)
  n <- nrow(X)
  ridge_base <- mean(colSums(X * X) + 1e-12)
  ridge <- max(1e-8, 1e-6 * ridge_base)

  if (p <= n) {
    XtX <- crossprod(X)
    return(solve(XtX + diag(ridge, p), crossprod(X, scores)))
  }

  XXt <- tcrossprod(X)
  t(X) %*% solve(XXt + diag(ridge, n), scores)
}

.msa_rescale_correspondences <- function(corr, cross_edge_weight) {
  if (is.null(corr) || !nrow(corr)) return(corr)
  out <- corr
  out$weight <- pmax(as.numeric(out$weight), 0) * as.numeric(cross_edge_weight)
  out
}

.msa_build_joint_adjacency <- function(within_graphs, correspondences, domain_sizes) {
  n_total <- sum(domain_sizes)
  A_within <- Matrix::bdiag(lapply(within_graphs, `[[`, "A"))

  if (is.null(correspondences) || !nrow(correspondences)) {
    A_joint <- methods::as(A_within, "dgCMatrix")
    Matrix::diag(A_joint) <- 0
    return(Matrix::drop0(A_joint))
  }

  offsets <- c(0L, cumsum(domain_sizes))
  gi <- offsets[correspondences$domain_i] + correspondences$index_i
  gj <- offsets[correspondences$domain_j] + correspondences$index_j
  w <- correspondences$weight

  A_cross <- Matrix::sparseMatrix(
    i = c(gi, gj),
    j = c(gj, gi),
    x = c(w, w),
    dims = c(n_total, n_total)
  )

  A_joint <- methods::as(A_within + A_cross, "dgCMatrix")
  Matrix::diag(A_joint) <- 0
  Matrix::drop0(A_joint)
}

.msa_build_joint_laplacian <- function(A_joint) {
  degree <- as.numeric(Matrix::rowSums(A_joint))
  degree[degree <= 0] <- 1
  list(
    laplacian = methods::as(Matrix::Diagonal(x = degree) - A_joint, "dgCMatrix"),
    degree = degree
  )
}

.msa_diffusion_eig_exact <- function(A, k) {
  A <- as.matrix((A + Matrix::t(A)) / 2)
  n <- nrow(A)
  if (k < 1L) {
    stop("Need at least one eigenpair.", call. = FALSE)
  }

  if (requireNamespace("RSpectra", quietly = TRUE) && k < (n - 1L)) {
    eig <- safe_compute(
      RSpectra::eigs_sym(A, k = k, which = "LA"),
      "Diffusion eigendecomposition failed in multiscale_manifold_align()."
    )
    vals <- as.numeric(eig$values)
    vecs <- as.matrix(eig$vectors)
  } else {
    eig <- eigen(A, symmetric = TRUE)
    vals <- as.numeric(eig$values[seq_len(k)])
    vecs <- as.matrix(eig$vectors[, seq_len(k), drop = FALSE])
  }

  ord <- order(vals, decreasing = TRUE)
  list(values = vals[ord], vectors = vecs[, ord, drop = FALSE])
}

.msa_diffusion_eig_randomized <- function(A, k, ctrl) {
  A <- as.matrix((A + Matrix::t(A)) / 2)
  n <- nrow(A)
  if (k < 1L) {
    stop("Need at least one eigenpair.", call. = FALSE)
  }

  ell <- min(n, max(k, as.integer(ctrl$randomized_sketch_size)) + as.integer(ctrl$randomized_oversample))
  set.seed(as.integer(ctrl$seed))
  Omega <- matrix(stats::rnorm(n * ell), n, ell)

  Y <- A %*% Omega
  if (as.integer(ctrl$randomized_power_iters) > 0L) {
    for (iter in seq_len(as.integer(ctrl$randomized_power_iters))) {
      Y <- A %*% qr.Q(qr(Y))
    }
  }

  Q <- qr.Q(qr(Y))
  B <- crossprod(Q, A %*% Q)
  B <- (B + t(B)) / 2
  eig <- eigen(B, symmetric = TRUE)
  ord <- order(eig$values, decreasing = TRUE)
  take <- seq_len(min(k, length(ord)))
  vals <- as.numeric(eig$values[ord][take])
  vecs <- Q %*% eig$vectors[, ord, drop = FALSE][, take, drop = FALSE]
  list(values = vals, vectors = vecs)
}

.msa_compress_scale_block <- function(block, ctrl) {
  block <- as.matrix(block)
  if (!nrow(block) || !ncol(block)) {
    return(matrix(0, nrow(block), 0L))
  }

  k <- min(nrow(block), ncol(block), as.integer(ctrl$max_basis_per_scale))
  if (k < 1L) {
    return(matrix(0, nrow(block), 0L))
  }

  sv <- if (requireNamespace("irlba", quietly = TRUE) && k < min(dim(block))) {
    irlba::irlba(block, nu = k, nv = 0)
  } else {
    s <- svd(block, nu = k, nv = 0)
    list(u = s$u[, seq_len(k), drop = FALSE], d = s$d[seq_len(k)])
  }

  keep <- which(is.finite(sv$d) & sv$d >= ctrl$svd_tol * max(sv$d))
  if (!length(keep)) keep <- 1L
  as.matrix(sv$u[, keep, drop = FALSE])
}

.msa_correspondence_globals <- function(corr, domain_sizes) {
  if (is.null(corr) || !nrow(corr)) return(NULL)
  offsets <- c(0L, cumsum(domain_sizes))
  list(
    gi = offsets[corr$domain_i] + corr$index_i,
    gj = offsets[corr$domain_j] + corr$index_j
  )
}

.msa_scale_supervision_score <- function(block, corr_global) {
  if (is.null(corr_global) || !ncol(block)) return(NA_real_)
  matched_d2 <- mean(rowSums((block[corr_global$gi, , drop = FALSE] - block[corr_global$gj, , drop = FALSE])^2))
  baseline <- mean(rowSums(block^2)) + 1e-12
  1 / (1 + matched_d2 / baseline)
}

.msa_normalize_positive <- function(x) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  xmax <- max(x, na.rm = TRUE)
  if (!is.finite(xmax) || xmax <= 0) {
    rep(1, length(x))
  } else {
    x / xmax
  }
}

.msa_build_diffusion_dictionary <- function(node_eig, eigvals, corr_df, domain_sizes, ctrl) {
  n_nodes <- nrow(node_eig)
  lambda <- pmin(pmax(as.numeric(eigvals), 0), 1)
  corr_global <- .msa_correspondence_globals(corr_df, domain_sizes)

  make_level <- function(weights) {
    if (!length(weights) || !any(is.finite(weights) & weights > 0)) {
      return(matrix(0, n_nodes, 0L))
    }
    wmax <- max(weights, na.rm = TRUE)
    active <- which(is.finite(weights) & weights > ctrl$eps_rank * wmax & weights > 0)
    if (!length(active)) {
      return(matrix(0, n_nodes, 0L))
    }
    sweep(node_eig[, active, drop = FALSE], 2, weights[active], `*`)
  }

  raw_blocks <- list()
  block_labels <- character(0)
  block_times <- numeric(0)
  scaling_levels <- list()
  wavelet_levels <- list()
  scaling_weights <- list()
  wavelet_weights <- list()
  level_dims <- integer(0)

  for (j in seq_len(as.integer(ctrl$max_levels))) {
    t_j <- 2^(j - 1L)
    h_j <- lambda^t_j
    h_next <- lambda^(2 * t_j)
    g_j <- sqrt(pmax(h_j^2 - h_next^2, 0))

    scaling_levels[[j]] <- make_level(h_j)
    wavelet_levels[[j]] <- make_level(g_j)
    scaling_weights[[j]] <- h_j
    wavelet_weights[[j]] <- g_j
    level_dims[j] <- ncol(scaling_levels[[j]])

    if (ncol(wavelet_levels[[j]]) > 0L) {
      raw_blocks[[length(raw_blocks) + 1L]] <- wavelet_levels[[j]]
      block_labels <- c(block_labels, paste0("wavelet_", j))
      block_times <- c(block_times, t_j)
    }

    if (ncol(scaling_levels[[j]]) <= 1L) {
      break
    }
  }

  if (length(scaling_levels)) {
    raw_blocks[[length(raw_blocks) + 1L]] <- scaling_levels[[length(scaling_levels)]]
    block_labels <- c(block_labels, paste0("scaling_", length(scaling_levels)))
    block_times <- c(block_times, 2^(length(scaling_levels) - 1L))
  }

  if (!length(raw_blocks)) {
    raw_blocks <- list(node_eig[, seq_len(min(ncol(node_eig), 1L)), drop = FALSE])
    block_labels <- "wavelet_1"
    block_times <- 1
  }

  compressed_blocks <- lapply(raw_blocks, .msa_compress_scale_block, ctrl = ctrl)
  energy_scores <- vapply(raw_blocks, function(B) {
    sqrt(sum(B^2) / max(1, nrow(B) * max(1, ncol(B))))
  }, numeric(1))
  supervision_scores <- vapply(compressed_blocks, .msa_scale_supervision_score, numeric(1), corr_global = corr_global)

  energy_norm <- .msa_normalize_positive(energy_scores)
  supervision_norm <- .msa_normalize_positive(supervision_scores)
  scale_weights <- switch(
    ctrl$scale_selection,
    energy = energy_norm,
    supervised = supervision_norm,
    hybrid = (1 - ctrl$supervision_weight) * energy_norm + ctrl$supervision_weight * supervision_norm
  )

  keep <- which(scale_weights >= max(scale_weights) * ctrl$min_scale_weight)
  if (!length(keep)) keep <- seq_along(compressed_blocks)
  selected_blocks <- Map(
    function(B, w) as.matrix(B) * sqrt(max(w, 1e-8)),
    compressed_blocks[keep],
    scale_weights[keep]
  )
  dictionary <- do.call(cbind, Filter(function(B) ncol(B) > 0L, selected_blocks))
  if (is.null(dictionary) || !ncol(dictionary)) {
    dictionary <- node_eig[, seq_len(min(ncol(node_eig), 1L)), drop = FALSE]
  }

  list(
    dictionary = as.matrix(dictionary),
    levels = scaling_levels,
    wavelet_levels = wavelet_levels,
    level_dims = as.integer(level_dims),
    scaling_weights = scaling_weights,
    wavelet_weights = wavelet_weights,
    scale_weights = scale_weights,
    scale_labels = block_labels,
    selected_scales = block_labels[keep],
    selected_scale_weights = scale_weights[keep],
    scale_times = block_times,
    energy_scores = energy_scores,
    supervision_scores = supervision_scores
  )
}

.msa_correspondence_margin <- function(scores, corr_df, domain_sizes) {
  if (is.null(corr_df) || !nrow(corr_df)) return(NA_real_)

  offsets <- c(0L, cumsum(domain_sizes))
  margins <- numeric(0)

  pair_keys <- paste(corr_df$domain_i, corr_df$domain_j, sep = "_")
  for (pair_key in unique(pair_keys)) {
    idx <- which(pair_keys == pair_key)
    pair_corr <- corr_df[idx, , drop = FALSE]
    di <- pair_corr$domain_i[1]
    dj <- pair_corr$domain_j[1]
    target_pool <- offsets[dj] + seq_len(domain_sizes[dj])
    query_global <- offsets[di] + pair_corr$index_i
    target_global <- offsets[dj] + pair_corr$index_j

    matched_by_query <- split(target_global, query_global)
    for (query_id in names(matched_by_query)) {
      gi <- as.integer(query_id)
      matched_targets <- unique(as.integer(matched_by_query[[query_id]]))
      query_row <- as.numeric(scores[gi, , drop = FALSE])
      matched_diff <- sweep(scores[matched_targets, , drop = FALSE], 2, query_row, `-`)
      matched_d2 <- min(rowSums(matched_diff^2))
      non_targets <- setdiff(target_pool, matched_targets)
      if (!length(non_targets)) next
      nonmatch_diff <- sweep(scores[non_targets, , drop = FALSE], 2, query_row, `-`)
      nonmatch_d2 <- min(rowSums(nonmatch_diff^2))
      margins <- c(margins, sqrt(nonmatch_d2) - sqrt(matched_d2))
    }
  }

  if (!length(margins)) NA_real_ else mean(margins)
}

.msa_tuning_grid <- function(ctrl) {
  cross_vals <- if (is.null(ctrl$candidate_cross_edge_weight)) {
    sort(unique(pmax(0.25, c(ctrl$cross_edge_weight / 2, ctrl$cross_edge_weight, ctrl$cross_edge_weight * 2))))
  } else {
    sort(unique(as.numeric(ctrl$candidate_cross_edge_weight)))
  }
  rank_vals <- if (is.null(ctrl$candidate_rank_per_domain)) {
    sort(unique(pmax(4L, c(ctrl$rank_per_domain %/% 2L, ctrl$rank_per_domain, ctrl$rank_per_domain * 2L))))
  } else {
    sort(unique(as.integer(ctrl$candidate_rank_per_domain)))
  }
  level_vals <- if (is.null(ctrl$candidate_max_levels)) {
    sort(unique(pmax(2L, c(ctrl$max_levels - 1L, ctrl$max_levels, ctrl$max_levels + 1L))))
  } else {
    sort(unique(as.integer(ctrl$candidate_max_levels)))
  }
  expand.grid(
    cross_edge_weight = cross_vals,
    rank_per_domain = rank_vals,
    max_levels = level_vals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

.msa_fit_core <- function(domains, within, corr_df, domain_names, ncomp, ctrl, proc) {
  domain_sizes <- vapply(domains, function(d) nrow(d$x), integer(1))
  feature_blocks <- feature_block_indices(domains)

  A_joint <- .msa_build_joint_adjacency(within_graphs = within, correspondences = corr_df, domain_sizes = domain_sizes)
  lap_info <- .msa_build_joint_laplacian(A_joint)
  L_joint <- lap_info$laplacian
  d_total <- lap_info$degree
  degree_splits <- split(d_total, rep(seq_along(domain_sizes), domain_sizes))

  svd_list <- lapply(seq_along(domains), function(i) {
    .msa_domain_svd(
      X = domains[[i]]$x,
      degree = degree_splits[[i]],
      rank_per_domain = ctrl$rank_per_domain,
      svd_tol = ctrl$svd_tol
    )
  })

  V_big <- Matrix::bdiag(lapply(svd_list, function(x) Matrix::Matrix(x$V, sparse = TRUE)))
  d_inv_sqrt <- 1 / sqrt(pmax(d_total, 1e-12))
  diffusion_op <- Matrix::Diagonal(x = d_inv_sqrt) %*% A_joint %*% Matrix::Diagonal(x = d_inv_sqrt)
  T_reduced <- Matrix::crossprod(V_big, diffusion_op %*% V_big)
  T_reduced <- Matrix::Matrix((T_reduced + Matrix::t(T_reduced)) / 2, sparse = FALSE)

  r <- nrow(T_reduced)
  k_req <- .msa_basis_rank(ctrl, r = r, ncomp = ncomp)
  if (k_req < 1L) {
    stop("Reduced operator rank is too small for diffusion eigendecomposition.", call. = FALSE)
  }

  spectral_fit <- if (identical(ctrl$backend, "randomized_dwt")) {
    .msa_diffusion_eig_randomized(T_reduced, k = k_req, ctrl = ctrl)
  } else {
    .msa_diffusion_eig_exact(T_reduced, k = k_req)
  }

  keep <- which(is.finite(spectral_fit$values) & spectral_fit$values > ctrl$eigen_tol)
  if (!length(keep)) {
    stop("No non-trivial diffusion eigenvalues were found above `eigen_tol`.", call. = FALSE)
  }

  eigvals <- pmin(pmax(as.numeric(spectral_fit$values[keep]), 0), 1)
  eigvecs <- spectral_fit$vectors[, keep, drop = FALSE]
  node_eig <- as.matrix(V_big %*% eigvecs)

  dictionary_info <- .msa_build_diffusion_dictionary(
    node_eig = node_eig,
    eigvals = eigvals,
    corr_df = corr_df,
    domain_sizes = domain_sizes,
    ctrl = ctrl
  )
  dictionary <- dictionary_info$dictionary

  A_dict <- crossprod(dictionary, as.matrix(L_joint %*% dictionary))
  B_dict <- crossprod(dictionary, sweep(dictionary, 1L, d_total, `*`))
  solved <- .ssma_solve_reduced_geigen(
    A_r = A_dict,
    B_r = B_dict,
    ncomp = as.integer(ncomp),
    eig_tol = ctrl$eigen_tol,
    regularization = ctrl$regularization
  )

  scores <- as.matrix(dictionary %*% solved$alpha)
  objective <- .msa_correspondence_margin(scores, corr_df = corr_df, domain_sizes = domain_sizes)

  domain_offsets <- cumsum(domain_sizes)
  domain_starts <- c(1L, head(domain_offsets, -1L) + 1L)
  loading_blocks <- lapply(seq_along(domains), function(i) {
    Xi <- as.matrix(domains[[i]]$x)
    Si <- scores[domain_starts[i]:domain_offsets[i], , drop = FALSE]
    .msa_fit_loadings(Xi, Si)
  })
  loadings <- do.call(rbind, loading_blocks)

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "multiscale_manifold_align",
    extras = list(
      backend = ctrl$backend,
      operator = "diffusion_symmetric",
      levels = dictionary_info$levels,
      wavelet_levels = dictionary_info$wavelet_levels,
      level_dims = dictionary_info$level_dims,
      multiscale_times = dictionary_info$scale_times,
      filter_bank = list(
        scaling = dictionary_info$scaling_weights,
        wavelet = dictionary_info$wavelet_weights
      ),
      reduced_eigenvalues = eigvals,
      generalized_eigenvalues = solved$eigenvalues,
      correspondences = corr_df,
      domain_names = domain_names,
      scale_weights = dictionary_info$scale_weights,
      selected_scales = dictionary_info$selected_scales,
      selected_scale_weights = dictionary_info$selected_scale_weights,
      scale_energy = dictionary_info$energy_scores,
      scale_supervision = dictionary_info$supervision_scores,
      dictionary_dims = ncol(dictionary),
      correspondence_margin = objective,
      control = ctrl
    )
  )
}

#' @rdname multiscale_manifold_align
#' @method multiscale_manifold_align hyperdesign
#' @param y Optional unquoted label column used to build cross-domain correspondences.
#' @param correspondences Optional data.frame with columns domain_i,index_i,domain_j,index_j[,weight].
#' @param preproc Optional preprocessing function/preprocessor applied per domain.
#' @param ncomp Number of target embedding components.
#' @param control Control settings from [multiscale_manifold_align_control()].
#' @param verbose Logical progress toggle.
#' @export
multiscale_manifold_align.hyperdesign <- function(data,
                                                  y = NULL,
                                                  correspondences = NULL,
                                                  preproc = multivarious::center(),
                                                  ncomp = 20L,
                                                  control = multiscale_manifold_align_control(),
                                                  verbose = TRUE,
                                                  ...) {
  ctrl <- resolve_multiscale_manifold_align_control(control)
  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  resolved <- resolve_hyperdesign(data)
  prepped <- .msa_preprocess_domains(resolved$domains, preproc = preproc)
  domains <- prepped$domains
  n_domains <- length(domains)
  domain_sizes <- vapply(domains, function(d) nrow(d$x), integer(1))
  n_total <- sum(domain_sizes)
  feature_blocks <- feature_block_indices(domains)
  proc <- if (all(vapply(prepped$preproc, is.null, logical(1)))) {
    NULL
  } else {
    multivarious::concat_pre_processors(prepped$preproc, feature_blocks)
  }

  if (verbose || isTRUE(ctrl$verbose)) {
    message("multiscale_manifold_align: domains=", n_domains, ", nodes=", n_total)
  }

  within <- lapply(domains, function(dom) .msa_build_within_graph(dom$x, knn = ctrl$knn, sigma = ctrl$sigma))

  y_quo <- if (missing(y)) rlang::quo(NULL) else rlang::enquo(y)
  label_list <- .msa_extract_labels(domains, y_quo)

  corr_df <- .msa_validate_correspondences(correspondences, domain_sizes = domain_sizes, default_weight = 1)
  if (!nrow(corr_df)) {
    if (!is.null(label_list)) {
      corr_df <- .msa_build_label_correspondences(
        labels_list = label_list,
        max_pairs_per_label = ctrl$max_pairs_per_label,
        weight = 1,
        seed = ctrl$seed
      )
    } else {
      corr_df <- .msa_build_mnn_correspondences(
        domains = domains,
        mnn_k = ctrl$mnn_k,
        max_pairs_per_pair = ctrl$max_unsup_pairs_per_pair,
        weight = 1
      )
    }
  }

  fit_for_ctrl <- function(ctrl_i) {
    .msa_fit_core(
      domains = domains,
      within = within,
      corr_df = .msa_rescale_correspondences(corr_df, ctrl_i$cross_edge_weight),
      domain_names = resolved$domain_names,
      ncomp = ncomp,
      ctrl = ctrl_i,
      proc = proc
    )
  }

  if (!isTRUE(ctrl$tune) || !nrow(corr_df)) {
    return(fit_for_ctrl(ctrl))
  }

  grid <- .msa_tuning_grid(ctrl)
  best_fit <- NULL
  best_score <- -Inf
  best_ctrl <- ctrl

  for (i in seq_len(nrow(grid))) {
    ctrl_i <- resolve_multiscale_manifold_align_control(utils::modifyList(
      ctrl,
      list(
        cross_edge_weight = grid$cross_edge_weight[i],
        rank_per_domain = grid$rank_per_domain[i],
        max_levels = grid$max_levels[i],
        tune = FALSE
      )
    ))

    fit_i <- tryCatch(fit_for_ctrl(ctrl_i), error = function(e) NULL)
    if (is.null(fit_i)) next
    score_i <- fit_i$correspondence_margin
    if (!is.finite(score_i)) next
    if (score_i > best_score) {
      best_score <- score_i
      best_fit <- fit_i
      best_ctrl <- ctrl_i
    }
  }

  if (is.null(best_fit)) {
    best_fit <- fit_for_ctrl(resolve_multiscale_manifold_align_control(utils::modifyList(ctrl, list(tune = FALSE))))
    best_score <- best_fit$correspondence_margin
    best_ctrl <- best_fit$control
  }

  best_fit$tuning <- list(
    enabled = TRUE,
    objective = best_score,
    candidates_evaluated = nrow(grid),
    selected = list(
      cross_edge_weight = best_ctrl$cross_edge_weight,
      rank_per_domain = best_ctrl$rank_per_domain,
      max_levels = best_ctrl$max_levels
    )
  )
  best_fit
}

#' @rdname multiscale_manifold_align
#' @method multiscale_manifold_align list
#' @export
multiscale_manifold_align.list <- function(data,
                                           y = NULL,
                                           correspondences = NULL,
                                           preproc = multivarious::center(),
                                           ncomp = 20L,
                                           control = multiscale_manifold_align_control(),
                                           verbose = TRUE,
                                           ...) {
  if (!length(data) || !all(vapply(data, is.matrix, logical(1)))) {
    stop("For list input, `data` must be a non-empty list of matrices.", call. = FALSE)
  }
  hd <- as_hyperdesign(data)
  if (missing(y)) {
    return(multiscale_manifold_align.hyperdesign(
      hd,
      correspondences = correspondences,
      preproc = preproc,
      ncomp = ncomp,
      control = control,
      verbose = verbose,
      ...
    ))
  }
  y_quo <- rlang::enquo(y)
  rlang::inject(multiscale_manifold_align.hyperdesign(
    hd,
    y = !!y_quo,
    correspondences = correspondences,
    preproc = preproc,
    ncomp = ncomp,
    control = control,
    verbose = verbose,
    ...
  ))
}

#' @rdname multiscale_manifold_align
#' @method multiscale_manifold_align default
#' @export
multiscale_manifold_align.default <- function(data, ...) {
  stop("No applicable method for multiscale_manifold_align().", call. = FALSE)
}
