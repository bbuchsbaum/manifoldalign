#' Control settings for serial-correlation handling in SSMA graph construction
#'
#' @param enabled Logical; enable serial-correlation adjustments when
#'   `serial_index` is provided.
#' @param row_whiten One of `"none"` or `"ar1"`. Controls optional row-wise
#'   whitening applied before within-domain kNN graph construction.
#' @param ar1_shrink Scalar in `[0, 1)` used to shrink estimated AR(1)
#'   coefficient toward zero.
#' @param ar1_clip Maximum absolute AR(1) coefficient allowed for stability.
#' @param lag_mode One of `"none"`, `"hard"`, `"exponential"`, or `"gaussian"`.
#'   Controls how short-lag neighbor edges are downweighted.
#' @param lag_exclusion Non-negative integer lag exclusion radius used when
#'   `lag_mode = "hard"`.
#' @param lag_scale Positive scalar lag scale used for soft lag kernels.
#' @param preserve_degree Logical; if `TRUE`, rescale filtered adjacency to keep
#'   row-degree magnitudes near pre-filter values.
#'
#' @return A list of class `ssma_serial_control`.
#' @export
ssma_serial_control <- function(
  enabled = FALSE,
  row_whiten = c("none", "ar1"),
  ar1_shrink = 0,
  ar1_clip = 0.98,
  lag_mode = c("none", "hard", "exponential", "gaussian"),
  lag_exclusion = 0L,
  lag_scale = 2,
  preserve_degree = TRUE
) {
  row_whiten <- match.arg(row_whiten)
  lag_mode <- match.arg(lag_mode)

  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
    stop("`enabled` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(ar1_shrink) || length(ar1_shrink) != 1L || is.na(ar1_shrink) || ar1_shrink < 0 || ar1_shrink >= 1) {
    stop("`ar1_shrink` must be in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(ar1_clip) || length(ar1_clip) != 1L || is.na(ar1_clip) || ar1_clip <= 0 || ar1_clip >= 1) {
    stop("`ar1_clip` must be in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(lag_exclusion) || length(lag_exclusion) != 1L || is.na(lag_exclusion) || lag_exclusion < 0) {
    stop("`lag_exclusion` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(lag_scale) || length(lag_scale) != 1L || is.na(lag_scale) || lag_scale <= 0) {
    stop("`lag_scale` must be a positive scalar.", call. = FALSE)
  }
  if (!is.logical(preserve_degree) || length(preserve_degree) != 1L || is.na(preserve_degree)) {
    stop("`preserve_degree` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      enabled = enabled,
      row_whiten = row_whiten,
      ar1_shrink = as.numeric(ar1_shrink),
      ar1_clip = as.numeric(ar1_clip),
      lag_mode = lag_mode,
      lag_exclusion = as.integer(lag_exclusion),
      lag_scale = as.numeric(lag_scale),
      preserve_degree = preserve_degree
    ),
    class = "ssma_serial_control"
  )
}

resolve_ssma_serial_control <- function(serial_control = NULL) {
  defaults <- unclass(ssma_serial_control())

  if (missing(serial_control) || is.null(serial_control)) {
    return(do.call(ssma_serial_control, defaults))
  }
  if (inherits(serial_control, "ssma_serial_control")) {
    merged <- utils::modifyList(defaults, unclass(serial_control))
    return(do.call(ssma_serial_control, merged))
  }
  if (!is.list(serial_control)) {
    stop("`serial_control` must be NULL, a named list, or ssma_serial_control().", call. = FALSE)
  }
  if (is.null(names(serial_control)) || any(names(serial_control) == "")) {
    stop("`serial_control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(serial_control), names(defaults))
  if (length(unknown)) {
    stop("Unknown serial control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  merged <- utils::modifyList(defaults, serial_control)
  do.call(ssma_serial_control, merged)
}

#' Control settings for `ssma_align()`
#'
#' @param knn Integer k used to build within-domain kNN graphs.
#' @param sigma Optional positive heat-kernel bandwidth for within-domain
#'   graph weights. If `NULL`, per-domain auto-tuning via [choose_sigma()] is
#'   used.
#' @param corr_edge_weight Non-negative scalar correspondence edge weight.
#' @param max_pairs_per_label Integer cap per shared label/domain pair when
#'   generating correspondences from labels.
#' @param mnn_k Integer k used for unsupervised MNN correspondence discovery.
#' @param max_unsup_pairs_per_pair Integer cap per domain pair for unsupervised
#'   correspondences.
#' @param rank_per_domain Integer rank cap for per-domain compression basis.
#' @param svd_tol Relative singular-value tolerance used in compression.
#' @param eig_tol Non-negative threshold for filtering trivial generalized
#'   eigenvalues.
#' @param regularization Positive ridge regularizer added to reduced `B` matrix
#'   before whitening.
#' @param solver One of `"auto"`, `"reduced"`, or `"operator"`.
#' @param operator_tol Positive convergence tolerance for operator eigensolve.
#' @param operator_maxiter Maximum ARPACK iterations for operator eigensolve.
#' @param operator_power_iter Number of power-iterations used to estimate
#'   operator spectral scale.
#' @param seed Integer random seed used in correspondence sampling.
#' @param verbose Logical progress toggle.
#' @param serial_control Serial-correlation handling settings from
#'   [ssma_serial_control()].
#'
#' @return A list of class `ssma_align_control`.
#' @export
ssma_align_control <- function(
  knn = 15L,
  sigma = NULL,
  corr_edge_weight = 1,
  max_pairs_per_label = 100L,
  mnn_k = 5L,
  max_unsup_pairs_per_pair = 200L,
  rank_per_domain = 128L,
  svd_tol = 1e-6,
  eig_tol = 1e-8,
  regularization = 1e-6,
  solver = c("auto", "reduced", "operator"),
  operator_tol = 1e-6,
  operator_maxiter = 2000L,
  operator_power_iter = 25L,
  seed = 1L,
  verbose = FALSE,
  serial_control = ssma_serial_control()
) {
  solver <- match.arg(solver)

  if (!is.numeric(knn) || length(knn) != 1L || is.na(knn) || knn < 1) {
    stop("`knn` must be >= 1.", call. = FALSE)
  }
  if (!is.null(sigma) && (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0)) {
    stop("`sigma` must be NULL or a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(corr_edge_weight) || length(corr_edge_weight) != 1L || is.na(corr_edge_weight) || corr_edge_weight < 0) {
    stop("`corr_edge_weight` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(max_pairs_per_label) || length(max_pairs_per_label) != 1L || is.na(max_pairs_per_label) || max_pairs_per_label < 1) {
    stop("`max_pairs_per_label` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(mnn_k) || length(mnn_k) != 1L || is.na(mnn_k) || mnn_k < 1) {
    stop("`mnn_k` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(max_unsup_pairs_per_pair) || length(max_unsup_pairs_per_pair) != 1L || is.na(max_unsup_pairs_per_pair) || max_unsup_pairs_per_pair < 1) {
    stop("`max_unsup_pairs_per_pair` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(rank_per_domain) || length(rank_per_domain) != 1L || is.na(rank_per_domain) || rank_per_domain < 1) {
    stop("`rank_per_domain` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(svd_tol) || length(svd_tol) != 1L || is.na(svd_tol) || svd_tol <= 0) {
    stop("`svd_tol` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(eig_tol) || length(eig_tol) != 1L || is.na(eig_tol) || eig_tol < 0) {
    stop("`eig_tol` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(regularization) || length(regularization) != 1L || is.na(regularization) || regularization <= 0) {
    stop("`regularization` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(operator_tol) || length(operator_tol) != 1L || is.na(operator_tol) || operator_tol <= 0) {
    stop("`operator_tol` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(operator_maxiter) || length(operator_maxiter) != 1L || is.na(operator_maxiter) || operator_maxiter < 1) {
    stop("`operator_maxiter` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(operator_power_iter) || length(operator_power_iter) != 1L || is.na(operator_power_iter) || operator_power_iter < 1) {
    stop("`operator_power_iter` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  serial_control <- resolve_ssma_serial_control(serial_control)

  structure(
    list(
      knn = as.integer(knn),
      sigma = if (is.null(sigma)) NULL else as.numeric(sigma),
      corr_edge_weight = as.numeric(corr_edge_weight),
      max_pairs_per_label = as.integer(max_pairs_per_label),
      mnn_k = as.integer(mnn_k),
      max_unsup_pairs_per_pair = as.integer(max_unsup_pairs_per_pair),
      rank_per_domain = as.integer(rank_per_domain),
      svd_tol = as.numeric(svd_tol),
      eig_tol = as.numeric(eig_tol),
      regularization = as.numeric(regularization),
      solver = solver,
      operator_tol = as.numeric(operator_tol),
      operator_maxiter = as.integer(operator_maxiter),
      operator_power_iter = as.integer(operator_power_iter),
      seed = as.integer(seed),
      verbose = verbose,
      serial_control = serial_control
    ),
    class = "ssma_align_control"
  )
}

resolve_ssma_align_control <- function(control = NULL) {
  defaults <- unclass(ssma_align_control())

  if (missing(control) || is.null(control)) {
    return(do.call(ssma_align_control, defaults))
  }
  if (inherits(control, "ssma_align_control")) {
    merged <- utils::modifyList(defaults, unclass(control))
    return(do.call(ssma_align_control, merged))
  }
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or ssma_align_control().", call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  if ("serial_control" %in% names(control)) {
    control$serial_control <- resolve_ssma_serial_control(control$serial_control)
  }

  merged <- utils::modifyList(defaults, control)
  do.call(ssma_align_control, merged)
}

.ssma_extract_labels <- function(domains, y_quo) {
  if (rlang::quo_is_null(y_quo)) {
    return(NULL)
  }

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

.ssma_build_label_correspondences <- function(labels_list, max_pairs_per_label, weight, seed) {
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

.ssma_build_mnn_correspondences <- function(domains, mnn_k, max_pairs_per_pair, weight) {
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

.ssma_validate_correspondences <- function(corr, domain_sizes, default_weight = 1) {
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

.ssma_build_joint_laplacian <- function(within_graphs, correspondences, domain_sizes) {
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

.ssma_resolve_serial_indices <- function(domains, serial_index = NULL) {
  n_domains <- length(domains)

  if (missing(serial_index) || is.null(serial_index)) {
    return(replicate(n_domains, NULL, simplify = FALSE))
  }

  if (is.character(serial_index) && length(serial_index) == 1L) {
    nm <- serial_index
    out <- lapply(seq_len(n_domains), function(i) {
      di <- domains[[i]]$design
      if (is.null(di) || !is.data.frame(di) || !(nm %in% names(di))) {
        stop("`serial_index` column `", nm, "` not found in domain ", i, ".", call. = FALSE)
      }
      as.numeric(di[[nm]])
    })
  } else if (is.list(serial_index)) {
    if (length(serial_index) != n_domains) {
      stop("When `serial_index` is a list, its length must match number of domains.", call. = FALSE)
    }
    out <- lapply(serial_index, as.numeric)
  } else if (is.numeric(serial_index) && n_domains == 1L) {
    out <- list(as.numeric(serial_index))
  } else {
    stop("`serial_index` must be NULL, a single design-column name, or a list with one index vector per domain.", call. = FALSE)
  }

  for (i in seq_len(n_domains)) {
    if (is.null(out[[i]])) next
    if (length(out[[i]]) != nrow(domains[[i]]$x)) {
      stop("serial_index length mismatch in domain ", i, ".", call. = FALSE)
    }
    if (any(!is.finite(out[[i]]))) {
      stop("serial_index for domain ", i, " contains non-finite values.", call. = FALSE)
    }
  }

  out
}

.ssma_ar1_whiten <- function(X, serial_idx, shrink = 0, clip = 0.98, min_n = 6L) {
  X <- as.matrix(X)
  n <- nrow(X)

  if (n < min_n) {
    return(list(X = X, phi = NA_real_, applied = FALSE))
  }

  ord <- order(serial_idx)
  Xo <- X[ord, , drop = FALSE]
  Xo <- scale(Xo, center = TRUE, scale = FALSE)

  prev <- Xo[-n, , drop = FALSE]
  nextv <- Xo[-1L, , drop = FALSE]
  denom <- colSums(prev * prev)
  numer <- colSums(prev * nextv)

  phi_j <- numer / pmax(denom, 1e-12)
  phi_j <- phi_j[is.finite(phi_j)]
  phi <- if (length(phi_j)) stats::median(phi_j) else 0

  phi <- (1 - shrink) * phi
  phi <- max(min(phi, clip), -clip)

  R <- Xo
  if (n > 1L) {
    R[1L, ] <- sqrt(max(1 - phi^2, 1e-8)) * Xo[1L, ]
    R[-1L, ] <- Xo[-1L, , drop = FALSE] - phi * Xo[-n, , drop = FALSE]
  }

  out <- X
  out[ord, ] <- R

  list(X = out, phi = phi, applied = TRUE)
}

.ssma_apply_serial_lag_filter <- function(A, serial_idx, serial_ctrl) {
  A0 <- methods::as(Matrix::forceSymmetric(A, uplo = "U"), "CsparseMatrix")
  edge_df <- Matrix::summary(A0)

  if (!nrow(edge_df) || serial_ctrl$lag_mode == "none") {
    return(list(
      A = A0,
      serial_info = list(masked_edges = 0L, total_edges = as.integer(nrow(edge_df)), masked_fraction = 0)
    ))
  }

  lag <- abs(serial_idx[edge_df$i] - serial_idx[edge_df$j])
  lag <- as.numeric(lag)

  if (serial_ctrl$lag_mode == "hard") {
    edge_w <- as.numeric(lag > serial_ctrl$lag_exclusion)
  } else if (serial_ctrl$lag_mode == "exponential") {
    edge_w <- 1 - exp(-lag / serial_ctrl$lag_scale)
  } else {
    edge_w <- 1 - exp(-(lag^2) / (2 * serial_ctrl$lag_scale^2))
  }

  edge_w <- pmax(pmin(edge_w, 1), 0)

  A1 <- Matrix::sparseMatrix(
    i = edge_df$i,
    j = edge_df$j,
    x = edge_df$x * edge_w,
    dims = dim(A0),
    symmetric = FALSE
  )

  A1 <- Matrix::forceSymmetric(Matrix::drop0(A1), uplo = "U")
  Matrix::diag(A1) <- 0

  if (isTRUE(serial_ctrl$preserve_degree)) {
    deg0 <- as.numeric(Matrix::rowSums(A0))
    deg1 <- as.numeric(Matrix::rowSums(A1))
    s <- sqrt(pmax(deg0, 0) / pmax(deg1, 1e-12))
    s[!is.finite(s)] <- 1
    S <- Matrix::Diagonal(x = s)
    A1 <- Matrix::drop0(S %*% A1 %*% S)
    A1 <- Matrix::forceSymmetric(A1, uplo = "U")
    Matrix::diag(A1) <- 0
  }

  masked_edges <- sum(edge_w < 1)
  total_edges <- nrow(edge_df)

  list(
    A = A1,
    serial_info = list(
      masked_edges = as.integer(masked_edges),
      total_edges = as.integer(total_edges),
      masked_fraction = if (total_edges > 0) as.numeric(masked_edges / total_edges) else 0
    )
  )
}

.ssma_build_within_graph <- function(X, knn, sigma = NULL, serial_idx = NULL, serial_ctrl = ssma_serial_control()) {
  X <- as.matrix(X)
  n <- nrow(X)

  serial_info <- list(
    enabled = isTRUE(serial_ctrl$enabled) && !is.null(serial_idx),
    row_whiten = serial_ctrl$row_whiten,
    ar1_phi = NA_real_,
    lag_mode = serial_ctrl$lag_mode,
    lag_exclusion = serial_ctrl$lag_exclusion,
    lag_scale = serial_ctrl$lag_scale,
    masked_edges = 0L,
    total_edges = 0L,
    masked_fraction = 0
  )

  if (n <= 1L) {
    A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    return(list(A = A, degree = numeric(n), serial = serial_info))
  }

  X_graph <- X
  if (serial_info$enabled && serial_ctrl$row_whiten == "ar1") {
    white <- .ssma_ar1_whiten(
      X_graph,
      serial_idx = serial_idx,
      shrink = serial_ctrl$ar1_shrink,
      clip = serial_ctrl$ar1_clip
    )
    X_graph <- white$X
    serial_info$ar1_phi <- white$phi
  }

  k_use <- min(as.integer(knn), n - 1L)
  sigma_use <- if (is.null(sigma)) choose_sigma(X_graph) else sigma

  gw <- safe_compute(
    neighborweights::graph_weights(
      X_graph,
      k = k_use,
      weight_mode = "heat",
      sigma = sigma_use,
      neighbor_mode = "knn"
    ),
    "Within-domain graph construction failed in ssma_align()."
  )

  A <- neighborweights::adjacency(gw)
  A <- (A + Matrix::t(A)) / 2
  Matrix::diag(A) <- 0
  A <- Matrix::drop0(methods::as(A, "CsparseMatrix"))

  if (serial_info$enabled && serial_ctrl$lag_mode != "none") {
    filtered <- .ssma_apply_serial_lag_filter(A, serial_idx = serial_idx, serial_ctrl = serial_ctrl)
    A <- filtered$A
    serial_info$masked_edges <- filtered$serial_info$masked_edges
    serial_info$total_edges <- filtered$serial_info$total_edges
    serial_info$masked_fraction <- filtered$serial_info$masked_fraction
  }

  list(A = A, degree = as.numeric(Matrix::rowSums(A)), serial = serial_info)
}

.ssma_domain_basis <- function(X, degree, rank_per_domain, svd_tol) {
  X <- as.matrix(X)
  dhalf <- sqrt(pmax(as.numeric(degree), 1e-12))
  Xw <- X * dhalf

  k <- min(as.integer(rank_per_domain), nrow(Xw), ncol(Xw))
  if (k < 1L) {
    return(list(P = matrix(0, ncol(X), 1L), T = matrix(0, nrow(X), 1L), s = 0))
  }

  sv <- if (requireNamespace("irlba", quietly = TRUE) && k < min(dim(Xw))) {
    irlba::irlba(Xw, nv = k, nu = 0)
  } else {
    s <- svd(Xw, nu = 0, nv = k)
    list(v = s$v[, seq_len(k), drop = FALSE], d = s$d[seq_len(k)])
  }

  dmax <- if (length(sv$d)) max(sv$d) else 0
  keep <- which(is.finite(sv$d) & sv$d >= svd_tol * max(dmax, 1e-12))
  if (!length(keep)) keep <- 1L

  P <- as.matrix(sv$v[, keep, drop = FALSE])
  T <- X %*% P

  list(P = P, T = T, s = as.numeric(sv$d[keep]))
}

.ssma_build_reduced_system <- function(domains, within_graphs, L_joint, rank_per_domain, svd_tol,
                                       assemble_dense = TRUE) {
  n_domains <- length(domains)
  domain_sizes <- vapply(domains, function(d) nrow(d$x), integer(1))
  degree_list <- lapply(within_graphs, function(g) pmax(as.numeric(g$degree), 1e-12))

  basis <- vector("list", n_domains)
  ranks <- integer(n_domains)

  for (i in seq_len(n_domains)) {
    bi <- .ssma_domain_basis(
      X = domains[[i]]$x,
      degree = within_graphs[[i]]$degree,
      rank_per_domain = rank_per_domain,
      svd_tol = svd_tol
    )
    basis[[i]] <- bi
    ranks[i] <- ncol(bi$P)
  }

  rank_offsets <- c(0L, cumsum(ranks))
  sample_offsets <- c(0L, cumsum(domain_sizes))
  r_total <- sum(ranks)

  B_blocks <- lapply(seq_len(n_domains), function(i) {
    Ti <- basis[[i]]$T
    di <- degree_list[[i]]
    crossprod(Ti, Ti * di)
  })

  A_r <- NULL
  B_r <- NULL
  if (isTRUE(assemble_dense)) {
    A_r <- matrix(0, nrow = r_total, ncol = r_total)
    B_r <- matrix(0, nrow = r_total, ncol = r_total)

    for (i in seq_len(n_domains)) {
      ri <- (rank_offsets[i] + 1L):rank_offsets[i + 1L]
      ni <- (sample_offsets[i] + 1L):sample_offsets[i + 1L]
      Ti <- basis[[i]]$T

      B_r[ri, ri] <- as.matrix(B_blocks[[i]])

      for (j in i:n_domains) {
        rj <- (rank_offsets[j] + 1L):rank_offsets[j + 1L]
        nj <- (sample_offsets[j] + 1L):sample_offsets[j + 1L]

        Lij <- L_joint[ni, nj, drop = FALSE]
        if (Matrix::nnzero(Lij) == 0) next

        Tj <- basis[[j]]$T
        block <- crossprod(Ti, Lij %*% Tj)
        A_r[ri, rj] <- A_r[ri, rj] + as.matrix(block)
        if (i != j) {
          A_r[rj, ri] <- A_r[rj, ri] + t(as.matrix(block))
        }
      }
    }
  }

  list(
    A_r = if (isTRUE(assemble_dense)) 0.5 * (A_r + t(A_r)) else NULL,
    B_r = if (isTRUE(assemble_dense)) 0.5 * (B_r + t(B_r)) else NULL,
    B_blocks = B_blocks,
    basis = basis,
    ranks = ranks,
    domain_sizes = domain_sizes,
    degree_list = degree_list
  )
}

.ssma_solve_reduced_geigen <- function(A_r, B_r, ncomp, eig_tol, regularization) {
  r <- nrow(A_r)
  if (r < 1L) {
    stop("Reduced system has zero rank.", call. = FALSE)
  }

  ridge <- max(regularization, 1e-10 * mean(diag(B_r) + 1e-12))
  B_reg <- B_r + diag(ridge, r)
  B_reg <- 0.5 * (B_reg + t(B_reg))
  A_sym <- 0.5 * (A_r + t(A_r))

  R <- tryCatch(chol(B_reg), error = function(e) NULL)
  if (is.null(R)) {
    ridge2 <- max(1e-8, 100 * ridge)
    B_reg <- B_reg + diag(ridge2, r)
    R <- tryCatch(chol(B_reg), error = function(e) NULL)
    if (is.null(R)) {
      stop("Reduced generalized eigenproblem failed: B matrix is not SPD after regularization.", call. = FALSE)
    }
    ridge <- ridge + ridge2
  }

  W <- backsolve(R, diag(r))
  S <- crossprod(W, A_sym %*% W)
  S <- 0.5 * (S + t(S))

  k_req <- min(max(as.integer(ncomp) + 5L, as.integer(ncomp)), r)

  if (requireNamespace("RSpectra", quietly = TRUE) && k_req < (r - 1L)) {
    eig <- RSpectra::eigs_sym(S, k = k_req, which = "SA")
    vals <- as.numeric(eig$values)
    vecs <- as.matrix(eig$vectors)
  } else {
    eig <- eigen(S, symmetric = TRUE)
    vals <- as.numeric(eig$values[seq_len(k_req)])
    vecs <- as.matrix(eig$vectors[, seq_len(k_req), drop = FALSE])
  }

  ord <- order(vals, decreasing = FALSE)
  vals <- vals[ord]
  vecs <- vecs[, ord, drop = FALSE]

  keep <- which(is.finite(vals) & vals > eig_tol)
  if (!length(keep)) {
    keep <- seq_len(min(as.integer(ncomp), length(vals)))
  }

  k_out <- min(as.integer(ncomp), length(keep))
  keep <- keep[seq_len(k_out)]

  alpha <- W %*% vecs[, keep, drop = FALSE]

  list(
    alpha = alpha,
    eigenvalues = vals[keep],
    eigenvalues_all = vals,
    ridge = ridge
  )
}

.ssma_as_matrix <- function(x, nrow) {
  if (is.matrix(x)) {
    return(list(mat = x, vec = FALSE))
  }
  list(mat = matrix(x, nrow = nrow), vec = TRUE)
}

.ssma_split_blocks <- function(X, sizes) {
  offs <- c(0L, cumsum(sizes))
  lapply(seq_along(sizes), function(i) {
    X[(offs[i] + 1L):offs[i + 1L], , drop = FALSE]
  })
}

.ssma_bind_blocks <- function(blocks, as_vector = FALSE) {
  out <- do.call(rbind, blocks)
  if (as_vector) drop(out) else out
}

.ssma_lambda_max_power <- function(op, n, n_iter = 30L, tol = 1e-6, seed = 1L) {
  set.seed(seed)
  v <- rnorm(n)
  v <- v / sqrt(sum(v * v))
  lambda_prev <- 0

  for (it in seq_len(n_iter)) {
    w <- op(v)
    lambda <- sum(v * w)
    norm_w <- sqrt(sum(w * w))
    if (norm_w <= 0 || !is.finite(norm_w)) break
    v <- w / norm_w
    if (abs(lambda - lambda_prev) <= tol * (abs(lambda_prev) + 1e-12)) break
    lambda_prev <- lambda
  }

  lambda_prev
}

.ssma_solve_operator_geigen <- function(reduced,
                                        L_joint,
                                        ncomp,
                                        eig_tol,
                                        regularization,
                                        operator_tol,
                                        operator_maxiter,
                                        operator_power_iter,
                                        seed) {
  ranks <- reduced$ranks
  basis <- reduced$basis
  B_blocks <- reduced$B_blocks
  domain_sizes <- reduced$domain_sizes
  r_total <- sum(ranks)

  if (r_total < 1L) {
    stop("Reduced system has zero rank.", call. = FALSE)
  }

  B_chol <- lapply(B_blocks, function(Bi) {
    Bi <- as.matrix(Bi)
    ridge_i <- max(regularization, 1e-10 * mean(diag(Bi) + 1e-12))
    Ri <- tryCatch(chol(Bi + diag(ridge_i, nrow(Bi))), error = function(e) NULL)
    if (is.null(Ri)) {
      ridge_i <- max(1e-6, 100 * ridge_i)
      Ri <- tryCatch(chol(Bi + diag(ridge_i, nrow(Bi))), error = function(e) NULL)
    }
    if (is.null(Ri)) {
      stop("Operator generalized eigenproblem failed: block B is not SPD after regularization.", call. = FALSE)
    }
    Ri
  })

  apply_inv_half <- function(x) {
    xm <- .ssma_as_matrix(x, r_total)
    xb <- .ssma_split_blocks(xm$mat, ranks)
    yb <- lapply(seq_along(xb), function(i) backsolve(B_chol[[i]], xb[[i]], upper.tri = TRUE))
    .ssma_bind_blocks(yb, as_vector = xm$vec)
  }

  apply_inv_half_t <- function(x) {
    xm <- .ssma_as_matrix(x, r_total)
    xb <- .ssma_split_blocks(xm$mat, ranks)
    yb <- lapply(seq_along(xb), function(i) forwardsolve(t(B_chol[[i]]), xb[[i]], upper.tri = FALSE))
    .ssma_bind_blocks(yb, as_vector = xm$vec)
  }

  A_op <- function(x, args = NULL) {
    xm <- .ssma_as_matrix(x, r_total)
    xb <- .ssma_split_blocks(xm$mat, ranks)

    u_blocks <- lapply(seq_along(xb), function(i) basis[[i]]$T %*% xb[[i]])
    u <- do.call(rbind, u_blocks)
    v <- as.matrix(L_joint %*% u)
    v_blocks <- .ssma_split_blocks(v, domain_sizes)
    y_blocks <- lapply(seq_along(v_blocks), function(i) crossprod(basis[[i]]$T, v_blocks[[i]]))

    .ssma_bind_blocks(y_blocks, as_vector = xm$vec)
  }

  S_op <- function(x, args = NULL) {
    y <- apply_inv_half(x)
    a <- A_op(y)
    apply_inv_half_t(a)
  }

  lam_max <- .ssma_lambda_max_power(
    op = function(v) as.numeric(S_op(v)),
    n = r_total,
    n_iter = as.integer(operator_power_iter),
    tol = operator_tol,
    seed = seed
  )
  if (!is.finite(lam_max) || lam_max <= 0) lam_max <- 1
  c_shift <- 1.05 * lam_max + 1e-8

  shifted_op <- function(x, args = NULL) {
    xm <- .ssma_as_matrix(x, r_total)
    y <- xm$mat - (1 / c_shift) * as.matrix(S_op(xm$mat))
    if (xm$vec) drop(y) else y
  }

  k_req <- min(max(as.integer(ncomp) + 5L, as.integer(ncomp)), max(1L, r_total - 1L))
  if (k_req < 1L) {
    stop("Insufficient reduced rank for operator eigendecomposition.", call. = FALSE)
  }

  if (requireNamespace("RSpectra", quietly = TRUE)) {
    eig <- RSpectra::eigs_sym(
      shifted_op,
      k = k_req,
      n = r_total,
      which = "LA",
      opts = list(maxitr = as.integer(operator_maxiter), tol = as.numeric(operator_tol))
    )
    vals <- c_shift * (1 - as.numeric(eig$values))
    vecs <- as.matrix(eig$vectors)
  } else {
    stop("solver='operator' requires RSpectra.", call. = FALSE)
  }

  ord <- order(vals, decreasing = FALSE)
  vals <- vals[ord]
  vecs <- vecs[, ord, drop = FALSE]

  keep <- which(is.finite(vals) & vals > eig_tol)
  if (!length(keep)) {
    keep <- seq_len(min(as.integer(ncomp), length(vals)))
  }
  k_out <- min(as.integer(ncomp), length(keep))
  keep <- keep[seq_len(k_out)]

  q <- vecs[, keep, drop = FALSE]
  alpha <- apply_inv_half(q)
  if (is.null(dim(alpha))) alpha <- matrix(alpha, ncol = k_out)

  list(
    alpha = as.matrix(alpha),
    eigenvalues = vals[keep],
    eigenvalues_all = vals,
    ridge = regularization
  )
}

#' Semi-Supervised Manifold Alignment (SSMA)
#'
#' Baseline SSMA implementation based on Wang et al. feature-level manifold
#' alignment (`Z^T L Z f = \lambda Z^T D Z f`) solved in a reduced-rank space.
#'
#' @param data A `hyperdesign` object or list of matrices.
#' @param ... Additional arguments passed to methods.
#' @export
ssma_align <- function(data, ...) {
  UseMethod("ssma_align")
}

#' @rdname ssma_align
#' @method ssma_align hyperdesign
#' @param y Optional unquoted label column used for label-based correspondences
#'   when explicit `correspondences` are not provided.
#' @param correspondences Optional data.frame with columns
#'   `domain_i,index_i,domain_j,index_j[,weight]`.
#' @param serial_index Optional serial/time index specification used for
#'   within-domain graph decontamination. Can be:
#'   1) `NULL`, 2) a single design-column name, or 3) a list with one numeric
#'   index vector per domain.
#' @param preproc Optional preprocessing function/preprocessor applied per
#'   domain. Default: `multivarious::center()`.
#' @param ncomp Number of aligned components.
#' @param control Control settings from [ssma_align_control()].
#' @export
ssma_align.hyperdesign <- function(
  data,
  y = NULL,
  correspondences = NULL,
  serial_index = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = ssma_align_control(),
  ...
) {
  ctrl <- resolve_ssma_align_control(control)
  verbose <- isTRUE(ctrl$verbose)

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }

  resolved <- resolve_hyperdesign(data)
  domains_in <- resolved$domains
  domain_names <- resolved$domain_names
  n_domains <- resolved$n_domains

  if (n_domains < 2L) {
    stop("ssma_align requires at least two domains.", call. = FALSE)
  }

  pre <- .apply_preproc_domains_spectral_mnn(domains_in, preproc = preproc)
  domains <- pre$pdata
  proc <- pre$proc
  names(domains) <- domain_names

  serial_idx <- .ssma_resolve_serial_indices(domains, serial_index = serial_index)

  if (verbose) {
    message("ssma_align: building within-domain graphs for ", n_domains, " domains")
  }

  within <- lapply(seq_len(n_domains), function(i) {
    .ssma_build_within_graph(
      X = domains[[i]]$x,
      knn = ctrl$knn,
      sigma = ctrl$sigma,
      serial_idx = serial_idx[[i]],
      serial_ctrl = ctrl$serial_control
    )
  })

  domain_sizes <- vapply(domains, function(d) nrow(d$x), integer(1))

  y_quo <- if (missing(y)) rlang::quo(NULL) else rlang::enquo(y)
  label_list <- .ssma_extract_labels(domains, y_quo)

  corr_df <- .ssma_validate_correspondences(
    correspondences,
    domain_sizes = domain_sizes,
    default_weight = ctrl$corr_edge_weight
  )

  if (!nrow(corr_df)) {
    if (!is.null(label_list)) {
      corr_df <- .ssma_build_label_correspondences(
        labels_list = label_list,
        max_pairs_per_label = ctrl$max_pairs_per_label,
        weight = ctrl$corr_edge_weight,
        seed = ctrl$seed
      )
    } else {
      corr_df <- .ssma_build_mnn_correspondences(
        domains = domains,
        mnn_k = ctrl$mnn_k,
        max_pairs_per_pair = ctrl$max_unsup_pairs_per_pair,
        weight = ctrl$corr_edge_weight
      )
    }
  }

  L_joint <- .ssma_build_joint_laplacian(
    within_graphs = within,
    correspondences = corr_df,
    domain_sizes = domain_sizes
  )

  if (verbose) {
    message("ssma_align: assembling reduced generalized eigenproblem")
  }

  if (ctrl$solver == "auto") {
    solver_used <- "reduced"
  } else {
    solver_used <- ctrl$solver
  }

  reduced <- .ssma_build_reduced_system(
    domains = domains,
    within_graphs = within,
    L_joint = L_joint,
    rank_per_domain = ctrl$rank_per_domain,
    svd_tol = ctrl$svd_tol,
    assemble_dense = (solver_used == "reduced")
  )

  solved <- if (solver_used == "reduced") {
    .ssma_solve_reduced_geigen(
      A_r = reduced$A_r,
      B_r = reduced$B_r,
      ncomp = as.integer(ncomp),
      eig_tol = ctrl$eig_tol,
      regularization = ctrl$regularization
    )
  } else if (solver_used == "operator") {
    .ssma_solve_operator_geigen(
      reduced = reduced,
      L_joint = L_joint,
      ncomp = as.integer(ncomp),
      eig_tol = ctrl$eig_tol,
      regularization = ctrl$regularization,
      operator_tol = ctrl$operator_tol,
      operator_maxiter = ctrl$operator_maxiter,
      operator_power_iter = ctrl$operator_power_iter,
      seed = ctrl$seed
    )
  } else {
    stop("Unknown solver: ", solver_used, call. = FALSE)
  }

  rank_offsets <- c(0L, cumsum(reduced$ranks))

  loadings_blocks <- vector("list", n_domains)
  scores_blocks <- vector("list", n_domains)
  for (i in seq_len(n_domains)) {
    ri <- (rank_offsets[i] + 1L):rank_offsets[i + 1L]
    Fi <- reduced$basis[[i]]$P %*% solved$alpha[ri, , drop = FALSE]
    Xi <- as.matrix(domains[[i]]$x)

    loadings_blocks[[i]] <- Fi
    scores_blocks[[i]] <- Xi %*% Fi
  }

  loadings <- do.call(rbind, loadings_blocks)
  scores <- do.call(rbind, scores_blocks)

  feature_blocks <- feature_block_indices(domains)

  serial_diag <- lapply(within, `[[`, "serial")
  names(serial_diag) <- domain_names

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "ssma_align",
    extras = list(
      generalized_eigenvalues = solved$eigenvalues,
      generalized_eigenvalues_all = solved$eigenvalues_all,
      correspondences = corr_df,
      serial = serial_diag,
      reduced_ranks = reduced$ranks,
      domain_names = domain_names,
      solver = solver_used,
      control = ctrl,
      regularization_applied = solved$ridge
    )
  )
}

#' @rdname ssma_align
#' @method ssma_align list
#' @export
ssma_align.list <- function(
  data,
  y = NULL,
  correspondences = NULL,
  serial_index = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = ssma_align_control(),
  ...
) {
  if (!length(data) || !all(vapply(data, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    stop("For list input, `data` must be a non-empty list of matrices/data.frames.", call. = FALSE)
  }

  hd <- as_hyperdesign(data)

  if (missing(y)) {
    return(ssma_align.hyperdesign(
      hd,
      correspondences = correspondences,
      serial_index = serial_index,
      preproc = preproc,
      ncomp = ncomp,
      control = control,
      ...
    ))
  }

  y_quo <- rlang::enquo(y)
  rlang::inject(ssma_align.hyperdesign(
    hd,
    y = !!y_quo,
    correspondences = correspondences,
    serial_index = serial_index,
    preproc = preproc,
    ncomp = ncomp,
    control = control,
    ...
  ))
}

#' @rdname ssma_align
#' @method ssma_align default
#' @export
ssma_align.default <- function(data, ...) {
  stop("No applicable method for ssma_align().", call. = FALSE)
}
