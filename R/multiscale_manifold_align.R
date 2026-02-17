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
  if (is.null(preproc)) return(domains)

  n_domains <- length(domains)
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
    preproc_list <- preproc
  } else {
    stop("`preproc` must be NULL, a function/preprocessor, or a list matching domain count.", call. = FALSE)
  }

  out <- domains
  for (i in seq_along(domains)) {
    Xi <- as.matrix(domains[[i]]$x)
    pre_i <- preproc_list[[i]]
    if (is.null(pre_i)) {
      out[[i]]$x <- Xi
      next
    }
    if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
      pre_i <- unserialize(serialize(pre_i, connection = NULL))
    }
    templ <- multivarious::prep(pre_i)
    out[[i]]$x <- multivarious::init_transform(templ, Xi)
  }
  out
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
  domains <- .msa_preprocess_domains(resolved$domains, preproc = preproc)
  n_domains <- length(domains)
  domain_sizes <- vapply(domains, function(d) nrow(d$x), integer(1))
  n_total <- sum(domain_sizes)

  if (verbose || isTRUE(ctrl$verbose)) {
    message("multiscale_manifold_align: domains=", n_domains, ", nodes=", n_total)
  }

  within <- lapply(domains, function(dom) .msa_build_within_graph(dom$x, knn = ctrl$knn, sigma = ctrl$sigma))
  d_within <- unlist(lapply(within, `[[`, "degree"), use.names = FALSE)
  d_within[d_within <= 0] <- 1

  y_quo <- if (missing(y)) rlang::quo(NULL) else rlang::enquo(y)
  label_list <- .msa_extract_labels(domains, y_quo)

  corr_df <- .msa_validate_correspondences(correspondences, domain_sizes = domain_sizes, default_weight = ctrl$cross_edge_weight)
  if (!nrow(corr_df)) {
    if (!is.null(label_list)) {
      corr_df <- .msa_build_label_correspondences(
        labels_list = label_list,
        max_pairs_per_label = ctrl$max_pairs_per_label,
        weight = ctrl$cross_edge_weight,
        seed = ctrl$seed
      )
    } else {
      corr_df <- .msa_build_mnn_correspondences(
        domains = domains,
        mnn_k = ctrl$mnn_k,
        max_pairs_per_pair = ctrl$max_unsup_pairs_per_pair,
        weight = ctrl$cross_edge_weight
      )
    }
  }

  L_joint <- .msa_build_joint_laplacian(within_graphs = within, correspondences = corr_df, domain_sizes = domain_sizes)

  svd_list <- lapply(seq_along(domains), function(i) {
    .msa_domain_svd(
      X = domains[[i]]$x,
      degree = within[[i]]$degree,
      rank_per_domain = ctrl$rank_per_domain,
      svd_tol = ctrl$svd_tol
    )
  })

  V_big <- Matrix::bdiag(lapply(svd_list, function(x) Matrix::Matrix(x$V, sparse = TRUE)))
  d_inv_sqrt <- 1 / sqrt(pmax(d_within, 1e-12))
  D_inv <- Matrix::Diagonal(x = d_inv_sqrt)
  DV <- D_inv %*% V_big
  T_reduced <- Matrix::crossprod(V_big, D_inv %*% (L_joint %*% DV))
  T_reduced <- Matrix::Matrix((T_reduced + Matrix::t(T_reduced)) / 2, sparse = FALSE)

  backend_used <- ctrl$backend
  if (ctrl$backend == "randomized_dwt") {
    warning("randomized_dwt backend not yet implemented; using spectral backend.", call. = FALSE)
    backend_used <- "spectral"
  }

  eigvals <- NULL
  eigvecs <- NULL
  r <- nrow(T_reduced)
  k_req <- min(max(as.integer(ncomp), 2L), as.integer(ctrl$eigen_k_max), max(1L, r - 1L))
  if (k_req < 1L) stop("Reduced operator rank is too small for eigendecomposition.", call. = FALSE)

  if (requireNamespace("RSpectra", quietly = TRUE) && k_req < (r - 1L)) {
    eig <- safe_compute(
      RSpectra::eigs_sym(T_reduced, k = k_req, which = "SA"),
      "Spectral decomposition failed in multiscale_manifold_align()."
    )
    eigvals <- as.numeric(eig$values)
    eigvecs <- as.matrix(eig$vectors)
  } else {
    eig <- eigen(as.matrix(T_reduced), symmetric = TRUE)
    ord <- order(eig$values, decreasing = FALSE)
    eigvals <- as.numeric(eig$values[ord][seq_len(k_req)])
    eigvecs <- as.matrix(eig$vectors[, ord, drop = FALSE][, seq_len(k_req), drop = FALSE])
  }

  ord <- order(eigvals, decreasing = FALSE)
  eigvals <- eigvals[ord]
  eigvecs <- eigvecs[, ord, drop = FALSE]

  keep <- which(is.finite(eigvals) & eigvals > ctrl$eigen_tol)
  if (!length(keep)) {
    stop("No non-trivial reduced eigenvalues were found above `eigen_tol`.", call. = FALSE)
  }

  eigvals_pos <- eigvals[keep]
  eigvecs_pos <- eigvecs[, keep, drop = FALSE]
  k_out <- min(as.integer(ncomp), ncol(eigvecs_pos))
  E_out <- eigvecs_pos[, seq_len(k_out), drop = FALSE]

  scores <- as.matrix(DV %*% E_out)

  level_dims <- integer(0)
  for (j in seq_len(ctrl$max_levels)) {
    t_j <- 2^(j - 1L)
    thr <- log(1 / ctrl$eps_rank) / t_j
    p_j <- sum(eigvals_pos <= thr)
    if (p_j < 1L) break
    level_dims <- c(level_dims, p_j)
  }
  if (!length(level_dims)) level_dims <- k_out
  level_dims <- as.integer(level_dims)
  level_dims <- pmin(level_dims, ncol(eigvecs_pos))
  levels <- lapply(level_dims, function(pj) eigvecs_pos[, seq_len(pj), drop = FALSE])

  feature_blocks <- feature_block_indices(domains)
  total_features <- sum(vapply(domains, function(d) ncol(d$x), integer(1)))
  domain_offsets <- cumsum(domain_sizes)
  domain_starts <- c(1L, head(domain_offsets, -1L) + 1L)
  loading_blocks <- lapply(seq_along(domains), function(i) {
    Xi <- as.matrix(domains[[i]]$x)
    Si <- scores[domain_starts[i]:domain_offsets[i], , drop = FALSE]
    XtX <- crossprod(Xi)
    ridge <- 1e-6 * mean(diag(XtX) + 1e-12)
    solve(XtX + diag(ridge, nrow = ncol(Xi)), crossprod(Xi, Si))
  })
  loadings <- do.call(rbind, loading_blocks)

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = NULL,
    feature_blocks = feature_blocks,
    subclass = "multiscale_manifold_align",
    extras = list(
      backend = backend_used,
      levels = levels,
      level_dims = level_dims,
      reduced_eigenvalues = eigvals_pos,
      correspondences = corr_df,
      domain_names = resolved$domain_names,
      control = ctrl
    )
  )
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
