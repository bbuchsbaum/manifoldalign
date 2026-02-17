#' Control settings for `spectral_mnn_align()`
#'
#' @param knn Integer k for within-domain kNN graph construction.
#' @param sigma Numeric bandwidth for heat-kernel affinities. Use `NULL` to
#'   auto-tune per domain via [choose_sigma()].
#' @param use_laplacian Logical; if TRUE (default) use the normalized Laplacian.
#' @param k_embed Optional integer number of eigenpairs used to build heat
#'   kernel signatures. If `NULL`, defaults to `max(3*ncomp, ncomp + 10)` inside
#'   the solver.
#' @param q_descriptors Integer number of heat kernel signature time steps.
#' @param time_mode Either `"auto"` (derive a log-spaced time range from the
#'   retained eigenvalues) or `"fixed"` (use `time_range`).
#' @param time_range Length-2 numeric vector giving min/max times when
#'   `time_mode="fixed"`.
#' @param time_scale Positive scalar multiplier applied to the time grid.
#' @param mnn_k Integer k for mutual nearest-neighbor (MNN) discovery in
#'   descriptor space.
#' @param ratio Ratio test threshold in (0,1]. Used to discard ambiguous
#'   nearest neighbors when `mnn_k >= 2`.
#' @param max_anchors Maximum number of anchor pairs used per domain alignment.
#' @param min_anchors Minimum number of anchors required to estimate a stable
#'   Procrustes rotation. If `NULL`, defaults to `max(10, 3*ncomp)` inside the
#'   solver.
#' @param max_pairs_per_label Maximum correspondences sampled per shared label
#'   when `y` is provided and explicit `correspondences` are absent.
#' @param dist_quantile Quantile of (directed) NN distances kept as candidates
#'   prior to enforcing mutuality and one-to-one pairing. Lower values are more
#'   selective.
#' @param refine_rounds Integer number of ICP-style refinement rounds. Each
#'   round updates anchors in the *aligned* embedding space and recomputes the
#'   Procrustes rotation.
#' @param force_SO Logical; if TRUE (default), enforce det(R)=+1.
#' @param store_basis Logical; if TRUE, store eigenpairs in the returned object.
#' @param store_descriptors Logical; if TRUE, store HKS descriptors in the result.
#' @param seed Integer seed for sampling label-based correspondences.
#' @param verbose Logical toggle for progress messages.
#'
#' @return A list of class `spectral_mnn_align_control`.
#' @examples
#' ctrl <- spectral_mnn_align_control(knn = 10, k_embed = 20)
#' str(ctrl)
#' @export
spectral_mnn_align_control <- function(
  knn = 15L,
  sigma = NULL,
  use_laplacian = TRUE,
  k_embed = NULL,
  q_descriptors = 16L,
  time_mode = c("auto", "fixed"),
  time_range = c(0.1, 50),
  time_scale = 1,
  mnn_k = 10L,
  ratio = 0.9,
  max_anchors = 512L,
  min_anchors = NULL,
  max_pairs_per_label = 100L,
  dist_quantile = 0.3,
  refine_rounds = 1L,
  force_SO = TRUE,
  store_basis = FALSE,
  store_descriptors = FALSE,
  seed = 1L,
  verbose = TRUE
) {
  time_mode <- match.arg(time_mode)

  if (!is.numeric(knn) || length(knn) != 1L || is.na(knn) || knn < 1) {
    stop("`knn` must be a positive scalar.", call. = FALSE)
  }
  if (!is.null(sigma) && (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0)) {
    stop("`sigma` must be NULL or a positive scalar.", call. = FALSE)
  }
  if (!is.logical(use_laplacian) || length(use_laplacian) != 1L || is.na(use_laplacian)) {
    stop("`use_laplacian` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(k_embed) && (!is.numeric(k_embed) || length(k_embed) != 1L || is.na(k_embed) || k_embed < 2)) {
    stop("`k_embed` must be NULL or an integer >= 2.", call. = FALSE)
  }
  if (!is.numeric(q_descriptors) || length(q_descriptors) != 1L || is.na(q_descriptors) || q_descriptors < 1) {
    stop("`q_descriptors` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(time_range) || length(time_range) != 2L || any(!is.finite(time_range)) || any(time_range <= 0) || time_range[2] <= time_range[1]) {
    stop("`time_range` must be a length-2 numeric vector with 0 < min < max.", call. = FALSE)
  }
  if (!is.numeric(time_scale) || length(time_scale) != 1L || is.na(time_scale) || time_scale <= 0) {
    stop("`time_scale` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(mnn_k) || length(mnn_k) != 1L || is.na(mnn_k) || mnn_k < 1) {
    stop("`mnn_k` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(ratio) || length(ratio) != 1L || is.na(ratio) || ratio <= 0 || ratio > 1) {
    stop("`ratio` must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(max_anchors) || length(max_anchors) != 1L || is.na(max_anchors) || max_anchors < 1) {
    stop("`max_anchors` must be >= 1.", call. = FALSE)
  }
  if (!is.null(min_anchors) && (!is.numeric(min_anchors) || length(min_anchors) != 1L || is.na(min_anchors) || min_anchors < 1)) {
    stop("`min_anchors` must be NULL or >= 1.", call. = FALSE)
  }
  if (!is.numeric(max_pairs_per_label) || length(max_pairs_per_label) != 1L || is.na(max_pairs_per_label) || max_pairs_per_label < 1) {
    stop("`max_pairs_per_label` must be >= 1.", call. = FALSE)
  }
  if (!is.numeric(dist_quantile) || length(dist_quantile) != 1L || is.na(dist_quantile) || dist_quantile <= 0 || dist_quantile >= 1) {
    stop("`dist_quantile` must be in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(refine_rounds) || length(refine_rounds) != 1L || is.na(refine_rounds) || refine_rounds < 0) {
    stop("`refine_rounds` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.logical(force_SO) || length(force_SO) != 1L || is.na(force_SO)) {
    stop("`force_SO` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(store_basis) || length(store_basis) != 1L || is.na(store_basis)) {
    stop("`store_basis` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(store_descriptors) || length(store_descriptors) != 1L || is.na(store_descriptors)) {
    stop("`store_descriptors` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      knn = as.integer(knn),
      sigma = if (is.null(sigma)) NULL else as.numeric(sigma),
      use_laplacian = use_laplacian,
      k_embed = if (is.null(k_embed)) NULL else as.integer(k_embed),
      q_descriptors = as.integer(q_descriptors),
      time_mode = time_mode,
      time_range = as.numeric(time_range),
      time_scale = as.numeric(time_scale),
      mnn_k = as.integer(mnn_k),
      ratio = as.numeric(ratio),
      max_anchors = as.integer(max_anchors),
      min_anchors = if (is.null(min_anchors)) NULL else as.integer(min_anchors),
      max_pairs_per_label = as.integer(max_pairs_per_label),
      dist_quantile = as.numeric(dist_quantile),
      refine_rounds = as.integer(refine_rounds),
      force_SO = force_SO,
      store_basis = store_basis,
      store_descriptors = store_descriptors,
      seed = as.integer(seed),
      verbose = verbose
    ),
    class = "spectral_mnn_align_control"
  )
}

resolve_spectral_mnn_align_control <- function(control) {
  defaults <- spectral_mnn_align_control()

  if (missing(control) || is.null(control)) {
    return(defaults)
  }
  if (inherits(control, "spectral_mnn_align_control")) {
    merged <- modifyList(defaults, control)
    return(do.call(spectral_mnn_align_control, merged))
  }
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or spectral_mnn_align_control().", call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  merged <- modifyList(defaults, control)
  do.call(spectral_mnn_align_control, merged)
}

.graph_laplacian_spectral_mnn <- function(A, normalized = TRUE) {
  deg <- Matrix::rowSums(A)
  if (normalized) {
    deg[deg == 0] <- 1
    D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(deg))
    L <- Matrix::Diagonal(nrow(A)) - D_inv_sqrt %*% A %*% D_inv_sqrt
  } else {
    L <- Matrix::Diagonal(x = deg) - A
  }
  methods::as(L, "CsparseMatrix")
}

.rspectra_v0_spectral_mnn <- function(n, seed) {
  # Deterministic (seeded) initial vector for RSpectra/ARPACK. Avoids
  # run-to-run variability from random starting vectors.
  idx <- seq_len(n)
  s <- as.numeric(seed)
  x <- sin((idx + s) * 0.113) + cos((idx + s) * 0.271) + 0.5 * sin((idx + s) * 0.019)
  x <- x - mean(x)
  nx <- sqrt(sum(x * x))
  if (!is.finite(nx) || nx < 1e-12) {
    x <- idx - mean(idx)
    nx <- sqrt(sum(x * x))
  }
  x / nx
}

.with_seed_spectral_mnn <- function(seed, expr) {
  if (is.null(seed) || !is.finite(seed)) {
    return(eval.parent(substitute(expr)))
  }
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    },
    add = TRUE
  )
  set.seed(as.integer(seed))
  eval.parent(substitute(expr))
}

.compute_basis_spectral_mnn <- function(X, k_embed, sigma, knn, use_laplacian, seed = 1L, eig_tol = 1e-12) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 3L) stop("Each domain must have at least 3 rows.", call. = FALSE)

  k_use <- min(as.integer(knn), n - 1L)
  if (k_use < 1L) stop("`knn` too large for domain size.", call. = FALSE)

  sigma_use <- if (is.null(sigma)) {
    # choose_sigma expects data; fail-soft with a conservative fallback
    tryCatch(.with_seed_spectral_mnn(seed, choose_sigma(X)), error = function(e) 0.5)
  } else {
    as.numeric(sigma)
  }

  graph_w <- safe_compute(
    neighborweights::graph_weights(X, k = k_use, weight_mode = "heat",
                                   sigma = sigma_use, neighbor_mode = "knn"),
    "Graph construction failed in spectral_mnn_align"
  )
  A <- neighborweights::adjacency(graph_w)
  A <- (A + Matrix::t(A)) / 2

  L <- .graph_laplacian_spectral_mnn(A, normalized = isTRUE(use_laplacian))

  k_embed_safe <- min(as.integer(k_embed), n - 1L)
  if (k_embed_safe < 2L) stop("`k_embed` too large for domain size.", call. = FALSE)

  # Request a few extra modes to survive multiple near-zero eigenvalues.
  k_request <- min(k_embed_safe + 4L, n - 1L)

  decomp <- if (requireNamespace("RSpectra", quietly = TRUE)) {
    v0 <- tryCatch(.rspectra_v0_spectral_mnn(n, seed = seed), error = function(e) NULL)
    opts <- list(tol = 1e-6, maxitr = 1000)
    if (!is.null(v0)) opts$initvec <- v0
    safe_compute(
      RSpectra::eigs_sym(L, k = k_request, which = "SA",
                         opts = opts),
      "Eigen-decomposition failed in spectral_mnn_align"
    )
  } else if (requireNamespace("PRIMME", quietly = TRUE)) {
    safe_compute(
      PRIMME::eigs_sym(L, NEig = k_request, which = "SA"),
      "Eigen-decomposition failed in spectral_mnn_align"
    )
  } else {
    eig_full <- eigen(as.matrix(L), symmetric = TRUE)
    list(values = eig_full$values[seq_len(k_request)],
         vectors = eig_full$vectors[, seq_len(k_request), drop = FALSE])
  }

  evals <- as.numeric(decomp$values)
  evecs <- as.matrix(decomp$vectors)

  keep <- which(evals > eig_tol & is.finite(evals))
  if (!length(keep)) {
    stop("No non-trivial Laplacian eigenvalues found; try increasing `knn` or checking connectivity.", call. = FALSE)
  }

  keep <- keep[order(evals[keep])]
  keep <- keep[seq_len(min(k_embed_safe, length(keep)))]

  list(
    values = evals[keep],
    vectors = evecs[, keep, drop = FALSE],
    sigma = sigma_use,
    knn = k_use
  )
}

.compute_hks_spectral_mnn <- function(basis, ctrl) {
  compute_hks_descriptors(
    basis = basis,
    q = ctrl$q_descriptors,
    time_mode = ctrl$time_mode,
    time_range = ctrl$time_range,
    spacing = "log",
    time_scale = ctrl$time_scale,
    normalize = "both"
  )
}

.pairs_from_correspondences_spectral_mnn <- function(correspondences, ref_idx, dom_idx, n_ref, n_dom) {
  if (!is.data.frame(correspondences)) {
    stop("`correspondences` must be a data.frame.", call. = FALSE)
  }
  req <- c("domain_i", "index_i", "domain_j", "index_j")
  if (!all(req %in% names(correspondences))) {
    stop("`correspondences` must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  corr <- correspondences
  # Keep only rows involving (ref, dom)
  keep <- (corr$domain_i == ref_idx & corr$domain_j == dom_idx) |
    (corr$domain_i == dom_idx & corr$domain_j == ref_idx)
  corr <- corr[keep, , drop = FALSE]
  if (!nrow(corr)) return(list(i = integer(0), j = integer(0)))

  i <- integer(nrow(corr))
  j <- integer(nrow(corr))
  for (r in seq_len(nrow(corr))) {
    if (corr$domain_i[r] == ref_idx) {
      i[r] <- corr$index_i[r]
      j[r] <- corr$index_j[r]
    } else {
      i[r] <- corr$index_j[r]
      j[r] <- corr$index_i[r]
    }
  }

  ok <- is.finite(i) & is.finite(j) & i >= 1 & i <= n_ref & j >= 1 & j <= n_dom
  list(i = as.integer(i[ok]), j = as.integer(j[ok]))
}

.pairs_from_labels_spectral_mnn <- function(labels_ref, labels_dom, max_pairs_per_label, seed = 1L) {
  set.seed(seed)
  if (is.null(labels_ref) || is.null(labels_dom)) return(list(i = integer(0), j = integer(0)))
  labs1 <- as.character(labels_ref)
  labs2 <- as.character(labels_dom)
  u <- intersect(unique(labs1[!is.na(labs1)]), unique(labs2[!is.na(labs2)]))
  if (!length(u)) return(list(i = integer(0), j = integer(0)))

  i_all <- integer(0)
  j_all <- integer(0)
  for (lab in u) {
    I <- which(labs1 == lab)
    J <- which(labs2 == lab)
    m <- min(length(I), length(J), as.integer(max_pairs_per_label))
    if (m < 1L) next
    I <- sample(I, m)
    J <- sample(J, m)
    i_all <- c(i_all, I)
    j_all <- c(j_all, J)
  }
  list(i = as.integer(i_all), j = as.integer(j_all))
}

.pairs_from_labels_mnn_spectral_mnn <- function(D_ref, D_dom, labels_ref, labels_dom, ctrl) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    return(list(i = integer(0), j = integer(0), d = numeric(0), mode = "labels_mnn_unavailable"))
  }

  D_ref <- scale(as.matrix(D_ref))
  D_dom <- scale(as.matrix(D_dom))
  labs1 <- as.character(labels_ref)
  labs2 <- as.character(labels_dom)

  shared <- intersect(unique(labs1[!is.na(labs1)]), unique(labs2[!is.na(labs2)]))
  if (!length(shared)) {
    return(list(i = integer(0), j = integer(0), d = numeric(0), mode = "labels_mnn_none"))
  }

  k <- as.integer(ctrl$mnn_k)
  ratio <- as.numeric(ctrl$ratio)
  q_keep <- as.numeric(ctrl$dist_quantile)
  max_pairs <- as.integer(ctrl$max_anchors)
  max_per_label <- as.integer(ctrl$max_pairs_per_label)

  # When labels are very coarse, be less selective to avoid anchor collapse
  # into a small subset of "easy" points.
  if (length(shared) <= 4L) {
    q_keep <- max(q_keep, 0.5)
    ratio <- max(ratio, 0.98)
  }

  out_i <- integer(0)
  out_j <- integer(0)
  out_d <- numeric(0)

  for (lab in shared) {
    I <- which(labs1 == lab)
    J <- which(labs2 == lab)
    if (!length(I) || !length(J)) next

    k12 <- min(k, length(I))
    k21 <- min(k, length(J))

    nn12 <- RANN::nn2(D_ref[I, , drop = FALSE], D_dom[J, , drop = FALSE], k = k12)
    nn21 <- RANN::nn2(D_dom[J, , drop = FALSE], D_ref[I, , drop = FALSE], k = k21)

    d1 <- nn12$nn.dists[, 1]
    if (!length(d1) || all(!is.finite(d1))) next

    thresh <- stats::quantile(d1, probs = q_keep, na.rm = TRUE, names = FALSE, type = 7)
    cand <- which(is.finite(d1) & d1 <= thresh)
    if (!length(cand)) next

    if (length(cand) > max_per_label * 20L) {
      cand <- cand[order(d1[cand])[seq_len(max_per_label * 20L)]]
    }

    i_tmp <- integer(0)
    j_tmp <- integer(0)
    d_tmp <- numeric(0)

    for (qq in cand) {
      ii_local <- nn12$nn.idx[qq, 1]
      if (!is.finite(ii_local) || ii_local < 1L) next
      # Mutuality in the label-restricted set
      if (!any(nn21$nn.idx[ii_local, ] == qq)) next

      if (k12 >= 2L) {
        d2 <- nn12$nn.dists[qq, 2]
        if (is.finite(d2) && d2 > 0 && (d1[qq] / d2) > ratio) next
      }

      i_tmp <- c(i_tmp, I[ii_local])
      j_tmp <- c(j_tmp, J[qq])
      d_tmp <- c(d_tmp, d1[qq])
    }

    if (!length(i_tmp)) next
    one <- .greedy_one_to_one_spectral_mnn(i_tmp, j_tmp, d_tmp, max_pairs = min(max_per_label, length(i_tmp)))
    out_i <- c(out_i, one$i)
    out_j <- c(out_j, one$j)
    out_d <- c(out_d, one$d)
  }

  if (!length(out_i)) {
    return(list(i = integer(0), j = integer(0), d = numeric(0), mode = "labels_mnn_empty"))
  }

  one_all <- .greedy_one_to_one_spectral_mnn(out_i, out_j, out_d, max_pairs = max_pairs)
  list(i = one_all$i, j = one_all$j, d = one_all$d, mode = "labels_mnn")
}

.anchor_weights_spectral_mnn <- function(d) {
  if (is.null(d) || !length(d)) return(NULL)
  d <- as.numeric(d)
  ok <- is.finite(d)
  if (!any(ok)) return(NULL)
  d_ok <- d[ok]
  s <- stats::median(d_ok)
  if (!is.finite(s) || s <= 0) {
    s <- mean(d_ok[d_ok > 0])
  }
  if (!is.finite(s) || s <= 0) return(NULL)
  w <- exp(-(d^2) / (2 * s^2))
  w[!ok] <- 0
  if (!is.finite(sum(w)) || sum(w) <= 0) return(NULL)
  w
}

.greedy_one_to_one_spectral_mnn <- function(i, j, d, max_pairs) {
  if (!length(i)) return(list(i = integer(0), j = integer(0), d = numeric(0)))
  ord <- order(d)
  i <- i[ord]; j <- j[ord]; d <- d[ord]
  used_i <- rep(FALSE, max(i))
  used_j <- rep(FALSE, max(j))
  out_i <- integer(0); out_j <- integer(0); out_d <- numeric(0)
  for (idx in seq_along(i)) {
    ii <- i[idx]; jj <- j[idx]
    if (used_i[ii] || used_j[jj]) next
    used_i[ii] <- TRUE; used_j[jj] <- TRUE
    out_i <- c(out_i, ii)
    out_j <- c(out_j, jj)
    out_d <- c(out_d, d[idx])
    if (length(out_i) >= max_pairs) break
  }
  list(i = out_i, j = out_j, d = out_d)
}

.mnn_pairs_spectral_mnn <- function(D_ref, D_dom, ctrl) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN is required for mutual nearest-neighbor discovery.", call. = FALSE)
  }
  D_ref <- scale(as.matrix(D_ref))
  D_dom <- scale(as.matrix(D_dom))
  k <- as.integer(ctrl$mnn_k)
  max_pairs <- as.integer(ctrl$max_anchors)
  q_keep <- as.numeric(ctrl$dist_quantile)
  ratio <- as.numeric(ctrl$ratio)

  nn12 <- RANN::nn2(D_dom, D_ref, k = k)
  nn21 <- RANN::nn2(D_ref, D_dom, k = k)

  j_of_i <- nn12$nn.idx[, 1]
  d1 <- nn12$nn.dists[, 1]

  # Candidate subset to keep loop fast.
  thresh <- stats::quantile(d1, probs = q_keep, na.rm = TRUE, names = FALSE, type = 7)
  cand <- which(is.finite(d1) & d1 <= thresh)
  if (length(cand) > max_pairs * 20L) {
    cand <- cand[order(d1[cand])[seq_len(max_pairs * 20L)]]
  }

  out_i <- integer(0)
  out_j <- integer(0)
  out_d <- numeric(0)

  for (ii in cand) {
    jj <- j_of_i[ii]
    if (jj < 1L) next
    # Mutual test: ii must appear among jj's kNNs back to ref
    if (!any(nn21$nn.idx[jj, ] == ii)) next
    if (k >= 2L) {
      d2 <- nn12$nn.dists[ii, 2]
      if (is.finite(d2) && d2 > 0 && (d1[ii] / d2) > ratio) next
    }
    out_i <- c(out_i, ii)
    out_j <- c(out_j, jj)
    out_d <- c(out_d, d1[ii])
  }

  one2one <- .greedy_one_to_one_spectral_mnn(out_i, out_j, out_d, max_pairs = max_pairs)
  list(i = one2one$i, j = one2one$j, d = one2one$d, mode = "hks_mnn")
}

.direct_nn_pairs_spectral_mnn <- function(D_ref, D_dom, ctrl) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("RANN is required for nearest-neighbor discovery.", call. = FALSE)
  }
  D_ref <- scale(as.matrix(D_ref))
  D_dom <- scale(as.matrix(D_dom))
  k <- as.integer(ctrl$mnn_k)
  max_pairs <- as.integer(ctrl$max_anchors)
  q_keep <- as.numeric(ctrl$dist_quantile)
  ratio <- as.numeric(ctrl$ratio)

  nn12 <- RANN::nn2(D_dom, D_ref, k = k)
  j_of_i <- nn12$nn.idx[, 1]
  d1 <- nn12$nn.dists[, 1]

  thresh <- stats::quantile(d1, probs = q_keep, na.rm = TRUE, names = FALSE, type = 7)
  cand <- which(is.finite(d1) & d1 <= thresh)
  if (!length(cand)) return(list(i = integer(0), j = integer(0), d = numeric(0), mode = "hks_nn"))
  if (length(cand) > max_pairs * 20L) {
    cand <- cand[order(d1[cand])[seq_len(max_pairs * 20L)]]
  }

  i <- cand
  j <- j_of_i[cand]
  d <- d1[cand]

  if (k >= 2L) {
    d2 <- nn12$nn.dists[cand, 2]
    ok <- !(is.finite(d2) & d2 > 0 & (d / d2) > ratio)
    i <- i[ok]; j <- j[ok]; d <- d[ok]
  }

  one2one <- .greedy_one_to_one_spectral_mnn(i, j, d, max_pairs = max_pairs)
  list(i = one2one$i, j = one2one$j, d = one2one$d, mode = "hks_nn")
}

.procrustes_rotation_spectral_mnn <- function(A, B, weights = NULL, force_SO = TRUE) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  if (nrow(A) != nrow(B)) stop("Procrustes: A and B must have same rows.", call. = FALSE)
  if (ncol(A) != ncol(B)) stop("Procrustes: A and B must have same columns.", call. = FALSE)

  Ac <- scale(A, center = TRUE, scale = FALSE)
  Bc <- scale(B, center = TRUE, scale = FALSE)

  if (!is.null(weights)) {
    w <- as.numeric(weights)
    if (length(w) == nrow(Ac) && any(is.finite(w) & w > 0)) {
      w[!is.finite(w) | w < 0] <- 0
      if (sum(w) > 0) {
        ws <- sqrt(w)
        Ac <- Ac * ws
        Bc <- Bc * ws
      }
    }
  }

  M <- crossprod(Bc, Ac)
  sv <- svd(M)
  R <- sv$u %*% t(sv$v)
  if (isTRUE(force_SO) && det(R) < 0) {
    D <- diag(ncol(R))
    D[ncol(R), ncol(R)] <- -1
    R <- sv$u %*% D %*% t(sv$v)
  }
  R
}

.apply_preproc_domains_spectral_mnn <- function(domains, preproc) {
  nb <- length(domains)
  pdata <- domains

  if (is.null(preproc)) {
    for (i in seq_len(nb)) pdata[[i]]$x <- as.matrix(pdata[[i]]$x)
    return(list(pdata = pdata, proc = NULL))
  }

  is_stateful_preproc <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")
  broadcast_preproc <- is.function(preproc) || is_stateful_preproc
  if (broadcast_preproc) {
    if (is_stateful_preproc) {
      preproc <- replicate(
        nb,
        unserialize(serialize(preproc, connection = NULL)),
        simplify = FALSE
      )
    } else {
      preproc <- rep(list(preproc), nb)
    }
  } else if (is.list(preproc) && length(preproc) != nb) {
    stop("Length of `preproc` list (", length(preproc), ") must match number of domains (", nb, ").",
         call. = FALSE)
  }

  proclist <- vector("list", nb)
  for (i in seq_len(nb)) {
    Xi <- as.matrix(domains[[i]]$x)
    pre_i <- preproc[[i]]

    # Functions are applied directly (avoid deprecated multivarious::prep()).
    if (is.function(pre_i)) {
      Xi_tf <- pre_i(Xi)
      pdata[[i]]$x <- Xi_tf
      proclist[[i]] <- pre_i
      next
    }

    if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
      pre_i <- unserialize(serialize(pre_i, connection = NULL))
    }

    if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
      ft <- multivarious::fit_transform(pre_i, Xi)
      pdata[[i]]$x <- ft$transformed
      proclist[[i]] <- ft$preproc
    } else {
      # Back-compat for older multivarious versions
      templ <- multivarious::prep(pre_i)
      Xi_tf <- multivarious::init_transform(templ, Xi)
      pdata[[i]]$x <- Xi_tf
      proc_attr <- attr(Xi_tf, "preproc")
      proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
    }
  }

  feature_sizes <- vapply(pdata, function(x) ncol(x$x), integer(1))
  feature_offsets <- cumsum(c(0L, feature_sizes))
  block_idx_list <- lapply(seq_along(feature_sizes), function(i) {
    if (feature_sizes[i] == 0L) integer(0) else seq.int(feature_offsets[i] + 1L, feature_offsets[i + 1L])
  })
  names(block_idx_list) <- names(pdata)

  proc <- multivarious::concat_pre_processors(proclist, block_idx_list)
  list(pdata = pdata, proc = proc)
}

#' Spectral MNN Alignment
#'
#' A fast, correspondence-optional aligner that:
#' 1) builds a within-domain kNN graph per domain,
#' 2) computes Laplacian eigenmaps,
#' 3) constructs intrinsic Heat Kernel Signatures (HKS),
#' 4) discovers anchors via mutual nearest neighbors (or optional supervision),
#' 5) aligns domains to a reference via Procrustes.
#'
#' @param data A `hyperdesign` object or list of matrices.
#' @param ... Additional arguments passed to methods.
#'
#' @return A `multiblock_biprojector` with subclass `spectral_mnn_align`.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 10, 5)
#' X2 <- matrix(rnorm(50), 10, 5)
#' data_list <- list(domain1 = X1, domain2 = X2)
#' result <- spectral_mnn_align(data_list, ncomp = 3)
#' dim(result$scores)
#' }
#' @export
spectral_mnn_align <- function(data, ...) {
  UseMethod("spectral_mnn_align")
}

#' @rdname spectral_mnn_align
#' @method spectral_mnn_align hyperdesign
#' @param y Optional unquoted column name in design tables for label-based anchors
#' @param correspondences Optional data.frame with explicit anchor pairs
#' @param preproc Preprocessing function or list. Default: multivarious::center()
#' @param ncomp Number of components to extract. Default: 20
#' @param ref_idx Reference domain index. Default: 1
#' @param control Control settings from spectral_mnn_align_control()
#' @export
spectral_mnn_align.hyperdesign <- function(
  data,
  y = NULL,
  correspondences = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  ref_idx = 1L,
  control = spectral_mnn_align_control(),
  ...
) {
  ctrl <- resolve_spectral_mnn_align_control(control)
  verbose <- isTRUE(ctrl$verbose)

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }

  resolved <- resolve_hyperdesign(data)
  domains <- resolved$domains
  domain_names <- resolved$domain_names
  n_domains <- resolved$n_domains

  if (n_domains < 2L) stop("spectral_mnn_align requires at least 2 domains.", call. = FALSE)
  if (!is.numeric(ref_idx) || length(ref_idx) != 1L || is.na(ref_idx) || ref_idx < 1L || ref_idx > n_domains) {
    stop("`ref_idx` must be an integer in [1, n_domains].", call. = FALSE)
  }
  ref_idx <- as.integer(ref_idx)

  if (verbose) message("spectral_mnn_align: preprocessing domains")
  pre <- .apply_preproc_domains_spectral_mnn(domains, preproc = preproc)
  pdata <- pre$pdata
  proc <- pre$proc
  names(pdata) <- domain_names

  # Extract labels if requested (avoid evaluating `y`; support unquoted column names)
  y_name <- NULL
  labels <- NULL
  if (!missing(y)) {
    y_quo <- rlang::enquo(y)
    if (!rlang::quo_is_null(y_quo)) {
      y_name <- rlang::as_name(rlang::ensym(y))
      labels <- vector("list", n_domains)
      for (i in seq_len(n_domains)) {
        di <- pdata[[i]]$design
        if (is.null(di) || !is.data.frame(di) || !(y_name %in% names(di))) {
          stop("Label column `", y_name, "` not found in domain ", i, ".", call. = FALSE)
        }
        labels[[i]] <- di[[y_name]]
      }
    }
  }

  k_embed <- ctrl$k_embed
  if (is.null(k_embed)) {
    k_embed <- max(3L * as.integer(ncomp), as.integer(ncomp) + 10L)
  }
  k_embed <- as.integer(k_embed)

  if (verbose) message("spectral_mnn_align: computing spectral bases (k_embed=", k_embed, ")")
  bases <- vector("list", n_domains)
  for (i in seq_len(n_domains)) {
    bases[[i]] <- .compute_basis_spectral_mnn(
      pdata[[i]]$x,
      k_embed = k_embed,
      sigma = ctrl$sigma,
      knn = ctrl$knn,
      use_laplacian = ctrl$use_laplacian,
      seed = as.integer(ctrl$seed) + as.integer(i) * 1009L
    )
  }

  if (verbose) message("spectral_mnn_align: computing HKS descriptors (q=", ctrl$q_descriptors, ")")
  desc <- vector("list", n_domains)
  for (i in seq_len(n_domains)) desc[[i]] <- .compute_hks_spectral_mnn(bases[[i]], ctrl)

  # Embeddings used for output/alignment
  k_avail <- vapply(bases, function(b) ncol(b$vectors), integer(1))
  k_use <- min(as.integer(ncomp), min(k_avail))
  if (k_use < as.integer(ncomp) && verbose) {
    message("spectral_mnn_align: reducing ncomp to ", k_use, " due to small graphs")
  }

  E <- lapply(bases, function(b) b$vectors[, seq_len(k_use), drop = FALSE])
  # Center embeddings per domain for stable Procrustes/loadings
  E <- lapply(E, function(M) scale(as.matrix(M), center = TRUE, scale = FALSE))

  # Estimate per-domain rotations into ref space
  rotations <- vector("list", n_domains)
  rotations[[ref_idx]] <- diag(k_use)
  anchors <- vector("list", n_domains)
  anchors[[ref_idx]] <- NULL

  min_anchors <- ctrl$min_anchors
  if (is.null(min_anchors)) min_anchors <- max(10L, 3L * k_use)
  min_anchors <- as.integer(min_anchors)

  for (g in seq_len(n_domains)) {
    if (g == ref_idx) next
    n_ref <- nrow(E[[ref_idx]])
    n_dom <- nrow(E[[g]])

    pair_info <- NULL
    mode <- NULL

    if (!is.null(correspondences)) {
      pair_info <- .pairs_from_correspondences_spectral_mnn(correspondences, ref_idx, g, n_ref, n_dom)
      mode <- "correspondences"
    } else if (!is.null(labels)) {
      pair_info <- .pairs_from_labels_mnn_spectral_mnn(desc[[ref_idx]], desc[[g]],
                                                       labels[[ref_idx]], labels[[g]], ctrl)
      mode <- pair_info$mode
      if (length(pair_info$i) < min_anchors) {
        # Fall back to purely unsupervised anchors if labels are too coarse.
        pair_info <- .mnn_pairs_spectral_mnn(desc[[ref_idx]], desc[[g]], ctrl)
        mode <- paste0("labels_fallback+", pair_info$mode)
        if (length(pair_info$i) < min_anchors) {
          pair_info <- .direct_nn_pairs_spectral_mnn(desc[[ref_idx]], desc[[g]], ctrl)
          mode <- paste0("labels_fallback+", pair_info$mode)
        }
      }
    } else {
      pair_info <- .mnn_pairs_spectral_mnn(desc[[ref_idx]], desc[[g]], ctrl)
      mode <- pair_info$mode
      if (length(pair_info$i) < min_anchors) {
        pair_info <- .direct_nn_pairs_spectral_mnn(desc[[ref_idx]], desc[[g]], ctrl)
        mode <- pair_info$mode
      }
    }

    if (length(pair_info$i) < max(2L, k_use)) {
      if (verbose) message("spectral_mnn_align: domain ", g, " has insufficient anchors (", length(pair_info$i), "); using identity rotation")
      rotations[[g]] <- diag(k_use)
      anchors[[g]] <- data.frame(i_ref = integer(0), i_dom = integer(0), dist = numeric(0), mode = rep(mode, 0L),
                                 stringsAsFactors = FALSE)
      next
    }

    # Fit initial rotation from anchors in embedding space
    A <- E[[ref_idx]][pair_info$i, , drop = FALSE]
    B <- E[[g]][pair_info$j, , drop = FALSE]
    w1 <- .anchor_weights_spectral_mnn(pair_info$d)
    Rg <- .procrustes_rotation_spectral_mnn(A, B, weights = w1, force_SO = ctrl$force_SO)

    # Optional refinement in aligned space (cheap ICP-style)
    rounds <- as.integer(ctrl$refine_rounds)
    if (rounds > 0L) {
      for (rr in seq_len(rounds)) {
        B_aligned <- E[[g]] %*% Rg
        # Update anchors directly in the aligned embedding space.
        if (is.null(correspondences) && !is.null(labels)) {
          pair2 <- .pairs_from_labels_mnn_spectral_mnn(
            E[[ref_idx]], B_aligned,
            labels[[ref_idx]], labels[[g]],
            ctrl
          )
          if (length(pair2$i) < min_anchors) {
            pair2 <- .mnn_pairs_spectral_mnn(E[[ref_idx]], B_aligned, ctrl)
          }
        } else {
          pair2 <- .mnn_pairs_spectral_mnn(E[[ref_idx]], B_aligned, ctrl)
        }
        if (length(pair2$i) < min_anchors) pair2 <- .direct_nn_pairs_spectral_mnn(E[[ref_idx]], B_aligned, ctrl)
        if (length(pair2$i) < max(2L, k_use)) break
        A2 <- E[[ref_idx]][pair2$i, , drop = FALSE]
        B2 <- E[[g]][pair2$j, , drop = FALSE]
        w2 <- .anchor_weights_spectral_mnn(pair2$d)
        Rg <- .procrustes_rotation_spectral_mnn(A2, B2, weights = w2, force_SO = ctrl$force_SO)
        pair_info <- pair2
        mode <- paste0(mode, "+refine")
      }
    }

    rotations[[g]] <- Rg
    dvec <- pair_info$d
    if (is.null(dvec) || length(dvec) != length(pair_info$i)) {
      dvec <- rep(NA_real_, length(pair_info$i))
    }
    anchors[[g]] <- data.frame(
      i_ref = as.integer(pair_info$i),
      i_dom = as.integer(pair_info$j),
      dist = as.numeric(dvec),
      mode = mode,
      stringsAsFactors = FALSE
    )
  }

  # Apply rotations
  E_aligned <- lapply(seq_len(n_domains), function(i) {
    Ei <- E[[i]]
    Ri <- rotations[[i]]
    if (is.null(Ri)) Ri <- diag(k_use)
    Ei %*% Ri
  })

  scores <- do.call(rbind, E_aligned)

  loadings <- do.call(rbind, lapply(seq_len(n_domains), function(i) {
    Xi <- as.matrix(pdata[[i]]$x)
    Matrix::crossprod(Xi, E_aligned[[i]])
  }))

  feature_blocks <- feature_block_indices(pdata)

  extras <- list(
    ref_idx = ref_idx,
    rotations = rotations,
    anchors = anchors,
    control = ctrl
  )
  if (isTRUE(ctrl$store_basis)) extras$basis <- bases
  if (isTRUE(ctrl$store_descriptors)) extras$descriptors <- desc

  new_alignment_result(
    scores = scores,
    loadings = as.matrix(loadings),
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "spectral_mnn_align",
    extras = extras
  )
}

#' @rdname spectral_mnn_align
#' @method spectral_mnn_align list
#' @export
spectral_mnn_align.list <- function(
  data,
  ...
) {
  if (!is.list(data) || length(data) < 2) {
    stop("Provide a list with at least two matrices.", call. = FALSE)
  }
  if (!all(vapply(data, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    stop("Every element in the list must be a numeric matrix.", call. = FALSE)
  }
  strata <- lapply(seq_along(data), function(i) {
    list(x = as.matrix(data[[i]]), design = data.frame(row_id = seq_len(nrow(data[[i]]))))
  })
  nm <- names(data)
  if (is.null(nm) || any(nm == "")) nm <- paste0("domain_", seq_along(strata))
  names(strata) <- nm
  class(strata) <- "hyperdesign"
  spectral_mnn_align.hyperdesign(strata, ...)
}

#' @rdname spectral_mnn_align
#' @method spectral_mnn_align default
#' @export
spectral_mnn_align.default <- function(data, ...) {
  stop("No applicable method: data must be a 'hyperdesign' or list of matrices.", call. = FALSE)
}
