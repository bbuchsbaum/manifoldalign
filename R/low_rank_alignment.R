#' Create a similarity function from a label similarity matrix
#'
#' Builds a function that maps a vector of labels to a sparse similarity matrix.
#' Unknown/NA labels receive all-zero rows/columns, and the result is guaranteed
#' symmetric with zero diagonal.
#'
#' This is useful as `simfun` for [`lowrank_align()`], where labels may include
#' missing values in semi-supervised settings.
#'
#' @param S A square similarity matrix whose row/column names define the label
#'   vocabulary. Only labels present in `rownames(S)` are assigned similarities.
#' @param na_value Deprecated; ignored (kept for backward compatibility). NA and
#'   unknown labels always map to zero similarity to preserve sparsity.
#' @return A function `f(labels)` that returns an `n x n` sparse similarity matrix.
#' @examples
#' S <- diag(3)
#' dimnames(S) <- list(letters[1:3], letters[1:3])
#' simfun <- createSimFun(S)
#' M <- simfun(c("a", NA, "c"))
#' Matrix::nnzero(M)
#' @export
createSimFun <- function(S, na_value = 0) {
  # Extract the labels from the similarity matrix
  label_names <- rownames(S)

  if (!isTRUE(all.equal(na_value, 0))) {
    warning("createSimFun ignores non-zero na_value to preserve sparsity; using 0")
  }

  # The function to return
  function(labels) {
    # Find indices of input labels in the label_names
    indices <- match(labels, label_names)
    n <- length(labels)

    # Sparse initialise; fill recognised labels only
    M <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(n, n))

    ok <- !is.na(indices)
    if (any(ok)) {
      M[ok, ok] <- S[indices[ok], indices[ok], drop = FALSE]
    }

    # Guarantee symmetry & zero diagonal
    M <- (M + Matrix::t(M)) / 2
    Matrix::diag(M) <- 0

    M
  }
}

#' @noRd
.compute_R_factors <- function(X, sv_thresh = 1) {
  svd_X <- svd(t(X))
  a_raw <- 1 - 1 / (svd_X$d^2)
  keep <- (svd_X$d > sv_thresh) & (a_raw > 0)

  if (!any(keep)) {
    return(list(U = matrix(0, nrow(X), 0),
                a = numeric(0)))
  }

  list(U = svd_X$v[, keep, drop = FALSE],
       a = a_raw[keep])
}

#' @noRd
compute_R <- function(X, sv_thresh = 1) {
  fac <- .compute_R_factors(X, sv_thresh = sv_thresh)

  if (length(fac$a) == 0) {
    return(Matrix::Matrix(0, nrow(X), nrow(X), sparse = TRUE))
  }

  R <- fac$U %*% (Matrix::Diagonal(x = fac$a)) %*% t(fac$U)
  Matrix::forceSymmetric(Matrix(R, sparse = TRUE))
}

#' @noRd
compute_rank_matrices <- function(strata, sv_thresh = 1) {
  purrr::map(strata, function(x) compute_R(x$x, sv_thresh = sv_thresh))
}

#' @noRd
.factors_to_R <- function(fac, n) {
  if (!length(fac$a)) {
    return(Matrix::Matrix(0, n, n, sparse = TRUE))
  }

  R <- fac$U %*% (Matrix::Diagonal(x = fac$a)) %*% t(fac$U)
  Matrix::forceSymmetric(Matrix(R, sparse = TRUE))
}

#' @noRd
.apply_IminusR <- function(xk, Uk, ak) {
  if (!length(ak)) {
    return(xk)
  }

  tUx <- as.numeric(crossprod(Uk, xk))
  xk - Uk %*% (ak * tUx)
}

#' @noRd
.make_M_op <- function(Rfac_list, sizes) {
  function(x, ...) {
    y <- x
    offset <- 0L

    for (i in seq_along(Rfac_list)) {
      n_i <- sizes[i]
      idx <- (offset + 1L):(offset + n_i)
      fac <- Rfac_list[[i]]

      tmp <- .apply_IminusR(x[idx], fac$U, fac$a)
      y[idx] <- .apply_IminusR(tmp, fac$U, fac$a)

      offset <- offset + n_i
    }

    y
  }
}

#' @noRd
.lambda_max_power <- function(op, n, n_iter = 30L, tol = 1e-6, seed = 1L) {
  v <- .with_local_seed(seed, stats::rnorm(n))
  v <- v / sqrt(sum(v * v))
  lambda_prev <- 0

  for (it in seq_len(n_iter)) {
    w <- op(v)
    lambda <- sum(v * w)
    norm_w <- sqrt(sum(w * w))

    if (norm_w == 0) {
      break
    }

    v <- w / norm_w

    if (abs(lambda - lambda_prev) <= tol * (abs(lambda_prev) + 1e-12)) {
      break
    }

    lambda_prev <- lambda
  }

  lambda_prev
}

# Keep only the entries of a symmetric similarity matrix that join samples
# from different domains (Eq. 12 of the LRA paper). Works on the stored
# entries, so the cost is O(nnz) rather than O(n^2).
#' @noRd
.lowrank_cross_domain_edges <- function(C, block_id) {
  n <- nrow(C)
  C <- methods::as(methods::as(methods::as(C, "dMatrix"), "generalMatrix"), "TsparseMatrix")
  keep <- block_id[C@i + 1L] != block_id[C@j + 1L]
  Matrix::sparseMatrix(
    i = C@i[keep] + 1L,
    j = C@j[keep] + 1L,
    x = C@x[keep],
    dims = c(n, n)
  )
}

# Extreme eigenpairs of a symmetric matrix: PRIMME when installed, otherwise a
# dense eigendecomposition for small problems or RSpectra. Returns the k
# smallest ("SA") or largest ("LA") eigenpairs.
#' @noRd
.lowrank_eigs_sym <- function(A, k, which = c("SA", "LA")) {
  which <- match.arg(which)
  n <- nrow(A)
  k <- as.integer(min(k, n))
  if (.primme_available()) {
    eig <- PRIMME::eigs_sym(A, NEig = k, which = which,
                            method = "PRIMME_DEFAULT_MIN_MATVECS")
    return(list(values = eig$values, vectors = eig$vectors))
  }
  if (n <= 2000L || k >= n - 1L) {
    e <- eigen(as.matrix(A), symmetric = TRUE)
    # eigen() returns eigenvalues in decreasing order.
    idx <- if (which == "SA") rev(seq_len(n))[seq_len(k)] else seq_len(k)
    return(list(values = e$values[idx], vectors = e$vectors[, idx, drop = FALSE]))
  }
  eig <- RSpectra::eigs_sym(A, k = k, which = which)
  list(values = eig$values, vectors = eig$vectors)
}

#' @noRd
.make_Z_op <- function(M_op, Ls, mu) {
  function(x, ...) {
    (1 - mu) * M_op(x) + 2 * mu * as.numeric(Ls %*% x)
  }
}

# Orthonormal basis Q (n x K, sparse) of the per-domain constant vectors for
# the domain index `block_id` (one entry per sample): column k is the
# indicator of domain k divided by sqrt(n_k). The embedding is solved on the
# orthogonal complement of span(Q).
#' @noRd
.lowrank_constant_basis <- function(block_id) {
  g <- match(block_id, unique(block_id))
  sizes <- tabulate(g)
  Matrix::sparseMatrix(
    i = seq_along(g),
    j = g,
    x = 1 / sqrt(sizes[g]),
    dims = c(length(g), length(sizes))
  )
}

# Projector P = I - Q Q' onto the complement of the per-domain constant
# vectors, as a function: it removes the per-domain mean of each column of x
# in O(n) per column, without forming P or Q.
#' @noRd
.lowrank_projector <- function(block_id) {
  g <- match(block_id, unique(block_id))
  sizes <- tabulate(g)
  function(x) {
    if (is.null(dim(x))) {
      return(x - (as.numeric(rowsum(x, g, reorder = FALSE)) / sizes)[g])
    }
    x <- as.matrix(x)
    means <- rowsum(x, g, reorder = FALSE) / sizes
    dimnames(means) <- NULL
    x - means[g, , drop = FALSE]
  }
}

#' @noRd
.lowrank_project_complement <- function(x, block_id) {
  .lowrank_projector(block_id)(x)
}

# Smallest eigenpairs of the symmetric matrix A restricted to the orthogonal
# complement of the per-domain constant vectors, i.e. of P A P on range(P).
# `sigma` must exceed the largest eigenvalue of A; the constraint directions
# are placed at `sigma` so that they are never among the smallest eigenpairs.
# PRIMME's own `ortho` argument is not used because it does not converge when
# the constrained span is not invariant under A.
#' @noRd
.lowrank_eigs_restricted <- function(A, k, block_id, sigma, dense_max = 2000L) {
  n <- nrow(A)
  n_free <- n - length(unique(block_id))
  k <- as.integer(min(k, n_free))
  P <- .lowrank_projector(block_id)
  deflated <- function(x, ...) {
    px <- P(x)
    apx <- A %*% px
    apx <- if (is.null(dim(x))) as.numeric(apx) else as.matrix(apx)
    P(apx) + sigma * (x - px)
  }

  if (.primme_available()) {
    eig <- PRIMME::eigs_sym(deflated, NEig = k, which = "SA", n = n, isreal = TRUE,
                            method = "PRIMME_DEFAULT_MIN_MATVECS")
    vals <- eig$values
    vecs <- eig$vectors
  } else if (n <= dense_max) {
    Qd <- as.matrix(.lowrank_constant_basis(block_id))
    N <- qr.Q(qr(Qd), complete = TRUE)[, -seq_len(ncol(Qd)), drop = FALSE]
    e <- eigen(crossprod(N, as.matrix(A %*% N)), symmetric = TRUE)
    # eigen() returns eigenvalues in decreasing order.
    idx <- rev(seq_len(n_free))[seq_len(k)]
    vals <- e$values[idx]
    vecs <- N %*% e$vectors[, idx, drop = FALSE]
  } else {
    eig <- RSpectra::eigs_sym(deflated, k = k, n = n, which = "SA")
    vals <- eig$values
    vecs <- eig$vectors
  }

  ord <- order(vals)
  list(values = vals[ord],
       vectors = P(vecs[, ord, drop = FALSE]))
}

#' @noRd
.make_shifted_op <- function(Z_op, n, c_bound) {
  function(x, ...) {
    x - (1 / c_bound) * Z_op(x)
  }
}

#' Low-rank Alignment
#'
#' Performs low-rank alignment using eigenvalue decomposition. Balances 
#' low-rank structure with similarity-based constraints. Supports 
#' semi-supervised learning with missing labels.
#'
#' @param data Input data object
#' @param y Variable name for labels (unquoted). Can contain NA values for 
#'   unlabeled samples in semi-supervised learning scenarios.
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Low-rank alignment (Boucher, Carey, Mahadevan & Dyar, 2015) builds the
#' matrix Z = (1-μ) * M + 2μ * L, where M = (I-R)ᵀ(I-R) is the block-diagonal
#' low-rank reconstruction term (R_k = U_k diag(1 - 1/s²) U_kᵀ per domain) and
#' L is the graph Laplacian of the cross-domain similarity matrix. The method
#' balances preserving low-rank structure (μ=0) with enforcing similarity
#' constraints (μ=1).
#'
#' **Estimand.** The embedding F (n x ncomp, rows are samples of all domains
#' stacked in order) minimises trace(Fᵀ Z F) subject to Fᵀ F = I and
#' Fᵀ e_k = 0 for every domain indicator vector e_k, i.e. the `ncomp`
#' smallest eigenvectors of Z restricted to the orthogonal complement of the
#' per-domain constant vectors. Each domain's coordinates therefore have zero
#' mean. As in Laplacian eigenmaps, this removes trivial solutions: with
#' centred features every R_k annihilates the domain constant, so the global
#' constant vector is an exact eigenvector of Z (eigenvalue 1-μ) and, without
#' the constraint, it is returned as a zero-variance component whenever the
#' label term dominates. Per-domain offsets carry no within-domain structure
#' and cannot be reproduced by the no-intercept linear projection fitted to
#' centred data. Both solvers solve this same restricted problem.
#'
#' **Semi-supervised Learning Support:**
#' The algorithm handles NA labels gracefully. Unlabeled samples:
#' - Still contribute to the low-rank structure term M through their data
#' - Do not participate in the similarity constraints (L term)
#' - Receive coordinates in the joint embedding space
#' - Are isolated nodes of the similarity graph. For μ < 1, M is positive
#'   definite, so they do not create zero eigenvalues; for μ = 1 the zero
#'   eigenvalues of Z (null modes of L) are skipped.
#'
#' @return The return value depends on the specific method. For hyperdesign 
#'   objects, returns a multiblock_biprojector object containing alignment 
#'   results, eigenvectors, preprocessing information, and metadata.
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' library(multidesign)
#' 
#' # Create synthetic data
#' set.seed(123)
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#' 
#' # Create similarity function (NA-tolerant)
#' S <- matrix(runif(10*10), 10, 10)
#' S <- abs(cor(S))
#' row.names(S) <- colnames(S) <- 1:10
#' simfun <- createSimFun(S)  # Handles NA labels automatically
#' 
#' # Create hyperdesign and run alignment
#' hd <- hyperdesign(list(d1, d2, d3))
#' result <- lowrank_align(hd, y, simfun=simfun)
#' 
#' # Semi-supervised learning with missing labels
#' d1_semi <- d1
#' d1_semi$design$y[1:3] <- NA  # Mark some samples as unlabeled
#' d2_semi <- d2
#' d2_semi$design$y[1:2] <- NA
#' hd_semi <- hyperdesign(list(d1_semi, d2_semi, d3))
#' result_semi <- lowrank_align(hd_semi, y, simfun=simfun)
#' }
#'
#' @export
lowrank_align <- function(data, y, ...) {
  UseMethod("lowrank_align")
}

#' @rdname lowrank_align
#' @method lowrank_align hyperdesign
#'
#' @param preproc Preprocessing (default: center()): a single preprocessor
#'   broadcast to all domains, a list with one entry per domain, or `NULL` for
#'   none. `NULL` entries are recorded as an identity transform.
#' @param ncomp Number of components to extract (default: 2)
#' @param simfun Function to compute similarity matrix from labels. Should 
#'   handle NA labels gracefully (e.g., created with createSimFun())
#' @param mu Balance parameter between low-rank (μ=0) and similarity (μ=1)
#'   terms (default: 0.5). Must be a single number in `[0, 1]`.
#' @param lambda Regularization parameter for glmnet. If NULL (default), uses 
#'   cross-validation to select optimal lambda via cv.glmnet. If specified, 
#'   uses the provided value directly.
#' @param scale_M Logical. If TRUE, scales M matrix to have similar eigenvalue 
#'   magnitude as L. This can improve numerical conditioning but changes the 
#'   mathematical objective. When enabled, consider adjusting the mu parameter 
#'   accordingly (default: FALSE)
#' @param n_cores Number of threads for PRIMME eigenvalue computations. If 
#'   NULL (default), uses system default. Set to 1 for reproducible results 
#'   across systems.
#' @param sv_thresh Singular value threshold used when forming R (default: 1).
#'   Only singular values strictly greater than `max(sv_thresh, 1)` are kept:
#'   the closed-form factor `1 - 1/s^2` of Eq. 11 of the original Low-Rank
#'   Alignment paper is non-positive for `s <= 1`, so values of `sv_thresh`
#'   below 1 behave exactly like `sv_thresh = 1`. Larger values discard more
#'   directions. Must be a finite non-negative number.
#' @param solver Eigen solver backend. `"explicit"` (default) forms the
#'   matrices and uses PRIMME when it is installed (otherwise a dense
#'   eigendecomposition for up to 2000 samples, or RSpectra), matching the
#'   original implementation. `"operator"` keeps low-rank factors and uses
#'   RSpectra with a matrix-vector operator to reduce memory and runtime.
#'   Both restrict the eigenproblem to the complement of the per-domain
#'   constant vectors by projecting inside the matrix-vector product (the
#'   projector is never formed), so they target the same embedding.
#'
#' @details
#' The scale_M parameter controls whether to apply eigenvalue-based scaling:
#' - scale_M = FALSE (default): Uses original formulation 
#'   Z = (1-μ) * M + 2μ * L
#' - scale_M = TRUE: Applies scaling M := M * (λ₁(L)/λ₁(M)), changing the 
#'   objective
#'
#' When scale_M = TRUE, the mu parameter no longer has its original 
#' mathematical meaning for balancing the two terms, as the relative scales 
#' have been artificially adjusted.
#'
#' For reproducibility across different systems, set n_cores = 1 to ensure 
#' deterministic results from PRIMME eigenvalue computations.
#'
#' **Handling NA Labels:**
#' Samples with NA labels are supported through the following mechanism:
#' - They contribute to the low-rank reconstruction term M = (I-R)ᵀ(I-R)
#' - They do not participate in similarity constraints (zero rows/columns in C)
#' - They are isolated nodes with zero degree; zero eigenvalues of Z (which
#'   occur only for μ = 1) are skipped
#' - Final embedding includes coordinates for all samples (labeled and unlabeled)
#'
#' Each domain needs at least two samples for its coordinates to carry
#' within-domain structure; a single-sample domain is embedded at the origin,
#' and an error is raised when no domain has two or more samples.
#'
#' @export
lowrank_align.hyperdesign <- function(data, y,
                                     preproc = center(),
                                     ncomp = 2,
                                     simfun,
                                     mu = .5,
                                     lambda = NULL,
                                     scale_M = FALSE,
                                     n_cores = NULL,
                                     sv_thresh = 1,
                                     solver = c("explicit", "operator"),
                                     ...) {

  solver <- match.arg(solver)
  if (!is.numeric(sv_thresh) || length(sv_thresh) != 1L || !is.finite(sv_thresh) ||
      sv_thresh < 0) {
    stop("`sv_thresh` must be a finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(mu) || length(mu) != 1L || !is.finite(mu) || mu < 0 || mu > 1) {
    stop("`mu` must be a single number in [0, 1].", call. = FALSE)
  }
  y <- rlang::enquo(y)

  # Extract labels, preserving NA values for semi-supervised learning
  labels <- unlist(purrr::map(data, function(x) {
    x$design |> select(!!y) |> pull(!!y)
  }))

  n_labeled <- sum(!is.na(labels))
  n_unlabeled <- sum(is.na(labels))

  if (n_labeled == 0) {
    stop("No labeled samples found. Low-rank alignment requires at least some labeled samples for similarity constraints.",
         call. = FALSE)
  }

  if (n_unlabeled > 0) {
    message("Semi-supervised low-rank alignment: ", n_labeled,
            " labeled samples, ", n_unlabeled, " unlabeled samples")
  }

  ninstances <- length(labels)

  if (solver == "explicit") {
    old_threads <- getOption("PRIMME_num_threads")
    if (!is.null(n_cores)) {
      options(PRIMME_num_threads = n_cores)
      on.exit(options(PRIMME_num_threads = old_threads), add = TRUE)
    }
  } else if (!is.null(n_cores)) {
    warning("n_cores is ignored for solver=\"operator\"; control BLAS threads instead.")
  }

  pdata <- data
  nb <- length(data)

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
    Xi <- as.matrix(data[[i]]$x)
    if (!is.null(preproc) && !is.null(preproc[[i]])) {
      pre_i <- preproc[[i]]
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }
      ft <- .preproc_fit_transform(pre_i, Xi)
      pdata[[i]]$x <- ft$transformed
      proclist[[i]] <- ft$preproc
    } else {
      # Record an identity transform; `proclist[[i]] <- NULL` would delete
      # the element and misalign the blocks.
      pdata[[i]]$x <- Xi
      proclist[i] <- list(.identity_preproc(Xi))
    }
  }
  names(proclist) <- names(pdata)

  sample_sizes <- vapply(pdata, function(x) nrow(x$x), integer(1))

  feature_sizes <- vapply(pdata, function(x) ncol(x$x), integer(1))
  feature_offsets <- cumsum(c(0L, feature_sizes))
  block_idx_list <- lapply(seq_along(feature_sizes), function(i) {
    if (feature_sizes[i] == 0L) {
      integer(0)
    } else {
      seq.int(feature_offsets[i] + 1L, feature_offsets[i + 1L])
    }
  })
  names(block_idx_list) <- names(pdata)

  proc <- multivarious::concat_pre_processors(proclist, block_idx_list)

  R_factors <- purrr::map(pdata, function(x) .compute_R_factors(x$x, sv_thresh = sv_thresh))

  # Use simfun to create similarity matrix (should handle NA labels)
  C_block <- simfun(labels)
  C_block <- Matrix::Matrix(C_block, sparse = TRUE)
  C_block[C_block < 0] <- 0
  C_block <- (C_block + Matrix::t(C_block)) / 2
  Matrix::diag(C_block) <- 0

  # Mask to cross-set edges only (Eq. 12): filter the stored entries by
  # domain id instead of materialising an n x n mask.
  block_id <- rep.int(seq_along(sample_sizes), sample_sizes)
  C_block <- .lowrank_cross_domain_edges(C_block, block_id)

  deg <- Matrix::rowSums(C_block)
  iso <- which(deg == 0)
  k_skip <- length(iso)

  if (k_skip > 0) {
    message("Detected ", k_skip, " isolated nodes (unlabeled samples); they are ",
            "embedded through the low-rank term only.")
  }

  Ls <- Matrix::Diagonal(x = deg) - C_block

  # The embedding is solved on the orthogonal complement of the per-domain
  # constant vectors (each domain's coordinates have zero mean). With centred
  # features, every domain constant satisfies R e_k = 0, so M e_k = e_k and the
  # global constant 1 is an exact eigenvector of Z (L 1 = 0) with eigenvalue
  # 1 - mu; without the constraint it is returned as a zero-variance component
  # whenever the label term dominates.
  n_free <- ninstances - length(unique(block_id))
  if (n_free < 1L) {
    stop("lowrank_align() needs at least one domain with two or more samples: ",
         "per-domain constant directions carry no within-domain structure.",
         call. = FALSE)
  }

  nev <- min(ncomp + k_skip + 4L, max(1L, n_free - 1L))

  if (solver == "explicit") {
    Rs <- purrr::map2(R_factors, sample_sizes, .factors_to_R)
    R_block <- Matrix::bdiag(Rs)

    R_i <- Matrix::Diagonal(ninstances) - R_block
    M <- crossprod(R_i, R_i)

    # lambda_max(M) <= 1 because 0 < a < 1; scaled below with M when scale_M.
    m_upper <- 1
    if (scale_M) {
      s_m <- .lowrank_eigs_sym(M, 1L, which = "LA")
      s_l <- .lowrank_eigs_sym(Ls, 1L, which = "LA")
      if (s_m$values[1] > 0 && is.finite(s_m$values[1]) &&
          s_l$values[1] > 0 && is.finite(s_l$values[1])) {
        e_ratio <- s_l$values[1] / s_m$values[1]
        M <- M * e_ratio
        m_upper <- e_ratio
        if (getOption("lowrank_align.verbose", FALSE)) {
          message("Applied M scaling with ratio: ", round(e_ratio, 4),
                  ". This changes the objective function - consider adjusting mu.")
        }
      }
    }

    Z <- (1 - mu) * M + (2 * mu * Ls)

    # Rigorous upper bounds on lambda_max(Z): lambda_max(M) <= m_upper and
    # lambda_max(L) <= 2 max(deg), or Gershgorin on Z itself. The constraint
    # directions are shifted just above the smaller bound.
    z_upper <- min((1 - mu) * m_upper + 4 * mu * max(deg),
                   max(Matrix::rowSums(abs(Z))))
    eig <- .lowrank_eigs_restricted(Z, nev, block_id, sigma = 1.01 * z_upper + 1e-8)
    evals <- eig$values
    vecs <- eig$vectors

  } else {
    if (!requireNamespace("RSpectra", quietly = TRUE)) {
      stop("RSpectra package is required for solver=\"operator\".", call. = FALSE)
    }

    M_op <- .make_M_op(R_factors, sample_sizes)
    lamM <- .lambda_max_power(M_op, n = ninstances, n_iter = 25L)
    lamM <- if (is.finite(lamM) && lamM > 0) lamM else 0

    if (Matrix::nnzero(Ls) > 0) {
      lamL <- RSpectra::eigs_sym(Ls, k = 1, which = "LM")$values[1]
    } else {
      lamL <- 0
    }
    lamL <- if (is.finite(lamL) && lamL > 0) lamL else 0

    if (scale_M && lamM > 0 && lamL > 0) {
      scale_ratio <- lamL / lamM
      orig_M_op <- M_op
      M_op <- function(x, ...) scale_ratio * orig_M_op(x, ...)
      lamM <- scale_ratio * lamM
      if (getOption("lowrank_align.verbose", FALSE)) {
        message("Applied M scaling with ratio: ", round(scale_ratio, 4),
                ". This changes the objective function - consider adjusting mu.")
      }
    }

    Z_op <- .make_Z_op(M_op, Ls, mu)
    c_bound <- (1 - mu) * lamM + 2 * mu * lamL

    if (!is.finite(c_bound) || c_bound <= 0) {
      c_bound <- max(1, lamM + lamL + 1e-6)
    }

    # P (I - Z / c) P with P = I - Q Q': the constraint directions get
    # eigenvalue 0 and the complement keeps 1 - lambda / c.
    P <- .lowrank_projector(block_id)
    B_shift <- .make_shifted_op(Z_op, ninstances, c_bound)
    B_op <- function(x, ...) P(B_shift(P(x)))

    eig <- RSpectra::eigs_sym(B_op, k = nev, n = ninstances, which = "LM")
    evals <- c_bound * (1 - eig$values)
    vecs <- eig$vectors
    # Discard any returned constraint direction (only possible when nev is
    # close to n_free) and remove round-off outside the complement.
    in_constraint <- colSums(vecs^2 - P(vecs)^2) > 0.5
    evals <- evals[!in_constraint]
    vecs <- P(vecs[, !in_constraint, drop = FALSE])
    vecs <- sweep(vecs, 2, sqrt(colSums(vecs^2)), "/")
  }

  keep <- which(evals > 1e-8)
  if (length(keep) < ncomp) {
    warning("Only ", length(keep), " non-zero eigenvalues found, but ",
            ncomp, " components requested. Using ", length(keep), " components.",
            call. = FALSE)
    ncomp <- length(keep)
  }

  keep <- keep[order(evals[keep])[seq_len(ncomp)]]
  v <- vecs[, keep, drop = FALSE]

  Xp <- Matrix::bdiag(lapply(pdata, function(x) x$x))

  if (is.null(lambda)) {
    cvfit <- glmnet::cv.glmnet(Xp, v, family = "mgaussian", alpha = 0,
                               intercept = FALSE)
    lambda <- cvfit$lambda.min
  }

  rfit <- glmnet::glmnet(Xp, v, family = "mgaussian", alpha = 0,
                          lambda = lambda, intercept = FALSE)

  coef_list <- stats::coef(rfit, s = lambda)
  cfs <- do.call(cbind, coef_list)
  cfs <- cfs[-1, , drop = FALSE]

  multivarious::multiblock_biprojector(
    v = cfs,
    s = v,
    sdev = apply(v, 2, sd),
    preproc = proc,
    block_indices = block_idx_list,
    labels = labels,
    mu = mu,
    classes = "lowrank_align"
  )
}
