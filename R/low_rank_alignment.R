#' @noRd
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
  keep <- svd_X$d > sv_thresh

  if (!any(keep)) {
    return(list(U = matrix(0, nrow(X), 0),
                a = numeric(0)))
  }

  list(U = svd_X$v[, keep, drop = FALSE],
       a = 1 - 1 / (svd_X$d[keep]^2))
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
  set.seed(seed)
  v <- rnorm(n)
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

#' @noRd
.make_Z_op <- function(M_op, Ls, mu) {
  function(x, ...) {
    (1 - mu) * M_op(x) + 2 * mu * as.numeric(Ls %*% x)
  }
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
#' Low-rank alignment optimizes the objective function Z = (1-μ) * M + 2μ * L 
#' where M captures low-rank structure and L is the graph Laplacian from 
#' similarity matrix. The method balances preserving low-rank structure 
#' (μ=0) with enforcing similarity constraints (μ=1).
#' 
#' **Semi-supervised Learning Support:**
#' The algorithm handles NA labels gracefully. Unlabeled samples:
#' - Still contribute to the low-rank structure term M through their data
#' - Do not participate in the similarity constraints (L term)
#' - Receive coordinates in the joint embedding space
#' - Create isolated nodes that produce zero eigenvalues (automatically 
#'   skipped)
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
#' @param preproc Preprocessing function (default: center())
#' @param ncomp Number of components to extract (default: 2)
#' @param simfun Function to compute similarity matrix from labels. Should 
#'   handle NA labels gracefully (e.g., created with createSimFun())
#' @param mu Balance parameter between low-rank (μ=0) and similarity (μ=1) 
#'   terms (default: 0.5)
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
#'   Values at or below the threshold are discarded, matching Eq. 11 of the
#'   original Low-Rank Alignment paper.
#' @param solver Eigen solver backend. `"explicit"` (default) forms the dense
#'   matrices and uses PRIMME, matching the original implementation. `"operator"`
#'   keeps low-rank factors and uses RSpectra with a matrix-vector operator to
#'   reduce memory and runtime.
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
#' - They create isolated nodes with zero degree, producing zero eigenvalues
#' - The algorithm automatically detects and skips these zero modes
#' - Final embedding includes coordinates for all samples (labeled and unlabeled)
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
  y <- rlang::enquo(y)

  # Extract labels, preserving NA values for semi-supervised learning
  labels <- unlist(purrr::map(data, function(x) {
    x$design %>% select(!!y) %>% pull(!!y)
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

  pdata <- multivarious::init_transform(data, preproc)
  proclist <- attr(pdata, "preproc")
  names(proclist) <- names(pdata)

  sample_sizes <- vapply(pdata, function(x) nrow(x$x), integer(1))
  sample_offsets <- cumsum(c(0L, sample_sizes))

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

  # Mask to cross-set edges only (Eq. 12)
  if (length(sample_sizes) > 1L) {
    C_mask <- Matrix::Matrix(0, ninstances, ninstances, sparse = TRUE)
    for (a in seq_along(sample_sizes)) {
      ia <- (sample_offsets[a] + 1L):sample_offsets[a + 1L]
      if (length(ia) == 0L) next
      for (b in seq_along(sample_sizes)) {
        if (a == b) next
        jb <- (sample_offsets[b] + 1L):sample_offsets[b + 1L]
        if (length(jb) == 0L) next
        C_mask[ia, jb] <- 1
      }
    }
    C_block <- C_block * C_mask
  } else {
    C_block <- Matrix::Matrix(0, ninstances, ninstances, sparse = TRUE)
  }

  deg <- Matrix::rowSums(C_block)
  iso <- which(deg == 0)
  k_skip <- length(iso)

  if (k_skip > 0) {
    message("Detected ", k_skip, " isolated nodes (unlabeled samples). Will skip corresponding zero eigenvalue modes.")
  }

  Ls <- Matrix::Diagonal(x = deg) - C_block

  nev <- min(ncomp + k_skip + 4L, max(1L, ninstances - 2L))

  if (solver == "explicit") {
    Rs <- purrr::map2(R_factors, sample_sizes, .factors_to_R)
    R_block <- Matrix::bdiag(Rs)

    R_i <- Matrix::Diagonal(ninstances) - R_block
    M <- crossprod(R_i, R_i)

    if (scale_M) {
      s_m <- PRIMME::eigs_sym(M, NEig = 1, which = "LA",
                              method = "PRIMME_DEFAULT_MIN_MATVECS")
      s_l <- PRIMME::eigs_sym(Ls, NEig = 1, which = "LA",
                              method = "PRIMME_DEFAULT_MIN_MATVECS")
      if (s_m$values[1] > 0 && is.finite(s_m$values[1]) &&
          s_l$values[1] > 0 && is.finite(s_l$values[1])) {
        e_ratio <- s_l$values[1] / s_m$values[1]
        M <- M * e_ratio
        if (getOption("lowrank_align.verbose", FALSE)) {
          message("Applied M scaling with ratio: ", round(e_ratio, 4),
                  ". This changes the objective function - consider adjusting mu.")
        }
      }
    }

    Z <- (1 - mu) * M + (2 * mu * Ls)

    eig <- PRIMME::eigs_sym(Z, NEig = nev, which = "SA",
                             method = "PRIMME_DEFAULT_MIN_MATVECS")
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

    B_op <- .make_shifted_op(Z_op, ninstances, c_bound)

    eig <- RSpectra::eigs_sym(B_op, k = nev, n = ninstances, which = "LM")
    evals <- c_bound * (1 - eig$values)
    vecs <- eig$vectors
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
