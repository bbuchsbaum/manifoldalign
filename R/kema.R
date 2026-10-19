#' @importFrom Matrix bdiag Diagonal rowSums sparseMatrix
#' @importFrom kernlab kernelMatrix
#' @importFrom cluster pam
#' @importFrom chk chk_number chk_range chk_logical chk_true
#' @importFrom rlang enquo !!
#' @importFrom purrr map
#' @importFrom coop cosine
#' @importFrom methods new is as
#' @import Matrix
#' @keywords internal
#' @noRd
NULL

#' Block indices for data list
#'
#' @param data_list List of data blocks, each with an `x` component whose rows
#'   are samples. Dense matrices, sparse `Matrix` objects and vectors are
#'   accepted; a vector counts as one sample per element.
#' @return Matrix with start and end indices for each block
#' @examples
#' data1 <- list(x = matrix(rnorm(20), 10, 2))
#' data2 <- list(x = matrix(rnorm(30), 15, 2))
#' data_list <- list(data1, data2)
#' indices <- block_indices(data_list)
#' indices
#' @export
block_indices <- function(data_list) {
  # Generate block indices for a list of data blocks
  if (is.list(data_list)) {
    n_samples <- vapply(data_list, function(x) {
      if (!is.null(dim(x$x))) nrow(x$x) else length(x$x)
    }, numeric(1))
  } else {
    stop("data_list must be a list")
  }

  end_indices <- cumsum(n_samples)
  start_indices <- c(1, end_indices[-length(end_indices)] + 1)

  cbind(start=start_indices, end=end_indices)
}

#' Cosine kernel function
#'
#' @return A kernel object for cosine similarity
#' @examples
#' kern <- coskern()
#' x <- c(1, 2, 3)
#' y <- c(4, 5, 6)
#' # kern(x, y)  # Compute cosine similarity
#' @export
coskern <- function() {
  rval <- function(x, y = NULL) {
    if (!is(x, "vector"))
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y))
      stop("y must be a vector")
    if (is.null(y)) {
      return(coop::cosine(x, x))
    }
    if (!length(x) == length(y))
      stop("number of dimension must be the same on both data points")
    coop::cosine(x, y)
  }

  return(new("kernel", .Data = rval, kpar = list()))
}

# ---------------------------------------------------------------------------
# KEMA kernels, landmarks and the dense-memory guard
# ---------------------------------------------------------------------------

#' Select REKEMA landmark rows for each domain
#'
#' @param n_per_block Integer vector of domain sizes.
#' @param sample_frac Landmark fraction; `1` selects every row (full KEMA).
#' @param seed `NULL` to draw from the session RNG, otherwise a seed used
#'   locally (the caller's RNG stream is restored).
#' @return A list with one sorted integer vector of landmark rows per domain,
#'   or `NULL` for full KEMA.
#' @keywords internal
#' @noRd
.kema_select_landmarks <- function(n_per_block, sample_frac, seed = NULL) {
  if (sample_frac >= 1) {
    return(NULL)
  }
  draw <- function() {
    lapply(n_per_block, function(n_i) {
      n_landmarks <- max(1L, as.integer(round(sample_frac * n_i)))
      if (n_landmarks >= n_i) {
        seq_len(n_i)
      } else {
        sort(sample.int(n_i, n_landmarks, replace = FALSE))
      }
    })
  }
  if (is.null(seed)) draw() else .with_local_seed(seed, draw())
}

#' Approximate peak memory of the dense parts of a KEMA fit
#'
#' Every KEMA backend materialises the dense kernel blocks (`n_i x r_i`, held
#' in dense and sparse-format copies), their eigen/SVD factors and dense
#' label-comparison matrices over the labelled rows. The dense quotient
#' backends (`full_exact`, `reduced_exact`) additionally form the two
#' `r_total x r_total` quotient matrices plus a Cholesky preconditioner or
#' dense eigensolver workspace of the same order. `operator_exact` applies the
#' quotient matrix-free: beyond the kernel factors it keeps only per-domain
#' `r_i x r_i` preconditioner blocks and `O((n_total + r_total) * block)`
#' LOBPCG workspace. The estimate is a coarse upper bound on those
#' allocations: about 60 bytes per kernel entry, 96 bytes per entry of an
#' `r_total x r_total` matrix for the dense quotient backends (40 bytes per
#' entry of each `r_i x r_i` block for `operator_exact`), and 32 bytes per pair
#' of labelled rows. Calibrated against the peak resident memory of complete
#' fits on the RAPID scaling fixture (R 4.5.3, eigencore 2.0.0, single-threaded
#' OpenBLAS): at 11,000 samples `operator_exact` peaked at 4.8 GB (estimate
#' 6.2 GB; the dense backends are estimated at 15 GB); at 5,500 samples
#' 0.93 GB (`operator_exact`, estimate 1.6 GB) and 2.9 GB (`reduced_exact`,
#' estimate 3.8 GB); at 2,200 samples 0.51 GB for `full_exact` (estimate
#' 0.61 GB).
#'
#' @param backend Backend whose footprint is estimated.
#' @param ncomp Requested number of components (sets the LOBPCG block width).
#' @keywords internal
#' @noRd
.kema_dense_memory_estimate <- function(n_per_block, r_per_block, n_labeled,
                                        backend = "full_exact",
                                        ncomp = 2L) {
  n_per_block <- as.numeric(n_per_block)
  r_per_block <- as.numeric(r_per_block)
  r_total <- sum(r_per_block)
  kernel_bytes <- 60 * sum(n_per_block * r_per_block)
  label_bytes <- 32 * as.numeric(n_labeled)^2
  quotient_bytes <- if (identical(backend, "operator_exact")) {
    block <- max(2 * ncomp, ncomp + 4)
    40 * sum(r_per_block^2) +
      8 * 12 * block * (sum(n_per_block) + r_total)
  } else {
    96 * r_total^2
  }
  kernel_bytes + quotient_bytes + label_bytes
}

#' Apply the KEMA dense-memory guard to the backend candidates
#'
#' @return A list with the per-candidate `estimates` and the `admissible`
#'   candidates (in their original order). Errors when no candidate fits.
#' @keywords internal
#' @noRd
.kema_check_dense_memory <- function(n_per_block, r_per_block, n_labeled,
                                     max_bytes,
                                     candidates = "full_exact",
                                     ncomp = 2L) {
  est <- vapply(candidates, function(b) {
    .kema_dense_memory_estimate(
      n_per_block, r_per_block, n_labeled,
      backend = b, ncomp = ncomp
    )
  }, numeric(1))
  names(est) <- candidates
  fits <- !is.finite(max_bytes) | est <= max_bytes
  if (!any(fits)) {
    cheapest <- candidates[which.min(est)]
    fmt <- function(x) {
      format(structure(x, class = "object_size"), units = "auto")
    }
    stop(
      "KEMA would allocate dense O(n^2) objects of roughly ",
      fmt(min(est)), " for its cheapest candidate backend '", cheapest,
      "' (", sum(n_per_block), " samples, ", sum(r_per_block),
      " kernel columns, ", n_labeled, " labelled rows), which exceeds ",
      "backend_control$max_dense_bytes = ", fmt(max_bytes), ". Every KEMA ",
      "backend forms the dense kernel blocks and their eigen factors ",
      "(O(sum n_i r_i)); 'full_exact' and 'reduced_exact' also form dense ",
      "r_total x r_total quotient matrices, which 'operator_exact' applies ",
      "matrix-free. Use sample_frac < 1 (REKEMA) to reduce the kernel rank, ",
      "or raise max_dense_bytes if the memory is available.",
      call. = FALSE
    )
  }
  list(estimates = est, admissible = candidates[fits])
}


#' Compute the KEMA kernel blocks
#'
#' @param strata Preprocessed domains (each with a dense `x`).
#' @param kernel A kernlab kernel.
#' @param landmarks `NULL` for full KEMA (square blocks `K(X_i, X_i)`), or a
#'   list of landmark row indices per domain for REKEMA (`K(X_i, X_i[L_i, ])`).
#' @return A list of `dgCMatrix` kernel blocks (dense content, no explicit
#'   zeros); square full-KEMA blocks are exactly symmetric.
#' @keywords internal
#' @noRd
compute_kernels <- function(strata, kernel, landmarks = NULL) {
  Ks <- vector("list", length(strata))

  for (i in seq_along(strata)) {
    X <- strata[[i]]$x

    if (!is.matrix(X) || nrow(X) == 0 || ncol(X) == 0) {
      stop("Data block ", i, " is not a valid matrix or is empty", call. = FALSE)
    }

    if (!is.null(landmarks)) {
      X_landmarks <- X[landmarks[[i]], , drop = FALSE]
      K <- safe_compute(
        kernlab::kernelMatrix(kernel, X, X_landmarks),
        paste0("Kernel computation failed for data block ", i, " in REKEMA mode")
      )
    } else {
      K <- safe_compute(
        kernlab::kernelMatrix(kernel, X),
        paste0("Full kernel computation failed for data block ", i)
      )
    }

    if (any(!is.finite(K))) {
      stop("Kernel matrix for block ", i, " contains non-finite values", call. = FALSE)
    }

    K <- as.matrix(K)
    if (is.null(landmarks)) {
      # Square full-KEMA blocks are symmetric in exact arithmetic; remove the
      # rounding asymmetry on the dense content before the single conversion.
      K <- (K + t(K)) / 2
    }
    Ks[[i]] <- as(K, "dgCMatrix")
    rm(K)
  }

  Ks
}

# ---------------------------------------------------------------------------
# KEMA kernel encoder (out-of-sample prediction)
# ---------------------------------------------------------------------------

#' Build the kernel encoder of a KEMA fit
#'
#' KEMA's embedding of a row `x` of domain `i` is `k(x, R_i) alpha_i`, where
#' `R_i` are the reference rows of the domain (every preprocessed training
#' row for full KEMA, the stored landmark rows for REKEMA) and `alpha_i` is the
#' domain's block of dual coefficients. Training scores are exactly this map
#' applied to the training rows.
#'
#' @keywords internal
#' @noRd
.kema_build_encoder <- function(pdata, proclist, kernel, alpha, landmarks,
                                domain_names) {
  n_ref <- if (is.null(landmarks)) {
    vapply(pdata, function(d) nrow(d$x), integer(1))
  } else {
    vapply(landmarks, length, integer(1))
  }
  ends <- cumsum(n_ref)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  alpha <- as.matrix(alpha)
  if (nrow(alpha) != sum(n_ref)) {
    stop("KEMA dual coefficients do not match the reference rows.", call. = FALSE)
  }

  reference <- lapply(seq_along(pdata), function(i) {
    X <- as.matrix(pdata[[i]]$x)
    if (is.null(landmarks)) X else X[landmarks[[i]], , drop = FALSE]
  })
  alpha_blocks <- lapply(seq_along(pdata), function(i) {
    alpha[seq.int(starts[i], ends[i]), , drop = FALSE]
  })
  names(reference) <- domain_names
  names(alpha_blocks) <- domain_names
  if (!is.null(landmarks)) names(landmarks) <- domain_names

  linear <- inherits(kernel, "vanillakernel")
  structure(
    list(
      kernel = kernel,
      reference = reference,
      alpha = alpha_blocks,
      landmarks = landmarks,
      preproc = proclist,
      domain_names = domain_names,
      n_features = vapply(pdata, function(d) ncol(d$x), integer(1)),
      form = if (is.null(landmarks)) "full" else "reduced",
      linear_kernel = linear
    ),
    class = "kema_encoder"
  )
}

#' Exact linear loadings of a KEMA fit, when they exist
#'
#' For the linear kernel `k(x, r) = x'r`, `k(x, R_i) alpha_i = x R_i' alpha_i`,
#' so `R_i' alpha_i` are exact feature loadings. For any other kernel no
#' feature-space loading reproduces the encoder, so the loadings are `NA`.
#'
#' @keywords internal
#' @noRd
.kema_encoder_loadings <- function(encoder) {
  ncomp <- ncol(encoder$alpha[[1L]])
  blocks <- lapply(seq_along(encoder$alpha), function(i) {
    if (isTRUE(encoder$linear_kernel)) {
      crossprod(encoder$reference[[i]], encoder$alpha[[i]])
    } else {
      matrix(NA_real_, encoder$n_features[[i]], ncomp)
    }
  })
  do.call(rbind, blocks)
}

#' @keywords internal
#' @noRd
.kema_block_index <- function(encoder, side) {
  n_blocks <- length(encoder$alpha)
  if (length(side) != 1L || is.na(side)) {
    stop("`side` must identify exactly one KEMA domain.", call. = FALSE)
  }
  if (is.character(side)) {
    idx <- match(side, encoder$domain_names)
    if (is.na(idx) && side %in% c("i", "j")) {
      idx <- if (identical(side, "i")) 1L else 2L
    }
    if (is.na(idx) || idx > n_blocks) {
      stop("Unknown KEMA domain `", side, "`; expected one of: ",
           paste(encoder$domain_names, collapse = ", "), ".", call. = FALSE)
    }
    return(as.integer(idx))
  }
  if (!is.numeric(side) || abs(side - round(side)) > 1e-8 ||
      side < 1 || side > n_blocks) {
    stop("`side` must be a KEMA domain name or an index in 1..", n_blocks, ".",
         call. = FALSE)
  }
  as.integer(side)
}

#' @keywords internal
#' @noRd
.kema_apply_block_preproc <- function(pre, X) {
  if (is.null(pre)) {
    return(X)
  }
  if (is.function(pre)) {
    return(as.matrix(pre(X)))
  }
  as.matrix(multivarious::transform(pre, X))
}

#' @keywords internal
#' @noRd
.kema_encode <- function(encoder, newX, block) {
  X <- newX
  if (is.data.frame(X) || methods::is(X, "Matrix")) {
    X <- as.matrix(X)
  } else if (is.null(dim(X))) {
    X <- matrix(X, nrow = 1L)
  }
  if (!is.numeric(X)) {
    stop("`newX` must be numeric.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  if (ncol(X) != encoder$n_features[[block]]) {
    stop("`newX` has ", ncol(X), " columns; KEMA domain `",
         encoder$domain_names[[block]], "` was fitted with ",
         encoder$n_features[[block]], ".", call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("`newX` contains non-finite values.", call. = FALSE)
  }
  Xp <- .kema_apply_block_preproc(encoder$preproc[[block]], X)
  K <- kernlab::kernelMatrix(encoder$kernel, Xp, encoder$reference[[block]])
  out <- as.matrix(K) %*% encoder$alpha[[block]]
  out <- matrix(as.numeric(out), nrow = nrow(X))
  rownames(out) <- rownames(X)
  out
}

#' Encode new observations with a fitted KEMA model
#'
#' KEMA is a kernel method: the embedding of a new row `x` of domain `i` is
#' `k(phi_i(x), R_i) alpha_i`, where `phi_i` is the domain's fitted
#' preprocessing, `R_i` are the domain's reference rows (all preprocessed
#' training rows for full KEMA, the stored landmark rows for REKEMA,
#' `sample_frac < 1`) and `alpha_i` are the domain's dual coefficients.
#' Re-encoding a domain's training rows reproduces the stored training scores,
#' and each row is encoded independently of the other rows in `newX`.
#'
#' @param fit_or_transform A KEMA fit returned by [kema()].
#' @param newX New observations (rows) with the same columns as the domain's
#'   training data. A numeric vector is treated as a single row.
#' @param side Domain name or index. For two-domain fits `"i"` and `"j"` are
#'   accepted as aliases of the first and second domain.
#' @param ... Not used; supplying extra arguments is an error.
#' @return A numeric matrix with one row per row of `newX` and `ncomp`
#'   columns.
#' @export
oos_predict.kema_fit <- function(fit_or_transform, newX, side, ...) {
  if (length(list(...))) {
    stop("oos_predict() for KEMA fits accepts no additional arguments.",
         call. = FALSE)
  }
  encoder <- fit_or_transform$encoder
  if (!inherits(encoder, "kema_encoder")) {
    stop("This KEMA fit has no kernel encoder; refit it with kema().",
         call. = FALSE)
  }
  block <- .kema_block_index(encoder, side)
  .kema_encode(encoder, newX, block)
}

#' Projection methods for KEMA fits
#'
#' A KEMA fit is a per-domain kernel encoder, not a linear projector on the
#' concatenated feature space. `project_block()` encodes rows of one domain
#' with the fitted kernel map (see \code{\link[=oos_predict.kema_fit]{oos_predict.kema_fit()}}); `project()` and
#' `partial_project()` are not defined for KEMA and error rather than return
#' a linear surrogate. The stored loadings `fit$v` are exact only for the
#' linear kernel (`kernlab::vanilladot()`); for every other kernel they are
#' `NA`.
#'
#' @param x A KEMA fit.
#' @param new_data New observations of one domain.
#' @param block Domain name or index.
#' @param least_squares Must be `FALSE`; least-squares partial projection is
#'   not defined for a kernel encoder.
#' @param ... Not used.
#' @return For `project_block()`, a numeric matrix of embeddings.
#' @name kema_projection
#' @exportS3Method multivarious::project_block
project_block.kema_fit <- function(x, new_data, block, least_squares = FALSE, ...) {
  if (isTRUE(least_squares)) {
    stop("least_squares projection is not defined for KEMA kernel encoders.",
         call. = FALSE)
  }
  oos_predict.kema_fit(x, new_data, side = block)
}

#' @rdname kema_projection
#' @exportS3Method multivarious::project
project.kema_fit <- function(x, new_data, ...) {
  stop(
    "project() is not defined for KEMA fits: KEMA encodes each domain with ",
    "its own kernel map k(x, R_i) alpha_i, which has no linear loading on the ",
    "concatenated features. Use oos_predict(fit, newX, side = <domain>) or ",
    "project_block(fit, new_data, block = <domain>).",
    call. = FALSE
  )
}

#' @rdname kema_projection
#' @exportS3Method multivarious::partial_project
partial_project.kema_fit <- function(x, new_data, ...) {
  stop(
    "partial_project() is not defined for KEMA kernel encoders. Use ",
    "oos_predict(fit, newX, side = <domain>) to encode complete rows of one domain.",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.reject_ignored_kema_arguments <- function(explicit_legacy, dots) {
  rejected <- names(explicit_legacy)[explicit_legacy]
  if (length(dots)) {
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep("...", length(dots))
    }
    dot_names[is.na(dot_names) | !nzchar(dot_names)] <- "..."
    rejected <- c(rejected, dot_names)
  }
  rejected <- unique(rejected)
  if (length(rejected)) {
    stop(
      "Unsupported KEMA argument(s): ", paste(rejected, collapse = ", "),
      ". These compatibility controls were withdrawn because the current ",
      "implementation did not use them; silently accepting them could ",
      "misrepresent the fitted objective.",
      call. = FALSE
    )
  }
}

#' @keywords internal
#' @noRd
.validate_kema_solver <- function(solver) {
  if (length(solver) != 1L || !is.character(solver) || is.na(solver)) {
    stop("solver must be the single value 'exact'.", call. = FALSE)
  }
  if (identical(solver, "regression")) {
    stop(
      "solver='regression' has been withdrawn: it never selected a distinct ",
      "regression implementation and silently ran the exact KEMA solver. ",
      "Use solver='exact' or omit solver.",
      call. = FALSE
    )
  }
  if (!identical(solver, "exact")) {
    stop("solver must be 'exact'.", call. = FALSE)
  }
  invisible(solver)
}

#' Kernel Manifold Alignment for Multidesign Data
#'
#' Performs Kernel Manifold Alignment on multidesign data structures.
#' Automatically splits data by subject variable and aligns domains.
#'
#' @param data A multidesign object containing data with subject groupings
#' @param y Name of the label variable to use for alignment (can contain NA
#'   for unlabeled samples)
#' @param subject Name of the subject variable that defines the domains/strata
#' @param preproc Preprocessing function to apply to the data (default:
#'   center())
#' @param ncomp Number of components to extract (default: 2)
#' @param knn Number of nearest neighbors for graph construction (default: 5)
#' @param sigma Scale parameter with two roles. It is always the heat-kernel
#'   bandwidth of the within-domain kNN graph
#'   (`adjoin::graph_weights(..., sigma = sigma)`). When `kernel = NULL` it is
#'   also the inverse width of the default kernel `kernlab::rbfdot(sigma)`,
#'   `k(x, y) = exp(-sigma * ||x - y||^2)`; these are different
#'   parameterisations, so pass `kernel` explicitly to set the kernel width
#'   independently of the graph. Default 0.73 (the hyperdesign method's `NULL`
#'   resolves to 0.73).
#' @param u Trade-off parameter between data geometry and class alignment
#'   (0-1, default: 0.5)
#' @param kernel A kernlab kernel object. The defaults differ by method: the
#'   multidesign method defaults to `coskern()`; the hyperdesign method's
#'   default `NULL` means `kernlab::rbfdot(sigma = sigma)`.
#' @param sample_frac Fraction of each domain's rows kept as REKEMA landmarks
#'   (default 1, full KEMA). The selected landmark rows are stored in
#'   `fit$encoder$landmarks`.
#' @param seed Optional seed for REKEMA landmark selection. `NULL` (default)
#'   draws landmarks from the session RNG; a value draws them reproducibly
#'   without disturbing the caller's RNG stream. No effect when
#'   `sample_frac = 1`.
#' @param use_laplacian Withdrawn compatibility argument. Supplying it now
#'   errors because previous releases silently ignored it.
#' @param solver Eigensolver contract. Only `"exact"` is supported. The former
#'   `"regression"` value never selected a distinct implementation and is now
#'   rejected.
#' @param backend Backend for the original eigensolver. One of `"auto"`,
#'   `"full_exact"`, `"reduced_exact"`, or `"operator_exact"`. All backends
#'   solve the same nullspace-deflated quotient pencil. `"full_exact"` and
#'   `"reduced_exact"` form the dense r-by-r quotient matrices (r the total
#'   kernel rank) and use a dense generalized solve for small quotients,
#'   otherwise preconditioned LOBPCG with a certified dense fallback.
#'   `"operator_exact"` applies the quotient matrix-free through the
#'   per-domain kernel factors (certified dense fallback only for very small
#'   quotients). Every backend still materialises the dense kernel blocks and
#'   their eigen factors; a memory guard (`backend_control$max_dense_bytes`)
#'   estimates each backend's footprint before allocation, skips candidates
#'   that exceed it (recorded in `fit$fidelity_history`) and errors when none
#'   fits.
#' @param backend_control Optional list controlling auto backend thresholds,
#'   fidelity checks and the dense-memory guard (passed through to
#'   `kema_orig()`).
#' @param dweight,rweight,simfun,disfun Withdrawn compatibility arguments.
#'   Supplying any of them errors because they were previously ignored.
#' @param lambda Strictly positive Tikhonov regularization for the generalized
#'   metric (default: 0.0001).
#' @param centre_kernel Withdrawn compatibility argument. Supplying it errors.
#' @param ... No additional arguments are accepted; unknown arguments error.
#'
#' @return A KEMA fit of class `c("kema_fit", "kema_orig",
#'   "multiblock_biprojector", ...)`. `fit$s` holds the training scores,
#'   `fit$alpha` the dual coefficients and `fit$encoder` the kernel encoder
#'   (kernel, preprocessing, reference or landmark rows, per-domain dual
#'   coefficients). Use [oos_predict()] (or `project_block()`) to embed new
#'   rows; generic `project()` errors because KEMA has no linear loading.
#'   `fit$v` is exact only for the linear kernel and is `NA` otherwise.
#'
#' @examples
#' \donttest{
#' # Example with multidesign data
#' library(multidesign)
#'
#' # Create synthetic multi-subject data
#' set.seed(123)
#' data_design <- expand.grid(
#'   subject = factor(1:4),
#'   condition = factor(c("A", "B")),
#'   trial = 1:10
#' )
#'
#' # Generate synthetic data matrix
#' n_obs <- nrow(data_design)
#' n_features <- 20
#' X <- matrix(rnorm(n_obs * n_features), n_obs, n_features)
#'
#' # Create multidesign object
#' md <- multidesign(X, data_design)
#'
#' # Run KEMA alignment across subjects
#' result <- kema(md, y = condition, subject = subject, ncomp = 2)
#'
#' # Semi-supervised learning with missing labels
#' data_design$condition[sample(nrow(data_design), 20)] <- NA
#' md_semi <- multidesign(X, data_design)
#' result_semi <- kema(md_semi, y = condition, subject = subject, ncomp = 2)
#' }
#'
#' @rdname kema
#' @method kema multidesign
#' @export
#' @import Matrix
kema.multidesign <- function(data, y,
                             subject,
                             preproc=center(),
                             ncomp=2,
                             knn=5,
                             sigma=.73, u=.5,
                             kernel=coskern(),
                             sample_frac=1,
                             use_laplacian=TRUE,
                             solver="exact",
                             backend="auto",
                             backend_control=NULL,
                             dweight=.1,
                             rweight=0,
                             simfun=adjoin::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE,
                             seed=NULL,
                             ...) {

  subject <- rlang::enquo(subject)
  y <- rlang::enquo(y)

  .reject_ignored_kema_arguments(
    c(
      use_laplacian = !missing(use_laplacian),
      dweight = !missing(dweight),
      rweight = !missing(rweight),
      simfun = !missing(simfun),
      disfun = !missing(disfun),
      centre_kernel = !missing(centre_kernel)
    ),
    list(...)
  )
  .validate_kema_solver(solver)

  if (is.null(sigma)) {
    sigma <- 0.73
  }
  if (is.null(kernel)) {
    kernel <- kernlab::rbfdot(sigma = sigma)
  }

  result <- kema_orig.multidesign(
    data = data,
    y = !!y,
    subject = !!subject,
    preproc = preproc,
    ncomp = ncomp,
    knn = knn,
    sigma = sigma,
    mu = u,
    kernel = kernel,
    sample_frac = sample_frac,
    lambda = lambda,
    backend = backend,
    backend_control = backend_control,
    seed = seed
  )

  # Preserve legacy class tag for downstream code using `kema()`.
  result$classes <- "kema"
  if (is.numeric(result$eigenvalues)) {
    result$eigenvalues <- list(
      values = as.numeric(result$eigenvalues),
      all_values = as.numeric(result$eigenvalues),
      solver = "exact"
    )
  }
  result$regression_quality <- NULL
  result$retry_info <- NULL
  result
}

#' Kernel Manifold Alignment (KEMA) for Hyperdesign Data
#'
#' Performs Kernel Manifold Alignment on hyperdesign data structures.
#' Projects data from multiple domains into a shared latent space while
#' preserving manifold structure and aligning same-class samples.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the label variable to use for alignment (can contain NA
#'   for unlabeled samples)
#' @param preproc Preprocessing function to apply to the data (default:
#'   center())
#' @param ncomp Number of components to extract (default: 2)
#' @param knn Number of nearest neighbors for graph construction (default: 5)
#' @param sigma Scale parameter with two roles. It is always the heat-kernel
#'   bandwidth of the within-domain kNN graph
#'   (`adjoin::graph_weights(..., sigma = sigma)`). When `kernel = NULL` it is
#'   also the inverse width of the default kernel `kernlab::rbfdot(sigma)`,
#'   `k(x, y) = exp(-sigma * ||x - y||^2)`; these are different
#'   parameterisations, so pass `kernel` explicitly to set the kernel width
#'   independently of the graph. Default 0.73 (the hyperdesign method's `NULL`
#'   resolves to 0.73).
#' @param u Trade-off parameter between data geometry and class alignment
#'   (0-1, default: 0.5)
#' @param kernel A kernlab kernel object. The defaults differ by method: the
#'   multidesign method defaults to `coskern()`; the hyperdesign method's
#'   default `NULL` means `kernlab::rbfdot(sigma = sigma)`.
#' @param sample_frac Fraction of each domain's rows kept as REKEMA landmarks
#'   (default 1, full KEMA). The selected landmark rows are stored in
#'   `fit$encoder$landmarks`.
#' @param seed Optional seed for REKEMA landmark selection. `NULL` (default)
#'   draws landmarks from the session RNG; a value draws them reproducibly
#'   without disturbing the caller's RNG stream. No effect when
#'   `sample_frac = 1`.
#' @param use_laplacian Withdrawn compatibility argument. Supplying it now
#'   errors because previous releases silently ignored it.
#' @param solver Eigensolver contract. Only `"exact"` is supported. The former
#'   `"regression"` value never selected a distinct implementation and is now
#'   rejected.
#' @param dweight,rweight,simfun,disfun Withdrawn compatibility arguments.
#'   Supplying any of them errors because they were previously ignored.
#' @param lambda Strictly positive Tikhonov regularization for the generalized
#'   metric (default: 0.0001).
#' @param centre_kernel Withdrawn compatibility argument. Supplying it errors.
#' @param ... No additional arguments are accepted; unknown arguments error.
#'
#' @return A KEMA fit of class `c("kema_fit", "kema_orig",
#'   "multiblock_biprojector", ...)`; see the multidesign method.
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' # Create synthetic multi-domain data
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
#' )
#' hd <- structure(list(domain1 = domain1, domain2 = domain2), class = "hyperdesign")
#'
#' # Run KEMA with default settings
#' result <- kema(hd, y = labels, ncomp = 2)
#'
#' # Embed new rows of domain 2 with the fitted kernel encoder
#' new_rows <- matrix(rnorm(10), 5, 2)
#' oos_predict(result, new_rows, side = "domain2")
#'
#' # Semi-supervised learning with missing labels
#' hd_semi <- hd
#' hd_semi$domain1$design$labels[1:10] <- NA  # Mark some samples as unlabeled
#' result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#'
#' # Use exact solver for highest accuracy
#' result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#'
#' # Use REKEMA for large datasets (reproducible landmarks)
#' result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2, seed = 1)
#' }
#'
#' @details
#' `kema()` now delegates to the paper-faithful `kema_orig()` backend and solves
#' the original generalized eigenproblems from Tuia & Camps-Valls (2016),
#' including the reduced-rank REKEMA form when `sample_frac < 1`.
#'
#' Legacy extension arguments remain temporarily in the function signature so
#' callers receive an actionable error. They are not silently ignored.
#'
#' @references
#' Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain
#' adaptation. PLoS ONE, 11(2), e0148655.
#'
#' @rdname kema
#' @method kema hyperdesign
#' @export
kema.hyperdesign <- function(data, y,
                             preproc=center(),
                             ncomp=2,
                             knn=5,
                             sigma=NULL,
                             u=.5,
                             kernel=NULL,
                             sample_frac=1,
                             use_laplacian=TRUE,
                             solver="exact",
                             backend="auto",
                             backend_control=NULL,
                             dweight=.1,
                             rweight=0,
                             simfun=adjoin::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE,
                             seed=NULL,
                             ...) {

  y <- rlang::enquo(y)

  .reject_ignored_kema_arguments(
    c(
      use_laplacian = !missing(use_laplacian),
      dweight = !missing(dweight),
      rweight = !missing(rweight),
      simfun = !missing(simfun),
      disfun = !missing(disfun),
      centre_kernel = !missing(centre_kernel)
    ),
    list(...)
  )
  .validate_kema_solver(solver)

  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_range(sample_frac, c(0,1))
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_range(u, c(0,1))

  if (is.null(sigma)) {
    sigma <- 0.73
  }
  if (is.null(kernel)) {
    kernel <- kernlab::rbfdot(sigma = sigma)
  }

  result <- kema_orig.hyperdesign(
    data = data,
    y = !!y,
    preproc = preproc,
    ncomp = ncomp,
    knn = knn,
    sigma = sigma,
    mu = u,
    kernel = kernel,
    sample_frac = sample_frac,
    lambda = lambda,
    backend = backend,
    backend_control = backend_control,
    seed = seed
  )

  # Preserve legacy class tag for downstream code using `kema()`.
  result$classes <- "kema"
  if (is.numeric(result$eigenvalues)) {
    result$eigenvalues <- list(
      values = as.numeric(result$eigenvalues),
      all_values = as.numeric(result$eigenvalues),
      solver = "exact"
    )
  }
  result$regression_quality <- NULL
  result$retry_info <- NULL
  result

}


#' Choose optimal sigma for RBF kernel
#'
#' Estimates a reasonable sigma parameter for RBF kernels using the median
#' pairwise distance heuristic, which often works well in practice.
#'
#' @param X Data matrix (samples x features)
#' @param sample_size Maximum number of pairs to sample for distance computation
#'   (default: 1000 for efficiency)
#' @return Suggested sigma value
#' @examples
#' X <- matrix(rnorm(100), 20, 5)
#' sigma <- choose_sigma(X)
#' sigma
#' @export
choose_sigma <- function(X, sample_size = 1000) {
  if (!is.matrix(X) && !methods::is(X, "Matrix")) {
    stop("X must be a matrix", call. = FALSE)
  }

  n <- nrow(X)
  if (n < 2) {
    return(1.0)  # Default fallback
  }

  # For large datasets, sample pairs to avoid O(n^2) computation
  if (n > sqrt(sample_size)) {
    indices <- sample.int(n, min(n, floor(sqrt(sample_size))))
    X_sample <- X[indices, , drop = FALSE]
  } else {
    X_sample <- X
  }

  # Compute pairwise distances efficiently
  distances <- as.matrix(dist(as.matrix(X_sample)))

  # Remove diagonal (zero distances)
  distances <- distances[upper.tri(distances)]

  if (length(distances) == 0) {
    return(1.0)  # Fallback
  }

  # Median distance heuristic
  median_dist <- median(distances, na.rm = TRUE)

  # Convert to RBF sigma: typical rule of thumb is sigma = median_dist / sqrt(2)
  sigma <- median_dist / sqrt(2)

  # Ensure reasonable bounds
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- 1.0
  }

  sigma
}

#' @rdname kema
#' @method kema default
#' @export
kema.default <- function(data, ...) {
  stop("kema() requires either a hyperdesign or multidesign object. ",
       "Got: ", paste(class(data), collapse = ", "),
       ". See ?kema for usage examples.", call. = FALSE)
}
