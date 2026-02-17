#' Control settings for `gpca_align.hyperdesign()`
#'
#' @param knn Optional integer k for k-NN sparsification (applied to within and between graphs)
#' @param knn_mode Whether to use a "symmetric" (union) or "mutual" (intersection) k-NN graph
#' @param balance Whether to apply label balancing to "none", "within", or "all" edges
#' @param balance_power Exponent for inverse-frequency label scaling (>= 0)
#' @param normalize Normalisation mode for the component metrics: "fro", "edges", or "degree"
#' @param verbose Logical toggle for informational messages
#' @param max_dense_elems Maximum dense elements permitted when forming the block matrix (set `NA` to auto-tune)
#' @param memory_fraction Fraction of available RAM to consider safe when auto-tuning densification
#' @param memory_safety_bytes Absolute byte cushion to subtract when auto-tuning densification
#'
#' @return A list of class `gpca_align_control`
#' @examples
#' ctrl <- gpca_align_control(knn = 15, normalize = "fro")
#' @export
gpca_align_control <- function(
  knn = NA_integer_,
  knn_mode = c("symmetric", "mutual"),
  balance = c("none", "within", "all"),
  balance_power = 1,
  normalize = c("fro", "edges", "degree"),
  verbose = TRUE,
  max_dense_elems = 2e8,
  memory_fraction = 0.5,
  memory_safety_bytes = 512 * 1024^2
) {
  knn_mode <- match.arg(knn_mode)
  balance <- match.arg(balance)
  normalize <- match.arg(normalize)

  if (!(is.na(knn) || (is.numeric(knn) && length(knn) == 1L && is.finite(knn) && knn >= 1))) {
    stop("`knn` must be NA or a positive finite number.", call. = FALSE)
  }

  if (!is.numeric(balance_power) || length(balance_power) != 1L || is.na(balance_power) || balance_power < 0) {
    stop("`balance_power` must be a non-negative scalar.", call. = FALSE)
  }

  if (!(is.na(max_dense_elems) || (is.numeric(max_dense_elems) && length(max_dense_elems) == 1L && max_dense_elems > 0))) {
    stop("`max_dense_elems` must be NA or a positive number.", call. = FALSE)
  }

  if (!is.numeric(memory_fraction) || length(memory_fraction) != 1L || is.na(memory_fraction) || memory_fraction <= 0 || memory_fraction > 1) {
    stop("`memory_fraction` must be in (0, 1].", call. = FALSE)
  }

  if (!is.numeric(memory_safety_bytes) || length(memory_safety_bytes) != 1L || is.na(memory_safety_bytes) || memory_safety_bytes < 0) {
    stop("`memory_safety_bytes` must be a non-negative scalar.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      knn = if (is.na(knn)) NA_integer_ else as.integer(knn),
      knn_mode = knn_mode,
      balance = balance,
      balance_power = balance_power,
      normalize = normalize,
      verbose = verbose,
      max_dense_elems = if (is.na(max_dense_elems)) NA_real_ else as.numeric(max_dense_elems),
      memory_fraction = memory_fraction,
      memory_safety_bytes = memory_safety_bytes
    ),
    class = "gpca_align_control"
  )
}

resolve_gpca_align_control <- function(control) {
  defaults <- gpca_align_control()

  if (missing(control) || is.null(control)) {
    return(defaults)
  }

  if (inherits(control, "gpca_align_control")) {
    merged <- modifyList(defaults, control)
    return(do.call(gpca_align_control, merged))
  }

  if (!is.list(control)) {
    stop("`control` must be a list or a gpca_align_control() result.", call. = FALSE)
  }

  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  merged <- modifyList(defaults, control)
  do.call(gpca_align_control, merged)
}

#' Generalized PCA Alignment
#'
#' @param data A `hyperdesign` object or other supported input.
#' @param ... Additional arguments passed to methods.
#' @export
gpca_align <- function(data, ...) {
  UseMethod("gpca_align")
}

#' Generalized PCA Alignment for Hyperdesign Data
#'
#' Performs GPCA alignment on multi-domain data by constructing label-based similarity
#' metrics and extracting aligned components via generalized eigenanalysis.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the label variable to use for alignment
#' @param preproc Preprocessing function to apply to the data (default `center()`)
#' @param ncomp Number of components to extract (default 2)
#' @param simfun Function to compute similarity between labels
#' @param u Trade-off parameter between within-domain and between-domain similarity (0-1)
#' @param lambda Ridge regularisation added to the metric matrix
#' @param row_metric_scale Scaling applied to the assembled row metric before GPCA
#' @param control A list created by [gpca_align_control()] that tunes sparsification,
#'   balancing, normalisation, and memory safeguards
#' @param ... Additional arguments passed to the method (unused).
#'
#' @return A multiblock_biprojector object with aligned scores, loadings, and preprocessing info
#'
#' @examples
#' \donttest{
#' library(multidesign)
#' set.seed(1)
#' X1 <- matrix(rnorm(40), 20, 2)
#' X2 <- matrix(rnorm(40), 20, 2)
#' labs <- factor(rep(c("A", "B"), each = 10))
#' d1 <- list(x = X1, design = data.frame(label = labs))
#' d2 <- list(x = X2, design = data.frame(label = labs))
#' hd <- structure(list(d1, d2), class = "hyperdesign")
#' result <- gpca_align(hd, y = label, ncomp = 2)
#' }
#'
#' @importFrom genpca genpca
#' @importFrom multivarious prep init_transform pass multiblock_biprojector
#' @export
gpca_align.hyperdesign <- function(
  data, y,
  preproc = center(),
  ncomp = 2,
  simfun = neighborweights::binary_label_matrix,
  u = 0.5,
  lambda = 0.1,
  row_metric_scale = 1,
  control = gpca_align_control(),
  ...
) {
  ## -------------------------- control --------------------------
  control <- resolve_gpca_align_control(control)
  normalize_mode <- control$normalize
  verbose        <- isTRUE(control$verbose)

  knn_k    <- control$knn
  knn_mode <- control$knn_mode

  balance     <- control$balance
  balance_pow <- control$balance_power

  max_dense_elems_opt <- control$max_dense_elems
  mem_fraction        <- control$memory_fraction
  mem_safety          <- control$memory_safety_bytes

  ## -------------------------- argument checks --------------------------
  yn <- rlang::as_name(rlang::ensym(y))

  if (!is.list(data) || length(data) == 0L) {
    stop("`data` must be a non-empty list of hyperdesign blocks.")
  }

  if (!is.numeric(u) || length(u) != 1L || is.na(u) || u < 0 || u > 1) {
    stop("`u` must be a scalar in [0, 1].")
  }

  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("`lambda` must be a nonnegative scalar.")
  }

  if (!is.numeric(row_metric_scale) || length(row_metric_scale) != 1L ||
      !is.finite(row_metric_scale) || row_metric_scale <= 0) {
    stop("`row_metric_scale` must be a positive scalar.")
  }

  if (!is.function(simfun)) {
    stop("`simfun` must be a function(labels) -> similarity matrix.")
  }

  ## -------------------------- collect labels & verify blocks --------------------------
  nb <- length(data)
  label_list <- vector("list", nb)
  ns <- integer(nb)

  for (i in seq_len(nb)) {
    blk <- data[[i]]
    if (!is.list(blk) || is.null(blk$x) || is.null(blk$design)) {
      stop("Each block must be a list with $x and $design.")
    }

    Xi <- blk$x
    Di <- blk$design
    if (!is.matrix(Xi) && !inherits(Xi, "Matrix")) {
      stop("Block ", i, ": $x must be a matrix or Matrix.")
    }

    if (!is.data.frame(Di)) {
      stop("Block ", i, ": $design must be a data.frame.")
    }

    if (!yn %in% names(Di)) {
      stop("Block ", i, ": label column `", yn, "` not found in $design.")
    }

    yi <- as.factor(Di[[yn]])
    if (nrow(Xi) != length(yi)) {
      stop("Block ", i, ": nrow($x) != length($design[[y]]).")
    }

    label_list[[i]] <- yi
    ns[i] <- length(yi)
  }

  ## combine labels across domains with unified levels
  labels <- factor(unlist(lapply(label_list, as.character), use.names = FALSE))

  n <- length(labels)
  if (n == 0L) {
    stop("No samples found.")
  }

  ## -------------------------- similarity from simfun (sparse, symmetric, >=0) --------------------------
  W_full <- tryCatch(simfun(labels), error = function(e) {
    stop("`simfun` failed: ", e$message)
  })

  if (!inherits(W_full, "Matrix") && !is.matrix(W_full)) {
    stop("`simfun` must return a base matrix or a Matrix::Matrix object.")
  }

  if (!all(dim(W_full) == c(n, n))) {
    stop("`simfun` returned a matrix of size ", paste(dim(W_full), collapse = "x"),
         "; expected ", n, "x", n, ".")
  }

  if (!inherits(W_full, "sparseMatrix")) {
    if (n > 200 && verbose) {
      message("Converting `simfun` result to sparse; consider returning sparse directly.")
    }
    W_full <- Matrix::Matrix(W_full, sparse = TRUE)
  }

  W_full <- Matrix::forceSymmetric(W_full, uplo = "U")
  if (any(W_full@x < 0)) {
    if (verbose) {
      message("Clipping negative similarities to zero (required for PSD metric D+W).")
    }
    W_full <- methods::as(W_full, "dgCMatrix")
    W_full@x <- pmax(W_full@x, 0)
    W_full <- Matrix::forceSymmetric(W_full, uplo = "U")
  }
  if (any(Matrix::diag(W_full) != 0)) {
    W_full <- W_full - Matrix::Diagonal(n, Matrix::diag(W_full))
  }
  W_full <- Matrix::drop0(W_full)

  ## -------------------------- split within / between edges --------------------------
  block_ids <- rep.int(seq_len(nb), ns)
  W_full_dgc <- methods::as(W_full, "dgCMatrix")
  col_ptr <- W_full_dgc@p
  rows <- W_full_dgc@i + 1L
  cols <- rep.int(seq_len(n), diff(col_ptr))
  row_blocks <- block_ids[rows]
  col_blocks <- block_ids[cols]
  keep_within <- row_blocks == col_blocks

  if (length(keep_within) != length(W_full_dgc@x)) {
    stop("Internal error: mismatch in within-mask sizing.")
  }

  W_within <- Matrix::sparseMatrix(
    i = rows[keep_within],
    j = cols[keep_within],
    x = W_full_dgc@x[keep_within],
    dims = dim(W_full_dgc),
    symmetric = FALSE
  )
  W_within <- Matrix::forceSymmetric(Matrix::drop0(W_within), uplo = "U")

  W_between <- Matrix::drop0(W_full_dgc - W_within)
  W_between <- Matrix::forceSymmetric(W_between, uplo = "U")

  ## -------------------------- (NEW) label balancing --------------------------
  if (balance != "none") {
    counts <- table(labels)
    mean_n <- mean(as.numeric(counts))
    lab_scale <- (mean_n / as.numeric(counts))^balance_pow
    lab_scale <- setNames(lab_scale, names(counts))
    s_vec <- sqrt(unname(lab_scale[as.character(labels)]))
    Sdiag <- Matrix::Diagonal(n, s_vec)

    W_within <- Matrix::drop0(Sdiag %*% W_within %*% Sdiag)

    if (balance == "all") {
      W_between <- Matrix::drop0(Sdiag %*% W_between %*% Sdiag)
    }
    if (verbose) {
      message("Applied label balancing (", balance, "), power = ", balance_pow, ".")
    }
  }

  ## -------------------------- (NEW) k-NN sparsifier --------------------------
  knn_sparsify <- function(W, k, mode = c("symmetric", "mutual")) {
    mode <- match.arg(mode)
    if (!is.finite(k) || k <= 0) {
      return(W)
    }
    if (!inherits(W, "dgCMatrix")) {
      W <- methods::as(W, "dgCMatrix")
    }

    deg <- Matrix::rowSums(W > 0)
    if (all(deg <= k)) {
      return(Matrix::forceSymmetric(W, uplo = "U"))
    }

    summary_df <- Matrix::summary(W)
    if (nrow(summary_df) == 0L) {
      return(Matrix::forceSymmetric(Matrix::drop0(W), uplo = "U"))
    }

    row_entries <- split(summary_df, summary_df$i)
    nrows <- nrow(W)
    ii_list <- vector("list", nrows)
    jj_list <- vector("list", nrows)
    xx_list <- vector("list", nrows)

    for (nm in names(row_entries)) {
      i <- as.integer(nm)
      entries <- row_entries[[nm]]
      js <- entries$j
      xs <- entries$x
      if (length(js) > k) {
        sel <- utils::head(order(xs, decreasing = TRUE), k)
        js <- js[sel]
        xs <- xs[sel]
      }
      ii_list[[i]] <- rep.int(i, length(js))
      jj_list[[i]] <- js
      xx_list[[i]] <- xs
    }

    ii <- unlist(ii_list, use.names = FALSE)
    jj <- unlist(jj_list, use.names = FALSE)
    xx <- unlist(xx_list, use.names = FALSE)

    if (!length(ii)) {
      return(Matrix::forceSymmetric(Matrix::drop0(W), uplo = "U"))
    }

    Wrow <- Matrix::sparseMatrix(i = ii, j = jj, x = xx, dims = dim(W))

    if (mode == "symmetric") {
      Ws <- pmax(Wrow, Matrix::t(Wrow))
    } else {
      Ws <- pmin(Wrow, Matrix::t(Wrow))
    }
    Matrix::forceSymmetric(Matrix::drop0(Ws), uplo = "U")
  }

  if (is.finite(knn_k) && knn_k >= 1) {
    if (verbose) {
      message("Applying k-NN sparsifier with k = ", knn_k, " (", knn_mode, ").")
    }
    W_within  <- knn_sparsify(W_within,  knn_k, knn_mode)
    W_between <- knn_sparsify(W_between, knn_k, knn_mode)
  }

  ## -------------------------- build PSD metrics --------------------------
  build_metric_psd <- function(W) {
    W <- Matrix::forceSymmetric(W, uplo = "U")
    if (any(Matrix::diag(W) != 0)) {
      W <- W - Matrix::Diagonal(n, Matrix::diag(W))
    }
    W <- Matrix::drop0(W)
    D <- Matrix::Diagonal(n, as.numeric(Matrix::rowSums(W)))
    Matrix::forceSymmetric(D + W, uplo = "U")
  }

  M_within  <- build_metric_psd(W_within)
  M_between <- build_metric_psd(W_between)

  ## -------------------------- normalize component metrics --------------------------
  norm_fro    <- function(A) Matrix::norm(A, "F")
  norm_edges  <- function(W) if (length(W@x)) sum(W@x) else 0
  norm_degree <- function(W) sum(Matrix::rowSums(W))

  scale_with <- switch(
    normalize_mode,
    "edges"  = max(norm_edges(W_within),  .Machine$double.eps),
    "degree" = max(norm_degree(W_within), .Machine$double.eps),
    max(norm_fro(M_within), .Machine$double.eps)
  )
  scale_between <- switch(
    normalize_mode,
    "edges"  = max(norm_edges(W_between),  .Machine$double.eps),
    "degree" = max(norm_degree(W_between), .Machine$double.eps),
    max(norm_fro(M_between), .Machine$double.eps)
  )

  M_within_n  <- M_within  / scale_with
  M_between_n <- M_between / scale_between

  M <- u * M_within_n + (1 - u) * M_between_n
  M <- Matrix::forceSymmetric(M, uplo = "U")
  if (lambda > 0) {
    M <- M + Matrix::Diagonal(n, lambda)
  }

  eig_la <- NA_real_
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    eig_la <- tryCatch(
      RSpectra::eigs_sym(M, 1L, which = "LA", opts = list(retvec = FALSE))$values[1],
      error = function(e) NA_real_
    )
  } else if (requireNamespace("PRIMME", quietly = TRUE)) {
    eig_la <- tryCatch(
      PRIMME::eigs_sym(M, NEig = 1L, which = "LA", method = "PRIMME_DEFAULT_MIN_MATVECS")$values[1],
      error = function(e) NA_real_
    )
  } else if (n <= 1200L) {
    eig_la <- tryCatch(
      max(eigen(as.matrix(M), symmetric = TRUE, only.values = TRUE)$values),
      error = function(e) NA_real_
    )
  }
  if (is.finite(eig_la) && eig_la > 0) {
    M <- M / eig_la
  }

  ## -------------------------- preprocessing per domain --------------------------
  pdata <- data
  proclist <- vector("list", nb)

  broadcast <- is.function(preproc) || inherits(preproc, "prepper")
  if (broadcast) {
    if (verbose) {
      message("Applying the same preprocessor to all ", nb, " domains.")
    }
    preproc <- rep(list(preproc), nb)
  } else if (is.list(preproc) && length(preproc) != nb) {
    stop("Length of `preproc` list (", length(preproc), ") must match number of domains (", nb, ").")
  }

  for (i in seq_len(nb)) {
    if (!is.null(preproc) && !is.null(preproc[[i]])) {
      templ <- multivarious::prep(preproc[[i]])
      Xi_tf <- multivarious::init_transform(templ, data[[i]]$x)
      pdata[[i]]$x <- Xi_tf
      proclist[[i]] <- attr(Xi_tf, "preproc") %||% templ
    } else {
      pdata[[i]]$x <- data[[i]]$x
      proclist[[i]] <- NULL
    }
  }
  names(proclist) <- names(data)

  sample_block_idx <- block_indices(pdata)
  sample_blocks_list <- split(sample_block_idx, row(sample_block_idx))
  names(sample_blocks_list) <- names(pdata)
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)

  feat_per_block <- vapply(pdata, function(b) ncol(b$x), integer(1))
  end_idx   <- cumsum(feat_per_block)
  start_idx <- c(1L, head(end_idx, -1L) + 1L)
  feature_block_idx <- lapply(seq_along(feat_per_block), function(i) start_idx[i]:end_idx[i])
  names(feature_block_idx) <- names(pdata)

  ## -------------------------- assemble X (block-diagonal) --------------------------
  X_block <- Matrix::bdiag(lapply(pdata, function(b) {
    if (inherits(b$x, "Matrix")) {
      b$x
    } else {
      Matrix::Matrix(b$x, sparse = TRUE)
    }
  }))

  max_dense_elems <- max_dense_elems_opt
  if (is.na(max_dense_elems_opt)) {
    avail_bytes <- NA_real_
    if (requireNamespace("ps", quietly = TRUE)) {
      sm <- tryCatch(ps::ps_system_memory(), error = function(e) NULL)
      if (!is.null(sm)) {
        avail_bytes <- if (!is.null(sm$available)) {
          sm$available
        } else if (!is.null(sm$free)) {
          sm$free
        } else {
          NA_real_
        }
      }
    }
    if (is.finite(avail_bytes)) {
      target <- max((avail_bytes - mem_safety) * mem_fraction, 100e6)
      max_dense_elems <- floor(target / 8)  # doubles
      if (verbose) {
        message(
          "Adaptive densify cap: ~",
          format(target / 1024^2, digits = 3),
          " MB => ",
          format(max_dense_elems, scientific = TRUE),
          " elements."
        )
      }
    } else {
      max_dense_elems <- 2e8
      if (verbose) {
        message(
          "Adaptive memory unavailable; using fallback cap of ",
          format(max_dense_elems, scientific = TRUE),
          " elements."
        )
      }
    }
  }

  dense_elems <- as.double(nrow(X_block)) * as.double(ncol(X_block))
  if (dense_elems > max_dense_elems) {
    stop(
      "Densifying X would allocate ~",
      format(dense_elems, scientific = TRUE),
      " elements (> cap). Reduce dimensionality, use stronger preprocessing, ",
      "increase k-NN sparsity on M (which may enable stronger preprocs), ",
      "or call gpca_align_control(max_dense_elems = ...) and pass it via the `control` argument."
    )
  }
  X_block_dense <- as.matrix(X_block)

  ## -------------------------- run GPCA --------------------------
  ret <- genpca::genpca(
    X_block_dense,
    M = row_metric_scale * M,
    ncomp = ncomp,
    preproc = multivarious::pass()
  )

  multivarious::multiblock_biprojector(
    v = ret$v,
    s = ret$s,
    sdev = ret$sdev,
    preproc = proc,
    block_indices = feature_block_idx,
    labels = labels,
    classes = "gpca_align"
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a
