#' Structural graph descriptors (HKS, diffusion coordinates)
#'
#' @description
#' These helpers compute **intrinsic** node descriptors from a within-domain
#' spectral basis (eigenvalues/eigenvectors). They are designed to be reusable
#' across multiple alignment methods (e.g., GRASP, spectral MNN, and as features
#' for UOT costs).
#'
#' The main entry points are:
#' - [hks_time_grid()] to construct a diffusion-time grid.
#' - [compute_hks_descriptors()] to compute Heat Kernel Signatures (HKS).
#' - [compute_diffusion_coordinates()] to compute diffusion coordinates at one
#'   or more times.
#'
#' @details
#' **Basis format**
#'
#' These functions accept either:
#' - a list with fields `values` (length K) and `vectors` (N x K), or
#' - explicit `values`/`vectors` arguments.
#'
#' For Laplacian-based bases, the first eigenvalue is typically 0 and
#' corresponds to a constant mode. By default, eigenpairs with
#' `value <= eigen_tol` are dropped.
#'
#' @return
#' See individual functions.
#'
#' @name structure_descriptors
NULL

#' Construct a diffusion-time grid for HKS/diffusion descriptors
#'
#' @param values Numeric vector of eigenvalues (typically Laplacian eigenvalues).
#' @param q Integer number of time points.
#' @param mode Either `"auto"` (derive a time range from the spectrum) or
#'   `"fixed"` (use `time_range`).
#' @param time_range Length-2 positive numeric vector giving min/max times when
#'   `mode="fixed"`.
#' @param spacing Either `"log"` (default; typical for HKS) or `"linear"`.
#' @param time_scale Positive scalar multiplier applied to the resulting times.
#' @param eigen_tol Positive threshold; eigenvalues `<= eigen_tol` are ignored
#'   when `mode="auto"`.
#'
#' @return Numeric vector of length `q`.
#' @export
hks_time_grid <- function(values,
                          q,
                          mode = c("auto", "fixed"),
                          time_range = c(0.1, 50),
                          spacing = c("log", "linear"),
                          time_scale = 1,
                          eigen_tol = 1e-12) {

  mode <- match.arg(mode)
  spacing <- match.arg(spacing)

  chk::chk_numeric(values)
  chk::chk_number(q)
  chk::chk_true(q >= 1)
  q <- as.integer(q)

  chk::chk_numeric(time_range)
  chk::chk_true(length(time_range) == 2L)
  chk::chk_true(all(is.finite(time_range)))
  chk::chk_true(all(time_range > 0))
  chk::chk_true(time_range[2] > time_range[1])

  chk::chk_number(time_scale)
  chk::chk_true(is.finite(time_scale))
  chk::chk_true(time_scale > 0)

  chk::chk_number(eigen_tol)
  chk::chk_true(is.finite(eigen_tol))
  chk::chk_true(eigen_tol > 0)

  vals <- as.numeric(values)
  vals <- vals[is.finite(vals) & vals > eigen_tol]

  if (mode == "auto") {
    if (!length(vals)) {
      mode <- "fixed"
    }
  }

  if (mode == "fixed") {
    if (spacing == "linear") {
      times <- seq(time_range[1], time_range[2], length.out = q)
    } else {
      times <- exp(seq(log(time_range[1]), log(time_range[2]), length.out = q))
    }
    return(times * time_scale)
  }

  lam_min <- min(vals)
  lam_max <- max(vals)

  # Standard HKS heuristic: times ~ [1/λ_max, 1/λ_min], clamped.
  t_min <- 1 / (lam_max + 1e-12)
  t_max <- 1 / (lam_min + 1e-12)
  t_min <- max(t_min, 1e-3)
  t_max <- min(t_max, 1e3)
  if (!(t_max > t_min)) {
    t_min <- min(t_min, t_max)
    t_max <- max(t_min * 10, t_max)
  }

  if (spacing == "linear") {
    times <- seq(t_min, t_max, length.out = q)
  } else {
    times <- exp(seq(log(t_min), log(t_max), length.out = q))
  }
  times * time_scale
}

.rspectra_v0_deterministic <- function(n, seed = 1L) {
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

#' Compute a (normalized) graph Laplacian from an adjacency matrix
#'
#' @param A Square adjacency/affinity matrix (dense or sparse). Values should be
#'   nonnegative; the diagonal is set to 0 by default.
#' @param normalized Logical; if TRUE (default), compute the normalized Laplacian
#'   \eqn{L = I - D^{-1/2} A D^{-1/2}}. If FALSE, compute the unnormalized
#'   Laplacian \eqn{L = D - A}.
#' @param symmetrize Logical; if TRUE (default), symmetrize via `(A + t(A))/2`.
#' @param zero_diagonal Logical; if TRUE (default), set `diag(A)=0` before
#'   forming the Laplacian.
#'
#' @return A sparse symmetric matrix (CsparseMatrix).
#' @export
compute_graph_laplacian <- function(A,
                                    normalized = TRUE,
                                    symmetrize = TRUE,
                                    zero_diagonal = TRUE) {

  chk::chk_logical(normalized)
  chk::chk_logical(symmetrize)
  chk::chk_logical(zero_diagonal)

  if (!inherits(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  } else {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }

  if (nrow(A) != ncol(A)) stop("`A` must be square.", call. = FALSE)

  if (isTRUE(symmetrize)) {
    A <- (A + Matrix::t(A)) / 2
  }
  if (isTRUE(zero_diagonal)) {
    Matrix::diag(A) <- 0
    A <- Matrix::drop0(A)
  }

  deg <- Matrix::rowSums(A)
  n <- nrow(A)

  if (isTRUE(normalized)) {
    deg_safe <- pmax(deg, 1)
    D_inv_sqrt <- Matrix::Diagonal(x = 1 / sqrt(deg_safe))
    iso <- deg == 0
    if (any(iso)) {
      idx <- which(iso)
      D_inv_sqrt[idx, idx] <- 0
    }
    L <- Matrix::Diagonal(n) - D_inv_sqrt %*% A %*% D_inv_sqrt
  } else {
    L <- Matrix::Diagonal(x = deg) - A
  }

  # Return a general sparse matrix (not a symmetric class) for compatibility
  # with downstream eigensolvers such as RSpectra.
  methods::as(Matrix::Matrix(L, sparse = TRUE), "dgCMatrix")
}

#' Compute a Laplacian spectral basis (eigenpairs) from an adjacency matrix
#'
#' @param A Square adjacency/affinity matrix (dense or sparse).
#' @param k Integer number of non-trivial eigenpairs to return.
#' @param normalized Logical; use the normalized Laplacian (default TRUE).
#' @param which Which eigenpairs to request (passed to the eigensolver). Default
#'   is `"SA"` (smallest algebraic) for Laplacian bases.
#' @param eig_tol Positive threshold; eigenvalues `<= eig_tol` are treated as
#'   trivial and dropped.
#' @param extra Integer number of extra eigenpairs requested to survive multiple
#'   near-zero modes on disconnected graphs.
#' @param seed Integer seed used to build a deterministic initial vector for
#'   RSpectra (helps reproducibility).
#' @param tol Numeric tolerance passed to RSpectra/PRIMME when available.
#' @param maxitr Maximum iterations passed to RSpectra when available.
#'
#' @return A list with fields `values` (length k) and `vectors` (n x k).
#' @export
compute_laplacian_basis <- function(A,
                                   k,
                                   normalized = TRUE,
                                   which = "SA",
                                   eig_tol = 1e-12,
                                   extra = 4L,
                                   seed = 1L,
                                   tol = 1e-6,
                                   maxitr = 1000) {

  chk::chk_number(k)
  chk::chk_true(k >= 1)
  k <- as.integer(k)

  chk::chk_logical(normalized)
  chk::chk_character(which)
  chk::chk_number(eig_tol)
  chk::chk_true(is.finite(eig_tol))
  chk::chk_true(eig_tol > 0)
  chk::chk_number(extra)
  chk::chk_true(extra >= 0)
  extra <- as.integer(extra)

  L <- compute_graph_laplacian(A, normalized = normalized, symmetrize = TRUE, zero_diagonal = TRUE)
  n <- nrow(L)
  if (n < 3L) stop("Graph must have at least 3 nodes.", call. = FALSE)

  k_request <- min(k + extra, n - 1L)
  if (k_request < 2L) stop("Not enough eigenpairs available.", call. = FALSE)

  decomp <- if (requireNamespace("RSpectra", quietly = TRUE)) {
    v0 <- tryCatch(.rspectra_v0_deterministic(n, seed = seed), error = function(e) NULL)
    opts <- list(tol = tol, maxitr = as.integer(maxitr))
    if (!is.null(v0)) opts$initvec <- v0
    safe_compute(
      RSpectra::eigs_sym(L, k = k_request, which = which, opts = opts),
      "Eigen-decomposition failed in compute_laplacian_basis()"
    )
  } else if (requireNamespace("PRIMME", quietly = TRUE)) {
    safe_compute(
      PRIMME::eigs_sym(L, NEig = k_request, which = which),
      "Eigen-decomposition failed in compute_laplacian_basis()"
    )
  } else {
    eig_full <- eigen(as.matrix(L), symmetric = TRUE)
    list(values = eig_full$values[seq_len(k_request)],
         vectors = eig_full$vectors[, seq_len(k_request), drop = FALSE])
  }

  evals <- as.numeric(decomp$values)
  evecs <- as.matrix(decomp$vectors)

  ord <- if (which %in% c("LA", "LM")) order(evals, decreasing = TRUE) else order(evals)
  evals <- evals[ord]
  evecs <- evecs[, ord, drop = FALSE]

  keep <- which(is.finite(evals) & evals > eig_tol)
  if (!length(keep)) stop("No non-trivial eigenvalues found.", call. = FALSE)
  keep <- keep[seq_len(min(k, length(keep)))]

  list(values = evals[keep], vectors = evecs[, keep, drop = FALSE])
}

#' Convenience: compute HKS descriptors from an adjacency matrix
#'
#' @param A Square adjacency/affinity matrix (dense or sparse).
#' @param k_embed Number of eigenpairs used for the descriptor basis.
#' @param q Number of HKS time steps.
#' @param use_normalized_laplacian Logical; use the normalized Laplacian
#'   (default TRUE).
#' @param store_basis Logical; if TRUE, return `list(descriptors=, basis=)`.
#' @param ... Passed to [compute_hks_descriptors()] (e.g., time grid settings).
#'
#' @return A numeric matrix (n x q), or a list when `store_basis=TRUE`.
#' @export
compute_hks_from_adjacency <- function(A,
                                      k_embed = 30L,
                                      q = 16L,
                                      use_normalized_laplacian = TRUE,
                                      store_basis = FALSE,
                                      ...) {

  basis <- compute_laplacian_basis(A, k = k_embed, normalized = use_normalized_laplacian)
  desc <- compute_hks_descriptors(basis = basis, q = q, ...)
  if (isTRUE(store_basis)) {
    return(list(descriptors = desc, basis = basis))
  }
  desc
}

#' Convenience: compute diffusion coordinates from an adjacency matrix
#'
#' @param A Square adjacency/affinity matrix (dense or sparse).
#' @param k_embed Number of eigenpairs used for the basis.
#' @param time Positive scalar or vector of times.
#' @param use_normalized_laplacian Logical; use the normalized Laplacian
#'   (default TRUE).
#' @param store_basis Logical; if TRUE, return `list(coords=, basis=)`.
#' @param ... Passed to [compute_diffusion_coordinates()].
#'
#' @return A numeric matrix (n x k) or list of matrices, or a list when
#'   `store_basis=TRUE`.
#' @export
compute_diffusion_coordinates_from_adjacency <- function(A,
                                                        k_embed = 30L,
                                                        time = 1,
                                                        use_normalized_laplacian = TRUE,
                                                        store_basis = FALSE,
                                                        ...) {
  basis <- compute_laplacian_basis(A, k = k_embed, normalized = use_normalized_laplacian)
  coords <- compute_diffusion_coordinates(basis = basis, time = time, ...)
  if (isTRUE(store_basis)) {
    return(list(coords = coords, basis = basis))
  }
  coords
}

.as_basis_list <- function(basis = NULL, values = NULL, vectors = NULL) {
  if (!is.null(basis)) {
    if (!is.list(basis) || is.null(basis$values) || is.null(basis$vectors)) {
      stop("`basis` must be a list with fields `values` and `vectors`.", call. = FALSE)
    }
    values <- basis$values
    vectors <- basis$vectors
  }
  if (is.null(values) || is.null(vectors)) {
    stop("Provide either `basis` or both `values` and `vectors`.", call. = FALSE)
  }
  list(values = as.numeric(values), vectors = vectors)
}

.normalize_descriptor_matrix <- function(desc, normalize = c("both", "column", "row", "none")) {
  normalize <- match.arg(normalize)
  desc <- as.matrix(desc)
  if (!nrow(desc) || !ncol(desc)) return(desc)

  if (normalize %in% c("column", "both")) {
    cn <- sqrt(colSums(desc * desc))
    cn[!is.finite(cn) | cn < 1e-12] <- 1
    desc <- sweep(desc, 2, cn, "/")
  }
  if (normalize %in% c("row", "both")) {
    rn <- sqrt(rowSums(desc * desc))
    rn[!is.finite(rn) | rn < 1e-12] <- 1
    desc <- sweep(desc, 1, rn, "/")
  }
  desc
}

.square_eigenvectors <- function(phi) {
  if (inherits(phi, "Matrix")) {
    phi_sq <- Matrix::Matrix(phi, sparse = TRUE)
    phi_sq@x <- phi_sq@x^2
    return(phi_sq)
  }
  phi <- as.matrix(phi)
  phi * phi
}

#' Compute Heat Kernel Signatures (HKS) from a spectral basis
#'
#' @param basis Optional basis list with fields `values` and `vectors`.
#' @param values Optional numeric vector of eigenvalues (overrides `basis$values`).
#' @param vectors Optional matrix (or sparse Matrix) of eigenvectors (overrides
#'   `basis$vectors`).
#' @param q Integer number of time points.
#' @param times Optional numeric time grid. If provided, `time_*` arguments are
#'   ignored.
#' @param time_mode Either `"auto"` (derive a time range from the spectrum) or
#'   `"fixed"` (use `time_range`).
#' @param time_range Length-2 positive numeric vector giving min/max times when
#'   `time_mode="fixed"`.
#' @param spacing Either `"log"` (default; typical for HKS) or `"linear"`.
#' @param time_scale Positive scalar multiplier applied to the time grid.
#' @param eigen_tol Positive threshold; eigenpairs with `value <= eigen_tol` are
#'   dropped.
#' @param normalize Descriptor normalization: `"both"` (default), `"column"`,
#'   `"row"`, or `"none"`.
#' @param clamp_exponent Lower clamp for the exponent (default `-745`), used to
#'   avoid denormal underflow in `exp(-λ t)`.
#'
#' @return A numeric matrix of dimension N x q.
#' @export
compute_hks_descriptors <- function(basis = NULL,
                                    values = NULL,
                                    vectors = NULL,
                                    q = 16L,
                                    times = NULL,
                                    time_mode = c("auto", "fixed"),
                                    time_range = c(0.1, 50),
                                    spacing = c("log", "linear"),
                                    time_scale = 1,
                                    eigen_tol = 1e-12,
                                    normalize = c("both", "column", "row", "none"),
                                    clamp_exponent = -745) {

  time_mode <- match.arg(time_mode)
  spacing <- match.arg(spacing)
  normalize <- match.arg(normalize)

  chk::chk_number(q)
  chk::chk_true(q >= 1)
  q <- as.integer(q)

  chk::chk_number(eigen_tol)
  chk::chk_true(is.finite(eigen_tol))
  chk::chk_true(eigen_tol > 0)

  chk::chk_number(clamp_exponent)
  chk::chk_true(is.finite(clamp_exponent))
  chk::chk_true(clamp_exponent < 0)

  b <- .as_basis_list(basis = basis, values = values, vectors = vectors)
  lambda_vals <- as.numeric(b$values)
  phi <- b$vectors

  if (is.null(dim(phi)) || nrow(phi) < 1L) {
    stop("`vectors` must be a matrix-like object with at least one row.", call. = FALSE)
  }

  keep <- which(is.finite(lambda_vals) & lambda_vals > eigen_tol)
  if (!length(keep)) {
    stop("No eigenvalues > eigen_tol; check graph connectivity or reduce eigen_tol.", call. = FALSE)
  }

  lambda_vals <- lambda_vals[keep]
  phi <- if (inherits(phi, "Matrix")) phi[, keep, drop = FALSE] else as.matrix(phi)[, keep, drop = FALSE]

  tgrid <- if (!is.null(times)) {
    times <- as.numeric(times)
    if (!length(times) || any(!is.finite(times)) || any(times <= 0)) {
      stop("`times` must be a positive finite numeric vector.", call. = FALSE)
    }
    times
  } else {
    hks_time_grid(lambda_vals, q = q, mode = time_mode, time_range = time_range,
                  spacing = spacing, time_scale = time_scale, eigen_tol = eigen_tol)
  }

  phi_sq <- .square_eigenvectors(phi)
  logH <- -outer(lambda_vals, tgrid)
  H <- exp(pmax(logH, clamp_exponent))
  desc <- phi_sq %*% H
  .normalize_descriptor_matrix(desc, normalize = normalize)
}

#' Compute diffusion coordinates from a spectral basis
#'
#' @description
#' Diffusion coordinates at time \eqn{t} are given by
#' \eqn{\Phi \mathrm{diag}(\exp(-\lambda t))}, where \eqn{(\Phi,\lambda)} are
#' Laplacian eigenvectors/eigenvalues.
#'
#' @param basis Optional basis list with fields `values` and `vectors`.
#' @param values Optional numeric vector of eigenvalues (overrides `basis$values`).
#' @param vectors Optional matrix (or sparse Matrix) of eigenvectors (overrides
#'   `basis$vectors`).
#' @param time Positive numeric scalar or vector of times.
#' @param k Optional integer number of coordinates to keep (defaults to all).
#' @param eigen_tol Positive threshold; eigenpairs with `value <= eigen_tol` are
#'   dropped.
#' @param normalize Optional normalization: `"none"` (default) or `"row"`.
#'
#' @return If `time` is length-1, a numeric matrix N x k. If `time` has length
#'   > 1, a list of matrices (one per time).
#' @export
compute_diffusion_coordinates <- function(basis = NULL,
                                         values = NULL,
                                         vectors = NULL,
                                         time = 1,
                                         k = NULL,
                                         eigen_tol = 1e-12,
                                         normalize = c("none", "row")) {

  normalize <- match.arg(normalize)
  chk::chk_number(eigen_tol)
  chk::chk_true(is.finite(eigen_tol))
  chk::chk_true(eigen_tol > 0)

  b <- .as_basis_list(basis = basis, values = values, vectors = vectors)
  lambda_vals <- as.numeric(b$values)
  phi <- b$vectors

  keep <- which(is.finite(lambda_vals) & lambda_vals > eigen_tol)
  if (!length(keep)) {
    stop("No eigenvalues > eigen_tol; check graph connectivity or reduce eigen_tol.", call. = FALSE)
  }

  lambda_vals <- lambda_vals[keep]
  phi <- if (inherits(phi, "Matrix")) phi[, keep, drop = FALSE] else as.matrix(phi)[, keep, drop = FALSE]

  if (!is.null(k)) {
    chk::chk_number(k)
    chk::chk_true(k >= 1)
    k <- min(as.integer(k), length(lambda_vals))
    lambda_vals <- lambda_vals[seq_len(k)]
    phi <- phi[, seq_len(k), drop = FALSE]
  }

  time <- as.numeric(time)
  if (!length(time) || any(!is.finite(time)) || any(time <= 0)) {
    stop("`time` must be a positive finite numeric scalar or vector.", call. = FALSE)
  }

  compute_one <- function(t) {
    w <- exp(-lambda_vals * t)
    coords <- if (inherits(phi, "Matrix")) phi %*% Matrix::Diagonal(x = w) else sweep(phi, 2, w, "*")
    .normalize_descriptor_matrix(coords, normalize = if (normalize == "row") "row" else "none")
  }

  if (length(time) == 1L) {
    return(compute_one(time))
  }
  lapply(time, compute_one)
}
