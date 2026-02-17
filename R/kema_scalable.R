#' Low-rank class graph helpers and matrix-free kernel operators for scalable KEMA
#'
#' These helpers avoid materializing dense class similarity matrices and block
#' diagonal kernel matrices when solving the KEMA generalized eigenproblem. They
#' expose compact representations that can be combined to build reduced Grams in
#' O(n + rC) memory, where r is the aggregate kernel rank and C the number of
#' classes.
#'
#' @param labels Character/factor vector with NA for unlabeled samples
#' @return List with sparse label factor matrix `C`, indicator vector `ell`, and
#'   diagonal entries for the same/different-class Laplacians.
#' @keywords internal
label_factors <- function(labels) {
  n <- length(labels)
  lab_mask <- !is.na(labels)

  if (!any(lab_mask)) {
    stop("label_factors: need at least one labeled sample", call. = FALSE)
  }

  labs <- factor(labels[lab_mask])
  C <- Matrix::sparseMatrix(
    i = which(lab_mask),
    j = as.integer(labs),
    x = 1,
    dims = c(n, nlevels(labs)),
    giveCsparse = TRUE
  )
  n_lab <- sum(lab_mask)
  n_per_c <- Matrix::colSums(C)
  d_s <- numeric(n)
  d_s[lab_mask] <- n_per_c[as.integer(labs)]

  d_d <- numeric(n)
  if (n_lab > 0) {
    d_d_vals <- n_lab - n_per_c[as.integer(labs)]
    d_d[lab_mask] <- d_d_vals
  }

  list(C = C, ell = as.numeric(lab_mask), d_s = d_s, d_d = d_d)
}

#' Build matrix-free operator for block-diagonal kernel matrix constructed from
#' pre-computed kernel blocks.
#'
#' @param Ks List of `Matrix` objects (each n_i x r_i) representing kernel
#'   evaluations for a block.
#' @return List with dimensions and forward/transpose operators.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' K2 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1, K2))
#' str(Zop)
#' }
#' @export
make_Zop_from_Ks <- function(Ks) {
  if (!length(Ks)) {
    stop("make_Zop_from_Ks: empty kernel list", call. = FALSE)
  }

  block_rows <- vapply(Ks, nrow, integer(1))
  block_cols <- vapply(Ks, ncol, integer(1))
  n <- sum(block_rows)
  r <- sum(block_cols)
  col_offsets <- c(0L, cumsum(block_cols))
  row_offsets <- c(0L, cumsum(block_rows))

  mv <- function(x) {
    X <- if (is.matrix(x)) x else matrix(x, nrow = r)
    out <- matrix(0, nrow = n, ncol = ncol(X))
    for (i in seq_along(Ks)) {
      idx_col <- (col_offsets[i] + 1L):col_offsets[i + 1L]
      idx_row <- (row_offsets[i] + 1L):row_offsets[i + 1L]
      block_mat <- as.matrix(Ks[[i]] %*% X[idx_col, , drop = FALSE])
      out[idx_row, ] <- block_mat
    }
    if (ncol(out) == 1L) drop(out) else out
  }

  tmv <- function(y) {
    Y <- if (is.matrix(y)) y else matrix(y, nrow = n)
    out <- matrix(0, nrow = r, ncol = ncol(Y))
    for (i in seq_along(Ks)) {
      idx_col <- (col_offsets[i] + 1L):col_offsets[i + 1L]
      idx_row <- (row_offsets[i] + 1L):row_offsets[i + 1L]
      block_mat <- Matrix::crossprod(Ks[[i]], Y[idx_row, , drop = FALSE])
      out[idx_col, ] <- as.matrix(block_mat)
    }
    if (ncol(out) == 1L) drop(out) else out
  }

  list(n = n, r = r, mv = mv, tmv = tmv, block_cols = block_cols)
}

#' Compute Z^T Z without forming Z explicitly.
#'
#' @param Zop Operator as returned by `make_Zop_from_Ks` or similar.
#' @param block_size Number of identity columns processed per batch.
#' @return Symmetric dense matrix of size r x r.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' K2 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1, K2))
#' G <- gram_ZZ(Zop)
#' }
#' @export
gram_ZZ <- function(Zop, block_size = 128L) {
  r <- Zop$r
  G <- matrix(0, r, r)
  eye <- diag(1, r)
  for (start in seq.int(1L, r, by = block_size)) {
    end <- min(start + block_size - 1L, r)
    cols <- eye[, start:end, drop = FALSE]
    Zcols <- Zop$mv(cols)
    G[, start:end] <- Zop$tmv(Zcols)
  }
  G <- 0.5 * (G + t(G))
  dense <- Matrix::Matrix(G, sparse = FALSE)
  methods::as(dense, "generalMatrix")
}

#' Generic helper for Z^T A Z using an operator.
#'
#' @param Zop Matrix-free operator.
#' @param A Sparse or diagonal matrix compatible with Z columns.
#' @param block_size Number of columns processed per batch.
#' @return Symmetric dense matrix.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1))
#' A <- Diagonal(10, 1)
#' G <- gram_Z_A_Z(Zop, A)
#' }
#' @export
gram_Z_A_Z <- function(Zop, A, block_size = 64L) {
  if (inherits(A, "diagonalMatrix")) {
    d <- Matrix::diag(A)
    return(gram_Z_diag_Z(Zop, d, block_size))
  }
  if (!methods::is(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  r <- Zop$r
  G <- matrix(0, r, r)
  eye <- diag(1, r)
  for (start in seq.int(1L, r, by = block_size)) {
    end <- min(start + block_size - 1L, r)
    cols <- eye[, start:end, drop = FALSE]
    Zcols <- Zop$mv(cols)
    AZ <- A %*% Zcols
    G[, start:end] <- Zop$tmv(AZ)
  }
  G <- 0.5 * (G + t(G))
  dense <- Matrix::Matrix(G, sparse = FALSE)
  methods::as(dense, "generalMatrix")
}

#' Fast path for diagonal A in Z^T A Z.
#'
#' @param Zop Operator for Z.
#' @param diag_entries Numeric vector of length n (diagonal of A).
#' @param block_size Batch size for columns.
#' @return Symmetric dense matrix.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1))
#' diag_vals <- rep(1, 10)
#' G <- gram_Z_diag_Z(Zop, diag_vals)
#' }
#' @export
gram_Z_diag_Z <- function(Zop, diag_entries, block_size = 128L) {
  r <- Zop$r
  n <- length(diag_entries)
  if (n != Zop$n) {
    stop("gram_Z_diag_Z: diagonal length does not match operator rows", call. = FALSE)
  }
  G <- matrix(0, r, r)
  eye <- diag(1, r)
  for (start in seq.int(1L, r, by = block_size)) {
    end <- min(start + block_size - 1L, r)
    cols <- eye[, start:end, drop = FALSE]
    Zcols <- Zop$mv(cols)
    Dcols <- Zcols * diag_entries
    G[, start:end] <- Zop$tmv(Dcols)
  }
  G <- 0.5 * (G + t(G))
  dense <- Matrix::Matrix(G, sparse = FALSE)
  methods::as(dense, "generalMatrix")
}

#' Apply Z^T to arbitrary multi-vectors without forming Z.
#'
#' @param Zop Operator for Z.
#' @param S Dense or sparse matrix with n rows.
#' @return Dense matrix with r rows.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1))
#' S <- matrix(rnorm(30), 10, 3)
#' result <- Zt_times(Zop, S)
#' }
#' @export
Zt_times <- function(Zop, S) {
  if (!is.matrix(S) && !methods::is(S, "Matrix")) {
    S <- matrix(S, nrow = Zop$n)
  }
  Zop$tmv(S)
}

#' Build Gram for same-class Laplacian using label factors.
#'
#' Z^T L_s Z = Z^T D_s Z - (Z^T C)(Z^T C)^T
#'
#' @param Zop Operator for kernel blocks.
#' @param F Label factor list from `label_factors()`.
#' @return Symmetric dense matrix of size r x r.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1))
#' # F <- label_factors(...)  # Requires label factors
#' # G <- gram_Ls(Zop, F)
#' }
gram_Ls <- function(Zop, F) {
  Ds <- Matrix::Diagonal(x = F$d_s)
  XtDX <- gram_Z_A_Z(Zop, Ds)
  ZtC <- Zt_times(Zop, F$C)
  XtDX - (ZtC %*% Matrix::t(ZtC))
}

#' Build Gram for between-class Laplacian using label factors.
#'
#' Z^T L_d Z = Z^T D_d Z + (Z^T C)(Z^T C)^T - (Z^T ell)(Z^T ell)^T
#'
#' @param Zop Operator for kernel blocks.
#' @param F Label factor list.
#' @return Symmetric dense matrix.
#' @examples
#' \donttest{
#' library(Matrix)
#' K1 <- Matrix(rnorm(20), 10, 2)
#' Zop <- make_Zop_from_Ks(list(K1))
#' # F <- label_factors(...)  # Requires label factors
#' # G <- gram_Ld(Zop, F)
#' }
gram_Ld <- function(Zop, F) {
  Dd <- Matrix::Diagonal(x = F$d_d)
  XtDX <- gram_Z_A_Z(Zop, Dd)
  ZtC <- Zt_times(Zop, F$C)
  Ztell <- Zt_times(Zop, matrix(F$ell, ncol = 1L))
  XtDX + (ZtC %*% Matrix::t(ZtC)) - (Ztell %*% Matrix::t(Ztell))
}

#' Simple pivoted Cholesky kernel rank selection.
#'
#' @param X Data matrix (samples x features).
#' @param kernel kernlab kernel object.
#' @param tol Residual tolerance.
#' @param max_rank Optional rank cap.
#' @param diag_fun Optional shortcut for diagonal entries of kernel matrix.
#' @return List with selected indices and factor matrix W (n x r).
#' @examples
#' \donttest{
#' X <- matrix(rnorm(50), 10, 5)
#' kernel <- kernlab::rbfdot(sigma = 0.1)
#' result <- pivoted_chol_kernel_block(X, kernel, tol = 1e-4, max_rank = 5)
#' str(result)
#' }
#' @export
pivoted_chol_kernel_block <- function(X, kernel, tol = 1e-6, max_rank = Inf,
                                      diag_fun = function(X) rep(1, nrow(X))) {
  n <- nrow(X)
  if (n == 0) {
    stop("pivoted_chol_kernel_block: empty input", call. = FALSE)
  }

  diag_residual <- as.numeric(diag_fun(X))
  if (length(diag_residual) != n) {
    stop("diag_fun returned incorrect length", call. = FALSE)
  }

  pivots <- integer(0)
  G <- matrix(0, n, 0)

  repeat {
    if (length(pivots) >= max_rank) break
    pivot_idx <- which.max(diag_residual)
    if (diag_residual[pivot_idx] <= tol) break

    k_col <- as.numeric(kernlab::kernelMatrix(kernel, X, X[pivot_idx, , drop = FALSE]))
    if (ncol(G) > 0) {
      projection <- G %*% Matrix::t(G[pivot_idx, , drop = FALSE])
      k_col <- k_col - as.numeric(projection)
    }
    gamma <- sqrt(max(diag_residual[pivot_idx], 0))
    if (gamma == 0) break
    new_col <- k_col / gamma
    G <- cbind(G, new_col)
    diag_residual <- pmax(diag_residual - new_col^2, 0)
    pivots <- c(pivots, pivot_idx)
  }

  list(indices = pivots, W = G)
}
