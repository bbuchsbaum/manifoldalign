# Cross-method contracts for the spectral graph aligners (CONE-Align, CONE
# multiple, GRASP, MMA, coupled diagonalization): fitted preprocessing,
# out-of-sample projection, score row identity and loadings.

block_rows <- function(n) {
  ends <- cumsum(n)
  Map(seq.int, c(1L, head(ends, -1L) + 1L), ends)
}

# Projection of new rows through a fitted multiblock result must replay the
# per-block centring fitted on that block's training rows.
expect_centred_projection <- function(res, X_train, block) {
  expect_false(is.null(res$preproc))
  newx <- X_train[seq_len(min(4L, nrow(X_train))), , drop = FALSE] + 0.25
  p <- multivarious::project_block(res, newx, block = block)
  v_b <- res$v[res$block_indices[[block]], , drop = FALSE]
  expect_equal(unname(as.matrix(p)),
               unname(sweep(newx, 2, colMeans(X_train)) %*% v_b),
               tolerance = 1e-8)
}

# X_i v_i is the least-squares fit of the aligned scores: exact when X_i has
# full row rank, otherwise the residual is orthogonal to the columns of X_i.
expect_regression_loadings <- function(res, X_list, centre = TRUE) {
  rows <- block_rows(vapply(X_list, nrow, integer(1)))
  for (i in seq_along(X_list)) {
    Xi <- if (centre) sweep(X_list[[i]], 2, colMeans(X_list[[i]])) else X_list[[i]]
    vi <- res$v[res$block_indices[[i]], , drop = FALSE]
    Si <- res$s[rows[[i]], , drop = FALSE]
    resid <- Si - Xi %*% vi
    expect_lt(max(abs(crossprod(Xi, resid))), 1e-8)
  }
}

test_that("CONE-Align keeps domain-2 rows in input order and reports the assignment", {
  fx <- generate_permuted_alignment(n = 16, d = 3, noise_sd = 0.02, seed = 1001)
  hd <- create_alignment_hyperdesign(fx$X1, fx$X2)
  res <- suppressWarnings(cone_align(hd, preproc = NULL, ncomp = 4, sigma = 0.9,
                                     lambda = 0.05, max_iter = 20, tol = 0.01))
  emb <- compute_test_embeddings(hd, ncomp = 4, sigma = 0.9)
  n <- nrow(fx$X1)

  expect_false(identical(res$assignment, seq_len(n)))
  expect_equal(res$s[n + seq_len(n), ], emb[[2]], tolerance = 1e-10)
  expect_equal(res$s[seq_len(n), ], emb[[1]] %*% res$rotation[[1]], tolerance = 1e-10)
  expect_true(is.logical(res$converged) && length(res$converged) == 1L)
  expect_true(res$iterations >= 1L)
})

test_that("CONE-Align stores the fitted preprocessor and projects new rows", {
  fx <- generate_permuted_alignment(n = 18, d = 3, noise_sd = 0.02, seed = 31)
  X1 <- fx$X1 + 5
  X2 <- fx$X2 - 2
  res <- suppressWarnings(cone_align(list(X1, X2), ncomp = 3))
  expect_centred_projection(res, X1, 1)
  expect_centred_projection(res, X2, 2)
  expect_regression_loadings(res, list(X1, X2))
})

test_that("CONE-Align loadings reproduce the aligned scores for full-row-rank features", {
  set.seed(5)
  n <- 12
  X1 <- matrix(rnorm(n * 30), n, 30)
  X2 <- X1[c(2:n, 1L), ] + matrix(rnorm(n * 30, sd = 0.01), n, 30)
  res <- suppressWarnings(cone_align(list(X1, X2), preproc = NULL, ncomp = 3))
  rows <- block_rows(c(n, n))
  for (i in 1:2) {
    Xi <- list(X1, X2)[[i]]
    vi <- res$v[res$block_indices[[i]], , drop = FALSE]
    expect_equal(Xi %*% vi, res$s[rows[[i]], , drop = FALSE], tolerance = 1e-8)
  }
})

test_that("CONE-Align rejects the unavailable auction solver and is quiet by default", {
  fx <- generate_permuted_alignment(n = 12, d = 3, noise_sd = 0.02, seed = 3)
  hd <- create_alignment_hyperdesign(fx$X1, fx$X2)
  expect_error(cone_align(hd, solver = "auction"), "not available")
  expect_error(cone_align(list(fx$X1, fx$X2), solver = "auction"), "not available")
  expect_error(cone_align_multiple(list(fx$X1, fx$X2, fx$X1), solver = "auction"),
               "not available")
  expect_error(cone_align(list(fx$X1, fx$X2), anchors = node_id), "hyperdesign")
  expect_no_message(suppressWarnings(cone_align(hd, preproc = NULL, ncomp = 3)))
  expect_message(suppressWarnings(cone_align(hd, preproc = NULL, ncomp = 3, verbose = TRUE)),
                 "iterations")
})

test_that("CONE multiple stores preprocessing and regression loadings", {
  set.seed(8)
  base <- matrix(rnorm(20 * 3), 20, 3)
  mats <- list(base + 3, base[sample(20), ] - 1, base[sample(20), ] + 0.5)
  res <- suppressWarnings(cone_align_multiple(mats, ncomp = 3))
  for (b in 1:3) expect_centred_projection(res, mats[[b]], b)
  expect_regression_loadings(res, mats)
  expect_length(res$converged, 3L)
})

test_that("MMA stores preprocessing, keeps row order and fits regression loadings", {
  set.seed(12)
  mats <- lapply(c(18L, 16L, 20L), function(n) matrix(rnorm(n * 3), n, 3) + 4)
  res <- suppressWarnings(mma_align_multiple(mats, ncomp = 3, em_max_iter = 10))
  for (b in 1:3) expect_centred_projection(res, mats[[b]], b)
  expect_regression_loadings(res, mats)
  expect_equal(nrow(res$s), sum(vapply(mats, nrow, integer(1))))
})

test_that("GRASP keeps both domains in input order with unequal node counts", {
  set.seed(14)
  X1 <- matrix(rnorm(14 * 3), 14, 3)
  X2 <- matrix(rnorm(12 * 3), 12, 3)
  hd <- structure(list(a = list(x = X1, design = data.frame(i = 1:14)),
                       b = list(x = X2, design = data.frame(i = 1:12))),
                  class = "hyperdesign")
  # ncomp exceeds min(n) - 1: bases must still share one width
  res <- suppressWarnings(grasp(hd, ncomp = 15, q_descriptors = 8))
  expect_equal(nrow(res$s), 14 + 12)
  expect_false(any(rowSums(abs(res$s)) == 0))
  expect_length(res$assignment, 14)
  expect_equal(sum(is.na(res$assignment)), 2L)
  ok <- res$assignment[!is.na(res$assignment)]
  expect_true(all(ok %in% 1:12) && !anyDuplicated(ok))
  expect_centred_projection(res, X1, 1)
  expect_centred_projection(res, X2, 2)
  expect_regression_loadings(res, list(X1, X2))
})

test_that("GRASP score rows do not depend on the order of domain-2 rows", {
  set.seed(15)
  X1 <- matrix(rnorm(20 * 3), 20, 3)
  X2 <- X1 + matrix(rnorm(20 * 3, sd = 0.01), 20, 3)
  hd <- create_alignment_hyperdesign(X1, X2)
  res <- suppressWarnings(grasp(hd, preproc = NULL, ncomp = 5, q_descriptors = 10))
  # block 2 of the scores is domain 2's own aligned basis, row for row
  n <- 20
  bases <- manifoldalign:::compute_grasp_basis(unclass(hd), 5, TRUE)
  embed2 <- as.matrix(bases[[2]]$vectors %*% res$rotation %*% res$mapping_matrix)
  expect_equal(res$s[n + seq_len(n), ], embed2, tolerance = 1e-8)
})

test_that("coupled diagonalization concatenates preprocessors by feature block", {
  set.seed(16)
  X1 <- matrix(rnorm(30 * 3), 30, 3) + 2
  X2 <- matrix(rnorm(30 * 4), 30, 4) - 1
  hd <- structure(list(a = list(x = X1, design = data.frame(i = 1:30)),
                       b = list(x = X2, design = data.frame(i = 1:30))),
                  class = "hyperdesign")
  res <- suppressWarnings(coupled_diagonalization(hd, ncomp = 3, ncomp_per_domain = 5,
                                                  max_iter = 20))
  expect_centred_projection(res, X1, 1)
  expect_centred_projection(res, X2, 2)
  full <- multivarious::transform(res$preproc, cbind(X1, X2))
  expect_equal(full, cbind(sweep(X1, 2, colMeans(X1)), sweep(X2, 2, colMeans(X2))),
               ignore_attr = TRUE)
})

test_that("graph Laplacians are symmetrised before the symmetric eigensolve", {
  A <- Matrix::sparseMatrix(i = c(1, 2, 3, 4), j = c(2, 3, 4, 1), x = c(1, 2, 3, 4),
                            dims = c(4, 4))
  for (norm in c(TRUE, FALSE)) {
    L <- manifoldalign:::graph_laplacian(A, normalized = norm)
    expect_true(Matrix::isSymmetric(L))
    As <- (A + Matrix::t(A)) / 2
    L_ref <- manifoldalign:::graph_laplacian(As, normalized = norm)
    expect_equal(as.matrix(L), as.matrix(L_ref))
  }
})
