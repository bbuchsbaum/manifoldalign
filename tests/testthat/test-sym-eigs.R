# Shared symmetric eigensolver used by the graph aligners: spectrum ends,
# optional PRIMME, and explicit handling of non-converged solves.

test_that("dense fallback returns the requested end of the spectrum", {
  L <- diag(c(5, 1, 3, 2, 4, 6, 7))
  sa <- manifoldalign:::.dense_sym_eigs(L, 3, which = "SA")
  la <- manifoldalign:::.dense_sym_eigs(L, 3, which = "LA")
  expect_equal(sa$values, c(1, 2, 3))
  expect_equal(la$values, c(7, 6, 5))
  expect_equal(abs(sa$vectors[2, 1]), 1)
})

# Sparse symmetric matrix with a well-separated spectrum (fast iterative solves).
sym_test_matrix <- function(n) {
  A <- Matrix::bandSparse(n, k = 1, diagonals = list(rep(0.05, n - 1)), symmetric = TRUE)
  as(Matrix::Diagonal(x = seq_len(n)) + A, "CsparseMatrix")
}

test_that(".sym_eigs agrees across backends and sorts its output", {
  L <- sym_test_matrix(120)
  ref <- sort(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
  rs <- manifoldalign:::.sym_eigs(L, 6, which = "SA", backend = "rspectra", tol = 1e-12)
  expect_identical(rs$backend, "rspectra")
  expect_equal(rs$values, ref[1:6], tolerance = 1e-8)
  la <- manifoldalign:::.sym_eigs(L, 4, which = "LA", backend = "rspectra", tol = 1e-12)
  expect_equal(la$values, rev(ref)[1:4], tolerance = 1e-8)
  # residual check on the returned eigenvectors
  expect_lt(max(abs(as.matrix(L %*% rs$vectors) - sweep(rs$vectors, 2, rs$values, "*"))), 1e-6)

  if (manifoldalign:::.primme_available()) {
    pr <- manifoldalign:::.sym_eigs(L, 6, which = "SA", backend = "primme")
    expect_identical(pr$backend, "primme")
    expect_equal(pr$values, ref[1:6], tolerance = 1e-6)
  }
})

test_that(".sym_eigs falls back when PRIMME is unavailable", {
  L <- sym_test_matrix(40)
  ref <- sort(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
  testthat::local_mocked_bindings(.primme_available = function() FALSE,
                                  .package = "manifoldalign")
  small <- manifoldalign:::.sym_eigs(L, 5, which = "SA", backend = "primme")
  expect_identical(small$backend, "dense")
  expect_equal(small$values, ref[1:5], tolerance = 1e-10)
  big <- manifoldalign:::.sym_eigs(sym_test_matrix(300), 5, which = "SA", backend = "primme",
                                   tol = 1e-10)
  expect_identical(big$backend, "rspectra")
})

test_that("non-converged iterative solves are recomputed densely or rejected", {
  L <- sym_test_matrix(200)
  ref <- sort(eigen(as.matrix(L), symmetric = TRUE, only.values = TRUE)$values)
  fb <- manifoldalign:::.sym_eigs(L, 8, which = "SA", backend = "rspectra",
                                  tol = 1e-14, maxitr = 1L)
  expect_identical(fb$backend, "dense_fallback")
  expect_equal(fb$values, ref[1:8], tolerance = 1e-10)
  expect_error(manifoldalign:::.sym_eigs(L, 8, which = "SA", backend = "rspectra",
                                         tol = 1e-14, maxitr = 1L,
                                         dense_fallback_max_n = 50L),
               "converged only")
})

test_that("CONE and MMA embeddings do not depend on PRIMME availability", {
  set.seed(19)
  X <- matrix(rnorm(30 * 3), 30, 3)
  dom <- list(x = X)
  e_default <- manifoldalign:::compute_embedding(dom, ncomp = 4, sigma = 0.9,
                                                 use_laplacian = TRUE, knn = NULL)
  testthat::local_mocked_bindings(.primme_available = function() FALSE,
                                  .package = "manifoldalign")
  e_fallback <- manifoldalign:::compute_embedding(dom, ncomp = 4, sigma = 0.9,
                                                  use_laplacian = TRUE, knn = NULL)
  # eigenvectors are defined up to sign
  sgn <- sign(colSums(e_default * e_fallback))
  expect_equal(sweep(e_fallback, 2, sgn, "*"), e_default, tolerance = 1e-5)

  strata <- list(a = list(x = X), b = list(x = X[30:1, ]), c = list(x = X + 0.01))
  m_fb <- manifoldalign:::compute_mma_embeddings(strata, ncomp = 3, sigma = 0.9, knn = NULL,
                                                 embedding = "ctd", normalize = "none",
                                                 target_var = 0.95)
  expect_length(m_fb, 3L)
  expect_true(all(vapply(m_fb, function(E) all(is.finite(E)), logical(1))))
})

test_that("compute_laplacian_basis accepts SM/LM synonyms and rejects others", {
  set.seed(2)
  A <- Matrix::rsparsematrix(30, 30, density = 0.2, rand.x = function(n) runif(n))
  A <- A + Matrix::t(A)
  sa <- compute_laplacian_basis(A, k = 4, which = "SA")
  sm <- compute_laplacian_basis(A, k = 4, which = "SM")
  expect_equal(sm$values, sa$values, tolerance = 1e-8)
  expect_true(all(diff(sa$values) >= -1e-12))
  expect_error(compute_laplacian_basis(A, k = 4, which = "BE"), "must be one of")
})

test_that("GRASP adjacency bases use the largest eigenvalues with or without PRIMME", {
  set.seed(23)
  X <- matrix(rnorm(25 * 3), 25, 3)
  strata <- list(list(x = X), list(x = X[25:1, ]))
  # k close to n selects the dense branch (PRIMME missed a repeated eigenvalue here)
  b_default <- manifoldalign:::compute_grasp_basis(strata, ncomp = 22, use_laplacian = FALSE)
  testthat::local_mocked_bindings(.primme_available = function() FALSE,
                                  .package = "manifoldalign")
  b_fallback <- manifoldalign:::compute_grasp_basis(strata, ncomp = 22, use_laplacian = FALSE)
  for (i in 1:2) {
    expect_true(all(diff(b_default[[i]]$values) <= 1e-10))
    expect_equal(b_fallback[[i]]$values, b_default[[i]]$values, tolerance = 1e-6)
  }
})
