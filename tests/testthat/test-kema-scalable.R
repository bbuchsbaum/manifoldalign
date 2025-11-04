test_that("label_factors captures class structure", {
  labels <- c("A", "B", NA, "B", "A")
  F <- label_factors(labels)

  expect_s4_class(F$C, "dgCMatrix")
  expect_equal(length(F$ell), length(labels))
  expect_equal(sum(F$ell), 4)
  expect_true(all(F$d_s[F$ell == 1] > 0))
  expect_true(all(F$d_d[F$ell == 1] >= 0))
})

test_that("low-rank Gram builders match explicit matrices", {
  set.seed(123)
  labels <- c("A", "A", "B", NA)
  F <- label_factors(labels)

  # Build simple kernel block and operator
  K <- Matrix::Matrix(matrix(runif(16), nrow = 4), sparse = TRUE)
  Ks <- list(K)
  Zop <- make_Zop_from_Ks(Ks)
  Z <- Matrix::bdiag(Ks)

  # Explicit Laplacians from factors for comparison
  C <- F$C
  Ds <- Matrix::Diagonal(x = F$d_s)
  Dd <- Matrix::Diagonal(x = F$d_d)
  ell_vec <- Matrix::Matrix(matrix(F$ell, ncol = 1L), sparse = TRUE)

  Ls_explicit <- Ds - C %*% Matrix::t(C)
  Ld_explicit <- Dd + C %*% Matrix::t(C) - ell_vec %*% Matrix::t(ell_vec)

  gram_ls_expected <- Matrix::crossprod(Z, Ls_explicit %*% Z)
  gram_ld_expected <- Matrix::crossprod(Z, Ld_explicit %*% Z)

  expect_equal(gram_Ls(Zop, F), Matrix::Matrix(gram_ls_expected, sparse = FALSE), tolerance = 1e-10)
  expect_equal(gram_Ld(Zop, F), Matrix::Matrix(gram_ld_expected, sparse = FALSE), tolerance = 1e-10)
})
