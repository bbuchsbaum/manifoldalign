# Contracts of gpca_align(): preprocessing lists and eigensolver fallback.

gpca_contract_hd <- function(X1, X2, labels) {
  multidesign::hyperdesign(list(
    domain1 = multidesign::multidesign(X1, data.frame(lbl = labels)),
    domain2 = multidesign::multidesign(X2, data.frame(lbl = labels))
  ))
}

test_that("gpca_align accepts preproc = NULL and NULL list entries", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("multidesign")
  set.seed(31)
  n <- 20
  labels <- factor(rep(c("A", "B"), length.out = n))
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 4), n, 4)
  X1 <- sweep(X1, 2, colMeans(X1))
  X2 <- sweep(X2, 2, colMeans(X2))
  hd <- gpca_contract_hd(X1, X2, labels)

  fit_center <- gpca_align(hd, y = lbl, ncomp = 2)
  # Regression: `proclist[[i]] <- NULL` shrank the list and `names<-` failed.
  fit_null <- gpca_align(hd, y = lbl, ncomp = 2, preproc = NULL)
  fit_mixed <- gpca_align(hd, y = lbl, ncomp = 2,
                          preproc = list(multivarious::center(), NULL))

  # Data are already centred, so the fits coincide (up to column sign).
  for (fit in list(fit_null, fit_mixed)) {
    expect_equal(abs(unname(fit$s)), abs(unname(fit_center$s)), tolerance = 1e-6)
    expect_equal(nrow(fit$v), ncol(X1) + ncol(X2))
  }

  # The stored preprocessor has one identity block per unprocessed domain.
  Xnew <- cbind(X1[1:3, ], X2[1:3, ] + 1)
  expect_equal(unname(multivarious::transform(fit_null$preproc, Xnew)), unname(Xnew))
  tf_mixed <- multivarious::transform(fit_mixed$preproc, Xnew)
  expect_equal(unname(tf_mixed[, 4:7]), unname(Xnew[, 4:7]))
})

test_that("gpca_align metric normalisation does not depend on PRIMME", {
  set.seed(2)
  B <- Matrix::rsparsematrix(30, 30, density = 0.2)
  M <- Matrix::forceSymmetric(Matrix::crossprod(B) + Matrix::Diagonal(30, 0.1))
  expected <- max(eigen(as.matrix(M), symmetric = TRUE, only.values = TRUE)$values)
  expect_equal(manifoldalign:::.gpca_largest_eigenvalue(M), expected, tolerance = 1e-8)

  # If RSpectra fails, the fallback chain still finds the largest eigenvalue:
  # PRIMME when installed, otherwise a dense eigendecomposition.
  local_mocked_bindings(eigs_sym = function(...) stop("mock failure"), .package = "RSpectra")
  if (requireNamespace("PRIMME", quietly = TRUE)) {
    expect_equal(manifoldalign:::.gpca_largest_eigenvalue(M), expected, tolerance = 1e-8)
  }
  local_mocked_bindings(.primme_available = function() FALSE)
  expect_equal(manifoldalign:::.gpca_largest_eigenvalue(M), expected, tolerance = 1e-8)
})

test_that("gpca_align handles a center-then-standardize preprocessing pipeline", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("multidesign")
  # multivarious::standardize() needs matrixStats, which multivarious only suggests.
  skip_if_not_installed("matrixStats")
  set.seed(42)
  n <- 10
  labels <- factor(rep(c("A", "B"), length.out = n))
  hd <- gpca_contract_hd(matrix(rnorm(n * 10), n, 10),
                         matrix(rnorm(n * 10, mean = 1), n, 10), labels)
  custom_preproc <- multivarious::standardize(multivarious::center())
  result_preproc <- gpca_align(hd, y = lbl, ncomp = 2, preproc = custom_preproc)
  expect_true(all(is.finite(result_preproc$s)))
})
