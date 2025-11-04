test_that("mma_align_multiple w2 signatures run and produce valid rotations", {
  skip_if_missing_cd_deps()
  set.seed(11)
  n <- 30; d <- 5
  t <- runif(n, 0, 2*pi)
  Xbase <- cbind(cos(t), sin(t), matrix(rnorm(n*(d-2), sd = 0.03), n, d-2))
  # Random orthogonal transforms
  orth <- function(p) { M <- matrix(rnorm(p*p), p, p); sv <- svd(M); sv$u %*% t(sv$v) }
  X1 <- Xbase
  X2 <- Xbase %*% orth(d)
  X3 <- Xbase %*% orth(d)

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 5,
                            sigma = 0.6,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            signature_method = "w2",
                            em_max_iter = 15,
                            em_tol = 1e-6)

  expect_true(inherits(res, "mma_align_multiple"))
  expect_equal(ncol(res$s), 5)

  # Rotations orthogonal
  for (R in res$rotations) {
    I5 <- crossprod(R)
    expect_true(max(abs(I5 - diag(ncol(R)))) < 1e-6)
  }
})

