test_that("mma_align_multiple works with and without hypersphere normalization", {
  skip_if_missing_cd_deps()
  set.seed(7)
  n <- 35; d <- 4
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- X1 %*% {M <- matrix(rnorm(d*d), d, d); svd(M)$u %*% t(svd(M)$v)} + matrix(rnorm(n * d, sd = 0.02), n, d)
  X3 <- X1 %*% {M <- matrix(rnorm(d*d), d, d); svd(M)$u %*% t(svd(M)$v)} + matrix(rnorm(n * d, sd = 0.02), n, d)

  res_h <- mma_align_multiple(list(X1, X2, X3),
                              ref_idx = 1,
                              ncomp = 5,
                              sigma = 0.6,
                              embedding = "ctd",
                              normalize = "hypersphere",
                              em_max_iter = 15,
                              final_assignment = "hungarian")

  res_n <- mma_align_multiple(list(X1, X2, X3),
                              ref_idx = 1,
                              ncomp = 5,
                              sigma = 0.6,
                              embedding = "ctd",
                              normalize = "none",
                              em_max_iter = 15,
                              final_assignment = "hungarian")

  expect_true(inherits(res_h, "mma_align_multiple"))
  expect_true(inherits(res_n, "mma_align_multiple"))
  expect_equal(res_h$K, res_n$K)

  # EM stats present; iterations should be positive for non-ref domains
  it_h <- vapply(res_h$em_stats[-1], function(z) as.integer(z$iterations), integer(1))
  it_n <- vapply(res_n$em_stats[-1], function(z) as.integer(z$iterations), integer(1))
  expect_true(all(it_h > 0))
  expect_true(all(it_n > 0))
})
