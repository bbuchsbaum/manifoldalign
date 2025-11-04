test_that("w2_cost_matrix_1d returns valid costs and signs", {
  skip_if_missing_cd_deps()
  set.seed(123)
  n <- 50; K <- 6
  Uref <- matrix(rnorm(n * K), n, K)
  Ug   <- matrix(rnorm(n * K), n, K)
  wr <- runif(n); wg <- runif(n)

  # Call C++ kernel
  out_cpp <- w2_cost_matrix_1d(Uref, Ug, wr, wg)
  C_cpp <- as.matrix(out_cpp$C)
  S_cpp <- as.matrix(out_cpp$sign)
  expect_equal(dim(C_cpp), c(K, K))
  expect_equal(dim(S_cpp), c(K, K))
  expect_true(all(is.finite(C_cpp)))
  expect_true(all(C_cpp >= -1e-12))
  expect_true(all(S_cpp %in% c(-1L, 1L)))
  # Self-consistency: identical inputs should give near-zero diagonal and +1 signs
  out_id <- w2_cost_matrix_1d(Uref, Uref, wr, wr)
  C_id <- as.matrix(out_id$C); S_id <- as.matrix(out_id$sign)
  expect_true(max(diag(C_id)) < 1e-12)
  expect_true(all(diag(S_id) == 1L))
})
