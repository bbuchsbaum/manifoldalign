library(testthat)

test_that("assess_ssma returns core diagnostics for reduced solver", {
  skip_if_not_installed("adjoin")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(20260302)
  res <- assess_ssma(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.02,
    structure = "ring",
    n_anchors = 4L,
    holdout_fraction = 0.25,
    n_reps = 1L,
    solvers = "reduced",
    serial = FALSE,
    ssma_knn = 6L,
    ssma_rank_per_domain = 8L,
    decode_k = 3L,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)

  required_cols <- c(
    "solver", "serial_enabled", "runtime_sec", "is_multiblock", "scores_finite",
    "top1_raw_nn", "topk_raw_nn", "paired_mse_raw", "eig_n_positive", "error"
  )
  expect_true(all(required_cols %in% names(res)))

  expect_identical(res$solver[[1]], "reduced")
  expect_false(res$serial_enabled[[1]])
  expect_true(is.finite(res$runtime_sec[[1]]) && res$runtime_sec[[1]] >= 0)
  expect_true(isTRUE(res$is_multiblock[[1]]))
  expect_true(isTRUE(res$scores_finite[[1]]))
  expect_true(is.na(res$error[[1]]))
  expect_true(is.finite(res$top1_raw_nn[[1]]))
  expect_true(is.finite(res$topk_raw_nn[[1]]))
  expect_true(is.finite(res$paired_mse_raw[[1]]))
})

test_that("assess_ssma reports subspace distance when both solvers succeed", {
  skip_if_not_installed("adjoin")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra")

  set.seed(20260303)
  res <- assess_ssma(
    sizes = 14L,
    d = 4L,
    noise_sd = 0.03,
    structure = "grid",
    n_anchors = 5L,
    n_reps = 1L,
    solvers = c("reduced", "operator"),
    serial = FALSE,
    ssma_knn = 7L,
    ssma_rank_per_domain = 10L,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  ok <- is.na(res$error)
  if (!all(c("reduced", "operator") %in% res$solver[ok])) {
    skip("Both solvers did not complete successfully in this environment.")
  }

  dvals <- res$subspace_dist_to_reduced[ok]
  expect_true(any(is.finite(dvals)))
  expect_true(all(dvals[is.finite(dvals)] >= 0))
})
