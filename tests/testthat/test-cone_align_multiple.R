skip_if_not_installed("manifoldalign")

test_that("cone_align_multiple aligns to the reference space", {
  set.seed(321)
  n <- 18
  d <- 4
  base <- matrix(rnorm(n * d), n, d)
  mats <- lapply(1:3, function(i) base + matrix(rnorm(n * d, sd = 0.05), n, d))
  hd <- create_multi_alignment_hyperdesign(mats)

  result <- manifoldalign::cone_align_multiple(hd, preproc = NULL, ref_idx = 2, ncomp = 4,
                                sigma = 0.9, lambda = 0.05, max_iter = 20, tol = 0.01)

  expect_s3_class(result, "cone_align_multiple")
  expect_s3_class(result, "multiblock_biprojector")

  obj <- unclass(result)
  expect_equal(nrow(obj$s), length(mats) * n)
  expect_equal(ncol(obj$s), 4)

  embeddings <- compute_test_embeddings(hd, ncomp = 4, sigma = 0.9)
  for (g in seq_along(embeddings)) {
    if (g == 2) {
      next
    }
    align_res <- manifoldalign:::`.align_two_embeddings`(
      embeddings[[2]], embeddings[[g]], solver = "linear",
      max_iter = 20, tol = 0.01, lambda = 0.05
    )
    expect_equal(sort(align_res[["P"]]), seq_len(n))
    expect_true(all(abs(crossprod(align_res[["Q"]]) - diag(ncol(align_res[["Q"]]))) < 1e-8))
  }

  scores <- obj$s
  for (g in seq_along(mats)) {
    rows <- ((g - 1) * n + 1):(g * n)
    row_norms <- sqrt(rowSums(scores[rows, , drop = FALSE]^2))
    expect_true(all(row_norms <= 1 + 1e-6))
    expect_true(all(row_norms > 0))
  }
})

test_that("cone_align_multiple improves alignment against known permutations", {
  set.seed(4044)
  n <- 20
  d <- 3
  base <- cbind(
    sin(seq_len(n) / n * 2 * pi),
    cos(seq_len(n) / n * 2 * pi),
    seq_len(n) / n
  )
  perm2 <- c(3:n, 1:2)
  perm3 <- c(n, seq_len(n - 1))
  map2 <- match(seq_len(n), perm2)
  map3 <- match(seq_len(n), perm3)
  mats <- list(
    base,
    base[perm2, , drop = FALSE] + matrix(rnorm(n * d, 0, 0.003), n, d),
    base[perm3, , drop = FALSE] + matrix(rnorm(n * d, 0, 0.003), n, d)
  )
  hd <- create_multi_alignment_hyperdesign(mats)

  embeddings <- compute_test_embeddings(hd, ncomp = 4, sigma = 0.8)

  identity_Q <- diag(ncol(embeddings[[1]]))

  align2 <- manifoldalign:::`.align_two_embeddings`(
    embeddings[[1]], embeddings[[2]], solver = "linear",
    max_iter = 25, tol = 1e-3, lambda = 0.02
  )
  initial_err2 <- alignment_error(embeddings[[1]], embeddings[[2]], seq_len(n), identity_Q)
  final_err2 <- alignment_error(embeddings[[1]], embeddings[[2]], align2[["P"]], align2[["Q"]])

  align3 <- manifoldalign:::`.align_two_embeddings`(
    embeddings[[1]], embeddings[[3]], solver = "linear",
    max_iter = 25, tol = 1e-3, lambda = 0.02
  )
  initial_err3 <- alignment_error(embeddings[[1]], embeddings[[3]], seq_len(n), identity_Q)
  final_err3 <- alignment_error(embeddings[[1]], embeddings[[3]], align3[["P"]], align3[["Q"]])

  expect_lt(final_err2, initial_err2)
  expect_lt(final_err3, initial_err3)
})

test_that("cone_align_multiple rotations match pairwise alignment", {
  set.seed(567)
  n <- 15
  d <- 3
  base <- matrix(rnorm(n * d), n, d)
  mats <- lapply(1:3, function(i) base + matrix(rnorm(n * d, sd = 0.03), n, d))
  hd <- create_multi_alignment_hyperdesign(mats)

  params <- list(preproc = NULL, ncomp = 3, sigma = 0.8, lambda = 0.05, max_iter = 15, tol = 0.02)
  ref_idx <- 1
  embeddings <- compute_test_embeddings(hd, ncomp = 3, sigma = 0.8)

  for (g in seq_along(mats)) {
    if (g == ref_idx) next
    align_pair <- manifoldalign:::`.align_two_embeddings`(
      embeddings[[ref_idx]], embeddings[[g]], solver = "linear",
      max_iter = params$max_iter, tol = params$tol, lambda = params$lambda
    )
    expect_equal(sort(align_pair[["P"]]), seq_len(n))
    expect_true(all(abs(crossprod(align_pair[["Q"]]) - diag(ncol(align_pair[["Q"]]))) < 1e-8))
  }
})

test_that("cone_align_multiple enforces equal node counts", {
  mats <- list(
    matrix(rnorm(5 * 3), 5, 3),
    matrix(rnorm(6 * 3), 6, 3),
    matrix(rnorm(5 * 3), 5, 3)
  )
  expect_error(manifoldalign::cone_align_multiple(mats, preproc = NULL, ncomp = 2),
               "equal node counts", fixed = TRUE)
})

test_that("cone_align_multiple stops when embedding computation fails", {
  n <- 8
  mats <- replicate(3, matrix(rnorm(n * 2), n, 2), simplify = FALSE)
  hd <- create_multi_alignment_hyperdesign(mats)

  with_mocked_bindings(
    compute_cone_embeddings = function(...) list(matrix(0, 1, 1), matrix(0, 1, 1), NULL),
    {
      expect_error(manifoldalign::cone_align_multiple(hd, preproc = NULL, ncomp = 3),
                   "Embedding computation failed", fixed = TRUE)
    }
  , .package = "manifoldalign")
})

test_that("cone_align_multiple list method mirrors hyperdesign method", {
  set.seed(902)
  mats <- lapply(1:3, function(i) matrix(rnorm(20 * 3), 20, 3))
  hd <- create_multi_alignment_hyperdesign(mats)
  params <- list(preproc = NULL, ncomp = 3, sigma = 1.1, lambda = 0.05, max_iter = 10, tol = 0.01)

  res_hd <- do.call(manifoldalign::cone_align_multiple, c(list(hd), params))
  res_list <- do.call(manifoldalign::cone_align_multiple, c(list(mats), params))

  expect_s3_class(res_hd, "cone_align_multiple")
  expect_s3_class(res_list, "cone_align_multiple")
  expect_equal(nrow(res_hd$s), nrow(res_list$s))
  expect_equal(ncol(res_hd$s), ncol(res_list$s))
})
