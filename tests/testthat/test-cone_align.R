skip_if_not_installed("manifoldalign")

nrow_block <- function(n) seq_len(n)

test_that("cone_align returns normalized aligned embeddings", {
  fixture <- generate_permuted_alignment(n = 14, d = 3, noise_sd = 0.02, seed = 1001)
  hd <- create_alignment_hyperdesign(fixture$X1, fixture$X2)
  result <- manifoldalign::cone_align(hd, preproc = NULL, ncomp = 5, sigma = 0.9,
                       lambda = 0.05, max_iter = 20, tol = 0.01)

  expect_s3_class(result, "cone_align")
  expect_s3_class(result, "multiblock_biprojector")
  expect_equal(nrow(result$s), 2 * nrow(fixture$X1))
  expect_equal(ncol(result$s), 5)

  embeddings <- compute_test_embeddings(hd, ncomp = 5, sigma = 0.9)
  iter_res <- suppressMessages(
    manifoldalign:::cone_align_iterate(embeddings, solver = "linear",
                                       max_iter = 20, tol = 0.01, lambda = 0.05)
  )

  P_vec <- iter_res[["P"]]
  Q_mat <- iter_res[["Q"]]
  expect_equal(sort(P_vec), seq_len(nrow(fixture$X1)))
  expect_equal(ncol(Q_mat), 5)
  expect_true(all(abs(crossprod(Q_mat) - diag(ncol(Q_mat))) < 1e-8))

  scores <- result$s
  n <- nrow(fixture$X1)
  scores_1 <- scores[nrow_block(n), , drop = FALSE]
  scores_2 <- scores[n + nrow_block(n), , drop = FALSE]
  row_norms_1 <- sqrt(rowSums(scores_1^2))
  row_norms_2 <- sqrt(rowSums(scores_2^2))

  expect_true(all(row_norms_1 <= 1 + 1e-6))
  expect_true(all(row_norms_2 <= 1 + 1e-6))
  expect_true(all(row_norms_1 > 0))
  expect_true(all(row_norms_2 > 0))
})

test_that("cone_align recovers the identity mapping on identical graphs", {
  fixture <- generate_identity_alignment(n = 20, d = 4, seed = 2024)
  hd <- create_alignment_hyperdesign(fixture$X1, fixture$X2)
  embeddings <- compute_test_embeddings(hd, ncomp = 5, sigma = 0.9)
  iter_res <- suppressMessages(
    manifoldalign:::cone_align_iterate(embeddings, solver = "linear",
                                       max_iter = 25, tol = 1e-3, lambda = 0.05)
  )

  expect_equal(iter_res[["P"]], seq_len(nrow(fixture$X1)))
})

test_that("cone_align iteration reduces alignment error", {
  set.seed(1307)
  n <- 16
  X1 <- matrix(rnorm(n * 4), n, 4)
  perm <- sample(n)
  X2 <- X1[perm, , drop = FALSE] + matrix(rnorm(n * 4, 0, 0.01), n)

  hd <- create_alignment_hyperdesign(X1, X2)
  embeddings <- compute_test_embeddings(hd, ncomp = 5, sigma = 0.9)

  initial_err <- alignment_error(embeddings[[1]], embeddings[[2]], seq_len(n), diag(ncol(embeddings[[1]])))
  iter_res <- suppressMessages(
    manifoldalign:::cone_align_iterate(embeddings, solver = "linear",
                                       max_iter = 20, tol = 1e-3, lambda = 0.05)
  )
  P_vec <- iter_res[["P"]]
  Q_mat <- iter_res[["Q"]]
  final_err <- alignment_error(embeddings[[1]], embeddings[[2]], P_vec, Q_mat)

  expect_lt(final_err, initial_err)
  expect_true(all(abs(crossprod(Q_mat) - diag(ncol(Q_mat))) < 1e-6))
})

test_that("solve_procrustes_cone returns the correct rotation", {
  Z1 <- matrix(rnorm(60), nrow = 20, ncol = 3)
  theta <- pi / 4
  Q_true <- matrix(c(cos(theta), -sin(theta), 0,
                     sin(theta),  cos(theta), 0,
                     0,           0,          1), 3, 3, byrow = TRUE)
  Z2 <- Z1 %*% Q_true

  Q_est <- manifoldalign:::solve_procrustes_cone(Z1, Z2, seq_len(nrow(Z1)), lambda = 0)
  expect_equal(Q_est, Q_true, tolerance = 1e-6)
  expect_equal(t(Q_est) %*% Q_est, diag(3), tolerance = 1e-8)
})

test_that("cone_align enforces equal node counts", {
  X1 <- matrix(rnorm(5 * 3), 5, 3)
  X2 <- matrix(rnorm(4 * 3), 4, 3)
  expect_error(manifoldalign::cone_align(list(X1, X2), preproc = NULL, ncomp = 2),
               "equal node counts", fixed = TRUE)
})

test_that("cone_align stops when embedding computation fails", {
  set.seed(11)
  n <- 6
  X <- matrix(rnorm(n * 2), n, 2)
  hd <- create_alignment_hyperdesign(X, X)

  with_mocked_bindings(
    compute_cone_embeddings = function(...) list(matrix(0, 1, 1), NULL),
    {
      expect_error(manifoldalign::cone_align(hd, preproc = NULL, ncomp = 2),
                   "Embedding computation failed", fixed = TRUE)
    }
  , .package = "manifoldalign")
})

test_that("cone_align parameter validation rejects invalid inputs", {
  n <- 8
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  hd <- create_alignment_hyperdesign(X1, X2)

  expect_error(manifoldalign::cone_align(hd, ncomp = 0), "ncomp > 0")
  expect_error(manifoldalign::cone_align(hd, sigma = 0), "sigma > 0")
  expect_error(manifoldalign::cone_align(hd, lambda = -1), "lambda >= 0")
  expect_error(manifoldalign::cone_align(hd, max_iter = 0), "max_iter > 0")
  expect_error(manifoldalign::cone_align(hd, tol = 0), "tol > 0")
  expect_error(manifoldalign::cone_align(hd, solver = "bogus"), "should be one of")
})

test_that("cone_align solver choices produce valid permutations", {
  set.seed(55)
  n <- 9
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- X1[sample(n), ]
  hd <- create_alignment_hyperdesign(X1, X2)

  embeddings <- compute_test_embeddings(hd, ncomp = 3, sigma = 0.9)

  res_linear <- suppressMessages(
    manifoldalign:::cone_align_iterate(embeddings, solver = "linear",
                                       max_iter = 10, tol = 0.01, lambda = 0.05)
  )
  res_auction <- suppressMessages(
    manifoldalign:::cone_align_iterate(embeddings, solver = "auction",
                                       max_iter = 10, tol = 0.01, lambda = 0.05)
  )

  expect_equal(sort(res_linear[["P"]]), seq_len(n))
  expect_equal(sort(res_auction[["P"]]), seq_len(n))
})
