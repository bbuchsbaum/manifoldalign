# Seeded helpers of multiscale_manifold_align() must be deterministic and must
# not disturb the caller's RNG stream.

test_that("label correspondence sampling restores the caller's RNG", {
  labels <- list(rep(c("A", "B", "C"), each = 6), rep(c("A", "B", "C"), times = 6))
  set.seed(77)
  expected <- runif(2)
  set.seed(77)
  c1 <- manifoldalign:::.msa_build_label_correspondences(labels, max_pairs_per_label = 3,
                                                         weight = 1, seed = 5)
  expect_equal(runif(2), expected)
  c2 <- manifoldalign:::.msa_build_label_correspondences(labels, max_pairs_per_label = 3,
                                                         weight = 1, seed = 5)
  expect_identical(c1, c2)
  expect_equal(nrow(c1), 9L)
  expect_true(all(labels[[1]][c1$index_i] == labels[[2]][c1$index_j]))
})

test_that("randomized diffusion eigensolver restores the caller's RNG", {
  set.seed(1)
  B <- matrix(rnorm(40 * 40), 40, 40)
  A <- crossprod(B) / 40
  ctrl <- list(seed = 3L, randomized_sketch_size = 10L, randomized_oversample = 5L,
               randomized_power_iters = 2L)
  set.seed(88)
  expected <- runif(2)
  set.seed(88)
  e1 <- manifoldalign:::.msa_diffusion_eig_randomized(A, k = 3L, ctrl = ctrl)
  expect_equal(runif(2), expected)
  e2 <- manifoldalign:::.msa_diffusion_eig_randomized(A, k = 3L, ctrl = ctrl)
  expect_identical(e1, e2)
  # Randomized sketch: approximate leading eigenvalues.
  expect_equal(e1$values, eigen(A, symmetric = TRUE)$values[1:3], tolerance = 1e-2)
})

test_that("the joint Laplacian builder used by the fit is the adjacency-based one", {
  A <- Matrix::sparseMatrix(i = c(1, 2, 2, 3), j = c(2, 1, 3, 2), x = c(1, 1, 2, 2), dims = c(4, 4))
  lap <- manifoldalign:::.msa_build_joint_laplacian(A)
  expect_named(lap, c("laplacian", "degree"))
  # Isolated node 4 gets unit degree; L = D - A.
  expect_equal(lap$degree, c(1, 3, 2, 1))
  expect_equal(as.matrix(lap$laplacian), diag(c(1, 3, 2, 1)) - as.matrix(A))
})
