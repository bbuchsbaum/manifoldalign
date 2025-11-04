context("mma_align_multiple: smoke and structure tests")

test_that("mma_align_multiple reference mode works with equal sizes", {
  set.seed(1)
  X1 <- matrix(rnorm(60), 20, 3)
  X2 <- matrix(rnorm(60), 20, 3)
  X3 <- matrix(rnorm(60), 20, 3)

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 5,
                            sigma = 0.7,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            em_max_iter = 10,
                            em_tol = 1e-5)

  expect_true(inherits(res, "mma_align_multiple"))
  expect_true(is.matrix(res$s))
  expect_equal(nrow(res$s), nrow(X1) + nrow(X2) + nrow(X3))
  expect_equal(ncol(res$s), 5)

  rotations <- res$rotations
  expect_equal(length(rotations), 3L)
  for (R in rotations) {
    expect_equal(dim(R), c(5,5))
    I5 <- crossprod(R)
    expect_true(max(abs(I5 - diag(5))) < 1e-6)
  }
})

test_that("mma_align_multiple consensus mode handles unequal sizes and returns template", {
  set.seed(2)
  X1 <- matrix(rnorm(40*3), 40, 3)
  X2 <- matrix(rnorm(35*3), 35, 3)
  X3 <- matrix(rnorm(30*3), 30, 3)

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 1,
                            ncomp = 4,
                            sigma = 0.6,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            match_to = "consensus",
                            consensus_centers = "min",
                            consensus_init = "ref",
                            em_max_iter = 10,
                            final_assignment = "hungarian")

  expect_true(inherits(res, "mma_align_multiple"))
  expect_true(is.matrix(res$s))
  expect_equal(ncol(res$s), 4)
  expect_equal(nrow(res$s), nrow(X1) + nrow(X2) + nrow(X3))

  # consensus template present
  C <- res$consensus
  expect_true(is.matrix(C))
  expect_equal(ncol(C), 4)
  # n_centers should be min sizes = 30
  expect_equal(nrow(C), 30)

  # posteriors shapes: list of (n_centers x n_i)
  alphas <- res$posteriors
  expect_equal(length(alphas), 3L)
  expect_equal(dim(alphas[[1]]), c(30, 40))
  expect_equal(dim(alphas[[2]]), c(30, 35))
  expect_equal(dim(alphas[[3]]), c(30, 30))

  # final assignment has length n_i (may include NAs for padding)
  assign <- res$assignment
  expect_equal(length(assign), 3L)
  expect_length(assign[[1]], 40)
  expect_length(assign[[2]], 35)
  expect_length(assign[[3]], 30)
})

test_that("mma_align_multiple auto-K (ncomp=NULL) selects reasonable K", {
  set.seed(3)
  X1 <- matrix(rnorm(25*3), 25, 3)
  X2 <- matrix(rnorm(22*3), 22, 3)
  X3 <- matrix(rnorm(20*3), 20, 3)

  res <- mma_align_multiple(list(X1, X2, X3),
                            ref_idx = 2,
                            ncomp = NULL,
                            target_var = 0.9,
                            sigma = 0.7,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            em_max_iter = 5)

  expect_true(inherits(res, "mma_align_multiple"))
  K <- res$K
  expect_true(K >= 2)
  expect_true(K <= 64)
  expect_equal(ncol(res$s), K)
})

