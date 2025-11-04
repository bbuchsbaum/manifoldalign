# tests/testthat/test-generalized_procrustes.R
# -----------------------------------------------------------------------------
# Rigorous unit tests for generalized_procrustes()

library(testthat)
library(manifoldalign)

set.seed(1)

rand_O <- function(d) {
  qr.Q(qr(matrix(rnorm(d * d), d)))
}

estimate_global_rotation <- function(O_hat, O_true) {
  M <- Reduce(`+`, Map(function(Oh, Ot) t(Ot) %*% Oh, O_hat, O_true))
  sv <- svd(M)
  sv$u %*% t(sv$v)
}

rotate_consensus <- function(A_est, Q) {
  Q %*% A_est
}

# -----------------------------------------------------------------------------
test_that("generalized_procrustes recovers template up to global rotation", {
  d <- 3; L <- 6; n <- 4
  template <- matrix(rnorm(d * L), d, L)
  O_true <- replicate(n, rand_O(d), simplify = FALSE)
  A_full <- lapply(O_true, function(O) O %*% template)
  labels_all <- replicate(n, seq_len(L), simplify = FALSE)

  gp <- generalized_procrustes(
    A_list = A_full,
    task_labels_list = labels_all,
    L = L,
    max_iter = 120,
    tol = 1e-7,
    verbose = FALSE
  )

  expect_true(gp$converged || gp$final_diff < 1e-6)
  lapply(gp$O_mats, function(O) {
    expect_equal(crossprod(O), diag(d), tolerance = 1e-6)
  })

  Q <- estimate_global_rotation(gp$O_mats, O_true)
  for (i in seq_len(n)) {
    expect_equal(gp$O_mats[[i]], O_true[[i]] %*% Q, tolerance = 1e-4)
  }

  A_rot <- rotate_consensus(gp$A_est, Q)
  expect_equal(A_rot, template, tolerance = 1e-4)
})

# -----------------------------------------------------------------------------
test_that("partial observations average correctly and mask unobserved tasks", {
  d <- 3; L <- 6; n <- 4
  template <- matrix(rnorm(d * L), d, L)
  O_true <- replicate(n, rand_O(d), simplify = FALSE)

  task_sets <- list(1:3, 2:4, 3:5, 4:5)  # task 6 left unobserved
  A_partial <- lapply(seq_len(n), function(i) O_true[[i]] %*% template[, task_sets[[i]], drop = FALSE])

  gp <- generalized_procrustes(
    A_list = A_partial,
    task_labels_list = task_sets,
    L = L,
    max_iter = 150,
    tol = 1e-7,
    verbose = FALSE
  )

  expect_true(gp$converged || gp$final_diff < 1e-6)

  Q <- estimate_global_rotation(gp$O_mats, O_true)
  A_rot <- rotate_consensus(gp$A_est, Q)

  expect_true(all(is.na(A_rot[, 6])))
  for (task in 1:5) {
    expect_equal(A_rot[, task], template[, task], tolerance = 5e-3)
  }
})

# -----------------------------------------------------------------------------
test_that("small additive noise is tolerated", {
  d <- 3; L <- 5; n <- 3
  template <- matrix(rnorm(d * L), d, L)
  O_true <- replicate(n, rand_O(d), simplify = FALSE)
  noise_sd <- 5e-4
  A_noise <- lapply(O_true, function(O) O %*% template + noise_sd * matrix(rnorm(d * L), d, L))
  labels <- replicate(n, seq_len(L), simplify = FALSE)

  gp <- generalized_procrustes(A_noise, labels, L, max_iter = 150, tol = 1e-6)
  expect_true(gp$converged)

  Q <- estimate_global_rotation(gp$O_mats, O_true)
  A_rot <- rotate_consensus(gp$A_est, Q)
  rel_err <- norm(A_rot - template, "F") / norm(template, "F")
  expect_lt(rel_err, 1e-1)
})

# -----------------------------------------------------------------------------
test_that("duplicate task labels trigger an informative error", {
  d <- 2
  A_bad <- list(matrix(rnorm(d * 2), d, 2), matrix(rnorm(d * 2), d, 2))
  tl_bad <- list(c(1, 1), c(2, 3))
  expect_error(generalized_procrustes(A_bad, tl_bad, L = 3), "duplicate task labels")
})

# -----------------------------------------------------------------------------
test_that("base::svd initialization path works", {
  d <- 3; L <- 5
  template <- matrix(rnorm(d * L), d, L)
  A_list <- lapply(seq_len(3), function(i) rand_O(d) %*% template)
  labels <- replicate(3, seq_len(L), simplify = FALSE)

  gp <- generalized_procrustes(A_list, labels, L, svd_method = "base", max_iter = 80, tol = 1e-6)
  expect_true(gp$converged)
})

# -----------------------------------------------------------------------------
test_that("iterations counter reflects early stopping", {
  d <- 3; L <- 4
  template <- matrix(rnorm(d * L), d, L)
  A_list <- lapply(seq_len(3), function(i) rand_O(d) %*% template)
  labels <- replicate(3, seq_len(L), simplify = FALSE)

  gp <- generalized_procrustes(A_list, labels, L, max_iter = 1, tol = 0)
  expect_false(gp$converged)
  expect_equal(gp$iterations, 1L)
})

# -----------------------------------------------------------------------------
test_that("tightness certificate runs and reports small residuals", {
  d <- 3; L <- 5
  template <- matrix(rnorm(d * L), d, L)
  A_list <- lapply(seq_len(3), function(i) rand_O(d) %*% template)
  labels <- replicate(3, seq_len(L), simplify = FALSE)

  out <- capture.output({
    gp <- generalized_procrustes(A_list, labels, L, max_iter = 80, tol = 1e-6, verbose = TRUE)
  })
  expect_true(any(grepl("Tightness certificate", out)))
  expect_true(gp$converged)
})

# -----------------------------------------------------------------------------
test_that("preallocated sparse triplets handle moderately large problems", {
  d <- 3; L <- 10; n <- 5
  template <- matrix(rnorm(d * L), d, L)
  O_true <- replicate(n, rand_O(d), simplify = FALSE)

  labels <- vector("list", n)
  A_list <- vector("list", n)
  for (i in seq_len(n)) {
    tasks_i <- sort(sample(L, size = sample(4:6, 1)))
    labels[[i]] <- tasks_i
    A_list[[i]] <- O_true[[i]] %*% template[, tasks_i, drop = FALSE]
  }

  t0 <- Sys.time()
  gp <- generalized_procrustes(A_list, labels, L, max_iter = 150, tol = 1e-5)
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")

  expect_true(gp$converged)
  expect_lt(elapsed, 5)
  lapply(gp$O_mats, function(O) expect_equal(crossprod(O), diag(d), tolerance = 1e-4))
})

# -----------------------------------------------------------------------------
test_that("list interface enforces equal task counts", {
  mats <- list(matrix(rnorm(6), 3, 2), matrix(rnorm(6), 3, 2))
  expect_no_error(generalized_procrustes(mats, L = NULL))

  mats_bad <- list(matrix(rnorm(6), 3, 2), matrix(rnorm(9), 3, 3))
  expect_error(generalized_procrustes(mats_bad), "same number of columns")
})

# -----------------------------------------------------------------------------
test_that("single subject is rejected", {
  d <- 3; L <- 5
  A_single <- list(matrix(rnorm(d * L), d, L))
  labels_single <- list(seq_len(L))
  expect_error(generalized_procrustes(A_single, labels_single, L), "at least 2 subjects")
})
