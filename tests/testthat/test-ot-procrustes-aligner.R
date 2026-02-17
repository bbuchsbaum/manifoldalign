# tests/testthat/test-ot-procrustes-aligner.R

library(testthat)

rand_orthogonal <- function(d, force_SO = FALSE) {
  if (d == 1L) return(matrix(1, 1, 1))
  A <- matrix(rnorm(d * d), d, d)
  sv <- svd(A)
  R <- sv$u %*% t(sv$v)
  if (isTRUE(force_SO) && det(R) < 0) {
    sv$u[, d] <- -sv$u[, d]
    R <- sv$u %*% t(sv$v)
  }
  R
}

test_that("ot_procrustes_aligner fits and returns an orthogonal transform", {
  set.seed(101)
  n <- 30
  d <- 4
  X <- matrix(rnorm(n * d), n, d)
  R_true <- rand_orthogonal(d)
  Y <- X %*% R_true

  algo <- ot_procrustes_aligner()
  fit <- fit_pair(
    algo,
    X,
    Y,
    epsilon0 = 1,
    epsilon_min = 0.05,
    decay = 0.9,
    max_iter = 60,
    tol = 1e-6,
    sinkhorn_max_iter = 2000,
    sinkhorn_tol = 1e-9,
    store_transport = FALSE
  )

  expect_s3_class(fit, "ot_procrustes_pair_fit")
  expect_true(is.matrix(fit$rotation))
  expect_equal(dim(fit$rotation), c(d, d))

  # Check orthogonality: R^T R ≈ I
  RtR <- t(fit$rotation) %*% fit$rotation
  expect_lt(norm(RtR - diag(d), type = "F"), 1e-6)
})

test_that("ot_procrustes_aligner recovers a known rotation on synthetic data", {
  set.seed(202)
  n <- 40
  d <- 3
  X <- matrix(rnorm(n * d), n, d)
  R_true <- rand_orthogonal(d)
  Y <- X %*% R_true

  algo <- ot_procrustes_aligner()
  fit <- fit_pair(
    algo,
    X,
    Y,
    epsilon0 = 1,
    epsilon_min = 1e-3,
    decay = 0.8,
    n_init = 12,
    seed = 202,
    max_iter = 80,
    tol = 1e-6,
    sinkhorn_max_iter = 2000,
    sinkhorn_tol = 1e-9,
    init = "identity",
    store_transport = FALSE
  )

  R_hat <- fit$rotation
  # Compare via relative rotation: should be close to identity.
  rel <- R_hat %*% t(R_true)
  expect_lt(norm(rel - diag(d), type = "F"), 0.2)
})

test_that("ot_procrustes_aligner can optionally store a transport plan", {
  set.seed(303)
  n <- 15
  d <- 2
  X <- matrix(rnorm(n * d), n, d)
  R_true <- rand_orthogonal(d)
  Y <- X %*% R_true

  algo <- ot_procrustes_aligner()
  fit <- fit_pair(
    algo,
    X,
    Y,
    epsilon0 = 0.5,
    epsilon_min = 0.1,
    decay = 1,
    max_iter = 5,
    sinkhorn_max_iter = 2000,
    store_transport = TRUE
  )

  P <- fit$transport
  expect_true(is.matrix(P))
  expect_equal(dim(P), c(n, n))
  expect_equal(rowSums(P), rep(1 / n, n), tolerance = 1e-3)
  expect_equal(colSums(P), rep(1 / n, n), tolerance = 1e-3)
})

test_that("ot_procrustes_aligner integrates with align_many()", {
  set.seed(404)
  n <- 25
  d <- 3
  X0 <- matrix(rnorm(n * d), n, d)
  R2 <- rand_orthogonal(d)
  R3 <- rand_orthogonal(d)
  domains <- list(
    X0,
    X0 %*% R2,
    X0 %*% R3
  )

  res <- align_many(
    domains,
    ot_procrustes_aligner(),
    graph = "complete",
    consensus = "sync",
    k = d,
    parallel = FALSE,
    epsilon0 = 1,
    epsilon_min = 1e-3,
    decay = 0.8,
    n_init = 12,
    seed = 404,
    max_iter = 80,
    sinkhorn_max_iter = 2000,
    store_transport = FALSE
  )

  expect_true(is.list(res$transforms))
  expect_length(res$transforms, 3)

  # Each transform should be orthogonal.
  for (i in seq_along(res$transforms)) {
    op <- as.matrix(res$transforms[[i]]$op)
    expect_equal(dim(op), c(d, d))
    expect_lt(norm(t(op) %*% op - diag(d), type = "F"), 1e-6)
  }

  # Embeddings should be in a shared space (up to synchronization gauge).
  Z1 <- res$embeddings[[1]]
  Z2 <- res$embeddings[[2]]
  Z3 <- res$embeddings[[3]]
  expect_true(is.matrix(Z1) && is.matrix(Z2) && is.matrix(Z3))
  expect_equal(dim(Z1), c(n, d))
  expect_equal(dim(Z2), c(n, d))
  expect_equal(dim(Z3), c(n, d))

  expect_lt(mean((Z1 - Z2)^2), 0.5)
  expect_lt(mean((Z1 - Z3)^2), 0.5)
})
