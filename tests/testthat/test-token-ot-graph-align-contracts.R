library(testthat)
library(Matrix)
library(manifoldalign)

make_graph_hyperdesign <- function(X1, X2, anchors1 = NULL, anchors2 = NULL) {
  n1 <- nrow(X1); n2 <- nrow(X2)
  d1 <- data.frame(node_id = seq_len(n1))
  d2 <- data.frame(node_id = seq_len(n2))
  if (!is.null(anchors1)) d1$anchors <- anchors1
  if (!is.null(anchors2)) d2$anchors <- anchors2
  hd <- list(domain1 = list(x = X1, design = d1),
             domain2 = list(x = X2, design = d2))
  class(hd) <- "hyperdesign"
  hd
}

test_that("token_ot_graph_align returns a valid sparse coupling with near-uniform marginals", {
  set.seed(20260224)
  n <- 24
  d <- 4
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  map12 <- match(seq_len(n), perm21)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.01), n, d)

  # Add a few anchors to reduce ambiguity
  anchor_idx <- seq(1, n, by = 6)
  a1 <- rep(NA_integer_, n)
  a2 <- rep(NA_integer_, n)
  for (k in seq_along(anchor_idx)) {
    i <- anchor_idx[[k]]
    j <- map12[[i]]
    a1[[i]] <- k
    a2[[j]] <- k
  }

  hd <- make_graph_hyperdesign(X1, X2, anchors1 = a1, anchors2 = a2)

  ctrl <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 8L,
    token_mode = "view_only",
    eps_token = 0.05,
    eps_node = 0.05,
    iters_node = 500L,
    projection = "greedy",
    anchor_weight = 2,
    anchor_penalty = 5,
    verbose = FALSE
  )

  fit <- token_ot_graph_align(hd, anchors = anchors, ncomp = 4, control = ctrl)

  expect_s3_class(fit, "token_ot_graph_align")
  expect_s3_class(fit, "multiblock_biprojector")
  expect_true(inherits(fit$transport_plan, "Matrix"))

  W <- fit$transport_plan
  expect_equal(dim(W), c(n, n))
  expect_true(all(is.finite(W@x)))
  expect_true(all(W@x >= 0))

  rsum <- as.numeric(rowSums(W))
  csum <- as.numeric(colSums(W))
  expect_equal(rsum, rep(1 / n, n), tolerance = 5e-2)
  expect_equal(csum, rep(1 / n, n), tolerance = 5e-2)
})

test_that("token_ot_graph_align coupling is equivariant to target relabeling (view_only)", {
  set.seed(77)
  n <- 20
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.01), n, d)

  hd1 <- make_graph_hyperdesign(X1, X2)
  ctrl <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 10L,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )
  fit1 <- token_ot_graph_align(hd1, ncomp = 3, control = ctrl)
  W1 <- as.matrix(fit1$transport_plan)

  q <- sample(n)
  X2p <- X2[q, , drop = FALSE]
  hd2 <- make_graph_hyperdesign(X1, X2p)
  fit2 <- token_ot_graph_align(hd2, ncomp = 3, control = ctrl)
  W2 <- as.matrix(fit2$transport_plan)

  expect_equal(W2, W1[, q, drop = FALSE], tolerance = 1e-6)
})

test_that("token_ot_graph_aligner produces a sparse perm transform usable by apply_transform()", {
  set.seed(123)
  n <- 18
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.02), n, d)

  algo <- token_ot_graph_aligner()
  ctrl <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 8L,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )
  pair_fit <- fit_pair(algo, X1, X2, ncomp = 3, control = ctrl)
  tr <- relative_transform(pair_fit, from = "i", to = "j")

  expect_true(inherits(tr$op, "Matrix"))

  Xhat <- apply_transform(tr, X1)
  expect_equal(dim(Xhat), c(n, d))
  expect_true(all(is.finite(Xhat)))
})

test_that("token_ot_graph_align runs in multilevel mode and returns hierarchy metadata", {
  set.seed(999)
  n <- 30
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.01), n, d)
  hd <- make_graph_hyperdesign(X1, X2)

  ctrl <- token_ot_graph_align_control(
    n_levels = 2L,
    prior_mode = "hard",
    views = "raw",
    candidate_k = 6L,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )

  fit <- token_ot_graph_align(hd, ncomp = 3, control = ctrl)
  expect_true(!is.null(fit$multilevel))
  expect_true(is.list(fit$multilevel$fits))
  expect_true(fit$multilevel$n_levels >= 1)
})

test_that("token_ot_graph_align supports Louvain coarsening in multilevel mode", {
  set.seed(1001)
  n <- 35
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  perm21 <- sample(n)
  X2 <- X1[perm21, , drop = FALSE] + matrix(rnorm(n * d, sd = 0.01), n, d)
  hd <- make_graph_hyperdesign(X1, X2)

  ctrl <- token_ot_graph_align_control(
    n_levels = 2L,
    coarsen_method = "louvain",
    prior_mode = "hard",
    views = "raw",
    candidate_k = 6L,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )

  fit <- token_ot_graph_align(hd, ncomp = 3, control = ctrl)
  expect_true(!is.null(fit$multilevel))
  expect_true(fit$multilevel$n_levels >= 1)
})

test_that("sparse Sinkhorn matches dense Sinkhorn when fully connected (oracle test)", {
  set.seed(20260224)
  n <- 8L
  m <- 9L
  C <- matrix(rexp(n * m, rate = 1), n, m)
  i_idx <- rep(seq_len(n), each = m)
  j_idx <- rep(seq_len(m), times = n)
  cost <- as.numeric(C[cbind(i_idx, j_idx)])

  eps <- 0.07
  iters <- 800L
  tol <- 1e-12

  sol <- manifoldalign:::.sparse_sinkhorn_items(
    i_idx = i_idx,
    j_idx = j_idx,
    cost = cost,
    n = n,
    m = m,
    eps = eps,
    iters = iters,
    tol = tol
  )

  W_sparse <- matrix(0, n, m)
  W_sparse[cbind(i_idx, j_idx)] <- sol$w

  C_shift <- C
  if (is.finite(sol$cost_shift) && sol$cost_shift != 0) C_shift <- C - sol$cost_shift
  C_scaled <- C_shift / sol$cost_scale

  a <- rep(1 / n, n)
  b <- rep(1 / m, m)
  W_dense <- sinkhorn_unified(
    C_scaled,
    a,
    b,
    epsilon = eps,
    max_iter = iters,
    tol = tol,
    stabilized = TRUE
  )

  expect_equal(W_sparse, W_dense, tolerance = 2e-4)
  expect_equal(rowSums(W_sparse), a, tolerance = 2e-4)
  expect_equal(colSums(W_sparse), b, tolerance = 2e-4)
})

test_that("ensure_cols prevents zero-incoming targets in sparse Sinkhorn", {
  set.seed(99)
  n <- 12
  d <- 3
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- matrix(0, n, d) # identical targets -> many missing columns when k=1
  hd <- make_graph_hyperdesign(X1, X2)

  ctrl_bad <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 1L,
    ensure_cols = FALSE,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )
  expect_error(
    token_ot_graph_align(hd, ncomp = 3, control = ctrl_bad),
    "zero incoming candidates"
  )

  ctrl_good <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 1L,
    ensure_cols = TRUE,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )
  fit <- token_ot_graph_align(hd, ncomp = 3, control = ctrl_good)
  expect_true(inherits(fit$transport_plan, "Matrix"))
  expect_equal(dim(fit$transport_plan), c(n, n))
})

test_that("token_ot_graph_align supports rectangular problems", {
  set.seed(101)
  n1 <- 18
  n2 <- 24
  d <- 4
  X1 <- matrix(rnorm(n1 * d), n1, d)
  X2 <- matrix(rnorm(n2 * d), n2, d)
  hd <- make_graph_hyperdesign(X1, X2)

  ctrl <- token_ot_graph_align_control(
    views = "raw",
    candidate_k = 10L,
    token_mode = "view_only",
    projection = "greedy",
    verbose = FALSE
  )
  fit <- token_ot_graph_align(hd, ncomp = 4, control = ctrl)
  W <- fit$transport_plan
  expect_equal(dim(W), c(n1, n2))
  expect_equal(as.numeric(rowSums(W)), rep(1 / n1, n1), tolerance = 5e-2)
  expect_equal(as.numeric(colSums(W)), rep(1 / n2, n2), tolerance = 5e-2)
})

test_that("benchmark_graph_alignment_methods runs (token_ot_graph only)", {
  set.seed(1)
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    n_reps = 1L,
    n_anchors = 3L,
    methods = "token_ot_graph",
    candidate_k = 6L,
    views = "raw",
    token_mode = "view_only",
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("method", "n", "rep", "runtime_sec", "top1_accuracy", "coverage") %in% names(res)))
  expect_true(all(res$method == "token_ot_graph_align"))
})

test_that("benchmark_graph_alignment_methods runs (fpgw only)", {
  set.seed(2)
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0,
    structure = "random",
    n_reps = 1L,
    methods = "fpgw",
    fpgw_omega1 = 1,
    fpgw_epsilon = 0.01,
    fpgw_max_iter = 60L,
    fpgw_inner_max_iter = 20L,
    fpgw_tol = 1e-7,
    verbose = FALSE
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(res$method == "fpgw"))
  expect_true(is.na(res$error))
  expect_equal(res$top1_accuracy, 1)
  expect_equal(res$coverage, 1)
})

test_that("benchmark_graph_alignment_methods runs (ssma only)", {
  skip_if_not_installed("adjoin")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(3)
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.02,
    structure = "ring",
    n_reps = 1L,
    n_anchors = 3L,
    methods = "ssma",
    ssma_solver = "reduced",
    ssma_knn = 6L,
    ssma_rank_per_domain = 8L,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(res$method == "ssma_align"))
  expect_true(all(is.na(res$error)))
  expect_true(all(is.finite(res$top1_accuracy)))
  expect_true(all(is.finite(res$coverage)))
})

test_that("benchmark SSMA procrustes decode is opt-in", {
  skip_if_not_installed("adjoin")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("clue")

  args <- list(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.02,
    structure = "ring",
    n_reps = 1L,
    n_anchors = 4L,
    methods = "ssma",
    ssma_solver = "reduced",
    ssma_knn = 6L,
    ssma_rank_per_domain = 8L,
    seed = 77L,
    verbose = FALSE
  )

  res_nn <- do.call(
    benchmark_graph_alignment_methods,
    c(args, list(decode_mode = "common_nn", ssma_procrustes = FALSE))
  )
  res_lap_default <- do.call(
    benchmark_graph_alignment_methods,
    c(args, list(decode_mode = "common_procrustes_lap", ssma_procrustes = FALSE))
  )
  res_lap_optin <- do.call(
    benchmark_graph_alignment_methods,
    c(args, list(decode_mode = "common_procrustes_lap", ssma_procrustes = TRUE))
  )

  expect_equal(res_lap_default$top1_accuracy, res_nn$top1_accuracy)
  expect_equal(res_lap_default$coverage, res_nn$coverage)
  expect_true(all(is.finite(res_lap_optin$top1_accuracy)))
  expect_true(all(is.finite(res_lap_optin$coverage)))
})

test_that("benchmark_graph_alignment_methods runs (lra only)", {
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("RSpectra")

  set.seed(4)
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.02,
    structure = "ring",
    n_reps = 1L,
    n_anchors = 3L,
    methods = "lra",
    lra_solver = "operator",
    lra_lambda = 0.01,
    lra_mu = 0.5,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(res$method == "lowrank_align"))
  expect_true(all(is.na(res$error)))
  expect_true(all(is.finite(res$top1_accuracy)))
  expect_true(all(is.finite(res$coverage)))
})

test_that("benchmark_graph_alignment_methods runs (gpca only)", {
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("genpca")
  skip_if_not_installed("PRIMME")

  set.seed(5)
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.02,
    structure = "ring",
    n_reps = 1L,
    n_anchors = 3L,
    methods = "gpca",
    gpca_u = 0.4,
    gpca_lambda = 0.1,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(res$method == "gpca_align"))
  expect_true(all(is.na(res$error)))
  expect_true(all(is.finite(res$top1_accuracy)))
  expect_true(all(is.finite(res$coverage)))
})
