rapid_prototype_relation <- function(...) {
  manifoldalign:::.rapid_prototype_relation(...)
}

rapid_relmatch_fit <- function(...) {
  manifoldalign:::.rapid_relmatch_fit(...)
}

test_that("sparse prototype accumulation matches dense Q transpose A Q", {
  set.seed(1401)
  n <- 17L
  K <- 5L
  A <- Matrix::rsparsematrix(n, n, density = 0.18)
  A@x <- abs(A@x)
  Matrix::diag(A) <- 0
  Q <- matrix(runif(n * K), nrow = n)
  Q[Q < 0.45] <- 0

  fit <- rapid_prototype_relation(
    A, Q, control = list(normalization = "mass", storage = "dense")
  )
  expected <- crossprod(Q, as.matrix(A) %*% Q)

  expect_equal(as.matrix(fit$raw), expected, tolerance = 1e-12)
  expect_equal(fit$graph, expected / sum(expected), tolerance = 1e-12)
  expect_identical(fit$diagnostics$status, "ok")
  expect_false(fit$diagnostics$dense_node_square_allocated)
})

test_that("node and prototype permutations preserve projected structure", {
  set.seed(1402)
  n <- 29L
  K <- 6L
  A <- Matrix::rsparsematrix(n, n, density = 0.12)
  A@x <- abs(A@x)
  Q <- matrix(runif(n * K), nrow = n)
  Q[Q < 0.6] <- 0
  node_permutation <- sample.int(n)
  prototype_permutation <- sample.int(K)

  original <- rapid_prototype_relation(A, Q, control = list(storage = "dense"))
  node_permuted <- rapid_prototype_relation(
    A[node_permutation, node_permutation],
    Q[node_permutation, , drop = FALSE],
    control = list(storage = "dense")
  )
  prototype_permuted <- rapid_prototype_relation(
    A,
    Q[, prototype_permutation, drop = FALSE],
    control = list(storage = "dense")
  )

  expect_equal(node_permuted$graph, original$graph, tolerance = 1e-12)
  expect_equal(
    prototype_permuted$graph,
    original$graph[prototype_permutation, prototype_permutation],
    tolerance = 1e-12
  )
  expect_equal(
    as.matrix(prototype_permuted$raw),
    as.matrix(original$raw)[prototype_permutation, prototype_permutation],
    tolerance = 1e-12
  )
})

test_that("isomorphic domains have zero loss and corruption increases it", {
  n <- 40L
  K <- 4L
  group <- rep(seq_len(K), each = n / K)
  next_group <- (group %% K) + 1L
  target <- vapply(seq_len(n), function(i) {
    which(group == next_group[[i]])[[1L]]
  }, integer(1))
  A <- Matrix::sparseMatrix(
    i = seq_len(n), j = target, x = 1, dims = c(n, n)
  )
  Q <- Matrix::sparseMatrix(
    i = seq_len(n), j = group, x = 1, dims = c(n, K)
  )
  permutation <- rev(seq_len(n))
  isomorphic <- rapid_relmatch_fit(
    list(first = list(topology = A),
         second = list(topology = A[permutation, permutation])),
    list(Q, Q[permutation, , drop = FALSE]),
    control = list(storage = "dense")
  )

  corrupted <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = rep(which(group == 1L)[[1L]], n),
    x = 1,
    dims = c(n, n)
  )
  damaged <- rapid_relmatch_fit(
    list(first = list(topology = A), second = list(topology = corrupted)),
    list(Q, Q),
    control = list(storage = "dense")
  )

  expect_lt(isomorphic$total_loss, 1e-12)
  expect_gt(damaged$total_loss, isomorphic$total_loss + 0.1)
  expect_equal(
    isomorphic$barycenters$topology,
    isomorphic$projected$first$topology$graph,
    tolerance = 1e-12
  )
})

test_that("robust relation gates suppress cross-domain corruption", {
  n <- 48L
  K <- 4L
  group <- rep(seq_len(K), each = n / K)
  Q <- Matrix::sparseMatrix(
    i = seq_len(n), j = group, x = 1, dims = c(n, K)
  )
  good <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = c(2:n, 1L),
    x = 1,
    dims = c(n, n)
  )
  bad <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = rep(1L, n),
    x = 1,
    dims = c(n, n)
  )
  fit <- rapid_relmatch_fit(
    list(
      d1 = list(good = good, unreliable = good),
      d2 = list(good = good, unreliable = bad)
    ),
    list(Q, Q),
    control = list(
      relation_gate = "robust",
      gate_temperature = 10,
      storage = "dense"
    )
  )

  expect_gt(fit$relation_weights[["good"]], 0.9)
  expect_lt(fit$relation_weights[["unreliable"]], 0.1)
  expect_lt(
    max(fit$losses[, "good"]),
    min(fit$losses[, "unreliable"])
  )
})

test_that("empty relations, missing relations, and unmatched mass are explicit", {
  n <- 12L
  K <- 3L
  empty <- Matrix::Matrix(0, n, n, sparse = TRUE)
  cycle <- Matrix::sparseMatrix(
    i = seq_len(n), j = c(2:n, 1L), x = 1, dims = c(n, n)
  )
  Q <- Matrix::sparseMatrix(
    i = seq_len(n), j = rep(seq_len(K), length.out = n), x = 1,
    dims = c(n, K)
  )
  zero_Q <- Matrix::Matrix(0, n, K, sparse = TRUE)

  empty_fit <- rapid_prototype_relation(empty, Q)
  unmatched_fit <- rapid_prototype_relation(cycle, zero_Q)
  combined <- rapid_relmatch_fit(
    list(d1 = list(topology = cycle, empty = empty),
         d2 = list(topology = cycle)),
    list(Q, Q)
  )

  expect_identical(empty_fit$diagnostics$status, "empty_relation")
  expect_identical(unmatched_fit$diagnostics$status, "unmatched_mass")
  expect_identical(
    combined$diagnostics$relation_status$d1[["empty"]],
    "empty_relation"
  )
  expect_identical(
    combined$diagnostics$relation_status$d2[["empty"]],
    "missing_relation"
  )
  expect_equal(combined$relation_weights[["empty"]], 0)
  expect_identical(dim(combined$barycenters$empty), c(K, K))
})

test_that("fixed relation weights are respected", {
  n <- 16L
  Q <- Matrix::Diagonal(n = n)
  A <- Matrix::sparseMatrix(
    i = seq_len(n), j = c(2:n, 1L), x = 1, dims = c(n, n)
  )
  fit <- rapid_relmatch_fit(
    list(d1 = list(a = A, b = A), d2 = list(a = A, b = A)),
    list(Q, Q),
    control = list(
      relation_gate = "fixed",
      relation_weights = c(a = 3, b = 1)
    )
  )

  expect_equal(unname(fit$relation_weights[c("a", "b")]), c(0.75, 0.25))
})

test_that("storage policy only densifies prototype-sized graphs", {
  n <- 20L
  K <- 5L
  Q <- Matrix::sparseMatrix(
    i = seq_len(n), j = rep(seq_len(K), length.out = n), x = 1,
    dims = c(n, K)
  )
  sparse_A <- Matrix::sparseMatrix(
    i = seq_len(n), j = c(2:n, 1L), x = 1, dims = c(n, n)
  )
  dense_A <- Matrix::Matrix(matrix(1, n, n) - diag(n), sparse = TRUE)

  sparse_fit <- rapid_prototype_relation(
    sparse_A, Q, control = list(storage = "sparse")
  )
  dense_fit <- rapid_prototype_relation(
    dense_A, Q, control = list(storage = "dense")
  )

  expect_s4_class(sparse_fit$graph, "dgCMatrix")
  expect_true(is.matrix(dense_fit$graph))
  expect_identical(dim(dense_fit$graph), c(K, K))
  expect_false(dense_fit$diagnostics$dense_node_square_allocated)
})

test_that("large sparse accumulation retains no node-square intermediate", {
  n <- 5000L
  K <- 20L
  q <- 2L
  A <- Matrix::sparseMatrix(
    i = seq_len(n), j = c(2:n, 1L), x = 1, dims = c(n, n)
  )
  Q <- Matrix::sparseMatrix(
    i = rep(seq_len(n), each = q),
    j = cbind(
      ((seq_len(n) - 1L) %% K) + 1L,
      (seq_len(n) %% K) + 1L
    ),
    x = 0.5,
    dims = c(n, K)
  )
  fit <- rapid_prototype_relation(A, Q, control = list(storage = "sparse"))

  expect_identical(dim(fit$graph), c(K, K))
  expect_identical(fit$diagnostics$relation_nnz, n)
  expect_identical(fit$diagnostics$q_max, q)
  expect_equal(fit$diagnostics$operation_bound, n * q^2)
  expect_lte(fit$diagnostics$intermediate_nnz, n * K)
  expect_false(fit$diagnostics$dense_node_square_allocated)
  expect_lt(as.numeric(object.size(fit)), n * n * 8)
})
