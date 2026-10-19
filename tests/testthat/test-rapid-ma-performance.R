# Exactness of the RAPID-MA performance paths. Each fast path is compared
# with the implementation it replaced (kept here as a reference) or with an
# independent brute-force oracle, and the default fit is pinned to values
# produced by the code before these changes (commit 2bf46b7).

ns <- asNamespace("manifoldalign")

# Pre-change reference implementations -----------------------------------

reference_force_dgC <- function(A) {
  if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  if (methods::is(A, "nMatrix")) A <- methods::as(A, "dMatrix")
  if (methods::is(A, "symmetricMatrix") ||
      methods::is(A, "triangularMatrix") ||
      methods::is(A, "diagonalMatrix")) {
    A <- methods::as(A, "generalMatrix")
  }
  A <- Matrix::drop0(Matrix::Matrix(A, sparse = TRUE))
  sm <- Matrix::summary(A)
  if (!nrow(sm)) {
    return(ns$.rapid_empty_relation(nrow(A), dimnames = dimnames(A)))
  }
  Matrix::sparseMatrix(
    i = as.integer(sm$i), j = as.integer(sm$j), x = as.numeric(sm$x),
    dims = dim(A), dimnames = dimnames(A)
  )
}

reference_mutual <- function(A) {
  if (!length(A@x)) return(A)
  At <- Matrix::t(A)
  mask <- (A != 0) * (At != 0)
  reference_force_dgC(Matrix::drop0(((A + At) / 2) * mask))
}

reference_scale_dense <- function(Z, drop_constant = FALSE) {
  Z <- as.matrix(Z)
  center <- colMeans(Z)
  center[!is.finite(center)] <- 0
  Zc <- sweep(Z, 2L, center, "-")
  denom <- max(nrow(Zc) - 1L, 1L)
  scale <- sqrt(colSums(Zc * Zc) / denom)
  active <- is.finite(scale) & scale > 1e-12
  scale[!active] <- 1
  out <- sweep(Zc, 2L, scale, "/")
  out[!is.finite(out)] <- 0
  if (isTRUE(drop_constant)) out <- out[, active, drop = FALSE]
  list(view = out, center = center, scale = scale, active = active)
}

# Canonical exact kNN: the `k` nearest other rows ordered by (distance,
# lexicographic row rank). Exact for integer-valued views.
oracle_knn <- function(X, k) {
  rank <- ns$.rapid_view_tie_rank(X)
  n <- nrow(X)
  index <- distance <- vector("list", n)
  for (i in seq_len(n)) {
    d2 <- colSums((t(X) - X[i, ])^2)
    d2[i] <- Inf
    o <- order(d2, rank)[seq_len(k)]
    index[[i]] <- o
    distance[[i]] <- sqrt(d2[o])
  }
  list(candidate = unlist(index), distance = unlist(distance))
}

random_relation <- function(n, density, seed, dimnames = NULL) {
  ns$.with_local_seed(seed, {
    A <- Matrix::rsparsematrix(n, n, density, rand.x = function(m) runif(m))
  })
  A <- methods::as(methods::as(A, "generalMatrix"), "CsparseMatrix")
  if (!is.null(dimnames)) dimnames(A) <- dimnames
  A
}

# kNN backends -------------------------------------------------------------

test_that("blocked exact kNN matches the kd-tree on continuous views", {
  for (d in c(5L, 10L, 24L)) {
    X <- ns$.with_local_seed(100L + d, matrix(stats::rnorm(700L * d), 700L, d))
    X[1:350, 1:3] <- X[1:350, 1:3] + 4
    for (k in c(1L, 7L, 15L)) {
      blocked <- ns$.rapid_self_knn(X, k, "blocked_exact")
      tree <- ns$.rapid_self_knn(X, k, "rann")
      expect_identical(blocked$query, tree$query)
      expect_identical(blocked$candidate, tree$candidate)
      expect_equal(blocked$distance, tree$distance, tolerance = 1e-14)
    }
  }
})

test_that("blocked exact kNN equals a brute-force oracle with exact ties", {
  # Integer grid coordinates and duplicated rows produce many exactly tied
  # distances, including at the k-th neighbour.
  grid <- as.matrix(expand.grid(a = 0:5, b = 0:4, c = 0:1, d = 0:1, e = 0:1))
  X <- rbind(grid, grid[seq(1, nrow(grid), by = 7L), ])
  storage.mode(X) <- "double"
  oracle <- oracle_knn(X, 6L)
  for (leaf in c(8L, 16L, 64L, 256L)) {
    found <- ns$rapid_exact_knn_cpp(X, ns$.rapid_view_tie_rank(X), 6L, leaf)
    expect_identical(found$index, oracle$candidate)
    expect_identical(found$distance, oracle$distance)
  }
})

test_that("blocked exact kNN is equivariant to row permutations", {
  # Distinct rows with many tied distances: lexicographic ranks of distinct
  # rows do not depend on row order, so the answer is equivariant.
  X <- as.matrix(expand.grid(0:4, 0:4, 0:2, 0:1, 0:1))
  storage.mode(X) <- "double"
  perm <- ns$.with_local_seed(5L, sample.int(nrow(X)))
  base <- ns$.rapid_self_knn(X, 9L, "blocked_exact")
  moved <- ns$.rapid_self_knn(X[perm, , drop = FALSE], 9L, "blocked_exact")
  # Row r of the permuted view is row perm[r] of X.
  mapped <- matrix(perm[moved$candidate], ncol = 9L, byrow = TRUE)
  expected <- matrix(base$candidate, ncol = 9L, byrow = TRUE)[perm, ]
  expect_identical(mapped, expected)
  expect_identical(
    matrix(moved$distance, ncol = 9L, byrow = TRUE),
    matrix(base$distance, ncol = 9L, byrow = TRUE)[perm, ]
  )
})

test_that("blocked exact kNN rejects invalid input", {
  X <- matrix(stats::rnorm(40), 8L, 5L)
  rank <- ns$.rapid_view_tie_rank(X)
  expect_error(ns$rapid_exact_knn_cpp(X, rank, 8L), "1..nrow")
  expect_error(ns$rapid_exact_knn_cpp(X, rank[-1L], 2L), "tie_rank")
  X[2L, 3L] <- NaN
  expect_error(ns$rapid_exact_knn_cpp(X, rank, 2L), "finite")
})

test_that("neighbour relations choose the backend by view width", {
  expect_identical(ns$.rapid_knn_backend(2L), "rann")
  expect_identical(ns$.rapid_knn_backend(4L), "rann")
  expect_identical(ns$.rapid_knn_backend(5L), "blocked_exact")
  X <- ns$.with_local_seed(3L, matrix(stats::rnorm(300L * 8L), 300L, 8L))
  built <- ns$.rapid_build_neighbor_relation(X, k = 6L, degree_cap = 16L)
  expect_identical(built$metadata$neighbor_backend, "blocked_exact")
  testthat::local_mocked_bindings(
    .rapid_knn_backend = function(n_columns) "rann",
    .package = "manifoldalign"
  )
  tree <- ns$.rapid_build_neighbor_relation(X, k = 6L, degree_cap = 16L)
  expect_identical(tree$metadata$neighbor_backend, "rann")
  expect_identical(built$A, tree$A)
  radius <- ns$.rapid_build_neighbor_relation(
    X, k = 6L, degree_cap = 16L, radius = 3
  )
  testthat::local_mocked_bindings(
    .rapid_knn_backend = function(n_columns) "blocked_exact",
    .package = "manifoldalign"
  )
  expect_identical(
    ns$.rapid_build_neighbor_relation(X, k = 6L, degree_cap = 16L, radius = 3)$A,
    radius$A
  )
})

# Sparse relation helpers ----------------------------------------------------

test_that("mutual relations and the dgC fast path match the matrix algebra", {
  names <- paste0("node", 1:40)
  cases <- list(
    random_relation(40L, 0.08, 1L),
    random_relation(40L, 0.3, 2L, dimnames = list(names, names)),
    random_relation(1L, 1, 3L)
  )
  with_zeros <- random_relation(30L, 0.2, 4L)
  with_zeros@x[c(1L, 5L, 9L)] <- 0
  cases[[length(cases) + 1L]] <- with_zeros
  symmetric <- random_relation(25L, 0.2, 6L)
  cases[[length(cases) + 1L]] <- methods::as(
    methods::as(symmetric + Matrix::t(symmetric), "generalMatrix"),
    "CsparseMatrix"
  )
  for (A in cases) {
    expect_identical(ns$.rapid_force_dgC(A), reference_force_dgC(A))
    expect_identical(ns$.rapid_mutual_relation(A), reference_mutual(A))
  }
  empty <- ns$.rapid_empty_relation(5L)
  expect_identical(ns$.rapid_mutual_relation(empty), empty)
})

test_that("column scaling matches the sweep implementation", {
  Z <- ns$.with_local_seed(8L, matrix(stats::rnorm(600), 100L, 6L))
  Z[, 3L] <- 2
  colnames(Z) <- letters[1:6]
  for (drop in c(FALSE, TRUE)) {
    expect_identical(
      ns$.rapid_scale_dense(Z, drop_constant = drop),
      reference_scale_dense(Z, drop_constant = drop)
    )
  }
  expect_identical(
    ns$.rapid_scale_dense(Z[0, , drop = FALSE])$view,
    reference_scale_dense(Z[0, , drop = FALSE])$view
  )
})

test_that("prototype relation uses a dense intermediate without changing values", {
  A <- random_relation(300L, 0.03, 9L)
  Q <- random_relation(300L, 0.02, 10L)[, 1:6]
  Q <- Q + Matrix::sparseMatrix(i = 1:300, j = rep(1:6, 50), x = 0.5,
                                dims = c(300L, 6L))
  Q <- methods::as(Q, "CsparseMatrix")
  fit <- ns$.rapid_prototype_relation(A, Q)
  raw <- reference_force_dgC(Matrix::crossprod(Q, A %*% Q))
  expect_identical(fit$raw, raw)
  expect_identical(fit$diagnostics$intermediate_nnz,
                   as.integer(Matrix::nnzero(A %*% Q)))
})

test_that("sparse UOT keeps every finite candidate when q covers the bank", {
  fixture <- rapid_ma_benchmark_fixture(
    n_common = 80L, extra_target = 10L, labels_per_class = 3L, seed = 2L
  )
  embedding <- ns$.with_local_seed(2L, matrix(stats::rnorm(80L * 4L), 80L, 4L))
  prototypes <- ns$.rapid_initialize_prototypes(
    list(embedding), list(fixture$labels[[1L]])
  )
  control <- ns$.rapid_uot_control(q = nrow(prototypes$embedding),
                                   hard_labels = TRUE, block_size = 32L)
  sparse <- ns$.rapid_uot_build_cost(
    embedding, prototypes, labels = fixture$labels[[1L]], control = control
  )
  dense <- ns$.rapid_uot_build_cost(
    embedding, prototypes, labels = fixture$labels[[1L]], control = control,
    mode = "dense"
  )
  expect_identical(sparse$prototype_index, dense$prototype_index)
  rows <- rep(seq_len(80L), diff(sparse$cost$row_ptr))
  cost <- dense$cost
  finite <- which(is.finite(cost), arr.ind = TRUE)
  finite <- finite[order(finite[, 1L], finite[, 2L]), , drop = FALSE]
  expect_identical(rows, as.integer(finite[, 1L]))
  expect_identical(sparse$cost$col_idx, as.integer(finite[, 2L]))
  expect_identical(sparse$cost$cost, cost[finite])
  # Hard labels make some candidates infinite, so rows differ in support.
  expect_gt(length(unique(diff(sparse$cost$row_ptr))), 1L)
})

# End-to-end pins -----------------------------------------------------------

fit_small_fixture <- function() {
  fixture <- rapid_ma_benchmark_fixture(
    n_common = 150L, extra_target = 30L, labels_per_class = 3L, seed = 11L
  )
  rapid_ma(
    fixture$data, labels = fixture$labels, positions = fixture$positions,
    attributes = fixture$attributes, ncomp = 6L,
    control = list(max_iter = 3L, min_iter = 1L, seed = 11L)
  )
}

test_that("default fit reproduces the values of the pre-optimization code", {
  fit <- fit_small_fixture()
  # Values printed with 17 significant digits from commit 2bf46b7.
  expect_equal(
    fit$objective_history$objective,
    c(0.8635941775972994, 0.51176922165744432, 0.41721012477992481,
      0.37049754394243123),
    tolerance = 1e-10
  )
  expect_equal(
    unname(vapply(fit$scores, sum, numeric(1))),
    c(3.6900118212374062, -10.075350528046529),
    tolerance = 1e-9
  )
  expect_equal(
    unname(vapply(fit$scores, function(x) sum(x^2), numeric(1))),
    c(140.39152647460475, 56.68455069309843),
    tolerance = 1e-10
  )
  expect_equal(
    unname(vapply(fit$couplings, function(x) sum(x@x * seq_along(x@x)), numeric(1))),
    c(452.56642261443341, 506.71083664457632),
    tolerance = 1e-10
  )
  relation_nnz <- unlist(lapply(fit$relations, function(state) {
    vapply(state$relations, Matrix::nnzero, numeric(1))
  }), use.names = FALSE)
  expect_identical(relation_nnz, c(1642, 846, 1856, 1846, 1030, 2246))
  relation_energy <- unlist(lapply(fit$relations, function(state) {
    vapply(state$relations, function(A) sum(A@x^2), numeric(1))
  }), use.names = FALSE)
  expect_equal(
    relation_energy,
    c(18.680903486550175, 36.42728468122386, 20.544879871472361,
      23.593472189789498, 42.557700239051599, 24.391160463586434),
    tolerance = 1e-10
  )
  expect_identical(
    unname(lapply(fit$predictions, function(p) as.vector(table(p)))),
    list(c(34L, 40L, 44L, 32L), c(28L, 53L, 54L, 45L))
  )
  backends <- unlist(lapply(fit$relations, function(state) {
    vapply(state$metadata, function(m) m$neighbor_backend, character(1))
  }), use.names = FALSE)
  expect_identical(
    backends,
    c("blocked_exact", "rann", "rann", "blocked_exact", "rann", "rann")
  )
})

test_that("default fit is identical under the kd-tree backend on continuous data", {
  fast <- fit_small_fixture()
  testthat::local_mocked_bindings(
    .rapid_knn_backend = function(n_columns) "rann",
    .package = "manifoldalign"
  )
  tree <- fit_small_fixture()
  for (field in c("scores", "couplings", "prediction_probabilities",
                  "predictions", "objective_history", "base_scores")) {
    expect_identical(fast[[field]], tree[[field]], info = field)
  }
  expect_identical(
    lapply(fast$relations, `[[`, "relations"),
    lapply(tree$relations, `[[`, "relations")
  )
})

test_that("blocked kNN relations keep linear retained state", {
  sizes <- c(400L, 800L)
  built <- lapply(sizes, function(n) {
    X <- ns$.with_local_seed(n, matrix(stats::rnorm(n * 12L), n, 12L))
    ns$.rapid_build_neighbor_relation(X, k = 10L, degree_cap = 10L,
                                      symmetrize = "none")
  })
  for (i in seq_along(sizes)) {
    expect_identical(built[[i]]$metadata$neighbor_backend, "blocked_exact")
    expect_equal(Matrix::nnzero(built[[i]]$A), sizes[[i]] * 10)
    expect_lte(max(diff(built[[i]]$A@p)), sizes[[i]])
    expect_identical(as.integer(table(Matrix::summary(built[[i]]$A)$i)),
                     rep(10L, sizes[[i]]))
  }
  retained <- vapply(built, function(x) as.numeric(utils::object.size(x)),
                     numeric(1))
  expect_lt(retained[[2L]] / retained[[1L]], 2.5)
})
