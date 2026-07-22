rapid_prepare_relations <- function(...) {
  manifoldalign:::.rapid_prepare_relations(...)
}

expect_sparse_relation <- function(A, n, degree_cap) {
  expect_true(inherits(A, "dgCMatrix"))
  expect_identical(dim(A), c(n, n))
  expect_true(all(is.finite(A@x)))
  expect_true(all(A@x >= 0))
  counts <- tabulate(Matrix::summary(A)$i, nbins = n)
  expect_lte(max(counts), degree_cap)
  expect_lte(Matrix::nnzero(A), n * degree_cap)
}

test_that("generated and custom relations stay sparse and degree capped", {
  set.seed(101)
  n <- 48L
  X <- matrix(rnorm(n * 80L), nrow = n)
  topology <- Matrix::sparseMatrix(
    i = c(rep(1L, n - 1L), 2:n),
    j = c(2:n, rep(1L, n - 1L)),
    x = 1,
    dims = c(n, n)
  )

  prepared <- rapid_prepare_relations(
    X,
    relations = list(topology = topology),
    build_feature = TRUE,
    feature_k = 12L,
    feature_sketch_dim = 16L,
    degree_cap = 5L
  )

  expect_setequal(names(prepared$relations), c("topology", "feature"))
  lapply(prepared$relations, expect_sparse_relation, n = n, degree_cap = 5L)
  expect_true(all(prepared$diagnostics$summary$max_out_degree <= 5L))
  expect_true(prepared$metadata$feature$view$projected)
  expect_identical(prepared$metadata$feature$view$output_dim, 16L)
  expect_identical(ncol(prepared$views$feature$view), 16L)
  expect_identical(prepared$metadata$topology$symmetrize, "mutual")
  expect_gt(prepared$metadata$topology$dropped_for_degree_cap, 0L)
})

test_that("relation construction is equivariant to node permutations", {
  set.seed(202)
  n <- 36L
  X <- matrix(rnorm(n * 9L), nrow = n)
  positions <- matrix(rnorm(n * 2L), nrow = n)
  attributes <- data.frame(
    elevation = rnorm(n),
    cover = paste0("class-", seq_len(n)),
    flag = rep(c(TRUE, FALSE), length.out = n),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = seq_len(n),
    to = c(2:n, 1L),
    weight = seq(0.5, 1.5, length.out = n)
  )
  permutation <- sample.int(n)

  original <- rapid_prepare_relations(
    X,
    relations = list(cycle = edges),
    positions = positions,
    attributes = attributes,
    build_feature = TRUE,
    feature_k = 7L,
    spatial_k = 7L,
    attribute_k = 7L,
    degree_cap = 7L,
    seed = 17L
  )

  inverse <- integer(n)
  inverse[permutation] <- seq_len(n)
  permuted_edges <- transform(
    edges,
    from = inverse[from],
    to = inverse[to]
  )
  permuted <- rapid_prepare_relations(
    X[permutation, , drop = FALSE],
    relations = list(cycle = permuted_edges),
    positions = positions[permutation, , drop = FALSE],
    attributes = attributes[permutation, , drop = FALSE],
    build_feature = TRUE,
    feature_k = 7L,
    spatial_k = 7L,
    attribute_k = 7L,
    degree_cap = 7L,
    seed = 17L
  )

  expect_setequal(names(permuted$relations), names(original$relations))
  for (name in names(original$relations)) {
    expected <- original$relations[[name]][permutation, permutation]
    expect_equal(
      as.matrix(permuted$relations[[name]]),
      as.matrix(expected),
      tolerance = 1e-9,
      info = name
    )
  }
})

test_that("relative spatial graphs ignore translation, rotation, and scale", {
  set.seed(303)
  n <- 50L
  X <- matrix(rnorm(n * 3L), nrow = n)
  positions <- matrix(rnorm(n * 2L), nrow = n)
  theta <- pi / 5
  rotation <- matrix(
    c(cos(theta), -sin(theta), sin(theta), cos(theta)),
    nrow = 2L,
    byrow = TRUE
  )
  transformed <- 11 * positions %*% rotation
  transformed <- sweep(transformed, 2L, c(140, -73), "+")

  original <- rapid_prepare_relations(
    X,
    relations = list(),
    positions = positions,
    build_feature = FALSE,
    spatial_mode = "relative",
    spatial_k = 9L,
    degree_cap = 9L
  )
  changed <- rapid_prepare_relations(
    X,
    relations = list(),
    positions = transformed,
    build_feature = FALSE,
    spatial_mode = "relative",
    spatial_k = 9L,
    degree_cap = 9L
  )

  expect_equal(
    as.matrix(changed$relations$spatial),
    as.matrix(original$relations$spatial),
    tolerance = 1e-9
  )
  coordinate_meta <- original$metadata$spatial$coordinates
  expect_true(coordinate_meta$translation_invariant)
  expect_true(coordinate_meta$rotation_invariant)
  expect_true(coordinate_meta$scale_invariant)

  common <- rapid_prepare_relations(
    X,
    relations = list(),
    positions = positions,
    build_feature = FALSE,
    spatial_mode = "common",
    spatial_radius = 2,
    spatial_sigma = 1,
    degree_cap = 6L
  )
  expect_false(common$metadata$spatial$coordinates$scale_invariant)
})

test_that("mixed attributes, missing positions, and isolates stay finite", {
  set.seed(404)
  n <- 30L
  X <- matrix(rnorm(n * 4L), nrow = n)
  positions <- matrix(rnorm(n * 2L), nrow = n)
  positions[c(3L, 19L), ] <- NA_real_
  attributes <- data.frame(
    numeric = c(NA_real_, Inf, rnorm(n - 2L)),
    category = c(NA_character_, rep(c("forest", "water"), length.out = n - 1L)),
    logical = rep(c(TRUE, FALSE, NA), length.out = n),
    constant = 1,
    all_missing = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  isolate_edges <- data.frame(
    from = c(1L, 2L, 4L),
    to = c(2L, 1L, 5L),
    weight = c(1, 1, 2)
  )

  prepared <- rapid_prepare_relations(
    X,
    relations = list(isolates = isolate_edges),
    positions = positions,
    attributes = attributes,
    build_feature = FALSE,
    spatial_k = 5L,
    attribute_k = 5L,
    normalization = "random_walk",
    degree_cap = 5L
  )

  expect_true(all(is.finite(prepared$views$position$view)))
  expect_true(all(is.finite(prepared$views$attribute$view)))
  expect_lte(ncol(prepared$views$attribute$view), 3L * 32L)
  expect_identical(
    unname(prepared$views$attribute$metadata$missing_by_column),
    c(2L, 1L, 10L, 0L, n)
  )
  expect_identical(prepared$views$attribute$metadata$n_missing_rows, n)
  expect_equal(prepared$metadata$spatial$missing_rows, c(3L, 19L))
  expect_gte(prepared$metadata$isolates$n_isolates, n - 5L)

  for (relation in prepared$relations) {
    expect_true(all(is.finite(relation@x)))
    row_sums <- as.numeric(Matrix::rowSums(relation))
    expect_true(all(row_sums == 0 | abs(row_sums - 1) < 1e-10))
  }
})

test_that("radius and dense relations cannot bypass sparsity guards", {
  set.seed(505)
  n <- 42L
  X <- matrix(rnorm(n * 2L), nrow = n)
  positions <- matrix(rnorm(n * 2L, sd = 0.01), nrow = n)

  radius <- rapid_prepare_relations(
    X,
    relations = list(),
    positions = positions,
    build_feature = FALSE,
    spatial_radius = 100,
    degree_cap = 3L
  )
  expect_sparse_relation(radius$relations$spatial, n, degree_cap = 3L)

  dense <- matrix(1, nrow = 6L, ncol = 6L)
  diag(dense) <- 0
  expect_error(
    rapid_prepare_relations(
      matrix(rnorm(12), nrow = 6L),
      relations = list(dense = dense),
      dense_max_n = 5L
    ),
    "Dense relation"
  )
  small <- rapid_prepare_relations(
    matrix(rnorm(12), nrow = 6L),
    relations = list(dense = dense),
    degree_cap = 2L,
    dense_max_n = 6L
  )
  expect_sparse_relation(small$relations$dense, 6L, degree_cap = 2L)
})

test_that("large sparse custom relations are processed without densification", {
  n <- 5000L
  chain <- Matrix::sparseMatrix(
    i = c(seq_len(n - 1L), 2:n),
    j = c(2:n, seq_len(n - 1L)),
    x = 1,
    dims = c(n, n)
  )
  prepared <- rapid_prepare_relations(
    matrix(seq_len(n), ncol = 1L),
    relations = list(chain = chain),
    build_feature = FALSE,
    degree_cap = 3L,
    dense_max_n = 10L
  )

  relation <- prepared$relations$chain
  expect_true(inherits(relation, "dgCMatrix"))
  expect_identical(dim(relation), c(n, n))
  expect_lte(Matrix::nnzero(relation), 2L * (n - 1L))
  expect_true(all(is.finite(relation@x)))
  expect_lt(as.numeric(object.size(relation)), 2e6)
})

test_that("directed edge relations retain direction under random-walk normalization", {
  n <- 8L
  edges <- data.frame(
    from = c(1L, 1L, 2L, 3L, 6L),
    to = c(2L, 3L, 4L, 5L, 7L),
    weight = c(1, 3, 2, 1, 4)
  )
  prepared <- rapid_prepare_relations(
    matrix(seq_len(n), ncol = 1L),
    relations = list(directed = edges),
    build_feature = FALSE,
    normalization = "random_walk",
    degree_cap = 3L
  )

  relation <- prepared$relations$directed
  row_sums <- as.numeric(Matrix::rowSums(relation))
  expect_true(prepared$metadata$directed$directed)
  expect_equal(relation[2L, 1L], 0)
  expect_equal(row_sums[row_sums > 0], rep(1, sum(row_sums > 0)))
})
