rapid_prepare_diffusion_relations <- function(...) {
  manifoldalign:::.rapid_prepare_relations(...)
}

rapid_diffusion_encode <- function(...) {
  manifoldalign:::.rapid_diffusion_encode(...)
}

directed_cycle <- function(n, offset = 1L) {
  Matrix::sparseMatrix(
    i = seq_len(n),
    j = ((seq_len(n) - 1L + offset) %% n) + 1L,
    x = 1,
    dims = c(n, n)
  )
}

test_that("diffusion is equivariant to node permutations", {
  set.seed(1101)
  n <- 47L
  X <- matrix(rnorm(n * 9L), nrow = n)
  B <- cbind(X[, 1L:4L], rowMeans(X[, 5L:9L, drop = FALSE]))
  topology <- directed_cycle(n)
  permutation <- sample.int(n)

  original_relations <- rapid_prepare_diffusion_relations(
    X,
    relations = list(topology = topology),
    build_feature = FALSE,
    normalization = "random_walk",
    custom_symmetrize = "none",
    degree_cap = 3L
  )
  permuted_relations <- rapid_prepare_diffusion_relations(
    X[permutation, , drop = FALSE],
    relations = list(topology = topology[permutation, permutation]),
    build_feature = FALSE,
    normalization = "random_walk",
    custom_symmetrize = "none",
    degree_cap = 3L
  )
  control <- list(
    steps = c(0L, 1L, 2L, 4L),
    output_dim = 7L,
    gate = "uniform",
    seed = 19L
  )

  original <- rapid_diffusion_encode(
    X, original_relations, control = control, precomputed = B
  )
  permuted <- rapid_diffusion_encode(
    X[permutation, , drop = FALSE],
    permuted_relations,
    control = control,
    precomputed = B[permutation, , drop = FALSE]
  )

  expect_equal(
    permuted$embedding,
    original$embedding[permutation, , drop = FALSE],
    tolerance = 1e-8
  )
  expect_equal(permuted$relation_weights, original$relation_weights)
})

test_that("stored diffusion powers agree with direct sparse multiplication", {
  set.seed(1102)
  n <- 23L
  X <- matrix(rnorm(n * 3L), nrow = n)
  B <- matrix(rnorm(n * 4L), nrow = n)
  relations <- rapid_prepare_diffusion_relations(
    X,
    relations = list(topology = directed_cycle(n)),
    build_feature = FALSE,
    normalization = "random_walk",
    custom_symmetrize = "none",
    degree_cap = 2L
  )
  fit <- rapid_diffusion_encode(
    X,
    relations,
    precomputed = B,
    control = list(
      steps = 0:3,
      output_dim = 5L,
      gate = "uniform",
      direction_mode = "forward",
      store_propagations = TRUE,
      block_size = 2L
    )
  )

  A <- relations$relations$topology
  expected <- fit$base
  for (step in 1:3) {
    expected <- as.matrix(A %*% expected)
    expect_equal(
      fit$propagations$topology[[as.character(step)]],
      expected,
      tolerance = 1e-12,
      info = paste("diffusion step", step)
    )
  }
  expect_identical(fit$diagnostics$sparse_multiplies, 3L)
})

test_that("energy and label gates suppress an adversarial relation", {
  n <- 60L
  labels <- rep(1:2, each = n / 2L)
  within_to <- unlist(lapply(split(seq_len(n), labels), function(index) {
    c(index[-1L], index[1L])
  }))
  good <- Matrix::sparseMatrix(
    i = seq_len(n), j = within_to, x = 1, dims = c(n, n)
  )
  bad <- Matrix::sparseMatrix(
    i = seq_len(n),
    j = c((n / 2L + 1L):n, 1L:(n / 2L)),
    x = 1,
    dims = c(n, n)
  )
  X <- cbind(
    ifelse(labels == 1L, -1, 1),
    ifelse(labels == 1L, 1, -1),
    seq_len(n) / n
  )
  relations <- rapid_prepare_diffusion_relations(
    X,
    relations = list(good = good, adversarial = bad),
    build_feature = FALSE,
    normalization = "random_walk",
    custom_symmetrize = "none",
    degree_cap = 4L
  )

  fit <- rapid_diffusion_encode(
    X,
    relations,
    labels = labels,
    precomputed = X,
    control = list(
      steps = 0:2,
      output_dim = 3L,
      gate_temperature = 3,
      label_gate_weight = 2
    )
  )

  expect_gt(fit$relation_weights[["good"]], 0.99)
  expect_lt(fit$relation_weights[["adversarial"]], 1e-4)
  expect_lt(
    fit$relation_quality$good$energy,
    fit$relation_quality$adversarial$energy
  )
  expect_equal(fit$relation_quality$good$label_agreement, 1)
  expect_equal(fit$relation_quality$adversarial$label_agreement, 0)
})

test_that("empty and disconnected relations produce finite bounded output", {
  X <- matrix(2, nrow = 11L, ncol = 5L)
  empty <- rapid_prepare_diffusion_relations(X, build_feature = FALSE)
  no_edges <- rapid_prepare_diffusion_relations(
    X,
    relations = list(disconnected = Matrix::Matrix(0, 11L, 11L, sparse = TRUE)),
    build_feature = FALSE
  )

  empty_fit <- rapid_diffusion_encode(
    X, empty, control = list(output_dim = 4L)
  )
  disconnected_fit <- rapid_diffusion_encode(
    X, no_edges, control = list(output_dim = 4L)
  )

  expect_identical(dim(empty_fit$signatures), c(11L, 0L))
  expect_true(all(is.finite(empty_fit$embedding)))
  expect_true(all(is.finite(disconnected_fit$embedding)))
  expect_equal(disconnected_fit$relation_weights[["disconnected"]], 0)
  expect_identical(disconnected_fit$diagnostics$sparse_multiplies, 0L)
})

test_that("callbacks and precomputed base embeddings are interchangeable", {
  set.seed(1103)
  n <- 31L
  X <- matrix(rnorm(n * 6L), nrow = n)
  B <- cbind(X[, 1L] + X[, 2L], X[, 3L] - X[, 4L], X[, 5L:6L])
  relations <- rapid_prepare_diffusion_relations(
    X,
    relations = list(cycle = directed_cycle(n)),
    build_feature = FALSE,
    custom_symmetrize = "none"
  )
  control <- list(output_dim = 4L, gate = "uniform", seed = 29L)

  precomputed <- rapid_diffusion_encode(
    X, relations, control = control, precomputed = B
  )
  callback <- rapid_diffusion_encode(
    X,
    relations,
    control = control,
    encoder = function(X, relation_state) {
      expect_identical(relation_state, relations)
      B
    }
  )

  expect_equal(callback$base, precomputed$base)
  expect_equal(callback$embedding, precomputed$embedding, tolerance = 1e-12)
  expect_identical(callback$encoder_fit$base$source, "callback")
  expect_identical(precomputed$encoder_fit$base$source, "precomputed")
})

test_that("directed relations create capped forward and reverse channels", {
  n <- 40L
  X <- matrix(seq_len(n * 3L), nrow = n)
  star <- Matrix::sparseMatrix(
    i = rep(1L, n - 1L), j = 2:n, x = 1, dims = c(n, n)
  )
  relations <- rapid_prepare_diffusion_relations(
    X,
    relations = list(star = star),
    build_feature = FALSE,
    custom_symmetrize = "none",
    normalization = "random_walk",
    degree_cap = 6L
  )
  fit <- rapid_diffusion_encode(
    X,
    relations,
    control = list(
      steps = c(0L, 1L),
      output_dim = 3L,
      gate = "uniform",
      reverse_degree_cap = 4L,
      store_propagations = TRUE
    )
  )

  expect_setequal(
    fit$diagnostics$channels,
    c("star__forward", "star__reverse")
  )
  expect_equal(unname(fit$channel_weights), c(0.5, 0.5))
  expect_lte(fit$diagnostics$channel_max_out_degree[["star__forward"]], 6L)
  expect_lte(fit$diagnostics$channel_max_out_degree[["star__reverse"]], 4L)
  expect_lte(fit$diagnostics$channel_nnz, 12L)
  expect_named(
    fit$propagations,
    c("star__forward", "star__reverse"),
    ignore.order = TRUE
  )
})

test_that("fixed-hop work and retained storage scale with sparse edges", {
  make_fit <- function(n) {
    X <- cbind(seq_len(n) / n, sin(seq_len(n) / 13), cos(seq_len(n) / 17))
    relations <- rapid_prepare_diffusion_relations(
      X,
      relations = list(chain = directed_cycle(n)),
      build_feature = FALSE,
      normalization = "random_walk",
      custom_symmetrize = "none",
      degree_cap = 2L
    )
    rapid_diffusion_encode(
      X,
      relations,
      precomputed = X,
      control = list(
        steps = c(0L, 1L, 2L, 4L),
        output_dim = 3L,
        gate = "uniform",
        direction_mode = "forward",
        store_propagations = FALSE
      )
    )
  }

  small <- make_fit(400L)
  large <- make_fit(800L)
  expect_identical(small$diagnostics$total_relation_nnz, 400L)
  expect_identical(large$diagnostics$total_relation_nnz, 800L)
  expect_identical(small$diagnostics$sparse_multiplies, 4L)
  expect_identical(large$diagnostics$sparse_multiplies, 4L)
  expect_lt(as.numeric(object.size(large)) / as.numeric(object.size(small)), 2.5)
  expect_false(small$diagnostics$dense_node_square_allocated)
  expect_false(large$diagnostics$dense_node_square_allocated)

  implementation <- paste(vapply(
    c(
      ".rapid_diffusion_encode",
      ".rapid_propagate_channels",
      ".rapid_compress_diffusion"
    ),
    function(name) paste(deparse(get(name, envir = asNamespace("manifoldalign"))),
                          collapse = "\n"),
    character(1)
  ), collapse = "\n")
  expect_false(grepl("eigen\\(|svd\\(|irlba|RSpectra|PRIMME", implementation))
})
