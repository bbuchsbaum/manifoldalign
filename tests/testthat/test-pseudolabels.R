context("pseudolabel diagnostics")

test_that("assign_pseudolabels yields consistent stats for symmetric storage", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("igraph")

  edges <- data.frame(
    i = c(1, 1, 2),
    j = c(2, 3, 3),
    x = c(0.9, 0.8, 0.7)
  )

  sim_sym <- Matrix::sparseMatrix(
    i = edges$i,
    j = edges$j,
    x = edges$x,
    symmetric = TRUE,
    dims = c(3, 3)
  )

  sim_gen <- as(sim_sym, "dgCMatrix")

  set.seed(123)
  res_sym <- assign_pseudolabels(
    sim_matrix = sim_sym,
    sim_threshold = 0.5,
    min_clusters = 1,
    max_clusters = 3,
    seed = 99,
    verbose = FALSE
  )

  set.seed(123)
  res_gen <- assign_pseudolabels(
    sim_matrix = sim_gen,
    sim_threshold = 0.5,
    min_clusters = 1,
    max_clusters = 3,
    seed = 99,
    verbose = FALSE
  )

  expect_equal(res_sym$cluster_info$avg_similarity, res_gen$cluster_info$avg_similarity)
  expect_equal(res_sym$cluster_info$observed_similarity, res_gen$cluster_info$observed_similarity)
  expect_equal(res_sym$cluster_info$edge_density, res_gen$cluster_info$edge_density)
})

test_that("assign_pseudolabels returns NA labels when no edges remain", {
  skip_if_not_installed("Matrix")

  sim_null <- Matrix::sparseMatrix(
    i = c(1, 2, 2, 3),
    j = c(2, 1, 3, 2),
    x = rep(0.2, 4),
    dims = c(3, 3)
  )

  res <- suppressWarnings(assign_pseudolabels(
    sim_matrix = sim_null,
    sim_threshold = 0.9,
    min_clusters = 1,
    max_clusters = 3,
    seed = 42,
    verbose = FALSE
  ))

  expect_equal(res$n_clusters, 0)
  expect_true(all(is.na(res$labels)))
})

test_that("high_sim_pseudolabels maps component indices back to data rows", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("RcppAnnoy")
  ann_ctor <- get0("AnnoyAngular", envir = asNamespace("RcppAnnoy"))
  if (is.null(ann_ctor) || !inherits(ann_ctor, "C++Class")) {
    skip("RcppAnnoy::AnnoyAngular constructor not available")
  }
  ann_instance <- try(new(ann_ctor, 2), silent = TRUE)
  if (inherits(ann_instance, "try-error")) {
    skip("Unable to instantiate AnnoyAngular index")
  }
  ann_instance$unload()

  strata <- list(
    list(
      x = matrix(c(1, 0, 0, 1), ncol = 2, byrow = TRUE),
      design = data.frame(domain = "A")
    ),
    list(
      x = matrix(c(1, 0, 0, 1), ncol = 2, byrow = TRUE),
      design = data.frame(domain = "B")
    )
  )

  set.seed(2024)
  labels <- high_sim_pseudolabels(
    strata = strata,
    k = 1,
    cos_thresh = 0.95,
    min_size = 2,
    ann_trees = 2,
    verbose = FALSE
  )

  expect_length(labels, 4)
  expect_true(!is.na(labels[1]))
  expect_true(!is.na(labels[2]))
  expect_equal(labels[1], labels[3])
  expect_equal(labels[2], labels[4])
  expect_false(labels[1] == labels[2])
})
