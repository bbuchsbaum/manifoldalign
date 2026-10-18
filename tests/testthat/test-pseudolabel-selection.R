# Pseudolabel representative selection and graph construction.

test_that("greedy diversity selection is farthest-first (uses max similarity)", {
  S <- matrix(0, 4, 4)
  S[1, 2] <- 0.9   # candidate 2 is a near-duplicate of candidate 1
  S[1, 4] <- 0.1
  S[3, 4] <- 0.1
  S <- S + t(S)
  S <- Matrix::Matrix(S, sparse = TRUE)
  sel <- manifoldalign:::.select_diverse_subset(S, 1:4, target_size = 3,
                                                diversity_weight = 0.9)
  # 1 (confidence tie -> first), then 3 (farthest from 1), then 4: candidate 2
  # is 0.9-similar to the selected 1 and must not be preferred.
  expect_identical(sel, c(1L, 3L, 4L))
})

test_that("greedy and farthest-first paths agree on pure diversity order", {
  set.seed(4)
  n <- 60
  P <- matrix(stats::runif(n * 2), n, 2)
  D <- as.matrix(dist(P))
  S <- Matrix::Matrix(exp(-D^2 / 0.05) * (D < 0.4), sparse = TRUE)
  Matrix::diag(S) <- 0
  greedy <- manifoldalign:::.select_diverse_subset(S, seq_len(n), 10, diversity_weight = 0.999)
  ffs <- manifoldalign:::.select_diverse_subset_advanced(S, seq_len(n), 10,
                                                         diversity_weight = 0.999)
  expect_identical(greedy, ffs)
})

test_that("threshold graph always has one vertex per sample", {
  # Highest-index samples have no edges: an edge-list graph would stop at 4.
  M <- Matrix::sparseMatrix(i = c(1, 3), j = c(2, 4), x = 1, dims = c(7, 7),
                            symmetric = TRUE)
  for (min_n in c(1L, 50000L)) {
    g <- manifoldalign:::.pseudolabel_threshold_graph(
      as(M, "generalMatrix"), edge_list_min_n = min_n
    )
    expect_equal(igraph::vcount(g), 7L)
    memb <- igraph::components(g)$membership
    expect_length(memb, 7L)
    expect_equal(memb[1], memb[2])
    expect_equal(memb[3], memb[4])
    expect_false(memb[1] == memb[3])
  }
})

test_that("assign_pseudolabels handles n >= 50000 with isolated trailing rows", {
  n <- 50000L
  sim <- Matrix::sparseMatrix(i = c(1, 3), j = c(2, 4), x = c(0.9, 0.95),
                              dims = c(n, n), symmetric = TRUE)
  res <- assign_pseudolabels(as(sim, "generalMatrix"), sim_threshold = 0.5,
                             min_clusters = 1, max_clusters = 10,
                             min_cluster_size = 2)
  expect_length(res$labels, n)
  expect_equal(res$n_clusters, 2)
  expect_true(all(!is.na(res$labels[1:4])))
  expect_true(all(is.na(res$labels[5:n])))
  expect_equal(as.character(res$labels[1]), as.character(res$labels[2]))
})

test_that("assign_pseudolabels seed is local", {
  syn <- create_synthetic_similarity_matrix(n_samples = 120, n_clusters = 12)
  set.seed(10)
  before <- .Random.seed
  a <- assign_pseudolabels(syn$sim_matrix, min_clusters = 2, max_clusters = 4,
                           seed = 5)
  expect_identical(.Random.seed, before)
  b <- assign_pseudolabels(syn$sim_matrix, min_clusters = 2, max_clusters = 4,
                           seed = 5)
  expect_identical(a$representatives, b$representatives)
})
