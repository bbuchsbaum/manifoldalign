test_that("build_dataset_graph_: complete graph has N choose 2 edges", {
  N <- 5
  E <- manifoldalign:::build_dataset_graph_(N, graph = "complete")
  expect_equal(nrow(E), N * (N - 1) / 2)
  expect_true(all(E$i < E$j))
})

test_that("build_dataset_graph_: MST returns N-1 edges and covers all nodes", {
  skip_on_cran()
  set.seed(1)
  N <- 6
  # Make synthetic domains with slight shifts
  domains <- lapply(seq_len(N), function(i) matrix(rnorm(100) + i * 0.2, 50, 2))
  E <- manifoldalign:::build_dataset_graph_(N, graph = "mst", domains = domains)
  expect_equal(nrow(E), N - 1)
  nodes <- sort(unique(c(E$i, E$j)))
  expect_equal(nodes, seq_len(N))
})

test_that("build_dataset_graph_: k-NN returns symmetric edges and reasonable count", {
  skip_on_cran()
  set.seed(2)
  N <- 7
  domains <- lapply(seq_len(N), function(i) matrix(rnorm(120) + i * 0.3, 40, 3))
  E <- manifoldalign:::build_dataset_graph_(N, graph = "knn", domains = domains, graph_opts = list(k = 2L))
  expect_true(nrow(E) >= N - 1)         # at least tree-like
  expect_true(all(E$i < E$j))           # upper-tri edges
  nodes <- sort(unique(c(E$i, E$j)))
  expect_equal(nodes, seq_len(N))
})

