test_that("grasp_multiset works with hyperdesign objects", {
  
  # Create test data with 3 domains
  set.seed(123)
  n_nodes <- 30
  n_features <- 5
  
  # Create three domains with some similarity
  base_data <- matrix(rnorm(n_nodes * n_features), n_nodes, n_features)
  
  # Create stratum-like structures directly
  domain1 <- list(
    x = base_data + rnorm(n_nodes * n_features, sd = 0.1),
    design = data.frame(node_id = 1:n_nodes)
  )
  
  domain2 <- list(
    x = base_data + rnorm(n_nodes * n_features, sd = 0.1),
    design = data.frame(node_id = 1:n_nodes)
  )
  
  domain3 <- list(
    x = base_data + rnorm(n_nodes * n_features, sd = 0.1),
    design = data.frame(node_id = 1:n_nodes)
  )
  
  # Create hyperdesign-like structure
  hd <- list(
    domain1 = domain1,
    domain2 = domain2,
    domain3 = domain3
  )
  class(hd) <- "hyperdesign"
  
  # Test basic functionality
  result <- grasp_multiset(hd, ncomp = 10, q_descriptors = 20, max_iter = 10)
  
  expect_s3_class(result, "grasp_multiset")
  
  # Check output structure
  expect_equal(length(result$embeddings), 3)
  expect_equal(length(result$permutations), 3)
  expect_equal(length(result$rotations), 3)
  expect_equal(length(result$mapping_diag), 3)
  expect_equal(nrow(result$embeddings[[1]]), n_nodes)
  expect_equal(ncol(result$embeddings[[1]]), 10)
  
  # Check that permutations are valid
  for (i in 1:3) {
    expect_equal(sort(result$permutations[[i]]), 1:n_nodes)
  }
  
  # Check that rotations are orthogonal
  for (i in 1:3) {
    Q <- result$rotations[[i]]
    QQt <- as.matrix(Q %*% Matrix::t(Q))
    expect_equal(QQt, diag(nrow(QQt)), tolerance = 1e-10)
  }
})

test_that("grasp_multiset works with list of matrices", {
  set.seed(456)
  
  # Create four similar matrices
  ms <- lapply(1:4, function(i) {
    matrix(rnorm(40 * 3), 40, 3)
  })
  
  result <- grasp_multiset(ms, ncomp = 8, q_descriptors = 15, max_iter = 10)
  
  expect_s3_class(result, "grasp_multiset")
  expect_equal(dim(result$embeddings[[1]]), c(40, 8))
  expect_equal(length(result$permutations), 4)
  expect_equal(result$anchor, 1L)
})

test_that("grasp_multiset rejects invalid inputs", {
  # Too few domains
  expect_error(
    grasp_multiset(list(matrix(1:10, 5, 2))),
    "at least two"
  )
  
  # Invalid input type
  expect_error(
    grasp_multiset(matrix(1:10, 5, 2)),
    "list"
  )
  
  # Non-matrix elements
  expect_error(
    grasp_multiset(list(1:10, 1:10, 1:10)),
    "matrix"
  )
})

test_that("grasp_multiset handles different anchor choices", {
  set.seed(789)
  
  # Create test data
  ms <- lapply(1:3, function(i) matrix(rnorm(25 * 3), 25, 3))
  
  # Test with different anchor indices
  result1 <- grasp_multiset(ms, anchor = 1, ncomp = 5, q_descriptors = 10, max_iter = 5)
  result2 <- grasp_multiset(ms, anchor = 2, ncomp = 5, q_descriptors = 10, max_iter = 5)
  result_mean <- grasp_multiset(ms, anchor = "mean", ncomp = 5, q_descriptors = 10, max_iter = 5)
  
  expect_equal(result1$anchor, 1)
  expect_equal(result2$anchor, 2)
  expect_equal(result_mean$anchor, "mean")
  
  # Results should be different with different anchors
  expect_false(all(result1$embeddings[[2]] == result2$embeddings[[2]]))
})

test_that("grasp_multiset handles convergence correctly", {
  set.seed(999)
  
  # Create easily alignable data for quick convergence
  base <- matrix(rnorm(20 * 3), 20, 3)
  ms <- lapply(1:3, function(i) base + rnorm(20 * 3, sd = 0.01))
  
  # Should converge with reasonable parameters
  result <- grasp_multiset(ms, ncomp = 3, q_descriptors = 10, 
                          max_iter = 50, tol = 0.01, lambda = 0.1)
  
  expect_true(result$converged || result$iterations < 50)
  
  # Should hit max iterations with very tight tolerance
  expect_warning(
    grasp_multiset(ms, ncomp = 3, q_descriptors = 10, 
                   max_iter = 2, tol = 1e-10),
    "maximum iterations"
  )
})

test_that("grasp_multiset produces reasonable alignment quality", {
  set.seed(1234)
  
  # Create graphs with known correspondence (shuffled versions)
  n <- 30
  base_graph <- matrix(rnorm(n * 4), n, 4)
  
  # Create permuted versions
  perm1 <- sample(n)
  perm2 <- sample(n)
  
  graphs <- list(
    base_graph,
    base_graph[perm1, ] + rnorm(n * 4, sd = 0.05),
    base_graph[perm2, ] + rnorm(n * 4, sd = 0.05)
  )
  
  result <- grasp_multiset(graphs, ncomp = 10, q_descriptors = 30, 
                          lambda = 0.2, max_iter = 30)
  
  # Check if we can recover some of the true correspondences
  # The first graph is unpermuted, so check alignment to it
  recovered1 <- result$permutations[[2]][perm1]
  recovered2 <- result$permutations[[3]][perm2]
  
  # Should recover at least some correct matches
  correct1 <- sum(recovered1 == 1:n)
  correct2 <- sum(recovered2 == 1:n)
  
  # With noise, we don't expect perfect recovery, but should do better than random
  expect_gt(correct1, n / 10)  # Better than 10% correct (random would be ~1/n)
  expect_gt(correct2, n / 10)
})

test_that("grasp_multiset handles solver options", {
  set.seed(5678)
  
  # Create larger graphs to test auction solver
  ms <- lapply(1:3, function(i) matrix(rnorm(100 * 5), 100, 5))
  
  result_auction <- grasp_multiset(ms, ncomp = 5, q_descriptors = 10, 
                                   solver = "auction", max_iter = 5)
  result_linear <- grasp_multiset(ms, ncomp = 5, q_descriptors = 10, 
                                  solver = "linear", max_iter = 5)
  
  expect_s3_class(result_auction, "grasp_multiset")
  expect_s3_class(result_linear, "grasp_multiset")
  
  # Both should produce valid permutations
  expect_equal(sort(result_auction$permutations[[2]]), 1:100)
  expect_equal(sort(result_linear$permutations[[2]]), 1:100)
}) 