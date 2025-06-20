test_that("cone_align_multiple works with hyperdesign objects", {
  
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
  result <- cone_align_multiple(hd, ncomp = 5, max_iter = 10)
  
  expect_s3_class(result, "cone_align_multiple")
  expect_s3_class(result, "multiblock_biprojector")
  
  # Check output structure
  expect_equal(length(result$assignment), 3)
  expect_equal(length(result$rotation), 3)
  expect_equal(nrow(result$s), n_nodes * 3)
  expect_equal(ncol(result$s), 5)
  
  # Check that assignments are valid permutations
  for (i in 1:3) {
    expect_equal(sort(result$assignment[[i]]), 1:n_nodes)
  }
  
  # Check that rotations are orthogonal
  for (i in 1:3) {
    Q <- result$rotation[[i]]
    expect_equal(nrow(Q), ncol(Q))
    expect_equal(Q %*% t(Q), diag(nrow(Q)), tolerance = 1e-10)
  }
})

test_that("cone_align_multiple works with list of matrices", {
  set.seed(456)
  
  # Create three similar matrices
  ms <- lapply(1:3, function(i) {
    matrix(rnorm(50 * 3), 50, 3)
  })
  
  result <- cone_align_multiple(ms, ncomp = 5, max_iter = 10)
  
  expect_s3_class(result, "cone_align_multiple")
  expect_equal(dim(result$s), c(150, 5))
  expect_equal(length(result$assignment), 3)
})

test_that("cone_align_multiple rejects invalid inputs", {
  # Too few domains
  expect_error(
    cone_align_multiple(list(matrix(1:10, 5, 2), matrix(1:10, 5, 2))),
    "at least three"
  )
  
  # Invalid input type
  expect_error(
    cone_align_multiple(matrix(1:10, 5, 2)),
    "list"
  )
  
  # Non-matrix elements
  expect_error(
    cone_align_multiple(list(1:10, 1:10, 1:10)),
    "matrix"
  )
})

test_that("cone_align_multiple handles different reference indices", {
  
  set.seed(789)
  
  # Create test data
  ms <- lapply(1:4, function(i) matrix(rnorm(30 * 3), 30, 3))
  
  # Test with different reference indices
  result1 <- cone_align_multiple(ms, ref_idx = 1, ncomp = 3, max_iter = 5)
  result2 <- cone_align_multiple(ms, ref_idx = 2, ncomp = 3, max_iter = 5)
  
  expect_equal(result1$ref_idx, 1)
  expect_equal(result2$ref_idx, 2)
  
  # Results should be different with different references
  expect_false(all(result1$assignment[[2]] == result2$assignment[[2]]))
})

test_that("cone_align_multiple convergence messages work", {
  
  set.seed(999)
  
  # Create easily alignable data for quick convergence
  base <- matrix(rnorm(20 * 3), 20, 3)
  ms <- lapply(1:3, function(i) base + rnorm(20 * 3, sd = 0.01))
  
  # Should converge quickly and print message
  expect_message(
    cone_align_multiple(ms, ncomp = 3, max_iter = 50, tol = 0.01),
    "Converged after"
  )
  
  # Should hit max iterations
  expect_message(
    cone_align_multiple(ms, ncomp = 3, max_iter = 1, tol = 1e-10),
    "Maximum iterations reached"
  )
}) 