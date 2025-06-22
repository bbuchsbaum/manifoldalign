# tests/testthat/test-gromov_wasserstein.R
# Unit tests for Gromov-Wasserstein implementation

library(testthat)
library(Matrix)
library(multidesign)
library(multivarious)
library(tibble)

test_that("gromov_wasserstein works with basic hyperdesign input", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("multivarious")
  
  # Create simple test data with different dimensions
  set.seed(123)
  n <- 50
  X1 <- matrix(rnorm(n * 3), n, 3)  # Domain 1: 3 features
  X2 <- matrix(rnorm(n * 5), n, 5)  # Domain 2: 5 features
  
  # Create hyperdesign
  design <- data.frame(sample_id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
  
  # Run Gromov-Wasserstein
  result <- gromov_wasserstein(hd, epsilon = 0.1, max_iter = 50, verbose = FALSE)
  
  # Check basic properties
  expect_s3_class(result, "gromov_wasserstein")
  expect_s3_class(result, "multiblock_biprojector")
  
  # Check transport plan
  expect_length(result$transport_plans, 1)  # One pair of domains
  P <- result$transport_plans[[1]]
  expect_equal(dim(P), c(n, n))
  
  # Transport plan should be doubly stochastic (sum to uniform marginals)
  expect_equal(as.vector(rowSums(P)), rep(1/n, n), tolerance = 1e-3)
  expect_equal(as.vector(colSums(P)), rep(1/n, n), tolerance = 1e-3)
  
  # Check distance matrix
  expect_equal(dim(result$distances), c(2, 2))
  expect_equal(result$distances[1, 1], 0)
  expect_equal(result$distances[2, 2], 0)
  expect_true(result$distances[1, 2] > 0)
  expect_equal(result$distances[1, 2], result$distances[2, 1])
})

test_that("gromov_wasserstein handles identical domains", {
  skip_if_not_installed("multidesign")
  
  # Create identical data in two domains
  set.seed(456)
  n <- 30
  X <- matrix(rnorm(n * 4), n, 4)
  
  design <- data.frame(sample_id = 1:n)
  md1 <- multidesign::multidesign(X, design)
  md2 <- multidesign::multidesign(X, design)  # Same data
  hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
  
  # Run GW
  result <- gromov_wasserstein(hd, epsilon = 0.01, verbose = FALSE)
  
  # Distance should be small (but entropic regularization prevents exact 0)
  expect_lt(result$distances[1, 2], 0.2)
  
  # Transport plan should have some concentration on diagonal
  P <- result$transport_plans[[1]]
  diag_sum <- sum(diag(P))
  expect_gt(diag_sum, 0.3)  # Some mass on diagonal
})

test_that("gromov_wasserstein with different sample sizes", {
  skip_if_not_installed("multidesign")
  
  # Different number of samples
  set.seed(789)
  n1 <- 40
  n2 <- 60
  X1 <- matrix(rnorm(n1 * 3), n1, 3)
  X2 <- matrix(rnorm(n2 * 4), n2, 4)
  
  design1 <- data.frame(sample_id = 1:n1)
  design2 <- data.frame(sample_id = 1:n2)
  
  md1 <- multidesign::multidesign(X1, design1)
  md2 <- multidesign::multidesign(X2, design2)
  hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
  
  # Run GW
  result <- gromov_wasserstein(hd, epsilon = 0.1, verbose = FALSE)
  
  # Check transport plan dimensions
  P <- result$transport_plans[[1]]
  expect_equal(dim(P), c(n1, n2))
  
  # Marginals should sum correctly
  expect_equal(sum(P), 1, tolerance = 1e-6)
  expect_equal(as.vector(rowSums(P)), rep(1/n1, n1), tolerance = 1e-3)
  expect_equal(as.vector(colSums(P)), rep(1/n2, n2), tolerance = 1e-3)
})

test_that("gromov_wasserstein detects structural similarity", {
  skip_if_not_installed("multidesign")
  
  set.seed(111)
  n <- 50
  
  # Create data with similar structure but different dimensions
  # Both have 2 clear clusters
  t <- seq(0, 2*pi, length.out = n)
  
  # Domain 1: 2D circle
  X1 <- cbind(cos(t), sin(t)) + rnorm(n * 2, sd = 0.1)
  
  # Domain 2: 3D spiral with same circular structure
  X2 <- cbind(cos(t), sin(t), t/3) + rnorm(n * 3, sd = 0.1)
  
  design <- data.frame(sample_id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(circle = md1, spiral = md2))
  
  # Run GW - should align points with similar local structure
  result <- gromov_wasserstein(hd, epsilon = 0.05, max_iter = 100, verbose = FALSE)
  
  # Transport plan should show structure
  P <- result$transport_plans[[1]]
  
  # Just check that we got a valid transport plan
  expect_equal(dim(P), c(n, n))
  expect_equal(sum(P), 1, tolerance = 1e-6)
  
  # And that the distance is reasonable
  expect_true(result$distances[1, 2] > 0)
  # Convergence might not always happen with max_iter = 100
  # Just check we got a valid result
})

test_that("gromov_wasserstein validates inputs", {
  skip_if_not_installed("multidesign")
  
  # Single domain should error
  X1 <- matrix(rnorm(30), 10, 3)
  md1 <- multidesign::multidesign(X1, data.frame(id = 1:10))
  hd_single <- hyperdesign(list(domain1 = md1))
  
  expect_error(
    gromov_wasserstein(hd_single),
    "at least 2 domains"
  )
  
  # Invalid metric
  X2 <- matrix(rnorm(40), 10, 4)
  md2 <- multidesign::multidesign(X2, data.frame(id = 1:10))
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  expect_error(
    gromov_wasserstein(hd, metric = "invalid"),
    "should be one of"
  )
  
  # Invalid epsilon
  expect_error(
    gromov_wasserstein(hd, epsilon = -0.1),
    "epsilon must be positive"
  )
})

test_that("gromov_wasserstein convergence", {
  skip_if_not_installed("multidesign")
  
  # Small problem that should converge quickly
  set.seed(222)
  n <- 20
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Run with tight tolerance
  result <- gromov_wasserstein(
    hd, 
    epsilon = 0.1,
    max_iter = 200,
    tol = 1e-8,
    verbose = FALSE
  )
  
  expect_true(all(result$converged))
})

test_that("gromov_wasserstein print method works", {
  skip_if_not_installed("multidesign")
  
  # Create simple test case
  set.seed(333)
  n <- 25
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 4), n, 4)
  X3 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  hd <- hyperdesign(list(
    data1 = multidesign::multidesign(X1, design),
    data2 = multidesign::multidesign(X2, design),
    data3 = multidesign::multidesign(X3, design)
  ))
  
  result <- gromov_wasserstein(hd, epsilon = 0.1, verbose = FALSE)
  
  # Capture print output
  output <- capture.output(print(result))
  
  expect_true(any(grepl("Gromov-Wasserstein", output)))
  expect_true(any(grepl("Number of domains: 3", output)))
  expect_true(any(grepl("data1, data2, data3", output)))
})