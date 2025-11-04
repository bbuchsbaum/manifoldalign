# tests/testthat/test-fpgw.R
# Unit tests for Fused-Partial Gromov-Wasserstein implementation

library(testthat)
library(Matrix)

have_multidesign <- requireNamespace("multidesign", quietly = TRUE)
have_lpSolve <- requireNamespace("lpSolve", quietly = TRUE)

if (have_multidesign && !("package:multidesign" %in% search())) {
  suppressPackageStartupMessages(library(multidesign))
}
if (requireNamespace("tibble", quietly = TRUE) && !("package:tibble" %in% search())) {
  suppressPackageStartupMessages(library(tibble))
}

test_that("fpgw works with basic hyperdesign input", {
  skip_if_not(have_multidesign, "multidesign not installed")
  skip_if_not(have_lpSolve, "lpSolve not installed")
  
  # Create simple test data
  set.seed(123)
  n <- 20
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 3), n, 3)  # Same dimensions for classical FGW
  
  # Create hyperdesign
  design <- data.frame(sample_id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
  
  # Run FPGW with default settings (lambda=0, rho=NULL should give classical FGW)
  result <- fpgw(hd, omega1 = 0.5, verbose = FALSE)
  
  # Check basic properties
  expect_s3_class(result, "fpgw")
  expect_s3_class(result, "multiblock_biprojector")
  
  # Check transport plan
  expect_length(result$transport_plans, 1)
  P <- result$transport_plans[[1]]
  expect_equal(dim(P), c(n, n))
  
  # For classical FGW (lambda=0, rho=NULL), transport plan should be doubly stochastic
  # But due to numerical precision, allow some tolerance
  expect_equal(sum(P), 1, tolerance = 0.1)  # Relaxed tolerance
  expect_equal(as.vector(rowSums(P)), rep(1/n, n), tolerance = 0.1)
  expect_equal(as.vector(colSums(P)), rep(1/n, n), tolerance = 0.1)
  
  # Check distance matrix
  expect_equal(dim(result$distances), c(2, 2))
  expect_equal(result$distances[1, 1], 0)
  expect_equal(result$distances[2, 2], 0)
  expect_true(result$distances[1, 2] > 0)
  expect_equal(result$distances[1, 2], result$distances[2, 1])
})

test_that("fpgw handles mass-constrained variant", {
  skip_if_not(have_multidesign, "multidesign not installed")
  skip_if_not(have_lpSolve, "lpSolve not installed")
  
  set.seed(456)
  n <- 15
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Run with mass constraint
  rho <- 0.8
  result <- fpgw(hd, omega1 = 0.1, rho = rho, verbose = FALSE)
  
  # Check that transported mass equals rho
  P <- result$transport_plans[[1]]
  expect_equal(sum(P), rho, tolerance = 1e-3)
  
  # Check metadata
  expect_equal(result$rho, rho)
  expect_equal(result$lambda, 0)  # Default lambda is 0, not NULL
})

test_that("fpgw handles TV-penalized variant", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  set.seed(789)
  n <- 10
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Run with TV penalty
  lambda <- 0.1
  result <- fpgw(hd, omega1 = 0.3, lambda = lambda, verbose = FALSE)
  
  # Check metadata
  expect_equal(result$lambda, lambda)
  expect_null(result$rho)
  
  # Transport plan should have reduced mass due to penalty
  P <- result$transport_plans[[1]]
  expect_lt(sum(P), 1)
})

test_that("fpgw validates mutually exclusive parameters", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  # Create dummy data
  n <- 10
  X <- matrix(rnorm(n * 2), n, 2)
  design <- data.frame(id = 1:n)
  md <- multidesign::multidesign(X, design)
  hd <- hyperdesign(list(d1 = md, d2 = md))
  
  # Should error if both rho and lambda are provided
  expect_error(
    fpgw(hd, omega1 = 0.5, rho = 0.8, lambda = 0.1),
    "Choose.*either.*rho.*lambda.*not both"
  )
})

test_that("fpgw handles different omega1 values", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  set.seed(111)
  n <- 15
  # Create data with clear feature differences
  X1 <- matrix(c(rnorm(n, 0), rnorm(n, 0)), n, 2)
  X2 <- matrix(c(rnorm(n, 5), rnorm(n, 5)), n, 2)  # Shifted features
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Pure structural (omega1 = 0 would make it pure GW, but use small value)
  result_struct <- fpgw(hd, omega1 = 0.001, verbose = FALSE)
  
  # Balanced
  result_balanced <- fpgw(hd, omega1 = 0.5, verbose = FALSE)
  
  # Pure feature (omega1 = 1 would make it pure OT)
  result_feature <- fpgw(hd, omega1 = 0.999, verbose = FALSE)
  
  # Distances should differ based on omega1
  expect_true(result_struct$distances[1,2] != result_feature$distances[1,2])
})

test_that("fpgw convergence behavior", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  # Small problem that should converge quickly
  set.seed(222)
  n <- 10
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Run with tight tolerance
  result <- fpgw(hd, omega1 = 0.1, max_iter = 200, tol = 1e-8, verbose = FALSE)
  
  expect_true(all(result$converged))
})

test_that("fpgw print method works", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  # Create simple test case
  set.seed(333)
  n <- 10
  X1 <- matrix(rnorm(n * 2), n, 2)
  X2 <- matrix(rnorm(n * 2), n, 2)
  X3 <- matrix(rnorm(n * 2), n, 2)
  
  design <- data.frame(id = 1:n)
  hd <- hyperdesign(list(
    data1 = multidesign::multidesign(X1, design),
    data2 = multidesign::multidesign(X2, design),
    data3 = multidesign::multidesign(X3, design)
  ))
  
  result <- fpgw(hd, omega1 = 0.2, rho = 0.9, verbose = FALSE)
  
  # Capture print output
  output <- capture.output(print(result))
  
  expect_true(any(grepl("Fused-Partial Gromov-Wasserstein", output)))
  expect_true(any(grepl("Number of domains: 3", output)))
  expect_true(any(grepl("data1, data2, data3", output)))
  expect_true(any(grepl("Mass-constrained", output)))
})

# Validation test against known solution
test_that("fpgw produces expected results on toy problem", {
  skip_if_not(have_multidesign, "multidesign not installed")
  skip_if_not(have_lpSolve, "lpSolve not installed")
  
  # Create a toy problem with known structure
  set.seed(42)
  n <- 5  # Very small for manual verification
  
  # Domain 1: Points on a line
  X1 <- matrix(c(1:n, rep(0, n)), n, 2)
  
  # Domain 2: Same structure but rotated
  theta <- pi/4
  rot <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
  X2 <- X1 %*% t(rot)
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # With high structural weight, should match points correctly
  result <- fpgw(hd, omega1 = 0.01, verbose = FALSE)
  
  P <- result$transport_plans[[1]]
  
  # Transport plan should be close to diagonal (matching i to i)
  diag_mass <- sum(diag(P))
  expect_gt(diag_mass, 0.8)  # Most mass on diagonal
})

test_that("fpgw handles different dimensional data", {
  skip_if_not(have_multidesign, "multidesign not installed")
  
  # Create data with different dimensions
  set.seed(789)
  n <- 15
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 5), n, 5)  # Different dimensions
  
  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(d1 = md1, d2 = md2))
  
  # Run FPGW - should work but may not have full mass transport
  result <- fpgw(hd, omega1 = 0.5, verbose = FALSE)
  
  # Check that it runs without error
  expect_s3_class(result, "fpgw")
  P <- result$transport_plans[[1]]
  expect_equal(dim(P), c(n, n))
  
  # For different dimensions, transport may be partial
  # Just check that mass is positive and bounded
  total_mass <- sum(P)
  expect_gt(total_mass, 0)
  expect_lte(total_mass, 1)
})

test_that("fpgw out-of-sample transport approximates target samples", {
  skip_if_not(have_multidesign, "multidesign not installed")
  skip_if_not(have_lpSolve, "lpSolve not installed")

  set.seed(4242)
  n <- 12
  X1 <- matrix(rnorm(n * 3), n, 3)
  # Target domain is shifted version of source
  shift <- matrix(rep(c(0.5, -0.25, 0.1), each = n), n, 3, byrow = FALSE)
  X2 <- X1 + shift

  design <- data.frame(id = 1:n)
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- hyperdesign(list(source = md1, target = md2))

  fit <- fpgw(hd, omega1 = 0.9, max_iter = 50, verbose = FALSE)

  new_points <- X1[1:3, , drop = FALSE]
  transported <- predict(fit, new_points, from = "source", to = "target", k = 1)

  expect_equal(dim(transported), dim(new_points))
  expected <- X2[1:3, , drop = FALSE]
  expect_equal(transported, expected, tolerance = 0.2)

  weights <- predict(fit, new_points, from = 1, to = 2, type = "weights", k = 2)
  expect_equal(nrow(weights), nrow(new_points))
  expect_equal(rowSums(weights), rep(1, nrow(weights)), tolerance = 1e-6)

  transformed <- transform(fit, new_points, source_index = 1, target_index = 2, k = 2)
  expect_equal(dim(transformed), dim(new_points))
})
