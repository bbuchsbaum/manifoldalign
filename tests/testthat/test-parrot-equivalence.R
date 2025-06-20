# test-parrot-equivalence.R - Test numerical equivalence between R and C++ implementations
# -------------------------------------------------------------------------
# This ensures that the C++ optimizations produce identical results to the R version
# -------------------------------------------------------------------------

library(testthat)
library(Matrix)
library(manifoldalign)
library(multidesign)
library(multivarious)
library(tibble)

# Loosen tolerance for R/C++ equivalence testing due to floating point differences
EQUIVALENCE_TOL <- 1e-5

# Helper to create a hyperdesign with anchors for testing
create_test_hyperdesign_anchors <- function(X1, X2, anchor_pairs) {
  
  # Create design data frames with a specified anchor column
  design1 <- data.frame(id = 1:nrow(X1), anchor_col = NA)
  design1$anchor_col[anchor_pairs[,1]] <- anchor_pairs[,1]
  
  design2 <- data.frame(id = 1:nrow(X2), anchor_col = NA)
  design2$anchor_col[anchor_pairs[,2]] <- anchor_pairs[,1]
  
  # Create multidesign objects using the correct package
  md1 <- multidesign::multidesign(X1, design1)
  md2 <- multidesign::multidesign(X2, design2)
  
  # Create a simple list that mimics a hyperdesign structure
  hd <- list(A = md1, B = md2)
  class(hd) <- "hyperdesign"
  attr(hd, "anchors") <- data.frame(
    id1 = anchor_pairs[,1], 
    id2 = anchor_pairs[,2]
  )
  hd
}

# Helper to force both implementations
run_both_implementations <- function(fun_name, ...) {
  # Get the dispatcher function
  fun_dispatch <- get(fun_name, envir = asNamespace("manifoldalign"))
  
  # Run R version
  result_r <- fun_dispatch(..., use_rcpp = FALSE)
  
  # Run C++ version
  result_cpp <- fun_dispatch(..., use_rcpp = TRUE)
  
  list(r = result_r, cpp = result_cpp)
}

test_that("Sinkhorn R and C++ produce equivalent results", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  sizes <- c(10, 25, 50)
  tau_values <- c(0.001, 0.01, 0.1)
  
  for (n in sizes) {
    for (tau in tau_values) {
      # Create test cost matrix
      C <- matrix(runif(n * n), n, n)
      C <- (C + t(C)) / 2  # Make symmetric for consistency
      
      # Run both implementations
      results <- run_both_implementations(
        "solve_sinkhorn_stabilized",
        C = C, tau = tau, max_iter = 100, tol = 1e-6
      )
      
      # Check equivalence
      expect_equal(
        results$r, results$cpp,
        tolerance = EQUIVALENCE_TOL,
        info = sprintf("Sinkhorn equivalence failed for n=%d, tau=%.3f", n, tau)
      )
      
      # Verify properties are preserved
      expect_equal(rowSums(results$cpp), rep(1, n), tolerance = 1e-6)
      expect_equal(colSums(results$cpp), rep(1, n), tolerance = 1e-6)
    }
  }
})

test_that("Edge gradient R and C++ produce equivalent results", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 20
  
  # Create test networks
  X1 <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  X2 <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  
  # Create sparse adjacency matrices
  A1 <- Matrix::rsparsematrix(n_nodes, n_nodes, nnz = 30, rand.x = NULL)
  A1 <- A1 + t(A1)  # Make symmetric
  A1@x <- rep(1, length(A1@x))  # Binary adjacency
  
  A2 <- Matrix::rsparsematrix(n_nodes, n_nodes, nnz = 30, rand.x = NULL)
  A2 <- A2 + t(A2)
  A2@x <- rep(1, length(A2@x))
  
  # Create transport plan
  S <- matrix(runif(n_nodes * n_nodes), n_nodes, n_nodes)
  S <- S / sum(S) * n_nodes  # Normalize
  
  # Package as networks
  networks <- list(
    list(adjacency = A1, features = X1),
    list(adjacency = A2, features = X2)
  )
  
  # Run both implementations
  results <- run_both_implementations(
    "compute_edge_gradient",
    S = S, networks = networks
  )
  
  # Check equivalence
  expect_equal(
    results$r, results$cpp,
    tolerance = EQUIVALENCE_TOL,
    info = "Edge gradient equivalence failed"
  )
})

test_that("Squared distance computation R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  sizes <- list(c(10, 15), c(50, 50), c(100, 80))
  
  for (dims in sizes) {
    n1 <- dims[1]
    n2 <- dims[2]
    d <- 5
    
    X1 <- matrix(rnorm(n1 * d), n1, d)
    X2 <- matrix(rnorm(n2 * d), n2, d)
    
    # R implementation
    X1_sq <- rowSums(X1^2)
    X2_sq <- rowSums(X2^2)
    D_r <- outer(X1_sq, X2_sq, "+") - 2 * X1 %*% t(X2)
    
    # C++ implementation
    D_cpp <- manifoldalign:::compute_squared_distances_cpp(X1, X2)
    
    # Check equivalence
    expect_equal(
      D_r, D_cpp,
      tolerance = EQUIVALENCE_TOL,
      info = sprintf("Distance computation failed for %dx%d", n1, n2)
    )
  }
})

test_that("RWR computation R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 25
  n_anchors <- 5
  
  # Create test network
  X <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  
  # Build adjacency matrix
  knn <- 5
  A <- Matrix::Matrix(0, n_nodes, n_nodes)
  for (i in 1:n_nodes) {
    dists <- sqrt(rowSums((X - matrix(X[i,], n_nodes, 3, byrow = TRUE))^2))
    nearest <- order(dists)[2:(knn+1)]
    A[i, nearest] <- 1
    A[nearest, i] <- 1
  }
  
  # Create transition matrix
  deg <- Matrix::rowSums(A)
  deg[deg == 0] <- 1
  W <- Matrix::Diagonal(x = 1/deg) %*% A
  WT <- Matrix::t(W)
  
  # Create restart vectors
  E <- Matrix::sparseMatrix(
    i = 1:n_anchors,
    j = 1:n_anchors,
    x = 1,
    dims = c(n_nodes, n_anchors)
  )
  
  # R implementation
  sigma <- 0.15
  R_r <- as.matrix(E) / n_nodes
  for (iter in 1:20) {
    R_old <- R_r
    R_r <- (1 - sigma) * as.matrix(WT %*% R_r) + sigma * as.matrix(E)
    if (max(abs(R_r - R_old)) < 1e-6) break
  }
  
  # C++ implementation
  R_cpp <- manifoldalign:::compute_rwr_vectorized_cpp(
    as.matrix(WT), as.matrix(E), sigma, 20, 1e-6
  )
  
  # Check equivalence
  expect_equal(
    R_r, R_cpp,
    tolerance = EQUIVALENCE_TOL,
    info = "RWR computation equivalence failed"
  )
})

test_that("Full PARROT pipeline R and C++ are equivalent", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  
  set.seed(123)
  n_nodes <- 15
  
  # Create simple aligned networks
  X1 <- matrix(rnorm(n_nodes * 2), n_nodes, 2)
  X2 <- X1 + matrix(rnorm(n_nodes * 2, 0, 0.1), n_nodes, 2)
  
  # Create hyperdesign with anchors
  hd <- create_test_hyperdesign_anchors(X1, X2, cbind(1:3, 1:3))
  
  # Save current option
  old_option <- getOption("manifoldalign.parrot.use_rcpp", default = FALSE)
  
  # Run with R implementation
  options(manifoldalign.parrot.use_rcpp = FALSE)
  result_r <- parrot(hd, anchors = anchor_col, tau = 0.1, max_iter = 20)
  
  # Run with C++ implementation
  options(manifoldalign.parrot.use_rcpp = TRUE)
  result_cpp <- parrot(hd, anchors = anchor_col, tau = 0.1, max_iter = 20)
  
  # Restore option
  options(manifoldalign.parrot.use_rcpp = old_option)
  
  # Check key outputs are equivalent
  expect_equal(
    result_r$transport_plan,
    result_cpp$transport_plan,
    tolerance = 1e-6,  # Slightly relaxed for full pipeline
    info = "Full PARROT pipeline equivalence failed"
  )
  
  # Check that both produce valid transport plans
  expect_true(all(result_cpp$transport_plan >= 0))
  expect_equal(rowSums(result_cpp$transport_plan), rep(1, n_nodes), tolerance = 1e-6)
})

test_that("Performance improvement is achieved with C++", {
  skip_if_not(requireNamespace("Rcpp", quietly = TRUE))
  skip_on_cran()  # Skip performance tests on CRAN
  
  set.seed(123)
  n <- 100
  C <- matrix(runif(n * n), n, n)
  
  # Time R implementation
  time_r <- system.time({
    for (i in 1:10) {
      result_r <- solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, 
                                           tol = 1e-6, use_rcpp = FALSE)
    }
  })
  
  # Time C++ implementation
  time_cpp <- system.time({
    for (i in 1:10) {
      result_cpp <- solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100,
                                             tol = 1e-6, use_rcpp = TRUE)
    }
  })
  
  speedup <- time_r[3] / time_cpp[3]
  message(sprintf("Sinkhorn speedup: %.1fx (R: %.3fs, C++: %.3fs)", 
                  speedup, time_r[3], time_cpp[3]))
  
  # Expect at least 2x speedup
  expect_gt(speedup, 2)
})