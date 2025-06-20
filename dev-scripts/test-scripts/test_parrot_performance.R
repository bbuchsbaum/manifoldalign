#!/usr/bin/env Rscript

# Performance test for PARROT C++ optimization

library(manifoldalign)
library(Matrix)

cat("PARROT C++ Performance Benchmarking\n")
cat("=====================================\n\n")

# Test 1: Sinkhorn solver performance
cat("1. Sinkhorn Solver Performance\n")
cat("------------------------------\n")

sizes <- c(50, 100, 200, 400)
for (n in sizes) {
  C <- matrix(runif(n * n), n, n)
  
  # R version
  time_r <- system.time({
    result_r <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, 
                                                          tol = 1e-6, use_cpp = FALSE)
  })[3]
  
  # C++ version
  time_cpp <- system.time({
    result_cpp <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, 
                                                            tol = 1e-6, use_cpp = TRUE)
  })[3]
  
  speedup <- time_r / time_cpp
  max_diff <- max(abs(result_r - result_cpp))
  
  cat(sprintf("  n=%d: R=%.3fs, C++=%.3fs, Speedup=%.1fx, MaxDiff=%.2e\n", 
              n, time_r, time_cpp, speedup, max_diff))
}

# Test 2: Edge gradient computation
cat("\n2. Edge Gradient Performance\n")
cat("----------------------------\n")

n1 <- 100
n2 <- 100
d <- 10

# Create test data
X1 <- matrix(rnorm(n1 * d), n1, d)
X2 <- matrix(rnorm(n2 * d), n2, d)

# Create sparse adjacency matrices (10% density)
density <- 0.1
A1 <- Matrix(rbinom(n1 * n1, 1, density), n1, n1, sparse = TRUE)
A2 <- Matrix(rbinom(n2 * n2, 1, density), n2, n2, sparse = TRUE)
diag(A1) <- 0
diag(A2) <- 0

# Create test transport plan
S <- matrix(1/(n1*n2), n1, n2)

# Create network structures for R version
networks <- list(
  list(adjacency = A1, features = X1),
  list(adjacency = A2, features = X2)
)

# Time R version
time_r <- system.time({
  grad_r <- manifoldalign:::compute_edge_gradient(S, networks, use_cpp = FALSE)
})[3]

# Time C++ version
time_cpp <- system.time({
  grad_cpp <- manifoldalign:::compute_edge_gradient(S, networks, use_cpp = TRUE)
})[3]

speedup <- time_r / time_cpp
max_diff <- max(abs(grad_r - grad_cpp))

cat(sprintf("  n=%d, d=%d: R=%.3fs, C++=%.3fs, Speedup=%.1fx, MaxDiff=%.2e\n", 
            n1, d, time_r, time_cpp, speedup, max_diff))

# Test 3: Squared distance computation
cat("\n3. Squared Distance Performance\n")
cat("-------------------------------\n")

sizes <- c(100, 200, 400)
d <- 50

for (n in sizes) {
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- matrix(rnorm(n * d), n, d)
  
  # R version
  time_r <- system.time({
    D_r <- manifoldalign:::compute_squared_distances(X1, X2, use_cpp = FALSE)
  })[3]
  
  # C++ version  
  time_cpp <- system.time({
    D_cpp <- manifoldalign:::compute_squared_distances(X1, X2, use_cpp = TRUE)
  })[3]
  
  speedup <- time_r / time_cpp
  max_diff <- max(abs(D_r - D_cpp))
  
  cat(sprintf("  n=%d, d=%d: R=%.3fs, C++=%.3fs, Speedup=%.1fx, MaxDiff=%.2e\n", 
              n, d, time_r, time_cpp, speedup, max_diff))
}

cat("\nPerformance benchmarking complete!\n")