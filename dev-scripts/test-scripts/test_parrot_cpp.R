#!/usr/bin/env Rscript

# Quick test of PARROT C++ optimization

library(manifoldalign)
library(Matrix)

# Source test helper
source("tests/testthat/test-parrot.R")

# Create small test case
set.seed(123)
n <- 50
C <- matrix(runif(n * n), n, n)

# Test R version - use dispatcher with flag
cat("Testing R implementation...\n")
manifoldalign::manifoldalign::set_parrot_use_rcpp(FALSE)
time_r <- system.time({
  result_r <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, tol = 1e-6)
})
cat("R time:", time_r[3], "seconds\n")

# Test C++ version - use dispatcher with flag
cat("\nTesting C++ implementation...\n")
manifoldalign::set_parrot_use_rcpp(TRUE)
time_cpp <- system.time({
  result_cpp <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, tol = 1e-6)
})
cat("C++ time:", time_cpp[3], "seconds\n")

# Check equivalence
max_diff <- max(abs(result_r - result_cpp))
cat("\nMax difference between R and C++:", max_diff, "\n")
cat("Speedup:", time_r[3] / time_cpp[3], "x\n")

# Check doubly stochastic property
cat("\nChecking doubly stochastic property:\n")
cat("Row sum error (C++):", max(abs(rowSums(result_cpp) - 1)), "\n")
cat("Col sum error (C++):", max(abs(colSums(result_cpp) - 1)), "\n")

# Test larger matrix
cat("\n\nTesting larger matrix (n=200)...\n")
n <- 200
C_large <- matrix(runif(n * n), n, n)

# Enable C++ globally
manifoldalign::set_parrot_use_rcpp(TRUE)
cat("Global C++ flag set to:", manifoldalign:::get_parrot_use_rcpp(), "\n")

# Time with C++
time_cpp_large <- system.time({
  result_cpp_large <- manifoldalign:::solve_sinkhorn_stabilized(C_large, tau = 0.01, max_iter = 100, tol = 1e-6)
})
cat("Large matrix C++ time:", time_cpp_large[3], "seconds\n")

# Disable C++ globally
manifoldalign::set_parrot_use_rcpp(FALSE)
time_r_large <- system.time({
  result_r_large <- manifoldalign:::solve_sinkhorn_stabilized(C_large, tau = 0.01, max_iter = 100, tol = 1e-6)
})
cat("Large matrix R time:", time_r_large[3], "seconds\n")
cat("Large matrix speedup:", time_r_large[3] / time_cpp_large[3], "x\n")

cat("\nOptimization test complete!\n")