#!/usr/bin/env Rscript

# Quick test of PARROT C++ optimization

library(manifoldalign)
library(Matrix)

# Create small test case
set.seed(123)
n <- 50
C <- matrix(runif(n * n), n, n)

# Test R version - use dispatcher with flag
cat("Testing R implementation...\n")
time_r <- system.time({
  result_r <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, tol = 1e-6, use_cpp = FALSE)
})
cat("R time:", time_r[3], "seconds\n")

# Test C++ version - use dispatcher with flag
cat("\nTesting C++ implementation...\n")
time_cpp <- system.time({
  result_cpp <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.01, max_iter = 100, tol = 1e-6, use_cpp = TRUE)
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

# Test with full PARROT
cat("\n\nTesting full PARROT with use_cpp flag...\n")

# Create test hyperdesign
X1 <- matrix(rnorm(30 * 3), 30, 3)
X2 <- X1 + matrix(rnorm(30 * 3, 0, 0.1), 30, 3)

# Create design data frames with anchors
design1 <- data.frame(
  node_id = 1:30,
  anchors = c(1:5, rep(NA, 25))
)
design2 <- data.frame(
  node_id = 1:30,
  anchors = c(1:5, rep(NA, 25))
)

# Create multidesign-like objects
domain1 <- list(x = X1, design = design1)
domain2 <- list(x = X2, design = design2)
hd <- list(domain1 = domain1, domain2 = domain2)
class(hd) <- c("hyperdesign", "list")

# Test with R
time_parrot_r <- system.time({
  result_parrot_r <- parrot(hd, anchors = anchors, tau = 0.1, max_iter = 20, use_cpp = FALSE)
})
cat("PARROT R time:", time_parrot_r[3], "seconds\n")

# Test with C++
time_parrot_cpp <- system.time({
  result_parrot_cpp <- parrot(hd, anchors = anchors, tau = 0.1, max_iter = 20, use_cpp = TRUE)
})
cat("PARROT C++ time:", time_parrot_cpp[3], "seconds\n")
cat("PARROT speedup:", time_parrot_r[3] / time_parrot_cpp[3], "x\n")

# Check transport plan equivalence
max_diff_parrot <- max(abs(result_parrot_r$transport_plan - result_parrot_cpp$transport_plan))
cat("\nMax difference in transport plans:", max_diff_parrot, "\n")

cat("\nOptimization test complete!\n")