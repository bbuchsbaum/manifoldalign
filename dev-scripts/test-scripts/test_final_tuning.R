#!/usr/bin/env Rscript

# Final tuning test to find optimal parameters
devtools::load_all()
library(Matrix)

set.seed(42)

# Create a very simple, structured isomorphic test case
create_simple_isomorphic <- function(n_nodes = 8) {
  # Simple circular graph
  adj_matrix <- Matrix::sparseMatrix(
    i = c(1:n_nodes, 1),
    j = c(2:n_nodes, 1, 1),
    x = 1,
    dims = c(n_nodes, n_nodes)
  )
  adj_matrix <- adj_matrix + Matrix::t(adj_matrix)
  
  # Simple features based on position
  angles <- 2 * pi * (1:n_nodes) / n_nodes
  X1 <- cbind(cos(angles), sin(angles), rnorm(n_nodes, 0, 0.01))
  
  # Simple permutation: reverse order
  perm <- n_nodes:1
  X2 <- X1[perm, ] + matrix(rnorm(n_nodes * 3, 0, 0.01), n_nodes, 3)
  
  # Create hyperdesign
  design1 <- data.frame(node_id = 1:n_nodes)
  design2 <- data.frame(node_id = 1:n_nodes)
  domain1 <- list(x = X1, design = design1)
  domain2 <- list(x = X2, design = design2)
  hd <- list(domain1 = domain1, domain2 = domain2)
  class(hd) <- c("hyperdesign", "list")
  
  list(hyperdesign = hd, true_permutation = perm, adj_matrix = adj_matrix)
}

cat("=== FINAL PARAMETER TUNING ===\n")

# Test different parameter combinations
param_grid <- expand.grid(
  ncomp = c(4, 6),
  q_desc = c(20, 40),
  sigma = c(0.8, 1.2),
  lambda = c(0.02, 0.05, 0.1)
)

n_nodes <- 8
test_case <- create_simple_isomorphic(n_nodes)

best_accuracy <- 0
best_params <- NULL

cat("Testing", nrow(param_grid), "parameter combinations...\n")

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  
  result <- tryCatch({
    bases <- manifoldalign:::compute_grasp_basis(test_case$hyperdesign, 
                                                ncomp = params$ncomp, use_laplacian = TRUE)
    descriptors <- manifoldalign:::compute_grasp_descriptors(bases, 
                                                            q_descriptors = params$q_desc, 
                                                            sigma = params$sigma)
    alignment <- manifoldalign:::align_grasp_bases(bases[[1]], bases[[2]], 
                                                  descriptors[[1]], descriptors[[2]], 
                                                  lambda = params$lambda)
    assignment <- manifoldalign:::compute_grasp_assignment(bases[[1]], bases[[2]], 
                                                          descriptors[[1]], descriptors[[2]], 
                                                          alignment$rotation, 
                                                          distance_method = "cosine", 
                                                          solver_method = "linear")
    
    accuracy <- sum(assignment$assignment == test_case$true_permutation) / n_nodes
    list(accuracy = accuracy, success = TRUE)
  }, error = function(e) {
    list(accuracy = 0, success = FALSE, error = e$message)
  })
  
  if (result$success && result$accuracy > best_accuracy) {
    best_accuracy <- result$accuracy
    best_params <- params
  }
  
  cat(sprintf("Combo %d: ncomp=%d, q=%d, sigma=%.1f, lambda=%.2f -> %.1f%%\n", 
             i, params$ncomp, params$q_desc, params$sigma, params$lambda, result$accuracy * 100))
}

cat("\nBest result:\n")
cat("Accuracy:", sprintf("%.1f%%", best_accuracy * 100), "\n")
cat("Parameters:", "\n")
print(best_params)

if (best_accuracy >= 0.7) {
  cat("\nSUCCESS: Found parameters that achieve ≥70% accuracy!\n")
} else {
  cat("\nNeed to further adjust parameters or test design...\n")
}