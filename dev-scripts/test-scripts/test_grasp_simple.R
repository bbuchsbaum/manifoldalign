#!/usr/bin/env Rscript

# Simple test for full GRASP functionality
devtools::load_all()

# Create minimal test case using the same pattern as the debug script
set.seed(123)
n_nodes <- 8
X1 <- matrix(rnorm(n_nodes * 2), n_nodes, 2)
X2 <- X1 + matrix(rnorm(n_nodes * 2, 0, 0.1), n_nodes, 2)

# Create simple hyperdesign structure (like debug script)
create_test_hyperdesign <- function(X1, X2) {
  design1 <- data.frame(node_id = 1:nrow(X1))
  design2 <- data.frame(node_id = 1:nrow(X2))
  
  domain1 <- list(x = X1, design = design1)
  domain2 <- list(x = X2, design = design2)
  
  hd <- list(domain1 = domain1, domain2 = domain2)
  class(hd) <- c("hyperdesign", "list")
  
  hd
}

hd <- create_test_hyperdesign(X1, X2)

cat("Testing full GRASP interface...\n")

tryCatch({
  # Test the low-level functions directly first (like debug script)
  cat("Testing GRASP step by step...\n")
  
  bases <- manifoldalign:::compute_grasp_basis(hd, ncomp = 5, use_laplacian = TRUE)
  cat("✓ Spectral basis computed\n")
  
  descriptors <- manifoldalign:::compute_grasp_descriptors(bases, q_descriptors = 10, sigma = 1.0)
  cat("✓ Descriptors computed\n")
  
  alignment <- manifoldalign:::align_grasp_bases(bases[[1]], bases[[2]], 
                                                descriptors[[1]], descriptors[[2]], 
                                                lambda = 0.1)
  cat("✓ Base alignment computed\n")
  
  assignment <- manifoldalign:::compute_grasp_assignment(bases[[1]], bases[[2]], 
                                                        descriptors[[1]], descriptors[[2]], 
                                                        alignment$rotation, 
                                                        distance_method = "cosine",
                                                        solver_method = "linear")
  cat("✓ Assignment computed\n")
  
  cat("SUCCESS: All GRASP components work!\n")
  cat("Assignment:", assignment$assignment, "\n")
  
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  traceback()
})