#!/usr/bin/env Rscript

# Basic test of C++ function availability

# Clean load of package
if ("manifoldalign" %in% loadedNamespaces()) {
  unloadNamespace("manifoldalign")
}

# Load the package
library(manifoldalign)

# Check if the C++ function is available
cat("Checking C++ function availability...\n")

# List available C entry points
dll_info <- getLoadedDLLs()[["manifoldalign"]]
if (!is.null(dll_info)) {
  cat("DLL path:", dll_info[["path"]], "\n")
}

# Try to call the C++ function directly
tryCatch({
  # Create test matrix
  n <- 10
  C <- matrix(runif(n * n), n, n)
  
  # Call the C++ function directly 
  result <- manifoldalign:::solve_sinkhorn_stabilized_cpp(C, tau = 0.1, max_iter = 100, tol = 1e-6)
  
  cat("SUCCESS: C++ function called successfully!\n")
  cat("Result dimensions:", dim(result), "\n")
  cat("Row sums:", round(rowSums(result), 6), "\n")
}, error = function(e) {
  cat("ERROR calling C++ function:", e$message, "\n")
})

# Now try with dispatcher
cat("\nTesting with dispatcher function...\n")
tryCatch({
  C <- matrix(runif(10 * 10), 10, 10)
  
  # Test R version
  result_r <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.1, max_iter = 100, tol = 1e-6, use_cpp = FALSE)
  cat("R version successful\n")
  
  # Test C++ version
  result_cpp <- manifoldalign:::solve_sinkhorn_stabilized(C, tau = 0.1, max_iter = 100, tol = 1e-6, use_cpp = TRUE)
  cat("C++ version successful\n")
  
  max_diff <- max(abs(result_r - result_cpp))
  cat("Max difference:", max_diff, "\n")
}, error = function(e) {
  cat("ERROR in dispatcher test:", e$message, "\n")
})