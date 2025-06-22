library(manifoldalign)
library(multidesign)

# Comprehensive FPGW Test Suite
# Based on examples from the implementation

set.seed(123)

# Helper function to create test data
create_test_hyperdesign <- function(n = 30, dims = c(3, 5)) {
  # Domain 1: 3D data with two clusters
  X1 <- matrix(rnorm(n * dims[1]), n, dims[1])
  X1[1:(n/2), ] <- X1[1:(n/2), ] + 2
  
  # Domain 2: 5D data with similar structure
  X2 <- matrix(rnorm(n * dims[2]), n, dims[2])
  X2[1:(n/2), ] <- X2[1:(n/2), ] + 2
  
  # Create design
  design <- data.frame(id = 1:n, cluster = rep(1:2, each = n/2))
  
  # Create hyperdesign
  hd <- hyperdesign(list(
    visual = multidesign(X1, design),
    semantic = multidesign(X2, design)
  ))
  
  return(hd)
}

cat("=== FPGW Comprehensive Test Suite ===\n\n")

# Test 1: Basic FPGW functionality
cat("Test 1: Basic FPGW with balanced transport\n")
tryCatch({
  hd <- create_test_hyperdesign(n = 20)
  
  # Run FPGW with different omega values
  result1 <- fpgw(hd, omega1 = 0.5, max_iter = 50, verbose = FALSE)
  
  cat("✓ FPGW object created successfully\n")
  cat("  - Number of domains:", length(result1$domain_names), "\n")
  cat("  - Distance matrix computed:", !is.null(result1$distances), "\n")
  cat("  - Transport plans available:", length(result1$transport_plans), "\n")
  cat("  - All converged:", all(result1$converged), "\n")
  
  # Check distance properties
  D <- result1$distances
  cat("  - Distance matrix symmetric:", all(abs(D - t(D)) < 1e-10), "\n")
  cat("  - Diagonal is zero:", all(diag(D) < 1e-10), "\n")
  cat("  - Non-negative:", all(D >= -1e-10), "\n")
  
}, error = function(e) {
  cat("✗ Error in basic FPGW:", e$message, "\n")
})

cat("\n")

# Test 2: Mass-constrained variant
cat("Test 2: Mass-constrained FPGW (partial transport)\n")
tryCatch({
  # Create data with outliers
  n_clean <- 20
  n_noisy <- 30
  X1 <- matrix(rnorm(n_clean * 3), n_clean, 3)
  X1[1:10, ] <- X1[1:10, ] + 2
  
  X2_clean <- matrix(rnorm(n_clean * 5), n_clean, 5)
  X2_clean[1:10, ] <- X2_clean[1:10, ] + 2
  X2_outliers <- matrix(rnorm(10 * 5, sd = 5), 10, 5)
  X2 <- rbind(X2_clean, X2_outliers)
  
  design1 <- data.frame(id = 1:n_clean)
  design2 <- data.frame(id = 1:n_noisy)
  
  hd_noisy <- hyperdesign(list(
    clean = multidesign(X1, design1),
    noisy = multidesign(X2, design2)
  ))
  
  # Test with 75% mass transport
  result2 <- fpgw(hd_noisy, omega1 = 0.3, rho = 0.75, max_iter = 50, verbose = FALSE)
  
  cat("✓ Mass-constrained FPGW completed\n")
  cat("  - Rho parameter:", result2$rho, "\n")
  
  # Check transported mass
  P <- result2$transport_plans[[1]]
  actual_mass <- sum(P)
  expected_mass <- 0.75 * min(1, 1)  # min of total masses
  cat("  - Expected mass:", expected_mass, "\n")
  cat("  - Actual transported mass:", actual_mass, "\n")
  cat("  - Mass constraint satisfied:", abs(actual_mass - expected_mass) < 0.1, "\n")
  
}, error = function(e) {
  cat("✗ Error in mass-constrained FPGW:", e$message, "\n")
})

cat("\n")

# Test 3: TV-penalized variant
cat("Test 3: TV-penalized FPGW\n")
tryCatch({
  hd <- create_test_hyperdesign(n = 16)
  
  # Run with TV penalty
  result3 <- fpgw(hd, omega1 = 0.2, lambda = 0.1, max_iter = 50, verbose = FALSE)
  
  cat("✓ TV-penalized FPGW completed\n")
  cat("  - Lambda parameter:", result3$lambda, "\n")
  
  # Check sparsity of transport plan
  P <- result3$transport_plans[[1]]
  sparsity <- sum(P > 1e-10) / length(P)
  cat("  - Transport plan sparsity:", round(sparsity, 3), "\n")
  cat("  - TV penalty active:", result3$lambda > 0, "\n")
  
}, error = function(e) {
  cat("✗ Error in TV-penalized FPGW:", e$message, "\n")
})

cat("\n")

# Test 4: Multi-domain alignment
cat("Test 4: Multi-domain alignment (3+ domains)\n")
tryCatch({
  n <- 20
  X1 <- matrix(rnorm(n * 3), n, 3)
  X1[1:10, ] <- X1[1:10, ] + 2
  
  X2 <- matrix(rnorm(n * 4), n, 4)
  X2[1:10, ] <- X2[1:10, ] + 2
  
  X3 <- matrix(rnorm(n * 5), n, 5)
  X3[11:20, ] <- X3[11:20, ] + 1.5
  
  design <- data.frame(id = 1:n)
  
  hd_multi <- hyperdesign(list(
    modality1 = multidesign(X1, design),
    modality2 = multidesign(X2, design),
    modality3 = multidesign(X3, design)
  ))
  
  result4 <- fpgw(hd_multi, omega1 = 0.2, max_iter = 30, verbose = FALSE)
  
  cat("✓ Multi-domain FPGW completed\n")
  cat("  - Number of domains:", length(result4$domain_names), "\n")
  cat("  - Number of pairwise distances:", sum(upper.tri(result4$distances)), "\n")
  cat("  - Distance matrix:\n")
  print(round(result4$distances, 3))
  
}, error = function(e) {
  cat("✗ Error in multi-domain FPGW:", e$message, "\n")
})

cat("\n")

# Test 5: Edge cases
cat("Test 5: Edge cases and parameter validation\n")

# Test 5a: Invalid parameters
tryCatch({
  hd <- create_test_hyperdesign(n = 10)
  result <- fpgw(hd, omega1 = 0.5, rho = 0.5, lambda = 0.1)
  cat("✗ Should have failed with both rho and lambda\n")
}, error = function(e) {
  cat("✓ Correctly rejected both rho and lambda:", substr(e$message, 1, 50), "...\n")
})

# Test 5b: Single domain
tryCatch({
  X1 <- matrix(rnorm(20 * 3), 20, 3)
  design <- data.frame(id = 1:20)
  hd_single <- hyperdesign(list(only = multidesign(X1, design)))
  result <- fpgw(hd_single)
  cat("✗ Should have failed with single domain\n")
}, error = function(e) {
  cat("✓ Correctly rejected single domain:", substr(e$message, 1, 50), "...\n")
})

# Test 5c: Extreme omega values
tryCatch({
  hd <- create_test_hyperdesign(n = 10)
  result_feat <- fpgw(hd, omega1 = 1.0, max_iter = 20)  # Pure feature
  result_struct <- fpgw(hd, omega1 = 0.0, max_iter = 20)  # Pure structure
  cat("✓ Extreme omega values handled correctly\n")
}, error = function(e) {
  cat("✗ Error with extreme omega:", e$message, "\n")
})

cat("\n")

# Test 6: Convergence behavior
cat("Test 6: Convergence behavior\n")
tryCatch({
  hd <- create_test_hyperdesign(n = 30)
  
  # Test with very few iterations
  result_few <- fpgw(hd, omega1 = 0.5, max_iter = 5, verbose = FALSE)
  
  # Test with many iterations
  result_many <- fpgw(hd, omega1 = 0.5, max_iter = 200, verbose = FALSE)
  
  cat("✓ Convergence tests completed\n")
  cat("  - Few iterations converged:", all(result_few$converged), "\n")
  cat("  - Many iterations converged:", all(result_many$converged), "\n")
  cat("  - Distance difference:", 
      abs(result_few$distances[1,2] - result_many$distances[1,2]), "\n")
  
}, error = function(e) {
  cat("✗ Error in convergence test:", e$message, "\n")
})

cat("\n")

# Test 7: Print and summary methods
cat("Test 7: S3 methods (print, predict)\n")
tryCatch({
  hd <- create_test_hyperdesign(n = 20)
  result <- fpgw(hd, omega1 = 0.3, max_iter = 30, verbose = FALSE)
  
  # Test print
  cat("✓ Print method output:\n")
  capture.output(print(result))
  cat("  - Print completed without error\n")
  
  # Test predict (if implemented)
  if (exists("predict.fpgw")) {
    new_data <- matrix(rnorm(5 * 3), 5, 3)
    pred <- predict(result, new_data, from = 1, to = 2)
    cat("✓ Predict method available\n")
  } else {
    cat("  - Predict method not yet implemented\n")
  }
  
}, error = function(e) {
  cat("✗ Error in S3 methods:", e$message, "\n")
})

cat("\n=== Test Summary ===\n")
cat("Tests check FPGW implementation against expected behavior.\n")
cat("Note: Some tests may fail due to the network simplex segfault issue.\n")
cat("Once that's fixed, all algorithmic components should pass.\n")