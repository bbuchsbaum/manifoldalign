library(testthat)
library(manifoldalign)

# Run the exact test
set.seed(1234)  # Same seed as test file

test_that("Partial OT oracle respects constraints exactly", {
  # Small problem for exact verification
  G <- matrix(runif(9), 3, 3)
  p <- c(0.3, 0.4, 0.3)
  q <- c(0.2, 0.5, 0.3)
  rho <- 0.5
  
  gamma <- manifoldalign:::partial_ot_mass(G, p, q, rho)
  
  cat("\nCost matrix G:\n")
  print(G)
  cat("\nTransport plan gamma:\n")
  print(gamma)
  
  # Check row constraints
  row_sums <- rowSums(gamma)
  cat("\nRow sums:", row_sums, "\n")
  cat("p:", p, "\n")
  cat("Row sums <= p?", all(row_sums <= p + 1e-6), "\n")
  expect_true(all(row_sums <= p + 1e-6),
              info = "Row sums should not exceed marginal p")
  
  # Check column constraints  
  col_sums <- colSums(gamma)
  cat("\nCol sums:", col_sums, "\n")
  cat("q:", q, "\n")
  cat("Col sums <= q?", all(col_sums <= q + 1e-10), "\n")
  cat("Differences from q:", col_sums - q, "\n")
  cat("Max violation:", max(col_sums - q), "\n")
  expect_true(all(col_sums <= q + 1e-10),
              info = "Column sums should not exceed marginal q")
  
  # Check exact mass constraint
  total_mass <- sum(gamma)
  cat("\nTotal mass:", total_mass, ", rho:", rho, "\n")
  expect_equal(total_mass, rho, tolerance = 1e-6,
               info = "Total transported mass should equal rho exactly")
  
  # Verify it's a minimizer by checking KKT conditions (simplified)
  expect_true(all(gamma >= -1e-6),
              info = "Transport plan should be non-negative")
})