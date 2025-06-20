#!/usr/bin/env Rscript
# Test a single PARROT case

library(testthat) 
library(manifoldalign)

# Source helpers
source("tests/testthat/test-parrot.R")

# Run just the transport plan test
test_that("PARROT transport plan has correct mathematical properties", {
  # Create test case
  n_nodes <- 25
  test_case <- create_aligned_networks(n_nodes, n_features = 3, noise_level = 0.1)
  
  # Run PARROT with tau = 0.01
  result <- parrot(test_case$hyperdesign, anchors = anchors, 
                   tau = 0.01, lambda = 0.1, max_iter = 50)
  
  transport_plan <- result$transport_plan
  
  # Check dimensions
  cat("\nTransport plan shape:", dim(transport_plan), "\n")
  
  # Check doubly stochastic property
  row_sums <- rowSums(transport_plan)
  col_sums <- colSums(transport_plan)
  
  cat("Row sums: min=", min(row_sums), "max=", max(row_sums), "mean=", mean(row_sums), "\n")
  cat("Col sums: min=", min(col_sums), "max=", max(col_sums), "mean=", mean(col_sums), "\n")
  cat("Total mass:", sum(transport_plan), "\n")
  
  # The actual test
  expect_true(all(abs(row_sums - 1) < 0.01), 
              info = paste("Row sums not 1"))
})