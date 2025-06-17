# test-gpca_align.R - Comprehensive tests for GPCA alignment
# -------------------------------------------------------------------------
# Test that gpca_align correctly performs Generalized PCA alignment on 
# multi-domain data, handling similarity matrices and trade-off parameters
# -------------------------------------------------------------------------

library(testthat)
library(Matrix)
library(dplyr)
library(multidesign)
library(multivarious)
library(manifoldalign)
library(neighborweights)

set.seed(42)

# Helper to create test hyperdesign objects
create_test_hyperdesign <- function(n_per_domain = 20, n_domains = 2, n_features = 5, 
                                    label_overlap = TRUE) {
  # Create domains with some systematic differences
  domains <- list()
  label_list <- list()
  
  for (i in 1:n_domains) {
    # Each domain has slightly shifted data
    X <- matrix(rnorm(n_per_domain * n_features, mean = i - 1), n_per_domain, n_features)
    
    if (label_overlap) {
      # Domains share same label structure
      labels <- rep(c("A", "B"), length.out = n_per_domain)
    } else {
      # Domains have different labels
      labels <- rep(paste0("D", i, c("A", "B")), length.out = n_per_domain)
    }
    
    # Create multidesign object
    md <- multidesign::multidesign(X, data.frame(lbl = factor(labels)))
    domains[[paste0("domain", i)]] <- md
    label_list[[i]] <- labels
  }
  
  # Return hyperdesign and label info
  list(
    hd = multidesign::hyperdesign(domains),
    labels = label_list
  )
}

# Test 1: Basic functionality and parameter handling
test_that("gpca_align performs basic alignment with proper structure", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("PRIMME")
  
  # Create test data with 3 domains
  test_data <- create_test_hyperdesign(n_per_domain = 30, n_domains = 3, n_features = 4)
  hd <- test_data$hd
  
  # Test basic alignment with default parameters
  result <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2)
  
  # Check result structure
  expect_s3_class(result, "multiblock_biprojector")
  expect_true("v" %in% names(result))  # Projection matrix
  expect_true("s" %in% names(result))  # Scores
  expect_true("sdev" %in% names(result))  # Standard deviations
  expect_true("block_indices" %in% names(result))
  expect_true("labels" %in% names(result))
  
  # Check dimensions
  total_samples <- sum(sapply(hd, function(x) nrow(x$x)))
  expect_equal(nrow(result$s), total_samples)
  expect_equal(ncol(result$s), 2)  # ncomp = 2
  
  # Check scores are finite
  expect_true(all(is.finite(result$s)))
  
  # Test with different u values (trade-off parameter)
  result_within <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, u = 0.9)  # Favor within-domain
  result_between <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, u = 0.1)  # Favor between-domain
  
  # Both should produce valid results
  expect_true(all(is.finite(result_within$s)))
  expect_true(all(is.finite(result_between$s)))
  
  # Results should differ based on u parameter
  expect_false(all(abs(result_within$s - result_between$s) < 1e-10))
  
  # Test with custom similarity function
  custom_simfun <- function(labels) {
    # Simple custom similarity: exact match = 1, otherwise = 0.1
    S <- outer(labels, labels, function(x, y) ifelse(x == y, 1, 0.1))
    Matrix::Matrix(S, sparse = TRUE)
  }
  
  result_custom <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, simfun = custom_simfun)
  expect_true(all(is.finite(result_custom$s)))
  
  # Test regularization parameter effect
  result_low_reg <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, lambda = 0.01)
  result_high_reg <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, lambda = 1.0)
  
  # Both should work but produce different results
  expect_true(all(is.finite(result_low_reg$s)))
  expect_true(all(is.finite(result_high_reg$s)))
  expect_false(all(abs(result_low_reg$s - result_high_reg$s) < 1e-10))
})

# Test 2: Multi-domain alignment quality and cross-domain relationships
test_that("gpca_align correctly handles cross-domain alignment", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("PRIMME")
  
  # Create test data where domains have clear correspondence
  # Domain 1: Original data
  # Domain 2: Rotated version
  # Domain 3: Scaled version
  set.seed(123)
  n_samples <- 40
  
  # Base data with clear cluster structure
  base_data <- rbind(
    matrix(rnorm(20 * 3, mean = c(-2, 0, 0)), 20, 3),  # Cluster A
    matrix(rnorm(20 * 3, mean = c(2, 0, 0)), 20, 3)    # Cluster B
  )
  labels <- rep(c("A", "B"), each = 20)
  
  # Create rotation matrix
  theta <- pi/6  # 30 degrees
  rot_mat <- matrix(c(cos(theta), -sin(theta), 0,
                      sin(theta), cos(theta), 0,
                      0, 0, 1), 3, 3)
  
  # Create three related domains
  domain1_data <- base_data
  domain2_data <- base_data %*% rot_mat  # Rotated
  domain3_data <- base_data * 1.5  # Scaled
  
  # Build hyperdesign
  md1 <- multidesign::multidesign(domain1_data, data.frame(lbl = factor(labels)))
  md2 <- multidesign::multidesign(domain2_data, data.frame(lbl = factor(labels)))
  md3 <- multidesign::multidesign(domain3_data, data.frame(lbl = factor(labels)))
  
  hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2, domain3 = md3))
  
  # Perform alignment favoring between-domain relationships
  result <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, u = 0.2)
  
  # Extract scores for each domain
  scores_d1 <- result$s[1:40, ]
  scores_d2 <- result$s[41:80, ]
  scores_d3 <- result$s[81:120, ]
  
  # Check that same-label samples across domains are aligned
  # Compute centroids for each label in each domain
  centroid_A_d1 <- colMeans(scores_d1[labels == "A", ])
  centroid_B_d1 <- colMeans(scores_d1[labels == "B", ])
  centroid_A_d2 <- colMeans(scores_d2[labels == "A", ])
  centroid_B_d2 <- colMeans(scores_d2[labels == "B", ])
  centroid_A_d3 <- colMeans(scores_d3[labels == "A", ])
  centroid_B_d3 <- colMeans(scores_d3[labels == "B", ])
  
  # Cross-domain distances for same labels should be smaller than different labels
  same_label_dist_A <- mean(c(
    dist(rbind(centroid_A_d1, centroid_A_d2)),
    dist(rbind(centroid_A_d1, centroid_A_d3)),
    dist(rbind(centroid_A_d2, centroid_A_d3))
  ))
  
  diff_label_dist <- mean(c(
    dist(rbind(centroid_A_d1, centroid_B_d2)),
    dist(rbind(centroid_A_d1, centroid_B_d3)),
    dist(rbind(centroid_B_d1, centroid_A_d2))
  ))
  
  # Same-label centroids should be closer across domains
  expect_lt(same_label_dist_A, diff_label_dist)
  
  # Test extreme u values
  # u = 1: Only within-domain (should preserve domain structure)
  result_within_only <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, u = 1.0)
  
  # u = 0: Only between-domain (should maximize cross-domain alignment)
  result_between_only <- gpca_align.hyperdesign(hd, y = lbl, ncomp = 2, u = 0.0)
  
  # Both should produce valid results
  expect_true(all(is.finite(result_within_only$s)))
  expect_true(all(is.finite(result_between_only$s)))
  
  # Check that between-domain alignment (u=0) produces better cross-domain alignment
  scores_between_d1 <- result_between_only$s[1:40, ]
  scores_between_d2 <- result_between_only$s[41:80, ]
  
  # Compute cross-domain correlation for matched samples
  cross_corr_balanced <- cor(scores_d1[,1], scores_d2[,1])  # u=0.2
  cross_corr_between <- cor(scores_between_d1[,1], scores_between_d2[,1])  # u=0
  
  # With u=0 (between-only), cross-domain correlation should be higher
  expect_gt(abs(cross_corr_between), abs(cross_corr_balanced) - 0.1)  # Allow some tolerance
})

# Test 3: Robustness to edge cases and numerical stability
test_that("gpca_align handles edge cases and maintains numerical stability", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("PRIMME")
  
  # Test 1: Small dataset with more components than samples per domain
  small_data <- create_test_hyperdesign(n_per_domain = 5, n_domains = 2, n_features = 10)
  hd_small <- small_data$hd
  
  # Should handle gracefully even when ncomp is large relative to samples
  expect_no_error({
    result_small <- gpca_align.hyperdesign(hd_small, y = lbl, ncomp = 3)
  })
  expect_true(all(is.finite(result_small$s)))
  
  # Test 2: Domains with no shared labels (no between-domain similarity)
  no_overlap_data <- create_test_hyperdesign(n_per_domain = 20, n_domains = 3, 
                                             n_features = 4, label_overlap = FALSE)
  hd_no_overlap <- no_overlap_data$hd
  
  # Should still work even with no cross-domain label matches
  expect_message({
    result_no_overlap <- gpca_align.hyperdesign(hd_no_overlap, y = lbl, ncomp = 2, u = 0.5)
  }, "M_between has zero norm")  # Should detect no between-domain similarities
  
  expect_true(all(is.finite(result_no_overlap$s)))
  
  # Test 3: Single domain (edge case)
  single_domain <- list(domain1 = multidesign::multidesign(
    matrix(rnorm(30 * 4), 30, 4),
    data.frame(lbl = factor(rep(c("A", "B", "C"), each = 10)))
  ))
  hd_single <- multidesign::hyperdesign(single_domain)
  
  # Should work with single domain (only within-domain similarity)
  expect_no_error({
    result_single <- gpca_align.hyperdesign(hd_single, y = lbl, ncomp = 2)
  })
  expect_true(all(is.finite(result_single$s)))
  
  # Test 4: Sparse similarity matrix with isolated samples
  sparse_simfun <- function(labels) {
    # Create very sparse similarity (some samples have no connections)
    n <- length(labels)
    S <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    
    # Only connect first few samples
    for (i in 1:min(5, n)) {
      for (j in 1:min(5, n)) {
        if (labels[i] == labels[j]) {
          S[i, j] <- 1
        }
      }
    }
    # Ensure diagonal is 1
    diag(S) <- 1
    S
  }
  
  test_data_sparse <- create_test_hyperdesign(n_per_domain = 20, n_domains = 2)
  hd_sparse <- test_data_sparse$hd
  
  # Should handle sparse similarities with regularization
  expect_no_error({
    result_sparse <- gpca_align.hyperdesign(hd_sparse, y = lbl, ncomp = 2, 
                                            simfun = sparse_simfun, lambda = 0.5)
  })
  expect_true(all(is.finite(result_sparse$s)))
  
  # Test 5: High-dimensional data (p >> n)
  high_dim_data <- create_test_hyperdesign(n_per_domain = 10, n_domains = 2, n_features = 50)
  hd_high_dim <- high_dim_data$hd
  
  # Should handle high-dimensional data
  expect_no_error({
    result_high_dim <- gpca_align.hyperdesign(hd_high_dim, y = lbl, ncomp = 2)
  })
  expect_true(all(is.finite(result_high_dim$s)))
  
  # Test 6: Perfect correlation between domains
  n <- 25
  X_base <- matrix(rnorm(n * 3), n, 3)
  labels_perfect <- rep(c("A", "B", "C", "D", "E"), each = 5)
  
  # Create perfectly correlated domains (just copies)
  md1_perfect <- multidesign::multidesign(X_base, data.frame(lbl = factor(labels_perfect)))
  md2_perfect <- multidesign::multidesign(X_base + 0.001, data.frame(lbl = factor(labels_perfect)))  # Tiny noise
  
  hd_perfect <- multidesign::hyperdesign(list(domain1 = md1_perfect, domain2 = md2_perfect))
  
  # Should handle near-perfect correlation
  expect_no_error({
    result_perfect <- gpca_align.hyperdesign(hd_perfect, y = lbl, ncomp = 2)
  })
  expect_true(all(is.finite(result_perfect$s)))
  
  # Test 7: Custom preprocessing
  custom_preproc <- multivarious::center() %>% multivarious::standardize()
  
  expect_no_error({
    result_preproc <- gpca_align.hyperdesign(hd_small, y = lbl, ncomp = 2, 
                                             preproc = custom_preproc)
  })
  expect_true(all(is.finite(result_preproc$s)))
  
  # Test 8: Verify PSD correction messages
  # Create a similarity function that might produce non-PSD matrix
  risky_simfun <- function(labels) {
    n <- length(labels)
    S <- matrix(0.5, n, n)  # All 0.5 except diagonal
    diag(S) <- 1
    # This creates a matrix that might have negative eigenvalues
    Matrix::Matrix(S - 0.6 * diag(n) + diag(n), sparse = TRUE)
  }
  
  expect_message({
    result_psd <- gpca_align.hyperdesign(hd_small, y = lbl, ncomp = 2, 
                                         simfun = risky_simfun, lambda = 0.001)
  }, "Applied PSD correction")
  
  expect_true(all(is.finite(result_psd$s)))
})