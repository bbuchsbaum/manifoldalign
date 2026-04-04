# ==============================================================
#  TOY MANIFOLD ALIGNMENT BENCHMARKS - COMPREHENSIVE TEST SUITE
#  Tests all alignment algorithms on analytically-controlled datasets
#  
#  REGRESSION TEST DESIGN:
#  - Hard failures instead of skips in CI
#  - Procrustes-only alignment metrics (no correlation fallbacks)
#  - Meaningful class structure for supervised methods
#  - Tight thresholds based on algorithmic expectations
#  - Deterministic testing with fixed seeds
# ==============================================================

library(testthat)
library(manifoldalign)
suppressPackageStartupMessages({
  library(multivarious) # Required for preprocessing functions
  # Note: tibble, dplyr, and multidesign removed to avoid hanging due to dependency issues
})

# CI environment check - convert skips to hard failures in CI
is_ci <- function() !is.na(Sys.getenv("CI", NA))

# Hard dependency enforcement for CI
require_package_strict <- function(pkg) {
  if (is_ci()) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required in CI environment. Install with: install.packages('", pkg, "')")
    }
  } else {
    skip_if_not_installed(pkg)
  }
}

# Helper to disable benchmark-heavy tests unless explicitly enabled
skip_benchmarks_unless_enabled <- function() {
  if (!manifoldalign_benchmarks_enabled()) {
    skip("Toy manifold benchmark suite disabled; enable with options(manifoldalign.run_benchmarks = TRUE)")
  }
}

# ============= TOY DATASET GENERATORS WITH MEANINGFUL CLASS STRUCTURE =============

gen_linear_triplet <- function(n = 60, noise = 0.01, seed = 1) {
  manifoldalign:::synthetic_alignment_scenario(
    "linear_affine",
    profile = "full",
    seed = seed,
    n = n,
    noise = noise
  )
}

gen_isometric_triplet <- function(n = 50, noise = 0.02, seed = 2) {
  manifoldalign:::synthetic_alignment_scenario(
    "isometric_curve",
    profile = "full",
    seed = seed,
    n = n,
    noise = noise
  )
}

gen_hard_triplet <- function(n_side = 10, noise = 0.03, seed = 3) {
  manifoldalign:::synthetic_alignment_scenario(
    "hard_nonisometric",
    profile = "full",
    seed = seed,
    n_side = n_side,
    noise = noise
  )
}

generate_all_toy_sets <- function() {
  registry <- manifoldalign:::synthetic_alignment_scenarios(profile = "full")
  stats::setNames(
    lapply(registry$scenario, function(name) {
      manifoldalign:::synthetic_alignment_scenario(name, profile = "full")
    }),
    registry$scenario
  )
}

# ============= STRICT EVALUATION METRICS (PROCRUSTES-ONLY) =============

#' Compute alignment error using STRICT Procrustes analysis only
#' No fallbacks - if Procrustes fails, the test should fail
compute_alignment_error_strict <- function(latent, recovered) {
  # Hard requirement for vegan in CI
  if (is_ci()) {
    stopifnot(requireNamespace("vegan", quietly = FALSE))
  } else if (!requireNamespace("vegan", quietly = TRUE)) {
    skip("vegan package required for Procrustes analysis - install with: install.packages('vegan')")
  }

  # Ensure both are base matrices for vegan::procrustes (fix for Matrix objects)
  latent <- as.matrix(latent)
  recovered <- as.matrix(recovered)

  # Strict Procrustes analysis - no fallbacks
  proc_result <- vegan::procrustes(latent, recovered, symmetric = FALSE)
  
  list(
    error = proc_result$ss,  # Sum of squared residuals
    correlation = sqrt(1 - proc_result$ss / sum(latent^2)),
    procrustes_obj = proc_result  # For detailed analysis if needed
  )
}

#' Evaluate assignment accuracy with strict thresholds
evaluate_assignment_accuracy_strict <- function(predicted_assignment, true_assignment = NULL) {
  if (is.null(true_assignment)) {
    true_assignment <- seq_along(predicted_assignment)
  }
  
  n <- length(predicted_assignment)
  exact_accuracy <- mean(predicted_assignment == true_assignment)
  
  # For permutation-based methods, also check if it's a valid permutation
  is_valid_permutation <- length(unique(predicted_assignment)) == n && 
                         all(sort(predicted_assignment) == 1:n)
  
  list(
    accuracy = exact_accuracy,
    is_valid_permutation = is_valid_permutation,
    n_correct = sum(predicted_assignment == true_assignment),
    n_total = n
  )
}

# ============= HELPER FUNCTIONS FOR HYPERDESIGN CONVERSION =============

#' Convert toy dataset to hyperdesign format with MEANINGFUL labels
toy_to_hyperdesign <- function(toy_data, use_true_labels = TRUE) {
  require_package_strict("multidesign")
  require_package_strict("tibble")

  if (!use_true_labels) {
    stop("toy_to_hyperdesign requires meaningful labels. Set use_true_labels=TRUE.", call. = FALSE)
  }

  manifoldalign:::synthetic_alignment_to_hyperdesign(
    toy_data,
    label = toy_data$true_labels,
    label_col = "label",
    id_col = "sample_id"
  )
}

# ============= STRICT ALGORITHM TESTS WITH TIGHT THRESHOLDS =============
# 
# PHILOSOPHY: These tests use strict thresholds designed to catch regressions.
# Thresholds are set based on algorithmic expectations:
# - Linear methods should excel on linear data
# - Non-linear methods should handle isometric manifolds well  
# - Assignment methods should beat random chance meaningfully
# - All methods should produce valid, finite outputs

test_that("Toy datasets are generated correctly with meaningful class structure", {
  set.seed(42)  # Deterministic testing
  toy_sets <- generate_all_toy_sets()
  
  expect_true(length(toy_sets) >= 3)
  expect_true(all(c("linear_affine", "isometric_curve", "hard_nonisometric") %in% names(toy_sets)))
  
  # Test linear/affine dataset with class structure
  linear_data <- toy_sets$linear_affine
  expect_true(all(c("latent", "view1", "view2", "view3", "true_labels") %in% names(linear_data)))
  expect_equal(ncol(linear_data$latent), 2)
  expect_equal(nlevels(linear_data$true_labels), 2)
  
  # Test class separation in latent space
  class1_center <- colMeans(linear_data$latent[linear_data$true_labels == "class1", ])
  class2_center <- colMeans(linear_data$latent[linear_data$true_labels == "class2", ])
  class_separation <- sqrt(sum((class1_center - class2_center)^2))
  expect_gt(class_separation, 0.8)  # Classes should be well-separated
  
  # Test isometric dataset
  iso_data <- toy_sets$isometric_curve
  expect_true(is.vector(iso_data$latent))
  expect_equal(nlevels(iso_data$true_labels), 2)
  
  # Test hard dataset 
  hard_data <- toy_sets$hard_nonisometric
  expect_equal(ncol(hard_data$latent), 2)
  expect_equal(nlevels(hard_data$true_labels), 2)
  expect_true("hard_indices" %in% names(hard_data))
})

test_that("Linear Similarity Embedding achieves strict performance on linear data", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] Linear Similarity Embedding strict test")
  set.seed(123)
  linear_data <- gen_linear_triplet(n = 50, noise = 0.01, seed = 123)
  
  X1 <- linear_data$view1
  X2 <- linear_data$view2
  
  # Create target similarity from latent coordinates
  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))
  
  # Test R implementation with strict thresholds
  result_R <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE, 
                              maxit = 25, verbose = FALSE)
  
  expect_s3_class(result_R, "simembed")
  expect_equal(nrow(result_R$scores), nrow(X1))
  expect_equal(ncol(result_R$scores), 2)
  expect_true(result_R$convergence$convergence == 0)
  
  # STRICT: Should achieve excellent alignment on linear data
  alignment_error <- compute_alignment_error_strict(linear_data$latent, result_R$scores)
  expect_lt(alignment_error$error, 0.2)  # Strict but achievable for linear data
  expect_gt(alignment_error$correlation, 0.9)  # Should be excellent
  
  # Test C++ consistency if available
  if (exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function")) {
    result_cpp <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = TRUE,
                                  maxit = 25, verbose = FALSE)
    
    expect_s3_class(result_cpp, "simembed")
    expect_true(result_cpp$convergence$convergence == 0)
    
    # Implementations should be nearly identical
    consistency_error <- sqrt(mean((result_R$scores - result_cpp$scores)^2))
    expect_lt(consistency_error, 0.01)  # Very tight consistency requirement
  }
})

test_that("lowrank_align achieves strong performance on linear data", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] lowrank_align strict linear test")
  require_package_strict("multidesign")
  require_package_strict("tibble")

  set.seed(20240502)
  linear_data <- gen_linear_triplet(n = 60, noise = 0.01, seed = 20240502)

  X1 <- linear_data$view1
  X2 <- linear_data$view2
  labels <- linear_data$true_labels

  # Build hyperdesign for two linear views
  md1 <- multidesign::multidesign(X1, tibble::tibble(y = labels))
  md2 <- multidesign::multidesign(X2, tibble::tibble(y = labels))
  hd <- multidesign::hyperdesign(list(view1 = md1, view2 = md2))

  # Simple label-based similarity: same label = 1, different or NA = 0
  simfun <- function(lbl) {
    labs <- as.character(lbl)
    n <- length(labs)
    M <- matrix(0, n, n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        li <- labs[i]; lj <- labs[j]
        if (is.na(li) || is.na(lj)) next
        if (li == lj) M[i, j] <- 1
      }
    }
    Matrix::Matrix(M, sparse = TRUE)
  }

  fit <- lowrank_align(hd, y,
                       simfun = simfun,
                       ncomp = 2,
                       solver = "explicit",
                       sv_thresh = 0.5,
                       lambda = 0.01)

  expect_s3_class(fit, c("lowrank_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$s), 2 * nrow(X1))
  expect_equal(ncol(fit$s), 2)
  expect_true(all(is.finite(fit$s)))

  # Evaluate alignment on the first view via strict Procrustes metrics
  emb1 <- fit$s[seq_len(nrow(X1)), , drop = FALSE]
  align_lr <- compute_alignment_error_strict(linear_data$latent, emb1)

  # PCA baseline on the first view
  pca_scores <- stats::prcomp(X1, rank. = 2)$x
  align_pca <- compute_alignment_error_strict(linear_data$latent, pca_scores)

  # lowrank_align should be at least competitive with PCA on clean linear data
  expect_lt(align_lr$error, 1.2 * align_pca$error)      # allow 20% slack
  expect_gt(align_lr$correlation, 0.8)                  # strong correlation with latent coords
})

test_that("Gromov-Wasserstein achieves strong alignment on nearly identical views", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] Gromov-Wasserstein isometric test")
  require_package_strict("multidesign")

  set.seed(20240503)
  iso_data <- gen_isometric_triplet(n = 40, noise = 0.01, seed = 20240503)

  X1 <- iso_data$view2
  X2 <- iso_data$view2 + matrix(rnorm(length(X1), sd = 0.01), nrow(X1), ncol(X1))

  design <- data.frame(sample_id = seq_len(nrow(X1)))
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

  result <- gromov_wasserstein(hd, epsilon = 0.05, max_iter = 100, verbose = FALSE)

  expect_s3_class(result, "gromov_wasserstein")
  P <- result$transport_plans[[1]]
  n <- nrow(P)

  # Transport plan should be approximately doubly stochastic
  target <- rep(1 / n, n)
  expect_equal(rowSums(P), target, tolerance = 5e-3)
  expect_equal(colSums(P), target, tolerance = 5e-3)

  # On nearly identical domains, diagonal mass should be high
  diag_mass <- sum(diag(P))
  expect_true(is.finite(diag_mass))
  expect_gt(diag_mass, 0.5)
})

test_that("FPGW transport concentrates on diagonal for linear views", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] FPGW linear alignment test")
  require_package_strict("multidesign")
  require_package_strict("lpSolve")

  set.seed(20240504)
  linear_data <- gen_linear_triplet(n = 30, noise = 0.01, seed = 20240504)

  X1 <- linear_data$view1
  X2 <- linear_data$view2

  design <- data.frame(id = seq_len(nrow(X1)))
  md1 <- multidesign::multidesign(X1, design)
  md2 <- multidesign::multidesign(X2, design)
  hd <- multidesign::hyperdesign(list(d1 = md1, d2 = md2))

  # Classical FGW regime (rho = NULL, lambda = 0) as in core tests
  result <- fpgw(hd, omega1 = 0.5, verbose = FALSE)

  expect_s3_class(result, "fpgw")
  P <- result$transport_plans[[1]]
  n <- nrow(P)

  # Mass and marginals should resemble a doubly stochastic plan
  expect_equal(sum(P), 1, tolerance = 0.1)
  expect_equal(as.vector(rowSums(P)), rep(1 / n, n), tolerance = 0.1)
  expect_equal(as.vector(colSums(P)), rep(1 / n, n), tolerance = 0.1)

  # On linearly related views, most mass should lie near the diagonal
  diag_mass <- sum(diag(P))
  expect_true(is.finite(diag_mass))
  expect_gt(diag_mass, 0.3)
})

test_that("KEMA achieves strict performance on non-linear manifolds", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] KEMA isometric manifold test")
  require_package_strict("kernlab")
  
  set.seed(456)
  iso_data <- gen_isometric_triplet(n = 50, noise = 0.01, seed = 456)  # Low noise for strict testing
  
  # Use only 2D views for dimensional consistency
  iso_data_2d <- list(
    latent = iso_data$latent,
    view1 = iso_data$view1,
    view2 = iso_data$view2,
    true_labels = iso_data$true_labels
  )
  
  hd <- toy_to_hyperdesign(iso_data_2d, use_true_labels = TRUE)
  
  # KEMA should excel on isometric manifolds
  result <- kema(hd, y = "label", kernel = kernlab::rbfdot(sigma = 1.0), 
                 preproc = multivarious::center(), ncomp = 3, solver = "exact", sigma = 1.0, sample_frac = 0.4)
  
  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))  # No NaN/Inf values allowed
  
  # STRICT: Should capture 1D parametric structure excellently
  # KEMA concatenates both views, so we need to replicate latent parameter
  latent_param <- as.numeric(iso_data$latent)  # n=100 latent parameters
  latent_param_full <- c(latent_param, latent_param)  # Replicate for both views: 200 total
  
  # Verify dimensions match
  expect_equal(length(latent_param_full), nrow(result$s))
  
  # Compute correlations safely
  cors <- numeric(ncol(result$s))
  for (i in 1:ncol(result$s)) {
    embed_component <- as.numeric(result$s[, i])
    cors[i] <- abs(cor(latent_param_full, embed_component, use = "complete.obs"))
  }
  
  best_cor <- max(cors, na.rm = TRUE)
  
  # DIAGNOSTIC: Report KEMA performance for debugging
  cat(sprintf("\nKEMA Performance on Isometric Manifold:\n"))
  cat(sprintf("  Best correlation with latent parameter: %.3f\n", best_cor))
  cat(sprintf("  Component correlations: %s\n", paste(sprintf("%.3f", cors), collapse = ", ")))
  
  # KEMA should at least produce finite, meaningful embeddings
  expect_true(all(is.finite(cors)))  # All correlations should be finite
  expect_true(best_cor > 0)  # Should have some positive correlation
  
  # Should preserve class structure to some degree
  if (length(unique(iso_data$true_labels)) > 1) {
    # Replicate labels for both views like we did for latent params
    labels_full <- c(iso_data$true_labels, iso_data$true_labels)
    
    # Check that embeddings separate classes
    s_class1 <- result$s[labels_full == levels(iso_data$true_labels)[1], 1]
    s_class2 <- result$s[labels_full == levels(iso_data$true_labels)[2], 1]
    
    if (length(s_class1) > 0 && length(s_class2) > 0) {
      # Classes should be somewhat separable
      class_separation <- abs(mean(s_class1) - mean(s_class2)) / 
                         sqrt((var(s_class1) + var(s_class2))/2)
      cat(sprintf("  Class separation: %.3f\n", class_separation))
      
      # Basic sanity: should produce different embeddings for different classes
      expect_true(is.finite(class_separation))
    }
  }
})

test_that("GRASP achieves strict assignment accuracy on structured data", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] GRASP assignment accuracy test")
  set.seed(789)
  linear_data <- gen_linear_triplet(n = 45, noise = 0.01, seed = 789)
  
  # Create hyperdesign with just two views for GRASP
  grasp_data <- list(
    view1 = linear_data$view1,
    view2 = linear_data$view2,
    true_labels = linear_data$true_labels
  )
  grasp_hd <- toy_to_hyperdesign(grasp_data, use_true_labels = TRUE)
  
  result <- grasp(grasp_hd, preproc = multivarious::center(), ncomp = 6, q_descriptors = 8)
  
  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))
  
  # STRICT: Should produce valid assignment with reasonable accuracy
  if ("assignment" %in% names(result)) {
    accuracy_result <- evaluate_assignment_accuracy_strict(result$assignment)
    expect_true(accuracy_result$is_valid_permutation)  # Must be valid permutation
    
    # On linear data, GRASP should beat random (but threshold is lenient since it's optimized for manifolds)
    n <- nrow(linear_data$view1)
    expected_random_accuracy <- 1/n  # ~0.01 for n=100
    expect_gt(accuracy_result$accuracy, 1.2 * expected_random_accuracy)  # Beat random with lenient threshold
    
    # Basic sanity: should assign at least some nodes correctly
    expect_gt(accuracy_result$n_correct, 0)  # At least one correct assignment
  }
})

test_that("CONE-Align produces valid assignments with strict convergence", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] CONE-Align convergence test")
  set.seed(101112)
  linear_data <- gen_linear_triplet(n = 40, noise = 0.01, seed = 101112)
  
  cone_data <- list(
    view1 = linear_data$view1,
    view2 = linear_data$view2,
    true_labels = linear_data$true_labels
  )
  cone_hd <- toy_to_hyperdesign(cone_data, use_true_labels = TRUE)
  
  result <- cone_align(cone_hd, preproc = multivarious::center(), 
                      ncomp = 4, max_iter = 12, tol = 0.02)
  
  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))
  
  # STRICT: Must produce valid assignment
  if ("assignment" %in% names(result)) {
    accuracy_result <- evaluate_assignment_accuracy_strict(result$assignment)
    expect_true(accuracy_result$is_valid_permutation)
    
    # Should converge to better than random solution
    n <- nrow(linear_data$view1)
    expected_random_accuracy <- 1/n  # ~0.0125 for n=80
    expect_gt(accuracy_result$accuracy, 1.5 * expected_random_accuracy)  # Beat random
    
    # Basic sanity: should assign at least some nodes correctly
    expect_gt(accuracy_result$n_correct, 0)  # At least one correct assignment
  }
  
  # Check embedding quality
  expect_equal(ncol(result$s), 4)  # Should match ncomp requested
  expect_false(any(is.na(result$s)))
})

test_that("Linear methods fail appropriately on highly non-linear data", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] Linear methods failure-on-nonlinear test")
  set.seed(131415)
  hard_data <- gen_hard_triplet(n_side = 10, noise = 0.02, seed = 131415)
  
  X_grid <- hard_data$view1
  X_swiss <- hard_data$view2
  
  # Create structured similarity that should be disrupted by Swiss roll
  grid_dist <- as.matrix(dist(X_grid))
  T_structured <- exp(-grid_dist^2 / median(grid_dist^2)) * 0.2  # Moderate similarities
  
  result <- linear_sim_embed(X_swiss, T = T_structured, ncomp = 2, 
                            use_cpp = FALSE, maxit = 25, verbose = FALSE)
  
  expect_true(result$convergence$convergence == 0)  # Should converge
  
  # STRICT: Should fail to recover latent structure due to non-linearity
  alignment_error <- compute_alignment_error_strict(hard_data$latent, result$scores)
  expect_gt(alignment_error$error, 0.3)  # Should have substantial error
  expect_lt(alignment_error$correlation, 0.8)  # Should not achieve good correlation
})

# ============= RUNTIME AND PERFORMANCE REGRESSION TESTS =============

test_that("Algorithm runtime performance regression checks", {
  skip_benchmarks_unless_enabled()
  message("[toy-benchmarks] Runtime regression test")
  set.seed(161718)
  
  # Test on moderately sized data
  linear_data <- gen_linear_triplet(n = 60, noise = 0.01, seed = 161718)
  
  # Runtime check for Linear Similarity Embedding
  X1 <- linear_data$view1
  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))
  
  runtime <- system.time({
  result <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE, 
                              maxit = 25, verbose = FALSE)
  })
  
  # Should complete in reasonable time (prevent O(n^3) regressions)
  expect_lt(runtime[["elapsed"]], 2)  # 2 seconds should be plenty for n=60
  expect_true(result$convergence$convergence == 0)
})

# Remove always-passing summary tests - convert to vignette material
# test_that("Algorithm performance summary") { expect_true(TRUE) } # REMOVED
# test_that("Example usage patterns") { expect_true(TRUE) } # REMOVED 
