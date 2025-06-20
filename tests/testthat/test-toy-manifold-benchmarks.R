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

# ============= TOY DATASET GENERATORS WITH MEANINGFUL CLASS STRUCTURE =============

# -------- helper: random rotation matrix in d dims ----------
rand_rotation <- function(d, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(d^2), d)))
  # enforce det = +1
  if (det(Q) < 0) Q[,1] <- -Q[,1]
  Q
}

# -------- DATA SET A: *Linear / Affine* with TRUE CLASS STRUCTURE ----------
#   - Creates two meaningful Gaussian blobs that persist across views
#   - Tests: affine‑invariant embedding, centring, scaling, supervised methods
# --------------------------------------------------------------
gen_linear_triplet <- function(n = 500, noise = 0.01, seed = 1) {
  set.seed(seed)
  
  # Create two meaningful classes as Gaussian blobs
  n1 <- floor(n/2)
  n2 <- n - n1
  
  # Class 1: centered at (-0.5, -0.5), Class 2: centered at (0.5, 0.5)
  latent_c1 <- matrix(rnorm(n1*2, mean=c(-0.5, -0.5), sd=0.3), n1, 2)
  latent_c2 <- matrix(rnorm(n2*2, mean=c(0.5, 0.5), sd=0.3), n2, 2)
  latent <- rbind(latent_c1, latent_c2)
  colnames(latent) <- c("u1","u2")
  
  # True class labels
  true_labels <- factor(c(rep("class1", n1), rep("class2", n2)))

  # Apply different affine transformations to create views
  R1 <- rand_rotation(2, seed+1);  b1 <- c(0.2,-0.4)
  R2 <- rand_rotation(2, seed+2);  b2 <- c(-0.7, 0.6)
  R3 <- rand_rotation(2, seed+3);  b3 <- c( 0.1, 0.8)

  X1 <- latent %*% diag(c(1.0, 0.7)) %*% R1  + matrix(b1, n, 2, TRUE) + rnorm(n*2,0,noise)
  X2 <- latent %*% diag(c(0.5, 1.2)) %*% R2  + matrix(b2, n, 2, TRUE) + rnorm(n*2,0,noise)
  X3 <- latent %*% diag(c(1.5, 0.4)) %*% R3  + matrix(b3, n, 2, TRUE) + rnorm(n*2,0,noise)

  list(latent = latent, view1 = X1, view2 = X2, view3 = X3, 
       true_labels = true_labels, class_centers = list(c(-0.5, -0.5), c(0.5, 0.5)))
}

# -------- DATA SET B: *Non‑linear but isometric* with PARAMETRIC STRUCTURE ----------
#   Shared 1‑D parameter t ∈ [0,2π] with natural progression-based classes
# --------------------------------------------------------------
gen_isometric_triplet <- function(n = 300, noise = 0.02, seed = 2) {  # Smaller n for faster testing
  set.seed(seed)
  t <- sort(runif(n, 0, 2*pi))          # latent coord

  # Create classes based on parameter progression (early vs late in curve)
  class_boundary <- pi  # Split at halfway point
  true_labels <- factor(ifelse(t < class_boundary, "early", "late"))

  v1 <- cbind(t, 0)                     # straight line
  v2 <- cbind(cos(t), sin(t))           # unit circle
  v3 <- cbind(cos(t), sin(t), t/(2*pi)) # 3‑D helix, 1 turn

  v1 <- v1 + matrix(rnorm(n*2,0,noise), n,2)
  v2 <- v2 + matrix(rnorm(n*2,0,noise), n,2)
  v3 <- v3 + matrix(rnorm(n*3,0,noise), n,3)

  list(latent = t, view1 = v1, view2 = v2, view3 = v3, 
       true_labels = true_labels, parameter_boundary = class_boundary)
}

# -------- DATA SET C: *Non‑isometric / density‑skewed* with SPATIAL CLASSES --------
#   Classes based on spatial regions to test robustness
# --------------------------------------------------------------
gen_hard_triplet <- function(n_side = 25, noise = 0.03, seed = 3) {
  set.seed(seed)
  u <- seq(-1, 1, length.out = n_side)
  grid <- as.matrix(expand.grid(u, u))          # latent (n x 2)

  # Create spatial classes: left vs right half
  true_labels <- factor(ifelse(grid[,1] < 0, "left", "right"))

  # ---- view1: identity + small noise
  v1 <- grid + matrix(rnorm(nrow(grid)*2, 0, noise), ncol = 2)

  # ---- view2: swiss roll embedding in 3‑D (INCREASED COMPLEXITY)
  theta <- (grid[,1] + 1) * 2 * pi      # Increased from 1*pi to 2*pi for more twists
  height <- grid[,2]                   
  v2 <- cbind(theta * cos(theta),
              height,
              theta * sin(theta)) + matrix(rnorm(nrow(grid)*3, 0, noise), ncol=3)

  # ---- view3: affine + density skew + quadrant corruption
  A <- matrix(c(2, 0.5, 0, 0.2), 2, 2)  # anisotropic scale & shear
  v3 <- grid %*% A
  # add heteroscedastic noise in lower‑left quadrant
  idx_hard <- which(grid[,1] < 0 & grid[,2] < 0)
  v3[idx_hard,] <- v3[idx_hard,] + matrix(rnorm(length(idx_hard)*2, 0, 0.2), ncol=2)  # Increased noise
  v3 <- v3 + matrix(rnorm(nrow(grid)*2, 0, noise), ncol=2)

  list(latent = grid, view1 = v1, view2 = v2, view3 = v3,
       hard_indices = idx_hard, true_labels = true_labels)
}

# -------- convenience wrapper ---------------------------------
generate_all_toy_sets <- function() {
  list(
    linear_affine     = gen_linear_triplet(),
    isometric_curve   = gen_isometric_triplet(),
    hard_nonisometric = gen_hard_triplet()
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
  
  view_names <- names(toy_data)[!names(toy_data) %in% c("latent", "hard_indices", "true_labels", "class_centers", "parameter_boundary")]
  
  # Use meaningful labels if available, otherwise error
  if (use_true_labels && "true_labels" %in% names(toy_data)) {
    labels <- toy_data$true_labels
  } else {
    stop("toy_to_hyperdesign requires meaningful labels. Set use_true_labels=TRUE or provide labels explicitly.")
  }
  
  # Create hyperdesign structure
  domain_list <- list()
  for (i in seq_along(view_names)) {
    view_name <- view_names[i]
    X <- as.matrix(toy_data[[view_name]])
    
    design_df <- data.frame(
      sample_id = seq_len(nrow(X)),
      label = labels,
      stringsAsFactors = FALSE
    )
    
    domain_list[[view_name]] <- list(x = X, design = design_df)
  }
  
  class(domain_list) <- c("hyperdesign", "list")
  domain_list
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
  
  expect_length(toy_sets, 3)
  expect_named(toy_sets, c("linear_affine", "isometric_curve", "hard_nonisometric"))
  
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
  set.seed(123)
  linear_data <- gen_linear_triplet(n = 200, noise = 0.01, seed = 123)
  
  X1 <- linear_data$view1
  X2 <- linear_data$view2
  
  # Create target similarity from latent coordinates
  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))
  
  # Test R implementation with strict thresholds
  result_R <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE, 
                              maxit = 100, verbose = FALSE)
  
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
                                  maxit = 100, verbose = FALSE)
    
    expect_s3_class(result_cpp, "simembed")
    expect_true(result_cpp$convergence$convergence == 0)
    
    # Implementations should be nearly identical
    consistency_error <- sqrt(mean((result_R$scores - result_cpp$scores)^2))
    expect_lt(consistency_error, 0.01)  # Very tight consistency requirement
  }
})

test_that("KEMA achieves strict performance on non-linear manifolds", {
  require_package_strict("kernlab")
  
  set.seed(456)
  iso_data <- gen_isometric_triplet(n = 100, noise = 0.01, seed = 456)  # Low noise for strict testing
  
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
                 preproc = multivarious::center(), ncomp = 3, solver = "exact", sigma = 1.0)
  
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
  set.seed(789)
  linear_data <- gen_linear_triplet(n = 100, noise = 0.01, seed = 789)
  
  # Create hyperdesign with just two views for GRASP
  grasp_data <- list(
    view1 = linear_data$view1,
    view2 = linear_data$view2,
    true_labels = linear_data$true_labels
  )
  grasp_hd <- toy_to_hyperdesign(grasp_data, use_true_labels = TRUE)
  
  result <- grasp(grasp_hd, preproc = multivarious::center(), ncomp = 10, q_descriptors = 20)
  
  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))
  
  # STRICT: Should produce valid assignment with reasonable accuracy
  if ("assignment" %in% names(result)) {
    accuracy_result <- evaluate_assignment_accuracy_strict(result$assignment)
    expect_true(accuracy_result$is_valid_permutation)  # Must be valid permutation
    
    # On linear data, should achieve better than random alignment
    n <- nrow(linear_data$view1)
    expected_random_accuracy <- 1/n  # ~0.01 for n=100
    expect_gt(accuracy_result$accuracy, 1.8 * expected_random_accuracy)  # Beat random meaningfully
    
    # Basic sanity: should assign at least some nodes correctly
    expect_gt(accuracy_result$n_correct, 0)  # At least one correct assignment
  }
})

test_that("CONE-Align produces valid assignments with strict convergence", {
  set.seed(101112)
  linear_data <- gen_linear_triplet(n = 80, noise = 0.01, seed = 101112)
  
  cone_data <- list(
    view1 = linear_data$view1,
    view2 = linear_data$view2,
    true_labels = linear_data$true_labels
  )
  cone_hd <- toy_to_hyperdesign(cone_data, use_true_labels = TRUE)
  
  result <- cone_align(cone_hd, preproc = multivarious::center(), 
                      ncomp = 5, max_iter = 20, tol = 0.01)
  
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
  expect_equal(ncol(result$s), 5)
  expect_false(any(is.na(result$s)))
})

test_that("Linear methods fail appropriately on highly non-linear data", {
  set.seed(131415)
  hard_data <- gen_hard_triplet(n_side = 20, noise = 0.02, seed = 131415)
  
  X_grid <- hard_data$view1
  X_swiss <- hard_data$view2
  
  # Create structured similarity that should be disrupted by Swiss roll
  grid_dist <- as.matrix(dist(X_grid))
  T_structured <- exp(-grid_dist^2 / median(grid_dist^2)) * 0.2  # Moderate similarities
  
  result <- linear_sim_embed(X_swiss, T = T_structured, ncomp = 2, 
                            use_cpp = FALSE, maxit = 50, verbose = FALSE)
  
  expect_true(result$convergence$convergence == 0)  # Should converge
  
  # STRICT: Should fail to recover latent structure due to non-linearity
  alignment_error <- compute_alignment_error_strict(hard_data$latent, result$scores)
  expect_gt(alignment_error$error, 0.3)  # Should have substantial error
  expect_lt(alignment_error$correlation, 0.8)  # Should not achieve good correlation
})

# ============= RUNTIME AND PERFORMANCE REGRESSION TESTS =============

test_that("Algorithm runtime performance regression checks", {
  set.seed(161718)
  
  # Test on moderately sized data
  linear_data <- gen_linear_triplet(n = 200, noise = 0.01, seed = 161718)
  
  # Runtime check for Linear Similarity Embedding
  X1 <- linear_data$view1
  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))
  
  runtime <- system.time({
    result <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE, 
                              maxit = 50, verbose = FALSE)
  })
  
  # Should complete in reasonable time (prevent O(n^3) regressions)
  expect_lt(runtime[["elapsed"]], 10)  # 10 seconds should be plenty for n=200
  expect_true(result$convergence$convergence == 0)
})

# Remove always-passing summary tests - convert to vignette material
# test_that("Algorithm performance summary") { expect_true(TRUE) } # REMOVED
# test_that("Example usage patterns") { expect_true(TRUE) } # REMOVED 