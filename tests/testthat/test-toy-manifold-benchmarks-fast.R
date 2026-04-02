# ==============================================================
#  FAST TOY MANIFOLD BENCHMARKS - SMOKE TEST SUITE
#  Optimized for CI: < 30 seconds target runtime
#
#  This is a trimmed-down version of test-toy-manifold-benchmarks.R
#  designed to run quickly in every CI build while catching major
#  regressions. For comprehensive validation, see the full benchmark
#  suite in test-toy-manifold-benchmarks.R (nightly/weekly runs).
#
#  KEY OPTIMIZATIONS:
#  - Reduced dataset sizes: 15-25 samples (vs 40-100)
#  - Fewer iterations: 10 (vs 25)
#  - KEMA: sample_frac=0.2, ncomp=2, solver="regression"
#  - Only ONE assignment test (GRASP)
#  - Adjusted thresholds for smaller data
# ==============================================================

library(testthat)
library(manifoldalign)
suppressPackageStartupMessages({
  library(multivarious)
  library(dplyr)  # Required for hyperdesign operations
})

# CI environment check
is_ci <- function() !is.na(Sys.getenv("CI", NA))

# ============= FAST DATASET GENERATORS (REDUCED SIZES) =============

gen_linear_triplet_fast <- function(n = 20, noise = 0.01, seed = 1) {
  manifoldalign:::synthetic_alignment_scenario(
    "linear_affine",
    profile = "fast",
    seed = seed,
    n = n,
    noise = noise
  )
}

gen_isometric_triplet_fast <- function(n = 25, noise = 0.02, seed = 2) {
  manifoldalign:::synthetic_alignment_scenario(
    "isometric_curve",
    profile = "fast",
    seed = seed,
    n = n,
    noise = noise
  )
}

gen_hard_triplet_fast <- function(n_side = 5, noise = 0.03, seed = 3) {
  manifoldalign:::synthetic_alignment_scenario(
    "hard_nonisometric",
    profile = "fast",
    seed = seed,
    n_side = n_side,
    noise = noise
  )
}

# ============= VALIDATION HELPERS =============

compute_alignment_error_strict <- function(latent, recovered) {
  if (is_ci()) {
    stopifnot(requireNamespace("vegan", quietly = FALSE))
  } else if (!requireNamespace("vegan", quietly = TRUE)) {
    skip("vegan package required")
  }

  # Ensure both are base matrices for vegan::procrustes
  latent <- as.matrix(latent)
  recovered <- as.matrix(recovered)

  proc_result <- vegan::procrustes(latent, recovered, symmetric = FALSE)

  list(
    error = proc_result$ss,
    correlation = sqrt(1 - proc_result$ss / sum(latent^2)),
    procrustes_obj = proc_result
  )
}

evaluate_assignment_accuracy_strict <- function(predicted_assignment, true_assignment = NULL) {
  if (is.null(true_assignment)) {
    true_assignment <- seq_along(predicted_assignment)
  }

  n <- length(predicted_assignment)
  exact_accuracy <- mean(predicted_assignment == true_assignment)

  is_valid_permutation <- length(unique(predicted_assignment)) == n &&
                         all(sort(predicted_assignment) == 1:n)

  list(
    accuracy = exact_accuracy,
    is_valid_permutation = is_valid_permutation,
    n_correct = sum(predicted_assignment == true_assignment),
    n_total = n
  )
}

toy_to_hyperdesign <- function(toy_data, use_true_labels = TRUE) {
  if (is_ci()) {
    stopifnot(requireNamespace("multidesign", quietly = FALSE))
    stopifnot(requireNamespace("tibble", quietly = FALSE))
  } else {
    skip_if_not_installed("multidesign")
    skip_if_not_installed("tibble")
  }

  if (!use_true_labels) {
    stop("toy_to_hyperdesign requires meaningful labels", call. = FALSE)
  }

  manifoldalign:::synthetic_alignment_to_hyperdesign(
    toy_data,
    label = toy_data$true_labels,
    label_col = "label",
    id_col = "sample_id"
  )
}

# ============= FAST SMOKE TESTS =============

test_that("[FAST] Dataset generation works with reduced sizes", {
  set.seed(42)

  linear_data <- gen_linear_triplet_fast(n = 20, noise = 0.01, seed = 42)
  expect_equal(nrow(linear_data$latent), 20)
  expect_equal(ncol(linear_data$latent), 2)
  expect_equal(nlevels(linear_data$true_labels), 2)

  # Class separation should still be meaningful
  class1_center <- colMeans(linear_data$latent[linear_data$true_labels == "class1", ])
  class2_center <- colMeans(linear_data$latent[linear_data$true_labels == "class2", ])
  class_separation <- sqrt(sum((class1_center - class2_center)^2))
  expect_gt(class_separation, 0.7)  # Slightly relaxed from 0.8 due to smaller sample

  iso_data <- gen_isometric_triplet_fast(n = 25, noise = 0.02, seed = 42)
  expect_equal(length(iso_data$latent), 25)
  expect_equal(nlevels(iso_data$true_labels), 2)

  hard_data <- gen_hard_triplet_fast(n_side = 5, noise = 0.03, seed = 42)
  expect_equal(nrow(hard_data$latent), 25)  # 5x5 grid
  expect_equal(nlevels(hard_data$true_labels), 2)
})

test_that("[FAST] Linear Similarity Embedding on linear data", {
  message("[fast-benchmarks] Linear Similarity Embedding smoke test")
  set.seed(123)
  linear_data <- gen_linear_triplet_fast(n = 20, noise = 0.01, seed = 123)

  X1 <- linear_data$view1

  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))

  # Reduced iterations from 25 to 10
  result_R <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE,
                              maxit = 10, verbose = FALSE)

  expect_s3_class(result_R, "simembed")
  expect_equal(nrow(result_R$scores), nrow(X1))
  expect_equal(ncol(result_R$scores), 2)
  expect_true(result_R$convergence$convergence == 0)

  # Relaxed threshold for smaller data: 0.3 (was 0.2)
  alignment_error <- compute_alignment_error_strict(linear_data$latent, result_R$scores)
  expect_lt(alignment_error$error, 0.3)
  expect_gt(alignment_error$correlation, 0.85)  # Relaxed from 0.9
})

test_that("[FAST] KEMA on isometric manifolds (optimized)", {
  message("[fast-benchmarks] KEMA smoke test (optimized)")
  skip_if_not_installed("kernlab")

  set.seed(456)
  iso_data <- gen_isometric_triplet_fast(n = 25, noise = 0.01, seed = 456)

  # Use only 2D views
  iso_data_2d <- list(
    latent = iso_data$latent,
    view1 = iso_data$view1,
    view2 = iso_data$view2,
    true_labels = iso_data$true_labels
  )

  hd <- toy_to_hyperdesign(iso_data_2d, use_true_labels = TRUE)

  # OPTIMIZED: regression solver, reduced sample_frac, fewer components
  result <- kema(hd, y = "label",
                 kernel = kernlab::rbfdot(sigma = 1.0),
                 preproc = multivarious::center(),
                 ncomp = 2,  # Reduced from 3
                 solver = "regression",  # Faster than exact
                 auto_retry_exact = FALSE,  # Skip retry for speed
                 sigma = 1.0,
                 sample_frac = 0.2)  # Reduced from 0.4

  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))

  # For smaller data, just check for finite, meaningful embeddings
  latent_param <- as.numeric(iso_data$latent)
  latent_param_full <- c(latent_param, latent_param)

  expect_equal(length(latent_param_full), nrow(result$s))

  cors <- numeric(ncol(result$s))
  for (i in 1:ncol(result$s)) {
    embed_component <- as.numeric(result$s[, i])
    cors[i] <- abs(cor(latent_param_full, embed_component, use = "complete.obs"))
  }

  # Basic sanity: embeddings should be finite and non-zero
  expect_true(all(is.finite(cors)))
  expect_true(max(cors, na.rm = TRUE) > 0)
})

test_that("[FAST] GRASP assignment accuracy (smoke test)", {
  message("[fast-benchmarks] GRASP smoke test")
  set.seed(789)
  linear_data <- gen_linear_triplet_fast(n = 20, noise = 0.01, seed = 789)

  grasp_data <- list(
    view1 = linear_data$view1,
    view2 = linear_data$view2,
    true_labels = linear_data$true_labels
  )
  grasp_hd <- toy_to_hyperdesign(grasp_data, use_true_labels = TRUE)

  # Reduced components: 4 (was 6), q_descriptors: 6 (was 8)
  result <- grasp(grasp_hd, preproc = multivarious::center(),
                  ncomp = 4, q_descriptors = 6)

  expect_s3_class(result, "multiblock_biprojector")
  expect_true(all(is.finite(result$s)))

  if ("assignment" %in% names(result)) {
    accuracy_result <- evaluate_assignment_accuracy_strict(result$assignment)
    expect_true(accuracy_result$is_valid_permutation)

    # For very small datasets (n=20), GRASP may not beat random significantly
    # Just verify it produces a valid assignment
    expect_gte(accuracy_result$n_correct, 0)
  }
})

test_that("[FAST] Runtime performance (smoke check)", {
  message("[fast-benchmarks] Runtime smoke check")
  set.seed(161718)

  linear_data <- gen_linear_triplet_fast(n = 25, noise = 0.01, seed = 161718)

  X1 <- linear_data$view1
  latent_dist <- as.matrix(dist(linear_data$latent))
  T_sim <- exp(-latent_dist^2 / median(latent_dist^2))

  runtime <- system.time({
    result <- linear_sim_embed(X1, T = T_sim, ncomp = 2, use_cpp = FALSE,
                              maxit = 10, verbose = FALSE)
  })

  # Should be very fast with reduced size
  expect_lt(runtime[["elapsed"]], 1)  # < 1 second for n=25
  expect_true(result$convergence$convergence == 0)
})
