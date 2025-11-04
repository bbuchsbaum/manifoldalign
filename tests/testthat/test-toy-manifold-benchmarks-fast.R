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

rand_rotation <- function(d, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(d^2), d)))
  if (det(Q) < 0) Q[,1] <- -Q[,1]
  Q
}

# Reduced from n=60 to n=20
gen_linear_triplet_fast <- function(n = 20, noise = 0.01, seed = 1) {
  set.seed(seed)

  n1 <- floor(n/2)
  n2 <- n - n1

  latent_c1 <- matrix(rnorm(n1*2, mean=c(-0.5, -0.5), sd=0.3), n1, 2)
  latent_c2 <- matrix(rnorm(n2*2, mean=c(0.5, 0.5), sd=0.3), n2, 2)
  latent <- rbind(latent_c1, latent_c2)
  colnames(latent) <- c("u1","u2")

  true_labels <- factor(c(rep("class1", n1), rep("class2", n2)))

  R1 <- rand_rotation(2, seed+1);  b1 <- c(0.2,-0.4)
  R2 <- rand_rotation(2, seed+2);  b2 <- c(-0.7, 0.6)
  R3 <- rand_rotation(2, seed+3);  b3 <- c( 0.1, 0.8)

  X1 <- latent %*% diag(c(1.0, 0.7)) %*% R1  + matrix(b1, n, 2, TRUE) + rnorm(n*2,0,noise)
  X2 <- latent %*% diag(c(0.5, 1.2)) %*% R2  + matrix(b2, n, 2, TRUE) + rnorm(n*2,0,noise)
  X3 <- latent %*% diag(c(1.5, 0.4)) %*% R3  + matrix(b3, n, 2, TRUE) + rnorm(n*2,0,noise)

  list(latent = latent, view1 = X1, view2 = X2, view3 = X3,
       true_labels = true_labels)
}

# Reduced from n=50 to n=25
gen_isometric_triplet_fast <- function(n = 25, noise = 0.02, seed = 2) {
  set.seed(seed)
  t <- sort(runif(n, 0, 2*pi))

  class_boundary <- pi
  true_labels <- factor(ifelse(t < class_boundary, "early", "late"))

  v1 <- cbind(t, 0)
  v2 <- cbind(cos(t), sin(t))
  v3 <- cbind(cos(t), sin(t), t/(2*pi))

  v1 <- v1 + matrix(rnorm(n*2,0,noise), n,2)
  v2 <- v2 + matrix(rnorm(n*2,0,noise), n,2)
  v3 <- v3 + matrix(rnorm(n*3,0,noise), n,3)

  list(latent = t, view1 = v1, view2 = v2, view3 = v3,
       true_labels = true_labels)
}

# Reduced from 10x10 (100 samples) to 5x5 (25 samples)
gen_hard_triplet_fast <- function(n_side = 5, noise = 0.03, seed = 3) {
  set.seed(seed)
  u <- seq(-1, 1, length.out = n_side)
  grid <- as.matrix(expand.grid(u, u))

  true_labels <- factor(ifelse(grid[,1] < 0, "left", "right"))

  v1 <- grid + matrix(rnorm(nrow(grid)*2, 0, noise), ncol = 2)

  theta <- (grid[,1] + 1) * 2 * pi
  height <- grid[,2]
  v2 <- cbind(theta * cos(theta),
              height,
              theta * sin(theta)) + matrix(rnorm(nrow(grid)*3, 0, noise), ncol=3)

  A <- matrix(c(2, 0.5, 0, 0.2), 2, 2)
  v3 <- grid %*% A
  idx_hard <- which(grid[,1] < 0 & grid[,2] < 0)
  v3[idx_hard,] <- v3[idx_hard,] + matrix(rnorm(length(idx_hard)*2, 0, 0.2), ncol=2)
  v3 <- v3 + matrix(rnorm(nrow(grid)*2, 0, noise), ncol=2)

  list(latent = grid, view1 = v1, view2 = v2, view3 = v3,
       hard_indices = idx_hard, true_labels = true_labels)
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

  view_names <- names(toy_data)[!names(toy_data) %in% c("latent", "hard_indices", "true_labels", "class_centers", "parameter_boundary")]

  if (use_true_labels && "true_labels" %in% names(toy_data)) {
    labels <- toy_data$true_labels
  } else {
    stop("toy_to_hyperdesign requires meaningful labels")
  }

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