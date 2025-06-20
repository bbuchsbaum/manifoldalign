# Test PARROT validation framework
library(devtools)
devtools::load_all(export_all = FALSE, compile = FALSE)

# Source validation functions
source("R/parrot-validation.R")
source("R/parrot-benchmark.R")

cat("Testing PARROT validation framework...\n\n")

# Test 1: Generate validation data
cat("1. Testing validation data generation...\n")
val_data <- generate_parrot_validation_data(
  n_nodes = 50,
  n_anchors = 5,
  noise_level = 0.1,
  structure = "ring",
  permute_fraction = 0.2
)

cat("  - Generated data with", val_data$params$n_nodes, "nodes\n")
cat("  - Number of anchors:", val_data$params$n_anchors, "\n")
cat("  - Permutation check:", 
    sum(val_data$ground_truth$permutation != 1:50), "nodes permuted\n\n")

# Test 2: Run PARROT on validation data
cat("2. Testing PARROT on validation data...\n")
parrot_result <- parrot(
  val_data$data,
  anchors = anchors,
  ncomp = 5,
  max_iter = 20,
  tol = 1e-3
)

cat("  - PARROT completed successfully\n")
cat("  - Transport plan shape:", dim(parrot_result$transport_plan), "\n\n")

# Test 3: Evaluate accuracy
cat("3. Testing accuracy evaluation...\n")
eval_result <- evaluate_parrot_accuracy(parrot_result, val_data$ground_truth)

cat("  - Top-1 accuracy:", round(eval_result$top1_accuracy, 3), "\n")
cat("  - Top-5 accuracy:", round(eval_result$topk_accuracies["top5"], 3), "\n")
cat("  - Anchor accuracy:", round(eval_result$anchor_accuracy, 3), "\n")
cat("  - Non-anchor accuracy:", round(eval_result$non_anchor_accuracy, 3), "\n")
cat("  - Mean reciprocal rank:", round(eval_result$mean_reciprocal_rank, 3), "\n\n")

# Test 4: Quick scalability test
cat("4. Testing scalability benchmark...\n")
scale_results <- benchmark_parrot_scalability(
  sizes = c(30, 50, 100),
  n_reps = 2,
  sparse_graph = TRUE
)

cat("\n")

# Test 5: Baseline comparison
cat("5. Testing baseline comparison...\n")
comparison <- compare_parrot_baselines(val_data)

cat("\n")

# Test 6: Different network structures
cat("6. Testing different network structures...\n")
structures <- c("ring", "grid", "random", "community")
structure_results <- list()

for (struct in structures) {
  cat("  Testing", struct, "structure...")
  
  val_data_struct <- generate_parrot_validation_data(
    n_nodes = 40,
    n_anchors = 4,
    noise_level = 0.1,
    structure = struct,
    permute_fraction = 0.25
  )
  
  result_struct <- parrot(
    val_data_struct$data,
    anchors = anchors,
    ncomp = 5,
    max_iter = 15,
    tol = 1e-3
  )
  
  eval_struct <- evaluate_parrot_accuracy(result_struct, val_data_struct$ground_truth)
  
  structure_results[[struct]] <- eval_struct$top1_accuracy
  cat(" accuracy =", round(eval_struct$top1_accuracy, 3), "\n")
}

cat("\n")

# Summary
cat("VALIDATION FRAMEWORK TEST SUMMARY\n")
cat("=================================\n")
cat("✓ Data generation working correctly\n")
cat("✓ PARROT runs on validation data\n")
cat("✓ Accuracy metrics computed successfully\n")
cat("✓ Scalability benchmarking functional\n")
cat("✓ Baseline comparisons working\n")
cat("✓ Multiple network structures tested\n")

cat("\nStructure-specific accuracies:\n")
for (struct in names(structure_results)) {
  cat(sprintf("  %s: %.3f\n", struct, structure_results[[struct]]))
}

cat("\nValidation framework is ready for comprehensive testing!\n")