# Run comprehensive PARROT validation
library(devtools)
devtools::load_all(export_all = FALSE, compile = FALSE)

# Source validation functions
source("R/parrot-validation.R")
source("R/parrot-benchmark.R")

cat("Running comprehensive PARROT validation...\n")
cat("==========================================\n\n")

# Create output directory
output_dir <- "parrot_validation_results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Part 1: Quick validation suite (reduced for demonstration)
cat("Part 1: Running validation suite (reduced version)...\n")

# Simplified configurations for quick run
configs <- expand.grid(
  n_nodes = c(50, 100),
  n_anchors = c(5, 10),
  noise_level = c(0.05, 0.15),
  structure = c("ring", "community"),
  permute_fraction = c(0.2, 0.4),
  stringsAsFactors = FALSE
)

# Just test default and high regularization
parrot_params <- list(
  default = list(sigma = 0.15, lambda = 0.1, tau = 0.01, alpha = 0.5, gamma = 0.1),
  high_reg = list(sigma = 0.15, lambda = 0.5, tau = 0.01, alpha = 0.5, gamma = 0.1)
)

# Run with just 2 replications for speed
results <- list()
result_idx <- 1

pb <- txtProgressBar(min = 0, max = nrow(configs) * length(parrot_params) * 2, style = 3)

for (config_idx in 1:nrow(configs)) {
  config <- configs[config_idx, ]
  
  for (param_name in names(parrot_params)) {
    params <- parrot_params[[param_name]]
    
    for (rep in 1:2) {
      setTxtProgressBar(pb, result_idx)
      
      # Generate validation data
      val_data <- generate_parrot_validation_data(
        n_nodes = config$n_nodes,
        n_anchors = config$n_anchors,
        noise_level = config$noise_level,
        structure = config$structure,
        permute_fraction = config$permute_fraction
      )
      
      # Run PARROT
      start_time <- Sys.time()
      
      tryCatch({
        parrot_result <- parrot(
          val_data$data,
          anchors = anchors,
          sigma = params$sigma,
          lambda = params$lambda,
          tau = params$tau,
          alpha = params$alpha,
          gamma = params$gamma,
          ncomp = 10,
          max_iter = 30,
          tol = 1e-3
        )
        
        runtime <- as.numeric(Sys.time() - start_time, units = "secs")
        
        # Evaluate accuracy
        eval_result <- evaluate_parrot_accuracy(parrot_result, val_data$ground_truth)
        
        # Store results
        results[[result_idx]] <- data.frame(
          config,
          param_set = param_name,
          replication = rep,
          runtime = runtime,
          top1_accuracy = eval_result$top1_accuracy,
          top5_accuracy = eval_result$topk_accuracies["top5"],
          anchor_accuracy = eval_result$anchor_accuracy,
          non_anchor_accuracy = eval_result$non_anchor_accuracy,
          mrr = eval_result$mean_reciprocal_rank,
          stringsAsFactors = FALSE
        )
        
      }, error = function(e) {
        results[[result_idx]] <- data.frame(
          config,
          param_set = param_name,
          replication = rep,
          runtime = NA,
          top1_accuracy = NA,
          top5_accuracy = NA,
          anchor_accuracy = NA,
          non_anchor_accuracy = NA,
          mrr = NA,
          error = as.character(e$message),
          stringsAsFactors = FALSE
        )
      })
      
      result_idx <- result_idx + 1
    }
  }
}

close(pb)

# Combine results
results_df <- do.call(rbind, results)

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write.csv(results_df, 
          file.path(output_dir, paste0("parrot_validation_", timestamp, ".csv")),
          row.names = FALSE)

cat("\n\nPart 2: Performance benchmarking...\n")

# Scalability test
scale_results <- benchmark_parrot_scalability(
  sizes = c(50, 100, 200),
  n_reps = 3
)

write.csv(scale_results,
          file.path(output_dir, paste0("parrot_scalability_", timestamp, ".csv")),
          row.names = FALSE)

# Part 3: Generate summary report
cat("\n\nGenerating summary report...\n")

report_file <- file.path(output_dir, paste0("parrot_validation_summary_", timestamp, ".txt"))
sink(report_file)

cat("PARROT VALIDATION SUMMARY REPORT\n")
cat("================================\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total configurations tested:", nrow(results_df), "\n\n")

# Overall accuracy
cat("OVERALL PERFORMANCE\n")
cat("-------------------\n")
successful <- results_df[!is.na(results_df$top1_accuracy), ]
cat(sprintf("Success rate: %.1f%%\n", 100 * nrow(successful) / nrow(results_df)))
cat(sprintf("Mean Top-1 Accuracy: %.3f (SD: %.3f)\n", 
            mean(successful$top1_accuracy), sd(successful$top1_accuracy)))
cat(sprintf("Mean Top-5 Accuracy: %.3f (SD: %.3f)\n",
            mean(successful$top5_accuracy), sd(successful$top5_accuracy)))
cat(sprintf("Mean Anchor Accuracy: %.3f\n", mean(successful$anchor_accuracy)))
cat(sprintf("Mean Non-Anchor Accuracy: %.3f\n", mean(successful$non_anchor_accuracy)))
cat(sprintf("Mean Runtime: %.2f seconds\n\n", mean(successful$runtime)))

# By structure
cat("ACCURACY BY NETWORK STRUCTURE\n")
cat("-----------------------------\n")
for (struct in unique(successful$structure)) {
  struct_data <- successful[successful$structure == struct, ]
  cat(sprintf("%s: %.3f (n=%d)\n", 
              struct, 
              mean(struct_data$top1_accuracy),
              nrow(struct_data)))
}
cat("\n")

# By noise level
cat("ACCURACY BY NOISE LEVEL\n")
cat("-----------------------\n")
for (noise in sort(unique(successful$noise_level))) {
  noise_data <- successful[successful$noise_level == noise, ]
  cat(sprintf("%.2f: %.3f (n=%d)\n", 
              noise,
              mean(noise_data$top1_accuracy),
              nrow(noise_data)))
}
cat("\n")

# By parameter set
cat("ACCURACY BY PARAMETER SET\n")
cat("-------------------------\n")
for (pset in unique(successful$param_set)) {
  pset_data <- successful[successful$param_set == pset, ]
  cat(sprintf("%s: %.3f (n=%d)\n", 
              pset,
              mean(pset_data$top1_accuracy),
              nrow(pset_data)))
}

# Scalability summary
cat("\n\nSCALABILITY ANALYSIS\n")
cat("--------------------\n")
scale_summary <- aggregate(total_time ~ n_nodes, data = scale_results, 
                          FUN = function(x) c(mean = mean(x), sd = sd(x)))
print(scale_summary)

sink()

cat("\n\nValidation complete!\n")
cat("Results saved to:", output_dir, "\n")
cat("Summary report:", report_file, "\n\n")

# Print key findings
cat("KEY FINDINGS:\n")
cat("=============\n")
cat(sprintf("1. Overall accuracy: %.1f%% top-1, %.1f%% top-5\n",
            100 * mean(successful$top1_accuracy),
            100 * mean(successful$top5_accuracy)))
cat(sprintf("2. Anchor nodes achieve %.1f%% accuracy\n",
            100 * mean(successful$anchor_accuracy)))
cat(sprintf("3. Best structure: %s\n",
            names(which.max(tapply(successful$top1_accuracy, 
                                 successful$structure, mean)))))
cat(sprintf("4. Performance scales approximately O(n^%.1f)\n",
            coef(lm(log(total_time) ~ log(n_nodes), 
                   data = aggregate(total_time ~ n_nodes, scale_results, mean)))[2]))