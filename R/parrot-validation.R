#' Validation Suite for PARROT Implementation
#'
#' Comprehensive validation functions to test PARROT accuracy and performance
#' against synthetic datasets with known ground truth.
#'
#' @name parrot-validation
#' @keywords internal
NULL

#' Generate Synthetic Network Alignment Data
#' 
#' Creates two networks with known correspondence for validating PARROT
#' 
#' @param n_nodes Number of nodes per network
#' @param n_anchors Number of anchor correspondences
#' @param noise_level Noise level for second network (0-1)
#' @param structure Type of network structure: "ring", "grid", "random", "community"
#' @param permute_fraction Fraction of nodes to permute in second network
#' @return List with two networks and ground truth alignment
#' @export
generate_parrot_validation_data <- function(n_nodes = 100, 
                                          n_anchors = 10,
                                          noise_level = 0.1,
                                          structure = c("ring", "grid", "random", "community"),
                                          permute_fraction = 0.3) {
  
  structure <- match.arg(structure)
  
  # Generate base network structure
  if (structure == "ring") {
    # Ring/circular structure
    angles <- 2 * pi * (1:n_nodes) / n_nodes
    X1 <- cbind(
      cos(angles),
      sin(angles),
      runif(n_nodes, -0.1, 0.1)  # Small z-variation
    )
  } else if (structure == "grid") {
    # Grid structure
    grid_size <- ceiling(sqrt(n_nodes))
    grid_coords <- expand.grid(
      x = 1:grid_size / grid_size,
      y = 1:grid_size / grid_size
    )[1:n_nodes, ]
    X1 <- cbind(
      grid_coords$x,
      grid_coords$y,
      runif(n_nodes, -0.1, 0.1)
    )
  } else if (structure == "random") {
    # Random points
    X1 <- matrix(rnorm(n_nodes * 3), n_nodes, 3)
  } else {  # community
    # Community structure
    n_communities <- 4
    community_size <- n_nodes / n_communities
    X1 <- matrix(0, n_nodes, 3)
    
    for (i in 1:n_communities) {
      idx <- ((i-1)*community_size + 1):(i*community_size)
      center <- c(cos(2*pi*i/n_communities), sin(2*pi*i/n_communities), 0)
      X1[idx, ] <- matrix(rnorm(length(idx) * 3, sd = 0.2), length(idx), 3) + 
                   matrix(center, length(idx), 3, byrow = TRUE)
    }
  }
  
  # Create permutation for second network
  n_permute <- floor(n_nodes * permute_fraction)
  permute_idx <- sample(n_nodes, n_permute)
  perm <- 1:n_nodes
  if (n_permute > 1) {
    perm[permute_idx] <- permute_idx[sample(n_permute)]
  }
  
  # Create second network with noise and permutation
  X2 <- X1[perm, ] + matrix(rnorm(n_nodes * 3, sd = noise_level), n_nodes, 3)
  
  # Select anchor nodes (ensuring they're not in permuted set)
  available_anchors <- setdiff(1:n_nodes, permute_idx)
  if (length(available_anchors) < n_anchors) {
    available_anchors <- 1:n_nodes  # Use all if not enough unpermuted
  }
  anchor_idx <- sort(sample(available_anchors, min(n_anchors, length(available_anchors))))
  
  # Create anchor vectors
  anchors1 <- rep(NA, n_nodes)
  anchors2 <- rep(NA, n_nodes)
  anchors1[anchor_idx] <- anchor_idx
  anchors2[perm[anchor_idx]] <- anchor_idx
  
  # Create hyperdesign structure
  domain1 <- list(
    x = X1,
    design = data.frame(
      node_id = 1:n_nodes,
      anchors = anchors1,
      true_labels = 1:n_nodes  # Ground truth node IDs
    )
  )
  
  domain2 <- list(
    x = X2,
    design = data.frame(
      node_id = 1:n_nodes,
      anchors = anchors2,
      true_labels = perm  # Permuted ground truth
    )
  )
  
  hd <- list(domain1 = domain1, domain2 = domain2)
  class(hd) <- c("hyperdesign", "list")
  
  list(
    data = hd,
    ground_truth = list(
      permutation = perm,
      inverse_perm = order(perm),
      anchor_indices = list(net1 = anchor_idx, net2 = perm[anchor_idx])
    ),
    params = list(
      n_nodes = n_nodes,
      n_anchors = n_anchors,
      noise_level = noise_level,
      structure = structure,
      permute_fraction = permute_fraction
    )
  )
}

#' Evaluate PARROT Alignment Accuracy
#' 
#' Computes various accuracy metrics for PARROT alignment results
#' 
#' @param parrot_result Result from parrot() function
#' @param ground_truth Ground truth from generate_parrot_validation_data
#' @param k Top-k accuracy levels to compute (default: c(1, 5, 10))
#' @return List of accuracy metrics
#' @export
evaluate_parrot_accuracy <- function(parrot_result, ground_truth, k = c(1, 5, 10)) {
  
  transport_plan <- parrot_result$transport_plan
  n1 <- nrow(transport_plan)
  n2 <- ncol(transport_plan)
  
  # Get predicted alignments (top-1 for each node in network 1)
  predicted_alignment <- apply(transport_plan, 1, which.max)
  
  # True alignment (inverse permutation maps network 2 back to network 1)
  true_alignment <- ground_truth$inverse_perm
  
  # Top-1 accuracy
  top1_accuracy <- mean(predicted_alignment == true_alignment)
  
  # Top-k accuracies
  topk_accuracies <- sapply(k, function(k_val) {
    # For each node, check if true match is in top-k predictions
    correct <- sapply(1:n1, function(i) {
      top_k_idx <- order(transport_plan[i, ], decreasing = TRUE)[1:min(k_val, n2)]
      true_alignment[i] %in% top_k_idx
    })
    mean(correct)
  })
  names(topk_accuracies) <- paste0("top", k)
  
  # Accuracy on anchor nodes only
  anchor_idx1 <- ground_truth$anchor_indices$net1
  anchor_accuracy <- mean(predicted_alignment[anchor_idx1] == true_alignment[anchor_idx1])
  
  # Accuracy on non-anchor nodes
  non_anchor_idx <- setdiff(1:n1, anchor_idx1)
  non_anchor_accuracy <- if (length(non_anchor_idx) > 0) {
    mean(predicted_alignment[non_anchor_idx] == true_alignment[non_anchor_idx])
  } else {
    NA
  }
  
  # Mean reciprocal rank (MRR)
  mrr <- mean(sapply(1:n1, function(i) {
    rank_of_true <- which(order(transport_plan[i, ], decreasing = TRUE) == true_alignment[i])
    1 / rank_of_true
  }))
  
  # Transport plan statistics
  transport_stats <- list(
    sparsity = mean(transport_plan < 1e-6),
    max_value = max(transport_plan),
    entropy = -sum(transport_plan * log(transport_plan + 1e-16)),
    row_sum_deviation = sd(rowSums(transport_plan)),
    col_sum_deviation = sd(colSums(transport_plan))
  )
  
  list(
    top1_accuracy = top1_accuracy,
    topk_accuracies = topk_accuracies,
    anchor_accuracy = anchor_accuracy,
    non_anchor_accuracy = non_anchor_accuracy,
    mean_reciprocal_rank = mrr,
    transport_stats = transport_stats,
    predicted_alignment = predicted_alignment,
    confusion_matrix = table(predicted = predicted_alignment, true = true_alignment)
  )
}

#' Run PARROT Validation Suite
#' 
#' Comprehensive validation across different network structures and parameters
#' 
#' @param output_dir Directory to save validation results
#' @param n_replications Number of replications per configuration
#' @return Data frame with validation results
#' @export
run_parrot_validation_suite <- function(output_dir = NULL, n_replications = 5) {
  
  # Validation configurations
  configs <- expand.grid(
    n_nodes = c(50, 100, 200),
    n_anchors = c(5, 10, 20),
    noise_level = c(0.05, 0.1, 0.2),
    structure = c("ring", "grid", "community"),
    permute_fraction = c(0.1, 0.3, 0.5),
    stringsAsFactors = FALSE
  )
  
  # PARROT parameter sets to test
  parrot_params <- list(
    default = list(sigma = 0.15, lambda = 0.1, tau = 0.01, alpha = 0.5, gamma = 0.1),
    high_reg = list(sigma = 0.15, lambda = 0.5, tau = 0.01, alpha = 0.5, gamma = 0.1),
    low_entropy = list(sigma = 0.15, lambda = 0.1, tau = 0.001, alpha = 0.5, gamma = 0.1),
    high_rwr = list(sigma = 0.3, lambda = 0.1, tau = 0.01, alpha = 0.3, gamma = 0.2)
  )
  
  results <- list()
  result_idx <- 1
  
  cat("Running PARROT validation suite...\n")
  pb <- txtProgressBar(min = 0, max = nrow(configs) * length(parrot_params) * n_replications, style = 3)
  
  for (config_idx in 1:nrow(configs)) {
    config <- configs[config_idx, ]
    
    for (param_name in names(parrot_params)) {
      params <- parrot_params[[param_name]]
      
      for (rep in 1:n_replications) {
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
            max_iter = 50,
            tol = 1e-4
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
            top10_accuracy = eval_result$topk_accuracies["top10"],
            anchor_accuracy = eval_result$anchor_accuracy,
            non_anchor_accuracy = eval_result$non_anchor_accuracy,
            mrr = eval_result$mean_reciprocal_rank,
            transport_sparsity = eval_result$transport_stats$sparsity,
            transport_entropy = eval_result$transport_stats$entropy,
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
            top10_accuracy = NA,
            anchor_accuracy = NA,
            non_anchor_accuracy = NA,
            mrr = NA,
            transport_sparsity = NA,
            transport_entropy = NA,
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
  
  # Save if output directory provided
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(results_df, 
              file.path(output_dir, paste0("parrot_validation_", timestamp, ".csv")),
              row.names = FALSE)
    
    # Generate summary report
    generate_parrot_validation_report(results_df, output_dir, timestamp)
  }
  
  results_df
}

#' Generate PARROT Validation Report
#' 
#' Creates a summary report from validation results
#' 
#' @param results_df Data frame from run_parrot_validation_suite
#' @param output_dir Directory to save report
#' @param timestamp Timestamp for file naming
#' @keywords internal
generate_parrot_validation_report <- function(results_df, output_dir, timestamp) {
  
  report_file <- file.path(output_dir, paste0("parrot_validation_report_", timestamp, ".txt"))
  
  sink(report_file)
  
  cat("PARROT VALIDATION REPORT\n")
  cat("========================\n\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Total configurations tested:", nrow(results_df), "\n\n")
  
  # Overall accuracy summary
  cat("OVERALL ACCURACY SUMMARY\n")
  cat("------------------------\n")
  cat(sprintf("Mean Top-1 Accuracy: %.3f (SD: %.3f)\n", 
              mean(results_df$top1_accuracy, na.rm = TRUE),
              sd(results_df$top1_accuracy, na.rm = TRUE)))
  cat(sprintf("Mean Top-5 Accuracy: %.3f (SD: %.3f)\n",
              mean(results_df$top5_accuracy, na.rm = TRUE),
              sd(results_df$top5_accuracy, na.rm = TRUE)))
  cat(sprintf("Mean MRR: %.3f (SD: %.3f)\n",
              mean(results_df$mrr, na.rm = TRUE),
              sd(results_df$mrr, na.rm = TRUE)))
  cat("\n")
  
  # Accuracy by structure type
  cat("ACCURACY BY NETWORK STRUCTURE\n")
  cat("-----------------------------\n")
  by_structure <- aggregate(
    cbind(top1_accuracy, anchor_accuracy, non_anchor_accuracy) ~ structure,
    data = results_df,
    FUN = function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
  )
  print(by_structure)
  cat("\n")
  
  # Accuracy by noise level
  cat("ACCURACY BY NOISE LEVEL\n")
  cat("-----------------------\n")
  by_noise <- aggregate(
    cbind(top1_accuracy, top5_accuracy) ~ noise_level,
    data = results_df,
    FUN = function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
  )
  print(by_noise)
  cat("\n")
  
  # Parameter set comparison
  cat("ACCURACY BY PARAMETER SET\n")
  cat("-------------------------\n")
  by_params <- aggregate(
    cbind(top1_accuracy, runtime) ~ param_set,
    data = results_df,
    FUN = function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
  )
  print(by_params)
  cat("\n")
  
  # Runtime analysis
  cat("RUNTIME ANALYSIS\n")
  cat("----------------\n")
  cat(sprintf("Mean runtime: %.2f seconds\n", mean(results_df$runtime, na.rm = TRUE)))
  runtime_by_size <- aggregate(
    runtime ~ n_nodes,
    data = results_df,
    FUN = function(x) sprintf("%.2f", mean(x, na.rm = TRUE))
  )
  cat("Runtime by network size:\n")
  print(runtime_by_size)
  
  sink()
  
  cat("\nValidation report saved to:", report_file, "\n")
}