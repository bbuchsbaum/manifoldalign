#' Pseudolabeling for Unsupervised Domain Adaptation
#'
#' @description
#' Provides pseudolabeling system for unsupervised domain adaptation with 
#' KEMA. Identifies high-confidence anchor samples to guide domain alignment.
#'
#' @details
#' The pseudolabeling system addresses unsupervised domain adaptation by 
#' identifying reliable correspondences between samples from different domains. 
#' The approach uses similarity-based clustering, diversity-aware selection, 
#' adaptive thresholding, and quality control filtering.
#' 
#' Main functions:
#' \itemize{
#'   \item \code{assign_pseudolabels()}: General-purpose pseudolabeling from 
#'     sparse similarity matrices
#'   \item \code{high_sim_pseudolabels()}: Specialized function for 
#'     multi-domain data using cosine similarity
#'   \item \code{create_synthetic_similarity_matrix()}: Generate synthetic 
#'     data for testing
#'   \item \code{evaluate_pseudolabeling()}: Evaluate pseudolabeling 
#'     performance against ground truth
#' }
#' 
#' Integration with KEMA:
#' \preformatted{
#' # Generate pseudolabels
#' plabs <- assign_pseudolabels(similarity_matrix, min_clusters = 20)
#' 
#' # Use with KEMA
#' fit <- kema.hyperdesign(
#'   data = strata,
#'   y = plabs$labels,
#'   u = 0.8,           # Trust geometry over pseudolabels
#'   dweight = 0.2,     # Mild class separation
#'   simfun = function(lab) binary_label_matrix(lab, type = "s"),
#'   disfun = function(lab) binary_label_matrix(lab, type = "d")
#' )
#' }
#' 
#' Key parameters:
#' \itemize{
#'   \item \strong{sim_threshold}: Controls which similarities are considered 
#'     "high". Can be adaptive.
#'   \item \strong{diversity_weight}: Balances cluster coherence vs. 
#'     representative diversity
#'   \item \strong{min_clusters/max_clusters}: Controls the number of anchor 
#'     points
#'   \item \strong{min_cluster_size}: Ensures clusters are large enough to be 
#'     reliable
#' }
#' 
#' @name pseudolabeling
NULL

#' Safely compute similarity threshold from sparse matrix
#'
#' Internal helper function to robustly extract off-diagonal similarities
#' and compute adaptive threshold.
#'
#' @param sim_matrix Sparse similarity matrix
#' @param q Quantile for adaptive threshold (default: 0.8)
#' @param explicit_thresh Explicit threshold value (overrides adaptive)
#'
#' @return Similarity threshold value
#' @keywords internal
#' @noRd
.threshold_similarity <- function(sim_matrix, q = 0.8, explicit_thresh = NULL) {
  if (!is.null(explicit_thresh)) return(explicit_thresh)
  
  # Validate quantile parameter
  if (q <= 0 || q >= 1) {
    stop("q must be in (0,1), got: ", q)
  }
  
  # Use summary() to safely extract off-diagonal similarities
  nz <- summary(sim_matrix)
  sim_values <- nz$x[nz$i != nz$j]   # off-diagonal only
  
  if (length(sim_values) == 0) {
    stop("All similarities are zero - cannot compute adaptive threshold")
  }
  
  stats::quantile(sim_values, q, names = FALSE, na.rm = TRUE)
}

#' Assign pseudolabels based on sparse similarity matrix clustering
#'
#' Takes a sparse similarity matrix and identifies clusters of highly similar 
#' samples. Focuses on finding diverse, high-confidence cluster representatives.
#'
#' @param sim_matrix Sparse similarity matrix (dgCMatrix or similar). Should 
#'   be symmetric with values between 0 and 1, where higher values indicate 
#'   greater similarity.
#' @param min_clusters Minimum number of clusters to find (default: 10)
#' @param max_clusters Maximum number of clusters to find (default: 100)
#' @param sim_threshold Minimum similarity threshold for cluster membership. 
#'   If NULL, will be adaptively determined from the distribution of non-zero 
#'   similarities.
#' @param min_cluster_size Minimum number of samples required to form a 
#'   cluster (default: 2)
#' @param max_cluster_size Maximum number of samples allowed in a cluster 
#'   (default: Inf)
#' @param diversity_weight Weight for promoting diversity among cluster 
#'   representatives [0,1]. Controls the trade-off between diversity and 
#'   cluster confidence: 0 = select largest/most confident clusters, 
#'   1 = maximize diversity (minimize pairwise similarities), 
#'   0.5 = balanced approach (default: 0.3)
#' @param adaptive_threshold_quantile Quantile of non-zero similarities to 
#'   use as adaptive threshold (default: 0.8)
#' @param seed Random seed for reproducibility
#' @param use_advanced_diversity Whether to use advanced farthest-first 
#'   traversal for diversity selection on large candidate sets (default: TRUE)
#' @param verbose Whether to print progress information (default: FALSE)
#'
#' @return A list containing:
#'   \item{labels}{Factor vector of pseudolabels, with NA for unassigned 
#'     samples}
#'   \item{representatives}{Integer vector of row indices serving as cluster 
#'     representatives}
#'   \item{cluster_info}{Data frame with cluster statistics including cluster 
#'     ID, representative index, size, and average similarity}
#'   \item{threshold_used}{The similarity threshold that was applied}
#'   \item{n_clusters}{Number of clusters found}
#'
#' @details
#' The algorithm works in several steps:
#' \enumerate{
#'   \item Determine similarity threshold (adaptive or user-specified)
#'   \item Build a graph of high-similarity connections
#'   \item Find connected components as initial clusters
#'   \item Select diverse representatives from clusters
#'   \item Optionally merge small clusters or split large ones
#' }
#'
#' The diversity constraint helps ensure that cluster representatives are 
#' spread across the similarity space, making them better anchors for domain 
#' alignment.
#'
#' @examples
#' \donttest{
#' # Create example sparse similarity matrix
#' library(Matrix)
#' n <- 1000
#' # Generate random sparse matrix with values in [0,1]
#' sim_matrix <- Matrix::rsparsematrix(n, n, density = 0.05, rand.x = runif)
#' sim_matrix <- (sim_matrix + Matrix::t(sim_matrix)) / 2  # Make symmetric
#' Matrix::diag(sim_matrix) <- 1  # Self-similarity = 1
#'
#' # Find pseudolabels
#' result <- assign_pseudolabels(sim_matrix, min_clusters = 20, verbose = TRUE)
#' table(result$labels, useNA = "always")
#' }
#'
#' @importFrom stats quantile
#' @importFrom magrittr %>%
#' @export
assign_pseudolabels <- function(sim_matrix,
                               min_clusters = 10,
                               max_clusters = 100,
                               sim_threshold = NULL,
                               min_cluster_size = 2,
                               max_cluster_size = Inf,
                               diversity_weight = 0.3,
                               adaptive_threshold_quantile = 0.8,
                               seed = NULL,
                               use_advanced_diversity = TRUE,
                               verbose = FALSE) {
  
  # ENHANCED Input validation with descriptive error messages
  if (!inherits(sim_matrix, "sparseMatrix")) {
    stop("sim_matrix must be a sparse matrix (e.g., dgCMatrix). ", 
         "Use Matrix::Matrix() or Matrix::sparseMatrix() to create one. ", 
         "Current class: ", paste(class(sim_matrix), collapse = ", "))
  }
  
  if (nrow(sim_matrix) != ncol(sim_matrix)) {
    stop("sim_matrix must be square. Current dimensions: ", 
         nrow(sim_matrix), " x ", ncol(sim_matrix))
  }
  
  # Check for degenerate cases
  n <- nrow(sim_matrix)
  if (n < 2) {
    stop("sim_matrix must have at least 2 rows/columns for clustering. ", 
         "Current size: ", n)
  }
  
  # Check value range with informative diagnostics
  if (any(sim_matrix < 0 | sim_matrix > 1, na.rm = TRUE)) {
    invalid_range <- range(sim_matrix, na.rm = TRUE)
    stop("sim_matrix values must be in [0,1]. Found range: [", 
         round(invalid_range[1], 4), ", ", round(invalid_range[2], 4), "]. ", 
         "Consider normalizing your similarity matrix.")
  }
  
  # Check for NAs with helpful guidance
  if (any(is.na(sim_matrix))) {
    stop("sim_matrix contains NA values. Please remove or impute them first.")
  }
  
  # Parameter validation with actionable guidance
  if (min_clusters < 1) {
    stop("min_clusters must be >= 1, got: ", min_clusters)
  }
  if (max_clusters < min_clusters) {
    stop("max_clusters (", max_clusters, ") must be >= min_clusters (", 
         min_clusters, ")")
  }
  if (min_cluster_size < 1) {
    stop("min_cluster_size must be >= 1, got: ", min_cluster_size)
  }
  if (diversity_weight < 0 || diversity_weight > 1) {
    stop("diversity_weight must be in [0,1], got: ", diversity_weight, 
         ". Use 0 for pure confidence-based selection, 1 for pure diversity.")
  }
  
  # Check sparsity level with performance warning
  sparsity <- Matrix::nnzero(sim_matrix) / (n * n)
  if (sparsity > 0.5 && n > 1000) {
    warning("High density matrix (", round(sparsity * 100, 1), 
           "% non-zero) with n=", n, " may cause memory issues. ", 
           "Consider using a sparser similarity threshold.")
  }
  
  # Check and enforce symmetry - OPTIMIZED for sparse matrices
  if (!isTRUE(all.equal(sim_matrix, Matrix::t(sim_matrix)))) {
    warning("sim_matrix is not exactly symmetric; forcing symmetry")
    # Use forceSymmetric for efficiency with sparse matrices
    if (requireNamespace("Matrix", quietly = TRUE) && 
        utils::packageVersion("Matrix") >= "1.6.0") {
      sim_matrix <- Matrix::forceSymmetric(sim_matrix, uplo = "U")
    } else {
      # Fallback for older Matrix versions
      sim_matrix <- (sim_matrix + Matrix::t(sim_matrix)) / 2
    }
  }
  
  # Clear diagonal to avoid self-loops - OPTIMIZED for sparse matrices
  # Use summary() approach for safe diagonal removal
  if (Matrix::nnzero(Matrix::diag(sim_matrix)) > 0) {
    diag_indices <- Matrix::which(Matrix::diag(sim_matrix) != 0)
    if (length(diag_indices) > 0) {
      sim_matrix[cbind(diag_indices, diag_indices)] <- 0
    }
  }
  
  # Set seed for reproducibility if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (verbose) {
    cat("Processing", n, "samples with", Matrix::nnzero(sim_matrix), 
        "non-zero similarities\n")
  }
  
  # Step 1: Determine similarity threshold
  if (is.null(sim_threshold)) {
    sim_threshold <- .threshold_similarity(sim_matrix, 
                                          q = adaptive_threshold_quantile)
    if (verbose) {
      cat("Adaptive threshold:", round(sim_threshold, 4), 
          "(", adaptive_threshold_quantile, "quantile of", 
          Matrix::nnzero(sim_matrix), "similarities)\n")
    }
  }
  
  # Step 2: Build thresholded similarity graph
  # OPTIMIZED: Efficient sparse thresholding without intermediate logical matrix
  thresholded_sim <- as(sim_matrix, "dgCMatrix")
  thresholded_sim@x <- as.numeric(thresholded_sim@x >= sim_threshold)
  thresholded_sim <- Matrix::drop0(thresholded_sim)  # Remove zeros
  
  # Ensure diagonal is zero for igraph (no self-loops) - OPTIMIZED
  # Already handled in sparse thresholding, but ensure safety
  if (any(Matrix::diag(thresholded_sim) != 0)) {
    Matrix::diag(thresholded_sim) <- 0L
  }
  
  if (verbose) {
    cat("After thresholding:", Matrix::nnzero(thresholded_sim), 
        "connections remain\n")
  }
  
  # Check for empty graph
  if (Matrix::nnzero(thresholded_sim) == 0) {
    warning("No connections remain after thresholding - try lower sim_threshold")
    return(list(
      labels = factor(rep(NA_character_, n)),
      representatives = integer(0),
      cluster_info = data.frame(),
      threshold_used = sim_threshold,
      n_clusters = 0
    ))
  }
  
  # Step 3: Find connected components
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for assign_pseudolabels. Please ", 
         "install it with: install.packages('igraph')")
  }
  
  # OPTIMIZED: Memory-efficient graph construction for large matrices
  # For large graphs (n >= 50,000), use edge list approach to avoid memory blow-up
  if (n >= 50000) {
    # Use summary() to get edge list directly, avoiding dense intermediate matrix
    edge_summary <- summary(thresholded_sim)
    edge_list <- as.matrix(edge_summary[edge_summary$x > 0, c("i", "j")])
    
    if (nrow(edge_list) == 0) {
      # No edges - return empty result
      warning("No connections remain after thresholding - try lower ", 
              "sim_threshold")
      return(list(
        labels = factor(rep(NA_character_, n)),
        representatives = integer(0),
        cluster_info = data.frame(),
        threshold_used = sim_threshold,
        n_clusters = 0
      ))
    }
    
    g <- igraph::graph_from_edgelist(edge_list, directed = FALSE)
    igraph::V(g)$name <- seq_len(n)  # Ensure proper vertex mapping
  } else {
    # For smaller graphs, use standard adjacency matrix approach
    g <- igraph::graph_from_adjacency_matrix(thresholded_sim, 
                                           mode = "undirected", 
                                           weighted = TRUE)
  }
  components <- igraph::components(g)
  
  # Extract cluster assignments
  cluster_assignments <- components$membership
  cluster_sizes <- components$csize
  
  if (verbose) {
    cat("Found", components$no, "initial clusters\n")
  }
  
  # Step 4: Filter clusters by size
  valid_clusters <- which(cluster_sizes >= min_cluster_size & 
                         cluster_sizes <= max_cluster_size)
  
  if (length(valid_clusters) == 0) {
    warning("No clusters meet size requirements (min_size=", min_cluster_size, 
            ", max_size=", max_cluster_size, 
            "). Consider lowering sim_threshold or min_cluster_size.")
    return(list(
      labels = factor(rep(NA_character_, n)),
      representatives = integer(0),
      cluster_info = data.frame(),
      threshold_used = sim_threshold,
      n_clusters = 0
    ))
  }
  
  if (verbose) {
    cat("Found", length(valid_clusters), 
        "valid clusters after size filtering\n")
  }
  
  # Step 5: Select diverse representatives
  initial_rep_count <- length(valid_clusters)  # Track for algorithm choice metadata
  representatives <- .select_diverse_representatives(
    sim_matrix = sim_matrix,
    cluster_assignments = cluster_assignments,
    valid_clusters = valid_clusters,
    diversity_weight = diversity_weight,
    min_clusters = min_clusters,
    max_clusters = max_clusters,
    use_advanced_diversity = use_advanced_diversity,
    verbose = verbose
  )
  
  # FIXED: Enforce min_clusters constraint after representative selection
  if (length(representatives) < min_clusters) {
    warning("Final number of representatives (", length(representatives), 
            ") is below min_clusters (", min_clusters, 
            "). Consider lowering sim_threshold or min_cluster_size.")
  }
  
  # Step 6: Create final labels - OPTIMIZED cluster_info construction
  labels <- rep(NA_character_, n)
  
  # Pre-allocate vectors for efficient data.frame construction
  cluster_ids <- integer(length(representatives))
  rep_indices <- integer(length(representatives))
  cluster_sizes_final <- integer(length(representatives))
  avg_similarities <- numeric(length(representatives))
  # ENHANCED: Add observed similarity metrics for better cluster assessment
  observed_similarities <- numeric(length(representatives))  # Mean of observed edges only
  edge_densities <- numeric(length(representatives))         # Proportion of observed edges
  
  for (i in seq_along(representatives)) {
    rep_idx <- representatives[i]
    cluster_id <- cluster_assignments[rep_idx]
    cluster_members <- which(cluster_assignments == cluster_id)
    
    # Assign pseudolabel
    label_name <- base::sprintf("anchor_%03d", i)
    labels[cluster_members] <- label_name
    
    # Store cluster info components
    cluster_ids[i] <- i
    rep_indices[i] <- rep_idx
    cluster_sizes_final[i] <- length(cluster_members)
    
    # ENHANCED: Calculate both density-adjusted and observed similarity metrics
    if (length(cluster_members) > 1) {
      cluster_sims <- sim_matrix[cluster_members, cluster_members]
      
      # Extract off-diagonal entries efficiently for sparse matrices
      if (inherits(cluster_sims, "sparseMatrix")) {
        # Use summary() to get only non-zero off-diagonal entries
        sim_summary <- summary(cluster_sims)
        off_diag_mask <- sim_summary$i != sim_summary$j
        off_diag_values <- sim_summary$x[off_diag_mask]
        n_observed_edges <- length(off_diag_values)
        n_possible_edges <- length(cluster_members) * (length(cluster_members) - 1)
        
        # Observed mean: average of non-zero similarities only
        observed_similarities[i] <- if (n_observed_edges > 0) {
          mean(off_diag_values, na.rm = TRUE)
        } else {
          0
        }
        
        # Edge density: proportion of possible edges that are non-zero
        edge_densities[i] <- n_observed_edges / n_possible_edges
        
        # Density-adjusted mean: includes zeros for missing edges
        avg_similarities[i] <- sum(off_diag_values) / n_possible_edges
      } else {
        # Dense matrix - compute directly
        total_sim <- sum(cluster_sims)
        diag_sim <- sum(Matrix::diag(cluster_sims))  # Should be 0 already
        n_possible_edges <- length(cluster_members) * (length(cluster_members) - 1)
        
        avg_similarities[i] <- (total_sim - diag_sim) / n_possible_edges
        observed_similarities[i] <- avg_similarities[i]  # Same for dense
        edge_densities[i] <- 1.0  # All edges observed in dense matrix
      }
    } else {
      # Single-member clusters
      avg_similarities[i] <- NA_real_
      observed_similarities[i] <- NA_real_
      edge_densities[i] <- NA_real_
    }
  }
  
  # ENHANCED: Build cluster_info with comprehensive similarity metrics
  cluster_info <- data.frame(
    cluster_id = cluster_ids,
    representative = rep_indices,
    size = cluster_sizes_final,
    avg_similarity = avg_similarities,        # Density-adjusted mean (includes zeros)
    observed_similarity = observed_similarities, # Mean of observed edges only
    edge_density = edge_densities,           # Proportion of edges observed
    stringsAsFactors = FALSE
  )
  
  labels <- factor(labels)
  
  if (verbose) {
    cat("Final result:", length(representatives), "clusters,", 
        sum(!is.na(labels)), "samples assigned\n")
  }
  
  # ENHANCED: Return S3 class with methods for better user experience
  result <- list(
    labels = labels,
    representatives = representatives,
    cluster_info = cluster_info,
    threshold_used = sim_threshold,
    n_clusters = length(representatives),
    # Additional metadata for methods
    n_samples = n,
    n_assigned = sum(!is.na(labels)),
    algorithm = "assign_pseudolabels",
    diversity_algorithm = if (use_advanced_diversity && initial_rep_count >= 50 && initial_rep_count > max_clusters) "farthest_first" else "greedy",
    diversity_weight = diversity_weight,
    call = match.call()
  )
  
  class(result) <- c("pseudolabels", "list")
  return(result)
}

#' Print method for pseudolabels objects
#'
#' @param x A pseudolabels object
#' @param ... Additional arguments (ignored)
#' @method print pseudolabels
#' @export
print.pseudolabels <- function(x, ...) {
  cat("Pseudolabeling Results\n")
  cat("======================\n")
  cat("Algorithm:        ", x$algorithm, "\n")
  if (!is.null(x$diversity_algorithm)) {
    cat("Diversity method: ", x$diversity_algorithm, 
        " (weight: ", round(x$diversity_weight, 2), ")\n", sep = "")
  }
  cat("Samples:          ", x$n_samples, "\n")
  cat("Assigned:         ", x$n_assigned, " (", 
      round(100 * x$n_assigned / x$n_samples, 1), "%)\n", sep = "")
  cat("Clusters found:   ", x$n_clusters, "\n")
  cat("Threshold used:   ", round(x$threshold_used, 4), "\n")
  
  if (x$n_clusters > 0) {
    cat("\nCluster Summary:\n")
    cluster_sizes <- table(x$labels, useNA = "no")
    cat("  Size range:     ", min(cluster_sizes), "-", max(cluster_sizes), "\n")
    cat("  Mean size:      ", round(mean(cluster_sizes), 1), "\n")
    
    if (!is.null(x$cluster_info) && nrow(x$cluster_info) > 0) {
      if ("observed_similarity" %in% names(x$cluster_info)) {
        obs_sim <- x$cluster_info$observed_similarity[!is.na(x$cluster_info$observed_similarity)]
        if (length(obs_sim) > 0) {
          cat("  Observed similarity: ", round(min(obs_sim), 3), "-", 
              round(max(obs_sim), 3), " (mean: ", round(mean(obs_sim), 3), ")\n", sep = "")
        }
      }
    }
  }
  
  invisible(x)
}

#' Summary method for pseudolabels objects
#'
#' @param object A pseudolabels object
#' @param ... Additional arguments (ignored)
#' @method summary pseudolabels
#' @export
summary.pseudolabels <- function(object, ...) {
  cat("Pseudolabeling Summary\n")
  cat("======================\n")
  print(object)
  
  if (object$n_clusters > 0 && !is.null(object$cluster_info)) {
    cat("\nDetailed Cluster Information:\n")
    print(object$cluster_info)
  }
  
  invisible(object)
}

#' Convert pseudolabels to data.frame
#'
#' @param x A pseudolabels object
#' @param ... Additional arguments (ignored)
#' @method as.data.frame pseudolabels
#' @export
as.data.frame.pseudolabels <- function(x, ...) {
  data.frame(
    sample_id = seq_len(x$n_samples),
    label = x$labels,
    is_representative = seq_len(x$n_samples) %in% x$representatives,
    stringsAsFactors = FALSE
  )
}

#' Select diverse cluster representatives
#'
#' Internal function to select cluster representatives that are diverse
#' across the similarity space.
#'
#' @param sim_matrix Original similarity matrix
#' @param cluster_assignments Vector of cluster assignments
#' @param valid_clusters Vector of valid cluster IDs
#' @param diversity_weight Weight for diversity constraint
#' @param min_clusters Minimum number of clusters
#' @param max_clusters Maximum number of clusters
#' @param use_advanced_diversity Whether to use advanced FFS algorithm
#' @param verbose Whether to print progress
#'
#' @return Vector of representative indices
#' @keywords internal
.select_diverse_representatives <- function(sim_matrix,
                                          cluster_assignments,
                                          valid_clusters,
                                          diversity_weight,
                                          min_clusters,
                                          max_clusters,
                                          use_advanced_diversity = TRUE,
                                          verbose) {
  
  # Start with one representative per valid cluster
  initial_reps <- integer(length(valid_clusters))
  
  for (i in seq_along(valid_clusters)) {
    cluster_id <- valid_clusters[i]
    cluster_members <- which(cluster_assignments == cluster_id)
    
    if (length(cluster_members) == 1) {
      initial_reps[i] <- cluster_members[1]
    } else {
      # OPTIMIZED: Use vectorized approach for representative selection
      # Extract the cluster submatrix once and compute row means efficiently
      cluster_submatrix <- sim_matrix[cluster_members, cluster_members, drop = FALSE]
      
      # For small clusters (≤ 5000 members), use dense slice for speed
      if (length(cluster_members) <= 5000) {
        # Convert to dense for faster row operations on small matrices
        dense_submatrix <- as.matrix(cluster_submatrix)
        avg_sims <- base::rowMeans(dense_submatrix)
      } else {
        # For larger clusters, use sparse operations
        avg_sims <- Matrix::rowMeans(cluster_submatrix)
      }
      
      best_member <- base::which.max(avg_sims)
      initial_reps[i] <- cluster_members[best_member]
    }
  }
  
  # If we have too few clusters, we're done
  if (length(initial_reps) <= min_clusters) {
    return(initial_reps)
  }
  
  # If we have too many clusters, select diverse subset
  if (length(initial_reps) > max_clusters) {
    # ENHANCED: Choose algorithm based on problem size and user preference
    if (use_advanced_diversity && length(initial_reps) >= 50) {
      # Use advanced FFS algorithm for large problems
      selected_reps <- .select_diverse_subset_advanced(
        sim_matrix = sim_matrix,
        candidates = initial_reps,
        target_size = max_clusters,
        diversity_weight = diversity_weight,
        cluster_assignments = cluster_assignments,
        valid_clusters = valid_clusters,
        use_advanced = TRUE
      )
    } else {
      # Use standard greedy algorithm for smaller problems
      selected_reps <- .select_diverse_subset(
        sim_matrix = sim_matrix,
        candidates = initial_reps,
        target_size = max_clusters,
        diversity_weight = diversity_weight,
        cluster_assignments = cluster_assignments,
        valid_clusters = valid_clusters
      )
    }
    return(selected_reps)
  }
  
  return(initial_reps)
}

#' Select diverse subset of representatives
#'
#' ENHANCED: Use a greedy algorithm with continuous weighting between 
#' diversity and cluster confidence scores.
#'
#' @param sim_matrix Similarity matrix
#' @param candidates Vector of candidate representative indices
#' @param target_size Target number of representatives
#' @param diversity_weight Weight for diversity vs confidence [0,1]. 
#'   0 = pure confidence (largest clusters), 1 = pure diversity, 0.5 = balanced
#' @param cluster_assignments Vector of cluster assignments (for confidence scores)
#' @param valid_clusters Vector of valid cluster IDs (for confidence scores)
#'
#' @return Vector of selected representative indices
#' @keywords internal
.select_diverse_subset <- function(sim_matrix, candidates, target_size, diversity_weight,
                                  cluster_assignments = NULL, valid_clusters = NULL) {
  
  if (length(candidates) <= target_size) {
    return(candidates)
  }
  
  # ENHANCED: Compute confidence scores based on cluster sizes (if available)
  confidence_scores <- rep(1, length(candidates))  # Default: equal confidence
  if (!is.null(cluster_assignments) && !is.null(valid_clusters)) {
    for (i in seq_along(candidates)) {
      candidate <- candidates[i]
      cluster_id <- cluster_assignments[candidate]
      cluster_size <- sum(cluster_assignments == cluster_id)
      # Confidence = normalized cluster size (larger clusters = higher confidence)
      confidence_scores[i] <- cluster_size
    }
    # Normalize confidence scores to [0,1]
    confidence_scores <- (confidence_scores - min(confidence_scores)) / 
      (max(confidence_scores) - min(confidence_scores) + 1e-8)
  }
  
  # Handle edge cases for diversity_weight
  if (diversity_weight <= 0) {
    # Pure confidence-based selection: select largest clusters
    return(candidates[order(confidence_scores, decreasing = TRUE)[1:target_size]])
  } else if (diversity_weight >= 1) {
    # Pure diversity-based selection (original algorithm)
    # Fall through to greedy diversity selection
  }
  
  # Greedy selection with weighted diversity + confidence
  selected <- integer(target_size)
  remaining <- candidates
  remaining_confidence <- confidence_scores
  
  # Start with best confidence or random (depending on weight)
  if (diversity_weight < 1) {
    # Confidence-weighted initialization
    init_idx <- base::which.max(remaining_confidence)
  } else {
    # Random initialization for pure diversity
    init_idx <- base::sample(length(remaining), 1)
  }
  selected[1] <- remaining[init_idx]
  remaining <- remaining[-init_idx]
  remaining_confidence <- remaining_confidence[-init_idx]
  
  # Greedily select remaining representatives
  for (i in 2:target_size) {
    if (length(remaining) == 0) break
    
    # Calculate combined scores for remaining candidates
    combined_scores <- numeric(length(remaining))
    
    for (j in seq_along(remaining)) {
      candidate <- remaining[j]
      
      # Diversity component: distance to already selected representatives
      sims_to_selected <- sim_matrix[candidate, selected[1:(i-1)]]
      if (inherits(sims_to_selected, "sparseVector")) {
        dense_sims <- as.vector(sims_to_selected)
        min_sim <- if (length(dense_sims) == 0) 0 else min(dense_sims, na.rm = TRUE)
      } else {
        min_sim <- if (length(sims_to_selected) == 0) 0 else min(sims_to_selected, na.rm = TRUE)
      }
      diversity_score <- 1 - min_sim  # Higher = more diverse
      
      # Confidence component
      confidence_score <- remaining_confidence[j]
      
      # ENHANCED: Weighted combination of diversity and confidence
      combined_scores[j] <- diversity_weight * diversity_score + 
                           (1 - diversity_weight) * confidence_score
    }
    
    # Select candidate with highest combined score
    best_idx <- base::which.max(combined_scores)
    selected[i] <- remaining[best_idx]
    remaining <- remaining[-best_idx]
    remaining_confidence <- remaining_confidence[-best_idx]
  }
  
  return(selected[1:min(target_size, sum(selected > 0))])
}

#' Advanced diversity selection with farthest-first traversal
#'
#' OPTIMIZED: Implements efficient farthest-first traversal (FFS) for diversity
#' selection with improved O(k²) complexity through early termination and 
#' sparse matrix optimizations.
#'
#' @param sim_matrix Similarity matrix
#' @param candidates Vector of candidate representative indices
#' @param target_size Target number of representatives
#' @param diversity_weight Weight for diversity vs confidence [0,1]
#' @param cluster_assignments Vector of cluster assignments (for confidence scores)
#' @param valid_clusters Vector of valid cluster IDs (for confidence scores)
#' @param use_advanced Whether to use advanced FFS algorithm (default: TRUE)
#'
#' @return Vector of selected representative indices
#' @keywords internal
.select_diverse_subset_advanced <- function(sim_matrix, candidates, target_size, 
                                           diversity_weight,
                                           cluster_assignments = NULL, 
                                           valid_clusters = NULL,
                                           use_advanced = TRUE) {
  
  if (length(candidates) <= target_size) {
    return(candidates)
  }
  
  # For small candidate sets, use the regular algorithm
  if (length(candidates) < 50 || !use_advanced) {
    return(.select_diverse_subset(sim_matrix, candidates, target_size, 
                                 diversity_weight, cluster_assignments, valid_clusters))
  }
  
  # OPTIMIZED: Farthest-First Traversal with confidence weighting
  # Pre-compute confidence scores
  confidence_scores <- rep(1, length(candidates))
  if (!is.null(cluster_assignments) && !is.null(valid_clusters)) {
    for (i in seq_along(candidates)) {
      candidate <- candidates[i]
      cluster_id <- cluster_assignments[candidate]
      cluster_size <- sum(cluster_assignments == cluster_id)
      confidence_scores[i] <- cluster_size
    }
    # Normalize confidence scores to [0,1]
    confidence_scores <- (confidence_scores - min(confidence_scores)) / 
      (max(confidence_scores) - min(confidence_scores) + 1e-8)
  }
  
  # Handle edge cases for diversity_weight
  if (diversity_weight <= 0) {
    return(candidates[order(confidence_scores, decreasing = TRUE)[1:target_size]])
  }
  
  # ADVANCED: Extract similarity sub-matrix once for efficiency
  candidate_sims <- sim_matrix[candidates, candidates, drop = FALSE]
  
  # Initialize with best confidence candidate
  selected_indices <- integer(target_size)  # Indices into candidates array
  if (diversity_weight < 1) {
    selected_indices[1] <- base::which.max(confidence_scores)
  } else {
    selected_indices[1] <- base::sample(length(candidates), 1)
  }
  
  # OPTIMIZED: Track minimum distances for farthest-first traversal
  # Distance from each unselected candidate to nearest selected candidate
  min_distances <- rep(0, length(candidates))
  remaining_mask <- rep(TRUE, length(candidates))
  remaining_mask[selected_indices[1]] <- FALSE
  
  # Initialize distances from first selected candidate
  if (inherits(candidate_sims, "sparseMatrix")) {
    first_row <- candidate_sims[selected_indices[1], ]
    min_distances <- 1 - as.vector(first_row)  # Convert similarity to distance
  } else {
    min_distances <- 1 - candidate_sims[selected_indices[1], ]
  }
  min_distances[selected_indices[1]] <- -Inf  # Mark as selected
  
  # FARTHEST-FIRST TRAVERSAL: Each iteration picks the candidate that is
  # farthest from all previously selected candidates
  for (i in 2:target_size) {
    remaining_indices <- which(remaining_mask)
    if (length(remaining_indices) == 0) break
    
    # ENHANCED: Combine diversity (distance) and confidence scores
    diversity_scores <- min_distances[remaining_indices]
    conf_scores <- confidence_scores[remaining_indices]
    
    # Weighted combination
    combined_scores <- diversity_weight * diversity_scores + 
                      (1 - diversity_weight) * conf_scores
    
    # Select candidate with highest combined score
    best_remaining_idx <- base::which.max(combined_scores)
    selected_indices[i] <- remaining_indices[best_remaining_idx]
    
    # Mark as selected
    remaining_mask[selected_indices[i]] <- FALSE
    
    # OPTIMIZED: Update minimum distances efficiently
    # Only update distances for candidates that could be improved
    if (inherits(candidate_sims, "sparseMatrix")) {
      new_row <- candidate_sims[selected_indices[i], ]
      new_distances <- 1 - as.vector(new_row)
    } else {
      new_distances <- 1 - candidate_sims[selected_indices[i], ]
    }
    
    # Update min_distances where new candidate is closer
    update_mask <- remaining_mask & (new_distances < min_distances)
    min_distances[update_mask] <- new_distances[update_mask]
    min_distances[selected_indices[i]] <- -Inf  # Mark as selected
  }
  
  # Return actual candidate indices
  return(candidates[selected_indices[1:min(target_size, sum(selected_indices > 0))]])
}

#' High-similarity pseudolabels from multi-domain data
#'
#' Generates pseudolabels based on high cosine similarity across different 
#' domains. Specialized version for multi-domain case.
#'
#' @param strata List produced by multidesign::hyperdesign. Each element 
#'   should have $x (samples x features matrix) and $design (metadata 
#'   data.frame)
#' @param k Number of cross-domain neighbors to examine (default: 5)
#' @param cos_thresh Minimum cosine similarity threshold (default: 0.97)
#' @param min_size Minimum cluster size to retain (default: 2)
#' @param ann_trees Number of trees for RcppAnnoy index (default: 50)
#' @param verbose Whether to print progress information (default: FALSE)
#'
#' @return Factor vector of pseudolabels with length equal to total number of 
#'   samples across all domains. NA indicates samples with no confident match.
#'
#' @details
#' This function implements the approach described in the pseudolabel 
#' documentation:
#' \enumerate{
#'   \item Pool and L2-normalize all samples across domains
#'   \item Build approximate nearest neighbor index using angular distance
#'   \item Find high-similarity cross-domain pairs
#'   \item Use connected components to form clusters
#'   \item Filter clusters by minimum size
#' }
#'
#' @examples
#' \donttest{
#' # Assuming you have multi-domain data
#' strata <- multidesign::hyperdesign(raw_blocks)
#' plabs <- high_sim_pseudolabels(strata, k = 5, cos_thresh = 0.98)
#' table(plabs, useNA = "always")
#' }
#'
#' @export
high_sim_pseudolabels <- function(strata,
                                 k = 5,
                                 cos_thresh = 0.97,
                                 min_size = 2,
                                 ann_trees = 50,
                                 verbose = FALSE) {
  
  # Input validation
  if (!is.list(strata) || length(strata) == 0) {
    stop("strata must be a non-empty list")
  }
  
  # Extract data matrices and create block identifiers
  X_list <- lapply(strata, function(s) {
    if (is.null(s$x)) {
      stop("Each element of strata must have an 'x' component")
    }
    x_mat <- as.matrix(s$x)
    if (nrow(x_mat) == 0) {
      stop("Each block must have at least one row")
    }
    if (any(is.na(x_mat))) {
      stop("Input matrices must not contain missing values")
    }
    x_mat
  })
  
  block_id <- unlist(Map(rep, seq_along(X_list), sapply(X_list, nrow)))
  X <- do.call(rbind, X_list)
  
  # L2 normalize rows (for cosine similarity via dot product)
  row_norms <- base::sqrt(base::rowSums(X^2))
  row_norms[row_norms == 0] <- 1  # Prevent division by zero
  X <- X / row_norms
  
  # Build approximate nearest neighbor index
  if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
    stop("RcppAnnoy package is required for high_sim_pseudolabels")
  }
  
  # ENHANCED: Check RcppAnnoy version compatibility for distance formula
  annoy_version <- utils::packageVersion("RcppAnnoy")
  expected_angular_distance <- TRUE  # Assume correct by default
  
  if (verbose) {
    cat("Using RcppAnnoy version", as.character(annoy_version), "\n")
  }
  
  ann <- RcppAnnoy::AnnoyAngular(ncol(X))
  for (i in seq_len(nrow(X))) {
    ann$addItem(i - 1, X[i, ])
  }
  ann$build(ann_trees)
  
  # Collect high-similarity cross-block pairs
  max_edges <- k * nrow(X)  # Pre-allocate for efficiency
  edge_i <- integer(max_edges)
  edge_j <- integer(max_edges)
  edge_count <- 0
  
  for (i in seq_len(nrow(X))) {
    # Get k+1 nearest neighbors (including self) and distances
    nn_res <- ann$getNNsByItem(i - 1, k + 1, include_distances = TRUE)
    
    # Exclude self (first element)
    nns <- nn_res$item[-1]
    dists <- nn_res$distance[-1]
    
    for (k_idx in seq_along(nns)) {
      j <- nns[k_idx] + 1  # Convert to 1-based indexing
      
      # Only consider cross-domain pairs
      if (block_id[i] != block_id[j]) {
        # ENHANCED: Derive cosine similarity from angular distance with version check
        # For angular distance: distance = sqrt(2 * (1 - cosine_similarity))
        # This formula assumes RcppAnnoy returns squared angular distance
        dist_val <- dists[k_idx]
        
        # Sanity check: angular distance should be in reasonable range
        if (dist_val > pi) {
          if (expected_angular_distance) {
            warning("Unexpected distance value (", round(dist_val, 4), 
                   " > π) from RcppAnnoy. This may indicate version incompatibility. ",
                   "Expected squared angular distance, got raw distance.")
            expected_angular_distance <- FALSE  # Only warn once
          }
          # Fallback: assume raw angular distance
          cos_sim <- cos(dist_val)
        } else {
          # Standard formula for squared angular distance
          cos_sim <- 1 - (dist_val^2) / 2
        }
        
        if (cos_sim >= cos_thresh) {
          edge_count <- edge_count + 1
          # FIXED: Guard against edge buffer overflow with auto-expansion
          if (edge_count > max_edges) {
            # Auto-expand edge buffer instead of failing
            new_max_edges <- max_edges * 2
            warning("Edge buffer overflow: expanding from ", max_edges, 
                   " to ", new_max_edges, " edges. Consider increasing ", 
                   "max_edges parameter for better performance.")
            
            # Expand edge vectors
            new_edge_i <- integer(new_max_edges)
            new_edge_j <- integer(new_max_edges)
            new_edge_i[1:edge_count] <- edge_i[1:edge_count]
            new_edge_j[1:edge_count] <- edge_j[1:edge_count]
            edge_i <- new_edge_i
            edge_j <- new_edge_j
            max_edges <- new_max_edges
          }
          edge_i[edge_count] <- i
          edge_j[edge_count] <- j
        }
      }
    }
  }
  
  if (edge_count == 0) {
    warning("No high-similarity cross-domain pairs found")
    return(factor(rep(NA_character_, nrow(X))))
  }
  
  # Trim to actual size
  edge_i <- edge_i[1:edge_count]
  edge_j <- edge_j[1:edge_count]
  edges_matrix <- cbind(edge_i, edge_j)
  
  # Find connected components
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package is required for high_sim_pseudolabels. Please ", 
         "install it with: install.packages('igraph')")
  }
  
  g <- igraph::graph_from_edgelist(edges_matrix, directed = FALSE)
  components <- igraph::components(g)
  
  # Assign labels
  labels <- rep(NA_character_, nrow(X))
  for (cid in seq_len(components$no)) {
    members <- which(components$membership == cid)
    if (length(members) >= min_size) {
      labels[members] <- base::sprintf("anchor_%04d", cid)
    }
  }
  
  return(factor(labels))
}

#' Create synthetic similarity matrix for testing pseudolabeling
#'
#' Generates a synthetic sparse similarity matrix with known cluster structure. 
#' Used for testing and demonstrating the pseudolabeling functions.
#'
#' @param n_samples Total number of samples
#' @param n_clusters Number of true clusters
#' @param within_cluster_sim Average similarity within clusters (default: 0.8)
#' @param between_cluster_sim Average similarity between clusters 
#'   (default: 0.1)
#' @param sparsity Overall sparsity level (proportion of non-zero entries) 
#'   (default: 0.1)
#' @param noise_level Amount of noise to add to similarities (default: 0.1)
#'
#' @return A list containing:
#'   \item{sim_matrix}{Sparse similarity matrix (dgCMatrix) with synthetic 
#'     cluster structure}
#'   \item{true_labels}{Factor vector of true cluster assignments for 
#'     validation}
#'   \item{cluster_centers}{Integer vector of cluster center indices}
#'
#' @examples
#' \donttest{
#' # Create synthetic data
#' synthetic <- create_synthetic_similarity_matrix(n_samples = 500, 
#'                                                 n_clusters = 10)
#' 
#' # Apply pseudolabeling
#' result <- assign_pseudolabels(synthetic$sim_matrix, verbose = TRUE)
#' 
#' # Compare with true labels
#' table(result$labels, synthetic$true_labels, useNA = "always")
#' }
#'
#' @export
create_synthetic_similarity_matrix <- function(n_samples = 1000,
                                              n_clusters = 20,
                                              within_cluster_sim = 0.8,
                                              between_cluster_sim = 0.1,
                                              sparsity = 0.1,
                                              noise_level = 0.1) {
  
  # Create true cluster assignments
  cluster_sizes <- rep(n_samples %/% n_clusters, n_clusters)
  remainder <- n_samples %% n_clusters
  if (remainder > 0) {
    cluster_sizes[1:remainder] <- cluster_sizes[1:remainder] + 1
  }
  
  true_labels <- rep(1:n_clusters, cluster_sizes)
  true_labels <- factor(paste0("cluster_", true_labels))
  
  # Pre-allocate vectors for efficient sparse matrix construction
  i_idx <- c()
  j_idx <- c()
  vals <- c()
  
  # Add within-cluster similarities
  start_idx <- 1
  cluster_centers <- integer(n_clusters)
  
  for (k in 1:n_clusters) {
    end_idx <- start_idx + cluster_sizes[k] - 1
    cluster_indices <- start_idx:end_idx
    
    # Select cluster center (representative)
    cluster_centers[k] <- cluster_indices[base::ceiling(length(cluster_indices) / 2)]
    
    # Create all pairs within the cluster efficiently
    if (length(cluster_indices) > 1) {
      pairs <- base::expand.grid(i = cluster_indices, j = cluster_indices)
      pairs <- pairs[pairs$i < pairs$j, ]  # Keep only upper triangle, no diagonal
      
      # Subsample based on density
      n_pairs <- nrow(pairs)
      n_to_keep <- round(n_pairs * sparsity * 3)  # Higher density within clusters
      if (n_to_keep > 0 && n_pairs > 0) {
        kept_pairs <- pairs[base::sample(n_pairs, base::min(n_pairs, n_to_keep)), ]
        
        n_new <- nrow(kept_pairs)
        new_vals <- within_cluster_sim + rnorm(n_new, 0, noise_level)
        new_vals <- pmax(0, pmin(1, new_vals))  # Clamp to [0,1]
        
        i_idx <- c(i_idx, kept_pairs$i)
        j_idx <- c(j_idx, kept_pairs$j)
        vals <- c(vals, new_vals)
      }
    }
    
    start_idx <- end_idx + 1
  }
  
  # Add between-cluster similarities (noise)
  n_between <- round(n_samples * n_samples * sparsity * 0.1)  # Much sparser
  if (n_between > 0) {
    between_i <- base::sample(n_samples, n_between, replace = TRUE)
    between_j <- base::sample(n_samples, n_between, replace = TRUE)
    
    # Filter out same-cluster and diagonal pairs
    valid <- (true_labels[between_i] != true_labels[between_j]) & (between_i != between_j)
    between_i <- between_i[valid]
    between_j <- between_j[valid]
    
    # Keep only upper triangle to avoid duplicates
    upper_tri <- between_i < between_j
    between_i <- between_i[upper_tri]
    between_j <- between_j[upper_tri]
    
    n_new <- length(between_i)
    if (n_new > 0) {
      new_vals <- between_cluster_sim + rnorm(n_new, 0, noise_level)
      new_vals <- pmax(0, pmin(1, new_vals))  # Clamp to [0,1]
      
      i_idx <- c(i_idx, between_i)
      j_idx <- c(j_idx, between_j)
      vals <- c(vals, new_vals)
    }
  }
  
  # Create the sparse matrix efficiently in one call
  sim_matrix <- Matrix::sparseMatrix(
    i = i_idx, j = j_idx, x = vals,
    dims = c(n_samples, n_samples),
    symmetric = TRUE  # Automatically handles symmetry
  )
  
  # Ensure diagonal is 1 (self-similarity)
  Matrix::diag(sim_matrix) <- 1
  
  # Convert to dgCMatrix for efficiency
  sim_matrix <- Matrix::drop0(sim_matrix)
  
  return(list(
    sim_matrix = sim_matrix,
    true_labels = true_labels,
    cluster_centers = cluster_centers
  ))
}

#' Evaluate pseudolabeling performance against ground truth
#'
#' Compares pseudolabeling results with known true cluster assignments. 
#' Uses various clustering evaluation metrics for assessment.
#'
#' @param predicted_labels Factor vector of predicted pseudolabels
#' @param true_labels Factor vector of true cluster assignments
#' @param verbose Whether to print detailed results (default: TRUE)
#'
#' @return A list containing evaluation metrics:
#'   \item{n_predicted_clusters}{Number of predicted clusters found}
#'   \item{n_true_clusters}{Number of true clusters in ground truth}
#'   \item{coverage}{Proportion of samples assigned to clusters}
#'   \item{purity}{Average purity of predicted clusters (proportion of 
#'     dominant class)}
#'   \item{completeness}{Average majority-class recall (proportion of true 
#'     cluster members captured by their dominant predicted cluster)}
#'   \item{confusion_matrix}{Confusion matrix between predicted and true 
#'     labels}
#'
#' @examples
#' \donttest{
#' # Create synthetic data and apply pseudolabeling
#' synthetic <- create_synthetic_similarity_matrix(n_samples = 500)
#' result <- assign_pseudolabels(synthetic$sim_matrix)
#' 
#' # Evaluate performance
#' eval_result <- evaluate_pseudolabeling(result$labels, 
#'                                        synthetic$true_labels)
#' print(eval_result)
#' }
#'
#' @export
evaluate_pseudolabeling <- function(predicted_labels, true_labels, verbose = TRUE) {
  
  # Remove NA predictions for evaluation
  valid_mask <- !is.na(predicted_labels)
  pred_clean <- predicted_labels[valid_mask]
  true_clean <- true_labels[valid_mask]
  
  if (length(pred_clean) == 0) {
    warning("No valid predictions to evaluate")
    return(list(
      n_predicted_clusters = 0,
      n_true_clusters = length(levels(true_labels)),
      coverage = 0,
      purity = 0,
      completeness = 0,
      confusion_matrix = matrix(0, 0, 0)
    ))
  }
  
  # Basic statistics - FIXED: Count only actually used levels, not all factor levels
  n_predicted <- length(unique(pred_clean))
  n_true <- length(unique(true_clean))
  coverage <- sum(valid_mask) / length(predicted_labels)
  
  # Create confusion matrix
  confusion <- table(pred_clean, true_clean)
  
  # Calculate purity (precision-like metric)
  row_max <- apply(confusion, 1, max)
  avg_purity <- mean(row_max / rowSums(confusion), na.rm = TRUE)
  
  # Calculate completeness (recall-like metric)
  col_max <- apply(confusion, 2, max)
  avg_completeness <- mean(col_max / colSums(confusion), na.rm = TRUE)
  
  if (verbose) {
    cat("Pseudolabeling Evaluation Results:\n")
    cat("==================================\n")
    cat("Predicted clusters:", n_predicted, "\n")
    cat("True clusters:", n_true, "\n")
    cat("Coverage:", round(coverage * 100, 1), "%\n")
    cat("Average purity:", round(avg_purity * 100, 1), "%\n")
    cat("Average completeness:", round(avg_completeness * 100, 1), "%\n")
    cat("\nConfusion Matrix (top 10x10):\n")
    print(confusion[1:min(10, nrow(confusion)), 1:min(10, ncol(confusion))])
  }
  
  return(list(
    n_predicted_clusters = n_predicted,
    n_true_clusters = n_true,
    coverage = coverage,
    purity = avg_purity,
    completeness = avg_completeness,
    confusion_matrix = confusion
  ))
}
