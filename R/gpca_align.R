#' Generalized PCA Alignment for Hyperdesign Data
#'
#' Performs Generalized Principal Component Analysis (GPCA) alignment on hyperdesign data structures.
#' This method aligns data from multiple domains using a similarity-based approach that balances
#' within-domain and between-domain relationships.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the label variable to use for alignment
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of components to extract (default: 2)
#' @param simfun Function to compute similarity between labels (default: neighborweights::binary_label_matrix)
#' @param u Trade-off parameter between within-domain and between-domain similarity (0-1, default: 0.5).
#'   Both similarity matrices are normalized to unit Frobenius norm before blending, so u=0.5 
#'   gives truly equal weight to within-domain and between-domain similarities, making this 
#'   parameter interpretable regardless of the scale of the similarity function.
#' @param lambda Regularization parameter for matrix conditioning (default: 0.1)
#' @param row_metric_scale Scaling factor for the row metric matrix M (default: 1). 
#'   Since GPCA is scale-invariant to M up to a constant, this parameter typically doesn't 
#'   affect results but can be adjusted for numerical reasons if needed.
#'
#' @return A multiblock_biprojector object containing the GPCA alignment
#' @importFrom genpca genpca
#' @export
gpca_align.hyperdesign <- function(data, y, 
                    preproc=center(), 
                    ncomp=2,
                    simfun=neighborweights::binary_label_matrix,
                    u=.5,
                    lambda=.1,
                    row_metric_scale=1) {
  
  y <- rlang::enquo(y)
  
  label_list <- purrr::map(data, function(x) x$design %>% select(!!y) %>% pull(!!y))
  labels <- factor(unlist(label_list))
  label_set <- levels(labels)
  
  # FIXED: Use simfun correctly following KEMA pattern
  # simfun should take the full label vector and return a similarity matrix
  # Compute overall similarity matrix using simfun
  M_full <- tryCatch({
    simfun(labels)
  }, error = function(e) {
    stop("Failed to compute similarity matrix using simfun. ",
         "simfun should take a label vector and return a similarity matrix. ",
         "Original error: ", e$message, call. = FALSE)
  })
  
  # Validate similarity matrix
  if (!methods::is(M_full, "Matrix") && !is.matrix(M_full)) {
    stop("simfun must return a matrix or Matrix object", call. = FALSE)
  }
  
  if (!all(dim(M_full) == c(length(labels), length(labels)))) {
    stop("simfun returned matrix with incorrect dimensions: ", 
         paste(dim(M_full), collapse = " x "), 
         ", expected: ", length(labels), " x ", length(labels), call. = FALSE)
  }
  
  # EFFICIENCY: Ensure sparse matrix for performance
  # simfun should return sparse matrices, but enforce this for efficiency
  if (!methods::is(M_full, "sparseMatrix")) {
    # Convert to sparse matrix if not already sparse
    M_full <- Matrix::Matrix(M_full, sparse = TRUE)
    
    # Inform user if conversion was needed (helps identify inefficient simfun)
    if (nrow(M_full) > 100) {  # Only warn for larger matrices where it matters
      message("Converted simfun result to sparse matrix for efficiency. ",
              "Consider using a simfun that returns sparse matrices directly.")
    }
  }
  
  # PSD-PRESERVING APPROACH: Instead of M_between = M_full - M_within (which can break PSD),
  # we construct M_within and M_between separately to ensure both are PSD
  
  # Compute within-domain similarity by masking the full similarity matrix
  # Create block diagonal mask for within-domain relationships
  domain_end_indices <- cumsum(sapply(label_list, length))
  domain_start_indices <- c(1, domain_end_indices[-length(domain_end_indices)] + 1)
  
  # Create within-domain similarity matrix (block diagonal of M_full)
  M_within <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                   dims = c(length(labels), length(labels)))
  
  for (i in seq_along(label_list)) {
    start_idx <- domain_start_indices[i]
    end_idx <- domain_end_indices[i]
    
    # Add block diagonal entries for this domain
    block_i <- seq(start_idx, end_idx)
    M_within[block_i, block_i] <- M_full[block_i, block_i]
  }
  
  # PSD-PRESERVING: Compute between-domain similarity using only off-diagonal blocks
  # This ensures M_between is PSD since it's built from PSD components
  M_between <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                    dims = c(length(labels), length(labels)))
  
  # Add off-diagonal blocks (between-domain similarities)
  for (i in seq_along(label_list)) {
    for (j in seq_along(label_list)) {
      if (i != j) {  # Only off-diagonal blocks
        start_i <- domain_start_indices[i]
        end_i <- domain_end_indices[i]
        start_j <- domain_start_indices[j]
        end_j <- domain_end_indices[j]
        
        block_i <- seq(start_i, end_i)
        block_j <- seq(start_j, end_j)
        
        # Copy the between-domain similarity block
        M_between[block_i, block_j] <- M_full[block_i, block_j]
      }
    }
  }
  
  # INTERPRETABILITY: Normalize M_within and M_between to comparable scales
  # This ensures that u behaves intuitively (u=0.5 means equal weight)
  # Compute Frobenius norms for both matrices
  M_within_norm <- Matrix::norm(M_within, "F")
  M_between_norm <- Matrix::norm(M_between, "F")
  
  # Normalize both matrices to unit Frobenius norm (if non-zero)
  # This ensures they contribute equally when u=0.5
  if (M_within_norm > 0) {
    M_within_normalized <- M_within / M_within_norm
  } else {
    M_within_normalized <- M_within  # Keep zero matrix as is
    message("M_within has zero norm - no within-domain similarities detected")
  }
  
  if (M_between_norm > 0) {
    M_between_normalized <- M_between / M_between_norm
  } else {
    M_between_normalized <- M_between  # Keep zero matrix as is
    message("M_between has zero norm - no between-domain similarities detected")
  }
  
  # PSD-PRESERVING: Combine using convex combination with normalized matrices
  # Now u=0.5 truly means equal weight between within and between domain similarities
  M <- u * M_within_normalized + (1 - u) * M_between_normalized
  
  # Report the relative scales for interpretability
  if (M_within_norm > 0 && M_between_norm > 0) {
    scale_ratio <- M_between_norm / M_within_norm
    message("Matrix normalization: M_between/M_within Frobenius norm ratio = ", 
            round(scale_ratio, 3), " (before normalization)")
  }
  
  # Add regularization
  M <- M + Matrix::Diagonal(x=rep(lambda, nrow(M)))
  
  # PSD GUARANTEE: Check and correct eigenvalues to ensure strict PSD
  # This is critical because genpca() requires PSD matrices
  min_eigenval <- tryCatch({
    # For large matrices, estimate minimum eigenvalue efficiently
    if (nrow(M) > 500) {
      # Use PRIMME to estimate smallest eigenvalue
      eig_result <- PRIMME::eigs_sym(M, NEig=1, which="SA", tol=1e-6)
      eig_result$values[1]
    } else {
      # For smaller matrices, compute exactly
      min(eigen(as.matrix(M), symmetric=TRUE, only.values=TRUE)$values)
    }
  }, error = function(e) {
    # If eigenvalue computation fails, assume we need regularization
    warning("Could not compute minimum eigenvalue, adding extra regularization: ", e$message, call. = FALSE)
    -1e-6  # Force regularization
  })
  
  # PSD CORRECTION: If minimum eigenvalue is negative or too small, add regularization
  if (min_eigenval < 1e-8) {
    correction <- abs(min_eigenval) + 1e-6
    M <- M + Matrix::Diagonal(x=rep(correction, nrow(M)))
    message("Applied PSD correction: added ", round(correction, 8), " to diagonal (min eigenvalue was ", 
            round(min_eigenval, 8), ")")
  }
  
  # Normalize by largest eigenvalue for numerical stability
  evm <- PRIMME::eigs_sym(M, NEig=1, which="LA", method='PRIMME_DEFAULT_MIN_MATVECS')
  M <- M / evm$values[1]
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  pdata <- multivarious::init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
  
  # Synchronize names between proclist and pdata
  names(proclist) <- names(pdata)
  
  # Get sample block indices using the exported function
  sample_block_idx <- block_indices(pdata)
  
  # CRITICAL: Use the 'split' method to create the list of block indices.
  # This preserves the 'start' and 'end' names within each list element,
  # which concat_pre_processors appears to require.
  sample_blocks_list <- split(sample_block_idx, row(sample_block_idx))
  
  # The names of the list elements themselves don't seem to matter, but let's
  # be consistent with the data block names for clarity.
  names(sample_blocks_list) <- names(pdata)
  
  # This call should now succeed.
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)
  
  # Compute feature block indices for multiblock_biprojector
  feat_per_block <- vapply(pdata, function(b) ncol(b$x), integer(1))
  end_idx   <- cumsum(feat_per_block)
  start_idx <- c(1L, head(end_idx, -1) + 1L)
  
  feature_block_idx <- lapply(seq_along(feat_per_block), function(i) {
    start_idx[i]:end_idx[i]
  })
  names(feature_block_idx) <- names(pdata)
  
  # Use feature block indices for multiblock_biprojector
  final_block_indices <- feature_block_idx
  
  X_block <- Matrix::bdiag(lapply(pdata, function(x) x$x))
  
  # genpca expects a regular matrix, not a sparse Matrix
  X_block_dense <- as.matrix(X_block)
  
  ret <- genpca::genpca(X_block_dense, M=row_metric_scale*M, ncomp=ncomp, preproc=multivarious::pass())

  
  multivarious::multiblock_biprojector(
    v=ret$v,
    s=ret$s,
    sdev=ret$sdev,
    preproc=proc,
    block_indices=final_block_indices,
    labels=labels,
    classes="gpca_align"
  )
}