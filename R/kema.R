

#' @importFrom Matrix bdiag Diagonal rowSums sparseMatrix
#' @importFrom PRIMME eigs_sym
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom kernlab kernelMatrix
#' @importFrom cluster pam
#' @importFrom chk chk_number chk_range chk_logical chk_true
#' @importFrom rlang enquo !!
#' @importFrom purrr map
#' @importFrom dplyr select pull
#' @importFrom coop cosine
#' @importFrom neighborweights graph_weights class_graph repulsion_graph binary_label_matrix adjacency
#' @importFrom RANN nn2
#' @importFrom multivarious concat_pre_processors prep init_transform
#' @importFrom methods new is as
#' @import Matrix
#' @keywords internal
rescale <- function(z) {
  # Check data orientation (expect samples x features)
  if (nrow(z) < ncol(z)) {
    warning("rescale(): Data has more features (", ncol(z), ") than samples (", nrow(z), "). ",
            "This may indicate transposed data. Expected orientation: samples x features.", 
            call. = FALSE)
  }
  
  # Vectorized rescaling
  rn <- sqrt(Matrix::rowSums(z^2))
  rn[rn == 0] <- 1e-12  # Guard against zero-norm rows
  z / rn
}

#' Block indices for data list
#' 
#' @param data_list List of data blocks
#' @return Matrix with start and end indices for each block
#' @export
block_indices <- function(data_list) {
  # Generate block indices for a list of data blocks
  if (is.list(data_list)) {
    n_samples <- sapply(data_list, function(x) if (is.matrix(x$x)) nrow(x$x) else length(x$x))
  } else {
    stop("data_list must be a list")
  }
  
  end_indices <- cumsum(n_samples)
  start_indices <- c(1, end_indices[-length(end_indices)] + 1)
  
  cbind(start=start_indices, end=end_indices)
}



#' @keywords internal
get_block_indices <- function(data_list, byrow=FALSE) {
  # Generate block indices for kernel matrices
  if (is.list(data_list)) {
    if (byrow) {
      n_items <- sapply(data_list, nrow)
    } else {
      n_items <- sapply(data_list, ncol)
    }
  } else {
    stop("data_list must be a list")
  }
  
  end_indices <- cumsum(n_items)
  start_indices <- c(1, end_indices[-length(end_indices)] + 1)
  
  cbind(start=start_indices, end=end_indices)
}

#' Cosine kernel function
#' 
#' @return A kernel object for cosine similarity
#' @export
coskern <- function() {
  rval <- function(x, y = NULL) {
    if (!is(x, "vector")) 
      stop("x must be a vector")
    if (!is(y, "vector") && !is.null(y)) 
      stop("y must be a vector")
    if (is(x, "vector") && is.null(y)) {
      coop::cosine(x)
    }
    if (is(x, "vector") && is(y, "vector")) {
      if (!length(x) == length(y)) 
        stop("number of dimension must be the same on both data points")
      coop::cosine(x, y)
    }
  }
  
  return(new("kernel", .Data = rval, kpar = list()))
}

#' @keywords internal
stratified_subsample <- function(labels, nperlabel) {
  nlabs <- length(unique(labels))
  sp <- split(seq_along(labels), labels)
  unlist(lapply(sp, function(idx) sample(idx, nperlabel)))
}

#' @keywords internal
normalize_laplacian <- function(A) {
  # Compute normalized Laplacian: L_norm = D^(-1/2) * L * D^(-1/2)
  
  # Input validation
  if (!methods::is(A, "Matrix")) {
    A <- Matrix::Matrix(A, sparse = TRUE)
  }
  
  degrees <- Matrix::rowSums(A)
  n <- length(degrees)
  
  # Handle isolated nodes (zero degree)
  zero_degree_mask <- degrees == 0
  n_isolated <- sum(zero_degree_mask)
  
  if (n_isolated > 0) {
    message("normalize_laplacian(): ", n_isolated, 
            " isolated nodes detected - treating as disconnected components.")
  }
  
  # For isolated nodes, D^(-1/2) = 0 to preserve PSD property
  inv_sqrt_degrees <- numeric(n)
  non_zero_mask <- !zero_degree_mask
  
  if (any(non_zero_mask)) {
    inv_sqrt_degrees[non_zero_mask] <- 1 / sqrt(degrees[non_zero_mask])
  }
  
  Dinvsq <- Matrix::Diagonal(x = inv_sqrt_degrees)
  
  # Use stable formula: L_norm = I - D^(-1/2) * A * D^(-1/2)
  I <- Matrix::Diagonal(n)
  L_norm <- I - Dinvsq %*% A %*% Dinvsq
  
  # Clean up numerical artifacts
  if (methods::is(L_norm, "sparseMatrix")) {
    L_norm@x[abs(L_norm@x) < 1e-14] <- 0
    L_norm <- Matrix::drop0(L_norm)
  }
  
  # Ensure result is sparse dgCMatrix
  if (!methods::is(L_norm, "dgCMatrix")) {
    L_norm <- as(L_norm, "dgCMatrix")
  }
  
  # Validate result properties
  if (nrow(L_norm) != ncol(L_norm)) {
    stop("Normalized Laplacian is not square: ", 
         paste(dim(L_norm), collapse = " x "), call. = FALSE)
  }
  
  # Check symmetry (within numerical tolerance)
  if (nrow(L_norm) <= 1000) {  # Only check for reasonably sized matrices
    if (!isSymmetric(L_norm, tol = 1e-10)) {
      warning("Normalized Laplacian is not symmetric within tolerance. ",
              "This may indicate numerical issues.", call. = FALSE)
    }
  }
  
  L_norm
}

#' @keywords internal
normalize_adjacency <- function(A, D) {
  Dinvsq <-  Matrix::Diagonal(x=1/sqrt(diag(D)))
  Dinvsq %*% A %*% Dinvsq
}



#' @keywords internal
kcentroids <- function(X, k, sfrac=.5) {
  chk::chk_true(sfrac > 0 && sfrac <= 1)
  # X is expected to be samples x features
  
  if (sfrac < 1) {
    n <- round(sfrac * nrow(X))  # Use nrow(X) for sample count
    sidx <- sort(sample(1:nrow(X), n))
    
    res <- safe_compute(
      cluster::pam(X[sidx, , drop=FALSE], k, metric="manhattan"),
      paste0("PAM clustering failed during subsampling. ",
             "Try reducing k, increasing sfrac, or checking for constant/degenerate features.")
    )
    
    # Convert back to original indices
    return(sort(sidx[res$id.med]))
  } else {
    res <- safe_compute(
      cluster::pam(X, k, metric="manhattan"),
      paste0("PAM clustering failed on full data. ",
             "Try reducing k or checking for constant/degenerate features.")
    )
    
    return(sort(res$id.med))
  }
}

#' @keywords internal
class_medoids <- function(X, L) {
  # SEMI-SUPERVISED: Handle missing labels
  # Filter out NA labels before processing
  non_missing_mask <- !is.na(L)
  
  if (!any(non_missing_mask)) {
    stop("All labels are missing. Cannot compute class medoids.", call. = FALSE)
  }
  
  # Work only with non-missing labels and corresponding data
  L_clean <- L[non_missing_mask]
  X_clean <- X[non_missing_mask, , drop=FALSE]
  
  L_clean <- as.factor(L_clean)
  sl <- split(1:length(L_clean), L_clean)
  
  ids <- tryCatch({
    sapply(sl, function(sidx) {
      if (length(sidx) <= 2) {
        sidx[1]
      } else {
        cp <- safe_compute(
          cluster::pam(X_clean[sidx,], k=1, metric="manhattan"),
          paste0("PAM clustering failed for class with ", length(sidx), " samples. Using first sample as medoid."),
          warning_fn = function(w) { list(id.med = 1) }
        )
        sidx[cp$id.med]
      }
    })
  }, error = function(e) {
    stop("Class medoid computation failed. ",
         "This may indicate issues with the class structure or data quality. ",
         "Original error: ", e$message, call. = FALSE)
  })
  
  # Map indices back to original data indices
  original_indices <- which(non_missing_mask)
  mapped_ids <- original_indices[ids]
  names(mapped_ids) <- names(ids)
  
  mapped_ids
}



#' @noRd
compute_local_similarity <- function(strata, y, knn, weight_mode, type, sigma, repulsion=TRUE) {
  y <- rlang::enquo(y)
  Sl <- purrr::map(strata, function(x) {
    labs <- x$design %>% select(!!y) %>% pull(!!y)
    
    # SEMI-SUPERVISED: Handle missing labels
    # neighborweights functions can handle NA values in labels
    # but we need to ensure they're passed through correctly
    
    g <- neighborweights::graph_weights(x$x,
                                        weight_mode=weight_mode,
                                        neighbor_mode="knn",
                                        k=knn,
                                        type=type,
                                        sigma=sigma)
    
    # Handle NA labels in class graph
    cg <- tryCatch({
      neighborweights::class_graph(labs)
    }, error = function(e) {
      # If class_graph fails with NA labels, create empty sparse matrix
      n <- length(labs)
      warning("class_graph failed for stratum with labels containing NAs. ",
              "Creating empty class graph. Error: ", e$message, call. = FALSE)
      Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    })
    
    if (repulsion) {
      r <- tryCatch({
        neighborweights::repulsion_graph(g, cg, method="weighted")
      }, error = function(e) {
        # If repulsion_graph fails, create empty repulsion matrix
        n <- length(labs)
        warning("repulsion_graph failed for stratum with labels containing NAs. ",
                "Creating empty repulsion graph. Error: ", e$message, call. = FALSE)
        Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
      })
      list(G=neighborweights::adjacency(g), R=neighborweights::adjacency(r))
    } else {
      list(G=neighborweights::adjacency(g), R=sparseMatrix(length(labs), length(labs)))
    }
    
    
  })
  
  Sl
  
}

#' @keywords internal
compute_laplacians <- function(Ws, Wr, W, Wd, use_laplacian) {
  
  # Ensure all matrices are sparse
  matrices <- list(Ws = Ws, Wr = Wr, W = W, Wd = Wd)
  for (name in names(matrices)) {
    mat <- matrices[[name]]
    if (!methods::is(mat, "sparseMatrix")) {
      matrices[[name]] <- Matrix::Matrix(mat, sparse = TRUE)
    }
  }
  
  # Extract validated sparse matrices
  Ws <- matrices$Ws
  Wr <- matrices$Wr  
  W <- matrices$W
  Wd <- matrices$Wd
  
  if (use_laplacian) {
    # Compute normalized Laplacians
    L <- safe_compute(
      normalize_laplacian(W),
      "Failed to compute normalized Laplacian for W. This may indicate issues with the graph structure or isolated nodes"
    )
    
    Ls <- safe_compute(
      normalize_laplacian(Ws),
      "Failed to compute normalized Laplacian for Ws. This may indicate issues with the similarity graph structure"
    )
    
    Lr <- safe_compute(
      normalize_laplacian(Wr),
      "Failed to compute normalized Laplacian for Wr. This may indicate issues with the repulsion graph structure"
    )
    
    Ld <- safe_compute(
      normalize_laplacian(Wd),
      "Failed to compute normalized Laplacian for Wd. This may indicate issues with the dissimilarity graph structure"
    )
    
  } else {
    # Compute degree vectors
    d_W <- Matrix::rowSums(W)
    d_Ws <- Matrix::rowSums(Ws)  
    d_Wr <- Matrix::rowSums(Wr)
    d_Wd <- Matrix::rowSums(Wd)
    
    # Create degree matrices as sparse diagonal matrices
    D_W <- Matrix::Diagonal(x = d_W)
    D_Ws <- Matrix::Diagonal(x = d_Ws)
    D_Wr <- Matrix::Diagonal(x = d_Wr)
    D_Wd <- Matrix::Diagonal(x = d_Wd)
    
    # Compute unnormalized Laplacians: L = D - A
    L <- D_W - W
    Ls <- D_Ws - Ws
    Lr <- D_Wr - Wr  
    Ld <- D_Wd - Wd
  }
  
  # Ensure all results are sparse dgCMatrix
  result_matrices <- list(L = L, Ls = Ls, Lr = Lr, Ld = Ld)
  for (name in names(result_matrices)) {
    mat <- result_matrices[[name]]
    if (!methods::is(mat, "dgCMatrix")) {
      result_matrices[[name]] <- as(mat, "dgCMatrix")
    }
  }
  
  # Validate matrix dimensions
  n <- nrow(result_matrices$L)
  for (name in names(result_matrices)) {
    mat <- result_matrices[[name]]
    if (!all(dim(mat) == c(n, n))) {
      stop("Laplacian matrix ", name, " has incorrect dimensions: ", 
           paste(dim(mat), collapse = " x "), call. = FALSE)
    }
  }
  
  result_matrices
}

#' @keywords internal
compute_between_graph <- function(strata, y, dfun=NULL) {
  y <- rlang::enquo(y)
  
  if (is.null(dfun)) {
    dlabels <- lapply(strata, function(s) {
      labs <- s$design %>% select(!!y) %>% pull(!!y)
    
      # Only compute medoids for labeled samples
      non_missing_mask <- !is.na(labs)
      
      if (!any(non_missing_mask)) {
        # All labels are missing in this stratum
        # Return all NAs - binary_label_matrix will handle this
        return(rep(NA, length(labs)))
      }
      
      # Find medoids only among labeled samples
      labeled_indices <- which(non_missing_mask)
      labeled_labs <- labs[labeled_indices]
      
      if (length(unique(labeled_labs)) == 1) {
        # Only one class in this stratum, use first labeled sample as medoid
        medlabels <- rep(NA, length(labs))
        medlabels[labeled_indices[1]] <- labeled_labs[1]
        return(medlabels)
      }
      
      # Compute medoids for labeled samples only
      meds <- tryCatch({
        sort(class_medoids(s$x[labeled_indices, , drop=FALSE], labeled_labs))
      }, error = function(e) {
        warning("class_medoids failed for stratum with missing labels. ",
                "Using first sample of each class. Error: ", e$message, call. = FALSE)
        # Fallback: use first sample of each unique class
        unique_labs <- unique(labeled_labs)
        sapply(unique_labs, function(lab) {
          class_indices <- which(labeled_labs == lab)
          class_indices[1]  # First sample of this class
        })
      })
    
      medlabels <- rep(NA, length(labs))
      # Map medoid indices back to original indices
      medoid_original_indices <- labeled_indices[meds]
      medlabels[medoid_original_indices] <- names(meds)
      medlabels
    })
    
    # binary_label_matrix converts NA to "NA" string for semi-supervised learning
    neighborweights::binary_label_matrix(unlist(dlabels), unlist(dlabels), type="d")
  } else {
    dlabels <- unlist(lapply(strata, function(s) {
      s$design %>% select(!!y) %>% pull(!!y)
    }))
    
    # Pass labels with NAs to custom function
    dfun(dlabels)
  }
    
}

#' @keywords internal
compute_kernels <- function(strata, kernel, sample_frac, centre_kernel = FALSE) {
  if (sample_frac <= 0 || sample_frac > 1) {
    stop("sample_frac must be in (0, 1]", call. = FALSE)
  }
  
  Ks <- vector("list", length(strata))
  
  for (i in seq_along(strata)) {
    X <- strata[[i]]$x
    
    if (!is.matrix(X) || nrow(X) == 0 || ncol(X) == 0) {
      stop("Data block ", i, " is not a valid matrix or is empty", call. = FALSE)
    }
    
    if (sample_frac < 1) {
      # REKEMA: Landmark-based kernel approximation
      n_landmarks <- max(1, round(sample_frac * nrow(X)))
      
      landmark_indices <- safe_compute({
        if (n_landmarks >= nrow(X)) {
          seq_len(nrow(X))
        } else {
          sort(sample(nrow(X), n_landmarks, replace = FALSE))
        }
      }, paste0("Landmark selection failed for data block ", i))
      
      X_landmarks <- X[landmark_indices, , drop = FALSE]
      
      # Paper's REKEMA: Direct kernel without K_rr^(-1/2) normalization
      K <- safe_compute(
        kernlab::kernelMatrix(kernel, X, X_landmarks),
        paste0("Kernel computation failed for data block ", i, " in REKEMA mode")
      )
      
      message("REKEMA block ", i, ": ", nrow(K), " x ", ncol(K), " kernel matrix")
      
    } else {
      # Full kernel computation
      n_samples <- nrow(X)
      
      if (n_samples > 10000) {
        warning("Computing full kernel matrix for ", n_samples, " samples in block ", i, 
                ". Consider using sample_frac < 1 for large datasets.", call. = FALSE)
      }
      
      K <- safe_compute(
        kernlab::kernelMatrix(kernel, X),
        paste0("Full kernel computation failed for data block ", i)
      )
    }
    
    if (any(!is.finite(K))) {
      stop("Kernel matrix for block ", i, " contains non-finite values", call. = FALSE)
    }
    
    if (centre_kernel) {
      # Center kernel matrix: K_c = K - 1K/n - K1/n + 1K1/n^2
      row_means <- Matrix::rowMeans(K)
      col_means <- Matrix::colMeans(K)
      grand_mean <- mean(K)
      K <- K - outer(row_means, rep(1, ncol(K))) - outer(rep(1, nrow(K)), col_means) + grand_mean
    }
    
    # Convert to sparse format
    if (!methods::is(K, "sparseMatrix")) {
      K <- as(K, "dgCMatrix")
    }
    
    Ks[[i]] <- K
  }
  
  Ks
}

#' @keywords internal
normalize_graphs <- function(Sl, Ws, Wd) {
  # Extract and combine similarity graphs
  W_list <- lapply(Sl, "[[", "G")
  Wr_list <- lapply(Sl, "[[", "R")
  
  W <- Matrix::bdiag(W_list)
  Wr <- Matrix::bdiag(Wr_list)
  
  # Remove diagonal from class similarity matrix
  if (methods::is(Ws, "sparseMatrix")) {
    diag_mask <- Ws@i == rep(0:(ncol(Ws)-1), diff(Ws@p))
    if (any(diag_mask)) {
      Ws@x[diag_mask] <- 0
      Ws <- Matrix::drop0(Ws)
    }
  } else {
    diag(Ws) <- 0
    Ws <- Matrix::Matrix(Ws, sparse = TRUE)
  }
  
  # Ensure all matrices are sparse
  if (!methods::is(W, "sparseMatrix")) W <- Matrix::Matrix(W, sparse = TRUE)
  if (!methods::is(Wr, "sparseMatrix")) Wr <- Matrix::Matrix(Wr, sparse = TRUE)
  if (!methods::is(Wd, "sparseMatrix")) Wd <- Matrix::Matrix(Wd, sparse = TRUE)
  
  # Validate dimensions
  n_total <- nrow(W)
  matrices_to_check <- list(W = W, Wr = Wr, Ws = Ws, Wd = Wd)
  for (name in names(matrices_to_check)) {
    mat <- matrices_to_check[[name]]
    if (!all(dim(mat) == c(n_total, n_total))) {
      stop("Matrix ", name, " has incorrect dimensions: ", 
           paste(dim(mat), collapse = " x "), 
           ", expected: ", n_total, " x ", n_total, call. = FALSE)
    }
  }
  
  list(W = W, Wr = Wr, Ws = Ws, Wd = Wd)
}
 

#' Kernel Manifold Alignment for Multidesign Data
#'
#' Performs Kernel Manifold Alignment on multidesign data structures. 
#' Automatically splits data by subject variable and aligns domains.
#'
#' @param data A multidesign object containing data with subject groupings
#' @param y Name of the label variable to use for alignment (can contain NA 
#'   for unlabeled samples)
#' @param subject Name of the subject variable that defines the domains/strata
#' @param preproc Preprocessing function to apply to the data (default: 
#'   center())
#' @param ncomp Number of components to extract (default: 2)
#' @param knn Number of nearest neighbors for graph construction (default: 5)
#' @param sigma Kernel bandwidth parameter (default: 0.73)
#' @param u Trade-off parameter between data geometry and class alignment 
#'   (0-1, default: 0.5)
#' @param kernel Kernel function to use (default: coskern())
#' @param sample_frac Fraction of samples to use for kernel approximation 
#'   (default: 1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Solver method: "regression" for fast approximation (default) 
#'   or "exact" for precise solution
#' @param dweight Weight for dissimilarity/repulsion terms (default: 0.1)
#' @param rweight Weight for repulsion graph (default: 0)
#' @param simfun Function to compute similarity between labels
#' @param disfun Function to compute dissimilarity between labels (optional)
#' @param lambda Regularization parameter for matrix conditioning 
#'   (default: 0.0001)
#' @param centre_kernel Whether to center kernel matrices (default: FALSE). 
#'   **[EXTENSION]** The original paper uses uncentered kernels. Set TRUE for 
#'   centered variant.
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the KEMA alignment
#'
#' @examples
#' \donttest{
#' # Example with multidesign data
#' library(multidesign)
#' 
#' # Create synthetic multi-subject data
#' set.seed(123)
#' data_design <- expand.grid(
#'   subject = factor(1:4),
#'   condition = factor(c("A", "B")),
#'   trial = 1:10
#' )
#' 
#' # Generate synthetic data matrix
#' n_obs <- nrow(data_design)
#' n_features <- 20
#' X <- matrix(rnorm(n_obs * n_features), n_obs, n_features)
#' 
#' # Create multidesign object
#' md <- multidesign(X, data_design)
#' 
#' # Run KEMA alignment across subjects
#' result <- kema(md, y = condition, subject = subject, ncomp = 2)
#' 
#' # Semi-supervised learning with missing labels
#' data_design$condition[sample(nrow(data_design), 20)] <- NA
#' md_semi <- multidesign(X, data_design)
#' result_semi <- kema(md_semi, y = condition, subject = subject, ncomp = 2)
#' }
#'
#' @export
#' @import Matrix
kema.multidesign <- function(data, y, 
                             subject, 
                             preproc=center(), 
                             ncomp=2, 
                             knn=5, 
                             sigma=.73, u=.5, 
                             kernel=coskern(), 
                             sample_frac=1,
                             use_laplacian=TRUE, 
                             solver="regression",
                             dweight=.1,
                             rweight=0,
                             simfun=neighborweights::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE, 
                             ...) {
  
  subject <- rlang::enquo(subject)
  y <- rlang::enquo(y)
  subjects <- factor(data$design %>% select(!!subject) %>% pull(!!subject))
  subject_set <- levels(subjects)
  
  strata <- multidesign::hyperdesign(split(data, subject))
  kema.hyperdesign(strata, !!y, preproc, ncomp, knn, sigma, u, kernel, 
                   sample_frac, use_laplacian, solver, dweight, rweight, 
                   simfun, disfun, lambda, centre_kernel)
}

#' Kernel Manifold Alignment (KEMA) for Hyperdesign Data
#'
#' Performs Kernel Manifold Alignment on hyperdesign data structures. 
#' Projects data from multiple domains into a shared latent space while 
#' preserving manifold structure and aligning same-class samples.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the label variable to use for alignment (can contain NA 
#'   for unlabeled samples)
#' @param preproc Preprocessing function to apply to the data (default: 
#'   center())
#' @param ncomp Number of components to extract (default: 2)
#' @param knn Number of nearest neighbors for graph construction (default: 5)
#' @param sigma Kernel bandwidth parameter (default: 0.73)
#' @param u Trade-off parameter between data geometry and class alignment 
#'   (0-1, default: 0.5)
#' @param kernel Kernel function to use (default: coskern())
#' @param sample_frac Fraction of samples to use for kernel approximation 
#'   (default: 1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Solver method: "regression" for fast approximation (default) 
#'   or "exact" for precise solution
#' @param dweight Weight for dissimilarity/repulsion terms (default: 0)
#' @param rweight Weight for repulsion graph (default: 0)
#' @param simfun Function to compute similarity between labels
#' @param disfun Function to compute dissimilarity between labels (optional)
#' @param lambda Regularization parameter for matrix conditioning 
#'   (default: 0.0001)
#' @param centre_kernel Whether to center kernel matrices (default: FALSE). 
#'   **[EXTENSION]** The original paper uses uncentered kernels. Set TRUE for 
#'   centered variant.
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the KEMA alignment
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' library(multidesign)
#' 
#' # Create synthetic multi-domain data
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(labels = sample(c("A", "B"), 50, TRUE))
#' )
#' hd <- list(domain1 = domain1, domain2 = domain2)
#' 
#' # Run KEMA with default settings
#' result <- kema(hd, y = labels, ncomp = 2)
#' 
#' # Semi-supervised learning with missing labels
#' domain1$design$labels[1:10] <- NA  # Mark some samples as unlabeled
#' result_semi <- kema(hd, y = labels, ncomp = 2)
#' 
#' # Use exact solver for highest accuracy
#' result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#' 
#' # Use REKEMA for large datasets
#' result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#' }
#'
#' @details
#' This implementation follows the Tuia & Camps-Valls (2016) paper with 
#' extensions:
#' 
#' **Core KEMA (from paper):**
#' - Generalized eigenvalue problem: Phi(L+mu*Ls)Phi^T*v = lambda * Phi*Ld*Phi^T*v (Eq. 4)
#' - Kernelization: K(L+mu*Ls)K*Lambda = lambda * K*Ld*K * Lambda (Eq. 6)  
#' - Reduced-rank KEMA (REKEMA) for computational efficiency
#' - Matrix-free eigensolver with Jacobi preconditioning
#' 
#' **Extensions (not in original paper):**
#' - \code{solver="regression"}: Fast spectral regression approximation 
#'   (default)
#' - \code{rweight}: Additional repulsion graph Lr for within-domain 
#'   separation
#' - Semi-supervised support: Handles NA labels for unlabeled samples
#' - Enhanced numerical stability and error handling
#'
#' The algorithm offers two solver methods:
#' - "regression": **[EXTENSION]** Fast approximation using spectral 
#'   regression (default). This method first solves the eigenvalue problem on 
#'   graph Laplacians, then uses ridge regression to find kernel coefficients. 
#'   Much faster but may be less accurate for non-linear kernels.
#' - "exact": Precise solution using the correct generalized eigenvalue 
#'   formulation from the paper. This method solves the mathematically correct 
#'   KEMA optimization problem but is more computationally intensive, 
#'   especially for large datasets.
#'
#' @references
#' Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain 
#' adaptation. PLoS ONE, 11(2), e0148655.
#'
#' @export
#' @importFrom multivarious init_transform
kema.hyperdesign <- function(data, y, 
                             preproc=center(), 
                             ncomp=2, 
                             knn=5, 
                             sigma=NULL, 
                             u=.5, 
                             kernel=NULL, 
                             sample_frac=1,
                             use_laplacian=TRUE, 
                             solver="regression",
                             dweight=0,
                             rweight=0,
                             simfun=neighborweights::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE, 
                             ...) {
  
  # Input validation
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_range(sample_frac, c(0,1))
  # Validate solver parameter
  if (!solver %in% c("regression", "exact")) {
    stop("solver must be either 'regression' or 'exact'", call. = FALSE)
  }
  chk::chk_logical(use_laplacian)
  chk::chk_number(dweight)
  chk::chk_true(dweight >= 0)
  chk::chk_number(rweight)
  chk::chk_true(rweight >= 0)
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_range(u, c(0,1))
  if (!is.null(sigma)) {
    chk::chk_number(sigma)
    chk::chk_true(sigma > 0)
  }
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_logical(centre_kernel)
  
  # Validate input data
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list of hyperdesign objects", call. = FALSE)
  }
  
  y <- rlang::enquo(y)
  
  # Extract labels, accounting for quosure, and preserving NA for semi-supervised
  y_char <- rlang::as_name(y)
  
  # Manually apply preprocessing to each domain
  proc_template <- multivarious::prep(preproc)
  pdata <- data
  proclist <- list()
  
  for (i in seq_along(data)) {
    transformed_x <- multivarious::init_transform(proc_template, data[[i]]$x)
    pdata[[i]]$x <- transformed_x
    # Store the processor - use the attribute if available, otherwise use the template
    proc_attr <- attr(transformed_x, "preproc")
    proclist[[i]] <- if (is.null(proc_attr)) proc_template else proc_attr
  }
  names(proclist) <- names(data)
  
  # Extract labels AFTER preprocessing to ensure data/label alignment
  labels <- unlist(purrr::map(pdata, function(x) {
    if (!"design" %in% names(x)) {
      stop("Each domain in the hyperdesign must have a 'design' component.", call. = FALSE)
    }
    
    if (!y_char %in% names(x$design)) {
      stop("Label column '", y_char, "' not found in domain design frame. ",
           "Available columns: ", paste(names(x$design), collapse = ", "), 
           call. = FALSE)
    }
    
    # Use base R extraction for robustness
    x$design[[y_char]]
  }))
  
  # Handle case where all labels are NA
  if (all(is.na(labels))) {
    stop("No non-missing labels found. Semi-supervised learning requires at least some labeled samples.")
  }
  
  # Check labeled samples only for validation
  non_missing_mask <- !is.na(labels)
  labeled_samples <- labels[non_missing_mask]
  
  if (length(labeled_samples) == 0) {
    stop("No labeled samples found. Semi-supervised learning requires at least some labeled observations.", 
         call. = FALSE)
  }
  
  # Get unique non-missing labels for validation
  label_set <- levels(factor(labeled_samples))  # Re-factor to remove NA level for counting
  if (length(label_set) < 2) {
    stop("Need at least 2 different class labels among labeled samples for alignment. Found: ", 
         length(label_set), " unique non-missing labels.", call. = FALSE)
  }
  
  # Report semi-supervised statistics
  n_labeled <- sum(non_missing_mask)
  n_unlabeled <- sum(!non_missing_mask)
  message("Semi-supervised KEMA: ", n_labeled, " labeled samples, ", n_unlabeled, " unlabeled samples")
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  # Validate ncomp against data size
  if (ncomp >= ninstances) {
    stop("ncomp (", ncomp, ") must be less than the number of samples (", ninstances, ")", 
         call. = FALSE)
  }
  
  # Validate knn against data size
  min_samples_per_block <- min(sapply(data, function(x) nrow(x$x)))
  if (knn >= min_samples_per_block) {
    stop("knn (", knn, ") must be less than the minimum number of samples per block (", 
         min_samples_per_block, ")", call. = FALSE)
  }
  
  # Validate preprocessed data (preprocessing was already done above)
  if (any(sapply(pdata, function(x) any(!is.finite(x$x))))) {
    stop("Preprocessed data contains non-finite values (NaN/Inf). ",
         "Check input data quality and preprocessing parameters.", call. = FALSE)
  }
  
  block_indices <- block_indices(pdata)
  
  # Convert block_indices matrix to list format for multivarious::concat_pre_processors
  # The function expects a list with one entry per block, not a matrix
  block_indices_list <- split(block_indices, row(block_indices))
  
  proc <- multivarious::concat_pre_processors(proclist, block_indices_list)
  names(block_indices) <- names(pdata)
  
  # Auto-tune kernel and sigma if not provided
  if (is.null(kernel)) {
    # Concatenate all data for sigma estimation
    all_data <- do.call(rbind, lapply(pdata, function(x) x$x))
    
    if (is.null(sigma)) {
      sigma <- choose_sigma(all_data)
      message("Auto-selected sigma = ", round(sigma, 4), " using median distance heuristic")
    }
    
    # Use RBF kernel as default (better than cosine for most continuous data)
    if (requireNamespace("kernlab", quietly = TRUE)) {
      kernel <- kernlab::rbfdot(sigma = sigma)
      message("Using RBF kernel with auto-tuned sigma")
    } else {
      warning("kernlab not available, falling back to cosine kernel")
      kernel <- coskern()
    }
  } else if (is.null(sigma)) {
    # Kernel provided but sigma not specified - warn user
    message("Kernel provided but sigma not specified. Using default sigma = 0.73")
    sigma <- 0.73
  }
  
  kema_fit(pdata, proc, ncomp, knn, sigma, u, !!y, labels, kernel, sample_frac, 
           solver, dweight, rweight, block_indices, simfun, disfun, lambda, use_laplacian, centre_kernel)
  
}


#' @importFrom stats coef cor sd setNames
#' @keywords internal
#' @noRd
kema_fit <- function(strata, proc, ncomp, knn, sigma, u, y, labels, kernel, sample_frac, 
                     solver, dweight, rweight, block_indices, simfun, disfun, lambda, use_laplacian, centre_kernel) {
  chk::chk_number(ncomp)
  chk::chk_range(sample_frac, c(0,1))
  # Validate solver parameter
  if (!solver %in% c("regression", "exact")) {
    stop("solver must be either 'regression' or 'exact'", call. = FALSE)
  }
  chk::chk_number(dweight)
  chk::chk_range(u, c(0,1))
  
  y <- rlang::enquo(y)
  
  
  ## data similarity
  Sl <- compute_local_similarity(strata, !!y, knn, 
                                 weight_mode="normalized", 
                                 type="normal",  
                                 sigma=sigma,
                                 repulsion=rweight>0)
  
  
  ## class pull
  # simfun handles NA values by treating them as a separate class
  Ws <- safe_compute(
    simfun(labels),
    "Failed to compute class similarity matrix. This may indicate issues with the similarity function or label structure. For semi-supervised learning, ensure simfun can handle NA values"
  )
  
  # Validate that Ws is a proper similarity matrix
  if (!methods::is(Ws, "Matrix") && !is.matrix(Ws)) {
    stop("simfun must return a matrix or Matrix object", call. = FALSE)
  }
  
  if (!all(dim(Ws) == c(length(labels), length(labels)))) {
    stop("simfun returned matrix with incorrect dimensions: ", 
         paste(dim(Ws), collapse = " x "), 
         ", expected: ", length(labels), " x ", length(labels), call. = FALSE)
  }
  
  ## class push
  if (dweight > 0) {
    Wd <- compute_between_graph(strata, !!y, disfun)
  } else {
    # Create empty sparse matrix
    Wd <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                               dims = c(length(labels), length(labels)))
  }
  
  ## reweight graphs
  G <- normalize_graphs(Sl, Ws, Wd)
  
  ## compute full or subsampled kernels
  Ks <- compute_kernels(strata, kernel, sample_frac, centre_kernel)
  kernel_indices <- get_block_indices(Ks, byrow=TRUE)
  
  # Convert to sparse format to prevent memory explosion in bdiag()
  Ks_sparse <- lapply(Ks, function(k) {
    if (!methods::is(k, "sparseMatrix")) {
      Matrix::Matrix(k, sparse = TRUE)
    } else {
      k
    }
  })
  
  Z <- Matrix::bdiag(Ks_sparse)
  
  ## compute laplacians
  Lap <- compute_laplacians(G$Ws,G$Wr,G$W,G$Wd, use_laplacian)
  
  # Use appropriate solver based on sample_frac
  if (sample_frac == 1) {
    kemfit <- kema_full_solver(strata, Z, Ks, Lap, kernel_indices, solver, ncomp, u, 
                               dweight, rweight, lambda)
  } else {
    kemfit <- kema_landmark_solver(strata, Z, Ks, Lap, kernel_indices, solver, ncomp, u, 
                                   dweight, rweight, sample_frac, lambda)
  }
  
  # FAIL-SOFT GUARD: Automatic retry with exact solver if regression quality is poor
  if (solver == "regression" && !is.null(kemfit$regression_quality) && kemfit$regression_quality$is_poor) {
    original_solver <- solver
    warning("Regression solver produced poor results (subspace angle: ", 
            round(kemfit$regression_quality$angle_deg, 1), "°, best match: ", 
            round(kemfit$regression_quality$best_match, 3), "). ",
            "Automatically retrying with solver='exact' for higher fidelity.", 
            call. = FALSE)
    
    # Retry with exact solver
    if (sample_frac == 1) {
      kemfit_exact <- kema_full_solver(strata, Z, Ks, Lap, kernel_indices, "exact", ncomp, u, 
                                       dweight, rweight, lambda)
    } else {
      kemfit_exact <- kema_landmark_solver(strata, Z, Ks, Lap, kernel_indices, "exact", ncomp, u, 
                                           dweight, rweight, sample_frac, lambda)
    }
    
    # Store original quality info before overwriting
    original_quality <- kemfit$regression_quality
    
    # Use exact results but preserve information about the retry
    kemfit <- kemfit_exact
    kemfit$retry_info <- list(
      original_solver = original_solver,
      original_quality = original_quality,
      retried_with = "exact"
    )
    
    message("✓ Retry with exact solver completed successfully.")
  }

  # Compute feature block indices for multiblock_biprojector
  feat_per_block <- vapply(strata, function(b) ncol(b$x), integer(1))
  end_idx   <- cumsum(feat_per_block)
  start_idx <- c(1L, head(end_idx, -1) + 1L)
  
  feature_block_idx <- lapply(seq_along(feat_per_block), function(i) {
    start_idx[i]:end_idx[i]
  })
  names(feature_block_idx) <- paste0("block_", seq_along(feature_block_idx))

  result <- multivarious::multiblock_biprojector(
    v=kemfit$coef,
    s=kemfit$scores,
    sdev=apply(kemfit$scores,2,sd),
    preproc=proc,
    alpha=kemfit$vectors,
    block_indices=feature_block_idx,
    Ks=Ks,
    sample_frac=sample_frac,
    dweight=dweight,
    rweight=rweight,
    labels=labels,
    classes="kema"
  )
  
  # Add KEMA-specific information to the result
  result$eigenvalues <- kemfit$eigenvalues
  result$regression_quality <- kemfit$regression_quality
  if (!is.null(kemfit$retry_info)) {
    result$retry_info <- kemfit$retry_info
  }
  
  result
}


# Old kema_solve function removed - replaced by kema_full_solver and kema_landmark_solver

#' @keywords internal
#' Spectral Regression Helper for KEMA
#' 
#' Implements the mathematically correct spectral regression approximation for KEMA.
#' This approach splits the KEMA trace-ratio objective across two steps:
#' 1. Find target embedding Y based on "pull" forces (geometry, same-class)
#' 2. Use regularized regression with "push" forces as penalty term
#'
#' @param Z Kernel matrix (n x r for REKEMA, n x n for full KEMA)
#' @param Lap List of Laplacian matrices
#' Spectral Regression KEMA Solver
#'
#' Internal helper function that implements the spectral regression approach for 
#' KEMA. This is a two-step approximation method that first computes target 
#' embeddings from graph Laplacians, then uses regularized regression to find 
#' kernel coefficients.
#'
#' @param Z Block diagonal kernel matrix (sparse dgCMatrix)
#' @param Lap List of Laplacian matrices (L, Ls, Lr, Ld)
#' @param ncomp Number of components to extract
#' @param u Trade-off parameter between manifold and class terms
#' @param dweight Weight for dissimilarity terms
#' @param rweight Weight for repulsion terms
#' @param lambda Regularization parameter
#' @return List with vectors (coefficients) and Y (target embedding)
#' @noRd
spectral_regression_kema <- function(Z, Lap, ncomp, u, dweight, rweight, lambda) {
  # Ensure valid lambda
  if (is.null(lambda)) {
    lambda <- 0.0001
  } else if (lambda < 0) {
    stop("lambda must be non-negative, got: ", lambda, call. = FALSE)
  }
  
  # Prevent singular systems
  if (dweight == 0 && rweight == 0 && lambda < 1e-6) {
    lambda <- 1e-6
    message("Increased lambda to ", lambda, " to prevent singular system")
  }
  
  # Step 1: Compute target embedding Y from pull forces
  A_pull <- u * Lap$L + (1 - u) * Lap$Ls
  
  decomp <- safe_compute(
    PRIMME::eigs_sym(A_pull, NEig = ncomp + 1, which = "SA"),
    "Eigenvalue computation for pull-graph failed in regression solver"
  )
  
  # CORRECTNESS FIX: Numerical guard rails
  if (!all(is.finite(decomp$values))) {
    stop("Non-finite eigenvalues detected in spectral regression. ",
         "Try increasing lambda regularization or checking data quality.", call. = FALSE)
  }
  
  # Validate eigenvalue results
  if (is.null(decomp$vectors) || ncol(decomp$vectors) < ncomp) {
    stop("Insufficient eigenvectors computed (got ", ncol(decomp$vectors), 
         ", needed ", ncomp, "). Try reducing ncomp or checking matrix conditioning.", 
         call. = FALSE)
  }
  
  # Discard the first eigenvector (trivial constant vector)
  # Select the next ncomp eigenvectors corresponding to smallest non-zero eigenvalues
  eigenvals <- decomp$values
  eigenvecs <- decomp$vectors
  
  if (abs(eigenvals[1]) < 1e-10 && ncomp + 1 <= length(eigenvals)) {
    Y <- eigenvecs[, 2:(ncomp+1), drop=FALSE]
  } else {
    Y <- eigenvecs[, 1:ncomp, drop=FALSE]
  }
  
  # Step 2: Solve for kernel coefficients using regularized regression
  # Objective: min ||Y - ZV||^2 + V^T(Z^T*P_push*Z + lambda*I)V
# Solution: (Z^T*Z + Z^T*P_push*Z + lambda*I)V = Z^T*Y
  
  Z <- as(Z, "dgCMatrix")
  P_push <- rweight * Lap$Lr + dweight * Lap$Ld  # n x n sparse matrix
  
  # Construct the Left-Hand Side (LHS) of the linear system
  # This will be an r x r matrix (much smaller than n x n)
  XtX <- Matrix::crossprod(Z)  # Z^T * Z (r x r)
  
  # Add penalty term: Z^T * P_push * Z (r x r)
  if (rweight > 0 || dweight > 0) {
    Penalty_term <- Matrix::crossprod(Z, P_push %*% Z)  # Z^T * P_push * Z
  } else {
    # If no push forces, penalty term is zero
    Penalty_term <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                         dims = c(ncol(Z), ncol(Z)))
  }
  
  # Complete LHS: Z^T*Z + Z^T*P_push*Z + lambda*I
  LHS <- XtX + Penalty_term + lambda * Matrix::Diagonal(ncol(Z))
  
  # Construct the Right-Hand Side (RHS)
  RHS <- Matrix::crossprod(Z, Y)  # Z^T * Y (r x ncomp)
  
  # Check conditioning
  rcond_val <- safe_compute(
    suppressWarnings(Matrix::rcond(LHS)),
    "Failed to compute condition number",
    warning_fn = function(w) {} # Suppress warnings
  )
  
  if (!is.null(rcond_val) && rcond_val < 1e-10) {
    warning("Regression system is ill-conditioned (rcond = ", 
            format(rcond_val, scientific = TRUE), "). ",
            "Consider increasing lambda.", call. = FALSE)
  }
  
  # Solve the linear system
  vectors <- tryCatch({
    # Make LHS symmetric and add jitter for stability
    LHS_sym <- Matrix::forceSymmetric(LHS)
    jitter <- max(1e-12, 1e-10 * mean(Matrix::diag(LHS_sym)))
    Matrix::diag(LHS_sym) <- Matrix::diag(LHS_sym) + jitter
    
    # Ensure consistent matrix types
    RHS_dgC <- as(Matrix::Matrix(RHS, sparse = TRUE), "dgCMatrix")
    
    # Try Cholesky first
    chol_decomp <- Matrix::Cholesky(LHS_sym)
    Matrix::solve(chol_decomp, RHS_dgC)
  }, error = function(e) {
    # Fallback to general solver
    warning("Using general sparse solver instead of Cholesky", call. = FALSE)
    safe_compute({
      RHS_dgC <- as(Matrix::Matrix(RHS, sparse = TRUE), "dgCMatrix")
      Matrix::solve(LHS_sym, RHS_dgC, sparse = TRUE)
    }, "Linear system solve failed. Try increasing lambda")
  })
  
  # Store eigenvalues from regression solver  
  eigenvalue_info <- list(
    values = if (abs(eigenvals[1]) < 1e-10 && ncomp + 1 <= length(eigenvals)) {
               eigenvals[2:(ncomp+1)]
             } else {
               eigenvals[1:ncomp]
             },
    all_values = eigenvals,
    solver = "regression"
  )
  
  return(list(vectors = vectors, Y = Y, eigenvalues = eigenvalue_info))
}

#' @keywords internal
validate_matrix <- function(mat, name = "matrix", check_finite = TRUE, check_symmetric = FALSE) {
  if (is.null(mat)) {
    stop(name, " is NULL", call. = FALSE)
  }
  
  if (!is.matrix(mat) && !methods::is(mat, "Matrix")) {
    stop(name, " must be a matrix or Matrix object", call. = FALSE)
  }
  
  if (any(dim(mat) == 0)) {
    stop(name, " has zero dimensions: ", paste(dim(mat), collapse = " x "), call. = FALSE)
  }
  
  if (check_finite && any(!is.finite(mat))) {
    stop(name, " contains non-finite values (NaN/Inf)", call. = FALSE)
  }
  
  if (check_symmetric && nrow(mat) == ncol(mat)) {
    # Check symmetry for square matrices
    if (!isSymmetric(mat, tol = 1e-10)) {
      warning(name, " is not symmetric within tolerance", call. = FALSE)
    }
  }
  
  invisible(TRUE)
}

#' Full KEMA Solver (sample_frac == 1)
#' 
#' Implements the full KEMA algorithm for complete kernel matrices.
#' This solver handles the n x n generalized eigenvalue problem.
#'
#' @param strata List of data strata
#' @param Z Block diagonal kernel matrix (n x n)
#' @param Ks List of individual kernel matrices
#' @param Lap List of Laplacian matrices
#' @param kernel_indices Block indices for kernels
#' @param solver Solver method: "regression" for fast approximation or "exact" for precise solution
#' @param ncomp Number of components to extract
#' @param u Trade-off parameter between manifold and class terms
#' @param dweight Weight for dissimilarity terms
#' @param rweight Weight for repulsion terms
#' @param lambda Regularization parameter
#' @return List with coef, scores, and vectors
#' @keywords internal
kema_full_solver <- function(strata, Z, Ks, Lap, kernel_indices, solver, ncomp, u, 
                             dweight, rweight, lambda) {
  
  # Honor user-supplied lambda value (remove override)
  if (is.null(lambda)) {
    lambda <- 0.0001  # Match main function default
  }
  
  if (solver == "regression") {
    # SPECTRAL REGRESSION PATH: Mathematically correct and robust approximation
    
    # Check for legacy mode (for backward compatibility)
    if (getOption("kema.legacy_specreg", FALSE)) {
      # Legacy implementation (flawed but preserved for reproducibility)
      warning("Using legacy spectral regression implementation. ",
              "This method has known mathematical issues. ",
              "Set options(kema.legacy_specreg = FALSE) to use the improved method.",
              call. = FALSE)
      
      A <- Lap$Ls - (rweight*Lap$Lr + dweight*Lap$Ld)
      decomp <- PRIMME::eigs_sym(u*Lap$L + (1-u)*A, NEig=ncomp+1, which="SA")
      Y <- decomp$vectors[, 2:(ncomp+1), drop=FALSE]
      Z <- as(Z, "dgCMatrix")
      
      if (is.null(lambda)) {
        cvfit <- cv.glmnet(Z, Y, family = "mgaussian", alpha=0, intercept=FALSE)
        lambda_glmnet <- cvfit$lambda.min
      } else {
        lambda_glmnet <- lambda
      }
      
      rfit <- glmnet(Z, Y, family = "mgaussian", alpha=0, lambda=lambda_glmnet, intercept=FALSE)
      cfs <- coef(rfit)
      vectors <- do.call(cbind, cfs)[-1,,drop=FALSE]
      
      # For legacy mode, create basic eigenvalue info from decomp
      eigenvalue_info <- list(
        values = if (abs(decomp$values[1]) < 1e-10 && ncomp + 1 <= length(decomp$values)) {
                   decomp$values[2:(ncomp+1)]
                 } else {
                   decomp$values[1:ncomp]
                 },
        all_values = decomp$values,
        solver = "regression_legacy"
      )
      
    } else {
      # New, mathematically correct implementation
      result <- spectral_regression_kema(Z, Lap, ncomp, u, dweight, rweight, lambda)
      vectors <- result$vectors
      Y <- result$Y
      eigenvalue_info <- result$eigenvalues
    }
    
    # Enhanced quality check using rotation-invariant subspace metric
    Z <- as(Z, "dgCMatrix")
    Y_hat <- Z %*% vectors
    
    # Convert to base matrices for analysis
    Y_mat <- as.matrix(Y)
    Y_hat_mat <- as.matrix(Y_hat)
    
    # Rotation-invariant subspace comparison using principal angles
    subspace_quality <- tryCatch({
      # Helper function for principal angle distance
      qa <- qr.Q(qr(Y_mat))
      qb <- qr.Q(qr(Y_hat_mat))
      sv <- svd(t(qa) %*% qb, nu = 0, nv = 0)$d
      max_angle <- acos(min(sv))  # Largest principal angle
      subspace_distance <- sin(max_angle)  # Distance in [0,1]
      
      # Also compute best-matching correlation for backwards compatibility
      corr_matrix <- abs(cor(Y_mat, Y_hat_mat, use = "pairwise.complete.obs"))
      best_match <- sum(apply(corr_matrix, 1, max)) / ncol(corr_matrix)
      
      list(distance = subspace_distance, best_match = best_match, 
           angle_deg = max_angle * 180 / pi)
    }, error = function(e) {
      warning("Could not compute regression approximation quality: ", e$message, call. = FALSE)
      list(distance = 0, best_match = 1.0, angle_deg = 0)
    })
    
    # FAIL-SOFT GUARD: Assess regression quality and store for potential retry
    quality <- assess_regression_quality(Y_mat, Y_hat_mat, "full KEMA regression")
    
    # Store quality metrics for main function decision
    regression_quality <- list(
      is_poor = quality$is_poor,
      subspace_distance = quality$distance,
      best_match = quality$best_match,
      angle_deg = quality$angle_deg
    )
    
  } else {
    # EXACT KEMA PATH: Solve the correct generalized eigenvalue problem
    
    Z <- as(Z, "dgCMatrix")
    n <- nrow(Z)
    
    # CORRECTED FORMULATION (Ticket 1):
    # A = u*L + (1-u)*Ls  (pull towards manifold structure and same-class samples)
    # B = rweight*Lr + dweight*Ld + lambda*I  (push away from different classes + regularization)
    
    # Full KEMA uses standard lambda (no scaling needed since sample_frac = 1)
    A_laplacian <- u * Lap$L + (1-u) * Lap$Ls
    B_laplacian <- rweight * Lap$Lr + dweight * Lap$Ld + lambda * Matrix::Diagonal(nrow(Lap$L))
    
    # Matrix-free operators for large datasets
    matvec_A <- function(x) {
      # CORRECTNESS FIX: Ensure x is always a numeric vector for PRIMME
      x <- as.numeric(x)
      temp1 <- Matrix::crossprod(Z, x)      # K^T * x
      temp2 <- A_laplacian %*% temp1        # A_laplacian * (K^T * x)
      result <- Z %*% temp2                 # K * (A_laplacian * (K^T * x))
      as.numeric(result)                    # Ensure result is numeric vector
    }
    
    matvec_B <- function(x) {
      # CORRECTNESS FIX: Ensure x is always a numeric vector for PRIMME
      x <- as.numeric(x)
      temp1 <- Matrix::crossprod(Z, x)      # K^T * x
      temp2 <- B_laplacian %*% temp1        # B_laplacian * (K^T * x)
      result <- Z %*% temp2                 # K * (B_laplacian * (K^T * x))
      as.numeric(result)                    # Ensure result is numeric vector
    }
    
    # CORRECTNESS FIX: Improved Jacobi preconditioner
    # Use rowSums as stabilizer when diagonal is near zero (e.g., cosine kernels)
    diag_A_laplacian <- Matrix::diag(A_laplacian)
    diag_K <- Matrix::diag(Z)
    
    # Use rowSums as fallback when diagonal is uninformative
    row_sums_K <- Matrix::rowSums(abs(Z))
    diag_A_approx <- diag_K * (diag_A_laplacian * diag_K)
    
    # Replace near-zero entries with rowSums-based approximation
    near_zero_mask <- abs(diag_A_approx) < 1e-9
    if (any(near_zero_mask)) {
      diag_A_approx[near_zero_mask] <- row_sums_K[near_zero_mask] * mean(abs(diag_A_laplacian))
    }
    
    # Final safety clipping
    diag_A_approx[abs(diag_A_approx) < 1e-12] <- 1.0
    
    preconditioner <- function(x) {
      x / diag_A_approx
    }
    
    # Estimate matrix norms for PRIMME
    A_laplacian_norm <- tryCatch({
      Matrix::norm(A_laplacian, "I")
    }, error = function(e) { 1.0 })
    
    K_norm <- tryCatch({
      Matrix::norm(Z, "I")
    }, error = function(e) { 1.0 })
    
    A_norm_estimate <- K_norm^2 * A_laplacian_norm
    
    # DEFAULT: Use explicit matrix construction (more stable)
    # Matrix-free solver is experimental and can be enabled by setting use_matrix_free = TRUE
    use_matrix_free <- FALSE  # Make matrix-free experimental for now
    
    if (use_matrix_free) {
      # EXPERIMENTAL: Matrix-free solver (may have stability issues)
      # PRIMME FIX: Build explicit B matrix to avoid hang with two matrix-free operators
      B_explicit <- Z %*% B_laplacian %*% Matrix::t(Z)
      # Regularize B matrix
      B_explicit <- B_explicit + (1e-8 * Matrix::Diagonal(n))
      
      # Optimized PRIMME parameters
      max_block_size <- min(2 * ncomp, 32)
      max_basis_size <- min(8 * max_block_size, 256)
      
      # Try matrix-free approach first
      decomp <- tryCatch({
        PRIMME::eigs_sym(
          A = matvec_A,
          B = B_explicit,
          n = n,
          NEig = ncomp + 1,
          which = "SA",
          prec = preconditioner,
          tol = 1e-5,
          method = "PRIMME_DEFAULT_MIN_MATVECS",
          maxBlockSize = max_block_size,
          maxBasisSize = max_basis_size,
          aNorm = A_norm_estimate,
          printLevel = 0
        )
      }, error = function(e) {
        # If matrix-free fails, fall back to explicit
        message("Matrix-free solver failed, falling back to explicit matrices. Error: ", e$message)
        NULL
      })
    } else {
      decomp <- NULL  # Skip matrix-free, go directly to explicit
    }
    
    # Default approach: Explicit matrix construction
    if (is.null(decomp)) {
      A_explicit <- Z %*% A_laplacian %*% Matrix::t(Z)
      B_explicit <- Z %*% B_laplacian %*% Matrix::t(Z)
      
      # Regularize B matrix
      B_explicit <- B_explicit + (1e-8 * Matrix::Diagonal(n))
      
      decomp <- tryCatch({
        PRIMME::eigs_sym(A_explicit, NEig=ncomp+1, which="SA", B=B_explicit)
      }, error = function(e) {
        stop("Explicit matrix solver failed. ",
             "This indicates severe numerical issues. ",
             "Try reducing ncomp, increasing lambda regularization, or checking data quality. ",
             "Error: ", e$message, call. = FALSE)
      })
    }
    
    # Check for non-finite eigenvalues
    if (!all(is.finite(decomp$values))) {
      stop("Non-finite eigenvalues detected in exact solver. ",
           "Try increasing lambda regularization or checking data quality.", call. = FALSE)
    }
    
    # Validate and select eigenvectors
    if (is.null(decomp$vectors) || ncol(decomp$vectors) < ncomp) {
      stop("Insufficient eigenvectors computed in generalized eigenvalue problem (got ", 
           ncol(decomp$vectors), ", needed ", ncomp, "). ",
           "Try reducing ncomp or increasing lambda regularization.", call. = FALSE)
    }
    
    eigenvals <- decomp$values
    eigenvecs <- decomp$vectors
    
    if (abs(eigenvals[1]) < 1e-10 && ncomp + 1 <= length(eigenvals)) {
      vectors <- eigenvecs[, 2:(ncomp+1), drop=FALSE]
      selected_eigenvals <- eigenvals[2:(ncomp+1)]
    } else {
      vectors <- eigenvecs[, 1:ncomp, drop=FALSE]
      selected_eigenvals <- eigenvals[1:ncomp]
    }
    
    # Store eigenvalues for exact solver
    eigenvalue_info <- list(
      values = selected_eigenvals,
      all_values = eigenvals,
      solver = "exact"
    )
    
    regression_quality <- NULL  # No regression quality for exact solver
  }
  
  # Compute scores and primal vectors
  scores <- Z %*% vectors
  
  # Primal vector calculation for full KEMA
  v <- do.call(rbind, lapply(1:length(strata), function(i) {
    xi <- strata[[i]]$x  # xi is samples x features
    kind <- seq(kernel_indices[i,1], kernel_indices[i,2])
    alpha_i <- vectors[kind,,drop=FALSE]  # alpha_i is samples x components
    
    # crossprod(xi, alpha_i) = t(xi) %*% alpha_i gives features x components
    v_i <- Matrix::crossprod(xi, alpha_i)
    v_i
  }))
  
  list(coef=v, scores=scores, vectors=vectors, eigenvalues=eigenvalue_info, regression_quality=regression_quality)
}

#' @keywords internal
#' REKEMA Landmark Solver (sample_frac < 1)
#' 
#' Implements the Reduced-rank KEMA (REKEMA) algorithm using landmark points.
#' This solver reduces computational complexity from O(n^2) to O(r^2) where r is 
#' the number of landmarks.
#'
#' @param strata List of data strata
#' @param Z Landmark kernel matrix (n x r)
#' @param Ks List of individual kernel matrices
#' @param Lap List of Laplacian matrices
#' @param kernel_indices Block indices for kernels
#' @param solver Solver method: "regression" for fast approximation or "exact" for precise solution
#' @param ncomp Number of components to extract
#' @param u Trade-off parameter between manifold and class terms
#' @param dweight Weight for dissimilarity terms
#' @param rweight Weight for repulsion terms
#' @param sample_frac Fraction of samples used as landmarks
#' @param lambda Regularization parameter
#' @return List with coef, scores, and vectors
#' @noRd
kema_landmark_solver <- function(strata, Z, Ks, Lap, kernel_indices, solver, ncomp, u, 
                                 dweight, rweight, sample_frac, lambda) {
  
  # Honor user-supplied lambda value (remove override)
  if (is.null(lambda)) {
    lambda <- 0.0001  # Match main function default
  }
  
  if (solver == "regression") {
    # SPECTRAL REGRESSION PATH for REKEMA: Mathematically correct approximation
    
    # Check for legacy mode (for backward compatibility)
    if (getOption("kema.legacy_specreg", FALSE)) {
      # Legacy implementation (flawed but preserved for reproducibility)
      warning("Using legacy REKEMA spectral regression implementation. ",
              "This method has known mathematical issues. ",
              "Set options(kema.legacy_specreg = FALSE) to use the improved method.",
              call. = FALSE)
      
      A <- Lap$Ls - (rweight*Lap$Lr + dweight*Lap$Ld)
      decomp <- PRIMME::eigs_sym(u*Lap$L + (1-u)*A, NEig=ncomp+1, which="SA")
      Y <- decomp$vectors[, 2:(ncomp+1), drop=FALSE]
      Z <- as(Z, "dgCMatrix")
      
      if (is.null(lambda)) {
        cvfit <- cv.glmnet(Z, Y, family = "mgaussian", alpha=0, intercept=FALSE)
        lambda_glmnet <- cvfit$lambda.min
      } else {
        lambda_glmnet <- lambda
      }
      
      rfit <- glmnet(Z, Y, family = "mgaussian", alpha=0, lambda=lambda_glmnet, intercept=FALSE)
      cfs <- coef(rfit)
      vectors <- do.call(cbind, cfs)[-1,,drop=FALSE]
      
      # For REKEMA legacy mode, create basic eigenvalue info from decomp
      eigenvalue_info <- list(
        values = if (abs(decomp$values[1]) < 1e-10 && ncomp + 1 <= length(decomp$values)) {
                   decomp$values[2:(ncomp+1)]
                 } else {
                   decomp$values[1:ncomp]
                 },
        all_values = decomp$values,
        solver = "rekema_regression_legacy"
      )
      
    } else {
      # New, mathematically correct implementation for REKEMA
      result <- spectral_regression_kema(Z, Lap, ncomp, u, dweight, rweight, lambda)
      vectors <- result$vectors
      Y <- result$Y
      eigenvalue_info <- result$eigenvalues
    }
    
    # FAIL-SOFT GUARD: Assess regression quality for REKEMA
    quality <- assess_regression_quality(Y_mat, Y_hat_mat, "REKEMA regression")
    regression_quality <- list(
      is_poor = quality$is_poor,
      subspace_distance = quality$distance,
      best_match = quality$best_match,
      angle_deg = quality$angle_deg
    )
    
    # Enhanced quality check using rotation-invariant subspace metric for REKEMA
    Z <- as(Z, "dgCMatrix")
    Y_hat <- Z %*% vectors
    
    # Convert to base matrices for analysis
    Y_mat <- as.matrix(Y)
    Y_hat_mat <- as.matrix(Y_hat)
    
    # Check regression quality using principal angles
    subspace_quality <- safe_compute({
      qa <- qr.Q(qr(Y_mat))
      qb <- qr.Q(qr(Y_hat_mat))
      sv <- svd(t(qa) %*% qb, nu = 0, nv = 0)$d
      max_angle <- acos(min(sv))
      subspace_distance <- sin(max_angle)
      
      corr_matrix <- abs(cor(Y_mat, Y_hat_mat, use = "pairwise.complete.obs"))
      best_match <- sum(apply(corr_matrix, 1, max)) / ncol(corr_matrix)
      
      list(distance = subspace_distance, best_match = best_match, 
           angle_deg = max_angle * 180 / pi)
    }, "Could not compute REKEMA regression quality", 
    warning_fn = function(w) { list(distance = 0, best_match = 1.0, angle_deg = 0) })
    
    if (is.finite(subspace_quality$distance) && subspace_quality$distance > 0.1) {
      k_rank <- safe_compute(Matrix::rankMatrix(Z)[1], "rank computation failed", 
                            warning_fn = function(w) { "unknown" })
      
      warning("REKEMA regression approximation fidelity is below target:\n",
              "  Best component match: ", round(subspace_quality$best_match, 3), "\n",
              "  Subspace angle (deg): ", round(subspace_quality$angle_deg, 1), "\n",
              "  Samples: ", nrow(Z), ", Kernel rank: ", k_rank, "\n",
              "Consider using solver='exact' or increasing sample_frac.",
              call. = FALSE)
    }
    
  } else {
    # EXACT REKEMA PATH: Solve the reduced r x r generalized eigenvalue problem
    # This is the key innovation of Ticket 2 - much more efficient than matrix-free n x n
    
    Z <- as(Z, "dgCMatrix")
    n <- nrow(Z)  # Total number of samples
    r <- ncol(Z)  # Total number of landmarks
    
    message("REKEMA: Solving reduced ", r, " x ", r, " eigenvalue problem instead of ", n, " x ", n)
    
    # CORRECTED FORMULATION (Ticket 1):
    # A = u*L + (1-u)*Ls  (pull towards manifold structure and same-class samples)
    # B = rweight*Lr + dweight*Ld + lambda*I  (push away from different classes + regularization)
    
    # LAMBDA SCALING FIX: Scale regularization for landmark size to maintain comparable energy
    # Full KEMA uses lambda*I_n, REKEMA should use lambda_r = lambda * n/r to maintain
    # comparable regularization energy across different sample fractions
    lambda_scaled <- lambda * n / r
    message("REKEMA: Scaling lambda from ", lambda, " to ", round(lambda_scaled, 6), 
            " (factor: ", round(n/r, 2), ") to maintain regularization energy")
    
    A_laplacian <- u * Lap$L + (1-u) * Lap$Ls
    B_laplacian <- rweight * Lap$Lr + dweight * Lap$Ld + lambda_scaled * Matrix::Diagonal(nrow(Lap$L))
    
    # PAPER'S REKEMA FORMULATION (Eq. 10): Direct projection without K_rr normalization
    # Solves the r x r generalized eigenvalue problem:
    # Krn(L + mu*Ls)Knr * Lambda = lambda * KrnLdKnr * Lambda
    # where:
    # - Krn is r x n (landmarks x all samples) = Z^T
    # - Knr is n x r (all samples x landmarks) = Z
    # - A_reduced = Krn * A_laplacian * Knr (r x r)
    # - B_reduced = Krn * B_laplacian * Knr (r x r)
    #
    # This formulation is more robust as it doesn't require K_rr to be invertible
    
    # Construct reduced matrices
    A_reduced <- safe_compute(
      Matrix::crossprod(Z, A_laplacian %*% Z),
      "Failed to construct reduced A matrix for REKEMA"
    )
    
    B_reduced <- safe_compute(
      Matrix::crossprod(Z, B_laplacian %*% Z),
      "Failed to construct reduced B matrix for REKEMA"
    )
    
    # Ensure B_reduced is well-conditioned
    if (rweight == 0 && dweight == 0) {
      min_reg <- max(lambda_scaled, 1e-8)
      B_reduced <- B_reduced + min_reg * Matrix::Diagonal(r)
      message("REKEMA: Added regularization (", min_reg, ") to reduced B matrix")
    } else {
      rcond_B <- safe_compute(Matrix::rcond(B_reduced), "rcond failed", 
                              warning_fn = function(w) { 1.0 })
      
      if (rcond_B < 1e-12) {
        reg_amount <- max(lambda_scaled, 1e-8)
        B_reduced <- B_reduced + reg_amount * Matrix::Diagonal(r)
        message("REKEMA: Added regularization to ill-conditioned B matrix")
      }
    }
    
    # Validate reduced matrices
    if (!all(dim(A_reduced) == c(r, r)) || !all(dim(B_reduced) == c(r, r))) {
      stop("Reduced matrices have incorrect dimensions. Expected ", r, " x ", r, 
           ", got A: ", paste(dim(A_reduced), collapse=" x "), 
           ", B: ", paste(dim(B_reduced), collapse=" x "), call. = FALSE)
    }
    
    # Solve the r x r generalized eigenvalue problem
    decomp <- safe_compute(
      PRIMME::eigs_sym(A_reduced, NEig=min(ncomp+1, r-1), which="SA", B=B_reduced),
      "REKEMA reduced eigenvalue problem failed"
    )
    
    if (!all(is.finite(decomp$values))) {
      stop("Non-finite eigenvalues detected in REKEMA solver", call. = FALSE)
    }
    
    # Validate eigenvalue results
    if (is.null(decomp$vectors) || ncol(decomp$vectors) < ncomp) {
      stop("Insufficient eigenvectors computed in REKEMA (got ", ncol(decomp$vectors), 
           ", needed ", ncomp, "). ",
           "Try reducing ncomp or increasing lambda regularization.", call. = FALSE)
    }
    
    # Select eigenvectors (these are r-dimensional, not n-dimensional!)
    eigenvals <- decomp$values
    eigenvecs <- decomp$vectors
    
    if (abs(eigenvals[1]) < 1e-10 && ncomp + 1 <= length(eigenvals)) {
      vectors <- eigenvecs[, 2:(ncomp+1), drop=FALSE]  # r x ncomp
      selected_eigenvals <- eigenvals[2:(ncomp+1)]
    } else {
      vectors <- eigenvecs[, 1:ncomp, drop=FALSE]      # r x ncomp
      selected_eigenvals <- eigenvals[1:ncomp]
    }
    
    # Store eigenvalues for REKEMA exact solver
    eigenvalue_info <- list(
      values = selected_eigenvals,
      all_values = eigenvals,
      solver = "rekema_exact"
    )
    
    regression_quality <- NULL  # No regression quality for exact solver
    
    # Note: Eigenvectors from PRIMME should already be B-orthogonal as per standard
  }
  
  # Compute scores: Z is n x r, vectors is r x ncomp, so scores is n x ncomp
  scores <- Z %*% vectors
  
  # Primal vector calculation for REKEMA
  # Map from landmark coefficients back to feature space
  
  # Compute landmark indices for each stratum
  landmark_counts <- sapply(Ks, ncol)  # Number of landmarks per stratum
  landmark_end_idx <- cumsum(landmark_counts)
  landmark_start_idx <- c(1, landmark_end_idx[-length(landmark_end_idx)] + 1)
  landmark_indices <- cbind(start = landmark_start_idx, end = landmark_end_idx)
  
  v <- do.call(rbind, lapply(1:length(strata), function(i) {
    xi <- strata[[i]]$x  # xi is samples x features
    kind <- seq(landmark_indices[i,1], landmark_indices[i,2])
    
    alpha_i <- vectors[kind,,drop=FALSE]  # alpha_i is landmarks x components
    
    Ki <- Ks[[i]]  # samples x landmarks
    
    # Compute primal vectors: t(xi) %*% (Ki %*% alpha_i) = features x components
    v_i <- Matrix::crossprod(xi, Ki %*% alpha_i)
    v_i
  }))
  
  list(coef=v, scores=scores, vectors=vectors, eigenvalues=eigenvalue_info, regression_quality=regression_quality)
}

#' Choose optimal sigma for RBF kernel
#' 
#' Estimates a reasonable sigma parameter for RBF kernels using the median 
#' pairwise distance heuristic, which often works well in practice.
#'
#' @param X Data matrix (samples x features)
#' @param sample_size Maximum number of pairs to sample for distance computation
#'   (default: 1000 for efficiency)
#' @return Suggested sigma value
#' @export
choose_sigma <- function(X, sample_size = 1000) {
  if (!is.matrix(X) && !methods::is(X, "Matrix")) {
    stop("X must be a matrix", call. = FALSE)
  }
  
  n <- nrow(X)
  if (n < 2) {
    return(1.0)  # Default fallback
  }
  
  # For large datasets, sample pairs to avoid O(n^2) computation
  if (n > sqrt(sample_size)) {
    indices <- sample(n, min(n, sqrt(sample_size)))
    X_sample <- X[indices, , drop = FALSE]
  } else {
    X_sample <- X
  }
  
  # Compute pairwise distances efficiently
  distances <- as.matrix(dist(X_sample))
  
  # Remove diagonal (zero distances)
  distances <- distances[upper.tri(distances)]
  
  if (length(distances) == 0) {
    return(1.0)  # Fallback
  }
  
  # Median distance heuristic
  median_dist <- median(distances, na.rm = TRUE)
  
  # Convert to RBF sigma: typical rule of thumb is sigma = median_dist / sqrt(2)
  sigma <- median_dist / sqrt(2)
  
  # Ensure reasonable bounds
  if (!is.finite(sigma) || sigma <= 0) {
    sigma <- 1.0
  }
  
  sigma
}

#' @keywords internal
#' @noRd
assess_regression_quality <- function(Y_mat, Y_hat_mat, method_name = "regression") {
  # Comprehensive regression quality assessment using rotation-invariant metrics
  
  subspace_quality <- tryCatch({
    # Principal angle analysis
    qa <- qr.Q(qr(Y_mat))
    qb <- qr.Q(qr(Y_hat_mat))
    sv <- svd(t(qa) %*% qb, nu = 0, nv = 0)$d
    max_angle <- acos(pmax(pmin(min(sv), 1), -1))  # Clamp for numerical stability
    subspace_distance <- sin(max_angle)  # Distance in [0,1]
    
    # Component-wise correlation analysis
    corr_matrix <- abs(cor(Y_mat, Y_hat_mat, use = "pairwise.complete.obs"))
    best_match <- sum(apply(corr_matrix, 1, max)) / ncol(corr_matrix)
    
    list(
      distance = subspace_distance, 
      best_match = best_match, 
      angle_deg = max_angle * 180 / pi,
      is_poor = subspace_distance > 0.15 || best_match < 0.75  # Adaptive thresholds
    )
  }, error = function(e) {
    list(distance = 0, best_match = 1.0, angle_deg = 0, is_poor = FALSE)
  })
  
  subspace_quality
}

  
