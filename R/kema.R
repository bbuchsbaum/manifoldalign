

#' @importFrom Matrix bdiag Diagonal rowSums sparseMatrix
#' @importFrom kernlab kernelMatrix
#' @importFrom cluster pam
#' @importFrom chk chk_number chk_range chk_logical chk_true
#' @importFrom rlang enquo !!
#' @importFrom purrr map
#' @importFrom coop cosine
#' @importFrom neighborweights graph_weights class_graph repulsion_graph binary_label_matrix adjacency
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
#' @examples
#' data1 <- list(x = matrix(rnorm(20), 10, 2))
#' data2 <- list(x = matrix(rnorm(30), 15, 2))
#' data_list <- list(data1, data2)
#' indices <- block_indices(data_list)
#' indices
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
#' @examples
#' kern <- coskern()
#' x <- c(1, 2, 3)
#' y <- c(4, 5, 6)
#' # kern(x, y)  # Compute cosine similarity
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
  y_name <- rlang::as_name(y)
  Sl <- purrr::map(strata, function(x) {
    labs <- x$design[[y_name]]
    
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
      list(
        G = neighborweights::adjacency(g),
        R = Matrix::sparseMatrix(
          i = integer(0),
          j = integer(0),
          x = numeric(0),
          dims = c(length(labs), length(labs))
        )
      )
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
  y_name <- rlang::as_name(y)
  
  if (is.null(dfun)) {
    dlabels <- lapply(strata, function(s) {
      labs <- s$design[[y_name]]
    
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
    
    # Sanitize labels before constructing between-class graph to avoid NA indices
    all_dlabels <- unlist(dlabels)
    all_dlabels[is.na(all_dlabels)] <- "UNLABELED"
    neighborweights::binary_label_matrix(all_dlabels, all_dlabels, type="d")
  } else {
    dlabels <- unlist(lapply(strata, function(s) {
      s$design[[y_name]]
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

#' Sanitize labels before computing similarity matrices
#' 
#' Ensures there are no NA values passed into simfun by replacing each missing
#' label with a unique placeholder. Using unique placeholders prevents all
#' unlabeled samples from being treated as the same class while still keeping
#' the input compatible with functions like neighborweights::binary_label_matrix,
#' which cannot accept NA values.
#'
#' @param labels Vector of labels that may contain NA entries
#' @return Vector with the same length and names as `labels` but with all
#'   missing entries replaced by unique placeholders.
#' @keywords internal
sanitize_labels_for_similarity <- function(labels) {
  if (!anyNA(labels)) {
    return(labels)
  }
  
  sanitized <- labels
  if (!is.character(sanitized)) {
    sanitized <- as.character(sanitized)
  }
  
  missing_idx <- which(is.na(sanitized))
  if (!length(missing_idx)) {
    return(sanitized)
  }
  
  sanitized[missing_idx] <- sprintf("__UNLABELED__%05d", missing_idx)
  if (!is.null(names(labels))) {
    names(sanitized) <- names(labels)
  }
  sanitized
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
#' @param use_laplacian Deprecated compatibility argument; ignored by current
#'   implementation.
#' @param solver Deprecated compatibility argument; currently accepted values are
#'   "regression" and "exact", but both route to the same original KEMA solver.
#' @param backend Backend for the original eigensolver. One of `"auto"`,
#'   `"full_exact"`, `"reduced_exact"`, or `"operator_exact"`.
#' @param backend_control Optional list controlling auto backend thresholds and
#'   fidelity checks (passed through to `kema_orig()`).
#' @param dweight Deprecated compatibility argument; ignored.
#' @param rweight Deprecated compatibility argument; ignored.
#' @param simfun Deprecated compatibility argument; ignored.
#' @param disfun Deprecated compatibility argument; ignored.
#' @param lambda Regularization parameter for matrix conditioning 
#'   (default: 0.0001)
#' @param centre_kernel Deprecated compatibility argument; ignored.
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
#' @rdname kema
#' @method kema multidesign
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
                             backend="auto",
                             backend_control=NULL,
                             dweight=.1,
                             rweight=0,
                             simfun=neighborweights::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE,
                             ...) {
  
  subject <- rlang::enquo(subject)
  y <- rlang::enquo(y)

  if (!solver %in% c("regression", "exact")) {
    stop("solver must be either 'regression' or 'exact'", call. = FALSE)
  }

  if (is.null(kernel)) {
    if (is.null(sigma)) {
      sigma <- 0.73
    }
    if (requireNamespace("kernlab", quietly = TRUE)) {
      kernel <- kernlab::rbfdot(sigma = sigma)
    } else {
      kernel <- coskern()
    }
  } else if (is.null(sigma)) {
    sigma <- 0.73
  }

  result <- kema_orig.multidesign(
    data = data,
    y = !!y,
    subject = !!subject,
    preproc = preproc,
    ncomp = ncomp,
    knn = knn,
    sigma = sigma,
    mu = u,
    kernel = kernel,
    sample_frac = sample_frac,
    lambda = lambda,
    backend = backend,
    backend_control = backend_control
  )

  # Preserve legacy class tag for downstream code using `kema()`.
  result$classes <- "kema"
  if (is.numeric(result$eigenvalues)) {
    legacy_solver <- if (identical(solver, "regression")) "exact_via_kema_orig" else "exact"
    result$eigenvalues <- list(
      values = as.numeric(result$eigenvalues),
      all_values = as.numeric(result$eigenvalues),
      solver = legacy_solver
    )
  }
  result$regression_quality <- NULL
  result$retry_info <- NULL
  result
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
#' @param use_laplacian Deprecated compatibility argument; ignored.
#' @param solver Deprecated compatibility argument; accepted values are
#'   `"regression"` and `"exact"`, but both currently route to the original
#'   KEMA solver.
#' @param dweight Deprecated compatibility argument; ignored.
#' @param rweight Deprecated compatibility argument; ignored.
#' @param simfun Deprecated compatibility argument; ignored.
#' @param disfun Deprecated compatibility argument; ignored.
#' @param lambda Regularization parameter for matrix conditioning 
#'   (default: 0.0001)
#' @param centre_kernel Deprecated compatibility argument; ignored.
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the KEMA alignment
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
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
#' hd <- structure(list(domain1 = domain1, domain2 = domain2), class = "hyperdesign")
#' 
#' # Run KEMA with default settings
#' result <- kema(hd, y = labels, ncomp = 2)
#' 
#' # Semi-supervised learning with missing labels
#' hd_semi <- hd
#' hd_semi$domain1$design$labels[1:10] <- NA  # Mark some samples as unlabeled
#' result_semi <- kema(hd_semi, y = labels, ncomp = 2)
#' 
#' # Use exact solver for highest accuracy
#' result_exact <- kema(hd, y = labels, solver = "exact", ncomp = 2)
#' 
#' # Use REKEMA for large datasets
#' result_rekema <- kema(hd, y = labels, sample_frac = 0.5, ncomp = 2)
#' }
#'
#' @details
#' `kema()` now delegates to the paper-faithful `kema_orig()` backend and solves
#' the original generalized eigenproblems from Tuia & Camps-Valls (2016),
#' including the reduced-rank REKEMA form when `sample_frac < 1`.
#'
#' Legacy extension arguments remain in the API for backward compatibility but
#' are ignored by the current implementation.
#'
#' @references
#' Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain 
#' adaptation. PLoS ONE, 11(2), e0148655.
#'
#' @rdname kema
#' @method kema hyperdesign
#' @export
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
                             backend="auto",
                             backend_control=NULL,
                             dweight=.1,
                             rweight=0,
                             simfun=neighborweights::binary_label_matrix,
                             disfun=NULL,
                             lambda=.0001,
                             centre_kernel=FALSE,
                             ...) {
  
  y <- rlang::enquo(y)

  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_range(sample_frac, c(0,1))
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_range(u, c(0,1))
  if (!solver %in% c("regression", "exact")) {
    stop("solver must be either 'regression' or 'exact'", call. = FALSE)
  }

  # Legacy compatibility: keep old arguments in signature but route to original KEMA.
  if (is.null(kernel)) {
    if (is.null(sigma)) {
      sigma <- 0.73
    }
    if (requireNamespace("kernlab", quietly = TRUE)) {
      kernel <- kernlab::rbfdot(sigma = sigma)
    } else {
      kernel <- coskern()
    }
  } else if (is.null(sigma)) {
    sigma <- 0.73
  }

  result <- kema_orig.hyperdesign(
    data = data,
    y = !!y,
    preproc = preproc,
    ncomp = ncomp,
    knn = knn,
    sigma = sigma,
    mu = u,
    kernel = kernel,
    sample_frac = sample_frac,
    lambda = lambda,
    backend = backend,
    backend_control = backend_control
  )

  # Preserve legacy class tag for downstream code using `kema()`.
  result$classes <- "kema"
  if (is.numeric(result$eigenvalues)) {
    legacy_solver <- if (identical(solver, "regression")) "exact_via_kema_orig" else "exact"
    result$eigenvalues <- list(
      values = as.numeric(result$eigenvalues),
      all_values = as.numeric(result$eigenvalues),
      solver = legacy_solver
    )
  }
  result$regression_quality <- NULL
  result$retry_info <- NULL
  result
  
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
#' @examples
#' X <- matrix(rnorm(100), 20, 5)
#' sigma <- choose_sigma(X)
#' sigma
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

#' @rdname kema
#' @method kema default
#' @export
kema.default <- function(data, ...) {
  stop("kema() requires either a hyperdesign or multidesign object. ",
       "Got: ", paste(class(data), collapse = ", "), 
       ". See ?kema for usage examples.", call. = FALSE)
}

  
