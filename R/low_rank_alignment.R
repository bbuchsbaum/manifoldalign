#' @noRd
createSimFun <- function(S, na_value = 0) {
  # Extract the labels from the similarity matrix
  label_names <- rownames(S)
  
  # The function to return
  function(labels) {
    # Find indices of input labels in the label_names
    indices <- match(labels, label_names)
    n <- length(labels)
    
    # Start with a zero matrix; fill only the recognised labels
    M <- Matrix::Matrix(na_value, n, n, sparse = TRUE)
    
    # Find which labels are recognized (not NA)
    ok <- which(!is.na(indices))
    if (length(ok) > 0) {
      M[ok, ok] <- S[indices[ok], indices[ok]]
    }
    
    # Guarantee symmetry & zero diagonal
    M <- (M + t(M)) / 2
    diag(M) <- 0
    
    return(M)
  }
}

#' @noRd
compute_R <- function(X) {
  svd_X <- svd(t(X))
  V1 <- svd_X$v[, svd_X$d > 1]
  S1_inv_sq <- Matrix::Diagonal(x=1 / (svd_X$d[svd_X$d > 1]^2))
  R <- V1 %*% (Matrix::Diagonal(ncol(V1)) - S1_inv_sq) %*% t(V1)
  
  # Return as sparse symmetric matrix for memory efficiency
  R <- Matrix::forceSymmetric(Matrix(R, sparse = TRUE))
  return(R)
}

#' @noRd
compute_rank_matrices <- function(strata) {
  Sl <- purrr::map(strata, function(x) {
    R <- compute_R(x$x)
  })
  
  Sl
}

#' Low-rank Alignment
#'
#' Performs low-rank alignment using eigenvalue decomposition. Balances 
#' low-rank structure with similarity-based constraints. Supports 
#' semi-supervised learning with missing labels.
#'
#' @param data Input data object
#' @param y Variable name for labels (unquoted). Can contain NA values for 
#'   unlabeled samples in semi-supervised learning scenarios.
#' @param ... Additional arguments passed to specific methods
#'
#' @details
#' Low-rank alignment optimizes the objective function Z = (1-μ) * M + 2μ * L 
#' where M captures low-rank structure and L is the graph Laplacian from 
#' similarity matrix. The method balances preserving low-rank structure 
#' (μ=0) with enforcing similarity constraints (μ=1).
#' 
#' **Semi-supervised Learning Support:**
#' The algorithm handles NA labels gracefully. Unlabeled samples:
#' - Still contribute to the low-rank structure term M through their data
#' - Do not participate in the similarity constraints (L term)
#' - Receive coordinates in the joint embedding space
#' - Create isolated nodes that produce zero eigenvalues (automatically 
#'   skipped)
#'
#' @return The return value depends on the specific method. For hyperdesign 
#'   objects, returns a multiblock_biprojector object containing alignment 
#'   results, eigenvectors, preprocessing information, and metadata.
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign data
#' library(multidesign)
#' 
#' # Create synthetic data
#' set.seed(123)
#' d1 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=1, run=rep(1:5, 2)))
#' d2 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=2, run=rep(1:5, 2)))
#' d3 <- multidesign(matrix(rnorm(10*20), 10, 20), 
#'                   data.frame(y=1:10, subject=3, run=rep(1:5, 2)))
#' 
#' # Create similarity function (NA-tolerant)
#' S <- matrix(runif(10*10), 10, 10)
#' S <- abs(cor(S))
#' row.names(S) <- colnames(S) <- 1:10
#' simfun <- createSimFun(S)  # Handles NA labels automatically
#' 
#' # Create hyperdesign and run alignment
#' hd <- hyperdesign(list(d1, d2, d3))
#' result <- lowrank_align(hd, y, simfun=simfun)
#' 
#' # Semi-supervised learning with missing labels
#' d1_semi <- d1
#' d1_semi$design$y[1:3] <- NA  # Mark some samples as unlabeled
#' d2_semi <- d2
#' d2_semi$design$y[1:2] <- NA
#' hd_semi <- hyperdesign(list(d1_semi, d2_semi, d3))
#' result_semi <- lowrank_align(hd_semi, y, simfun=simfun)
#' }
#'
#' @export
lowrank_align <- function(data, y, ...) {
  UseMethod("lowrank_align")
}

#' @rdname lowrank_align
#' @method lowrank_align hyperdesign
#'
#' @param preproc Preprocessing function (default: center())
#' @param ncomp Number of components to extract (default: 2)
#' @param simfun Function to compute similarity matrix from labels. Should 
#'   handle NA labels gracefully (e.g., created with createSimFun())
#' @param mu Balance parameter between low-rank (μ=0) and similarity (μ=1) 
#'   terms (default: 0.5)
#' @param lambda Regularization parameter for glmnet. If NULL (default), uses 
#'   cross-validation to select optimal lambda via cv.glmnet. If specified, 
#'   uses the provided value directly.
#' @param scale_M Logical. If TRUE, scales M matrix to have similar eigenvalue 
#'   magnitude as L. This can improve numerical conditioning but changes the 
#'   mathematical objective. When enabled, consider adjusting the mu parameter 
#'   accordingly (default: FALSE)
#' @param n_cores Number of threads for PRIMME eigenvalue computations. If 
#'   NULL (default), uses system default. Set to 1 for reproducible results 
#'   across systems.
#'
#' @details
#' The scale_M parameter controls whether to apply eigenvalue-based scaling:
#' - scale_M = FALSE (default): Uses original formulation 
#'   Z = (1-μ) * M + 2μ * L
#' - scale_M = TRUE: Applies scaling M := M * (λ₁(L)/λ₁(M)), changing the 
#'   objective
#'
#' When scale_M = TRUE, the mu parameter no longer has its original 
#' mathematical meaning for balancing the two terms, as the relative scales 
#' have been artificially adjusted.
#'
#' For reproducibility across different systems, set n_cores = 1 to ensure 
#' deterministic results from PRIMME eigenvalue computations.
#'
#' **Handling NA Labels:**
#' Samples with NA labels are supported through the following mechanism:
#' - They contribute to the low-rank reconstruction term M = (I-R)ᵀ(I-R)
#' - They do not participate in similarity constraints (zero rows/columns in C)
#' - They create isolated nodes with zero degree, producing zero eigenvalues
#' - The algorithm automatically detects and skips these zero modes
#' - Final embedding includes coordinates for all samples (labeled and unlabeled)
#'
#' @export
lowrank_align.hyperdesign <- function(data, y, 
                                     preproc=center(), 
                                     ncomp=2,
                                     simfun,
                                     mu=.5,
                                     lambda=NULL,
                                     scale_M=FALSE,
                                     n_cores=NULL,
                                     ...) {
  
  y <- rlang::enquo(y)
  
  # Extract labels, preserving NA values for semi-supervised learning
  labels <- unlist(purrr::map(data, function(x) {
    x$design %>% select(!!y) %>% pull(!!y)
  }))
  
  # Count labeled vs unlabeled samples
  n_labeled <- sum(!is.na(labels))
  n_unlabeled <- sum(is.na(labels))
  
  if (n_labeled == 0) {
    stop("No labeled samples found. Low-rank alignment requires at least ", 
         "some labeled samples for similarity constraints.", call. = FALSE)
  }
  
  if (n_unlabeled > 0) {
    message("Semi-supervised low-rank alignment: ", n_labeled, 
            " labeled samples, ", n_unlabeled, " unlabeled samples")
  }
  
  ninstances <- length(labels)
  nsets <- length(data)
  
  # Set PRIMME threading for reproducibility
  old_threads <- getOption("PRIMME_num_threads")
  if (!is.null(n_cores)) {
    options(PRIMME_num_threads = n_cores)
    on.exit(options(PRIMME_num_threads = old_threads), add = TRUE)
  }
  
  pdata <- multivarious::init_transform(data, preproc) 
  proclist <- attr(pdata, "preproc")
  
  names(proclist) <- names(pdata)
  
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, block_indices)
  names(block_indices) <- names(pdata)
  
  Rs <- compute_rank_matrices(pdata)
  R_block <- Matrix::bdiag(Rs)  # Rs are now sparse symmetric matrices
  
  # Use simfun to create similarity matrix (should handle NA labels)
  C_block <- simfun(labels)
  C_block[C_block < 0] <- 0
  C_block <- Matrix(C_block, sparse=TRUE)
  C_block <- (C_block + t(C_block))/2
  
  diag(C_block) <- 0
  
  # CRITICAL: Identify isolated nodes (NA labels create zero-degree nodes)
  deg <- Matrix::rowSums(C_block)
  iso <- which(deg == 0)  # isolated / NA-labelled nodes
  k_skip <- length(iso)   # how many extra zero modes we expect
  
  if (k_skip > 0) {
    message("Detected ", k_skip, " isolated nodes (unlabeled samples). ",
            "Will skip corresponding zero eigenvalue modes.")
  }
  
  R_i <- (Matrix::Diagonal(nrow(R_block)) - R_block)
  M <- crossprod(R_i, R_i)
  
  Ds <- Matrix::Diagonal(x=Matrix::rowSums(C_block))
  Ls <- Ds - C_block
  
  # Optional scaling to improve numerical conditioning
  # WARNING: This changes the mathematical objective unless mu is also adjusted
  if (scale_M) {
    s_m <- PRIMME::eigs_sym(M, NEig=1, which="LA", 
                           method='PRIMME_DEFAULT_MIN_MATVECS')
    s_l <- PRIMME::eigs_sym(Ls, NEig=1, which="LA", 
                           method='PRIMME_DEFAULT_MIN_MATVECS')
    e_ratio <- s_l$values[1] / s_m$values[1] 
    M <- M * e_ratio
    
    if (getOption("lowrank_align.verbose", FALSE)) {
      message("Applied M scaling with ratio: ", round(e_ratio, 4), 
              ". This changes the objective function - consider adjusting ", 
              "mu parameter.")
    }
  }
  
  # Target equation: Z = (1-μ) * M + 2μ * L (Eq. 20)
  Z <- (1-mu) * M + (2 * mu * Ls)

  # CRITICAL: Solve eigenvalue problem, skipping zero modes from isolated nodes
  nev <- ncomp + k_skip + 1  # ask for extra eigenvectors to account for zero modes
  eigen_result <- PRIMME::eigs_sym(Z, NEig=nev, which="SA", 
                                  method='PRIMME_DEFAULT_MIN_MATVECS')
  
  # Drop exact-zero eigenvalues and select informative modes
  keep <- which(eigen_result$values > 1e-8)
  if (length(keep) < ncomp) {
    warning("Only ", length(keep), " non-zero eigenvalues found, but ", 
            ncomp, " components requested. Using ", length(keep), 
            " components.", call. = FALSE)
    ncomp <- length(keep)
  }
  
  # Select the first ncomp non-zero modes
  keep <- keep[1:ncomp]
  v <- eigen_result$vectors[, keep, drop=FALSE]
  
  Xp <- Matrix::bdiag(lapply(pdata, function(x) x$x))
  
  # Post-projection: Map block-diagonal data matrix to eigenvectors via 
  # ridge regression
  if (is.null(lambda)) {
    # Use cross-validation to select optimal lambda
    cvfit <- cv.glmnet(Xp, v, family = "mgaussian", alpha=0, 
                       intercept=FALSE)
    lambda <- cvfit$lambda.min
  } 
  
  # Fit ridge regression model
  rfit <- glmnet(Xp, v, family = "mgaussian", alpha=0, lambda=lambda, 
                 intercept=FALSE)
  
  # Extract coefficients properly: use exact=TRUE and remove intercept row
  cfs <- coef(rfit, exact=TRUE)
  cfs <- do.call(cbind, cfs)[-1, , drop=FALSE]  # Remove intercept row
  
  multivarious::multiblock_biprojector(
    v=cfs,
    s=v,
    sdev=apply(v,2,sd),
    preproc=proc,
    block_indices=block_indices,
    labels=labels,
    mu=mu,
    classes="lowrank_align"
  )
}