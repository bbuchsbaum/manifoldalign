#' Graph Alignment by Spectral Corresponding Functions (GRASP) - Corrected Version
#'
#' Performs GRASP alignment on hyperdesign data structures. 
#' Projects graph data from multiple domains into aligned spectral spaces.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param preproc Preprocessing function to apply to the data (default: center)
#' @param ncomp Number of spectral components to extract (default: 30)
#' @param q_descriptors Number of diffusion-time descriptors (default: 100)
#' @param sigma Diffusion parameter for descriptor computation (default: 0.73)
#' @param lambda Regularization parameter for base alignment (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment method: "linear" for exact assignment (default)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the GRASP alignment
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' library(tibble)
#' 
#' # Create synthetic graph domains (node coordinates for graph construction)
#' set.seed(123)
#' X1 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 1
#' X2 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 2
#' 
#' # Create design data frames with node indices (GRASP doesn't use features)
#' design1 <- tibble(node_id = 1:50)
#' design2 <- tibble(node_id = 1:50)
#' 
#' # Create multidesign objects
#' md1 <- multidesign(X1, design1)
#' md2 <- multidesign(X2, design2)
#' 
#' # Create hyperdesign from multidesign objects
#' hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
#' 
#' # Run GRASP alignment (purely spectral, uses graph structure from X)
#' result <- grasp(hd, ncomp = 20, q_descriptors = 50)
#' }
#'
#' @export
#' @importFrom multivarious init_transform prep concat_pre_processors
#' @importFrom manifoldalign block_indices
#' @importFrom stats dist
grasp.hyperdesign <- function(data, 
                             preproc = multivarious::center(), 
                             ncomp = 30,
                             q_descriptors = 100,
                             sigma = 0.73,
                             lambda = 0.1,
                             use_laplacian = TRUE,
                             solver = "linear",
                             ...) {
  
  # Input validation (following package patterns + your suggestions)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_logical(use_laplacian)
  
  # Validate solver parameter
  if (!solver %in% c("linear", "hungarian", "auction")) {
    stop("solver must be 'linear', 'hungarian', or 'auction'", call. = FALSE)
  }
  
  # Validate input data - should be a hyperdesign object containing 2 domains
  chk::chk_s3_class(data, "hyperdesign")
  
  if (!is.list(data) || length(data) == 0) {
    stop("hyperdesign data must contain domains", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("GRASP currently supports exactly 2 domains. For multi-graph ", 
         "alignment with 3+ domains, use grasp_multiset()", call. = FALSE)
  }
  
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
  
  # Block indices computation
  block_indices <- block_indices(pdata)
  
  # This is the correct concatenated preprocessor
  full_proc <- multivarious::concat_pre_processors(proclist, 
                                             split(block_indices, row(block_indices)))
  names(block_indices) <- names(pdata)
  
  # Call GRASP fitting function
  grasp_fit(pdata, full_proc, ncomp, q_descriptors, sigma, lambda, 
             use_laplacian, solver, block_indices)
}

#' @keywords internal
grasp_fit <- function(strata, proc, ncomp, q_descriptors, sigma, lambda, 
                     use_laplacian, solver, block_indices) {
  
  # Block A: Spectral basis construction
  bases <- compute_grasp_basis(strata, ncomp, use_laplacian)
  
  # Block B: Multi-scale descriptors  
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)
  
  # Block C: Base alignment
  alignment_result <- align_grasp_bases(bases[[1]], bases[[2]], 
                                       descriptors[[1]], descriptors[[2]], lambda)
  
  # Blocks D & E: Functional map and assignment
  assignment_result <- compute_grasp_assignment(bases[[1]], bases[[2]], 
                                               descriptors[[1]], descriptors[[2]],
                                               alignment_result$rotation,
                                               distance_method = "cosine",
                                               solver_method = solver)
  
  # Compute scores (following package pattern)
  embed1 <- bases[[1]]$vectors
  embed2 <- bases[[2]]$vectors %*% alignment_result$rotation %*% assignment_result$mapping_matrix
  scores <- rbind(embed1, embed2)
  
  # Primal vectors computation (following KEMA pattern)
  v <- do.call(rbind, lapply(seq_along(strata), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- bases[[i]]$vectors
    Matrix::crossprod(xi, alpha_i)
  }))
  
  # Return multiblock_biprojector with proper class structure for S3 dispatch
  result <- multivarious::multiblock_biprojector(
    v = v,
    s = scores,
    sdev = apply(scores, 2, sd),
    preproc = proc,
    block_indices = block_indices,
    assignment = assignment_result$assignment,
    rotation = alignment_result$rotation,
    mapping_matrix = assignment_result$mapping_matrix,
    classes = "grasp"
  )
  
  # Add grasp class to outer object for S3 dispatch
  structure(result, class = c("grasp", class(result)))
}

#' Block A: Spectral Basis Construction - Corrected
#' 
#' @param strata List of data domains 
#' @param ncomp Number of spectral components
#' @param use_laplacian Whether to use Laplacian normalization
#' @return List of spectral bases for each domain
#' @keywords internal
#' @importFrom RSpectra eigs_sym
compute_grasp_basis <- function(strata, ncomp, use_laplacian = TRUE) {
  # Validate input
  if (!is.list(strata) || length(strata) == 0) {
    stop("strata must be a non-empty list", call. = FALSE)
  }
  
  bases <- lapply(strata, function(stratum) {
    # Ensure stratum$x is a matrix for neighborweights
    if (!is.matrix(stratum$x)) {
      stratum$x <- as.matrix(stratum$x)
      if (!is.matrix(stratum$x)) {
        stop("Failed to convert domain data to matrix format", call. = FALSE)
      }
    }
    
    # Validate dimensions
    if (nrow(stratum$x) < 3) {
      stop("GRASP requires at least 3 samples per domain", call. = FALSE)
    }
    
    # Construct graph using package patterns (following KEMA pattern)
    graph_weights <- safe_compute(
      neighborweights::graph_weights(
        stratum$x,
        weight_mode = "normalized",
        neighbor_mode = "knn", 
        k = min(max(3, nrow(stratum$x) %/% 3), nrow(stratum$x) - 1),  # Adaptive k
        type = "normal",
        sigma = 0.5  # More conservative sigma for better locality
      ),
      "Graph construction failed"
    )
    
    # Extract adjacency matrix
    adj_matrix <- neighborweights::adjacency(graph_weights)
    
    # Apply Laplacian normalization if requested
    if (use_laplacian) {
      # Improved isolated node handling (your suggestion)
      degrees <- Matrix::rowSums(adj_matrix)
      isolated_nodes <- degrees == 0
      
      if (any(isolated_nodes)) {
        warning("Found ", sum(isolated_nodes), " isolated nodes. They will be handled with zero eigenvalues.")
        # Keep isolated nodes with zero row/col rather than perturbing
        degrees[isolated_nodes] <- 1  # Temporary for sqrt computation
      }
      
      D_inv_sqrt <- Matrix::Diagonal(n = length(degrees), x = 1 / sqrt(degrees))
      if (any(isolated_nodes)) {
        isolated_indices <- which(isolated_nodes)
        D_inv_sqrt[isolated_indices, isolated_indices] <- 0  # Zero out isolated nodes
      }
      
      L <- Matrix::Diagonal(nrow(adj_matrix)) - D_inv_sqrt %*% adj_matrix %*% D_inv_sqrt
    } else {
      L <- adj_matrix
    }
    
    # Eigenvalue computation - improved selection (your suggestion)
    k_to_request <- min(ncomp + 1, nrow(L) - 1)
    
    # Use RSpectra by default with PRIMME as fallback (your suggestion)
    decomp <- safe_compute({
      if (requireNamespace("RSpectra", quietly = TRUE)) {
        RSpectra::eigs_sym(L, k = k_to_request, which = "SA")
      } else {
        PRIMME::eigs_sym(L, NEig = k_to_request, which = "SA")
      }
    }, "Eigenvalue computation failed for GRASP basis")
    
    # Improved eigenvector selection (your suggestion)
    if (use_laplacian) {
      # For Laplacian, take indices 2:(ncomp+1) to skip the constant vector
      if (k_to_request < ncomp + 1) {
        warning("Requested ", ncomp, " components but only ", k_to_request - 1, " available")
        ncomp_actual <- k_to_request - 1
      } else {
        ncomp_actual <- ncomp
      }
      
      vectors <- as(decomp$vectors[, 2:(ncomp_actual + 1), drop = FALSE], "CsparseMatrix")
      values <- decomp$values[2:(ncomp_actual + 1)]
    } else {
      # For adjacency matrix, use threshold as fallback
      non_trivial_idx <- which(abs(decomp$values) > 1e-10)
      ncomp_actual <- min(ncomp, length(non_trivial_idx))
      
      vectors <- as(decomp$vectors[, non_trivial_idx[1:ncomp_actual], drop = FALSE], "CsparseMatrix")
      values <- decomp$values[non_trivial_idx[1:ncomp_actual]]
    }
    
    list(vectors = vectors, values = values)
  })
  
  bases
}

#' Block B: Multi-scale Descriptors - Corrected
#' 
#' @param bases List of spectral bases
#' @param q_descriptors Number of descriptor functions
#' @param sigma Diffusion parameter
#' @return List of descriptor matrices
#' @keywords internal
compute_grasp_descriptors <- function(bases, q_descriptors, sigma = 0.73) {
  # Validation
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  
  # Generate time steps
  time_steps <- seq(0.1, 50, length.out = q_descriptors) * sigma
  
  descriptors <- lapply(bases, function(basis) {
    phi <- basis$vectors
    lambda_vals <- basis$values
    
    # Improved sparse matrix handling (your suggestion)
    # Keep phi^2 sparse
    phi_sq <- phi
    phi_sq@x <- phi_sq@x^2
    
    # Vectorized computation with numerical stability (your suggestion)
    # Compute in log-space to prevent underflow
    logH <- -outer(lambda_vals, time_steps)  # k x q
    H <- exp(pmax(logH, -745))  # double-precision floor
    desc_matrix <- phi_sq %*% H  # n x q (sparse x dense)
    
    # Efficient column normalization (your suggestion)
    col_norms <- sqrt(Matrix::colSums(desc_matrix^2))
    desc_matrix <- desc_matrix %*% Matrix::Diagonal(n = length(col_norms), x = 1 / pmax(col_norms, 1e-12))
    
    as(desc_matrix, "CsparseMatrix")
  })
  
  descriptors
}

#' Block C: Base Alignment - Corrected
#' 
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis  
#' @param desc1 First domain's descriptors
#' @param desc2 Second domain's descriptors
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return List with rotation matrix and convergence info
#' @keywords internal
align_grasp_bases <- function(basis1, basis2, desc1, desc2, lambda = 0.1, 
                             max_iter = 200, tol = 1e-8) {
  # Validation
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_number(max_iter)
  chk::chk_number(tol)
  
  # Extract matrices
  phi1 <- basis1$vectors
  phi2 <- basis2$vectors
  lambda2_vals <- basis2$values
  
  k <- ncol(phi1)
  M <- Matrix::Diagonal(k)  # Initialize as identity
  
  # Pre-compute constant parts
  F_T_Phi1 <- Matrix::crossprod(desc1, phi1)
  G_T_Phi2 <- Matrix::crossprod(desc2, phi2)
  Lambda2 <- Matrix::Diagonal(n = length(lambda2_vals), x = lambda2_vals)
  
  # Corrected objective function (your suggestion - remove λ₁ subtraction)
  objective_fn <- function(M) {
    aligned_lambda2_matrix <- Matrix::crossprod(M, Lambda2 %*% M)
    # off() is sum of squares of off-diagonal elements (independent of λ₁)
    diag_vals <- Matrix::diag(aligned_lambda2_matrix)
    off_mat <- aligned_lambda2_matrix - Matrix::Diagonal(n = length(diag_vals), x = diag_vals)
    off_penalty <- Matrix::norm(off_mat, "F")^2
    
    desc_penalty <- Matrix::norm(F_T_Phi1 - G_T_Phi2 %*% M, "F")^2
    return(off_penalty + lambda * desc_penalty)
  }

  # Corrected gradient function
  gradient_fn <- function(M) {
    A <- Lambda2 %*% M
    AtMA <- Matrix::crossprod(M, A)
    
    # Corrected gradient of off-diagonal penalty
    diag_AtMA <- Matrix::diag(AtMA)
    grad_off <- 2 * (A %*% AtMA - A %*% Matrix::Diagonal(n = length(diag_AtMA), x = diag_AtMA))
    
    # Gradient of descriptor penalty (optimized for sparse operations)
    # G_T_Phi2 is q×k, M is k×k, result should be k×k
    residual <- F_T_Phi1 - G_T_Phi2 %*% M  # q×k
    grad_desc <- -2 * Matrix::crossprod(G_T_Phi2, residual)  # sparse-friendly
    
    euc_grad <- grad_off + lambda * grad_desc
    
    # Project to tangent space of Stiefel manifold
    riemannian_grad <- euc_grad - M %*% Matrix::crossprod(M, euc_grad)
    return(riemannian_grad)
  }

  # Improved optimization with relative tolerance (your suggestion)
  prev_obj <- objective_fn(M)
  
  for (iter in seq_len(max_iter)) {
    grad <- gradient_fn(M)
    
    # Backtracking line search with smaller initial step
    step_size <- 0.5
    for (ls_iter in 1:8) {
      M_new_euc <- M - step_size * grad
      
      # SVD retraction (your suggestion - more stable than QR for small k)
      M_new_mat <- as.matrix(M_new_euc)
      if (any(!is.finite(M_new_mat))) {
        stop("Non-finite values in matrix during SVD retraction", call. = FALSE)
      }
      svd_res <- svd(M_new_mat)
      M_new <- as(svd_res$u %*% t(svd_res$v), "CsparseMatrix")
      
      new_obj <- objective_fn(M_new)
      if (new_obj < prev_obj) break
      step_size <- step_size * 0.5
    }
    
    M <- M_new
    
    # Improved convergence check - require minimum iterations
    rel_change <- abs(prev_obj - new_obj) / max(1, prev_obj)
    if (iter > 5 && rel_change < tol) {  # Minimum 5 iterations
      break
    }
    prev_obj <- new_obj
  }
  
  list(rotation = M, iterations = iter, converged = iter < max_iter)
}

#' Blocks D & E: Functional Map and Assignment - Corrected
#' 
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis
#' @param desc1 First domain's descriptors  
#' @param desc2 Second domain's descriptors
#' @param M Rotation matrix from base alignment
#' @param distance_method Distance method: "cosine" or "euclidean"
#' @param solver_method Assignment solver: "linear", "hungarian", or "auction"
#' @return List with assignment and related matrices
#' @keywords internal
compute_grasp_assignment <- function(basis1, basis2, desc1, desc2, M, 
                                   distance_method = "cosine", 
                                   solver_method = "linear") {
  # Validation
  if (!distance_method %in% c("cosine", "euclidean")) {
    stop("distance_method must be 'cosine' or 'euclidean'", call. = FALSE)
  }
  
  if (!solver_method %in% c("linear", "hungarian", "auction")) {
    stop("solver_method must be 'linear', 'hungarian', or 'auction'", call. = FALSE)
  }
  
  # Block D: Functional correspondence with ridge regularization (your suggestion)
  phi1 <- basis1$vectors
  phi2_aligned <- basis2$vectors %*% M
  
  # Solve F^T*Phi1 ≈ G^T*Phi2_aligned * C for diagonal C
  A_prime <- Matrix::crossprod(desc2, phi2_aligned)  # (q, k)
  B_prime <- Matrix::crossprod(desc1, phi1)         # (q, k)
  
  # Ridge regularization for numerical stability (your suggestion)
  c_diag <- Matrix::colSums(A_prime * B_prime) / (Matrix::colSums(A_prime^2) + 1e-9)
  C <- Matrix::Diagonal(n = length(c_diag), x = c_diag)
  
  # Block E: Final embeddings and assignment
  embed1 <- as.matrix(phi1)
  embed2 <- as.matrix(phi2_aligned %*% C)
  
  # Improved distance computation (cosine similarity)
  if (distance_method == "cosine") {
    # Robust cosine computation without external dependencies
    embed1_norms <- sqrt(rowSums(embed1^2) + 1e-12)
    embed2_norms <- sqrt(rowSums(embed2^2) + 1e-12)
    embed1_norm <- embed1 / embed1_norms
    embed2_norm <- embed2 / embed2_norms
    
    sim_matrix <- embed1_norm %*% t(embed2_norm)
    cost_matrix <- 1 - sim_matrix
    cost_matrix[cost_matrix < 0] <- 0  # Clamp negative values to preserve matrix structure
  } else {
    # Efficient Euclidean distance computation
    embed1_sq_norms <- rowSums(embed1^2)
    embed2_sq_norms <- rowSums(embed2^2)
    cross_product <- embed1 %*% t(embed2)
    
    cost_matrix <- sqrt(outer(embed1_sq_norms, embed2_sq_norms, "+") - 2 * cross_product)
    cost_matrix[cost_matrix < 0] <- 0  # Clamp negative values to preserve matrix structure
  }
  
  # Assignment solver with auction option (your suggestion)
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }
  
  # Assignment solver selection
  if (nrow(cost_matrix) >= 5000 && solver_method == "auction") {
    assignment <- clue::solve_LSAP(cost_matrix, method = "auction")
  } else if (solver_method == "hungarian") {
    assignment <- clue::solve_LSAP(cost_matrix, method = "dense")  # default Hungarian
  } else {
    assignment <- clue::solve_LSAP(cost_matrix)  # alias for "dense"
  }
  
  list(assignment = as.integer(assignment), 
       cost_matrix = cost_matrix,
       mapping_matrix = C)
}

#' @export
grasp <- function(data, ...) {
  UseMethod("grasp")
}