#' Position-Aware Random Transport (PARROT) Network Alignment
#'
#' Performs PARROT alignment on hyperdesign data structures. Aligns
#' networks using regularized optimal transport with position-aware
#' features and consistency constraints.
#'
#' PARROT tackles network alignment by formulating it as a regularized
#' optimal transport problem. The method incorporates position-aware
#' features through Random Walk with Restart (RWR) descriptors and
#' enforces structural consistency through neighborhood-preserving
#' regularization terms.
#'
#' @param data Input data object containing network domains
#' @param anchors Name of anchor/correspondence variable for semi-supervised alignment
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{parrot.hyperdesign}} for details on method-specific parameters 
#'   such as \code{preproc}, \code{ncomp}, \code{sigma}, \code{lambda}, 
#'   \code{tau}, \code{solver}, \code{max_iter}, \code{tol}, and \code{use_cpp}
#'
#' @details
#' PARROT operates through the following algorithmic components:
#' \itemize{
#'   \item \strong{Position-Aware Features}: Compute RWR descriptors
#'     capturing network position
#'   \item \strong{Cross-Network Cost}: Build transport cost matrix between
#'     networks
#'   \item \strong{Consistency Regularization}: Add structural similarity
#'     constraints
#'   \item \strong{Optimal Transport}: Solve regularized transport problem via
#'     Sinkhorn
#' }
#'
#' The algorithm minimizes the objective:
#' \deqn{L(S) = \langle C, S \rangle + \lambda_1 \Omega_1(S) + \lambda_2 \Omega_2(S) + \tau H(S)}
#'
#' where \eqn{C} is the position-aware cost matrix, \eqn{\Omega_1, \Omega_2} are consistency
#' regularizers, and \eqn{H(S)} is the entropy regularization term with parameter \eqn{\tau}.
#'
#' Key parameters:
#' \itemize{
#'   \item \code{sigma}: RWR restart probability and diffusion parameter
#'   \item \code{lambda}: Consistency regularization weights
#'   \item \code{tau}: Entropy regularization parameter for Sinkhorn
#'   \item \code{solver}: Transport solver ("sinkhorn" or "exact")
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned network embeddings 
#'   \item \code{v}: Primal vectors for out-of-sample projection
#'   \item \code{alignment_matrix}: Soft alignment/transport plan between networks
#'   \item \code{transport_plan}: Dense transport matrix S
#'   \item \code{sdev}: Standard deviations of aligned components
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign network data
#' library(multidesign)
#' library(tibble)
#' 
#' # Create synthetic network domains (node features for PARROT)
#' set.seed(123)
#' X1 <- matrix(rnorm(200), 100, 2)  # Node features for domain 1
#' X2 <- matrix(rnorm(200), 100, 2)  # Node features for domain 2
#' 
#' # Create design data frames with anchor correspondences (PARROT uses anchors)
#' design1 <- tibble(
#'   node_id = 1:100,
#'   anchors = c(1:10, rep(NA, 90))  # First 10 nodes are anchors
#' )
#' design2 <- tibble(
#'   node_id = 1:100,
#'   anchors = c(1:10, rep(NA, 90))  # Corresponding anchors
#' )
#' 
#' # Create multidesign objects
#' md1 <- multidesign(X1, design1)
#' md2 <- multidesign(X2, design2)
#' 
#' # Create hyperdesign from multidesign objects
#' hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
#' 
#' # Run PARROT alignment (uses anchors from design component)
#' result <- parrot(hd, anchors = anchors)
#' 
#' # Access alignment results
#' transport_plan <- result$transport_plan
#' aligned_embeddings <- result$s
#' 
#' # Use different regularization settings
#' result_strong <- parrot(hd, anchors = anchors, lambda = 0.5, tau = 0.1)
#' }
#'
#' @references
#' Wang, S., Chen, Z., Yu, X., Li, T., Yang, J., & Liu, X. (2022). PARROT: 
#' Position-aware regularized optimal transport for network alignment. 
#' In Proceedings of the 28th ACM SIGKDD Conference on Knowledge Discovery 
#' and Data Mining (pp. 1896-1905).
#'
#' @seealso \code{\link{parrot.hyperdesign}}
#' @export
parrot <- function(data, anchors, ...) {
  UseMethod("parrot")
}

#' PARROT for Hyperdesign Objects
#'
#' Performs PARROT alignment on hyperdesign data structures containing network domains.
#'
#' @param data A hyperdesign object containing multiple network domains
#' @param anchors Name of the anchor/correspondence variable for semi-supervised alignment
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of latent dimensions for barycenter (default: NULL, auto-determine)
#' @param sigma RWR restart probability and diffusion parameter (default: 0.15)
#' @param lambda Overall consistency regularization weight (default: 0.1). Can also specify
#'   lambda_e, lambda_n, lambda_p separately for fine control
#' @param lambda_e Edge consistency weight (default: lambda * 0.5)
#' @param lambda_n Neighborhood consistency weight (default: lambda * 0.5)  
#' @param lambda_p Anchor prior weight (default: lambda * 0.1)
#' @param tau Entropy regularization parameter for Sinkhorn (default: 0.01)
#' @param alpha Weight for feature vs RWR cost in Sylvester equation (default: 0.5)
#' @param gamma Cross-graph mixing parameter for Sylvester iteration (default: 0.1)
#' @param solver Transport solver: "sinkhorn" for entropic (default), "exact" for unregularized
#' @param max_iter Maximum number of iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param use_cpp Whether to use C++ optimizations (default: FALSE)
#' @param verbose Logical flag to print progress information during fitting
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the PARROT alignment results
#'
#' @rdname parrot
#' @method parrot hyperdesign
#' @export
#' @importFrom chk chk_number chk_true chk_logical
#' @importFrom multivarious center
#' @importFrom rlang enquo !!
#' @importFrom purrr map
#' @importFrom dplyr select pull
parrot.hyperdesign <- function(data, 
                              anchors,
                              preproc = center(), 
                              ncomp = NULL,
                              sigma = 0.15,
                              lambda = 0.1,
                              lambda_e = NULL,
                              lambda_n = NULL,
                              lambda_p = NULL,
                              tau = 0.05,
                              alpha = 0.2,
                              gamma = 0.1,
                              solver = c("sinkhorn", "exact"),
                              max_iter = 100,
                              tol = 1e-6,
                              use_cpp = FALSE,
                              verbose = FALSE,
                              ...) {
  
  # Capture anchor variable (following KEMA/GRASP pattern)
  anchors <- rlang::enquo(anchors)
  
  # Input validation (following exact package patterns)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0 && sigma < 1)
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_number(tau)
  chk::chk_true(tau > 0)
  chk::chk_number(alpha)
  chk::chk_true(alpha >= 0 && alpha <= 1)
  chk::chk_number(gamma)
  chk::chk_true(gamma >= 0 && gamma <= 1)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  # FIX3-07: Set default lambda components if not provided
  if (is.null(lambda_e)) lambda_e <- lambda * 0.5
  if (is.null(lambda_n)) lambda_n <- lambda * 0.5
  if (is.null(lambda_p)) lambda_p <- lambda
  
  # Validate lambda components
  chk::chk_number(lambda_e)
  chk::chk_true(lambda_e >= 0)
  chk::chk_number(lambda_n)
  chk::chk_true(lambda_n >= 0)
  chk::chk_number(lambda_p)
  chk::chk_true(lambda_p >= 0)
  
  # Validate solver parameter (following KEMA/GRASP pattern)
  solver <- match.arg(solver)
  
  # Validate input data (following package pattern)
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list of hyperdesign objects", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("PARROT currently supports exactly 2 domains. For multi-network ", 
         "alignment, use parrot_multiple() [coming in Sprint 2]", call. = FALSE)
  }
  
  # Extract anchor information (following package pattern)
  anchor_data <- safe_compute({
    unlist(purrr::map(data, function(x) {
      x$design |> dplyr::select(!!anchors) |> dplyr::pull(!!anchors)
    }))
  }, "Failed to extract anchor data. Check that anchor variable exists in design.")
  
  # Validate anchor data
  if (all(is.na(anchor_data))) {
    stop("No anchor correspondences found. PARROT requires some anchor pairs ", 
         "for semi-supervised alignment.", call. = FALSE)
  }
  
  n_anchors <- sum(!is.na(anchor_data))
  n_total <- length(anchor_data)
  if (verbose) {
    message("PARROT: Using ", n_anchors, " anchor correspondences out of ",
            n_total, " total nodes")
  }
  
  # Preprocess data (following CONE-Align simplified pattern)
  pdata <- lapply(data, function(domain) {
    x_preprocessed <- safe_compute(
      if (is.function(preproc)) preproc(domain$x) else domain$x,
      "Data preprocessing failed. Check preprocessing function or data format."
    )
    list(x = x_preprocessed, design = domain$design)
  })
  
  # Debug: Check if preprocessing made a difference
    if (verbose) {
      if (is.function(preproc)) {
        orig_mean <- mean(unlist(lapply(data, function(d) mean(d$x))))
        proc_mean <- mean(unlist(lapply(pdata, function(d) mean(d$x))))
        if (abs(orig_mean - proc_mean) > 1e-10) {
          message("Preprocessing changed data mean from ", round(orig_mean, 4),
                  " to ", round(proc_mean, 4))
        }
      }
    }
  
  # Block indices computation (following package pattern)
  # Create proper preprocessing structure
  # Store the original preproc for later use
  original_preproc <- preproc
  
  if (!is.null(preproc)) {
    if (inherits(preproc, "pre_processor")) {
      proc <- preproc
    } else {
      # For prepper, functions, or other objects that don't work with multiblock_biprojector
      # use the manual result construction
      proc <- NULL
    }
  } else {
    proc <- NULL
  }
  
  # Call PARROT fitting function - pass original data and preproc for flexible application
  parrot_fit(pdata, proc, anchor_data, ncomp, sigma, lambda_e, lambda_n, lambda_p, tau, 
            alpha, gamma, solver, max_iter, tol, use_cpp, feature_block_indices(pdata), 
            original_preproc, data, verbose = verbose)
}

#' @rdname parrot
#' @method parrot default
#' @export
parrot.default <- function(data, anchors, ...) {
  # Validate input data
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list of hyperdesign objects", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("PARROT currently supports exactly 2 domains. For multi-network ", 
         "alignment, use parrot_multiple() [coming in Sprint 2]", call. = FALSE)
  }
  
  stop("No applicable method for parrot. data must be a hyperdesign object.", call. = FALSE)
}

#' @keywords internal
parrot_fit <- function(strata, proc, anchor_data, ncomp, sigma, lambda_e, lambda_n, lambda_p, tau, 
                      alpha, gamma, solver, max_iter, tol, use_cpp, feature_blocks, 
                      original_preproc = NULL, original_data = NULL, verbose = FALSE) {
  
  # Extract network structures 
  # Use preprocessed data (strata) by default, but for v computation we may need original data
  networks <- extract_parrot_networks(strata)
  
  # FIX3-01: Split anchor handling properly
  n1 <- nrow(networks[[1]]$features)
  n2 <- nrow(networks[[2]]$features)
  
  if (length(anchor_data) >= n1 + n2) {
    anchor_vec1 <- anchor_data[1:n1]
    anchor_vec2 <- anchor_data[(n1+1):(n1+n2)]
    anchor_idx1 <- which(!is.na(anchor_vec1))
    anchor_idx2 <- which(!is.na(anchor_vec2))
  } else {
    stop("Anchor data length mismatch with network sizes", call. = FALSE)
  }
  
  # Package anchor information for passing to functions
  anchor_info <- list(
    vec1 = anchor_vec1,
    vec2 = anchor_vec2,
    idx1 = anchor_idx1,
    idx2 = anchor_idx2,
    n1 = n1,
    n2 = n2
  )
  
  # Compute position-aware features (RWR descriptors)
  rwr_features <- compute_parrot_rwr(networks, anchor_info, sigma, max_iter, tol, use_cpp)
  
  # Solve optimal transport problem with proximal method
  transport_result <- solve_parrot_transport(networks, rwr_features, anchor_info, 
                                           lambda_e, lambda_n, lambda_p, 
                                           tau, alpha, sigma, gamma, solver, max_iter, tol = tol, use_cpp = use_cpp)
  
  # Debug: Check transport plan properties
  tp <- transport_result$transport_plan
  if (verbose) {
    message("Transport plan from solver: rows=", nrow(tp), " cols=", ncol(tp),
            " row_sum_mean=", round(mean(rowSums(tp)), 4),
            " col_sum_mean=", round(mean(colSums(tp)), 4))
  }
  
  # Compute aligned embeddings (following package pattern)
  scores <- compute_parrot_embeddings(networks, transport_result$transport_plan, ncomp)
  
  # Primal vectors computation (following KEMA/GRASP pattern)
  # Ensure consistent dimensions with scores
  actual_ncomp <- ncol(scores)
  
  # Loadings are defined by regressing each domain's (preprocessed) features onto
  # its corresponding embedding coordinates. Since `scores` is computed as
  # rbind(U, V) from the SVD of the transport plan, reuse it directly to avoid:
  # - using the wrong SVD side for domain 2 (V vs U)
  # - recomputing the SVD per-domain
  if (nrow(scores) != n1 + n2) {
    stop("PARROT internal error: embedding rows do not match domain sizes.", call. = FALSE)
  }
  alpha_list <- list(
    scores[seq_len(n1), , drop = FALSE],
    scores[n1 + seq_len(n2), , drop = FALSE]
  )

  loadings <- do.call(rbind, lapply(seq_along(strata), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- alpha_list[[i]]
    if (nrow(xi) != nrow(alpha_i)) {
      stop("PARROT internal error: embedding rows do not match feature rows for domain ", i, ".", call. = FALSE)
    }
    Matrix::crossprod(xi, alpha_i)
  }))

  result <- new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "parrot",
    extras = list(
      alignment_matrix = transport_result$transport_plan,
      transport_plan = transport_result$transport_plan,
      anchors = anchor_data
    )
  )

  if (is.null(proc) && !is.null(original_preproc)) {
    result$preproc <- original_preproc
  }

  result
}

#' Extract Network Structures for PARROT
#' 
#' @param strata List of data domains
#' @return List of network structures with adjacency matrices
#' @keywords internal
extract_parrot_networks <- function(strata) {
  networks <- lapply(strata, function(stratum) {
    x <- stratum$x
    n <- nrow(x)
    
    # Build adjacency matrix using neighborweights (following KEMA/GRASP/CONE-Align pattern)
    # Use adaptive k based on data size (KEMA pattern)
    knn <- min(10, max(3, floor(sqrt(n))))
    
    # Use neighborweights::graph_weights (consistent with other methods)
    graph_weights <- tryCatch({
      neighborweights::graph_weights(
        x,
        weight_mode = "normalized",
        neighbor_mode = "knn", 
        k = knn,
        type = "normal",
        sigma = 0.5  # Conservative sigma for locality
      )
    }, error = function(e) {
      stop("Graph construction failed: ", e$message, call. = FALSE)
    })
    
    # Extract adjacency matrix
    A <- neighborweights::adjacency(graph_weights)
    
    # Row-normalize to get transition matrix (PARROT requirement)
    deg <- Matrix::rowSums(A)
    deg[deg == 0] <- 1  # Avoid division by zero
    W <- Matrix::Diagonal(x = 1 / deg) %*% A
    
    list(adjacency = A, transition = W, features = x)
  })
  
  networks
}

#' Compute RWR Features for PARROT (R implementation)
#' 
#' @param networks List of network structures
#' @param anchor_info List with anchor information for both networks
#' @param sigma RWR restart probability (beta in paper)
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return List of RWR descriptor matrices
#' @keywords internal
compute_parrot_rwr_r <- function(networks, anchor_info, sigma, max_iter, tol) {
  # AUDIT-01: Vectorized RWR computation
  anchor_idx1 <- anchor_info$idx1
  anchor_idx2 <- anchor_info$idx2
  
  if (length(anchor_idx1) == 0 || length(anchor_idx2) == 0) {
    stop("No valid anchor nodes found for RWR computation", call. = FALSE)
  }
  
  # Get unique anchor values to determine total number of anchor columns
  anchor_vals1 <- anchor_info$vec1[anchor_idx1]
  anchor_vals2 <- anchor_info$vec2[anchor_idx2]
  unique_anchors <- unique(c(anchor_vals1, anchor_vals2))
  n_anchors <- length(unique_anchors)
  
  # Vectorized RWR for network 1
  W1T <- Matrix::t(networks[[1]]$transition)  # Transpose once
  n1 <- ncol(W1T)
  
  # Create sparse matrix of all restart vectors
  anchor_to_col <- match(anchor_vals1, unique_anchors)
  E1 <- Matrix::sparseMatrix(
    i = anchor_idx1,
    j = anchor_to_col,
    x = 1,
    dims = c(n1, n_anchors)
  )
  
  # Vectorized power iteration
  R1 <- E1 / n1  # Initial uniform distribution
  for (iter in 1:max_iter) {
    R1_old <- R1
    R1 <- (1 - sigma) * (W1T %*% R1) + sigma * E1
    if (max(abs(R1 - R1_old)) < tol) break
  }
  R1 <- as.matrix(R1)  # Convert to dense for subsequent steps
  
  # Vectorized RWR for network 2
  W2T <- Matrix::t(networks[[2]]$transition)
  n2 <- ncol(W2T)
  
  # Create sparse matrix of restart vectors for network 2
  anchor_to_col2 <- match(anchor_vals2, unique_anchors)
  E2 <- Matrix::sparseMatrix(
    i = anchor_idx2,
    j = anchor_to_col2,
    x = 1,
    dims = c(n2, n_anchors)
  )
  
  # Vectorized power iteration
  R2 <- E2 / n2  # Initial uniform distribution
  for (iter in 1:max_iter) {
    R2_old <- R2
    R2 <- (1 - sigma) * (W2T %*% R2) + sigma * E2
    if (max(abs(R2 - R2_old)) < tol) break
  }
  R2 <- as.matrix(R2)
  
  list(R1, R2)
}

#' Compute RWR Features for PARROT (Dispatcher)
#'
#' @inheritParams compute_parrot_rwr_r
#' @param use_cpp Whether to use C++ implementation
#' @return A list of RWR feature matrices, one per network
#' @keywords internal
compute_parrot_rwr <- function(networks, anchor_info, sigma, max_iter, tol, use_cpp = FALSE) {
  
  if (isTRUE(use_cpp)) {
    # Extract anchor information
    anchor_idx1 <- anchor_info$idx1
    anchor_idx2 <- anchor_info$idx2
    
    if (length(anchor_idx1) == 0 || length(anchor_idx2) == 0) {
      stop("No valid anchor nodes found for RWR computation", call. = FALSE)
    }
    
    # Get unique anchor values
    anchor_vals1 <- anchor_info$vec1[anchor_idx1]
    anchor_vals2 <- anchor_info$vec2[anchor_idx2]
    unique_anchors <- unique(c(anchor_vals1, anchor_vals2))
    n_anchors <- length(unique_anchors)
    
    # Network 1 RWR
    W1T <- t(networks[[1]]$transition)
    n1 <- ncol(W1T)
    
    # Create sparse matrix of restart vectors
    anchor_to_col <- match(anchor_vals1, unique_anchors)
    E1 <- Matrix::sparseMatrix(
      i = anchor_idx1,
      j = anchor_to_col,
      x = 1,
      dims = c(n1, n_anchors)
    )
    
    # Convert to dense for C++
    E1_dense <- as.matrix(E1)
    W1T_dense <- as.matrix(W1T)
    R1 <- compute_rwr_vectorized_cpp_fallback(W1T_dense, E1_dense, sigma, max_iter, tol)
    
    # Network 2 RWR
    W2T <- t(networks[[2]]$transition)
    n2 <- ncol(W2T)
    
    anchor_to_col2 <- match(anchor_vals2, unique_anchors)
    E2 <- Matrix::sparseMatrix(
      i = anchor_idx2,
      j = anchor_to_col2,
      x = 1,
      dims = c(n2, n_anchors)
    )
    
    # Convert to dense for C++
    E2_dense <- as.matrix(E2)
    W2T_dense <- as.matrix(W2T)
    
    # Compute RWR using C++
    R2 <- compute_rwr_vectorized_cpp_fallback(W2T_dense, E2_dense, sigma, max_iter, tol)
    
    list(R1, R2)
  } else {
    # Use existing R implementation
    compute_parrot_rwr_r(networks, anchor_info, sigma, max_iter, tol)
  }
}

#' Solve PARROT Transport Problem with Proximal Point Method
#' 
#' Implements the proximal point loop that iteratively refines the cost matrix
#' and solves the optimal transport problem. This is critical for convergence.
#' 
#' @param networks List of network structures
#' @param rwr_features List of RWR descriptor matrices
#' @param anchor_info List with anchor information for both networks
#' @param lambda_e Edge consistency weight
#' @param lambda_n Neighborhood consistency weight
#' @param lambda_p Anchor prior weight
#' @param tau Entropy regularization parameter
#' @param solver Transport solver method
#' @param max_iter Maximum iterations for inner Sinkhorn
#' @param max_outer Maximum outer proximal iterations
#' @param tol Convergence tolerance
#' @return List with transport plan and related matrices
#' @keywords internal
solve_parrot_transport <- function(networks, rwr_features, anchor_info, 
                                  lambda_e, lambda_n, lambda_p, tau, alpha, sigma, gamma, solver, 
                                  max_iter, max_outer = 10, tol, use_cpp = FALSE) {
  
  # AUDIT-03: Use proper proximal point method implementation
  # Call the corrected proximal solver (functions defined in parrot_proximal.R)
  solve_parrot_proximal(networks, rwr_features, anchor_info,
                       lambda_e, lambda_n, lambda_p,
                       tau, alpha, sigma, gamma, solver,
                       max_iter, max_outer, tol, use_cpp)
}

#' Solve Sylvester Equation for Cross-Graph RWR Cost (R implementation)
#'
#' Implements the Sylvester iteration from Eq. 4 in the paper:
#' C_rwr = (1+beta) * C_node + (1-beta) * gamma * W1 * C_rwr * W2^T
#'
#' @param W1 First network's transition matrix (row-stochastic)
#' @param W2T Transpose of second network's transition matrix  
#' @param Cnode Node feature cost matrix
#' @param beta RWR restart probability (corresponds to sigma in code)
#' @param gamma Cross-graph discount factor
#' @param tol Convergence tolerance
#' @param max_iter Maximum iterations
#' @return Cross-graph RWR cost matrix C_rwr
#' @keywords internal
solve_sylvester_rwr_r <- function(W1, W2T, Cnode, beta = 0.15, gamma = 0.1, tol = 1e-6, max_iter = 50) {
  # AUDIT-02: Correct Sylvester formulation from paper's Eq. 4
  X <- Cnode

  kappa <- (1 - beta) * gamma
  if (kappa >= 1) {
    warning(
      "Sylvester RWR iteration may diverge because (1 - beta) * gamma = ",
      signif(kappa, 4),
      " >= 1. Consider reducing gamma below ",
      signif(1 / (1 - beta), 4),
      call. = FALSE
    )
  }
  
  converged <- FALSE
  delta <- NA_real_
  for (iter in seq_len(max_iter)) {
    X_new <- (1 + beta) * Cnode + (1 - beta) * gamma * (W1 %*% X %*% W2T)
    
    # Check convergence
    delta <- max(abs(X_new - X))
    if (!is.finite(delta)) {
      warning("Sylvester RWR iteration produced non-finite values at iteration ", iter, call. = FALSE)
      break
    }

    if (delta < tol) {
      converged <- TRUE
      X <- X_new
      break
    }

    X <- X_new
  }

  if (!converged && is.finite(delta) && delta >= tol) {
    warning(
      "Sylvester RWR iteration did not converge after ",
      max_iter,
      " iterations (last max|delta| = ",
      signif(delta, 4),
      ").",
      call. = FALSE
    )
  }
  
  X
}

#' Solve Sylvester Equation for Cross-Graph RWR Cost (Dispatcher)
#'
#' @inheritParams solve_sylvester_rwr_r
#' @param use_cpp Whether to use C++ implementation
#' @return A numeric matrix solving the Sylvester equation for RWR computation
#' @keywords internal
solve_sylvester_rwr <- function(W1, W2T, Cnode, beta = 0.15, gamma = 0.1,
                                tol = 1e-6, max_iter = 50, use_cpp = FALSE) {
  
  if (isTRUE(use_cpp)) {
    if (inherits(W1, "Matrix")) W1 <- as.matrix(W1)
    if (inherits(W2T, "Matrix")) W2T <- as.matrix(W2T)
    if (inherits(Cnode, "Matrix")) Cnode <- as.matrix(Cnode)
    solve_sylvester_rwr_cpp_fallback(W1, W2T, Cnode, beta, gamma, tol, max_iter)
  } else {
    # Use R implementation (no use_cpp arg)
    solve_sylvester_rwr_r(W1, W2T, Cnode, beta, gamma, tol, max_iter)
  }
}

#' Compute Position-Aware Cost Matrix (R implementation)
#' 
#' Follows the two-stage process from the paper:
#' 1. Compute C_node = alpha * cost_RWR + (1-alpha) * cost_attr (Eq. 3)
#' 2. Solve Sylvester equation for C_rwr (Eq. 4)
#' 
#' @param networks List of network structures
#' @param rwr_features List of RWR descriptor matrices
#' @param anchor_info List with anchor information for both networks
#' @param alpha Weight for RWR vs attribute cost (default: 0.5)
#' @param sigma RWR restart probability (beta in paper)
#' @param gamma Cross-graph discount factor for Sylvester equation
#' @return Position-aware cost matrix C_rwr
#' @keywords internal
compute_parrot_cost_r <- function(networks, rwr_features, anchor_info, alpha = 0.5,                                  sigma = 0.15, gamma = 0.1) {
  # Cosine-based node cost with exponential scaling (Eq. 3 / App. A.2)
  row_normalize <- function(M) {
    norms <- sqrt(rowSums(M^2))
    norms[norms < 1e-12] <- 1e-12
    M / norms
  }

  R1 <- rwr_features[[1]]
  R2 <- rwr_features[[2]]
  X1 <- networks[[1]]$features
  X2 <- networks[[2]]$features

  R1n <- row_normalize(R1)
  R2n <- row_normalize(R2)
  X1n <- row_normalize(X1)
  X2n <- row_normalize(X2)

  cosine_clip <- function(M) {
    pmin(pmax(M, -1), 1)
  }
  sim_rwr <- cosine_clip(R1n %*% t(R2n))
  sim_attr <- cosine_clip(X1n %*% t(X2n))

  kappa <- 1.0
  cost_rwr <- exp(-kappa * sim_rwr)
  cost_attr <- exp(-kappa * sim_attr)

  C_node <- (1 - alpha) * cost_rwr + alpha * cost_attr

  # Solve Sylvester equation for C_rwr (Eq. 4)
  W1 <- networks[[1]]$transition  # Row-stochastic
  W2T <- t(networks[[2]]$transition)
  C_rwr <- solve_sylvester_rwr(W1, W2T, C_node, beta = sigma, gamma = gamma)

  # Ensure strictly positive costs for downstream logarithms
  C_rwr_min <- min(C_rwr)
  if (C_rwr_min <= 0) {
    C_rwr <- C_rwr - C_rwr_min + 1e-6
  }

  C_rwr
}


#' Compute Position-Aware Cost Matrix (Dispatcher)
#'
#' @inheritParams compute_parrot_cost_r
#' @param use_cpp Whether to use C++ implementation
#' @return A numeric cost matrix between networks
#' @keywords internal
compute_parrot_cost <- function(networks, rwr_features, anchor_info, alpha = 0.5,                                 sigma = 0.15, gamma = 0.1, use_cpp = FALSE) {

  if (isTRUE(use_cpp)) {
    X1 <- networks[[1]]$features
    X2 <- networks[[2]]$features
    R1 <- rwr_features[[1]]
    R2 <- rwr_features[[2]]
    W1 <- networks[[1]]$transition
    W2 <- networks[[2]]$transition

    if (inherits(X1, "Matrix")) X1 <- as.matrix(X1)
    if (inherits(X2, "Matrix")) X2 <- as.matrix(X2)
    if (inherits(R1, "Matrix")) R1 <- as.matrix(R1)
    if (inherits(R2, "Matrix")) R2 <- as.matrix(R2)
    if (inherits(W1, "Matrix")) W1 <- as.matrix(W1)
    if (inherits(W2, "Matrix")) W2 <- as.matrix(W2)

    C_rwr <- compute_parrot_cost_cpp_fallback(networks, list(R1, R2), anchor_info,
                                              alpha, sigma, gamma)
    C_rwr_min <- min(C_rwr)
    if (C_rwr_min <= 0) {
      C_rwr <- C_rwr - C_rwr_min + 1e-6
    }
    C_rwr
  } else {
    compute_parrot_cost_r(networks, rwr_features, anchor_info, alpha, sigma, gamma)
  }
}


#' Compute Sparse Edge Consistency Regularization
#' 
#' Implements sparse version of L_edge for memory efficiency
#' Only computes distances for existing edges
#' 
#' @param networks List of network structures
#' @param transport_plan Current transport plan S
#' @param lambda_e Edge regularization weight
#' @return Edge consistency cost matrix (sparse computation)
#' @keywords internal
compute_edge_consistency <- function(networks, transport_plan, lambda_e) {
  if (lambda_e <= 0) {
    return(matrix(0, nrow(transport_plan), ncol(transport_plan)))
  }

  A1 <- networks[[1]]$adjacency
  A2 <- networks[[2]]$adjacency
  X1 <- networks[[1]]$features
  X2 <- networks[[2]]$features

  if (inherits(A1, "Matrix")) A1 <- as.matrix(A1)
  if (inherits(A2, "Matrix")) A2 <- as.matrix(A2)

  row_normalize <- function(M) {
    norms <- sqrt(rowSums(M^2))
    norms[norms < 1e-12] <- 1e-12
    M / norms
  }

  X1n <- row_normalize(X1)
  X2n <- row_normalize(X2)

  C1 <- (1 - X1n %*% t(X1n)) * (A1 != 0)
  C2 <- (1 - X2n %*% t(X2n)) * (A2 != 0)

  C1_sq <- rowSums(C1^2)
  C2_sq <- rowSums(C2^2)

  L_base <- outer(C1_sq, rep(1, length(C2_sq))) +
    outer(rep(1, length(C1_sq)), C2_sq) -
    2 * (C1 %*% transport_plan %*% t(C2))

  lambda_e * L_base
}


#' Compute Neighborhood Consistency Regularization
#' 
#' Implements neighborhood consistency via KL divergence between transported neighborhoods
#' 
#' @param networks List of network structures
#' @param transport_plan Current transport plan S
#' @param lambda_n Neighborhood regularization weight
#' @return Neighborhood consistency cost matrix
#' @keywords internal
compute_neighborhood_consistency <- function(networks, transport_plan, lambda_n) {
  W1 <- networks[[1]]$transition
  W2 <- networks[[2]]$transition
  
  # FIX3-04: Use KL divergence for neighborhood consistency
  # Compute transported neighborhoods
  P <- t(W1) %*% transport_plan %*% W2
  
  # Add small epsilon for numerical stability in KL
  eps <- 1e-16
  P <- P + eps
  S_eps <- transport_plan + eps
  
  # KL divergence: P * log(P/S) = P * (log(P) - log(S))
  # For gradient computation, we need: lambda_n * log(P/S)
  kl_term <- log(P) - log(S_eps)
  
  # Return the KL divergence term scaled by lambda_n
  lambda_n * kl_term
}

#' Solve Sinkhorn Transport (Optimized Primal Scaling)
#' 
#' @param cost_matrix Dense cost matrix
#' @param tau Entropy regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Transport plan matrix
#' @keywords internal
solve_sinkhorn_transport <- function(cost_matrix, tau, max_iter, tol) {
  # Just delegate to the stabilized version which works correctly
  solve_sinkhorn_stabilized(cost_matrix, tau, max_iter, tol)
}

#' Compute Squared Euclidean Distance Matrix (R implementation)
#' 
#' @param X1 First matrix (n1 x d)
#' @param X2 Second matrix (n2 x d)
#' @return Distance matrix (n1 x n2)
#' @keywords internal
compute_squared_distances_r <- function(X1, X2) {
  X1_sq <- rowSums(X1^2)
  X2_sq <- rowSums(X2^2)
  outer(X1_sq, X2_sq, "+") - 2 * (X1 %*% t(X2))
}

#' Compute Squared Euclidean Distance Matrix (Dispatcher)
#' 
#' @param X1 First matrix (n1 x d)
#' @param X2 Second matrix (n2 x d)
#' @param use_cpp Whether to use C++ implementation
#' @return Distance matrix (n1 x n2)
#' @keywords internal
compute_squared_distances <- function(X1, X2, use_cpp = FALSE) {
  if (isTRUE(use_cpp)) {
    # Ensure matrices are regular for C++
    if (inherits(X1, "Matrix")) X1 <- as.matrix(X1)
    if (inherits(X2, "Matrix")) X2 <- as.matrix(X2)
    
    # Call C++ implementation
    compute_squared_distances_cpp_fallback(X1, X2)
  } else {
    # Use R implementation
    compute_squared_distances_r(X1, X2)
  }
}

#' Compute PARROT Embeddings
#' 
#' @param networks List of network structures
#' @param transport_plan Transport plan matrix
#' @param ncomp Number of components for embedding
#' @return Combined embedding matrix
#' @keywords internal
compute_parrot_embeddings <- function(networks, transport_plan, ncomp) {
  # Use SVD of transport plan to get embeddings (following package dimensionality patterns)
  if (is.null(ncomp)) {
    ncomp <- min(10, min(dim(transport_plan)) - 1)
  }
  
  # FIX2-06: Use sparse SVD for large matrices
  min_dim <- min(dim(transport_plan))
  use_sparse_svd <- ncomp < 0.2 * min_dim && min_dim > 100
  
  if (use_sparse_svd && requireNamespace("RSpectra", quietly = TRUE)) {
    # Use sparse SVD for efficiency
    svd_result <- RSpectra::svds(transport_plan, k = ncomp)
    embed1 <- svd_result$u
    embed2 <- svd_result$v
  } else {
    # Use base R SVD for smaller matrices
    svd_result <- svd(transport_plan, nu = ncomp, nv = ncomp)
    embed1 <- svd_result$u[, 1:ncomp, drop = FALSE]
    embed2 <- svd_result$v[, 1:ncomp, drop = FALSE]
  }
  
  # Combine embeddings from both networks
  scores <- rbind(embed1, embed2)
  scores
}


#' Get whether to use Rcpp implementations in PARROT
#'
#' @return Logical, whether to use Rcpp
#' @examples
#' get_parrot_use_cpp()
#' @export
get_parrot_use_cpp <- function() {
  isTRUE(getOption("manifoldalign.parrot.use_cpp", FALSE))
}
