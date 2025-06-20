#' Multi-Graph GRASP: Graph Alignment by Spectral Corresponding Functions
#'
#' Performs GRASP alignment on three or more graph domains simultaneously. Extends
#' the pairwise GRASP algorithm to handle multiple graphs through joint diagonalization
#' and shared latent basis alignment.
#'
#' @description
#' The multi-graph GRASP algorithm works by:
#' \enumerate{
#'   \item Computing spectral bases and descriptors for each graph
#'   \item Joint diagonalization: Finding rotations that simultaneously align all graphs
#'   \item Shared latent basis: Projecting all graphs to a common anchor space
#'   \item Computing assignments between graphs (typically to a reference)
#' }
#'
#' The method minimizes:
#' \deqn{\sum_{s=1}^{m} \mathrm{off}((M^{(s)})^T \Lambda^{(s)} M^{(s)}) + 
#'       \lambda \sum_{s<t} ||F^{(s)T} \Phi^{(s)} - F^{(t)T} \Phi^{(t)} M^{(t)} (M^{(s)})^T||_F^2}
#'
#' @param data Input data containing three or more graph domains. Can be:
#'   \itemize{
#'     \item A \code{hyperdesign} object with 3+ domains
#'     \item A list of 3+ matrices (nodes x features)
#'   }
#' @param ncomp Number of spectral components to extract (default: 30)
#' @param q_descriptors Number of diffusion-time descriptors (default: 100)
#' @param sigma Diffusion parameter for descriptor computation (default: 0.73)
#' @param lambda Regularization parameter for base alignment (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param anchor Index of reference graph or "mean" for barycenter (default: 1)
#' @param solver Assignment method: "auction" (default) or "linear"
#' @param max_iter Maximum iterations for joint alignment (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param ... Additional arguments (currently unused)
#'
#' @return A \code{grasp_multiset} object containing:
#' \itemize{
#'   \item \code{embeddings}: List of aligned embeddings for each graph
#'   \item \code{permutations}: List of permutation matrices to reference
#'   \item \code{rotations}: List of orthogonal rotation matrices
#'   \item \code{mapping_diag}: List of diagonal mapping matrices
#'   \item \code{anchor}: The anchor index or "mean"
#'   \item \code{iterations}: Number of iterations until convergence
#' }
#'
#' @examples
#' \donttest{
#' # Example with multiple synthetic graphs
#' set.seed(123)
#' graphs <- lapply(1:4, function(i) matrix(rnorm(50*3), 50, 3))
#' result <- grasp_multiset(graphs, ncomp = 20, q_descriptors = 50)
#' 
#' # Access results
#' length(result$embeddings)  # 4 aligned embeddings
#' dim(result$embeddings[[1]])  # 50 x 20
#' 
#' # Use mean as anchor instead of first graph
#' result_mean <- grasp_multiset(graphs, anchor = "mean")
#' }
#'
#' @references
#' Ovsjanikov, M., Ben-Chen, M., Solomon, J., Butscher, A., & Guibas, L. (2012). 
#' Functional maps: a flexible representation of maps between shapes. ACM 
#' Transactions on Graphics, 31(4), 1-11.
#'
#' @seealso \code{\link{grasp}} for pairwise alignment
#' @export
grasp_multiset <- function(data, ...) {
  UseMethod("grasp_multiset")
}

#' GRASP Multiset for Hyperdesign Objects
#'
#' @param data A hyperdesign object containing three or more graph domains
#' @param ncomp Number of spectral components (default: 30)
#' @param q_descriptors Number of descriptors (default: 100)
#' @param sigma Diffusion parameter (default: 0.73)
#' @param lambda Regularization parameter (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param anchor Index of reference graph (1-based) or "mean" (default: 1)
#' @param solver Assignment solver: "auction" or "linear" (default: "auction")
#' @param max_iter Maximum iterations for joint optimization (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param ... Additional arguments (currently unused)
#'
#' @return A grasp_multiset object
#' @export
#' @importFrom chk chk_s3_class chk_number chk_true chk_logical
#' @importFrom Matrix Diagonal crossprod norm
#' @importFrom clue solve_LSAP
grasp_multiset.hyperdesign <- function(data,
                                     ncomp = 30,
                                     q_descriptors = 100,
                                     sigma = 0.73,
                                     lambda = 0.1,
                                     use_laplacian = TRUE,
                                     anchor = 1L,
                                     solver = c("auction", "linear"),
                                     max_iter = 100,
                                     tol = 1e-6,
                                     ...) {
  
  # Validate input
  chk::chk_s3_class(data, "hyperdesign")
  
  domains <- unclass(data)
  m <- length(domains)
  if (m < 2) {
    stop("Need at least two domains for multi-graph alignment", call. = FALSE)
  }
  
  if (m == 2) {
    message("Only 2 domains provided. Consider using grasp() for pairwise alignment.")
  }
  
  # Validate parameters
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_logical(use_laplacian)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  # Validate anchor
  if (is.numeric(anchor)) {
    chk::chk_number(anchor)
    chk::chk_true(anchor >= 1 && anchor <= m)
    anchor <- as.integer(anchor)
  } else if (!identical(anchor, "mean")) {
    stop("anchor must be a domain index (1 to ", m, ") or 'mean'", call. = FALSE)
  }
  
  solver <- match.arg(solver)
  
  # Blocks A & B: Use existing helpers from grasp.R
  bases <- compute_grasp_basis(domains, ncomp, use_laplacian)
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)
  
  # New: Joint basis alignment
  alignment_result <- align_bases_multiset(bases, descriptors, lambda, max_iter, tol)
  rotations <- alignment_result$rotations
  
  # New: Diagonal maps to a latent anchor
  latent_coef <- build_anchor(descriptors, bases, rotations, anchor)
  C_list <- lapply(seq_len(m), function(s) {
    solve_diagonal_map(descriptors[[s]], bases[[s]], rotations[[s]], latent_coef)
  })
  
  # Embeddings in latent space
  embeddings <- lapply(seq_len(m), function(s) {
    as.matrix(bases[[s]]$vectors %*% rotations[[s]] %*% C_list[[s]])
  })
  
  # Assignments - compute permutations to anchor
  if (is.numeric(anchor)) {
    anchor_idx <- anchor
  } else {
    # For mean anchor, use first domain as reference for assignments
    anchor_idx <- 1L
  }
  
  permutations <- lapply(seq_len(m), function(s) {
    if (s == anchor_idx) {
      # Identity permutation for anchor
      seq_len(nrow(embeddings[[s]]))
    } else {
      solve_assignment_multiset(embeddings[[anchor_idx]], embeddings[[s]], 
                               distance = "cosine", solver = solver)
    }
  })
  
  # Return grasp_multiset object
  structure(
    list(
      embeddings = embeddings,
      permutations = permutations,
      rotations = rotations,
      mapping_diag = C_list,
      anchor = anchor,
      iterations = alignment_result$iterations,
      converged = alignment_result$converged,
      n_domains = m
    ),
    class = "grasp_multiset"
  )
}

#' GRASP Multiset for List of Matrices
#'
#' @param data A list containing three or more matrices (nodes x features)
#' @param ... Additional arguments passed to hyperdesign method
#'
#' @export
grasp_multiset.list <- function(data, ...) {
  
  # Validate input
  if (!is.list(data) || length(data) < 2) {
    stop("Provide a list with at least two matrices.", call. = FALSE)
  }
  
  if (!all(vapply(data, is.matrix, logical(1)))) {
    stop("Every element in the list must be a numeric matrix.", call. = FALSE)
  }
  
  # Validate matrix dimensions
  nrows <- vapply(data, nrow, integer(1))
  ncols <- vapply(data, ncol, integer(1))
  
  if (any(nrows < 3)) {
    stop("Each matrix must have at least 3 rows (nodes)", call. = FALSE)
  }
  
  if (any(ncols < 1)) {
    stop("Each matrix must have at least 1 column (feature)", call. = FALSE)
  }
  
  # Convert to stratum format (following cone_align pattern)
  strata <- lapply(seq_along(data), function(i) {
    list(
      x = data[[i]],
      design = data.frame(node_id = seq_len(nrow(data[[i]])))
    )
  })
  names(strata) <- paste0("domain_", seq_along(strata))
  
  # Wrap in hyperdesign-like structure
  class(strata) <- "hyperdesign"
  
  # Call hyperdesign method
  grasp_multiset.hyperdesign(strata, ...)
}

#' @export
grasp_multiset.default <- function(data, ...) {
  stop("No applicable method for grasp_multiset. data must be a hyperdesign ",
       "object or list of matrices with 2+ elements.", call. = FALSE)
}

#' Joint Basis Alignment for Multiple Graphs
#'
#' Performs joint diagonalization to find rotations that simultaneously
#' align all graph spectral bases.
#'
#' @param bases List of spectral bases from compute_grasp_basis
#' @param descriptors List of descriptor matrices
#' @param lambda Regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#'
#' @return List with rotations and convergence info
#' @keywords internal
align_bases_multiset <- function(bases, descriptors, lambda, 
                                max_iter = 100, tol = 1e-6) {
  
  m <- length(bases)
  k <- ncol(bases[[1]]$vectors)
  
  # Initialize rotations as identity matrices
  Ms <- replicate(m, Matrix::Diagonal(k), simplify = FALSE)
  
  # Pre-compute F^T * Phi for every graph
  FtPhi <- Map(function(d, b) Matrix::crossprod(d, b$vectors),
               descriptors, bases)
  
  # Objective function (closure uses outer environment)
  compute_objective <- function() {
    obj_total <- 0
    
    for (s in seq_len(m)) {
      # Off-diagonal penalty for graph s
      Lambda_s <- Matrix::Diagonal(x = bases[[s]]$values)
      aligned <- Matrix::t(Ms[[s]]) %*% Lambda_s %*% Ms[[s]]
      diag_aligned <- Matrix::diag(aligned)
      off_diag <- aligned - Matrix::Diagonal(x = diag_aligned)
      off_penalty <- Matrix::norm(off_diag, "F")^2
      
      # Descriptor alignment penalty
      desc_penalty <- 0
      for (t in setdiff(seq_len(m), s)) {
        residual <- FtPhi[[s]] - FtPhi[[t]] %*% Ms[[t]] %*% Matrix::t(Ms[[s]])
        desc_penalty <- desc_penalty + Matrix::norm(residual, "F")^2
      }
      
      obj_total <- obj_total + off_penalty + lambda * desc_penalty
    }
    
    obj_total
  }
  
  # Main optimization loop
  prev_obj <- compute_objective()
  
  for (iter in seq_len(max_iter)) {
    # Update each rotation matrix
    for (s in seq_len(m)) {
      # Compute gradient for M_s
      grad <- gradient_single_multiset(bases[[s]], FtPhi[[s]], FtPhi, Ms, s, lambda)
      
      # Line search with retraction
      step_size <- 0.5
      for (ls_iter in 1:8) {
        M_new_euc <- Ms[[s]] - step_size * grad
        
        # SVD retraction to Stiefel manifold
        M_new_mat <- as.matrix(M_new_euc)
        if (any(!is.finite(M_new_mat))) {
          warning("Non-finite values in rotation matrix, reducing step size")
          step_size <- step_size * 0.5
          next
        }
        
        svd_res <- svd(M_new_mat)
        M_new <- Matrix::Matrix(svd_res$u %*% t(svd_res$v), sparse = TRUE)
        
        # Store old M_s, update, check objective
        M_old <- Ms[[s]]
        Ms[[s]] <- M_new
        new_obj <- compute_objective()
        
        if (new_obj < prev_obj) {
          break
        } else {
          Ms[[s]] <- M_old
          step_size <- step_size * 0.5
        }
      }
    }
    
    # Check convergence
    cur_obj <- compute_objective()
    rel_change <- abs(prev_obj - cur_obj) / max(1, prev_obj)
    
    if (iter > 5 && rel_change < tol) {
      message("Joint alignment converged after ", iter, " iterations")
      return(list(rotations = Ms, iterations = iter, converged = TRUE))
    }
    
    prev_obj <- cur_obj
  }
  
  warning("Joint alignment reached maximum iterations without convergence")
  list(rotations = Ms, iterations = max_iter, converged = FALSE)
}

#' Gradient for Single Rotation in Multi-Graph Setting
#'
#' Computes the gradient of the objective with respect to one rotation matrix
#' in the context of multiple graphs.
#'
#' @param basis_s Spectral basis for graph s
#' @param FtPhi_s Descriptor-basis product for graph s  
#' @param FtPhi_all List of all descriptor-basis products
#' @param Ms List of all rotation matrices
#' @param s Index of current graph
#' @param lambda Regularization parameter
#'
#' @return Gradient matrix
#' @keywords internal
gradient_single_multiset <- function(basis_s, FtPhi_s, FtPhi_all, Ms, s, lambda) {
  
  m <- length(Ms)
  Lambda_s <- Matrix::Diagonal(x = basis_s$values)
  
  # Gradient of off-diagonal penalty (same as pairwise)
  A_s <- Lambda_s %*% Ms[[s]]
  AtMA <- Matrix::t(Ms[[s]]) %*% A_s
  diag_AtMA <- Matrix::diag(AtMA)
  grad_off <- 2 * (A_s %*% AtMA - A_s %*% Matrix::Diagonal(x = diag_AtMA))
  
  # Gradient of descriptor penalty (sum over all other graphs)
  grad_desc <- Matrix::Matrix(0, nrow = nrow(Ms[[s]]), ncol = ncol(Ms[[s]]), sparse = TRUE)
  
  for (t in setdiff(seq_len(m), s)) {
    # Derivative w.r.t M_s of ||F_s^T Phi_s - F_t^T Phi_t M_t M_s^T||^2
    residual <- FtPhi_s - FtPhi_all[[t]] %*% Ms[[t]] %*% Matrix::t(Ms[[s]])
    grad_desc <- grad_desc - 2 * Matrix::t(Ms[[t]]) %*% Matrix::t(FtPhi_all[[t]]) %*% residual
  }
  
  # Total Euclidean gradient
  euc_grad <- grad_off + lambda * grad_desc
  
  # Project to tangent space of Stiefel manifold
  riemannian_grad <- euc_grad - Ms[[s]] %*% Matrix::t(Ms[[s]]) %*% euc_grad
  
  riemannian_grad
}

#' Build Anchor Representation
#'
#' Constructs the anchor coefficient matrix either from a specific graph
#' or as the mean of all aligned representations.
#'
#' @param desc List of descriptor matrices
#' @param bases List of spectral bases
#' @param Mlist List of rotation matrices
#' @param anchor Index of anchor graph or "mean"
#'
#' @return Anchor coefficient matrix
#' @keywords internal
build_anchor <- function(desc, bases, Mlist, anchor) {
  if (length(anchor) == 1L && is.numeric(anchor)) {
    # Use specific graph as anchor
    s <- anchor
    Matrix::crossprod(desc[[s]], bases[[s]]$vectors %*% Mlist[[s]])
  } else {
    # Use mean as anchor
    anchor_list <- Map(function(d, b, M) {
      Matrix::crossprod(d, b$vectors %*% M)
    }, desc, bases, Mlist)
    
    Reduce(`+`, anchor_list) / length(bases)
  }
}

#' Solve Diagonal Map to Anchor
#'
#' Finds the diagonal matrix C that best maps from a graph's aligned
#' representation to the anchor representation.
#'
#' @param desc Descriptor matrix for graph
#' @param basis Spectral basis for graph
#' @param M Rotation matrix for graph
#' @param Z Anchor coefficient matrix
#'
#' @return Diagonal mapping matrix
#' @keywords internal
solve_diagonal_map <- function(desc, basis, M, Z) {
  A <- Matrix::crossprod(desc, basis$vectors %*% M)  # q x k
  
  # Ridge regularized solution for diagonal elements
  diag_vals <- Matrix::colSums(A * Z) / (Matrix::colSums(A^2) + 1e-9)
  
  Matrix::Diagonal(x = diag_vals)
}

#' Solve Assignment for Multi-Graph GRASP
#'
#' Computes node assignment between two embedded graphs using
#' cosine or Euclidean distance.
#'
#' @param E1 First embedding matrix
#' @param E2 Second embedding matrix
#' @param distance Distance metric: "cosine" or "euclidean"
#' @param solver Assignment solver: "auction" or "linear"
#'
#' @return Assignment vector
#' @keywords internal
solve_assignment_multiset <- function(E1, E2, distance = "cosine", solver = "auction") {
  
  # Compute cost matrix
  if (distance == "cosine") {
    # Normalize rows for cosine similarity
    E1_norms <- sqrt(rowSums(E1^2) + 1e-12)
    E2_norms <- sqrt(rowSums(E2^2) + 1e-12)
    E1_norm <- E1 / E1_norms
    E2_norm <- E2 / E2_norms
    
    sim_matrix <- E1_norm %*% t(E2_norm)
    cost <- 1 - sim_matrix
    cost[cost < 0] <- 0
  } else {
    # Euclidean distance
    E1_sq <- rowSums(E1^2)
    E2_sq <- rowSums(E2^2)
    cross <- E1 %*% t(E2)
    
    cost <- sqrt(outer(E1_sq, E2_sq, "+") - 2 * cross)
    cost[cost < 0] <- 0
  }
  
  # Solve assignment (following existing grasp.R pattern)
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }
  
  # Assignment solver selection
  # Note: clue::solve_LSAP may not accept method parameter in all versions
  assignment <- tryCatch({
    if (nrow(cost) >= 1000 && solver == "auction") {
      clue::solve_LSAP(cost, method = "auction")
    } else {
      clue::solve_LSAP(cost)  # Default solver
    }
  }, error = function(e) {
    # Fallback if method parameter is not supported
    clue::solve_LSAP(cost)
  })
  
  as.integer(assignment)
} 