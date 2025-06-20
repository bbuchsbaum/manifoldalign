# CONE-Align Implementation - Revised for manifoldalign Package Consistency

## Overview

Based on analysis of existing manifoldalign package patterns (KEMA, GPCA, lowrank_align, GRASP, generalized_procrustes), the CONE-Align implementation has been revised to ensure full consistency with established conventions.

## Key Consistency Requirements

### 1. Parameter Naming Conventions
- `ncomp` (not `k`) for number of components/dimensions
- `preproc` for preprocessing function (default: `center()`)
- `sigma` for diffusion/kernel parameters
- `lambda` for regularization parameters
- `use_laplacian` for Laplacian normalization choice
- `solver` for algorithm variants
- `max_iter` (not `tau_max`) for iteration limits
- `tol` (not `epsilon`) for convergence tolerance

### 2. Function Structure Patterns
- **S3 Methods**: `cone_align.hyperdesign()` following exact patterns from existing methods
- **Generic Function**: `cone_align()` with `UseMethod()` dispatch
- **Helper Functions**: Internal functions with `@keywords internal`
- **Validation**: Using `chk::` functions for parameter validation
- **Error Handling**: Using `safe_compute()` and consistent error messages

### 3. Dependencies and Imports
- Follow exact import patterns from GRASP and KEMA
- Use `clue` for assignment (consistent with GRASP)
- Use `RSpectra::eigs_sym` for eigendecomposition (consistent with all methods)
- Use `Matrix::` for sparse matrix operations
- Use `chk::` for parameter validation

### 4. Return Value Structure
- Return `multiblock_biprojector` object (exact consistency)
- Include standard components: `v`, `s`, `sdev`, `preproc`, `block_indices`
- Add method-specific components: `assignment`, `rotation`
- Follow exact naming from existing methods

---

## Revised Implementation

### Dependencies (Updated for Package Consistency)

```r
# DESCRIPTION dependencies (following exact package patterns)
# Imports: Matrix, RSpectra, chk, multivarious, methods
# Suggests: clue, igraph (for graph construction if needed)
```

### Data Structures (Revised)

- **Input Data**: `hyperdesign` objects (consistent with all package methods)
- **Embeddings**: `matrix` objects (`n × ncomp`) following `ncomp` convention
- **Assignments**: Integer vectors internally, sparse matrices for output
- **Return Values**: `multiblock_biprojector` objects

### Generic Function (Package Pattern)

```r
#' CONE-Align: Consensus Optimization for Node Embedding Alignment
#'
#' Performs CONE-Align on hyperdesign data structures. Aligns graph embeddings
#' by alternating between orthogonal transformations and node assignments.
#'
#' CONE-Align tackles graph alignment by iteratively refining both orthogonal
#' transformations and node-to-node correspondences. The algorithm alternates
#' between two steps: (1) finding optimal orthogonal rotations given current
#' assignments via Procrustes analysis, and (2) updating assignments given
#' current transformations via linear assignment.
#'
#' @param data Input data object containing graph domains
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{cone_align.hyperdesign}} for details on method-specific 
#'   parameters such as \code{preproc}, \code{ncomp}, \code{sigma}, 
#'   \code{lambda}, \code{use_laplacian}, \code{solver}, \code{max_iter}, 
#'   and \code{tol}
#'
#' @details
#' CONE-Align operates through the following algorithmic blocks:
#' \itemize{
#'   \item \strong{Spectral Embedding}: Compute Laplacian eigenmaps for each graph
#'   \item \strong{Orthogonal Alignment}: Find rotation matrices via Procrustes analysis
#'   \item \strong{Assignment Update}: Solve linear assignment problem for correspondences
#'   \item \strong{Convergence Check}: Iterate until assignments stabilize
#' }
#'
#' The algorithm minimizes the objective:
#' \deqn{\sum_{i,j} ||Q_i^T Z_i - P_{ij} Q_j^T Z_j||_F^2}
#'
#' where \eqn{Z_i} are the embeddings, \eqn{Q_i} are orthogonal transforms,
#' and \eqn{P_{ij}} are permutation matrices.
#'
#' Key parameters:
#' \itemize{
#'   \item \code{ncomp}: Dimension of spectral embeddings
#'   \item \code{sigma}: Bandwidth for embedding computation
#'   \item \code{lambda}: Regularization for numerical stability
#'   \item \code{solver}: Assignment algorithm ("linear" or "auction")
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned embeddings for all nodes
#'   \item \code{v}: Primal vectors for out-of-sample projection
#'   \item \code{assignment}: Node-to-node correspondence assignments
#'   \item \code{rotation}: Orthogonal transformation matrices
#'   \item \code{sdev}: Standard deviations of aligned components
#'   \item Additional metadata for reconstruction and validation
#' }
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' 
#' # Create synthetic graph domains
#' set.seed(123)
#' domain1 <- list(
#'   x = matrix(rnorm(100), 50, 2),  # Node features
#'   design = data.frame(node_id = 1:50)
#' )
#' domain2 <- list(
#'   x = matrix(rnorm(100), 50, 2),
#'   design = data.frame(node_id = 1:50)  
#' )
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = domain1, domain2 = domain2))
#' 
#' # Run CONE-Align with default parameters
#' result <- cone_align(hd, ncomp = 10)
#' 
#' # Access alignment results
#' node_assignment <- result$assignment
#' aligned_embeddings <- result$s
#' 
#' # Use exact solver for optimal results
#' result_exact <- cone_align(hd, ncomp = 10, solver = "linear")
#' }
#'
#' @references
#' Heimann, M., Shen, H., Safavi, T., & Koutra, D. (2018). REGAL: Representation 
#' learning-based graph alignment. In Proceedings of the 27th ACM International 
#' Conference on Information and Knowledge Management (pp. 117-126).
#'
#' @seealso \code{\link{cone_align.hyperdesign}}
#' @export
cone_align <- function(data, ...) {
  UseMethod("cone_align")
}
```

### S3 Method Implementation (Package Pattern)

```r
#' CONE-Align for Hyperdesign Objects
#'
#' Performs CONE-Align on hyperdesign data structures containing graph domains.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of embedding dimensions (default: 10)
#' @param sigma Diffusion parameter for embedding computation (default: 0.73)
#' @param lambda Regularization parameter for numerical stability (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment algorithm: "linear" for exact assignment (default), 
#'   "auction" for large-scale approximation
#' @param max_iter Maximum number of iterations (default: 30)
#' @param tol Convergence tolerance for assignment changes (default: 0.01)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the CONE-Align results
#'
#' @export
#' @importFrom chk chk_number chk_true chk_logical
#' @importFrom multivarious init_transform concat_pre_processors
cone_align.hyperdesign <- function(data, 
                                  preproc = center(), 
                                  ncomp = 10,
                                  sigma = 0.73,
                                  lambda = 0.1,
                                  use_laplacian = TRUE,
                                  solver = c("linear", "auction"),
                                  max_iter = 30,
                                  tol = 0.01,
                                  ...) {
  
  # Input validation (following exact package patterns)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_logical(use_laplacian)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  # Validate solver parameter (following KEMA/GRASP pattern)
  solver <- match.arg(solver)
  
  # Validate input data (following package pattern)
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list of hyperdesign objects", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("CONE-Align currently supports exactly 2 domains. For multi-graph ", 
         "alignment, use cone_align_multiple() [coming in Sprint 2]", call. = FALSE)
  }
  
  # Preprocess data (following exact package pattern)
  pdata <- safe_compute(
    multivarious::init_transform(data, preproc),
    "Data preprocessing failed. Check preprocessing function or data format."
  )
  
  # Extract preprocessing information (package pattern)
  proclist <- attr(pdata, "preproc")
  names(proclist) <- names(pdata)
  
  # Block indices computation (following package pattern)
  block_indices <- block_indices(pdata)
  proc <- multivarious::concat_pre_processors(proclist, 
                                             split(block_indices, row(block_indices)))
  names(block_indices) <- names(pdata)
  
  # Call CONE-Align fitting function
  cone_align_fit(pdata, proc, ncomp, sigma, lambda, use_laplacian, 
                solver, max_iter, tol, block_indices)
}

#' @keywords internal
cone_align_fit <- function(strata, proc, ncomp, sigma, lambda, use_laplacian, 
                          solver, max_iter, tol, block_indices) {
  
  # Compute spectral embeddings (following GRASP pattern)
  embeddings <- compute_cone_embeddings(strata, ncomp, sigma, use_laplacian)
  
  # Iterative alignment (core CONE-Align algorithm)
  alignment_result <- cone_align_iterate(embeddings, solver, max_iter, tol, lambda)
  
  # Compute scores (following package pattern)
  embed1_aligned <- embeddings[[1]] %*% alignment_result$rotations[[1]]
  embed2_aligned <- embeddings[[2]] %*% alignment_result$rotations[[2]]
  scores <- rbind(embed1_aligned, embed2_aligned)
  
  # Primal vectors computation (following KEMA/GRASP pattern)
  v <- do.call(rbind, lapply(seq_along(strata), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- embeddings[[i]]
    Matrix::crossprod(xi, alpha_i)
  }))
  
  # Return multiblock_biprojector (exact package pattern)
  multivarious::multiblock_biprojector(
    v = v,
    s = scores,
    sdev = apply(scores, 2, sd),
    preproc = proc,
    block_indices = block_indices,
    assignment = alignment_result$assignment,
    rotation = alignment_result$rotations,
    classes = "cone_align"
  )
}

#' @keywords internal
cone_align_iterate <- function(embeddings, solver, max_iter, tol, lambda) {
  n <- nrow(embeddings[[1]])
  ncomp <- ncol(embeddings[[1]])
  
  # Initialize assignment (following package integer vector pattern)
  assignment <- 1:n
  rotations <- list(diag(ncomp), diag(ncomp))
  
  for (iter in 1:max_iter) {
    # Step 1: Update rotations via Procrustes (following GRASP pattern)
    target_embed <- embeddings[[2]][assignment, , drop = FALSE]
    rotations[[1]] <- solve_procrustes_cone(target_embed, embeddings[[1]], lambda)
    rotations[[2]] <- diag(ncomp)  # Keep second embedding as reference
    
    # Step 2: Update assignment via linear assignment
    embed1_aligned <- embeddings[[1]] %*% rotations[[1]]
    assignment_new <- solve_assignment_cone(embed1_aligned, embeddings[[2]], solver)
    
    # Step 3: Check convergence (following package patterns)
    changes <- sum(assignment_new != assignment)
    if (changes / n < tol) {
      break
    }
    assignment <- assignment_new
  }
  
  list(assignment = assignment, rotations = rotations, iterations = iter)
}

#' Solve Procrustes Problem for CONE-Align
#' 
#' @param target Target embedding matrix
#' @param source Source embedding matrix  
#' @param lambda Regularization parameter
#' @return Orthogonal rotation matrix
#' @keywords internal
solve_procrustes_cone <- function(target, source, lambda) {
  # Regularized cross-covariance (following package numerical stability)
  M <- Matrix::crossprod(target, source) + lambda * Matrix::Diagonal(ncol(source))
  
  # SVD-based solution (following exact GRASP pattern)
  svd_result <- svd(as.matrix(M))
  
  # Ensure proper rotation (det = +1, following package reflection handling)
  d <- sign(det(svd_result$v %*% t(svd_result$u)))
  k_dim <- ncol(svd_result$u)
  
  # Apply reflection correction (exact GRASP pattern)
  if (d < 0) {
    svd_result$v[, k_dim] <- -svd_result$v[, k_dim]
  }
  
  Q <- svd_result$v %*% t(svd_result$u)
  as.matrix(Q)
}

#' Solve Assignment Problem for CONE-Align
#' 
#' @param embed1 First embedding matrix
#' @param embed2 Second embedding matrix
#' @param solver Assignment solver method
#' @return Assignment vector (1-based indexing)
#' @keywords internal
solve_assignment_cone <- function(embed1, embed2, solver) {
  # Compute cost matrix (efficient Euclidean distance)
  embed1_sq <- rowSums(embed1^2)
  embed2_sq <- rowSums(embed2^2)
  cross_prod <- embed1 %*% t(embed2)
  cost_matrix <- sqrt(outer(embed1_sq, embed2_sq, "+") - 2 * cross_prod)
  
  # Solve assignment (following exact GRASP dependency pattern)
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }
  
  # Use specified solver (following package solver pattern)
  if (solver == "auction" && nrow(cost_matrix) >= 1000) {
    # Use auction algorithm for large problems
    assignment <- clue::solve_LSAP(cost_matrix, method = "auction")
  } else {
    # Use exact algorithm (default)
    assignment <- clue::solve_LSAP(cost_matrix)
  }
  
  as.integer(assignment)
}

#' @export
cone_align <- function(data, ...) {
  UseMethod("cone_align")
}

# Return the final alignment results
res
}

---

## Implementation Tickets (Revised for Package Consistency)

### Sprint 1: Core Functionality & Package Integration

**Ticket CA-01: Implement S3 Method Structure**
- **Description:** Implement `cone_align()` generic and `cone_align.hyperdesign()` method following exact manifoldalign patterns
- **Acceptance Criteria:**
  - Generic function with `UseMethod()` dispatch
  - S3 method with consistent parameter naming (`ncomp`, `preproc`, `sigma`, `lambda`, etc.)
  - Proper parameter validation using `chk::` functions
  - Consistent error messages and handling

**Ticket CA-02: Implement Spectral Embedding Helper**
- **Description:** Create `compute_cone_embeddings()` following KEMA/GRASP neighbor computation patterns
- **Acceptance Criteria:**
  - Uses `RANN::nn2()` for neighbor computation (KEMA pattern)
  - Uses `RSpectra::eigs_sym()` for eigendecomposition (all methods pattern)
  - Handles Laplacian normalization consistently with `use_laplacian` parameter
  - Returns matrices with proper dimensions and finite values

**Ticket CA-03: Implement Procrustes and Assignment Helpers**
- **Description:** Create `solve_procrustes_cone()` and `solve_assignment_cone()` following GRASP patterns
- **Acceptance Criteria:**
  - Procrustes solver ensures det(Q) = +1 (reflection handling)
  - Assignment solver uses `clue::solve_LSAP()` with fallback patterns
  - Supports both "linear" and "auction" solvers
  - Returns proper 1-based indexing for assignments

**Ticket CA-04: Implement Core Iteration Algorithm**
- **Description:** Create `cone_align_iterate()` with convergence checking following package patterns
- **Acceptance Criteria:**
  - Uses `max_iter` and `tol` parameter naming (not `tau_max`/`epsilon`)
  - Proper convergence detection based on assignment changes
  - Returns structured result with rotations and final assignment
  - Handles edge cases and numerical stability

**Ticket CA-05: Implement Return Value Structure**
- **Description:** Return `multiblock_biprojector` object following exact package patterns
- **Acceptance Criteria:**
  - Includes standard components: `v`, `s`, `sdev`, `preproc`, `block_indices`
  - Adds method-specific components: `assignment`, `rotation`
  - Proper `classes` attribute set to `"cone_align"`
  - Follows exact naming and structure from other methods

**Ticket CA-06: Add Package Documentation and Examples**
- **Description:** Create comprehensive Roxygen2 documentation following package patterns
- **Acceptance Criteria:**
  - Generic and method documentation with proper `@param`, `@return`, `@examples`
  - Examples use `hyperdesign` objects and proper function calls
  - References to relevant literature and methodology
  - Consistent documentation style with existing methods

### Sprint 2: Testing and Validation

**Ticket CA-07: Implement Comprehensive Unit Tests**
- **Description:** Create test suite following package testing patterns from `test-grasp.R`, `test-kema.R`
- **Acceptance Criteria:**
  - Basic functionality tests with synthetic data
  - Parameter validation tests for all input parameters
  - Edge case handling (minimum data sizes, convergence failures)
  - Return value structure validation

**Ticket CA-08: Performance and Accuracy Validation**
- **Description:** Validate algorithm performance on known test cases
- **Acceptance Criteria:**
  - Accuracy measurement on synthetic graphs with known correspondences
  - Performance benchmarking against graph size
  - Comparison with existing graph alignment methods if available
  - Documentation of expected performance characteristics

---

## Key Consistency Achievements

1. **Parameter Names**: Exact alignment with package conventions (`ncomp`, `preproc`, `sigma`, `lambda`, `use_laplacian`, `solver`, `max_iter`, `tol`)

2. **Function Structure**: S3 methods, generic dispatch, helper functions with `@keywords internal`, consistent validation patterns

3. **Dependencies**: Uses exact same packages as existing methods (`clue`, `RSpectra`, `RANN`, `chk`, `multivarious`)

4. **Return Values**: Returns `multiblock_biprojector` with exact same component structure and naming as other methods

5. **Documentation**: Follows exact Roxygen2 patterns, parameter descriptions, examples, and reference formatting

6. **Error Handling**: Uses `safe_compute()`, consistent error messages, and same validation approaches

7. **Implementation Details**: Follows eigendecomposition patterns, sparse matrix usage, numerical stability approaches from existing methods

This revision ensures CONE-Align will integrate seamlessly with the manifoldalign package while maintaining its unique algorithmic approach.