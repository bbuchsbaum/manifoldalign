#' Generalized Orthogonal Procrustes with Partial Tasks
#'
#' @description
#' Extends Generalized Orthogonal Procrustes to handle partial task 
#' observations. Uses efficient sparse matrix operations and robust 
#' initialization.
#'
#' @param A_list A list of length \eqn{n}, where each element \eqn{A_i} is a 
#'   \eqn{d x li} numeric matrix representing the observed data for subject 
#'   \eqn{i}. Each subject might have different \eqn{li < L}.
#' @param task_labels_list A list of length \eqn{n}, where 
#'   \eqn{task_labels_list[[i]]} is an integer (or numeric) vector of length 
#'   \eqn{li}, specifying which global tasks each column of \eqn{A_i} 
#'   corresponds to (in \eqn{1..L}).
#' @param L Integer. The total number of distinct tasks across all subjects.
#' @param max_iter Integer. Maximum number of GPM iterations. Default 100.
#' @param tol Numeric tolerance for convergence. Default \code{1e-6}.
#' @param tol_type Character. Type of tolerance ("relative" or "absolute"). 
#'   Note: The relative tolerance is calculated as `sqrt(diff_sq) / norm0`, 
#'   where `norm0` is the Frobenius norm of the initial `O_mats` (plus epsilon 
#'   for stability). The absolute tolerance uses `sqrt(diff_sq) / (sqrt(n*d) + 1)`.
#'   Default "relative".
#' @param verbose Logical. If \code{TRUE}, prints iteration info. Default 
#'   \code{FALSE}.
#' @param svd_method Backend for the initial SVD (one of \code{"irlba"}, 
#'   \code{"base"}). \code{"irlba"} is recommended and used by default, with 
#'   fallback to random initialization if needed. \code{"base"} forces the use 
#'   of dense SVD, which can be slow/memory-intensive. Default \code{"irlba"}.
#' @param svd_opts List of extra arguments passed to \code{irlba::irlba} or 
#'   \code{base::svd}.
#'
#' @details
#' \strong{Sparse Construction of D}:
#' We form a large sparse matrix \eqn{D \in R^{(n d) \times L}} using 
#' \code{Matrix::sparseMatrix}. Triplet indices (row, column, value) are 
#' collected efficiently.
#'
#' \strong{Robust Initialization}:
#' Initialization uses the top \eqn{d} left singular vectors of \eqn{D}, typically 
#' computed via \code{irlba::irlba}. If \code{irlba} fails or returns a rank 
#' deficient result (less than \eqn{d} vectors), it falls back to initializing 
#' with random orthogonal matrices to ensure robustness.
#'
#' \strong{Vectorized GPM Iteration}:
#' The core iteration avoids \eqn{O(n^2)} complexity by computing the update 
#' term \eqn{CO} using two efficient sparse matrix multiplications:
#' \enumerate{
#'   \item \eqn{Y <- D^T * stack(O_{old})}
#'   \item \eqn{CO_{stack} <- D * Y}
#' }
#' where \eqn{stack(O)} reshapes the \eqn{d \times d \times n} array of 
#' rotation matrices into an \eqn{(nd) \times d} matrix. The resulting 
#' \eqn{(nd) \times d} matrix \eqn{CO_{stack}} is then unstacked back into 
#' a \eqn{d \times d \times n} array.
#'
#' \strong{Efficient Projection}:
#' Each \eqn{d \times d} slice of the unstacked update is projected onto the 
#' orthogonal group O(d). An efficient projection method (\code{.projO}) is used, 
#' switching between an eigenvalue-based method for small \eqn{d} and an 
#' SVD-based method for larger \eqn{d}.
#'
#' \strong{Convergence Check}:
#' Convergence is assessed by comparing the Frobenius norm of the difference 
#' between consecutive projected iterates \eqn{O_{new}} and 
#' \eqn{O_{old}}. For relative tolerance, this difference is scaled by the 
#' Frobenius norm of the initial \eqn{O} matrices (calculated once).
#'
#' \strong{Final Consensus A_est}:
#' After convergence, the consensus matrix \eqn{A_{est}} is computed. 
#' For each task \eqn{j}, we average the transformed data 
#' \eqn{O_i^T A_i(*, j_i)} only over subjects \eqn{i} in \eqn{S_j} that observed 
#' task \eqn{j}. Task counts (number of subjects observing each task \eqn{j}) are 
#' computed by directly iterating over \code{task_labels_list}, ensuring correctness 
#' even if input matrices contain explicit zeros.
#'
#' @return A list with:
#' \item{O_mats}{List of \eqn{n} orthogonal \eqn{d \times d} matrices.}
#' \item{A_est}{A \eqn{d \times L} consensus matrix. Columns corresponding to 
#'   tasks observed by no subjects will contain \code{NA}.}
#' \item{iterations}{Number of iterations performed.}
#' \item{converged}{Logical indicating whether the GPM converged.}
#' \item{final_diff}{The final difference value (`sqrt(diff_sq)`) used for the 
#'   convergence check.}
#'
#' @examples
#' \donttest{
#' # Simple example with 2 subjects and 3 tasks
#' A1 <- matrix(rnorm(12), 3, 4)  # Subject 1: 3 features x 4 tasks
#' A2 <- matrix(rnorm(9), 3, 3)   # Subject 2: 3 features x 3 tasks
#' 
#' # Task labels (partial overlap)
#' task_labels1 <- c(1, 2, 3, 4)  # Subject 1 observes tasks 1,2,3,4
#' task_labels2 <- c(1, 3, 5)     # Subject 2 observes tasks 1,3,5
#' 
#' # Run alignment
#' result <- generalized_procrustes(
#'   A_list = list(A1, A2),
#'   task_labels_list = list(task_labels1, task_labels2),
#'   L = 5,  # Total of 5 unique tasks
#'   max_iter = 50,
#'   verbose = FALSE
#' )
#' 
#' # Check results
#' print(result$converged)
#' print(dim(result$A_est))  # Should be 3x5
#' }
#'
#' @importFrom Matrix sparseMatrix crossprod t
#' @importFrom stats rnorm
#' @importFrom rlang enquo !! as_name quo_is_symbol eval_tidy
#' @importFrom purrr map
#' @importFrom dplyr select pull
#' @importFrom multivarious init_transform center
#' @importFrom irlba irlba
#' @export
generalized_procrustes <- function(
    A_list,
    task_labels_list,
    L,
    max_iter   = 100,
    tol        = 1e-6,
    tol_type   = c("relative", "absolute"),
    verbose    = FALSE,
    svd_method = c("irlba", "base"),
    svd_opts   = list()
){
  tol_type   <- match.arg(tol_type)
  svd_method <- match.arg(svd_method)

  # --- Helper Functions --- Internal to this function
  # EFFICIENCY FIX: Vectorized stack/unstack using aperm for O(1) reshaping
  # Stack (d x d x n) array into (nd x d) matrix
  .stack_rot <- function(Oarr) {          # (d x d x n) -> (nd x d)
    # Vectorized approach: reshape using aperm (array permutation)
    matrix(aperm(Oarr, c(3, 1, 2)), nrow = dim(Oarr)[3] * dim(Oarr)[1], ncol = dim(Oarr)[2])
  }

  # Unstack (nd x d) matrix back to (d x d x n) array  
  .unstack <- function(M, d_, n_) {         # inverse of .stack_rot
    # Vectorized approach: reshape and permute back
    aperm(array(M, c(n_, d_, d_)), c(2, 3, 1))
  }

  # Efficient projection onto O(d)
  .projO <- function(A) {                 # => nearest orthogonal matrix
    d_ <- ncol(A)
    if(d_ == 0) return(A) # Handle edge case
    # Check for NaNs/Infs which cause svd/eigen to fail
    if(any(!is.finite(A))) {
        warning("Non-finite values encountered in matrix projection; returning identity.", call. = FALSE)
        return(diag(d_))
    }
    if (d_ < 8) {                    # Small d: Polar decomposition via eigen is often faster
      # A = UP => U = A P^-1 = A (A^T A)^(-1/2)
      # Use crossprod for A^T A
      AtA <- crossprod(A)
      # Ensure symmetry for eigen, although crossprod should guarantee it
      # Add small epsilon for numerical stability? Maybe not needed if AtA is well-behaved.
      ev <- eigen(AtA, symmetric = TRUE)
      # Check for negative eigenvalues due to numerical issues?
      ev$values[ev$values < 0] <- 0 # Clamp small negative values
      # P^(-1/2) = V D^(-1/2) V^T
      # ROBUSTNESS FIX: Guard against overflow with stronger clamping for singular AtA
      P_inv_sqrt <- ev$vectors %*% diag(1/sqrt(pmax(ev$values, 1e-8)), nrow=d_) %*% t(ev$vectors)
      O <- A %*% P_inv_sqrt
    } else { # Larger d: SVD-based projection U V^T is generally robust
      # A = U S V^T => U V^T is the orthogonal part
      # svd() might be faster than svds() from Rsvd/PRIMME for moderate d
      # We only need U and V
      sv <- svd(A, nu = d_, nv = d_)
      O <- sv$u %*% t(sv$v)
    }
    # Optional check for orthogonality? Or assume it's close enough.
    # sum((crossprod(O) - diag(d_))^2) should be small
    O
  }
  # --- End Helper Functions ---

  n <- length(A_list)
  # if(n < 1) stop("Must have at least 1 subject.") # Allow n=1? GPM needs n>=2.
  if(n < 2) stop("Must have at least 2 subjects for Procrustes alignment.")
  if(length(task_labels_list) != n){
    stop("task_labels_list must have the same length as A_list.")
  }

  # --- Input Validation & Dimension Check ---
  d <- NA_integer_
  # total_nnz calculation is okay, but we build triplets differently now
  for(i in seq_len(n)){
    Ai <- A_list[[i]]
    if (!is.matrix(Ai) || !is.numeric(Ai)) stop("A_list elements must be numeric matrices.")
    li <- ncol(Ai)
    current_d <- nrow(Ai)

    if (is.na(d)) {
      d <- current_d
    } else if (current_d != d) {
      stop("All subject matrices must have the same row dimension d.")
    }
    if (d <= 0) stop("Dimension d must be positive.")

    tLabels <- task_labels_list[[i]]
    if (length(tLabels) != li) {
      stop(sprintf("task_labels_list[[%d]] length (%d) must match ncol(A_list[[%d]]) (%d)",
                   i, length(tLabels), i, li))
    }
    if (li > 0) { # Only check labels if there are columns
        if (!is.numeric(tLabels) || any(tLabels < 1) || any(tLabels > L) || any(tLabels != floor(tLabels))) {
           stop(sprintf("task_labels_list[[%d]] must contain integers between 1 and L=%d.", i, L))
        }
        if (any(duplicated(tLabels))) {
            stop(sprintf("task_labels_list[[%d]] contains duplicate task labels.", i))
        }
    }
  }
  if (is.na(d)) stop("Cannot determine dimension d (list might be empty or contain non-matrices).")


  # ---- 1. Build sparse D (nd x L) & Calculate Subject Counts per Task ----
  
  # EFFICIENCY FIX: Pre-calculate total non-zeros to avoid quadratic concatenation
  total_nnz <- sum(vapply(A_list, function(Ai) nrow(Ai) * ncol(Ai), integer(1)))
  
  # Pre-allocate triplet vectors (MAJOR PERFORMANCE IMPROVEMENT)
  row_inds <- integer(total_nnz)
  col_inds <- integer(total_nnz)
  vals     <- numeric(total_nnz)
  
  row_offset <- 0L
  current_idx <- 1L  # Current position in pre-allocated vectors
  
  # --- Correctly count subjects observing each task --- VITAL for averaging!
  subject_counts_per_task <- integer(L) # Initialize counts to zero

  for (i in seq_len(n)) {
    Ai <- A_list[[i]]
    li <- ncol(Ai)
    tLabels <- task_labels_list[[i]] # Get labels even if li=0 (should be empty)

    # Increment count for tasks observed by this subject
    if(length(tLabels) > 0) {
        subject_counts_per_task[tLabels] <- subject_counts_per_task[tLabels] + 1L
    }

    # Build sparse matrix triplets (only if data exists)
    if (li > 0) {
        nnz_i <- d * li
        end_idx <- current_idx + nnz_i - 1L
        
        # Fill pre-allocated vectors (no concatenation = O(1) amortized)
        row_inds[current_idx:end_idx] <- row_offset + rep(seq_len(d), times = li)
        col_inds[current_idx:end_idx] <- rep(tLabels, each = d)
        vals[current_idx:end_idx] <- as.vector(Ai)  # Ensure vector for assignment
        
        current_idx <- end_idx + 1L
    }
    row_offset <- row_offset + d
  }
  
  # Trim vectors to actual used length (in case some matrices were empty)
  if (current_idx <= total_nnz) {
    row_inds <- row_inds[1:(current_idx-1L)]
    col_inds <- col_inds[1:(current_idx-1L)]
    vals <- vals[1:(current_idx-1L)]
  }

  # Check if any data was actually added to build D
  if (length(vals) == 0) {
      # Should not happen if n>=2 and validation passed, unless all A_list[[i]] were empty
      warning("Constructed sparse matrix D is empty. Check input A_list.")
      # Return something reasonable? Or let sparseMatrix handle it?
      # Let's construct it, subsequent steps might fail gracefully or produce NAs.
      # An empty D might cause issues with SVD/irlba though.
      # If no subjects observed any tasks, subject_counts_per_task is all zero,
      # A_est will be all NAs, O_mats will be random/identity based on init.
  }

  D <- Matrix::sparseMatrix(
    i = row_inds,
    j = col_inds,
    x = vals,
    dims = c(n*d, L)
    # Consider dims = c(n*d, L), index1=TRUE? Default is 1-based.
    # giveCsparse = ? Performance varies.
  )

  # ---- 2. Initialize Rotations O_i with SVD of D (Robustly) ----
  if (verbose) message("Initializing O_mats via SVD of D...")

  O_old <- array(0.0, dim = c(d, d, n))

  # Handle case where D might be empty or have zero columns/rows relevant to SVD
  svd_successful <- FALSE
  if (nrow(D) > 0 && ncol(D) > 0 && Matrix::nnzero(D) > 0) {
      if (svd_method == "base") {
          if (verbose) message("Using base::svd for initialization (can be slow/memory-intensive).")
          # Ensure D is not excessively large before trying to make dense
          # Heuristic: if estimated dense size > 1GB? Or just let it try?
          tryCatch({
              M_dense <- as.matrix(D) # Warning: Potential memory issue!
              args <- c(list(x = M_dense, nu = d, nv = 0), svd_opts)
              sres <- do.call(svd, args)
              U_init <- sres$u
              if (ncol(U_init) < d) stop("base::svd returned rank < d")
              svd_successful <- TRUE
          }, error = function(e) {
              warning("base::svd failed during initialization: ", e$message, ". Falling back to random orthogonal matrices.", call. = FALSE)
              U_init <- NULL # Signal failure
          })
      } else { # Default: try irlba, fallback to random
          if (verbose) message("Attempting irlba::irlba for initialization...")
          if (!requireNamespace("irlba", quietly = TRUE)) {
            stop("Package 'irlba' needed for svd_method='irlba'. Please install it.", call. = FALSE)
          }
          args <- c(list(A = D, nu = d, nv = 0), svd_opts)
          args$center <- NULL; args$scale <- NULL # Avoid conflicts
          U_init <- tryCatch({
              sres <- do.call(irlba::irlba, args)
              if (ncol(sres$u) < d) stop("irlba returned rank < d")
              svd_successful <- TRUE
              sres$u
          }, error = function(e) {
              warning("irlba::irlba failed: ", e$message, ". Falling back to random orthogonal matrices.", call. = FALSE)
              NULL # Signal failure
          })
      }
  } else {
      if (verbose) message("Sparse matrix D is empty or effectively zero; cannot perform SVD.")
      U_init <- NULL # Ensure fallback occurs
  }

  # If SVD failed, was rank deficient, or D was empty, use random orthogonal matrices
  if (!svd_successful) {
    if (verbose) message("Using random orthogonal initialization.")
    # Generate n*d x d standard normal matrix, then get Q from QR
    rand_mat <- matrix(stats::rnorm(n * d * d), nrow = n * d, ncol = d)
    # Ensure QR doesn't fail on zero matrix etc.
    qr_decomp <- tryCatch(qr(rand_mat), error = function(e) NULL)
    if(is.null(qr_decomp)) stop("QR decomposition failed during random initialization.")
    U_init <- qr.Q(qr_decomp)

    # Ensure dimensions are correct after QR
    if (is.null(U_init) || nrow(U_init) != n*d || ncol(U_init) < d) {
        stop("Random initialization failed to produce a valid matrix.")
    }
    # Take only the first d columns if QR produced more
    if (ncol(U_init) > d) U_init <- U_init[, 1:d, drop = FALSE]
  }

  # Project each subject's block onto O(d) using the efficient projector
  for (i in seq_len(n)) {
    row_idx <- ((i-1)*d + 1):(i*d)
    # Check if U_init has enough rows (should be n*d)
    if(max(row_idx) > nrow(U_init)) stop("Internal error: U_init row dimension mismatch.")
    Ui <- U_init[row_idx, , drop = FALSE]
    if (ncol(Ui) != d) stop(sprintf("Internal error: Initial block for subject %d has %d columns, expected %d.", i, ncol(Ui), d))
    O_old[,,i] <- .projO(Ui)
  }

  # ---- 3. Iterative GPM (Vectorized & Corrected Convergence) ----
  if (verbose) message("Starting GPM iterations...")

  # Calculate initial norm for relative tolerance denominator (once)
  # ROBUSTNESS FIX: Add floor to prevent hypersensitive relative test if random fallback produces singular O_old
  norm0 <- max(sqrt(sum(O_old^2)), 1e-8)
  converged  <- FALSE
  final_diff <- NA_real_
  iter <- 0 # Initialize iter outside loop

  for (iter in seq_len(max_iter)) {

    # Stack current rotations
    O_stack_old <- .stack_rot(O_old)

    # Compute Y = D^T * O_stack_old  (L x d)
    # Use Matrix::crossprod for efficient t(sparse) %*% dense
    Y <- Matrix::crossprod(D, O_stack_old)

    # Compute CO_stack = D * Y  ((nd) x d)
    CO_stack <- D %*% Y

    # Unstack CO_stack into an array (d x d x n)
    if(!isTRUE(all.equal(dim(CO_stack), c(n*d, d)))) {
        stop(sprintf("Internal error: CO_stack dimension mismatch. Iter %d. Expected (%d, %d), got (%d, %d)",
                     iter, n*d, d, nrow(CO_stack), ncol(CO_stack)))
    }
    CO_arr <- .unstack(CO_stack, d, n) # Unprojected updates

    # --- Project each CO_i onto O(d) --- Store in O_new
    O_new <- array(0.0, dim = c(d, d, n)) # Initialize storage for new projections
    for (i in seq_len(n)) {
      O_new[,,i] <- .projO(CO_arr[,,i])
    }

    # --- Calculate difference between NEW projection and OLD projection ---
    # Ensure dimensions match before subtraction
    if(!isTRUE(all.equal(dim(O_new), dim(O_old)))) stop("Internal error: O_new/O_old dimension mismatch.")
    diff_sq <- sum((O_new - O_old)^2)
    diff_val <- sqrt(diff_sq)
    final_diff <- diff_val # Store last difference value

    # --- Check Convergence ---
    converged <- FALSE
    if (tol_type == "relative") {
      rel_val <- diff_val / norm0 # Use frozen norm0
      if (!is.na(rel_val) && rel_val < tol) converged <- TRUE
      # Avoid string formatting cost if not verbose
      if (verbose && (iter %% 5L == 0L || converged || iter == 1L || iter == max_iter)) {
          message(sprintf("Iter %3d: Diff = %.6e, Rel Diff = %.6e", iter, diff_val, rel_val))
      }
    } else { # Absolute tolerance, scaled by sqrt(n*d)
      abs_val_scaled <- diff_val / (sqrt(n * d) + 1)
      if (!is.na(abs_val_scaled) && abs_val_scaled < tol) converged <- TRUE
      # Avoid string formatting cost if not verbose
      if (verbose && (iter %% 5L == 0L || converged || iter == 1L || iter == max_iter)) {
          message(sprintf("Iter %3d: Diff = %.6e, Scaled Abs Diff = %.6e", iter, diff_val, abs_val_scaled))
      }
    }

    # Update O_old for the next iteration
    O_old <- O_new

    if (converged) {
      if (verbose) message("Convergence criteria met.")
      break
    }

  } # End GPM iteration loop

  if (iter == max_iter && !converged && verbose) {
      message("Maximum iterations reached without meeting convergence criteria.")
  }

  O_final <- O_old # Use the last computed O (which has been projected)

  # ---- OPTIONAL: Tightness Certificate Check (Lambda-C test from Ling 2024, eqs. 3.15-3.17) ----
  # This provides a certificate of global optimality beyond just convergence
  if (verbose) {
    tryCatch({
      # Compute C = D D^T (correlation matrix)
      C_matrix <- Matrix::tcrossprod(D)  # D %*% t(D) but more efficient
      
      # Stack final rotations for eigenvalue computation
      O_stack_final <- .stack_rot(O_final)
      
          # Compute Lambda as the eigenvalues of the final generalized eigenvalue problem
    # For GPM: Lambda O = C O, so Lambda = O^T C O (Rayleigh quotients)
      CO_final <- C_matrix %*% O_stack_final
      lambda_vals <- colSums(O_stack_final * CO_final)  # Diagonal of O^T C O
      
      # Tightness check: ||CO - Lambda*O||_F should be small for global optimum
      lambda_O <- O_stack_final %*% diag(lambda_vals, nrow = d)
      tightness_error <- norm(CO_final - lambda_O, "F") / max(norm(CO_final, "F"), 1e-12)
      
      message(sprintf("Tightness certificate: ||CO - Lambda*O||_F / ||CO||_F = %.2e", tightness_error))
      
      if (tightness_error < sqrt(tol)) {
              message("(checkmark) Global optimality certificate satisfied (tight relaxation)")
    } else {
      message("(warning) Tightness error above threshold - solution may be local optimum")
      }
    }, error = function(e) {
      message("Tightness certificate check failed: ", e$message)
    })
  }

  # ---- 4. Compute Final Consensus A_est (Vectorized approach) ----
  if (verbose) message("Computing final consensus A_est...")
  A_est <- matrix(0.0, nrow = d, ncol = L) # Initialize with zeros

  for (i in seq_len(n)) {
    Ai <- A_list[[i]]
    li <- ncol(Ai)
    if (li > 0) {
        tLabels <- task_labels_list[[i]]
        # Transform Ai: O_i^T * Ai  (d x li)
        Ai_transformed <- crossprod(O_final[,,i], Ai)
        # Add to the correct columns in A_est
        A_est[, tLabels] <- A_est[, tLabels] + Ai_transformed
    }
  }

  # Average by the number of subjects who contributed to each task
  # Use the correctly computed subject_counts_per_task
  valid_tasks_idx <- which(subject_counts_per_task > 0)
  if (length(valid_tasks_idx) > 0) {
      valid_counts <- subject_counts_per_task[valid_tasks_idx]
      # Divide using matrix subsetting and recycling 'each=d'
      A_est[, valid_tasks_idx] <- A_est[, valid_tasks_idx] / rep(valid_counts, each = d)
  }

  # Set columns with zero counts to NA
  invalid_tasks_idx <- which(subject_counts_per_task == 0)
  if (length(invalid_tasks_idx) > 0) {
    A_est[, invalid_tasks_idx] <- NA_real_
  }

  # ---- 5. Prepare Results ----
  res <- list(
    O_mats     = lapply(seq_len(n), function(i_) O_final[,,i_]), # Convert array back to list
    A_est      = A_est,
    iterations = iter,
    converged  = converged,
    final_diff = final_diff
  )

  if (verbose) message("Generalized Procrustes finished.")
  return(res)
}

#' @rdname generalized_procrustes
#' @method generalized_procrustes hyperdesign
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the task/label variable to use for alignment
#' @param preproc Preprocessing function to apply to the data (default: 
#'   \code{center()})
#' @param max_iter Maximum number of iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param tol_type Type of tolerance check (default: "relative")
#' @param verbose Whether to print progress messages (default: FALSE)
#' @param svd_method SVD method to use (default: "irlba")
#' @param svd_opts Options for SVD method (default: empty list)
#' @param ... Additional arguments (currently unused)
#' @examples
#' \donttest{
#' # Create example hyperdesign data
#' library(multidesign)
#' 
#' # Domain 1: 10 features x 5 tasks
#' d1_data <- matrix(rnorm(50), 10, 5)
#' d1_design <- data.frame(task = factor(c("A", "B", "C", "D", "E")))
#' d1 <- multidesign(d1_data, d1_design)
#' 
#' # Domain 2: 10 features x 4 tasks (partial overlap)
#' d2_data <- matrix(rnorm(40), 10, 4) 
#' d2_design <- data.frame(task = factor(c("A", "C", "D", "F")))
#' d2 <- multidesign(d2_data, d2_design)
#' 
#' # Create hyperdesign
#' hd <- hyperdesign(list(domain1 = d1, domain2 = d2))
#' 
#' # Perform alignment
#' result <- generalized_procrustes(hd, task)
#' 
#' # Access results
#' print(result$converged)
#' print(dim(result$A_est))  # 10 features x 6 tasks
#' }
#' @export
generalized_procrustes.hyperdesign <- function(data, y,
                                               preproc = center(),
                                               max_iter = 100,
                                               tol = 1e-6,
                                               tol_type = c("relative", "absolute"),
                                               verbose = FALSE,
                                               svd_method = c("irlba", "base"),
                                               svd_opts = list(),
                                               ...) {
  
  # Input validation
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty hyperdesign object", call. = FALSE)
  }
  
  tol_type <- match.arg(tol_type)
  svd_method <- match.arg(svd_method)
  
  # Use tidy evaluation for the task variable
  y <- rlang::enquo(y)
  
  # Handle both quoted strings and unquoted symbols
  y_char <- tryCatch({
    # If y is a string, use it directly
    if (rlang::quo_is_symbol(y)) {
      rlang::as_name(y)
    } else {
      # If it's a string literal, extract the string
      rlang::eval_tidy(y)
    }
  }, error = function(e) {
    # Fallback: try to convert to name
    rlang::as_name(y)
  })
  
  # Apply preprocessing to get data matrices
  if (verbose) message("Applying preprocessing...")
  pdata <- multivarious::init_transform(data, preproc)
  
  # Validate preprocessed data
  if (any(sapply(pdata, function(x) any(!is.finite(x$x))))) {
    stop("Preprocessed data contains non-finite values (NaN/Inf). ",
         "Check input data quality and preprocessing parameters.", call. = FALSE)
  }
  
  # Extract task labels from each domain AFTER preprocessing to ensure consistency
  task_labels_list <- purrr::map(pdata, function(x) {
    if (!"design" %in% names(x)) {
      stop("Each domain must have a 'design' component", call. = FALSE)
    }
    
    # Validate that the column exists
    if (!y_char %in% names(x$design)) {
      stop("Task column '", y_char, "' not found in domain design frame.", call. = FALSE)
    }
    
    # Extract using base R approach for robustness
    x$design[[y_char]]
  })
  
  # Extract data matrices (transpose to get features x samples for Procrustes)
  A_list <- lapply(pdata, function(x) {
    # Procrustes expects features x tasks, so transpose the samples x features matrix
    t(x$x)
  })
  
  # Get all unique tasks and create mapping
  all_tasks <- unique(unlist(task_labels_list))
  if (length(all_tasks) == 0) {
    stop("No tasks found in the data", call. = FALSE)
  }
  
  # Convert task labels to integers for generalized_procrustes
  task_to_int <- setNames(seq_along(all_tasks), all_tasks)
  L <- length(all_tasks)
  
  # Convert task labels to integer lists
  task_labels_int_list <- lapply(task_labels_list, function(tasks) {
    task_to_int[as.character(tasks)]
  })
  
  if (verbose) {
    message("Data extraction summary:")
    message("  Original domains: ", length(data))
    message("  Preprocessed domains: ", length(pdata))
    message("  Task labels list length: ", length(task_labels_list))
    message("  A_list length: ", length(A_list))
    message("  Task labels int list length: ", length(task_labels_int_list))
    for (i in seq_along(A_list)) {
      message("  Domain ", i, " (", names(A_list)[i], "): ", 
              nrow(A_list[[i]]), " features x ", ncol(A_list[[i]]), " tasks")
      message("    Task labels: ", paste(task_labels_list[[i]], collapse = ", "))
      message("    Task labels int: ", paste(task_labels_int_list[[i]], collapse = ", "))
    }
  }
  
  # Validate dimensions
  d <- nrow(A_list[[1]])
  if (any(sapply(A_list, nrow) != d)) {
    stop("All domains must have the same number of features after preprocessing", call. = FALSE)
  }
  
  if (verbose) {
    message("Running Generalized Procrustes with:")
    message("  Domains: ", length(A_list))
    message("  Features: ", d)
    message("  Total tasks: ", L)
    message("  Tasks per domain: ", paste(sapply(A_list, ncol), collapse = ", "))
  }
  
  # Call the core generalized_procrustes function
  result <- generalized_procrustes(
    A_list = A_list,
    task_labels_list = task_labels_int_list,
    L = L,
    max_iter = max_iter,
    tol = tol,
    tol_type = tol_type,
    verbose = verbose,
    svd_method = svd_method,
    svd_opts = svd_opts
  )
  
  # Add domain names and task labels to result for interpretability
  result$domains <- names(data)
  result$task_labels <- all_tasks
  result$task_mapping <- task_to_int
  
  # Add preprocessing information
  result$preproc <- attr(pdata, "preproc")
  
  if (verbose) message("Generalized Procrustes alignment completed.")
  
  return(result)
}
