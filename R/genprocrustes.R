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
#' \strong{Optimality Certificate}:
#' The iteration maximises \eqn{tr(O^T C O)} with \eqn{C = D D^T}. At the
#' returned rotations, \eqn{\Lambda_i = (C O)_i O_i^T} must be symmetric
#' (first-order stationarity), and \eqn{\Lambda - C \succeq 0} with
#' \eqn{\Lambda = blkdiag(\Lambda_i)} certifies global optimality (Ling, 2024).
#' The smallest eigenvalue of \eqn{\Lambda - C} is computed densely when
#' \eqn{n d \le 1000}; larger problems report stationarity only.
#'
#' \strong{Translations}:
#' The model is rotation-only: \eqn{A_i \approx O_i A_{est}[, tasks_i]}. Inputs
#' must therefore share a common origin. In particular, centring each subject
#' over a different subset of tasks introduces subject-specific offsets that
#' the rotations cannot absorb and biases the consensus; centre all subjects
#' on a common set of tasks instead (the hyperdesign method does this by
#' default, see \code{translation}).
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
#' \item{final_diff}{The Frobenius norm of the change in the rotations at the
#'   last iteration (`sqrt(diff_sq)`), i.e. the value used for the final
#'   convergence check.}
#' \item{certificate}{A list with \code{status} (\code{"computed"},
#'   \code{"stationarity_only"}, \code{"skipped"} or \code{"failed"}),
#'   \code{stationarity_residual} (relative asymmetry of the
#'   \eqn{\Lambda_i}), \code{min_eigenvalue} and
#'   \code{relative_min_eigenvalue} of \eqn{\Lambda - C} (relative to
#'   \eqn{\lambda_{max}(C)}), \code{tol}, and \code{certified} (\code{TRUE}
#'   when both quantities are within \code{tol}; \code{NA} when the
#'   eigenvalue test was not run).}
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
#' @importFrom multivarious center
#' @importFrom irlba irlba
#' @keywords internal
.generalized_procrustes_impl <- function(
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
  .stack_rot <- function(Oarr) {
    d <- dim(Oarr)[1]; n <- dim(Oarr)[3]
    M <- aperm(Oarr, c(1, 3, 2))  # (d, n, d)
    dim(M) <- c(d * n, d)
    M
  }

  .unstack <- function(M, d_, n_) {
    A <- array(M, c(d_, n_, d_))
    aperm(A, c(1, 3, 2))
  }

  # Efficient projection onto O(d)
  .projO <- function(A) {
    d_ <- ncol(A)
    if (d_ == 0) return(A)
    if (any(!is.finite(A))) {
      warning("Non-finite values encountered in matrix projection; returning identity.", call. = FALSE)
      return(diag(d_))
    }
    if (sum(A * A) < 1e-14) return(diag(d_))
    sv <- svd(A, nu = d_, nv = d_)
    sv$u %*% t(sv$v)
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
              sres <- suppressWarnings(do.call(irlba::irlba, args))
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

  if (Matrix::nnzero(D) == 0) {
    warning("Constructed sparse matrix D has no non-zero entries; returning identity rotations.", call. = FALSE)
    O_final <- array(0.0, dim = c(d, d, n))
    for (i in seq_len(n)) O_final[,,i] <- diag(d)
    return(list(
      O_mats     = lapply(seq_len(n), function(i_) O_final[,,i_]),
      A_est      = matrix(NA_real_, d, L),
      iterations = 0,
      converged  = TRUE,
      final_diff = 0,
      certificate = list(status = "skipped", message = "D has no non-zero entries",
                         stationarity_residual = NA_real_, min_eigenvalue = NA_real_,
                         relative_min_eigenvalue = NA_real_, certified = NA,
                         tol = NA_real_)
    ))
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

  # ---- Optimality certificate (Lambda - C test, Ling 2024) ----
  certificate <- tryCatch(
    .gpa_certificate(D, O_final),
    error = function(e) {
      list(status = "failed", message = conditionMessage(e),
           stationarity_residual = NA_real_, min_eigenvalue = NA_real_,
           relative_min_eigenvalue = NA_real_, certified = NA, tol = NA_real_)
    }
  )
  if (verbose) {
    message(sprintf(
      "Optimality certificate (%s): stationarity residual = %.2e, relative min eigenvalue of Lambda - C = %.2e, certified = %s",
      certificate$status, certificate$stationarity_residual,
      certificate$relative_min_eigenvalue, format(certificate$certified)
    ))
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
  O_mats_list <- lapply(seq_len(n), function(i_) {
    matrix(O_final[, , i_], nrow = d, ncol = d)
  })

  res <- list(
    O_mats      = O_mats_list,
    A_est       = A_est,
    iterations  = iter,
    converged   = converged,
    final_diff  = final_diff,
    certificate = certificate
  )
  if (verbose) message("Generalized Procrustes finished.")
  return(res)
}

#' @rdname generalized_procrustes
#' @method generalized_procrustes hyperdesign
#' @param data A hyperdesign object containing multiple data domains
#' @param y Name of the task/label variable to use for alignment
#' @param preproc Preprocessing applied per domain (default:
#'   \code{center()}): a single preprocessor broadcast to all domains, a list
#'   with one entry per domain, or \code{NULL} for none. The fitted
#'   preprocessors are returned in \code{preproc} (one entry per domain,
#'   \code{NULL} where none was applied).
#' @param max_iter Maximum number of iterations (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param tol_type Type of tolerance check (default: "relative")
#' @param verbose Whether to print progress messages (default: FALSE)
#' @param svd_method SVD method to use (default: "irlba")
#' @param svd_opts Options for SVD method (default: empty list)
#' @param translation How per-domain translations are handled after
#'   preprocessing. \code{"shared_tasks"} (default) subtracts from every domain
#'   the mean of its rows on the tasks observed by \emph{all} domains, so that
#'   domains observing different task subsets share a common origin (per-domain
#'   centring over different subsets would otherwise bias the consensus). With
#'   complete designs and \code{center()} preprocessing this is a no-op. If no
#'   task is shared by all domains, a warning is issued and no translation is
#'   applied. \code{"none"} uses the preprocessed data as is (rotation-only
#'   model). The subtracted vectors are returned in \code{translations}.
#' @param ... Additional arguments (currently unused)
#' @examples
#' \donttest{
#' # Create example hyperdesign data
#' library(multidesign)
#' 
#' # Domain 1: 5 tasks x 10 features
#' d1_data <- matrix(rnorm(50), 5, 10)
#' d1_design <- data.frame(task = factor(c("A", "B", "C", "D", "E")))
#' d1 <- multidesign(d1_data, d1_design)
#' 
#' # Domain 2: 4 tasks x 10 features (partial overlap)
#' d2_data <- matrix(rnorm(40), 4, 10) 
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
#' print(dim(result$A_est))  # 10 features x 6 total tasks
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
                                               translation = c("shared_tasks", "none"),
                                               ...) {
  
  # Input validation
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty hyperdesign object", call. = FALSE)
  }
  
  tol_type <- match.arg(tol_type)
  svd_method <- match.arg(svd_method)
  translation <- match.arg(translation)
  
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
  pdata <- unclass(data)
  nb <- length(pdata)

  if (is.null(preproc)) {
    # `proclist[[i]] <- NULL` would delete elements; keep one NULL per domain.
    proclist <- vector("list", nb)
    for (i in seq_len(nb)) {
      pdata[[i]]$x <- as.matrix(pdata[[i]]$x)
    }
  } else {
    is_preproc_obj <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")

    if (!is.list(preproc) || is_preproc_obj) {
      preproc <- rep(list(preproc), nb)
    } else if (length(preproc) == 1L) {
      preproc <- rep(preproc, nb)
    } else if (length(preproc) != nb) {
      stop(
        "Length of `preproc` list (",
        length(preproc),
        ") must match number of domains (",
        nb,
        ").",
        call. = FALSE
      )
    }

    proclist <- vector("list", nb)
    for (i in seq_len(nb)) {
      Xi <- as.matrix(pdata[[i]]$x)
      pre_i <- preproc[[i]]

      if (is.null(pre_i)) {
        pdata[[i]]$x <- Xi
        proclist[i] <- list(NULL)
        next
      }

      # Functions are applied directly (avoid recipes-based preprocessing).
      if (is.function(pre_i)) {
        Xi_tf <- pre_i(Xi)
        pdata[[i]]$x <- Xi_tf
        proclist[[i]] <- pre_i
        next
      }

      # Clone stateful preprocessors before fitting.
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }

      if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
        ft <- multivarious::fit_transform(pre_i, Xi)
        pdata[[i]]$x <- ft$transformed
        proclist[[i]] <- ft$preproc
      } else {
        templ <- multivarious::prep(pre_i)
        Xi_tf <- multivarious::init_transform(templ, Xi)
        pdata[[i]]$x <- Xi_tf
        proc_attr <- attr(Xi_tf, "preproc")
        proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
      }
    }
  }

  names(pdata) <- names(data)
  class(pdata) <- class(data)
  names(proclist) <- names(pdata)
  attr(pdata, "preproc") <- proclist
  
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
  
  # Translation step: centre every domain on the tasks shared by all domains.
  translations <- NULL
  if (translation == "shared_tasks") {
    centred <- .gpa_center_on_shared_tasks(pdata, task_labels_list)
    if (is.null(centred$translations)) {
      warning("No task is observed by every domain, so per-domain translations ",
              "cannot be estimated from shared tasks; the data are used as ",
              "preprocessed. If domains were centred over different task ",
              "subsets the consensus may be biased.", call. = FALSE)
    } else {
      pdata <- centred$pdata
      translations <- centred$translations
    }
  }

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
  result <- .generalized_procrustes_impl(
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
  result$preproc <- proclist
  result$translations <- translations
  
  if (verbose) message("Generalized Procrustes alignment completed.")
  
  return(result)
}

#' @export
generalized_procrustes.default <- function(data = NULL,
                                           task_labels_list = NULL,
                                           L = NULL,
                                           A_list = NULL,
                                           ...) {
  if (is.null(data)) {
    if (!is.null(A_list)) {
      data <- A_list
    } else {
      stop("generalized_procrustes() requires either a hyperdesign object or a list of matrices.",
           " Supply `data` or `A_list`.", call. = FALSE)
    }
  }

  if (!(is.list(data) && all(vapply(data, is.matrix, logical(1L))))) {
    stop("generalized_procrustes() requires either a hyperdesign object or a list of matrices. ",
         "Got: ", paste(class(data), collapse = ", "),
         ". See ?generalized_procrustes for usage examples.", call. = FALSE)
  }

  n <- length(data)

  if (!is.null(task_labels_list)) {
    if (!is.list(task_labels_list) || length(task_labels_list) != n) {
      stop("task_labels_list must be a list matching the length of data", call. = FALSE)
    }
    if (is.null(L)) {
      all_labels <- unlist(task_labels_list, use.names = FALSE)
      if (length(all_labels) == 0) {
        stop("Cannot infer L from empty task_labels_list. Please provide L explicitly.", call. = FALSE)
      }
      L <- max(all_labels)
    }
    .generalized_procrustes_impl(
      A_list = data,
      task_labels_list = task_labels_list,
      L = L,
      ...
    )
  } else {
    Ls <- sapply(data, ncol)
    if (length(unique(Ls)) != 1L) {
      stop("For a raw list of matrices, all subjects must have the same number of columns (tasks).", call. = FALSE)
    }
    if (is.null(L)) {
      L <- Ls[1L]
    }
    task_labels_list <- replicate(n, seq_len(L), simplify = FALSE)
    .generalized_procrustes_impl(
      A_list = data,
      task_labels_list = task_labels_list,
      L = L,
      ...
    )
  }
}

# Optimality certificate for generalized orthogonal Procrustes.
#
# The GPM maximises f(O) = tr(O^T C O) over block-orthogonal O, with
# C = D D^T and O the (n d x d) stack of the rotations. At a stationary point
# (C O)_i = Lambda_i O_i with Lambda_i = (C O)_i O_i^T symmetric. If, in
# addition, Lambda - C is positive semidefinite (Lambda = blkdiag(Lambda_i)),
# O is a global maximiser (dual certificate of the SDP relaxation; Ling 2024).
#
# Returns a list with `status`, `stationarity_residual`, `min_eigenvalue`,
# `relative_min_eigenvalue` (divided by lambda_max(C)), `tol` and `certified`.
#' @noRd
.gpa_certificate <- function(D, O_arr, tol = 1e-6, max_dense = 1000L) {
  d <- dim(O_arr)[1]
  n <- dim(O_arr)[3]
  nd <- n * d
  O_stack <- aperm(O_arr, c(1, 3, 2))
  dim(O_stack) <- c(nd, d)

  CO <- as.matrix(D %*% Matrix::crossprod(D, O_stack))
  co_norm <- max(sqrt(sum(CO * CO)), .Machine$double.eps)

  Lambda <- vector("list", n)
  asym2 <- 0
  for (i in seq_len(n)) {
    rows <- ((i - 1L) * d + 1L):(i * d)
    Li <- CO[rows, , drop = FALSE] %*% t(matrix(O_arr[, , i], d, d))
    asym2 <- asym2 + sum((Li - t(Li))^2)
    Lambda[[i]] <- (Li + t(Li)) / 2
  }
  stationarity <- sqrt(asym2) / co_norm

  # lambda_max(C) is the squared spectral norm of D; use the smaller Gram matrix.
  G <- if (nrow(D) <= ncol(D)) Matrix::tcrossprod(D) else Matrix::crossprod(D)
  lam_max_C <- max(eigen(as.matrix(G), symmetric = TRUE, only.values = TRUE)$values)
  lam_max_C <- max(lam_max_C, .Machine$double.eps)

  if (nd > max_dense) {
    return(list(
      status = "stationarity_only",
      stationarity_residual = stationarity,
      min_eigenvalue = NA_real_,
      relative_min_eigenvalue = NA_real_,
      tol = tol,
      certified = NA
    ))
  }

  S <- as.matrix(Matrix::bdiag(Lambda)) - as.matrix(Matrix::tcrossprod(D))
  S <- (S + t(S)) / 2
  min_eig <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  rel_min <- min_eig / lam_max_C

  list(
    status = "computed",
    stationarity_residual = stationarity,
    min_eigenvalue = min_eig,
    relative_min_eigenvalue = rel_min,
    tol = tol,
    certified = isTRUE(stationarity <= tol && rel_min >= -tol)
  )
}

# Subtract from each domain the mean of its rows on the tasks observed by every
# domain. Under x_ij = R_i a_j + t_i this removes t_i + R_i mean_S(a) with the
# same task set S for every domain, so the rotation-only Procrustes model
# applies to the result. Returns the shifted data and the subtracted vectors
# (NULL translations when no task is shared by all domains).
#' @noRd
.gpa_center_on_shared_tasks <- function(pdata, task_labels_list) {
  task_sets <- lapply(task_labels_list, function(t) unique(as.character(t)))
  shared <- Reduce(intersect, task_sets)
  if (!length(shared)) {
    return(list(pdata = pdata, translations = NULL, shared_tasks = character(0)))
  }
  translations <- vector("list", length(pdata))
  for (i in seq_along(pdata)) {
    X <- as.matrix(pdata[[i]]$x)
    rows <- as.character(task_labels_list[[i]]) %in% shared
    shift <- colMeans(X[rows, , drop = FALSE])
    pdata[[i]]$x <- sweep(X, 2, shift, `-`)
    translations[[i]] <- shift
  }
  names(translations) <- names(pdata)
  list(pdata = pdata, translations = translations, shared_tasks = shared)
}
