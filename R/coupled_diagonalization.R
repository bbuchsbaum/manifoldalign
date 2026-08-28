#' Coupled Diagonalization for Hyperdesign Objects
#'
#' Performs Coupled Diagonalization alignment on hyperdesign data
#' structures containing multiple modalities. Finds coupled eigenvectors
#' that jointly diagonalize graph
#' Laplacians while respecting cross-modal correspondences.
#'
#' @param data A hyperdesign object containing multiple data domains
#' @param correspondence Optional list of correspondence matrices
#'   between modalities. If NULL (default), assumes full correspondence
#'   based on sample order.
#'   Each matrix should be n_i x n_j sparse matrix indicating correspondences 
#'   between modality i and j.
#' @param preproc Preprocessing function to apply to the data (default: center())
#' @param ncomp Number of coupled eigenvectors to extract (default: 10)
#' @param ncomp_per_domain Number of eigenvectors per domain to compute before 
#'   coupling (default: 20). Must be >= ncomp.
#' @param mu_coupling Coupling weight controlling alignment vs diagonalization 
#'   trade-off (default: 1). Higher values emphasize cross-modal alignment.
#' @param alpha_match Optional weight for matching the diagonal of 
#'   \code{A_i^T Lambda_i A_i} to target eigenvalues (default: 0).
#' @param lambda_target Optional list of length \code{m} providing target eigenvalue
#'   vectors for each modality (only the first \code{ncomp} entries are used).
#' @param knn Number of nearest neighbors for graph construction (default: 10)
#' @param sigma Gaussian kernel bandwidth for graph weights. If NULL (default), 
#'   uses median of nearest neighbor distances.
#' @param spectral_cache Optional precomputed spectral cache to reuse
#'   Laplacian subspaces across repeated runs (for example
#'   \code{previous_fit$diagnostics$spectral_cache}). When supplied,
#'   graph/Laplacian/eigendecomposition steps are skipped.
#' @param max_iter Maximum iterations for Stiefel manifold optimization
#'   (default: 200)
#' @param step_size Base step size for gradient descent (default: 0.3)
#' @param tol Convergence tolerance for relative cost change (default: 1e-6)
#' @param verbose Whether to print convergence information (default: FALSE)
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing:
#' \itemize{
#'   \item \code{s}: Coupled eigenvectors (scores) concatenated across domains
#'   \item \code{v}: Projection matrices for out-of-sample extension  
#'   \item \code{coupled_bases}: List of coupled bases V_i for each modality
#'   \item \code{sdev}: Standard deviations of coupled components
#'   \item \code{eigenvalues}: List of eigenvalues for each domain
#'   \item \code{diagnostics}: List with per-iteration cost diagnostics
#'   \item \code{converged}: Logical indicating convergence
#'   \item \code{iterations}: Number of iterations performed
#'   \item \code{final_cost}: Final objective value
#' }
#'
#' @details
#' The algorithm minimizes the objective
#' \preformatted{
#' F = sum_i ||A_i^T Lambda_i A_i - diag(A_i^T Lambda_i A_i)||_F^2 +
#'     mu_c sum_{i<j} ||F_i^T U_i A_i - F_j^T U_j A_j||_F^2
#' }
#' 
#' It uses projected gradient descent on the Stiefel manifold to maintain
#' orthogonality of the coupling matrices A_i.
#' When \code{alpha_match > 0}, the optional match term nudges the diagonal
#' entries of \code{A_i^T Lambda_i A_i} toward user-specified target eigenvalues
#' (without adding extra off-diagonal penalty).
#'
#' @examples
#' \donttest{
#' library(multidesign)
#' # Create synthetic multi-modal data
#' X1 <- matrix(rnorm(150), 50, 3)
#' X2 <- matrix(rnorm(200), 50, 4)
#'
#' # Create hyperdesign
#' hd <- hyperdesign(list(
#'   domain1 = multidesign(X1, data.frame(id = 1:50)),
#'   domain2 = multidesign(X2, data.frame(id = 1:50))
#' ))
#'
#' # Run coupled diagonalization
#' result <- coupled_diagonalization(hd, ncomp = 3)
#' }
#'
#' @references
#' Eynard et al. (2015). Multimodal manifold analysis by simultaneous 
#' diagonalization of Laplacians. IEEE TPAMI, 37(12), 2505-2517.
#'
#' @rdname coupled_diagonalization
#' @method coupled_diagonalization hyperdesign
#' @export
#' @importFrom Matrix sparseMatrix Diagonal crossprod colSums rowSums bdiag
#' @importFrom adjoin graph_weights
#' @importFrom multivarious center init_transform multiblock_biprojector
#' @importFrom chk chk_number chk_true
coupled_diagonalization.hyperdesign <- function(data,
                                              correspondence = NULL,
                                              preproc = center(),
                                              ncomp = 10,
                                              ncomp_per_domain = 20,
                                              mu_coupling = 1,
                                              alpha_match = 0,
                                              lambda_target = NULL,
                                              knn = 10,
                                              sigma = NULL,
                                              spectral_cache = NULL,
                                              max_iter = 200,
                                              step_size = 0.3,
                                              tol = 1e-6,
                                              verbose = FALSE,
                                              ...) {
  
  # Input validation
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(ncomp_per_domain)
  chk::chk_true(ncomp_per_domain >= ncomp)
  chk::chk_number(mu_coupling)
  chk::chk_true(mu_coupling >= 0)
  chk::chk_number(alpha_match)
  chk::chk_true(alpha_match >= 0)
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(step_size)
  chk::chk_true(step_size > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  # Extract number of domains
  m <- length(data)
  if (m < 2) {
    stop("Coupled diagonalization requires at least 2 domains", call. = FALSE)
  }
  user_lambda_target_supplied <- !is.null(lambda_target)
  if (!is.null(lambda_target)) {
    if (!is.list(lambda_target) || length(lambda_target) != m) {
      stop("lambda_target must be a list with one element per domain", call. = FALSE)
    }
    lambda_target <- lapply(lambda_target, function(x) {
      if (is.matrix(x)) {
        x <- diag(x)
      }
      as.numeric(x)
    })
  }
  
  # Apply preprocessing to each domain
  pdata <- data
  nb <- m

  if (is.null(preproc)) {
    proclist <- rep(list(NULL), nb)
    for (i in seq_len(nb)) {
      pdata[[i]]$x <- as.matrix(pdata[[i]]$x)
    }
  } else {
    is_stateful_preproc <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")
    broadcast_preproc <- is.function(preproc) || is_stateful_preproc
    if (broadcast_preproc) {
      if (is_stateful_preproc) {
        preproc <- replicate(
          nb,
          unserialize(serialize(preproc, connection = NULL)),
          simplify = FALSE
        )
      } else {
        preproc <- rep(list(preproc), nb)
      }
    } else if (is.list(preproc) && length(preproc) != nb) {
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
      Xi <- as.matrix(data[[i]]$x)
      pre_i <- preproc[[i]]

      if (is.null(pre_i)) {
        pdata[[i]]$x <- Xi
        proclist[[i]] <- NULL
        next
      }

      # Functions are applied directly (avoid deprecated multivarious::prep()).
      if (is.function(pre_i)) {
        Xi_tf <- pre_i(Xi)
        pdata[[i]]$x <- Xi_tf
        proclist[[i]] <- pre_i
        next
      }

      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }

      if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
        ft <- multivarious::fit_transform(pre_i, Xi)
        pdata[[i]]$x <- ft$transformed
        proclist[[i]] <- ft$preproc
      } else {
        # Back-compat for older multivarious versions
        templ <- multivarious::prep(pre_i)
        Xi_tf <- multivarious::init_transform(templ, Xi)
        pdata[[i]]$x <- Xi_tf
        proc_attr <- attr(Xi_tf, "preproc")
        proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
      }
    }
  }

  names(proclist) <- names(pdata)
  
  # Get sample block indices using the exported function
  sample_block_idx <- block_indices(pdata)
  
  sample_blocks_list <- lapply(seq_len(nrow(sample_block_idx)), function(i) {
    out <- as.integer(sample_block_idx[i, ])
    names(out) <- colnames(sample_block_idx)
    out
  })
  names(sample_blocks_list) <- names(pdata)
  
  # This call should now succeed.
  proc <- multivarious::concat_pre_processors(proclist, sample_blocks_list)
  
  feature_blocks <- feature_block_indices(pdata)
  
  # Extract data matrices
  X_list <- lapply(pdata, function(x) x$x)
  
  # Build correspondence matrices if not provided
  if (is.null(correspondence)) {
    # Assume full correspondence based on sample order
    n_samples <- sapply(X_list, nrow)
    
    if (length(unique(n_samples)) > 1) {
      warning("Domains have different numbers of samples. Using identity ",
              "correspondence for matching samples only.", call. = FALSE)
    }
    
    # Create identity correspondence matrices
    n_min <- min(n_samples)
    correspondence <- lapply(seq_len(m), function(i) {
      n_i <- n_samples[i]
      # Identity matrix truncated to minimum size
      F_i <- Matrix::sparseMatrix(i = seq_len(n_min), 
                                  j = seq_len(n_min), 
                                  x = 1, 
                                  dims = c(n_i, n_min))
      F_i
    })
  }
  
  # Validate correspondence matrices
  if (!is.list(correspondence) || length(correspondence) != m) {
    stop("correspondence must be a list of length ", m, " (one per domain)", 
         call. = FALSE)
  }
  corr_cols <- vapply(correspondence, ncol, integer(1))
  if (length(unique(corr_cols)) != 1) {
    stop("All correspondence matrices must share the same number of columns", call. = FALSE)
  }
  for (i in seq_len(m)) {
    if (nrow(correspondence[[i]]) != nrow(X_list[[i]])) {
      stop("correspondence[[", i, "]] has ", nrow(correspondence[[i]]),
           " rows but domain ", i, " has ", nrow(X_list[[i]]), " samples", call. = FALSE)
    }
  }
  
  used_spectral_cache <- !is.null(spectral_cache)
  eigendecomps <- vector("list", m)
  
  if (!used_spectral_cache) {
    # Build graph Laplacians using adjoin
    if (verbose) message("Building graph Laplacians...")
    
    for (i in seq_len(m)) {
      domain_name <- if (!is.null(names(data)) && length(names(data)) >= i) {
        names(data)[i]
      } else {
        paste0("domain_", i)
      }
      # Use adjoin for consistency with package
      W_graph <- adjoin::graph_weights(X_list[[i]],
                                              k = knn,
                                              weight_mode = "heat",
                                              sigma = sigma)
      W <- adjoin::adjacency(W_graph)
      # Symmetrize weights in case the builder returned a directed graph
      W <- 0.5 * (W + Matrix::t(W))
      Matrix::diag(W) <- 0
      
      # Compute normalized Laplacian
      d <- Matrix::rowSums(W)
      d[d == 0] <- 1
      Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
      L <- Matrix::forceSymmetric(Matrix::Diagonal(nrow(W)) - Dm12 %*% W %*% Dm12)
      
      # Compute eigendecomposition
      # Request one additional eigenpair so we can drop the trivial λ≈0 mode
      n_eigs_total <- min(ncomp_per_domain + 1L, nrow(L))
      eig <- tryCatch({
        if (requireNamespace("RSpectra", quietly = TRUE) && n_eigs_total < nrow(L)) {
          RSpectra::eigs_sym(L, n_eigs_total, which = "SM")
        } else if (n_eigs_total < nrow(L) / 2 && nrow(L) > 50) {
          PRIMME::eigs_sym(L, NEig = n_eigs_total, which = "SA",
                          method = "PRIMME_DEFAULT_MIN_MATVECS")
        } else {
          eig_full <- eigen(as.matrix(L), symmetric = TRUE)
          idx <- order(eig_full$values)[seq_len(n_eigs_total)]
          list(
            vectors = eig_full$vectors[, idx, drop = FALSE],
            values = eig_full$values[idx]
          )
        }
      }, error = function(e) {
        if (verbose) message("Sparse eigensolver failed, using base eigen: ", e$message)
        eig_full <- eigen(as.matrix(L), symmetric = TRUE)
        idx <- order(eig_full$values)[seq_len(n_eigs_total)]
        list(
          vectors = eig_full$vectors[, idx, drop = FALSE],
          values = eig_full$values[idx]
        )
      })
      
      ord <- order(eig$values)
      vals <- eig$values[ord]
      vecs <- eig$vectors[, ord, drop = FALSE]
      non_trivial <- vals > 1e-12
      if (!any(non_trivial)) {
        stop("No non-trivial Laplacian eigenvalues found for domain ", domain_name,
             call. = FALSE)
      }
      keep <- head(which(non_trivial), ncomp_per_domain)
      if (length(keep) < ncomp_per_domain) {
        warning("Domain ", domain_name, ": only ", length(keep),
                " non-trivial eigenvectors available; reducing ncomp_per_domain.",
                call. = FALSE)
      }
      keep <- keep[keep <= ncol(vecs)]
      U_i <- vecs[, keep, drop = FALSE]
      lambda_i <- vals[keep]
      # Normalize spectra so diagonalization terms are comparable across domains
      lambda_ref_idx <- max(1L, min(length(lambda_i), ncomp))
      lambda_ref <- lambda_i[lambda_ref_idx]
      scale_factor <- max(lambda_ref, 1e-8)
      lambda_i <- lambda_i / scale_factor
      attr(lambda_i, "scale_factor") <- scale_factor
      
      eigendecomps[[i]] <- list(
        vectors = U_i,
        values = lambda_i
      )
    }
  } else {
    if (verbose) message("Using cached spectral decomposition.")
    parsed_cache <- parse_spectral_cache_cd(spectral_cache, m)
    Ubar_cache <- parsed_cache$Ubar
    Lambda_cache <- parsed_cache$Lambda
    min_cached <- min(vapply(Ubar_cache, ncol, integer(1)))
    if (min_cached < ncomp) {
      stop("spectral_cache provides only ", min_cached,
           " components per domain, but ncomp = ", ncomp, ".", call. = FALSE)
    }
    if (min_cached < ncomp_per_domain) {
      warning("spectral_cache provides only ", min_cached,
              " components per domain. Reducing ncomp_per_domain.",
              call. = FALSE)
      ncomp_per_domain <- min_cached
    }
    for (i in seq_len(m)) {
      vals_i <- Lambda_cache[[i]][seq_len(ncomp_per_domain)]
      sf_i <- attr(Lambda_cache[[i]], "scale_factor")
      if (!is.null(sf_i) && is.finite(sf_i) && sf_i > 0) {
        attr(vals_i, "scale_factor") <- as.numeric(sf_i)
      }
      eigendecomps[[i]] <- list(
        vectors = Ubar_cache[[i]][, seq_len(ncomp_per_domain), drop = FALSE],
        values = vals_i
      )
    }
  }
  
  # Extract eigenvectors and eigenvalues
  Ubar <- lapply(eigendecomps, function(x) x$vectors)
  Lambda <- lapply(eigendecomps, function(x) x$values)
  
  # Check that we have enough eigenvectors for the requested components
  min_eigenvecs <- min(sapply(eigendecomps, function(x) ncol(x$vectors)))
  if (min_eigenvecs < ncomp) {
    warning("Requested ", ncomp, " components but only ", min_eigenvecs,
            " eigenvectors available. Reducing ncomp.", call. = FALSE)
    ncomp <- min_eigenvecs
  }
  
  scale_factors <- vapply(Lambda, function(vals) {
    sf <- attr(vals, "scale_factor")
    if (is.null(sf) || !is.finite(sf) || sf <= 0) 1 else as.numeric(sf)
  }, numeric(1))
  
  if (alpha_match > 0) {
    if (is.null(lambda_target)) {
      lambda_target <- lapply(seq_len(m), function(i) {
        vals <- Lambda[[i]]
        if (length(vals) == 0) numeric(0) else vals[seq_len(min(ncomp, length(vals)))]
      })
    } else if (user_lambda_target_supplied) {
      # Targets are interpreted in the raw (pre-normalization) eigenvalue scale.
      lambda_target <- Map(function(target_vec, sf) {
        if (is.null(target_vec)) return(NULL)
        as.numeric(target_vec) / max(sf, 1e-8)
      }, lambda_target, scale_factors)
    }
    lambda_target <- lapply(lambda_target, normalize_lambda_target_cd, ncomp = ncomp)
  } else {
    lambda_target <- NULL
  }
  
  # Precompute F_i^T U_i
  FiUbar <- lapply(seq_len(m), function(i) crossprod(correspondence[[i]], Ubar[[i]]))
  coupling_sizes <- vapply(FiUbar, nrow, integer(1))
  coupling_scale <- sqrt(pmax(1, coupling_sizes))
  FiUbar <- Map(function(mat, s) {
    if (is.null(mat) || length(mat) == 0L || s <= 0) return(mat)
    mat / s
  }, FiUbar, coupling_scale)
  
  coupling_gram <- NULL
  if (mu_coupling > 0 && m > 1) {
    coupling_gram <- compute_coupling_gram_cd(FiUbar)
  }
  
  # Initialize coupling matrices A_i on Stiefel manifold
  A <- lapply(seq_len(m), function(i) {
    n_vecs <- ncol(eigendecomps[[i]]$vectors)
    M <- matrix(0, n_vecs, ncomp)
    M[seq_len(ncomp), ] <- diag(ncomp)
    M
  })
  
  # Initialize objective state
  domain_components <- lapply(seq_len(m), function(i) {
    compute_domain_cost_components_cd(
      A[[i]],
      Lambda[[i]],
      if (alpha_match > 0) lambda_target[[i]] else NULL
    )
  })
  diag_costs <- vapply(domain_components, function(x) x$diagonal, numeric(1))
  match_costs <- vapply(domain_components, function(x) x$match, numeric(1))
  diag_sum <- sum(diag_costs)
  match_sum <- sum(match_costs)
  
  gram_sum <- vector("list", m)
  coupling_diag_terms <- numeric(m)
  cross_sum <- 0
  coupling_cost <- 0
  if (!is.null(coupling_gram)) {
    gram_sum <- compute_gram_sum_cd(coupling_gram, A)
    coupling_diag_terms <- vapply(seq_len(m), function(i) {
      sum(A[[i]] * (coupling_gram[[i]][[i]] %*% A[[i]]))
    }, numeric(1))
    cross_sum <- sum(vapply(seq_len(m), function(i) {
      sum(A[[i]] * gram_sum[[i]])
    }, numeric(1)))
    coupling_cost <- m * sum(coupling_diag_terms) - cross_sum
  }
  
  # Optimization on Stiefel manifold with backtracking line search
  converged <- FALSE
  beta_backtrack <- 0.5
  armijo_c <- 1e-4
  max_backtracks <- 12L
  current_step_size <- step_size
  
  cost_current <- diag_sum + alpha_match * match_sum + mu_coupling * coupling_cost
  best_cost <- cost_current
  best_A <- A
  best_components <- list(diagonal = diag_sum, coupling = coupling_cost, match = match_sum)
  
  cost_history <- rep(NA_real_, max_iter + 1L)
  component_history <- vector("list", max_iter + 1L)
  cost_history[1] <- cost_current
  component_history[[1]] <- list(diagonal = diag_sum, coupling = coupling_cost, match = match_sum)
  
  for (iter in seq_len(max_iter)) {
    cost_before_iter <- cost_current
    A_prev <- A
    any_update <- FALSE
    
    for (i in seq_len(m)) {
      grad_i <- compute_gradient_cd(
        idx = i,
        A = A,
        Lambda = Lambda,
        FiUbar = FiUbar,
        mu_coupling = mu_coupling,
        FiUbarA = NULL,
        sumB = NULL,
        lambda_target = lambda_target,
        alpha_match = alpha_match,
        coupling_gram = coupling_gram,
        gram_sum = gram_sum
      )
      if (any(!is.finite(grad_i))) {
        warning("Invalid gradient at iteration ", iter, " for domain ", i,
                ". Skipping update.", call. = FALSE)
        next
      }
      grad_i <- project_stiefel_gradient_cd(A[[i]], grad_i)
      gnorm <- sqrt(sum(grad_i * grad_i))
      if (gnorm < 1e-12) next
      
      direction <- grad_i / gnorm
      step_try <- current_step_size
      accepted <- FALSE
      old_Ai <- A[[i]]
      
      for (ls in seq_len(max_backtracks)) {
        candidate <- old_Ai - step_try * direction
        candidate <- qr.Q(qr(candidate))[, seq_len(ncomp), drop = FALSE]
        overlap <- crossprod(A_prev[[i]], candidate)
        align <- sign(diag(overlap))
        align[align == 0] <- 1
        candidate <- sweep(candidate, 2L, align, `*`)
        
        domain_trial <- compute_domain_cost_components_cd(
          candidate,
          Lambda[[i]],
          if (alpha_match > 0) lambda_target[[i]] else NULL
        )
        diag_sum_trial <- diag_sum - diag_costs[i] + domain_trial$diagonal
        match_sum_trial <- match_sum - match_costs[i] + domain_trial$match
        
        if (!is.null(coupling_gram)) {
          deltaA <- candidate - old_Ai
          gram_sum_trial <- vector("list", m)
          for (j in seq_len(m)) {
            gram_sum_trial[[j]] <- gram_sum[[j]] + coupling_gram[[j]][[i]] %*% deltaA
          }
          
          coupling_diag_terms_trial <- coupling_diag_terms
          coupling_diag_terms_trial[i] <- sum(candidate * (coupling_gram[[i]][[i]] %*% candidate))
          cross_sum_trial <- 0
          for (j in seq_len(m)) {
            Aj_trial <- if (j == i) candidate else A[[j]]
            cross_sum_trial <- cross_sum_trial + sum(Aj_trial * gram_sum_trial[[j]])
          }
          coupling_cost_trial <- m * sum(coupling_diag_terms_trial) - cross_sum_trial
        } else {
          gram_sum_trial <- gram_sum
          coupling_diag_terms_trial <- coupling_diag_terms
          cross_sum_trial <- cross_sum
          coupling_cost_trial <- 0
        }
        
        cost_trial <- diag_sum_trial +
          alpha_match * match_sum_trial +
          mu_coupling * coupling_cost_trial
        
        # Armijo condition for normalized steepest-descent direction.
        if (is.finite(cost_trial) &&
            cost_trial <= cost_current - armijo_c * step_try * gnorm) {
          A[[i]] <- candidate
          diag_costs[i] <- domain_trial$diagonal
          match_costs[i] <- domain_trial$match
          diag_sum <- diag_sum_trial
          match_sum <- match_sum_trial
          gram_sum <- gram_sum_trial
          coupling_diag_terms <- coupling_diag_terms_trial
          cross_sum <- cross_sum_trial
          coupling_cost <- coupling_cost_trial
          cost_current <- cost_trial
          accepted <- TRUE
          any_update <- TRUE
          break
        }
        step_try <- step_try * beta_backtrack
      }
      
      if (!accepted && verbose) {
        message("Backtracking rejected domain ", i, " at iteration ", iter)
      }
    }
    
    if (!any_update) {
      if (verbose) message("No successful updates at iteration ", iter)
      if (iter > 1) {
        converged <- TRUE
        break
      }
    }
    
    if (cost_current < best_cost) {
      best_cost <- cost_current
      best_A <- A
      best_components <- list(diagonal = diag_sum, coupling = coupling_cost, match = match_sum)
    }
    
    cost_history[iter + 1L] <- cost_current
    component_history[[iter + 1L]] <- list(
      diagonal = diag_sum,
      coupling = coupling_cost,
      match = match_sum
    )
    
    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %3d   cost = %.6e", iter, cost_current))
    }
    
    rel_change <- abs(cost_before_iter - cost_current) / (abs(cost_before_iter) + 1e-9)
    if (rel_change < tol) {
      converged <- TRUE
      if (verbose) message("Converged at iteration ", iter)
      break
    }
    
    if (cost_current > cost_before_iter) {
      current_step_size <- max(current_step_size * beta_backtrack, 1e-6)
    } else {
      current_step_size <- min(current_step_size * 1.05, step_size * 5)
    }
  }
  
  if (!converged && verbose) {
    warning("Did not converge within ", max_iter, " iterations", call. = FALSE)
  }
  
  if (best_cost < cost_current) {
    A <- best_A
    cost_current <- best_cost
    diag_sum <- best_components$diagonal
    coupling_cost <- best_components$coupling
    match_sum <- best_components$match
  }
  
  final_cost <- cost_current
  iterations_run <- if (exists("iter")) iter else max_iter
  n_history <- min(iterations_run + 1L, length(cost_history))
  cost_history_valid <- cost_history[seq_len(n_history)]
  component_history_valid <- component_history[seq_len(n_history)]
  component_history_valid <- Filter(function(x) !is.null(x), component_history_valid)
  if (length(cost_history_valid) > 0) {
    names(cost_history_valid) <- paste0("iter_", seq_along(cost_history_valid) - 1L)
  }
  if (length(component_history_valid) > 0) {
    names(component_history_valid) <- paste0("iter_", seq_along(component_history_valid) - 1L)
  }

  # Compute final coupled bases V_i = U_i A_i
  V <- lapply(seq_len(m), function(i) Ubar[[i]] %*% A[[i]])
  names(V) <- names(data)
  
  # Concatenate scores for multiblock structure
  s <- do.call(rbind, V)
  
  # Prepare v matrix for projections using a ridge-regularized blockwise fit
  ridge_penalty <- 1e-3
  v_blocks <- vector("list", m)
  for (i in seq_len(m)) {
    Xi <- X_list[[i]]
    Vi <- V[[i]]
    if (ncol(Xi) == 0) {
      v_blocks[[i]] <- matrix(0, 0, ncol(Vi))
      next
    }
    gram <- as.matrix(crossprod(Xi))
    rhs <- as.matrix(crossprod(Xi, Vi))
    penalty <- ridge_penalty * diag(ncol(Xi))
    v_blocks[[i]] <- solve(gram + penalty, rhs)
  }
  v <- do.call(rbind, v_blocks)
  
  diagnostics <- list(
    cost_history = cost_history_valid,
    component_history = component_history_valid,
    coupling_scale = coupling_scale,
    lambda_scale_factors = scale_factors,
    lambda_target_scaled = lambda_target,
    used_spectral_cache = used_spectral_cache,
    spectral_cache = list(
      subspace_bases = Ubar,
      eigenvalues = Lambda
    ),
    stiefel_factors = A,
    subspace_bases = Ubar
  )

  new_alignment_result(
    scores = s,
    loadings = v,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "coupled_diagonalization",
    extras = list(
      coupled_bases = V,
      eigenvalues = Lambda,
      converged = converged,
      iterations = iterations_run,
      final_cost = final_cost,
      diagnostics = diagnostics
    )
  )
}

parse_spectral_cache_cd <- function(spectral_cache, m) {
  if (!is.list(spectral_cache)) {
    stop("spectral_cache must be a list.", call. = FALSE)
  }
  Ubar <- spectral_cache$subspace_bases
  if (is.null(Ubar)) Ubar <- spectral_cache$Ubar
  if (is.null(Ubar)) Ubar <- spectral_cache$vectors
  Lambda <- spectral_cache$eigenvalues
  if (is.null(Lambda)) Lambda <- spectral_cache$Lambda
  if (is.null(Lambda)) Lambda <- spectral_cache$values
  
  if (!is.list(Ubar) || !is.list(Lambda) || length(Ubar) != m || length(Lambda) != m) {
    stop("spectral_cache must provide list entries `subspace_bases`/`Ubar` and ",
         "`eigenvalues`/`Lambda` with one element per domain.", call. = FALSE)
  }
  
  Ubar_out <- vector("list", m)
  Lambda_out <- vector("list", m)
  for (i in seq_len(m)) {
    Ui <- as.matrix(Ubar[[i]])
    li_raw <- Lambda[[i]]
    if (is.matrix(li_raw)) {
      li_raw <- diag(li_raw)
    }
    li <- as.numeric(li_raw)
    if (ncol(Ui) == 0L || length(li) == 0L) {
      stop("spectral_cache domain ", i, " has empty basis/eigenvalues.", call. = FALSE)
    }
    if (length(li) < ncol(Ui)) {
      stop("spectral_cache domain ", i, " has fewer eigenvalues (", length(li),
           ") than basis vectors (", ncol(Ui), ").", call. = FALSE)
    }
    li <- li[seq_len(ncol(Ui))]
    sf <- attr(Lambda[[i]], "scale_factor")
    if (is.null(sf) || !is.finite(sf) || sf <= 0) sf <- 1
    attr(li, "scale_factor") <- as.numeric(sf)
    Ubar_out[[i]] <- Ui
    Lambda_out[[i]] <- li
  }
  
  list(Ubar = Ubar_out, Lambda = Lambda_out)
}

normalize_lambda_target_cd <- function(target_vec, ncomp) {
  if (is.null(target_vec)) return(NULL)
  if (is.matrix(target_vec)) target_vec <- diag(target_vec)
  target_vec <- as.numeric(target_vec)
  if (length(target_vec) == 0L) return(rep(0, ncomp))
  if (length(target_vec) < ncomp) {
    target_vec <- c(target_vec, rep(tail(target_vec, 1L), ncomp - length(target_vec)))
  }
  target_vec[seq_len(ncomp)]
}

compute_domain_cost_components_cd <- function(Ai, lambda_i, target_vec = NULL) {
  lambda_i <- pmax(lambda_i, 1e-8)
  LiAi <- Ai * lambda_i
  Mi <- crossprod(Ai, LiAi)
  off_diag <- Mi - diag(diag(Mi))
  diagonal_cost <- sum(off_diag * off_diag)
  match_cost <- 0
  if (!is.null(target_vec)) {
    target_vec <- normalize_lambda_target_cd(target_vec, ncol(Ai))
    diag_residual <- diag(Mi) - target_vec
    match_cost <- sum(diag_residual * diag_residual)
  }
  list(diagonal = diagonal_cost, match = match_cost)
}

compute_coupling_gram_cd <- function(FiUbar) {
  m <- length(FiUbar)
  gram <- vector("list", m)
  for (i in seq_len(m)) {
    gram[[i]] <- vector("list", m)
    for (j in seq_len(m)) {
      gram[[i]][[j]] <- crossprod(FiUbar[[i]], FiUbar[[j]])
    }
  }
  gram
}

compute_gram_sum_cd <- function(coupling_gram, A) {
  m <- length(A)
  gram_sum <- vector("list", m)
  for (i in seq_len(m)) {
    acc <- matrix(0, nrow(A[[i]]), ncol(A[[i]]))
    for (j in seq_len(m)) {
      acc <- acc + coupling_gram[[i]][[j]] %*% A[[j]]
    }
    gram_sum[[i]] <- acc
  }
  gram_sum
}

project_stiefel_gradient_cd <- function(Ai, G) {
  AtG <- crossprod(Ai, G)
  G - Ai %*% (0.5 * (AtG + t(AtG)))
}

# Helper function: compute gradient for domain i
compute_gradient_cd <- function(idx,
                                A,
                                Lambda,
                                FiUbar,
                                mu_coupling,
                                FiUbarA = NULL,
                                sumB = NULL,
                                lambda_target = NULL,
                                alpha_match = 0,
                                coupling_gram = NULL,
                                gram_sum = NULL) {
  m <- length(A)
  Ai <- A[[idx]]
  lambda_i <- pmax(Lambda[[idx]], 1e-8)
  LiAi <- Ai * lambda_i
  Mi <- crossprod(Ai, LiAi)
  S <- Mi - diag(diag(Mi))
  
  # Gradient of off-diagonalization term
  grad_diag <- 4 * (LiAi %*% S)
  
  grad_match <- 0
  if (alpha_match > 0 && !is.null(lambda_target) && !is.null(lambda_target[[idx]])) {
    target_vec <- normalize_lambda_target_cd(lambda_target[[idx]], ncol(Ai))
    diag_residual <- diag(Mi) - target_vec
    grad_match <- 4 * sweep(LiAi, 2L, diag_residual, `*`)
  }
  
  grad_coupling <- 0
  if (mu_coupling > 0 && m > 1) {
    if (!is.null(coupling_gram)) {
      sumGA_i <- if (!is.null(gram_sum)) {
        gram_sum[[idx]]
      } else {
        acc <- matrix(0, nrow(Ai), ncol(Ai))
        for (j in seq_len(m)) {
          acc <- acc + coupling_gram[[idx]][[j]] %*% A[[j]]
        }
        acc
      }
      grad_coupling <- 2 * mu_coupling * (
        m * (coupling_gram[[idx]][[idx]] %*% Ai) - sumGA_i
      )
    } else {
      Bi <- if (!is.null(FiUbarA)) FiUbarA[[idx]] else NULL
      if (!is.null(Bi) && !is.null(sumB) && !is.null(FiUbar[[idx]])) {
        grad_coupling <- 2 * mu_coupling * (
          m * crossprod(FiUbar[[idx]], Bi) - crossprod(FiUbar[[idx]], sumB)
        )
      }
    }
  }
  
  grad_diag + grad_coupling + alpha_match * grad_match
}

# Helper function: compute objective cost
compute_cost_cd <- function(A,
                            Lambda,
                            FiUbar,
                            mu_coupling,
                            lambda_target = NULL,
                            alpha_match = 0,
                            coupling_gram = NULL) {
  m <- length(A)
  cost_d <- 0
  cost_match <- 0
  
  for (i in seq_len(m)) {
    Ai <- A[[i]]
    lambda_i <- pmax(Lambda[[i]], 1e-8)
    LiAi <- Ai * lambda_i
    Mi <- crossprod(Ai, LiAi)
    off_diag <- Mi - diag(diag(Mi))
    cost_d <- cost_d + sum(off_diag * off_diag)
    if (alpha_match > 0 && !is.null(lambda_target) && !is.null(lambda_target[[i]])) {
      target_vec <- normalize_lambda_target_cd(lambda_target[[i]], ncol(Ai))
      diag_residual <- diag(Mi) - target_vec
      cost_match <- cost_match + sum(diag_residual * diag_residual)
    }
  }
  
  cost_c <- 0
  if (mu_coupling > 0 && m > 1) {
    if (!is.null(coupling_gram)) {
      gram_sum <- compute_gram_sum_cd(coupling_gram, A)
      diag_terms <- vapply(seq_len(m), function(i) {
        sum(A[[i]] * (coupling_gram[[i]][[i]] %*% A[[i]]))
      }, numeric(1))
      cross_sum <- sum(vapply(seq_len(m), function(i) {
        sum(A[[i]] * gram_sum[[i]])
      }, numeric(1)))
      cost_c <- m * sum(diag_terms) - cross_sum
    } else {
      B <- vector("list", m)
      norms_B <- numeric(m)
      for (i in seq_len(m)) {
        if (!is.null(FiUbar[[i]]) && length(FiUbar[[i]]) > 0L) {
          B[[i]] <- FiUbar[[i]] %*% A[[i]]
          norms_B[i] <- sum(B[[i]] * B[[i]])
        } else {
          B[[i]] <- NULL
        }
      }
      non_null_B <- Filter(function(x) !is.null(x), B)
      if (length(non_null_B) > 0) {
        sumB <- non_null_B[[1]]
        if (length(non_null_B) > 1) {
          for (idx in 2:length(non_null_B)) {
            sumB <- sumB + non_null_B[[idx]]
          }
        }
        cost_c <- m * sum(norms_B) - sum(sumB * sumB)
      }
    }
  }
  
  total <- cost_d + mu_coupling * cost_c + alpha_match * cost_match
  attr(total, "components") <- list(
    diagonal = cost_d,
    coupling = cost_c,
    match = cost_match
  )
  total
}

#' @rdname coupled_diagonalization
#' @method coupled_diagonalization default
#' @export
coupled_diagonalization.default <- function(data, ...) {
  stop("No applicable method for coupled_diagonalization. ",
       "data must be a hyperdesign object.", call. = FALSE)
}

# For backward compatibility - keep the original function as internal
#' @keywords internal
coupled_diagonalisation_internal <- function(
    X_list,                # list of data matrices (modality i: n_i × d_i)
    F_list,                # list of correspondence matrices (n_i × q)
    k        = 10,         # # joint eigenvectors wanted
    kprime   = 20,         # # low-frequency eigenvectors per modality (≥k)
    mu_c     = 1,          # coupling weight
    knn      = 10,         # k-NN for graph
    sigma    = NULL,       # Gaussian σ (NULL => median NN dist)
    maxIter  = 200,        # outer iterations
    step0    = 0.3,        # base step size
    tol      = 1e-6) {     # relative cost tolerance
  
  # This is the original implementation preserved for reference
  # Users should use the new coupled_diagonalization() function instead
  
  m <- length(X_list)
  stopifnot(length(F_list) == m, kprime >= k)
  
  # Build Laplacians and get first k′ eigenpairs
  Ubar    <- vector("list", m)
  Lambda  <- vector("list", m)
  FiUbar  <- vector("list", m)
  
  for (i in seq_len(m)) {
    # Use FNN for backward compatibility in internal function
    nn <- FNN::get.knn(X_list[[i]], k = knn)
    n  <- nrow(X_list[[i]])
    i_idx  <- rep(seq_len(n), each = knn)
    j_idx  <- as.vector(nn$nn.index)
    d2 <- as.vector(nn$nn.dist)^2
    
    sigma_i <- sigma
    if (is.null(sigma_i)) {
      sigma_i <- median(sqrt(d2))
    }
    
    w  <- exp(-d2 / (2 * sigma_i^2))
    W  <- Matrix::sparseMatrix(i = c(i_idx, j_idx), j = c(j_idx, i_idx), 
                              x = c(w, w), dims = c(n, n))
    diag(W) <- 0
    
    # Laplacian
    d <- rowSums(W)
    d[d == 0] <- 1
    Dm12 <- Matrix::Diagonal(x = 1 / sqrt(d))
    L <- Matrix::Diagonal(n = nrow(W)) - Dm12 %*% W %*% Dm12
    
    eig <- RSpectra::eigs_sym(L, kprime, which = "SM")
    Ubar[[i]]   <- eig$vectors
    Lambda[[i]] <- eig$values
    FiUbar[[i]] <- crossprod(F_list[[i]], Ubar[[i]])
  }
  
  # Initialize A_i
  A <- lapply(seq_len(m), function(i) {
    M <- matrix(0, kprime, k)
    M[seq_len(k), ] <- diag(k)
    M
  })
  
  # Cost function
  cost_function <- function(A) {
    cost_d <- 0
    cost_c <- 0
    for (i in seq_len(m)) {
      Ai <- A[[i]]
      Li_diag <- diag(Lambda[[i]])
      Mi <- crossprod(Ai, Li_diag %*% Ai)
      cost_d <- cost_d + sum((Mi - diag(diag(Mi)))^2)
    }
    
    for (i in 1:(m - 1)) {
      for (j in (i + 1):m) {
        diff <- FiUbar[[i]] %*% A[[i]] - FiUbar[[j]] %*% A[[j]]
        cost_c <- cost_c + sum(diff * diff)
      }
    }
    cost_d + mu_c * cost_c
  }
  
  # Gradient
  gradient_block <- function(idx, A) {
    Ai      <- A[[idx]]
    Li_diag <- diag(Lambda[[idx]])
    Mi      <- crossprod(Ai, Li_diag %*% Ai)
    S       <- Mi - diag(diag(Mi))
    
    grad_diag <- 4 * (Li_diag %*% Ai %*% S)
    
    accum <- matrix(0, nrow = kprime, ncol = k)
    for (j in seq_len(m)) {
      if (j != idx) {
        diff   <- FiUbar[[idx]] %*% Ai - FiUbar[[j]] %*% A[[j]]
        accum  <- accum + crossprod(FiUbar[[idx]], diff)
      }
    }
    grad_coupling <- 2 * mu_c * accum
    
    grad_diag + grad_coupling
  }
  
  # Optimization
  cost_old <- Inf
  for (iter in seq_len(maxIter)) {
    for (i in seq_len(m)) {
      grad_i <- gradient_block(i, A)
      gnorm  <- sqrt(sum(grad_i * grad_i))
      eta    <- step0 / max(gnorm, 1e-9)
      A_new  <- A[[i]] - eta * grad_i
      A[[i]] <- qr.Q(qr(A_new))[, seq_len(k)]
    }
    
    cost_new <- cost_function(A)
    
    if ((cost_old - cost_new) / (abs(cost_old) + 1e-9) < tol) {
      message("Convergence criterion met.")
      break
    }
    cost_old <- cost_new
  }
  
  # Final coupled bases
  V <- lapply(seq_len(m), function(i) Ubar[[i]] %*% A[[i]])
  names(V) <- paste0("modality_", seq_len(m))
  invisible(V)
}
