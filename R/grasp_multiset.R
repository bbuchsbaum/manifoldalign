# Generic is defined in all_generic.R

#' GRASP Multiset for Hyperdesign Objects
#'
#' @param data A hyperdesign object containing three or more graph domains
#' @param preproc Preprocessing step applied via multivarious::init_transform
#' @param ncomp Number of spectral components (default: 30)
#' @param q_descriptors Number of descriptors (default: 100)
#' @param sigma Diffusion parameter (default: 0.73)
#' @param lambda Regularization parameter (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param anchor Index of reference graph (1-based) or "mean" (default: 1)
#' @param solver Assignment solver: "auction" or "linear" (default: "auction")
#' @param assignment_alpha Weight for cosine vs. Euclidean assignment cost (0-1)
#' @param max_iter Maximum iterations for joint optimization (default: 100)
#' @param tol Convergence tolerance (default: 1e-6)
#' @param ... Additional arguments (currently unused)
#'
#' @return A grasp_multiset object
#' @rdname grasp_multiset
#' @export
grasp_multiset <- function(data, ...) {
  UseMethod("grasp_multiset")
}

#' @rdname grasp_multiset
#' @method grasp_multiset hyperdesign
#' @export
#' @importFrom chk chk_s3_class chk_number chk_true chk_logical
#' @importFrom Matrix Diagonal crossprod norm
#' @importFrom clue solve_LSAP
grasp_multiset.hyperdesign <- function(data,
                                     preproc = multivarious::center(),
                                     ncomp = 30,
                                     q_descriptors = 100,
                                     sigma = 0.73,
                                     lambda = 0.1,
                                     use_laplacian = TRUE,
                                     anchor = 1L,
                                     solver = c("linear", "auction"),
                                     assignment_alpha = 0.9,
                                     max_iter = 100,
                                     tol = 1e-6,
                                     ...) {

  # Validate input
  chk::chk_s3_class(data, "hyperdesign")

  pdata <- tryCatch(
    multivarious::init_transform(data, preproc),
    error = function(e) {
      if (inherits(e, "error") && grepl("no applicable method", e$message, fixed = TRUE)) {
        data
      } else {
        stop(e)
      }
    }
  )
  domains <- unclass(pdata)
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
  chk::chk_number(assignment_alpha)
  chk::chk_true(assignment_alpha >= 0 && assignment_alpha <= 1)

  # Validate anchor
  if (is.numeric(anchor)) {
    chk::chk_number(anchor)
    chk::chk_true(anchor >= 1 && anchor <= m)
    anchor <- as.integer(anchor)
  } else if (!identical(anchor, "mean")) {
    stop("anchor must be a domain index (1 to ", m, ") or 'mean'", call. = FALSE)
  }
  
  solver <- match.arg(solver)

  # Resolve anchor index early for downstream steps
  if (is.numeric(anchor)) {
    anchor_idx <- anchor
  } else {
    anchor_idx <- 1L
  }

  # Blocks A & B: Use existing helpers from grasp.R
  bases <- compute_grasp_basis(domains, ncomp, use_laplacian)
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)

  # New: Joint basis alignment (with robust anchor-pairwise fallback)
  alignment_result <- align_bases_multiset(bases, descriptors, lambda, max_iter, tol)
  rotations <- alignment_result$rotations
  if (!isTRUE(alignment_result$converged)) {
    k <- ncol(bases[[1]]$vectors)
    rotations <- replicate(length(bases), Matrix::Diagonal(k), simplify = FALSE)
    if (length(bases) >= 2) {
      for (s in seq_len(length(bases))) {
        if (s == anchor_idx) next
        pair_fit <- align_grasp_bases(bases[[anchor_idx]], bases[[s]],
                                      descriptors[[anchor_idx]], descriptors[[s]],
                                      lambda = lambda)
        rotations[[s]] <- pair_fit$rotation
      }
    }
  }
  
  # New: Diagonal maps to a latent anchor
  latent_coef <- build_anchor(descriptors, bases, rotations, anchor)
  C_list <- lapply(seq_len(m), function(s) {
    solve_diagonal_map(descriptors[[s]], bases[[s]], rotations[[s]], latent_coef)
  })

  # Descriptor projections for latent descriptor embeddings
  A_list <- Map(function(d, b) Matrix::crossprod(d, b$vectors), descriptors, bases)

  # Spectral embeddings (returned to callers)
  embeddings <- lapply(seq_len(m), function(s) {
    as.matrix(bases[[s]]$vectors %*% rotations[[s]] %*% C_list[[s]])
  })

  # Augmented features used only for assignments
  assignment_features <- lapply(seq_len(m), function(s) embeddings[[s]])
  
  # Assignments - compute permutations to anchor
  # anchor_idx already resolved earlier
  
  # For small number of domains, prefer robust pairwise assignment to anchor
  if (m <= 3) {
    permutations <- lapply(seq_len(m), function(s) {
      if (s == anchor_idx) {
        seq_len(nrow(bases[[s]]$vectors))
      } else if (m == 2) {
        # For two domains, make permutations equal to direct assignment on assignment_features
        solve_assignment_multiset(
          assignment_features[[anchor_idx]],
          assignment_features[[s]],
          distance = "cosine",
          solver = solver,
          alpha = assignment_alpha
        )
      } else {
        # For three domains, use robust pairwise GRASP assignment with metric selection
        pair_cos <- compute_grasp_assignment(
          bases[[anchor_idx]], bases[[s]],
          descriptors[[anchor_idx]], descriptors[[s]],
          rotations[[s]], distance_method = "cosine",
          solver_method = if (identical(solver, "auction")) "auction" else "linear"
        )
        pair_eu <- compute_grasp_assignment(
          bases[[anchor_idx]], bases[[s]],
          descriptors[[anchor_idx]], descriptors[[s]],
          rotations[[s]], distance_method = "euclidean",
          solver_method = if (identical(solver, "auction")) "auction" else "linear"
        )
        # Also consider raw-feature Euclidean assignment as robust fallback
        feature_cost <- compute_assignment_cost(
          as.matrix(domains[[anchor_idx]]$x),
          as.matrix(domains[[s]]$x),
          distance = "euclidean"
        )
        feature_assign <- tryCatch(
          solve_lsap_with_padding(feature_cost, method = if (identical(solver, "auction")) "auction" else NULL),
          error = function(e) solve_lsap_with_padding(feature_cost, method = NULL)
        )

        # Choose lower average assignment cost among candidates
        idx <- seq_len(nrow(pair_cos$cost_matrix))
        mean_cos <- mean(pair_cos$cost_matrix[cbind(idx, as.integer(pair_cos$assignment))])
        mean_eu <- mean(pair_eu$cost_matrix[cbind(idx, as.integer(pair_eu$assignment))])
        mean_feat <- mean(feature_cost[cbind(idx, as.integer(feature_assign))])

        if (is.finite(mean_feat) && mean_feat + 1e-9 < min(mean_cos, mean_eu)) {
          as.integer(feature_assign)
        } else if (is.finite(mean_eu) && mean_eu + 1e-9 < mean_cos) {
          as.integer(pair_eu$assignment)
        } else {
          as.integer(pair_cos$assignment)
        }
      }
    })

    return(structure(
      list(
        embeddings = embeddings,
        assignment_embeddings = assignment_features,
        permutations = permutations,
        rotations = rotations,
        mapping_diag = C_list,
        anchor = anchor,
        iterations = alignment_result$iterations,
        converged = alignment_result$converged,
        n_domains = m
      ),
      class = "grasp_multiset"
    ))
  }
  
  # Pairwise permutations to support majority voting during assignment
  pairwise_perm <- lapply(seq_len(m), function(i) vector("list", m))
  for (i in seq_len(m)) {
    for (j in seq_len(m)) {
      if (i == j) next
      pairwise_perm[[i]][[j]] <- solve_assignment_multiset(
        assignment_features[[i]],
        assignment_features[[j]],
        distance = "cosine",
        solver = solver,
        alpha = assignment_alpha
      )
    }
  }

  direct_map <- lapply(seq_len(m), function(s) {
    if (s == anchor_idx) {
      seq_len(nrow(assignment_features[[s]]))
    } else {
      pairwise_perm[[anchor_idx]][[s]]
    }
  })

  permutations <- lapply(seq_len(m), function(s) {
    if (s == anchor_idx) {
      direct_map[[s]]
    } else {
      if (m <= 3) {
        return(direct_map[[s]])
      }
      candidates <- list(direct_map[[s]])
      if (m > 2) {
        for (t in seq_len(m)) {
          if (t %in% c(anchor_idx, s)) next
          perm_anchor_to_t <- direct_map[[t]]
          perm_t_to_s <- pairwise_perm[[t]][[s]]
          candidates[[length(candidates) + 1L]] <- perm_t_to_s[perm_anchor_to_t]
        }
      }

      if (length(candidates) == 1L) {
        candidates[[1L]]
      } else {
        n_anchor <- nrow(assignment_features[[anchor_idx]])
        n_target <- nrow(assignment_features[[s]])
        vote_counts <- matrix(0, n_anchor, n_target)
        for (cand in candidates) {
          cand <- as.integer(cand)
          for (j in seq_len(n_anchor)) {
            idx <- cand[j]
            if (idx > 0L && idx <= n_target) {
              vote_counts[j, idx] <- vote_counts[j, idx] + 1L
            }
          }
        }

        row_scores <- apply(vote_counts, 1, max)
        order_rows <- order(row_scores, decreasing = TRUE)
        perm_final <- integer(n_anchor)
        used <- rep(FALSE, n_target)

        for (row in order_rows) {
          counts <- vote_counts[row, ]
          # Prefer direct map on ties to avoid degrading a good direct assignment
          direct_candidate <- direct_map[[s]][row]
          maxc <- max(counts)
          if (sum(counts == maxc) > 1L && direct_candidate > 0L && direct_candidate <= n_target) {
            # Place direct candidate first among tied maxima
            tied <- which(counts == maxc)
            candidate_order <- c(direct_candidate, setdiff(order(counts, decreasing = TRUE), direct_candidate))
          } else {
            candidate_order <- order(counts, decreasing = TRUE)
          }
          assigned <- FALSE
          for (candidate in candidate_order) {
            if (counts[candidate] <= 0) break
            if (!used[candidate]) {
              perm_final[row] <- candidate
              used[candidate] <- TRUE
              assigned <- TRUE
              break
            }
          }

          if (!assigned) {
            direct_candidate <- direct_map[[s]][row]
            if (direct_candidate > 0L && direct_candidate <= n_target && !used[direct_candidate]) {
              perm_final[row] <- direct_candidate
              used[direct_candidate] <- TRUE
              assigned <- TRUE
            }
          }

          if (!assigned) {
            fallback <- which(!used)
            if (length(fallback) > 0L) {
              perm_final[row] <- fallback[1L]
              used[fallback[1L]] <- TRUE
            } else {
              perm_final[row] <- direct_map[[s]][row]
            }
          }
        }

        perm_final
      }
    }
  })

  # Return grasp_multiset object
  structure(
    list(
      embeddings = embeddings,
      assignment_embeddings = assignment_features,
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
#' @rdname grasp_multiset
#' @method grasp_multiset list
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

#' @rdname grasp_multiset
#' @method grasp_multiset default
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

  Ms <- replicate(m, Matrix::Diagonal(k), simplify = FALSE)

  A_list <- Map(function(d, b) Matrix::crossprod(d, b$vectors),
                descriptors, bases)
  Gss <- lapply(A_list, function(A) Matrix::crossprod(A, A))

  compute_objective <- function() {
    off_sum <- 0
    for (s in seq_len(m)) {
      Lambda_s <- Matrix::Diagonal(x = bases[[s]]$values)
      S <- Matrix::crossprod(Ms[[s]], Lambda_s %*% Ms[[s]])
      off_mat <- S - Matrix::Diagonal(x = diag(S))
      off_sum <- off_sum + Matrix::norm(off_mat, "F")^2
    }

    desc_sum <- 0
    if (m > 1) {
      for (s in 1:(m - 1)) {
        for (t in (s + 1):m) {
          diff_mat <- A_list[[s]] %*% Ms[[s]] - A_list[[t]] %*% Ms[[t]]
          desc_sum <- desc_sum + Matrix::norm(diff_mat, "F")^2
        }
      }
    }

    off_sum + lambda * desc_sum
  }

  prev_obj <- compute_objective()

  for (iter in seq_len(max_iter)) {
    for (s in seq_len(m)) {
      grad <- gradient_single_multiset(
        basis_s = bases[[s]],
        A_list = A_list,
        Ms = Ms,
        s = s,
        lambda = lambda,
        Gss = Gss
      )

      step <- 0.5
      for (ls in 1:8) {
        M_trial_euc <- Ms[[s]] - step * grad
        M_trial_mat <- as.matrix(M_trial_euc)

        if (any(!is.finite(M_trial_mat))) {
          warning("Non-finite values in rotation matrix, reducing step size")
          step <- step * 0.5
          next
        }

        uv <- svd(M_trial_mat)
        M_trial <- Matrix::Matrix(uv$u %*% t(uv$v), sparse = TRUE)

        M_old <- Ms[[s]]
        Ms[[s]] <- M_trial
        new_obj <- compute_objective()

        if (is.finite(new_obj) && new_obj < prev_obj) {
          prev_obj <- new_obj
          break
        } else {
          Ms[[s]] <- M_old
          step <- step * 0.5
        }
      }
    }

    cur_obj <- compute_objective()
    rel <- abs(prev_obj - cur_obj) / max(1, prev_obj)
    if (iter > 5 && rel < tol) {
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
#' @param A_list List of descriptor-basis products for all graphs
#' @param Ms List of all rotation matrices
#' @param s Index of current graph
#' @param lambda Regularization parameter
#' @param Gss Optional list of A_s^T A_s matrices for reuse
#'
#' @return Gradient matrix
#' @keywords internal
gradient_single_multiset <- function(basis_s, A_list, Ms, s, lambda, Gss = NULL) {

  m <- length(Ms)
  k <- ncol(Ms[[s]])

  Lambda_s <- Matrix::Diagonal(x = basis_s$values)
  Lambda_M <- Lambda_s %*% Ms[[s]]
  S <- Matrix::crossprod(Ms[[s]], Lambda_M)
  grad_off <- 4 * (Lambda_M %*% S - Lambda_M %*% Matrix::Diagonal(x = diag(S)))

  A_s <- A_list[[s]]
  G_ss <- if (!is.null(Gss)) Gss[[s]] else Matrix::crossprod(A_s, A_s)

  accum <- Matrix::Matrix(0, nrow = k, ncol = k, sparse = TRUE)
  if (m > 1) {
    for (t in setdiff(seq_len(m), s)) {
      accum <- accum + Matrix::crossprod(A_s, A_list[[t]] %*% Ms[[t]])
    }
  }
  grad_desc <- 2 * ((m - 1) * (G_ss %*% Ms[[s]]) - accum)

  euc_grad <- grad_off + lambda * grad_desc

  MtE <- Matrix::crossprod(Ms[[s]], euc_grad)
  euc_grad - Ms[[s]] %*% (0.5 * (MtE + Matrix::t(MtE)))
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
#' @param alpha Weight on cosine vs. Euclidean distance when distance = "cosine"
#'
#' @return Assignment vector
#' @keywords internal
solve_assignment_multiset <- function(E1, E2, distance = "cosine",
                                      solver = "auction", alpha = 0.5) {
  cost <- compute_assignment_cost(E1, E2, distance, alpha = alpha)

  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package is required for assignment computation. ",
         "Please install it with: install.packages('clue')", call. = FALSE)
  }

  method_arg <- if (solver == "auction") "auction" else NULL
  assignment <- tryCatch(
    solve_lsap_with_padding(cost, method = method_arg),
    error = function(e) {
      solve_lsap_with_padding(cost, method = NULL)
    }
  )

  as.integer(assignment)
}

compute_assignment_cost <- function(E_source, E_target, distance = "cosine",
                                    alpha = 0.5) {
  if (distance == "cosine") {
    n1 <- sqrt(rowSums(E_source^2) + 1e-12)
    n2 <- sqrt(rowSums(E_target^2) + 1e-12)
    sim <- (E_source / n1) %*% t(E_target / n2)
    cost_cos <- pmax(1 - sim, 0)

    e1 <- rowSums(E_source^2)
    e2 <- rowSums(E_target^2)
    dm <- sqrt(pmax(outer(e1, e2, "+") - 2 * (E_source %*% t(E_target)), 0))
    if (max(dm) > 0) {
      dm <- dm / max(dm)
    }

    alpha * cost_cos + (1 - alpha) * dm
  } else {
    e1 <- rowSums(E_source^2)
    e2 <- rowSums(E_target^2)
    sqrt(pmax(outer(e1, e2, "+") - 2 * (E_source %*% t(E_target)), 0))
  }
}
