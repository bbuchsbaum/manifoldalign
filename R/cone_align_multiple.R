# Generic is defined in all_generic.R

#' CONE-Align Multiple for Hyperdesign Objects
#'
#' @param data A hyperdesign object containing three or more graph domains
#' @param ref_idx Which domain acts as initial reference (default: 1)
#' @inheritParams cone_align.hyperdesign
#'
#' @rdname cone_align_multiple
#' @method cone_align_multiple hyperdesign
#' @export
#' @importFrom chk chk_number chk_true chk_logical
#' @importFrom multivarious center concat_pre_processors
cone_align_multiple.hyperdesign <- function(data,
                                          ref_idx = 1L,
                                          preproc = center(),
                                          ncomp = 10,
                                          sigma = 0.73,
                                          lambda = 0.1,
                                          use_laplacian = TRUE,
                                          solver = c("linear", "auction"),
                                          max_iter = 30,
                                          tol = 0.01,
                                          knn = NULL,
                                          ...) {
  
  # Validate input
  if (!inherits(data, "hyperdesign")) {
    stop("data must be a 'hyperdesign' object", call. = FALSE)
  }
  
  domains <- unclass(data)
  if (length(domains) < 3) {
    stop("cone_align_multiple() requires three or more domains; ",
         "use cone_align() for the pairwise case.", call. = FALSE)
  }
  
  # Validate ref_idx
  chk::chk_number(ref_idx)
  chk::chk_true(ref_idx >= 1 && ref_idx <= length(domains))
  
  # Validate other parameters (same as pairwise)
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
  
  solver <- match.arg(solver)
  
  # Preprocess each domain (identical to pairwise code)
  pdata <- lapply(domains, function(domain) {
    x_preprocessed <- safe_compute(
      if (is.function(preproc)) preproc(domain$x) else domain$x,
      "Data preprocessing failed. Check preprocessing function or data format."
    )
    
    # Add data scaling warning
    if (is.function(preproc) && !inherits(preproc, "scaler")) {
      warning("Consider scaling data before embedding as sigma parameter is ",
              "sensitive to data units", call. = FALSE)
    }
    
    list(x = x_preprocessed, design = domain$design)
  })
  names(pdata) <- names(domains)
  
  # Block indices computation
  block_indices <- compute_block_indices(pdata)
  
  # Create preprocessing structure
  proc <- if (is.function(preproc)) {
    proclist <- setNames(rep(list(preproc), length(pdata)), names(pdata))
    block_indices_list <- split(block_indices, row(block_indices))
    
    safe_compute(
      multivarious::concat_pre_processors(proclist, block_indices_list),
      "Preprocessing concatenation failed"
    )
  } else {
    NULL
  }
  
  # Call multi-graph fitting function
  cone_align_multiple_fit(pdata, proc, ref_idx, ncomp, sigma, lambda,
                         use_laplacian, solver, max_iter, tol, knn,
                         feature_block_indices(pdata))
}

#' CONE-Align Multiple for List of Matrices
#'
#' @param data A list containing three or more matrices (nodes x features)
#' @param ... Additional arguments passed to hyperdesign method
#'
#' @rdname cone_align_multiple
#' @method cone_align_multiple list
#' @export
cone_align_multiple.list <- function(data, ...) {
  
  # Validate input
  if (!is.list(data) || length(data) < 3) {
    stop("Provide a list with at least three matrices.", call. = FALSE)
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
  
  # Convert to stratum format directly (following cone_align.list pattern)
  # This avoids the multidesign API issues
  strata <- lapply(seq_along(data), function(i) {
    list(
      x = data[[i]],
      design = data.frame(node_id = seq_len(nrow(data[[i]])))
    )
  })
  names(strata) <- paste0("domain_", seq_along(strata))
  
  # Wrap in a simple class that mimics hyperdesign structure
  class(strata) <- "hyperdesign"
  
  # Call hyperdesign method
  cone_align_multiple.hyperdesign(strata, ...)
}

#' @rdname cone_align_multiple
#' @method cone_align_multiple default
#' @export
cone_align_multiple.default <- function(data, ...) {
  stop("No applicable method for cone_align_multiple. data must be a hyperdesign ",
       "object or list of matrices with 3+ elements.", call. = FALSE)
}

#' Core Multi-Graph CONE-Align Fitting
#'
#' Internal function that implements the iterative multi-graph
#' alignment algorithm.
#' Uses the helper functions from the pairwise implementation.
#'
#' @param strata List of preprocessed data domains
#' @param proc Preprocessing object
#' @param ref_idx Reference domain index
#' @param ncomp Number of components
#' @param sigma Diffusion parameter
#' @param lambda Regularization parameter
#' @param use_laplacian Whether to use Laplacian normalization
#' @param solver Assignment solver
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param knn Number of nearest neighbors
#' @param feature_blocks Feature block index structure
#'
#' @return multiblock_biprojector object with multi-graph alignment results
#' @keywords internal

# --------------------------------------------------------------------
# helper: pairwise alignment that returns P/Q ready for multi-graph use
# --------------------------------------------------------------------
.align_two_embeddings <- function(Eref, Eg, solver, max_iter, tol, lambda) {
  pair_res <- cone_align_iterate(
                list(Eref, Eg),      # what cone_align_iterate expects
                solver   = solver,
                max_iter = max_iter,
                tol      = tol,
                lambda   = lambda)

  list(P = pair_res$P,           # changed from assignment
       Q = pair_res$Q,           # direct mapping
       iterations = pair_res$iterations)
}

cone_align_multiple_fit <- function(strata, proc, ref_idx,
                                   ncomp, sigma, lambda, use_laplacian,
                                   solver, max_iter, tol, knn,
                                   feature_blocks) {
  
  m <- length(strata)
  
  # Step 1: Compute embeddings for all graphs
  embeddings <- compute_cone_embeddings(strata, ncomp, sigma, use_laplacian, knn)
  
  # Check for NULL embeddings (embedding computation failures)
  if (any(sapply(embeddings, is.null))) {
    stop("Embedding computation failed for one or more graphs.", call. = FALSE)
  }
  
  # -----------------------------------------------
  # Step 2: align each graph to the reference graph
  # -----------------------------------------------
  alignments <- lapply(seq_len(m), function(g) {
    if (g == ref_idx) {
      list(P = seq_len(nrow(embeddings[[g]])),
           Q = diag(ncol(embeddings[[g]])),
           iterations = 0L)
    } else {
      .align_two_embeddings(embeddings[[ref_idx]], embeddings[[g]],
                            solver, max_iter, tol, lambda)
    }
  })

  # Step 3: Extract results and construct the final object
  P <- lapply(alignments, `[[`, "P")
  Q <- lapply(alignments, `[[`, "Q")
  
  # Apply the transformations to get the final aligned embeddings
  s_list <- lapply(seq_len(m), function(g) {
    if (is.null(alignments[[g]])) {
      return(NULL)
    }

    Q_g <- alignments[[g]]$Q

    if (g == ref_idx) {
      embeddings[[g]] %*% Q_g
    } else {
      embeddings[[g]] %*% t(Q_g)
    }
  })
  
  # Filter out NULLs before stacking
  s_list <- s_list[!sapply(s_list, is.null)]
  
  # Stack the aligned embeddings
  s <- do.call(rbind, s_list)

  # Compute primal vectors (following cone_align pattern)
  loadings <- do.call(rbind, lapply(seq_len(m), function(i) {
    xi <- strata[[i]]$x
    alpha_i <- embeddings[[i]]
    Matrix::crossprod(xi, alpha_i)
  }))

  # Check for convergence messages
  iterations <- vapply(alignments, function(x) x$iterations %||% 0, integer(1))
  if (any(iterations >= max_iter)) {
      message("Maximum iterations reached in one or more pairwise alignments.")
  }

  # Construct the final multiblock_biprojector object
  rotations_to_ref <- lapply(seq_len(m), function(g) {
    if (g == ref_idx) {
      diag(ncol(embeddings[[g]]))
    } else {
      t(Q[[g]])
    }
  })
  
  new_alignment_result(
    scores = s,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "cone_align_multiple",
    extras = list(
      assignment = P,
      rotation = rotations_to_ref,
      n_domains = m,
      ref_idx = ref_idx,
      iterations = max(iterations)
    )
  )
}
