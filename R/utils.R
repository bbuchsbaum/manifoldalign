#' Utility functions for manifoldalign package
#' 
#' This file contains shared utility functions used across multiple modules.
#' 
#' @section Usage of safe_compute:
#' 
#' The \code{safe_compute} function should be used for operations that:
#' \itemize{
#'   \item May fail due to external dependencies or numerical issues
#'   \item Should stop execution with a clear error message (not provide fallbacks)
#'   \item Benefit from consistent error reporting across the package
#' }
#' 
#' Do NOT use \code{safe_compute} for:
#' \itemize{
#'   \item Operations that need custom fallback behavior
#'   \item Cases where you want to continue execution after an error
#'   \item Simple parameter validation (use chk:: functions instead)
#' }
#' 
#' Examples of good candidates:
#' \itemize{
#'   \item Eigenvalue computations: \code{safe_compute(PRIMME::eigs_sym(...), "Eigenvalue computation failed")}
#'   \item Graph construction: \code{safe_compute(neighborweights::graph_weights(...), "Graph construction failed")}
#'   \item Matrix operations: \code{safe_compute(solve(A), "Matrix solve failed")}
#' }

#' Safe computation wrapper with enhanced error handling
#'
#' Wraps expressions in tryCatch with informative error messages and optional
#' warning handling. Provides consistent error reporting across the package.
#'
#' @param expr Expression to evaluate
#' @param error_msg Custom error message to display on failure
#' @param warning_fn Optional function to handle warnings (default: standard warning)
#' @return Result of expression evaluation
#' @keywords internal
safe_compute <- function(expr, error_msg, warning_fn = NULL) {
  tryCatch(
    expr,
    error = function(e) {
      stop(error_msg, " Original error: ", e$message, call. = FALSE)
    },
    warning = function(w) {
      if (!is.null(warning_fn)) {
        warning_fn(w)
      } else {
        warning(w)
      }
      restart <- findRestart("muffleWarning")
      if (!is.null(restart)) {
        invokeRestart(restart)
      }
    }
  )
}

#' Resolve a hyperdesign object into a normalized domain list
#'
#' Guarantees a consistent list-of-lists with entries `x` and `design` while
#' preserving original names for downstream use.
#'
#' @param data Hyperdesign object produced by multidesign::hyperdesign() or a
#'   compatible structure used in tests.
#' @return List with `domains`, `domain_names`, and `n_domains`.
#' @keywords internal
resolve_hyperdesign <- function(data) {
  if (!inherits(data, "hyperdesign")) {
    stop("Expected an object of class 'hyperdesign'", call. = FALSE)
  }

  # Preferred structure from multidesign::hyperdesign stores blocks.
  if (!is.null(data$blocks)) {
    domains_raw <- data$blocks
  } else {
    domains_raw <- unclass(data)
  }

  if (!length(domains_raw)) {
    stop("Hyperdesign object contains no domains", call. = FALSE)
  }

  domain_names <- names(domains_raw)
  if (is.null(domain_names) || any(domain_names == "")) {
    domain_names <- paste0("domain", seq_along(domains_raw))
  }

  domains <- lapply(domains_raw, function(dom) {
    if (inherits(dom, "multidesign")) {
      x <- dom$x
      design <- dom$design
    } else if (is.list(dom) && !is.null(dom$x)) {
      x <- dom$x
      design <- dom$design
    } else if (is.matrix(dom) || is.data.frame(dom)) {
      x <- dom
      design <- NULL
    } else {
      stop("Unsupported domain format in hyperdesign", call. = FALSE)
    }

    if (!is.null(x) && !is.matrix(x)) {
      x <- as.matrix(x)
    }
    if (is.null(x)) {
      stop("Domain is missing numeric data matrix 'x'", call. = FALSE)
    }
    list(x = x, design = design)
  })

  list(domains = domains, domain_names = domain_names, n_domains = length(domains))
}

#' Compute feature block indices for concatenated matrices
#'
#' @param domains List of domain objects with element `$x`
#' @return Named list of integer sequences matching the column ranges per block
#' @keywords internal
feature_block_indices <- function(domains) {
  stopifnot(is.list(domains), length(domains) > 0)

  features_per_block <- vapply(domains, function(block) ncol(block$x), integer(1))
  end_idx <- cumsum(features_per_block)
  start_idx <- c(1L, head(end_idx, -1) + 1L)

  out <- lapply(seq_along(features_per_block), function(i) {
    seq.int(start_idx[i], end_idx[i])
  })
  names(out) <- names(domains)
  out
}

#' Construct an alignment result with consistent S3 metadata
#'
#' @param scores Matrix of concatenated scores (samples x components)
#' @param loadings Matrix of loadings/primal vectors
#' @param preproc Pre-processing object or NULL
#' @param feature_blocks Feature block indices (list or matrix)
#' @param subclass Character vector of additional classes to prepend
#' @param extras Named list of extra slots to attach to the result
#' @return Object inheriting from multiblock_biprojector and subclass
#' @keywords internal
new_alignment_result <- function(scores,
                                 loadings,
                                 preproc = NULL,
                                 feature_blocks,
                                 subclass,
                                 extras = list()) {
  if (!is.matrix(scores)) {
    scores <- as.matrix(scores)
  }
  if (!is.matrix(loadings)) {
    loadings <- as.matrix(loadings)
  }

  sdev <- if (ncol(scores) > 0) apply(scores, 2, stats::sd) else numeric(0)

  if (is.null(preproc)) {
    result <- c(list(
      v = loadings,
      s = scores,
      sdev = sdev,
      preproc = NULL,
      block_indices = feature_blocks
    ), extras)
    class(result) <- unique(c(subclass, "multiblock_biprojector"))
  } else {
    result <- multivarious::multiblock_biprojector(
      v = loadings,
      s = scores,
      sdev = sdev,
      preproc = preproc,
      block_indices = feature_blocks
    )
    for (nm in names(extras)) {
      result[[nm]] <- extras[[nm]]
    }
    class(result) <- unique(c(subclass, class(result)))
  }

  result
}
