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
      invokeRestart("muffleWarning")
    }
  )
} 