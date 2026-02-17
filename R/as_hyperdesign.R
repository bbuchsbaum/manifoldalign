#' Construct a lightweight hyperdesign from matrices and labels
#'
#' Builds a minimal 'hyperdesign' object (list of domains with fields `x` and
#' `design`) from a list of matrices and optional label vectors. This avoids
#' requiring the multidesign package for simple use cases.
#'
#' @param X_list List of matrices (samples x features) per domain
#' @param labels Optional list of label vectors aligned to rows of each matrix
#' @param label_name Name for the label column in each domain's design
#' @param domain_names Optional names for the domains
#' @return An object of class 'hyperdesign' suitable for kema()/generalized_procrustes()
#' @examples
#' # Create hyperdesign from list of matrices
#' X1 <- matrix(rnorm(100), 50, 2)
#' X2 <- matrix(rnorm(100), 50, 2)
#' X3 <- matrix(rnorm(100), 50, 2)
#' hd <- as_hyperdesign(list(X1, X2, X3))
#'
#' # With labels
#' labels <- list(
#'   rep(c("A", "B"), each = 25),
#'   rep(c("A", "B"), each = 25),
#'   rep(c("A", "B"), each = 25)
#' )
#' hd_labeled <- as_hyperdesign(list(X1, X2, X3), labels = labels)
#' @export
as_hyperdesign <- function(X_list, labels = NULL, label_name = "label", domain_names = NULL) {
  if (!is.list(X_list) || length(X_list) == 0) {
    stop("X_list must be a non-empty list of matrices", call. = FALSE)
  }
  n <- length(X_list)
  if (!is.null(labels)) {
    if (!is.list(labels) || length(labels) != n) {
      stop("labels must be NULL or a list of the same length as X_list", call. = FALSE)
    }
  }
  if (!is.null(domain_names)) {
    if (length(domain_names) != n) stop("domain_names must have length equal to X_list")
  } else {
    domain_names <- paste0("domain", seq_len(n))
  }

  strata <- lapply(seq_len(n), function(i) {
    Xi <- X_list[[i]]
    if (!is.matrix(Xi)) Xi <- as.matrix(Xi)
    if (!is.numeric(Xi)) stop("All X_list elements must be numeric matrices", call. = FALSE)
    di <- if (!is.null(labels)) labels[[i]] else NULL
    if (!is.null(di)) {
      if (length(di) != nrow(Xi)) stop("labels[[", i, "]] length must match nrow(X_list[[", i, "]])", call. = FALSE)
      design <- data.frame(setNames(list(di), label_name))
    } else {
      design <- data.frame()
    }
    list(x = Xi, design = design)
  })
  names(strata) <- domain_names
  class(strata) <- "hyperdesign"
  strata
}

