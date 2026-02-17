#' Additional methods for FPGW objects
#'
#' @description
#' These methods integrate FPGW results with the multiblock_biprojector
#' framework from the multivarious package.
#'
#' @return NULL (documentation page only).
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' # Use methods: summary(res), plot(res), coef(res)
#' }
#' @name fpgw-methods
NULL

#' Transform new data using FPGW alignment
#'
#' @param _data An fpgw object
#' @param newdata New data to transform
#' @param source_index Index of source domain
#' @param target_index Index of target domain
#' @param ... Additional arguments
#'
#' @return Transformed samples in the target domain feature space.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' X_new <- matrix(rnorm(15), 5, 3)
#' transformed <- transform(res, X_new, source_index = 1, target_index = 2)
#' }
#'
#' @export
transform.fpgw <- function(`_data`, newdata, source_index = 1, target_index = 2, ...) {
  predict(`_data`, newdata, from = source_index, to = target_index, type = "transport", ...)
}

#' Extract coefficients from FPGW
#'
#' @param object An fpgw object
#' @param type Type of coefficients to extract
#' @param ... Additional arguments
#'
#' @return Coefficients (transport plans or parameters)
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' transport_plans <- coef(res, type = "transport")
#' params <- coef(res, type = "parameters")
#' }
#'
#' @export
coef.fpgw <- function(object, type = c("transport", "parameters"), ...) {
  type <- match.arg(type)

  if (type == "transport") {
    return(object$transport_plans)
  } else {
    # Return key parameters
    params <- list(
      omega1 = object$omega1,
      lambda = object$lambda,
      rho = object$rho,
      n_domains = length(object$domain_names)
    )
    return(params)
  }
}

#' Summary method for FPGW
#'
#' @param object An fpgw object
#' @param ... Additional arguments
#'
#' @return The fpgw object, invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' summary(res)
#' }
#'
#' @method summary fpgw
#' @export
summary.fpgw <- function(object, ...) {
  cat("Fused-Partial Gromov-Wasserstein Summary\n")
  cat("========================================\n\n")
  
  # Basic info
  cat("Number of domains:", length(object$domain_names), "\n")
  cat("Domain names:", paste(object$domain_names, collapse = ", "), "\n")
  cat("Samples per domain:", paste(object$n_samples, collapse = ", "), "\n\n")
  
  # Parameters
  cat("Parameters:\n")
  cat("  Feature weight (omega1):", object$omega1, "\n")
  cat("  Structure weight (omega2):", 1 - object$omega1, "\n")
  
  if (!is.null(object$rho)) {
    cat("  Variant: Mass-constrained (rho =", object$rho, ")\n")
  } else if (!is.null(object$lambda) && object$lambda > 0) {
    cat("  Variant: TV-penalized (lambda =", object$lambda, ")\n")
  } else {
    cat("  Variant: Classical FGW\n")
  }
  
  cat("\nDistance Matrix:\n")
  print(round(object$distances, 4))
  
  # Transport statistics
  cat("\nTransport Plan Statistics:\n")
  for (i in seq_along(object$transport_plans)) {
    P <- object$transport_plans[[i]]
    cat(sprintf("  Plan %d: mass = %.4f, sparsity = %.2f%%\n",
                i, sum(P), 100 * mean(P < 1e-8)))
  }
  
  # Convergence
  if (!all(object$converged)) {
    cat("\nWarning: Some optimizations did not converge\n")
  }
  
  invisible(object)
}

#' Plot method for FPGW
#'
#' @param x An fpgw object
#' @param which Which plot to create (1: distance matrix, 2: transport plan)
#' @param pair For transport plan, which pair to plot (default: 1)
#' @param ... Additional arguments passed to plotting functions
#'
#' @return NULL, invisibly. Called for its side effect of producing a plot.
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' plot(res, which = 1)
#' plot(res, which = 2, pair = 1)
#' }
#'
#' @export
plot.fpgw <- function(x, which = 1, pair = 1, ...) {
  if (which == 1) {
    # Plot distance matrix as heatmap
    image(x$distances, 
          main = "FPGW Distance Matrix",
          xlab = "Domain", ylab = "Domain",
          col = heat.colors(100),
          axes = FALSE, ...)
    axis(1, at = seq(0, 1, length.out = length(x$domain_names)), 
         labels = x$domain_names)
    axis(2, at = seq(0, 1, length.out = length(x$domain_names)), 
         labels = x$domain_names)
  } else if (which == 2) {
    # Plot transport plan
    if (pair > length(x$transport_plans)) {
      stop("Invalid pair index", call. = FALSE)
    }
    
    P <- x$transport_plans[[pair]]
    image(P, 
          main = paste("Transport Plan", pair),
          xlab = "Source", ylab = "Target",
          col = heat.colors(100), ...)
  }
  invisible(NULL)
}

#' Extract fitted values
#'
#' @param object An fpgw object
#' @param ... Additional arguments
#'
#' @return List of transport plans
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' plans <- fitted(res)
#' }
#'
#' @export
fitted.fpgw <- function(object, ...) {
  object$transport_plans
}

#' Extract residuals (not applicable for FPGW)
#'
#' @param object An fpgw object
#' @param ... Additional arguments
#'
#' @return NULL with a warning
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' library(multidesign)
#' md1 <- multidesign(X1, data.frame(id = 1:10))
#' md2 <- multidesign(X2, data.frame(id = 1:10))
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' res <- fpgw(hd)
#' residuals(res)
#' }
#'
#' @export
residuals.fpgw <- function(object, ...) {
  warning("Residuals are not defined for FPGW objects", call. = FALSE)
  return(NULL)
}
