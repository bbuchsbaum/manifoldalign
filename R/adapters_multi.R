#' Multi-view aligner adapters (native delegates)
#'
#' These adapters declare native multi-domain support and delegate
#' to existing multi-view functions in this package.

#' KEMA aligner (native multi-view)
#' @export
kema_aligner <- function() {
  new_aligner("kema", group = "kernel", supports_multi = TRUE)
}

#' @export
fit_many.kema_aligner <- function(algo, domains, ...) {
  # Prefer native KEMA on hyperdesign objects
  if (inherits(domains, "hyperdesign")) {
    return(kema(domains, ...))
  }
  stop("kema_aligner: fit_many expects a 'hyperdesign' object as 'domains'. ",
       "Construct one via multidesign::hyperdesign(), with a label column for 'y'.")
}

#' Generalized Procrustes aligner (native multi-view)
#' @export
gpa_aligner <- function() {
  new_aligner("gpa", group = "O", supports_multi = TRUE)
}

#' @export
fit_many.gpa_aligner <- function(algo, domains, ...) {
  if (inherits(domains, "hyperdesign")) {
    return(generalized_procrustes(domains, ...))
  }
  stop("gpa_aligner: fit_many expects a 'hyperdesign' object as 'domains'. ",
       "See generalized_procrustes.hyperdesign for required inputs (task labels).")
}

#' CONE-Align multiple aligner (native multi-view)
#' @export
cone_multi_aligner <- function() {
  new_aligner("cone_multi", group = "O", supports_multi = TRUE)
}

#' @export
fit_many.cone_multi_aligner <- function(algo, domains, ...) {
  if (inherits(domains, "hyperdesign")) {
    return(cone_align_multiple(domains, ...))
  }
  if (is.list(domains) && all(vapply(domains, is.matrix, logical(1)))) {
    return(cone_align_multiple(domains, ...))
  }
  stop("cone_multi_aligner: fit_many expects a hyperdesign or list of matrices.")
}

