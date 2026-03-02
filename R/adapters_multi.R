#' Multi-view aligner adapters (native delegates)
#'
#' These adapters declare native multi-domain support and delegate
#' to existing multi-view functions in this package.

#' KEMA aligner (native multi-view)
#'
#' @return An object of class `c("kema_aligner", "aligner")` representing a
#'   KEMA alignment algorithm descriptor.
#' @examples
#' algo <- kema_aligner()
#' aligner_capabilities(algo)
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
#'
#' @return An object of class `c("gpa_aligner", "aligner")` representing a
#'   Generalized Procrustes alignment algorithm descriptor.
#' @examples
#' algo <- gpa_aligner()
#' aligner_capabilities(algo)
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
#'
#' @return An object of class `c("cone_multi_aligner", "aligner")` representing
#'   a CONE-Align multiple alignment algorithm descriptor.
#' @examples
#' algo <- cone_multi_aligner()
#' aligner_capabilities(algo)
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

#' Global geometry aligner (native multi-view)
#'
#' @return An object of class `c("global_geo_aligner", "aligner")` representing
#'   a global geometry alignment algorithm descriptor.
#' @examples
#' algo <- global_geo_aligner()
#' class(algo)
#' @export
global_geo_aligner <- function() {
  new_aligner("global_geo", group = "GL", supports_multi = TRUE)
}

#' @export
fit_many.global_geo_aligner <- function(algo, domains, ...) {
  if (inherits(domains, "hyperdesign")) {
    return(global_geo_align(domains, ...))
  }
  if (is.list(domains) && all(vapply(domains, is.matrix, logical(1)))) {
    return(global_geo_align(domains, ...))
  }
  stop("global_geo_aligner: fit_many expects a hyperdesign or list of matrices.")
}

#' Spectral MNN aligner (native multi-view)
#'
#' @return An object of class `c("spectral_mnn_aligner", "aligner")` representing
#'   a spectral MNN alignment algorithm descriptor.
#' @examples
#' algo <- spectral_mnn_aligner()
#' aligner_capabilities(algo)
#' @export
spectral_mnn_aligner <- function() {
  new_aligner("spectral_mnn", group = "O", supports_multi = TRUE)
}

#' @export
fit_many.spectral_mnn_aligner <- function(algo, domains, ...) {
  if (inherits(domains, "hyperdesign")) {
    return(spectral_mnn_align(domains, ...))
  }
  if (is.list(domains) && all(vapply(domains, is.matrix, logical(1)))) {
    return(spectral_mnn_align(domains, ...))
  }
  stop("spectral_mnn_aligner: fit_many expects a hyperdesign or list of matrices.")
}

#' SSMA aligner (native multi-view)
#'
#' @return An object of class `c("ssma_aligner", "aligner")` representing
#'   a semi-supervised manifold alignment (SSMA) algorithm descriptor.
#' @examples
#' algo <- ssma_aligner()
#' aligner_capabilities(algo)
#' @export
ssma_aligner <- function() {
  new_aligner("ssma", group = "GL", supports_multi = TRUE)
}

#' @export
fit_many.ssma_aligner <- function(algo, domains, ...) {
  if (inherits(domains, "hyperdesign")) {
    return(ssma_align(domains, ...))
  }
  if (is.list(domains) && all(vapply(domains, is.matrix, logical(1)))) {
    return(ssma_align(domains, ...))
  }
  stop("ssma_aligner: fit_many expects a hyperdesign or list of matrices.")
}
