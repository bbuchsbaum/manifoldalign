# Multi-view aligner adapters (native delegates)
#
# These adapters declare native multi-domain support and delegate
# to existing multi-view functions in this package. Their `fit_many()`
# methods translate the generic latent dimension `k` (as passed by
# [align_many()]) into the delegate's `ncomp` argument; passing both `k` and
# a different `ncomp` is an error.

# Resolve the generic `k` and the method-specific `ncomp` into one value.
.native_ncomp <- function(k, ncomp, aligner) {
  if (!is.null(k) && !is.null(ncomp) &&
      !isTRUE(all.equal(as.numeric(k), as.numeric(ncomp)))) {
    stop(aligner, ": `k` (", paste(k, collapse = ", "), ") and `ncomp` (",
         paste(ncomp, collapse = ", "), ") disagree; supply one of them.",
         call. = FALSE)
  }
  if (is.null(ncomp)) k else ncomp
}

#' KEMA aligner (native multi-view)
#'
#' The capability record declares a nonlinear (kernel) encoder whose
#' out-of-sample projection is evaluated from the fitted object alone, with no
#' new-domain calibration and no inverse map.
#'
#' @return An object of class `c("kema_aligner", "aligner")` representing a
#'   KEMA alignment algorithm descriptor.
#' @examples
#' algo <- kema_aligner()
#' aligner_capabilities(algo)
#' @export
kema_aligner <- function() {
  new_aligner("kema", group = "kernel", supports_multi = TRUE,
              encoder = "nonlinear",
              supports_predict = TRUE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = FALSE)
}

#' @export
fit_many.kema_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  ncomp <- .native_ncomp(k, ncomp, "kema_aligner")
  # Prefer native KEMA on hyperdesign objects
  if (inherits(domains, "hyperdesign")) {
    if (is.null(ncomp)) return(kema(domains, ...))
    return(kema(domains, ncomp = ncomp, ...))
  }
  stop("kema_aligner: fit_many expects a 'hyperdesign' object as 'domains'. ",
       "Construct one via multidesign::hyperdesign(), with a label column for 'y'.")
}

#' Generalized Procrustes aligner (native multi-view)
#'
#' Generalized Procrustes rotates every domain's feature space onto a common
#' consensus, so the output dimension is the shared feature dimension and no
#' latent dimension `k` can be chosen. The capability record declares a
#' linear, invertible (orthogonal) map without an out-of-sample prediction
#' method or new-domain calibration.
#'
#' @return An object of class `c("gpa_aligner", "aligner")` representing a
#'   Generalized Procrustes alignment algorithm descriptor.
#' @examples
#' algo <- gpa_aligner()
#' aligner_capabilities(algo)
#' @export
gpa_aligner <- function() {
  new_aligner("gpa", group = "O", supports_multi = TRUE,
              encoder = "linear",
              supports_predict = FALSE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = TRUE)
}

#' @export
fit_many.gpa_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  if (!is.null(k) || !is.null(ncomp)) {
    stop("gpa_aligner: generalized Procrustes keeps the shared feature ",
         "dimension and has no latent dimension; do not supply `k`/`ncomp`.",
         call. = FALSE)
  }
  if (inherits(domains, "hyperdesign")) {
    return(generalized_procrustes(domains, ...))
  }
  stop("gpa_aligner: fit_many expects a 'hyperdesign' object as 'domains'. ",
       "See generalized_procrustes.hyperdesign for required inputs (task labels).")
}

#' CONE-Align multiple aligner (native multi-view)
#'
#' CONE-Align embeds each domain's neighbourhood graph and matches nodes, so
#' it is transductive: the capability record declares no feature encoder and
#' no out-of-sample prediction.
#'
#' @return An object of class `c("cone_multi_aligner", "aligner")` representing
#'   a CONE-Align multiple alignment algorithm descriptor.
#' @examples
#' algo <- cone_multi_aligner()
#' aligner_capabilities(algo)
#' @export
cone_multi_aligner <- function() {
  new_aligner("cone_multi", group = "O", supports_multi = TRUE,
              encoder = "none",
              supports_predict = FALSE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = FALSE)
}

#' @export
fit_many.cone_multi_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  ncomp <- .native_ncomp(k, ncomp, "cone_multi_aligner")
  if (inherits(domains, "hyperdesign") ||
      (is.list(domains) && all(vapply(domains, is.matrix, logical(1))))) {
    if (is.null(ncomp)) return(cone_align_multiple(domains, ...))
    return(cone_align_multiple(domains, ncomp = ncomp, ...))
  }
  stop("cone_multi_aligner: fit_many expects a hyperdesign or list of matrices.")
}

#' Global geometry aligner (native multi-view)
#'
#' The capability record declares per-domain linear maps with out-of-sample
#' projection of new rows, no new-domain calibration and no inverse map.
#'
#' @return An object of class `c("global_geo_aligner", "aligner")` representing
#'   a global geometry alignment algorithm descriptor.
#' @examples
#' algo <- global_geo_aligner()
#' class(algo)
#' @export
global_geo_aligner <- function() {
  new_aligner("global_geo", group = "GL", supports_multi = TRUE,
              encoder = "linear",
              supports_predict = TRUE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = FALSE)
}

#' @export
fit_many.global_geo_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  ncomp <- .native_ncomp(k, ncomp, "global_geo_aligner")
  if (inherits(domains, "hyperdesign") ||
      (is.list(domains) && all(vapply(domains, is.matrix, logical(1))))) {
    if (is.null(ncomp)) return(global_geo_align(domains, ...))
    return(global_geo_align(domains, ncomp = ncomp, ...))
  }
  stop("global_geo_aligner: fit_many expects a hyperdesign or list of matrices.")
}

#' Spectral MNN aligner (native multi-view)
#'
#' Spectral MNN alignment embeds the joint graph spectrally and is
#' transductive: the capability record declares no feature encoder and no
#' validated out-of-sample prediction.
#'
#' @return An object of class `c("spectral_mnn_aligner", "aligner")` representing
#'   a spectral MNN alignment algorithm descriptor.
#' @examples
#' algo <- spectral_mnn_aligner()
#' aligner_capabilities(algo)
#' @export
spectral_mnn_aligner <- function() {
  new_aligner("spectral_mnn", group = "O", supports_multi = TRUE,
              encoder = "none",
              supports_predict = FALSE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = FALSE)
}

#' @export
fit_many.spectral_mnn_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  ncomp <- .native_ncomp(k, ncomp, "spectral_mnn_aligner")
  if (inherits(domains, "hyperdesign") ||
      (is.list(domains) && all(vapply(domains, is.matrix, logical(1))))) {
    if (is.null(ncomp)) return(spectral_mnn_align(domains, ...))
    return(spectral_mnn_align(domains, ncomp = ncomp, ...))
  }
  stop("spectral_mnn_aligner: fit_many expects a hyperdesign or list of matrices.")
}

#' SSMA aligner (native multi-view)
#'
#' The capability record declares per-domain linear projections (scores are
#' `X_i F_i`) with out-of-sample projection of new rows, no new-domain
#' calibration and no inverse map.
#'
#' @return An object of class `c("ssma_aligner", "aligner")` representing
#'   a semi-supervised manifold alignment (SSMA) algorithm descriptor.
#' @examples
#' algo <- ssma_aligner()
#' aligner_capabilities(algo)
#' @export
ssma_aligner <- function() {
  new_aligner("ssma", group = "GL", supports_multi = TRUE,
              encoder = "linear",
              supports_predict = TRUE,
              supports_calibration = FALSE,
              requires_auxiliary_at_predict = FALSE,
              invertible = FALSE)
}

#' @export
fit_many.ssma_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  ncomp <- .native_ncomp(k, ncomp, "ssma_aligner")
  if (inherits(domains, "hyperdesign") ||
      (is.list(domains) && all(vapply(domains, is.matrix, logical(1))))) {
    if (is.null(ncomp)) return(ssma_align(domains, ...))
    return(ssma_align(domains, ncomp = ncomp, ...))
  }
  stop("ssma_aligner: fit_many expects a hyperdesign or list of matrices.")
}
