#' Control settings for multiscale manifold alignment backends
#'
#' @param enabled Logical; enable multiscale mode.
#' @param backend Character backend name. One of `"spectral"` or `"randomized_dwt"`.
#' @param operator_mode Character operator mode. Currently `"compressed_graph_core"`.
#' @param rank_per_domain Integer per-domain truncation rank for compressed operators.
#' @param svd_tol Numeric positive tolerance for truncated SVD decisions.
#' @param eigen_tol Numeric positive tolerance for non-trivial eigenvalues.
#' @param max_levels Integer maximum hierarchy depth.
#' @param eps_rank Numeric relative truncation threshold in `(0, 1]`.
#' @param eigen_solver Character eigensolver selector. One of `"lobpcg"` or `"lanczos"`.
#' @param eigen_k_max Integer cap on eigenvectors for spectral hierarchy.
#' @param knn Integer nearest-neighbor graph degree for within-domain graphs.
#' @param sigma Optional positive scalar graph bandwidth. If `NULL`, auto-tuned per domain.
#' @param cross_edge_weight Numeric non-negative weight for cross-domain correspondence edges.
#' @param mnn_k Integer nearest-neighbor count for unsupervised MNN correspondences.
#' @param max_pairs_per_label Integer cap per label/pair when building label correspondences.
#' @param max_unsup_pairs_per_pair Integer cap per domain pair in unsupervised MNN mode.
#' @param randomized_sketch_size Integer sketch size for randomized backend.
#' @param randomized_oversample Integer oversampling budget for randomized backend.
#' @param randomized_power_iters Integer number of subspace power iterations.
#' @param shift_invert_delta Numeric positive regularization for shift-invert style operators.
#' @param cg_tol Numeric positive tolerance for iterative solves.
#' @param cg_max_iter Integer maximum iterations for iterative solves.
#' @param seed Integer random seed.
#' @param verbose Logical; print backend progress messages.
#' @param scale_selection How to weight scales when assembling the multiscale
#'   dictionary. One of `"energy"`, `"supervised"`, or `"hybrid"`.
#' @param supervision_weight Weight in `[0, 1]` used when
#'   `scale_selection = "hybrid"`.
#' @param min_scale_weight Relative threshold in `[0, 1)` for keeping a scale in
#'   the final dictionary.
#' @param max_basis_per_scale Integer cap on the number of basis vectors retained
#'   from each scale after compression.
#' @param regularization Positive ridge regularization for the final dictionary
#'   generalized eigenproblem.
#' @param tune Logical; whether to tune `cross_edge_weight`, `rank_per_domain`,
#'   and `max_levels` over small candidate grids before the final fit.
#' @param candidate_cross_edge_weight Optional numeric vector of candidate
#'   cross-edge weights used when `tune = TRUE`.
#' @param candidate_rank_per_domain Optional integer vector of candidate
#'   `rank_per_domain` values used when `tune = TRUE`.
#' @param candidate_max_levels Optional integer vector of candidate `max_levels`
#'   values used when `tune = TRUE`.
#'
#' @return A list with class `multiscale_manifold_align_control`.
#' @export
multiscale_manifold_align_control <- function(
  enabled = FALSE,
  backend = c("spectral", "randomized_dwt"),
  operator_mode = c("compressed_graph_core"),
  rank_per_domain = 128L,
  svd_tol = 1e-6,
  eigen_tol = 1e-8,
  max_levels = 8L,
  eps_rank = 1e-3,
  eigen_solver = c("lobpcg", "lanczos"),
  eigen_k_max = 512L,
  knn = 12L,
  sigma = NULL,
  cross_edge_weight = 1,
  mnn_k = 5L,
  max_pairs_per_label = 100L,
  max_unsup_pairs_per_pair = 200L,
  randomized_sketch_size = 256L,
  randomized_oversample = 16L,
  randomized_power_iters = 1L,
  shift_invert_delta = 1e-3,
  cg_tol = 1e-6,
  cg_max_iter = 500L,
  scale_selection = c("energy", "supervised", "hybrid"),
  supervision_weight = 0.5,
  min_scale_weight = 0.05,
  max_basis_per_scale = 32L,
  regularization = 1e-6,
  tune = FALSE,
  candidate_cross_edge_weight = NULL,
  candidate_rank_per_domain = NULL,
  candidate_max_levels = NULL,
  seed = 1L,
  verbose = FALSE
) {
  backend <- match.arg(backend)
  operator_mode <- match.arg(operator_mode)
  eigen_solver <- match.arg(eigen_solver)
  scale_selection <- match.arg(scale_selection)

  if (!is.logical(enabled) || length(enabled) != 1L || is.na(enabled)) {
    stop("`enabled` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(rank_per_domain) || length(rank_per_domain) != 1L || is.na(rank_per_domain) || rank_per_domain < 1) {
    stop("`rank_per_domain` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(svd_tol) || length(svd_tol) != 1L || is.na(svd_tol) || svd_tol <= 0) {
    stop("`svd_tol` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(eigen_tol) || length(eigen_tol) != 1L || is.na(eigen_tol) || eigen_tol <= 0) {
    stop("`eigen_tol` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(max_levels) || length(max_levels) != 1L || is.na(max_levels) || max_levels < 1) {
    stop("`max_levels` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(eps_rank) || length(eps_rank) != 1L || is.na(eps_rank) || eps_rank <= 0 || eps_rank > 1) {
    stop("`eps_rank` must be in (0, 1].", call. = FALSE)
  }
  if (!is.numeric(eigen_k_max) || length(eigen_k_max) != 1L || is.na(eigen_k_max) || eigen_k_max < 1) {
    stop("`eigen_k_max` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(knn) || length(knn) != 1L || is.na(knn) || knn < 1) {
    stop("`knn` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.null(sigma) && (!is.numeric(sigma) || length(sigma) != 1L || is.na(sigma) || sigma <= 0)) {
    stop("`sigma` must be NULL or a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(cross_edge_weight) || length(cross_edge_weight) != 1L || is.na(cross_edge_weight) || cross_edge_weight < 0) {
    stop("`cross_edge_weight` must be a non-negative scalar.", call. = FALSE)
  }
  if (!is.numeric(mnn_k) || length(mnn_k) != 1L || is.na(mnn_k) || mnn_k < 1) {
    stop("`mnn_k` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(max_pairs_per_label) || length(max_pairs_per_label) != 1L || is.na(max_pairs_per_label) || max_pairs_per_label < 1) {
    stop("`max_pairs_per_label` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(max_unsup_pairs_per_pair) || length(max_unsup_pairs_per_pair) != 1L || is.na(max_unsup_pairs_per_pair) || max_unsup_pairs_per_pair < 1) {
    stop("`max_unsup_pairs_per_pair` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(randomized_sketch_size) || length(randomized_sketch_size) != 1L || is.na(randomized_sketch_size) || randomized_sketch_size < 1) {
    stop("`randomized_sketch_size` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(randomized_oversample) || length(randomized_oversample) != 1L || is.na(randomized_oversample) || randomized_oversample < 0) {
    stop("`randomized_oversample` must be an integer >= 0.", call. = FALSE)
  }
  if (!is.numeric(randomized_power_iters) || length(randomized_power_iters) != 1L || is.na(randomized_power_iters) || randomized_power_iters < 0) {
    stop("`randomized_power_iters` must be an integer >= 0.", call. = FALSE)
  }
  if (!is.numeric(shift_invert_delta) || length(shift_invert_delta) != 1L || is.na(shift_invert_delta) || shift_invert_delta <= 0) {
    stop("`shift_invert_delta` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(cg_tol) || length(cg_tol) != 1L || is.na(cg_tol) || cg_tol <= 0) {
    stop("`cg_tol` must be a positive scalar.", call. = FALSE)
  }
  if (!is.numeric(cg_max_iter) || length(cg_max_iter) != 1L || is.na(cg_max_iter) || cg_max_iter < 1) {
    stop("`cg_max_iter` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(supervision_weight) || length(supervision_weight) != 1L || is.na(supervision_weight) ||
      supervision_weight < 0 || supervision_weight > 1) {
    stop("`supervision_weight` must be in [0, 1].", call. = FALSE)
  }
  if (!is.numeric(min_scale_weight) || length(min_scale_weight) != 1L || is.na(min_scale_weight) ||
      min_scale_weight < 0 || min_scale_weight >= 1) {
    stop("`min_scale_weight` must be in [0, 1).", call. = FALSE)
  }
  if (!is.numeric(max_basis_per_scale) || length(max_basis_per_scale) != 1L || is.na(max_basis_per_scale) ||
      max_basis_per_scale < 1) {
    stop("`max_basis_per_scale` must be an integer >= 1.", call. = FALSE)
  }
  if (!is.numeric(regularization) || length(regularization) != 1L || is.na(regularization) || regularization <= 0) {
    stop("`regularization` must be a positive scalar.", call. = FALSE)
  }
  if (!is.logical(tune) || length(tune) != 1L || is.na(tune)) {
    stop("`tune` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(candidate_cross_edge_weight) &&
      (!is.numeric(candidate_cross_edge_weight) || !length(candidate_cross_edge_weight) ||
       any(!is.finite(candidate_cross_edge_weight)) || any(candidate_cross_edge_weight <= 0))) {
    stop("`candidate_cross_edge_weight` must be NULL or a positive numeric vector.", call. = FALSE)
  }
  if (!is.null(candidate_rank_per_domain) &&
      (!is.numeric(candidate_rank_per_domain) || !length(candidate_rank_per_domain) ||
       any(!is.finite(candidate_rank_per_domain)) || any(candidate_rank_per_domain < 1))) {
    stop("`candidate_rank_per_domain` must be NULL or a positive integer vector.", call. = FALSE)
  }
  if (!is.null(candidate_max_levels) &&
      (!is.numeric(candidate_max_levels) || !length(candidate_max_levels) ||
       any(!is.finite(candidate_max_levels)) || any(candidate_max_levels < 1))) {
    stop("`candidate_max_levels` must be NULL or an integer vector with values >= 1.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be a numeric scalar.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      enabled = enabled,
      backend = backend,
      operator_mode = operator_mode,
      rank_per_domain = as.integer(rank_per_domain),
      svd_tol = as.numeric(svd_tol),
      eigen_tol = as.numeric(eigen_tol),
      max_levels = as.integer(max_levels),
      eps_rank = as.numeric(eps_rank),
      eigen_solver = eigen_solver,
      eigen_k_max = as.integer(eigen_k_max),
      knn = as.integer(knn),
      sigma = if (is.null(sigma)) NULL else as.numeric(sigma),
      cross_edge_weight = as.numeric(cross_edge_weight),
      mnn_k = as.integer(mnn_k),
      max_pairs_per_label = as.integer(max_pairs_per_label),
      max_unsup_pairs_per_pair = as.integer(max_unsup_pairs_per_pair),
      randomized_sketch_size = as.integer(randomized_sketch_size),
      randomized_oversample = as.integer(randomized_oversample),
      randomized_power_iters = as.integer(randomized_power_iters),
      shift_invert_delta = as.numeric(shift_invert_delta),
      cg_tol = as.numeric(cg_tol),
      cg_max_iter = as.integer(cg_max_iter),
      scale_selection = scale_selection,
      supervision_weight = as.numeric(supervision_weight),
      min_scale_weight = as.numeric(min_scale_weight),
      max_basis_per_scale = as.integer(max_basis_per_scale),
      regularization = as.numeric(regularization),
      tune = tune,
      candidate_cross_edge_weight = if (is.null(candidate_cross_edge_weight)) NULL else as.numeric(candidate_cross_edge_weight),
      candidate_rank_per_domain = if (is.null(candidate_rank_per_domain)) NULL else as.integer(candidate_rank_per_domain),
      candidate_max_levels = if (is.null(candidate_max_levels)) NULL else as.integer(candidate_max_levels),
      seed = as.integer(seed),
      verbose = verbose
    ),
    class = "multiscale_manifold_align_control"
  )
}

#' @noRd
resolve_multiscale_manifold_align_control <- function(control = NULL) {
  defaults <- multiscale_manifold_align_control()

  if (missing(control) || is.null(control)) {
    return(defaults)
  }

  if (inherits(control, "multiscale_manifold_align_control")) {
    merged <- utils::modifyList(defaults, control)
    return(do.call(multiscale_manifold_align_control, merged))
  }

  if (!is.list(control)) {
    stop("`multiscale_control` must be NULL, a named list, or multiscale_manifold_align_control().", call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`multiscale_control` must be a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown multiscale control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  }

  merged <- utils::modifyList(defaults, control)
  do.call(multiscale_manifold_align_control, merged)
}
