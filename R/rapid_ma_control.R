# RAPID-MA solver control --------------------------------------------------

.rapid_relation_control <- function(
  build_feature = TRUE,
  feature_k = 15L,
  feature_sketch_dim = NULL,
  spatial_mode = c("relative", "common", "ignore"),
  spatial_k = 8L,
  spatial_radius = NULL,
  spatial_sigma = NULL,
  attribute_mode = c("relation", "shared", "ignore"),
  attribute_k = 15L,
  attribute_hash_dim = 32L,
  normalization = c("random_walk", "symmetric"),
  weight_mode = c("heat", "binary"),
  degree_cap = 32L,
  custom_symmetrize = c("preserve", "mutual", "none"),
  dense_max_n = 2000L
) {
  spatial_mode <- match.arg(spatial_mode)
  attribute_mode <- match.arg(attribute_mode)
  normalization <- match.arg(normalization)
  weight_mode <- match.arg(weight_mode)
  custom_symmetrize <- match.arg(custom_symmetrize)
  if (length(build_feature) != 1L || !is.logical(build_feature) ||
      is.na(build_feature)) {
    stop("`build_feature` must be TRUE or FALSE.", call. = FALSE)
  }
  feature_k <- .rapid_int_scalar(feature_k, "feature_k")
  spatial_k <- .rapid_int_scalar(spatial_k, "spatial_k")
  attribute_k <- .rapid_int_scalar(attribute_k, "attribute_k")
  attribute_hash_dim <- .rapid_int_scalar(
    attribute_hash_dim, "attribute_hash_dim"
  )
  degree_cap <- .rapid_int_scalar(degree_cap, "degree_cap")
  dense_max_n <- .rapid_int_scalar(dense_max_n, "dense_max_n")
  if (!is.null(feature_sketch_dim)) {
    feature_sketch_dim <- .rapid_int_scalar(
      feature_sketch_dim, "feature_sketch_dim"
    )
  }
  for (name in c("spatial_radius", "spatial_sigma")) {
    value <- get(name)
    if (!is.null(value)) {
      value <- .rapid_number_scalar(value, name, lower = 0, strict = TRUE)
      assign(name, value)
    }
  }
  if (degree_cap > 4096L) {
    stop("`degree_cap` must be <= 4096.", call. = FALSE)
  }
  structure(
    list(
      build_feature = build_feature,
      feature_k = feature_k,
      feature_sketch_dim = feature_sketch_dim,
      spatial_mode = spatial_mode,
      spatial_k = spatial_k,
      spatial_radius = spatial_radius,
      spatial_sigma = spatial_sigma,
      attribute_mode = attribute_mode,
      attribute_k = attribute_k,
      attribute_hash_dim = attribute_hash_dim,
      normalization = normalization,
      weight_mode = weight_mode,
      degree_cap = degree_cap,
      custom_symmetrize = custom_symmetrize,
      dense_max_n = dense_max_n
    ),
    class = "rapid_ma_relation_control"
  )
}

.rapid_resolve_relation_control <- function(control = NULL) {
  defaults <- unclass(.rapid_relation_control())
  if (is.null(control)) return(do.call(.rapid_relation_control, defaults))
  if (inherits(control, "rapid_ma_relation_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`relation` must be NULL, a named list, or rapid_ma_relation_control.",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown relation control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  do.call(.rapid_relation_control, utils::modifyList(defaults, control))
}

#' Control RAPID-MA optimization and bounded sparse components
#'
#' @param max_iter Maximum accepted alternating updates.
#' @param min_iter Minimum updates before convergence can be declared.
#' @param tol Relative objective improvement tolerance.
#' @param objective_tolerance Absolute tolerance for monotone acceptance.
#' @param max_backtracks Maximum damped retries after an objective increase.
#' @param backtrack_factor Multiplicative step reduction in `(0, 1)`.
#' @param map_type Shared-coordinate map: `"ridge"`, `"orthogonal"`, or
#'   `"identity"`.
#' @param map_ridge Positive ridge penalty for small latent map solves.
#' @param map_step Initial map update fraction in `(0, 1]`.
#' @param prototype_step Initial prototype update fraction in `(0, 1]`.
#' @param uot_weight,structure_weight,supervised_weight Nonnegative objective
#'   component weights.
#' @param relation,diffusion,prototype,uot,relmatch Named nested control lists.
#' @param seed Nonnegative deterministic seed.
#' @param verbose Emit iteration progress.
#' @return A `rapid_ma_control` object.
#' @export
rapid_ma_control <- function(
  max_iter = 8L,
  min_iter = 1L,
  tol = 1e-4,
  objective_tolerance = 1e-8,
  max_backtracks = 4L,
  backtrack_factor = 0.5,
  map_type = c("ridge", "orthogonal", "identity"),
  map_ridge = 0.05,
  map_step = 0.75,
  prototype_step = 1,
  uot_weight = 1,
  structure_weight = 0.25,
  supervised_weight = 0.5,
  relation = NULL,
  diffusion = NULL,
  prototype = NULL,
  uot = NULL,
  relmatch = NULL,
  seed = 1L,
  verbose = FALSE
) {
  map_type <- match.arg(map_type)
  max_iter <- .rapid_int_scalar(max_iter, "max_iter")
  min_iter <- .rapid_int_scalar(min_iter, "min_iter", min_value = 0L)
  max_backtracks <- .rapid_int_scalar(
    max_backtracks, "max_backtracks", min_value = 0L
  )
  seed <- .rapid_int_scalar(seed, "seed", min_value = 0L)
  tol <- .rapid_number_scalar(tol, "tol", lower = 0)
  objective_tolerance <- .rapid_number_scalar(
    objective_tolerance, "objective_tolerance", lower = 0
  )
  backtrack_factor <- .rapid_number_scalar(
    backtrack_factor, "backtrack_factor", lower = 0, strict = TRUE
  )
  if (backtrack_factor >= 1) {
    stop("`backtrack_factor` must be < 1.", call. = FALSE)
  }
  map_ridge <- .rapid_number_scalar(map_ridge, "map_ridge", lower = 0,
                                    strict = TRUE)
  map_step <- .rapid_number_scalar(map_step, "map_step", lower = 0,
                                   strict = TRUE)
  prototype_step <- .rapid_number_scalar(
    prototype_step, "prototype_step", lower = 0, strict = TRUE
  )
  if (map_step > 1 || prototype_step > 1) {
    stop("`map_step` and `prototype_step` must be <= 1.", call. = FALSE)
  }
  weights <- c(
    uot_weight = uot_weight,
    structure_weight = structure_weight,
    supervised_weight = supervised_weight
  )
  for (name in names(weights)) {
    weights[[name]] <- .rapid_number_scalar(weights[[name]], name, lower = 0)
  }
  if (sum(weights) <= 0) {
    stop("At least one core objective weight must be positive.", call. = FALSE)
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  if (min_iter > max_iter) {
    stop("`min_iter` cannot exceed `max_iter`.", call. = FALSE)
  }

  structure(
    list(
      max_iter = max_iter,
      min_iter = min_iter,
      tol = tol,
      objective_tolerance = objective_tolerance,
      max_backtracks = max_backtracks,
      backtrack_factor = backtrack_factor,
      map_type = map_type,
      map_ridge = map_ridge,
      map_step = map_step,
      prototype_step = prototype_step,
      uot_weight = weights[["uot_weight"]],
      structure_weight = weights[["structure_weight"]],
      supervised_weight = weights[["supervised_weight"]],
      relation = .rapid_resolve_relation_control(relation),
      diffusion = .rapid_resolve_diffusion_control(diffusion),
      prototype = .rapid_resolve_prototype_control(prototype),
      uot = .rapid_resolve_uot_control(uot),
      relmatch = .rapid_resolve_relmatch_control(relmatch),
      seed = seed,
      verbose = verbose
    ),
    class = "rapid_ma_control"
  )
}

.rapid_resolve_control <- function(control = NULL) {
  defaults <- unclass(rapid_ma_control())
  if (is.null(control)) return(do.call(rapid_ma_control, defaults))
  if (inherits(control, "rapid_ma_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`control` must be NULL, a named list, or rapid_ma_control().",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown RAPID-MA control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  merged <- utils::modifyList(defaults, control)
  do.call(rapid_ma_control, merged)
}
