# RAPID-MA solver control --------------------------------------------------

# Resolve one RAPID-MA settings object. Every path is revalidated by the
# constructor: a plain named list, and also an already classed object (which
# may have been edited after construction). Fields are merged one level deep
# into the constructor defaults, so nested component lists reach their own
# resolver unmerged and are validated there; unknown fields are an error.
# A `NULL` field means "use the default", matching `utils::modifyList()`.
.rapid_resolve_settings <- function(control, constructor, class, label,
                                    type_message) {
  if (is.null(control)) return(constructor())
  if (inherits(control, class)) control <- unclass(control)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(!nzchar(names(control)))))) {
    stop(type_message, call. = FALSE)
  }
  defaults <- unclass(constructor())
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown ", label, " control field(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(names(control))) {
    stop("Duplicated ", label, " control field(s): ",
         paste(unique(names(control)[duplicated(names(control))]),
               collapse = ", "), ".", call. = FALSE)
  }
  control <- control[!vapply(control, is.null, logical(1))]
  arguments <- defaults[setdiff(names(defaults), names(control))]
  arguments <- arguments[!vapply(arguments, is.null, logical(1))]
  do.call(constructor, c(arguments, control))
}

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
  .rapid_resolve_settings(
    control, .rapid_relation_control, "rapid_ma_relation_control", "relation",
    "`relation` must be NULL, a named list, or rapid_ma_relation_control."
  )
}

#' Control RAPID-MA Optimization and Bounded Sparse Components
#'
#' Validates the alternating solver settings and resolves all nested sparse
#' graph, diffusion, prototype, transport, and relation-matching controls.
#'
#' @details
#' Each nested component (`relation`, `diffusion`, `prototype`, `uot`,
#' `relmatch`) accepts `NULL` (defaults) or a named list of fields. Every path,
#' including a plain named list passed as `control` to [rapid_ma()] and a
#' previously constructed control object that was edited afterwards, is
#' revalidated field by field; unknown field names are an error.
#'
#' Fields of the nested components:
#' * `relation`: `build_feature`, `feature_k`, `feature_sketch_dim`,
#'   `spatial_mode` (`"relative"`, `"common"`, `"ignore"`), `spatial_k`,
#'   `spatial_radius`, `spatial_sigma`, `attribute_mode` (`"relation"`,
#'   `"shared"`, `"ignore"`), `attribute_k`, `attribute_hash_dim`,
#'   `normalization` (`"random_walk"`, `"symmetric"`), `weight_mode`
#'   (`"heat"`, `"binary"`), `degree_cap`, `custom_symmetrize`, `dense_max_n`.
#'   When `feature_sketch_dim` is `NULL` the feature sketch width is
#'   `min(ncol(X), max(128, 4 * ncomp))`. Generated neighbour graphs are
#'   exact. Views with at least five columns use a blocked exact search that
#'   breaks ties at the `k`-th distance by lexicographic row order, so the
#'   graph does not depend on row order; narrower views (for example 2-D
#'   coordinates) use the RANN kd-tree. The fitted
#'   `relations[[domain]]$metadata[[relation]]$neighbor_backend` records which
#'   search was used.
#' * `diffusion`: `steps`, `sketch_dim`, `output_dim`, `gate` (`"energy"`,
#'   `"uniform"`, `"fixed"`), `relation_weights`, `gate_temperature`,
#'   `label_gate_weight`, `include_position`, `include_attribute`,
#'   `include_signatures`, `direction_mode` (`"bidirectional"`, `"forward"`),
#'   `reverse_degree_cap`, `block_size`, `ridge`, `store_propagations`, `seed`.
#'   `output_dim` and `seed` are overwritten by [rapid_ma()] with `ncomp` and
#'   the top-level `seed`, and `sketch_dim` is raised to at least the feature
#'   sketch width and `ncomp` (it is a floor, not a cap). In
#'   `"bidirectional"` mode a reverse propagation channel is added only for
#'   relations that are directed before normalization.
#' * `prototype`: `prototypes_per_class`, `unknown_prototypes`,
#'   `unknown_level`, `embedding_weight`, `structure_weight`,
#'   `attribute_weight`, `min_mass`, `class_smoothing`, `max_recoveries`,
#'   `seed`. Prototype initialization is deterministic farthest-point seeding
#'   in canonical row order and draws no random numbers; `seed` is only
#'   recorded in the prototype diagnostics.
#' * `uot`: `q`, `epsilon`, `rho_source`, `rho_target`, `latent_weight`,
#'   `structure_weight`, `attribute_weight`, `position_weight`, `label_weight`,
#'   `label_confidence_threshold`, `hard_labels`, `anchor_weight`,
#'   `hard_anchors`, `block_size`, `ensure_prototype_coverage`, `cost_cap`,
#'   `max_iter`, `tol`, `require_convergence`. [rapid_ma()] folds position
#'   coordinates into the structure cost view (weighted by
#'   `structure_weight`); `position_weight` applies only to a separate
#'   position view, which [rapid_ma()] never supplies, so a non-default value
#'   triggers a warning there.
#' * `relmatch`: `normalization` (`"mass"`, `"frobenius"`, `"symmetric"`),
#'   `storage`, `sparse_threshold`, `relation_gate` (`"robust"`, `"uniform"`,
#'   `"fixed"`), `relation_weights`, `gate_temperature`, `min_mass`,
#'   `dense_coupling_max_entries`.
#'
#' @param max_iter Maximum accepted alternating updates.
#' @param min_iter Minimum updates before convergence can be declared.
#' @param tol Relative objective improvement tolerance.
#' @param objective_tolerance Absolute tolerance for monotone acceptance.
#' @param max_backtracks Maximum damped retries after an objective increase.
#' @param backtrack_factor Multiplicative step reduction in `(0, 1)`.
#' @param map_type Shared-coordinate map: `"ridge"`, `"orthogonal"`, or
#'   `"identity"`. Orthogonal maps stay orthogonal under damped updates: the
#'   blended coefficient matrix is replaced by its polar factor.
#' @param map_ridge Positive ridge penalty for small latent map solves.
#' @param map_step,prototype_step Initial map and prototype update fractions in
#'   `(0, 1]`.
#' @param uot_weight,structure_weight,supervised_weight Nonnegative objective
#'   component weights.
#' @param relation,diffusion,prototype,uot,relmatch Optional named nested
#'   control lists for the corresponding RAPID-MA components; see Details.
#' @param seed Nonnegative deterministic seed for hashed sketches.
#' @param verbose Whether to emit iteration progress.
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
  .rapid_resolve_settings(
    control, rapid_ma_control, "rapid_ma_control", "RAPID-MA",
    "`control` must be NULL, a named list, or rapid_ma_control()."
  )
}
