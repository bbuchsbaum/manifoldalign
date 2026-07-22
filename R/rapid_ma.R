# RAPID-MA scalable structure-aware manifold alignment --------------------

.rapid_named_domain_argument <- function(x, domain_names, name, default = NULL) {
  M <- length(domain_names)
  if (is.null(x)) return(rep(list(default), M))
  if (!is.list(x) || length(x) != M) {
    stop("`", name, "` must be a list with one entry per domain.",
         call. = FALSE)
  }
  if (!is.null(names(x)) && all(domain_names %in% names(x))) {
    x <- x[domain_names]
  }
  unname(x)
}

.rapid_validate_domain_data <- function(data, domain_names) {
  if (!is.list(data) || length(data) < 2L) {
    stop("RAPID-MA requires at least two domains.", call. = FALSE)
  }
  lapply(seq_along(data), function(i) {
    X <- data[[i]]
    if (is.data.frame(X)) X <- as.matrix(X)
    tryCatch(
      {
        .rapid_validate_numeric_matrix(X, "features")
        X <- as.matrix(X)
        if (any(!is.finite(X))) {
          stop("features contain non-finite values", call. = FALSE)
        }
        X
      },
      error = function(e) {
        stop("Domain `", domain_names[[i]], "`: ", conditionMessage(e),
             call. = FALSE)
      }
    )
  })
}

.rapid_validate_domain_labels <- function(labels, sizes, domain_names) {
  if (is.null(labels)) return(lapply(sizes, function(n) rep(NA_character_, n)))
  if (!is.list(labels) || length(labels) != length(sizes)) {
    stop("`labels` must be a list with one vector per domain.", call. = FALSE)
  }
  lapply(seq_along(labels), function(i) {
    if (length(labels[[i]]) != sizes[[i]]) {
      stop("Domain `", domain_names[[i]], "`: labels must have one value per row.",
           call. = FALSE)
    }
    as.character(labels[[i]])
  })
}

.rapid_pad_embedding <- function(Z, width) {
  Z <- as.matrix(Z)
  if (ncol(Z) > width) Z <- Z[, seq_len(width), drop = FALSE]
  if (ncol(Z) < width) {
    Z <- cbind(Z, matrix(0, nrow(Z), width - ncol(Z)))
  }
  colnames(Z) <- paste0("rapid_component_", seq_len(width))
  Z
}

.rapid_align_named_views <- function(views, prefix) {
  prepared <- lapply(seq_along(views), function(i) {
    view <- views[[i]]
    if (is.null(view)) return(NULL)
    view <- as.matrix(view)
    if (!ncol(view)) return(NULL)
    names <- colnames(view)
    if (is.null(names) || any(!nzchar(names)) || anyDuplicated(names)) {
      names <- paste0(prefix, "_", seq_len(ncol(view)))
    } else {
      names <- paste0(prefix, "__", names)
    }
    colnames(view) <- names
    view
  })
  columns <- unique(unlist(lapply(prepared, colnames), use.names = FALSE))
  if (!length(columns)) return(NULL)
  columns <- sort(columns, method = "radix")
  lapply(seq_along(prepared), function(i) {
    n <- if (is.null(prepared[[i]])) nrow(views[[i]]) else nrow(prepared[[i]])
    out <- matrix(0, n, length(columns), dimnames = list(NULL, columns))
    if (!is.null(prepared[[i]])) {
      out[, colnames(prepared[[i]])] <- prepared[[i]]
    }
    out
  })
}

.rapid_domain_preprocess <- function(X, labels, relations, positions, attributes,
                                     ncomp, control, domain_name) {
  relation_control <- control$relation
  feature_width <- relation_control$feature_sketch_dim
  if (is.null(feature_width)) {
    feature_width <- min(ncol(X), max(128L, 4L * ncomp))
  }
  relation_state <- tryCatch(
    .rapid_prepare_relations(
      X,
      relations = relations,
      positions = positions,
      attributes = attributes,
      build_feature = relation_control$build_feature,
      feature_k = relation_control$feature_k,
      feature_sketch_dim = feature_width,
      spatial_mode = relation_control$spatial_mode,
      spatial_k = relation_control$spatial_k,
      spatial_radius = relation_control$spatial_radius,
      spatial_sigma = relation_control$spatial_sigma,
      attribute_mode = relation_control$attribute_mode,
      attribute_k = relation_control$attribute_k,
      attribute_hash_dim = relation_control$attribute_hash_dim,
      normalization = relation_control$normalization,
      weight_mode = relation_control$weight_mode,
      degree_cap = relation_control$degree_cap,
      custom_symmetrize = relation_control$custom_symmetrize,
      dense_max_n = relation_control$dense_max_n,
      seed = control$seed
    ),
    error = function(e) {
      stop("Domain `", domain_name, "` relation preparation: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  diffusion_control <- unclass(control$diffusion)
  diffusion_control$output_dim <- ncomp
  diffusion_control$sketch_dim <- max(
    diffusion_control$sketch_dim, feature_width, ncomp
  )
  diffusion_control$seed <- control$seed
  diffusion <- tryCatch(
    .rapid_diffusion_encode(
      X,
      relation_state,
      control = diffusion_control,
      labels = labels
    ),
    error = function(e) {
      stop("Domain `", domain_name, "` diffusion encoding: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  list(
    relations = relation_state,
    diffusion = diffusion,
    score = .rapid_pad_embedding(diffusion$embedding, ncomp),
    structure = diffusion$signatures,
    position = .rapid_named_view(relation_state$views$position),
    attribute = .rapid_named_view(relation_state$views$attribute),
    feature_width = as.integer(feature_width)
  )
}

.rapid_identity_map <- function(width) {
  list(
    coef = diag(width),
    intercept = rep(0, width),
    type = "identity",
    status = "initial",
    effective_rows = 0L
  )
}

.rapid_apply_latent_map <- function(Z, map) {
  sweep(Z %*% map$coef, 2L, map$intercept, "+")
}

.rapid_fit_latent_map <- function(Z, target, weight, type, ridge) {
  d <- ncol(Z)
  weight <- as.numeric(weight)
  positive <- is.finite(weight) & weight > max(weight, na.rm = TRUE) * 1e-10
  if (identical(type, "identity") || sum(positive) < 2L) {
    out <- .rapid_identity_map(d)
    out$status <- if (identical(type, "identity")) "disabled" else "insufficient_mass"
    out$effective_rows <- as.integer(sum(positive))
    return(out)
  }
  Zp <- Z[positive, , drop = FALSE]
  Tp <- target[positive, , drop = FALSE]
  wp <- weight[positive]
  wp <- wp / mean(wp)

  if (identical(type, "orthogonal")) {
    total <- sum(wp)
    center_z <- colSums(Zp * wp) / total
    center_t <- colSums(Tp * wp) / total
    Zc <- sweep(Zp, 2L, center_z, "-") * sqrt(wp)
    Tc <- sweep(Tp, 2L, center_t, "-") * sqrt(wp)
    decomposition <- svd(crossprod(Zc, Tc), nu = d, nv = d)
    coefficient <- decomposition$u %*% t(decomposition$v)
    intercept <- center_t - as.numeric(center_z %*% coefficient)
    return(list(
      coef = coefficient,
      intercept = intercept,
      type = type,
      status = "ok",
      effective_rows = as.integer(sum(positive))
    ))
  }

  design <- cbind(Zp, intercept = 1)
  weighted_design <- design * sqrt(wp)
  weighted_target <- Tp * sqrt(wp)
  penalty <- diag(c(rep(ridge, d), 0), d + 1L)
  system <- crossprod(weighted_design) + penalty
  rhs <- crossprod(weighted_design, weighted_target)
  coefficient <- tryCatch(
    solve(system, rhs),
    error = function(...) qr.solve(system, rhs, tol = 1e-10)
  )
  list(
    coef = coefficient[seq_len(d), , drop = FALSE],
    intercept = as.numeric(coefficient[d + 1L, ]),
    type = type,
    status = "ok",
    effective_rows = as.integer(sum(positive))
  )
}

.rapid_blend_map <- function(old, new, step) {
  list(
    coef = (1 - step) * old$coef + step * new$coef,
    intercept = (1 - step) * old$intercept + step * new$intercept,
    type = new$type,
    status = new$status,
    effective_rows = new$effective_rows
  )
}

.rapid_barycentric_target <- function(uot_fit, prototypes, fallback) {
  mass <- uot_fit$row_mass
  target <- as.matrix(uot_fit$coupling %*% prototypes$embedding)
  positive <- mass > max(mass, 0) * 1e-12
  if (any(positive)) {
    target[positive, ] <- target[positive, , drop = FALSE] / mass[positive]
  }
  if (any(!positive)) target[!positive, ] <- fallback[!positive, , drop = FALSE]
  list(target = target, weight = mass)
}

.rapid_blend_prototypes <- function(old, new, step) {
  out <- new
  for (name in c("embedding", "structure", "attribute")) {
    out[[name]] <- (1 - step) * old[[name]] + step * new[[name]]
  }
  out$class_prob <- (1 - step) * old$class_prob + step * new$class_prob
  out$class_prob <- out$class_prob / pmax(rowSums(out$class_prob), 1e-15)
  out$mass <- (1 - step) * old$mass + step * new$mass
  out$active <- old$active | new$active
  out$prior <- old$prior
  out$prototype_class <- old$prototype_class
  out$fixed_class <- old$fixed_class
  out$class_levels <- old$class_levels
  out
}

.rapid_generalized_kl <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y) || any(x < 0) || any(y < 0)) return(Inf)
  positive <- x > 0
  if (any(positive & y <= 0)) return(Inf)
  sum(ifelse(positive, x * log(x / pmax(y, .Machine$double.xmin)), 0) - x + y)
}

.rapid_uot_objective <- function(fit, prototypes, control) {
  Q <- fit$coupling
  sm <- Matrix::summary(Q)
  beta_global <- numeric(nrow(prototypes$embedding))
  beta_global[fit$prototype_index] <- fit$beta
  marginal_source <- .rapid_generalized_kl(fit$row_mass, fit$alpha)
  marginal_target <- .rapid_generalized_kl(fit$prototype_mass, beta_global)
  if (nrow(sm)) {
    reference <- fit$alpha[sm$i] * beta_global[sm$j]
    entropy <- sum(sm$x * log(sm$x / pmax(reference, .Machine$double.xmin))) -
      sum(sm$x) + sum(fit$alpha) * sum(beta_global)
  } else {
    entropy <- sum(fit$alpha) * sum(beta_global)
  }
  as.numeric(
    fit$diagnostics$transport_cost +
      control$rho_source * marginal_source +
      control$rho_target * marginal_target +
      control$epsilon * entropy
  )
}

.rapid_prediction_state <- function(couplings, prototypes, labels) {
  class_prior <- as.numeric(crossprod(prototypes$prior, prototypes$class_prob))
  names(class_prior) <- prototypes$class_levels
  probabilities <- predictions <- confidence <- vector("list", length(couplings))
  supervised_loss <- numeric(length(couplings))
  for (m in seq_along(couplings)) {
    Q <- couplings[[m]]$coupling
    mass <- couplings[[m]]$row_mass
    probability <- as.matrix(Q %*% prototypes$class_prob)
    positive <- mass > 0
    if (any(positive)) {
      probability[positive, ] <- probability[positive, , drop = FALSE] / mass[positive]
    }
    if (any(!positive)) probability[!positive, ] <- rep(class_prior, each = sum(!positive))
    probability <- probability / pmax(rowSums(probability), 1e-15)
    colnames(probability) <- prototypes$class_levels
    choice <- max.col(probability, ties.method = "first")
    predictions[[m]] <- prototypes$class_levels[choice]
    confidence[[m]] <- probability[cbind(seq_len(nrow(probability)), choice)]
    probabilities[[m]] <- probability

    observed <- which(!is.na(labels[[m]]))
    if (length(observed)) {
      class_index <- match(labels[[m]][observed], prototypes$class_levels)
      valid <- !is.na(class_index)
      observed <- observed[valid]
      class_index <- class_index[valid]
      if (length(observed)) {
        counts <- table(labels[[m]][observed])
        weight <- 1 / as.numeric(counts[labels[[m]][observed]])
        weight <- weight / mean(weight)
        truth_probability <- probability[cbind(observed, class_index)]
        supervised_loss[[m]] <- mean(
          weight * -log(pmax(truth_probability, 1e-15))
        )
      }
    }
  }
  list(
    probabilities = probabilities,
    predictions = predictions,
    confidence = confidence,
    supervised_loss = mean(supervised_loss)
  )
}

.rapid_evaluate_state <- function(scores, prototypes, relation_states,
                                  structures, attributes, labels,
                                  diffusion_fits, control, domain_names,
                                  warm_start = NULL) {
  uot_fits <- vector("list", length(scores))
  for (m in seq_along(scores)) {
    uot_fits[[m]] <- tryCatch(
      .rapid_uot_solve(
        scores[[m]],
        prototypes,
        structure = if (is.null(structures)) NULL else structures[[m]],
        attributes = if (is.null(attributes)) NULL else attributes[[m]],
        labels = labels[[m]],
        control = control$uot,
        warm_start = if (is.null(warm_start)) NULL else warm_start[[m]]
      ),
      error = function(e) {
        stop("Domain `", domain_names[[m]], "` transport: ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  names(uot_fits) <- domain_names
  has_relations <- any(vapply(
    relation_states, function(x) length(x$relations) > 0L, logical(1)
  ))
  relation_fit <- if (has_relations) {
    .rapid_relmatch_fit(
      relation_states,
      uot_fits,
      domain_relation_weights = stats::setNames(
        lapply(diffusion_fits, `[[`, "relation_weights"), domain_names
      ),
      control = control$relmatch
    )
  } else {
    NULL
  }
  prediction <- .rapid_prediction_state(uot_fits, prototypes, labels)
  uot_objective <- mean(vapply(
    uot_fits, .rapid_uot_objective, numeric(1),
    prototypes = prototypes, control = control$uot
  ))
  relation_objective <- if (is.null(relation_fit)) 0 else relation_fit$total_loss
  components <- c(
    uot = uot_objective,
    structure = relation_objective,
    supervised = prediction$supervised_loss
  )
  objective <- control$uot_weight * components[["uot"]] +
    control$structure_weight * components[["structure"]] +
    control$supervised_weight * components[["supervised"]]
  if (!is.finite(objective)) stop("RAPID-MA objective is non-finite.", call. = FALSE)
  list(
    objective = as.numeric(objective),
    components = components,
    uot = uot_fits,
    relation = relation_fit,
    prediction = prediction
  )
}

.rapid_history_row <- function(iteration, evaluation, accepted, step, backtracks,
                               status) {
  data.frame(
    iteration = as.integer(iteration),
    objective = evaluation$objective,
    uot = evaluation$components[["uot"]],
    structure = evaluation$components[["structure"]],
    supervised = evaluation$components[["supervised"]],
    accepted = accepted,
    step = step,
    backtracks = as.integer(backtracks),
    status = status,
    stringsAsFactors = FALSE
  )
}

.rapid_ma_fit <- function(data, labels, relations, positions, attributes,
                          ncomp, control, domain_names, call) {
  control <- .rapid_resolve_control(control)
  ncomp <- .rapid_int_scalar(ncomp, "ncomp")
  X <- .rapid_validate_domain_data(data, domain_names)
  sizes <- vapply(X, nrow, integer(1))
  names(sizes) <- domain_names
  labels <- .rapid_validate_domain_labels(labels, sizes, domain_names)
  relations <- .rapid_named_domain_argument(relations, domain_names, "relations")
  positions <- .rapid_named_domain_argument(positions, domain_names, "positions")
  attributes <- .rapid_named_domain_argument(attributes, domain_names, "attributes")

  preprocessed <- vector("list", length(X))
  for (m in seq_along(X)) {
    if (isTRUE(control$verbose)) {
      message("rapid_ma: preprocessing domain `", domain_names[[m]], "`")
    }
    preprocessed[[m]] <- .rapid_domain_preprocess(
      X[[m]], labels[[m]], relations[[m]], positions[[m]], attributes[[m]],
      ncomp, control, domain_names[[m]]
    )
  }
  names(preprocessed) <- domain_names
  relation_states <- lapply(preprocessed, `[[`, "relations")
  diffusion_fits <- lapply(preprocessed, `[[`, "diffusion")
  base_scores <- lapply(preprocessed, `[[`, "score")
  names(base_scores) <- domain_names
  oos_features <- lapply(seq_along(X), function(m) {
    existing <- relation_states[[m]]$views$feature
    if (!is.null(existing) && !is.null(existing$view)) return(existing)
    .rapid_feature_view(
      X[[m]],
      max_dim = preprocessed[[m]]$feature_width,
      seed = control$seed
    )
  })
  names(oos_features) <- domain_names

  structural_raw <- lapply(preprocessed, function(x) {
    signature <- x$structure
    position <- x$position
    if (!is.null(position) && ncol(position)) {
      colnames(position) <- paste0("coordinate_", seq_len(ncol(position)))
      signature <- cbind(signature, position)
    }
    signature
  })
  structures <- .rapid_align_named_views(structural_raw, "structure")
  attribute_raw <- lapply(seq_along(preprocessed), function(i) {
    view <- preprocessed[[i]]$attribute
    if (is.null(view)) matrix(numeric(0), nrow(base_scores[[i]]), 0L) else view
  })
  attributes_aligned <- .rapid_align_named_views(attribute_raw, "attribute")

  prototypes <- .rapid_initialize_prototypes(
    base_scores,
    labels,
    structures = structures,
    attributes = attributes_aligned,
    control = control$prototype
  )
  maps <- lapply(base_scores, function(x) .rapid_identity_map(ncol(x)))
  names(maps) <- domain_names
  scores <- base_scores
  evaluation <- .rapid_evaluate_state(
    scores, prototypes, relation_states, structures, attributes_aligned,
    labels, diffusion_fits, control, domain_names
  )
  history <- list(.rapid_history_row(
    0L, evaluation, TRUE, 0, 0L, "initial"
  ))
  status <- "max_iter"
  accepted_iterations <- 0L

  for (iteration in seq_len(control$max_iter)) {
    current_Q <- lapply(evaluation$uot, `[[`, "coupling")
    prototype_target <- .rapid_update_prototypes(
      prototypes,
      current_Q,
      scores,
      labels,
      structures = structures,
      attributes = attributes_aligned,
      control = control$prototype
    )
    fitted_maps <- vector("list", length(scores))
    for (m in seq_along(scores)) {
      target <- .rapid_barycentric_target(
        evaluation$uot[[m]], prototype_target, scores[[m]]
      )
      fitted_maps[[m]] <- .rapid_fit_latent_map(
        base_scores[[m]], target$target, target$weight,
        type = control$map_type, ridge = control$map_ridge
      )
    }
    names(fitted_maps) <- domain_names

    accepted <- FALSE
    candidate_error <- NULL
    accepted_step <- 0
    accepted_backtracks <- 0L
    for (backtrack in 0:control$max_backtracks) {
      damping <- control$backtrack_factor^backtrack
      map_step <- control$map_step * damping
      prototype_step <- control$prototype_step * damping
      candidate_maps <- Map(
        .rapid_blend_map, maps, fitted_maps,
        MoreArgs = list(step = map_step)
      )
      candidate_scores <- Map(.rapid_apply_latent_map, base_scores, candidate_maps)
      candidate_prototypes <- .rapid_update_prototypes(
        prototype_target,
        current_Q,
        candidate_scores,
        labels,
        structures = structures,
        attributes = attributes_aligned,
        control = control$prototype
      )
      candidate_prototypes <- .rapid_blend_prototypes(
        prototypes, candidate_prototypes, prototype_step
      )
      candidate_evaluation <- tryCatch(
        .rapid_evaluate_state(
          candidate_scores, candidate_prototypes, relation_states,
          structures, attributes_aligned, labels, diffusion_fits,
          control, domain_names, warm_start = evaluation$uot
        ),
        error = function(e) e
      )
      if (inherits(candidate_evaluation, "error")) {
        candidate_error <- conditionMessage(candidate_evaluation)
        next
      }
      if (candidate_evaluation$objective <=
          evaluation$objective + control$objective_tolerance) {
        previous_objective <- evaluation$objective
        scores <- candidate_scores
        maps <- candidate_maps
        prototypes <- candidate_prototypes
        evaluation <- candidate_evaluation
        accepted <- TRUE
        accepted_step <- map_step
        accepted_backtracks <- backtrack
        accepted_iterations <- accepted_iterations + 1L
        improvement <- previous_objective - evaluation$objective
        relative_improvement <- improvement / max(abs(previous_objective), 1e-12)
        history[[length(history) + 1L]] <- .rapid_history_row(
          iteration, evaluation, TRUE, accepted_step, accepted_backtracks,
          "accepted"
        )
        if (isTRUE(control$verbose)) {
          message(
            "rapid_ma: iteration ", iteration,
            " objective=", format(evaluation$objective, digits = 7),
            " relative improvement=", format(relative_improvement, digits = 4),
            " backtracks=", backtrack
          )
        }
        if (iteration >= control$min_iter && relative_improvement <= control$tol) {
          status <- "converged"
        }
        break
      }
    }
    if (!accepted) {
      status <- "stalled"
      history[[length(history) + 1L]] <- .rapid_history_row(
        iteration, evaluation, FALSE, 0, control$max_backtracks,
        if (is.null(candidate_error)) "rolled_back" else "candidate_error"
      )
      break
    }
    if (identical(status, "converged")) break
  }

  history <- do.call(rbind, history)
  rownames(history) <- NULL
  couplings <- lapply(evaluation$uot, `[[`, "coupling")
  names(couplings) <- domain_names
  probabilities <- evaluation$prediction$probabilities
  predictions <- evaluation$prediction$predictions
  confidence <- evaluation$prediction$confidence
  names(probabilities) <- names(predictions) <- names(confidence) <- domain_names
  dropped_mass <- vapply(
    evaluation$uot, function(x) x$diagnostics$dropped_fraction, numeric(1)
  )

  structure(
    list(
      call = call,
      scores = scores,
      s = do.call(rbind, scores),
      prototypes = prototypes,
      couplings = couplings,
      transport = evaluation$uot,
      relation_match = evaluation$relation,
      relation_barycenters = if (is.null(evaluation$relation)) {
        list()
      } else {
        evaluation$relation$barycenters
      },
      objective = evaluation$objective,
      objective_components = evaluation$components,
      objective_history = history,
      predictions = predictions,
      prediction_probabilities = probabilities,
      prediction_confidence = confidence,
      dropped_mass = dropped_mass,
      base_scores = base_scores,
      domain_maps = maps,
      relations = relation_states,
      diffusion = diffusion_fits,
      preprocessing = list(
        feature_width = vapply(preprocessed, `[[`, integer(1), "feature_width"),
        structures = structures,
        attributes = attributes_aligned,
        oos_feature_views = lapply(oos_features, `[[`, "view"),
        oos_feature_metadata = lapply(oos_features, `[[`, "metadata"),
        input_dims = stats::setNames(vapply(X, ncol, integer(1)), domain_names),
        embedding_only = TRUE,
        loadings_available = FALSE
      ),
      labels = labels,
      domain_names = domain_names,
      domain_sizes = sizes,
      ncomp = as.integer(ncomp),
      control = control,
      seed = control$seed,
      convergence = list(
        status = status,
        accepted_iterations = accepted_iterations,
        attempted_iterations = max(history$iteration),
        monotone_accepted = all(diff(history$objective[history$accepted]) <=
                                  control$objective_tolerance),
        last_candidate_error = candidate_error
      )
    ),
    class = c("rapid_ma", "manifold_alignment")
  )
}

#' Scalable structure-aware manifold alignment
#'
#' @param data A list of numeric domain matrices or a `hyperdesign`.
#' @param ... Method-specific arguments.
#' @return A `rapid_ma` embedding result with sparse couplings and diagnostics.
#' @export
rapid_ma <- function(data, ...) {
  UseMethod("rapid_ma")
}

#' @rdname rapid_ma
#' @param labels Optional list of partially observed label vectors.
#' @param relations Optional list of custom named sparse relations per domain.
#' @param positions Optional list of spatial coordinate matrices per domain.
#' @param attributes Optional list of mixed attribute data.frames per domain.
#' @param ncomp Shared embedding dimension.
#' @param control Settings from [rapid_ma_control()].
#' @export
rapid_ma.list <- function(
  data,
  labels = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  ...
) {
  domain_names <- names(data)
  if (is.null(domain_names) || any(!nzchar(domain_names))) {
    domain_names <- paste0("domain", seq_along(data))
  }
  if (!is.null(labels) && !is.null(names(labels)) &&
      all(domain_names %in% names(labels))) {
    labels <- unname(labels[domain_names])
  }
  .rapid_ma_fit(
    data, labels, relations, positions, attributes, ncomp, control,
    domain_names, match.call()
  )
}

#' @rdname rapid_ma
#' @param y Optional unquoted or string label column in each domain design.
#' @method rapid_ma hyperdesign
#' @export
rapid_ma.hyperdesign <- function(
  data,
  y = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  ...
) {
  resolved <- resolve_hyperdesign(data)
  y_quo <- rlang::enquo(y)
  labels <- .ssma_extract_labels(resolved$domains, y_quo)
  matrices <- lapply(resolved$domains, `[[`, "x")
  .rapid_ma_fit(
    matrices, labels, relations, positions, attributes, ncomp, control,
    resolved$domain_names, match.call()
  )
}

#' @export
print.rapid_ma <- function(x, ...) {
  cat("RAPID-MA structure-aware alignment\n")
  cat("  domains: ", length(x$domain_names), " (",
      paste(x$domain_sizes, collapse = ", "), ")\n", sep = "")
  cat("  components: ", x$ncomp, "\n", sep = "")
  cat("  objective: ", format(x$objective, digits = 7), "\n", sep = "")
  cat("  status: ", x$convergence$status, "\n", sep = "")
  invisible(x)
}
