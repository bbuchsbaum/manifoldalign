# RAPID-MA adapter, OOS interpolation, and bounded fine matching ------------

#' Construct a RAPID-MA aligner descriptor
#'
#' @return An object of class `c("rapid_ma_aligner", "aligner")`.
#' @export
rapid_ma_aligner <- function() {
  new_aligner("rapid_ma", group = "perm", supports_multi = TRUE)
}

.rapid_adapter_domains <- function(x, default = c("i", "j")) {
  if (is.null(x)) return(stats::setNames(vector("list", length(default)), default))
  if (!is.list(x) || length(x) != length(default)) {
    stop("Pairwise domain arguments must be lists of length two.", call. = FALSE)
  }
  if (!is.null(names(x)) && all(default %in% names(x))) x <- x[default]
  stats::setNames(unname(x), default)
}

.rapid_adapter_side <- function(fit, side) {
  if (length(side) != 1L || is.na(side)) {
    stop("`side` must identify exactly one fitted domain.", call. = FALSE)
  }
  if (is.character(side)) {
    if (side %in% c("i", "j")) side <- if (side == "i") 1L else 2L
    else side <- match(side, fit$domain_names)
  }
  if (!is.numeric(side) || !is.finite(side) ||
      abs(side - round(side)) > 1e-8) {
    stop("`side` must be 'i', 'j', a domain name, or a domain index.",
         call. = FALSE)
  }
  side <- as.integer(side)
  if (side < 1L || side > length(fit$domain_names)) {
    stop("`side` does not identify a fitted domain.", call. = FALSE)
  }
  side
}

#' @export
fit_many.rapid_ma_aligner <- function(algo, domains, k = NULL, ncomp = NULL, ...) {
  if (is.null(ncomp)) ncomp <- if (is.null(k)) 16L else k
  rapid_ma(domains, ncomp = ncomp, ...)
}

#' Fit RAPID-MA to a pair of domains
#'
#' @param algo A `rapid_ma_aligner` descriptor.
#' @param X_i,X_j Numeric domain matrices.
#' @param links Optional one-to-one anchor pairs or paired anchor-ID vectors.
#' @param labels,relations,positions,attributes Optional two-domain inputs.
#' @param ncomp Shared embedding dimension.
#' @param control Settings from [rapid_ma_control()].
#' @param fine_match Whether to compute bounded node-level matches.
#' @param match_control Named arguments forwarded to [rapid_ma_match()].
#' @param ... Unused.
#' @return A `rapid_ma_pair_fit` containing the shared core fit and fine match.
#' @export
fit_pair.rapid_ma_aligner <- function(
  algo,
  X_i,
  X_j,
  links = NULL,
  labels = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  fine_match = TRUE,
  match_control = list(),
  ...
) {
  domains <- list(i = as.matrix(X_i), j = as.matrix(X_j))
  if (!is.null(labels)) labels <- .rapid_adapter_domains(labels)
  relations <- .rapid_adapter_domains(relations)
  positions <- .rapid_adapter_domains(positions)
  attributes <- .rapid_adapter_domains(attributes)
  core <- rapid_ma(
    domains,
    labels = labels,
    relations = relations,
    positions = positions,
    attributes = attributes,
    ncomp = ncomp,
    control = control
  )
  matching <- NULL
  if (isTRUE(fine_match)) {
    if (!is.list(match_control) ||
        (!is.null(names(match_control)) && any(!nzchar(names(match_control)))) ||
        (length(match_control) && is.null(names(match_control)))) {
      stop("`match_control` must be a named list.", call. = FALSE)
    }
    matching <- do.call(
      rapid_ma_match,
      c(list(fit = core, from = "i", to = "j", anchors = links), match_control)
    )
  }
  structure(
    list(
      core = core,
      matching = matching,
      objective = core$objective,
      n1 = nrow(X_i),
      n2 = nrow(X_j),
      k = core$ncomp,
      anchors = links
    ),
    class = "rapid_ma_pair_fit"
  )
}

.rapid_apply_oos_feature_view <- function(newX, metadata) {
  .rapid_validate_numeric_matrix(newX, "newX")
  X <- as.matrix(newX)
  if (any(!is.finite(X))) stop("`newX` contains non-finite values.", call. = FALSE)
  if (ncol(X) != metadata$original_dim) {
    stop("`newX` has ", ncol(X), " columns but this domain expects ",
         metadata$original_dim, ".", call. = FALSE)
  }
  if (isTRUE(metadata$projected)) {
    projection <- Matrix::sparseMatrix(
      i = seq_len(metadata$original_dim),
      j = metadata$buckets,
      x = metadata$signs,
      dims = c(metadata$original_dim, metadata$output_dim)
    )
    X <- as.matrix(X %*% projection)
  }
  X <- sweep(X, 2L, metadata$center, "-")
  X <- sweep(X, 2L, metadata$scale, "/")
  X[!is.finite(X)] <- 0
  X
}

.rapid_apply_oos_position_view <- function(new_positions, metadata, n_rows) {
  .rapid_validate_numeric_matrix(
    new_positions, "new_positions", n_rows = n_rows
  )
  position <- as.matrix(new_positions)
  input_dim <- metadata$input_dim
  if (is.null(input_dim)) input_dim <- length(metadata$center)
  if (ncol(position) != input_dim) {
    stop(
      "`new_positions` has ", ncol(position),
      " columns but this domain expects ", input_dim, ".",
      call. = FALSE
    )
  }
  valid <- rowSums(!is.finite(position)) == 0L
  if (identical(metadata$mode, "relative") && any(valid)) {
    position[valid, ] <- sweep(
      position[valid, , drop = FALSE], 2L, metadata$center, "-"
    ) / metadata$scale
  }
  position[!is.finite(position)] <- 0
  if (isTRUE(metadata$missing_indicator)) {
    position <- cbind(position, position_missing = as.numeric(!valid))
  }
  position
}

.rapid_apply_oos_attribute_view <- function(new_attributes, metadata, n_rows) {
  .rapid_validate_attributes(new_attributes, n_rows)
  expected <- metadata$column_names
  if (is.null(expected)) {
    stop(
      "The fitted attribute encoder predates OOS attribute support; refit RAPID-MA.",
      call. = FALSE
    )
  }
  supplied <- names(new_attributes)
  if (!identical(supplied, expected)) {
    can_reorder <- !anyDuplicated(supplied) && !anyDuplicated(expected) &&
      length(supplied) == length(expected) && setequal(supplied, expected)
    if (!can_reorder) {
      stop(
        "`new_attributes` must contain the fitted columns: ",
        paste(expected, collapse = ", "), ".",
        call. = FALSE
      )
    }
    new_attributes <- new_attributes[, expected, drop = FALSE]
  }

  hash_dim <- metadata$hash_dim
  seed <- metadata$seed
  n <- nrow(new_attributes)
  numeric_hash <- matrix(0, n, hash_dim)
  categorical_hash <- matrix(0, n, hash_dim)
  missing_hash <- matrix(0, n, hash_dim)
  encoding_names <- metadata$encoding_names
  if (is.null(encoding_names)) {
    encoding_names <- vapply(seq_along(expected), function(j) {
      if (nzchar(expected[[j]])) expected[[j]] else paste0("attribute_", j)
    }, character(1))
  }

  for (j in seq_along(new_attributes)) {
    name <- encoding_names[[j]]
    value <- new_attributes[[j]]
    missing <- is.na(value)
    numeric_value <- is.numeric(value) && !is.factor(value) && !is.logical(value)
    if (identical(metadata$column_types[[j]], "numeric") && !numeric_value) {
      stop("`new_attributes$", expected[[j]], "` must be numeric.", call. = FALSE)
    }
    if (identical(metadata$column_types[[j]], "categorical") && numeric_value) {
      stop(
        "`new_attributes$", expected[[j]], "` must be categorical or logical.",
        call. = FALSE
      )
    }
    if (numeric_value) {
      x <- as.numeric(value)
      missing <- missing | !is.finite(x)
      stats <- metadata$numeric_stats[[name]]
      z <- (x - stats$center) / stats$scale
      z[missing | !is.finite(z)] <- 0
      z <- pmax(pmin(z, 10), -10)
      key <- paste0("numeric:", name)
      bucket <- .rapid_hash_bucket(key, hash_dim, seed)
      sign <- .rapid_hash_sign(key, seed)
      numeric_hash[, bucket] <- numeric_hash[, bucket] + sign * z
    } else {
      x <- as.character(value)
      rows <- which(!missing)
      if (length(rows)) {
        values <- enc2utf8(x[rows])
        unique_values <- unique(values)
        keys <- paste0("categorical:", name, "=", unique_values)
        buckets <- vapply(
          keys, .rapid_hash_bucket, integer(1), width = hash_dim, seed = seed
        )
        signs <- vapply(keys, .rapid_hash_sign, numeric(1), seed = seed)
        value_index <- match(values, unique_values)
        cells <- cbind(rows, buckets[value_index])
        categorical_hash[cells] <-
          categorical_hash[cells] + signs[value_index]
      }
    }
    if (any(missing)) {
      key <- paste0("missing:", name)
      bucket <- .rapid_hash_bucket(key, hash_dim, seed)
      sign <- .rapid_hash_sign(key, seed)
      rows <- which(missing)
      cells <- cbind(rows, rep.int(bucket, length(rows)))
      missing_hash[cells] <- missing_hash[cells] + sign
    }
  }

  raw <- cbind(numeric_hash, categorical_hash, missing_hash)
  if (is.null(metadata$center) || is.null(metadata$scale)) {
    stop(
      "The fitted attribute encoder predates OOS attribute support; refit RAPID-MA.",
      call. = FALSE
    )
  }
  encoded <- sweep(raw, 2L, metadata$center, "-")
  encoded <- sweep(encoded, 2L, metadata$scale, "/")
  encoded[!is.finite(encoded)] <- 0
  encoded <- encoded[, metadata$active_columns, drop = FALSE]
  if (ncol(encoded) != metadata$output_dim) {
    stop("The fitted OOS attribute encoder is internally inconsistent.", call. = FALSE)
  }
  encoded
}

.rapid_oos_view_weights <- function(view_weights) {
  allowed <- c("feature", "position", "attribute")
  if (is.null(view_weights)) view_weights <- c(feature = 1)
  if (!is.numeric(view_weights) || !length(view_weights) ||
      is.null(names(view_weights)) || any(!nzchar(names(view_weights))) ||
      anyDuplicated(names(view_weights)) ||
      any(!is.finite(view_weights)) || any(view_weights < 0)) {
    stop(
      "`view_weights` must be a named, finite, nonnegative numeric vector.",
      call. = FALSE
    )
  }
  unknown <- setdiff(names(view_weights), allowed)
  if (length(unknown)) {
    stop(
      "Unknown OOS views: ", paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  output <- stats::setNames(numeric(length(allowed)), allowed)
  output[names(view_weights)] <- view_weights
  if (sum(output) <= 0) {
    stop("At least one `view_weights` value must be positive.", call. = FALSE)
  }
  output / sum(output)
}

.rapid_oos_weights <- function(training, query, k, zero_tolerance) {
  k <- min(.rapid_int_scalar(k, "k"), nrow(training))
  zero_tolerance <- .rapid_number_scalar(
    zero_tolerance, "zero_tolerance", lower = 0
  )
  neighbours <- RANN::nn2(data = training, query = query, k = k)
  distance <- neighbours$nn.dists
  index <- neighbours$nn.idx
  weight <- matrix(0, nrow(distance), ncol(distance))
  for (row in seq_len(nrow(distance))) {
    exact <- is.finite(distance[row, ]) & distance[row, ] <= zero_tolerance
    if (any(exact)) {
      weight[row, exact] <- 1 / sum(exact)
    } else {
      finite <- distance[row, is.finite(distance[row, ])]
      scale <- if (length(finite)) stats::median(finite) else 1
      if (!is.finite(scale) || scale <= 1e-12) scale <- max(finite, 1)
      value <- exp(pmax(-0.5 * (distance[row, ] / scale)^2, -700))
      value[!is.finite(value)] <- 0
      if (sum(value) == 0) value[[which.min(distance[row, ])]] <- 1
      weight[row, ] <- value / sum(value)
    }
  }
  interpolation <- Matrix::sparseMatrix(
    i = rep(seq_len(nrow(index)), each = ncol(index)),
    j = as.integer(c(t(index))),
    x = as.numeric(c(t(weight))),
    dims = c(nrow(index), nrow(training))
  )
  list(
    interpolation = interpolation,
    index = index,
    distance = distance,
    weight = weight
  )
}

#' Out-of-sample RAPID-MA embedding and prediction
#'
#' @param fit_or_transform A fitted `rapid_ma` or `rapid_ma_pair_fit`.
#' @param newX New rows in the original feature space of `side`.
#' @param side Domain name, index, or pairwise alias `"i"`/`"j"`.
#' @param type Return embedding coordinates, class probabilities, class labels,
#'   or all diagnostics.
#' @param k Number of bounded interpolation neighbours.
#' @param zero_tolerance Distance treated as an exact in-sample reapplication.
#' @param new_positions Optional position matrix for structure-aware OOS
#'   interpolation.
#' @param new_attributes Optional attribute data frame for structure-aware OOS
#'   interpolation.
#' @param view_weights Named nonnegative weights over `feature`, `position`, and
#'   `attribute`. The default is feature-only for backward compatibility.
#' @param ... Unused.
#' @return A matrix/vector, or a diagnostic list when `type = "all"`.
#' @export
oos_predict.rapid_ma <- function(
  fit_or_transform,
  newX,
  side,
  type = c("embedding", "probabilities", "class", "all"),
  k = 8L,
  zero_tolerance = sqrt(.Machine$double.eps),
  new_positions = NULL,
  new_attributes = NULL,
  view_weights = c(feature = 1),
  ...
) {
  type <- match.arg(type)
  m <- .rapid_adapter_side(fit_or_transform, side)
  metadata <- fit_or_transform$preprocessing$oos_feature_metadata[[m]]
  training <- fit_or_transform$preprocessing$oos_feature_views[[m]]
  if (is.null(metadata) || is.null(training)) {
    stop("The fit does not contain RAPID-MA OOS preprocessing state.",
         call. = FALSE)
  }
  feature_query <- .rapid_apply_oos_feature_view(newX, metadata)
  feature_neighbours <- .rapid_oos_weights(
    training, feature_query, k, zero_tolerance
  )
  weights <- .rapid_oos_view_weights(view_weights)
  view_state <- list()
  if (weights[["feature"]] > 0) view_state$feature <- feature_neighbours

  retained_views <- fit_or_transform$relations[[m]]$views
  if (weights[["position"]] > 0) {
    position_state <- retained_views$position
    if (is.null(new_positions)) {
      stop(
        "Positive position weight requires `new_positions`.", call. = FALSE
      )
    }
    if (is.null(position_state$view) || !ncol(position_state$view)) {
      stop("This domain has no retained position OOS view.", call. = FALSE)
    }
    position_query <- .rapid_apply_oos_position_view(
      new_positions, position_state$metadata, nrow(feature_query)
    )
    if (ncol(position_query) != ncol(position_state$view)) {
      stop("The fitted OOS position encoder is internally inconsistent.",
           call. = FALSE)
    }
    view_state$position <- .rapid_oos_weights(
      position_state$view, position_query, k, zero_tolerance
    )
  }
  if (weights[["attribute"]] > 0) {
    attribute_state <- retained_views$attribute
    if (is.null(new_attributes)) {
      stop(
        "Positive attribute weight requires `new_attributes`.", call. = FALSE
      )
    }
    if (is.null(attribute_state$view) || !ncol(attribute_state$view)) {
      stop("This domain has no retained attribute OOS view.", call. = FALSE)
    }
    attribute_query <- .rapid_apply_oos_attribute_view(
      new_attributes, attribute_state$metadata, nrow(feature_query)
    )
    view_state$attribute <- .rapid_oos_weights(
      attribute_state$view, attribute_query, k, zero_tolerance
    )
  }

  active_weights <- weights[names(view_state)]
  active_weights <- active_weights / sum(active_weights)
  embedding <- matrix(
    0, nrow(feature_query), ncol(fit_or_transform$scores[[m]])
  )
  probability <- matrix(
    0, nrow(feature_query),
    ncol(fit_or_transform$prediction_probabilities[[m]])
  )
  for (view in names(view_state)) {
    interpolation <- view_state[[view]]$interpolation
    embedding <- embedding + active_weights[[view]] * as.matrix(
      interpolation %*% fit_or_transform$scores[[m]]
    )
    probability <- probability + active_weights[[view]] * as.matrix(
      interpolation %*% fit_or_transform$prediction_probabilities[[m]]
    )
  }

  # An exact feature match identifies a retained row. Preserve the historical
  # exact-reapplication contract even when auxiliary views contain duplicates.
  exact_feature <- rowSums(
    feature_neighbours$distance <= zero_tolerance
  ) > 0L
  if (weights[["feature"]] > 0 && any(exact_feature)) {
    embedding[exact_feature, ] <- as.matrix(
      feature_neighbours$interpolation[exact_feature, , drop = FALSE] %*%
        fit_or_transform$scores[[m]]
    )
    probability[exact_feature, ] <- as.matrix(
      feature_neighbours$interpolation[exact_feature, , drop = FALSE] %*%
        fit_or_transform$prediction_probabilities[[m]]
    )
  }
  colnames(embedding) <- colnames(fit_or_transform$scores[[m]])
  probability <- probability / pmax(rowSums(probability), 1e-15)
  colnames(probability) <- fit_or_transform$prototypes$class_levels
  choice <- max.col(probability, ties.method = "first")
  predicted <- colnames(probability)[choice]
  confidence <- probability[cbind(seq_len(nrow(probability)), choice)]

  if (type == "embedding") return(embedding)
  if (type == "probabilities") return(probability)
  if (type == "class") return(predicted)
  primary <- if ("feature" %in% names(view_state)) {
    view_state$feature
  } else {
    view_state[[1L]]
  }
  list(
    embedding = embedding,
    probabilities = probability,
    class = predicted,
    confidence = confidence,
    neighbours = primary$index,
    distances = primary$distance,
    weights = primary$weight,
    view_neighbours = lapply(view_state, `[[`, "index"),
    view_distances = lapply(view_state, `[[`, "distance"),
    view_interpolation_weights = lapply(view_state, `[[`, "weight"),
    view_weights = active_weights,
    views_used = names(view_state),
    side = fit_or_transform$domain_names[[m]],
    exact_reapplications = if (weights[["feature"]] > 0) {
      exact_feature
    } else {
      rowSums(primary$distance <= zero_tolerance) > 0L
    }
  )
}

#' @export
oos_predict.rapid_ma_pair_fit <- function(
  fit_or_transform,
  newX,
  side = c("i", "j"),
  ...
) {
  oos_predict(fit_or_transform$core, newX, side = side[[1L]], ...)
}

.rapid_match_anchor_pairs <- function(anchors, n_source, n_target) {
  empty <- matrix(integer(0), 0L, 2L,
                  dimnames = list(NULL, c("source", "target")))
  if (is.null(anchors)) return(empty)
  pairs <- NULL
  if (is.data.frame(anchors) || is.matrix(anchors)) {
    if (ncol(anchors) != 2L) stop("`anchors` must have two columns.", call. = FALSE)
    pairs <- as.matrix(anchors[, 1:2, drop = FALSE])
  } else if (is.list(anchors)) {
    if (!is.null(anchors$source) && !is.null(anchors$target)) {
      pairs <- cbind(anchors$source, anchors$target)
    } else if (!is.null(anchors$i) && !is.null(anchors$j)) {
      pairs <- cbind(anchors$i, anchors$j)
    } else if (!is.null(anchors$vec1) && !is.null(anchors$vec2)) {
      left <- anchors$vec1
      right <- anchors$vec2
      if (length(left) != n_source || length(right) != n_target) {
        stop("Anchor-ID vectors must match the two domain sizes.", call. = FALSE)
      }
      shared <- intersect(unique(left[!is.na(left)]), unique(right[!is.na(right)]))
      if (anyDuplicated(left[!is.na(left)]) || anyDuplicated(right[!is.na(right)])) {
        stop("Non-missing anchor IDs must be unique within each domain.",
             call. = FALSE)
      }
      pairs <- cbind(match(shared, left), match(shared, right))
    } else if (length(anchors) == 2L) {
      pairs <- cbind(anchors[[1L]], anchors[[2L]])
    }
  }
  if (is.null(pairs)) {
    stop("Unsupported `anchors`; use a two-column object, paired indices, or vec1/vec2 IDs.",
         call. = FALSE)
  }
  if (!length(pairs)) return(empty)
  storage.mode(pairs) <- "double"
  if (any(!is.finite(pairs)) || any(abs(pairs - round(pairs)) > 1e-8)) {
    stop("Anchor indices must be finite integers.", call. = FALSE)
  }
  pairs <- matrix(as.integer(pairs), ncol = 2L,
                  dimnames = list(NULL, c("source", "target")))
  if (any(pairs[, 1L] < 1L | pairs[, 1L] > n_source) ||
      any(pairs[, 2L] < 1L | pairs[, 2L] > n_target)) {
    stop("Anchor indices are outside the domain bounds.", call. = FALSE)
  }
  if (anyDuplicated(pairs[, 1L]) || anyDuplicated(pairs[, 2L])) {
    stop("Anchors must be one-to-one.", call. = FALSE)
  }
  pairs
}

.rapid_match_sparse_rows <- function(Q, top = NULL) {
  sm <- Matrix::summary(Q)
  output <- vector("list", nrow(Q))
  if (!nrow(sm)) return(output)
  grouped <- split(seq_len(nrow(sm)), sm$i)
  for (name in names(grouped)) {
    index <- grouped[[name]]
    ordering <- order(-sm$x[index], sm$j[index], method = "radix")
    if (!is.null(top)) ordering <- head(ordering, top)
    index <- index[ordering]
    value <- sm$x[index]
    norm <- sqrt(sum(value * value))
    if (is.finite(norm) && norm > 0) value <- value / norm
    output[[as.integer(name)]] <- list(index = as.integer(sm$j[index]), value = value)
  }
  output
}

.rapid_match_buckets <- function(rows, K, cap) {
  buckets <- vector("list", K)
  row_lengths <- vapply(rows, function(row) {
    if (is.null(row)) 0L else length(row$index)
  }, integer(1))
  if (!sum(row_lengths)) {
    return(lapply(buckets, function(...) {
      matrix(numeric(0), 0L, 2L,
             dimnames = list(NULL, c("node", "weight")))
    }))
  }
  node <- rep.int(seq_along(rows), row_lengths)
  prototype <- unlist(lapply(rows, function(row) {
    if (is.null(row)) integer(0) else row$index
  }), use.names = FALSE)
  weight <- unlist(lapply(rows, function(row) {
    if (is.null(row)) numeric(0) else row$value
  }), use.names = FALSE)
  ordering <- order(prototype, -weight, node, method = "radix")
  node <- node[ordering]
  prototype <- prototype[ordering]
  weight <- weight[ordering]
  grouped <- split(seq_along(prototype), prototype)
  for (name in names(grouped)) {
    index <- head(grouped[[name]], cap)
    buckets[[as.integer(name)]] <- cbind(
      node = node[index], weight = weight[index]
    )
  }
  lapply(buckets, function(bucket) {
    if (is.null(bucket)) {
      matrix(numeric(0), 0L, 2L,
             dimnames = list(NULL, c("node", "weight")))
    } else {
      matrix(bucket, ncol = 2L, dimnames = list(NULL, c("node", "weight")))
    }
  })
}

.rapid_match_similarity <- function(left, right) {
  if (is.null(left) || is.null(right)) return(0)
  shared <- intersect(left$index, right$index)
  if (!length(shared)) return(0)
  sum(left$value[match(shared, left$index)] *
        right$value[match(shared, right$index)])
}

.rapid_match_scale <- function(x) {
  positive <- x[is.finite(x) & x > 0]
  value <- if (length(positive)) stats::median(positive) else 1
  if (!is.finite(value) || value <= 1e-12) 1 else value
}

.rapid_match_view <- function(fit, source, target, view) {
  values <- switch(
    view,
    latent = fit$scores,
    structure = fit$preprocessing$structures,
    attribute = fit$preprocessing$attributes,
    position = lapply(fit$relations, function(relation) {
      state <- relation$views$position
      if (is.null(state)) NULL else state$view
    }),
    NULL
  )
  if (is.null(values) || length(values) < max(source, target) ||
      is.null(values[[source]]) || is.null(values[[target]])) {
    return(NULL)
  }
  left <- as.matrix(values[[source]])
  right <- as.matrix(values[[target]])
  if (!ncol(left) || ncol(left) != ncol(right) ||
      nrow(left) != fit$domain_sizes[[source]] ||
      nrow(right) != fit$domain_sizes[[target]]) {
    return(NULL)
  }
  left[!is.finite(left)] <- 0
  right[!is.finite(right)] <- 0
  list(source = left, target = right)
}

.rapid_match_view_candidates <- function(view, k) {
  k <- min(k, nrow(view$target))
  neighbours <- RANN::nn2(
    data = view$target, query = view$source, k = k,
    treetype = "kd", searchtype = "standard", eps = 0
  )
  distance <- neighbours$nn.dists
  index <- neighbours$nn.idx
  if (is.null(dim(distance))) distance <- matrix(distance, ncol = k)
  if (is.null(dim(index))) index <- matrix(index, ncol = k)
  scale <- .rapid_match_scale(as.numeric(distance))
  strength <- exp(pmax(-distance / scale, -700))
  strength[!is.finite(strength)] <- 0
  list(index = index, distance = distance, strength = strength)
}

#' Bounded multi-view matching for a RAPID-MA fit
#'
#' @param fit A fitted `rapid_ma` object.
#' @param from,to Source and target domains.
#' @param anchors Optional one-to-one anchor pairs or paired anchor-ID vectors.
#' @param prototype_per_node Number of strongest prototype buckets queried.
#' @param candidate_cap Maximum retained target candidates per source node.
#' @param prototype_bucket_cap Maximum target nodes retained per prototype.
#' @param candidate_views Retained comparable views that contribute bounded
#'   nearest-neighbour candidates in addition to prototype buckets.
#' @param view_candidate_k Number of candidates requested from each comparable
#'   view.
#' @param assignment Either one-to-one greedy assignment or independent
#'   retrieval. Pair transforms use the one-to-one default.
#' @param latent_weight,structure_weight,position_weight,attribute_weight,prototype_weight
#'   Nonnegative reranking weights.
#' @return A `rapid_ma_matching` object with a sparse assignment operator,
#'   match table, coverage, confidence, and unmatched-node diagnostics.
#' @export
rapid_ma_match <- function(
  fit,
  from = 1L,
  to = 2L,
  anchors = NULL,
  prototype_per_node = 3L,
  candidate_cap = 32L,
  prototype_bucket_cap = 128L,
  candidate_views = c("position", "attribute", "latent"),
  view_candidate_k = 12L,
  assignment = c("one_to_one", "independent"),
  latent_weight = 1,
  structure_weight = 0.25,
  position_weight = 1,
  attribute_weight = 0.15,
  prototype_weight = 0.25
) {
  if (!inherits(fit, "rapid_ma")) {
    stop("`fit` must be a fitted rapid_ma object.", call. = FALSE)
  }
  source <- .rapid_adapter_side(fit, from)
  target <- .rapid_adapter_side(fit, to)
  if (source == target) stop("`from` and `to` must differ.", call. = FALSE)
  prototype_per_node <- .rapid_int_scalar(
    prototype_per_node, "prototype_per_node"
  )
  candidate_cap <- .rapid_int_scalar(candidate_cap, "candidate_cap")
  prototype_bucket_cap <- .rapid_int_scalar(
    prototype_bucket_cap, "prototype_bucket_cap"
  )
  view_candidate_k <- .rapid_int_scalar(view_candidate_k, "view_candidate_k")
  assignment_mode <- match.arg(assignment)
  allowed_candidate_views <- c("position", "attribute", "structure", "latent")
  if (!is.character(candidate_views) || anyNA(candidate_views) ||
      any(!nzchar(candidate_views)) || anyDuplicated(candidate_views) ||
      any(!candidate_views %in% allowed_candidate_views)) {
    stop(
      "`candidate_views` must be a unique subset of: ",
      paste(allowed_candidate_views, collapse = ", "), ".",
      call. = FALSE
    )
  }
  view_weights <- c(
    latent = latent_weight,
    structure = structure_weight,
    position = position_weight,
    attribute = attribute_weight,
    prototype = prototype_weight
  )
  for (name in names(view_weights)) {
    view_weights[[name]] <- .rapid_number_scalar(
      view_weights[[name]], paste0(name, "_weight"), lower = 0
    )
  }
  if (sum(view_weights) <= 0) {
    stop("At least one matching weight must be positive.", call. = FALSE)
  }

  Q_source <- fit$couplings[[source]]
  Q_target <- fit$couplings[[target]]
  source_rows <- .rapid_match_sparse_rows(Q_source, top = prototype_per_node)
  target_rows <- .rapid_match_sparse_rows(Q_target)
  buckets <- .rapid_match_buckets(
    target_rows, ncol(Q_target), prototype_bucket_cap
  )
  candidate_state <- list()
  for (view in candidate_views) {
    matrices <- .rapid_match_view(fit, source, target, view)
    if (!is.null(matrices)) {
      candidate_state[[view]] <- .rapid_match_view_candidates(
        matrices, view_candidate_k
      )
    }
  }
  n_source <- nrow(Q_source)
  n_target <- nrow(Q_target)
  anchor_pairs <- .rapid_match_anchor_pairs(anchors, n_source, n_target)
  anchor_by_source <- stats::setNames(anchor_pairs[, 2L], anchor_pairs[, 1L])
  edge_target_rows <- vector("list", n_source)
  edge_prototype_rows <- vector("list", n_source)

  for (i in seq_len(n_source)) {
    row <- source_rows[[i]]
    candidates <- strengths <- numeric(0)
    if (!is.null(row)) {
      for (j in seq_along(row$index)) {
        bucket <- buckets[[row$index[[j]]]]
        if (!nrow(bucket)) next
        candidates <- c(candidates, as.integer(bucket[, "node"]))
        strengths <- c(strengths, row$value[[j]] * bucket[, "weight"])
      }
    }
    if (length(candidate_state)) {
      for (view in names(candidate_state)) {
        state <- candidate_state[[view]]
        candidates <- c(candidates, as.integer(state$index[i, ]))
        strengths <- c(strengths, as.numeric(state$strength[i, ]))
      }
    }
    if (as.character(i) %in% names(anchor_by_source)) {
      candidates <- c(candidates, anchor_by_source[[as.character(i)]])
      strengths <- c(strengths, Inf)
    }
    if (!length(candidates)) next
    aggregate_strength <- tapply(strengths, candidates, max)
    candidate <- as.integer(names(aggregate_strength))
    ordering <- order(-aggregate_strength, candidate, method = "radix")
    candidate <- head(candidate[ordering], candidate_cap)
    edge_target_rows[[i]] <- candidate
    edge_prototype_rows[[i]] <- vapply(candidate, function(j) {
        .rapid_match_similarity(source_rows[[i]], target_rows[[j]])
      }, numeric(1))
  }

  edge_lengths <- lengths(edge_target_rows)
  edge_source <- rep.int(seq_len(n_source), edge_lengths)
  edge_target <- as.integer(unlist(edge_target_rows, use.names = FALSE))
  edge_prototype <- as.numeric(unlist(
    edge_prototype_rows, use.names = FALSE
  ))

  represented <- unique(edge_source)
  missing_source <- setdiff(seq_len(n_source), represented)
  if (length(missing_source)) {
    k_fallback <- min(candidate_cap, n_target)
    fallback <- RANN::nn2(
      data = fit$scores[[target]],
      query = fit$scores[[source]][missing_source, , drop = FALSE],
      k = k_fallback
    )
    edge_source <- c(
      edge_source,
      rep(missing_source, each = k_fallback)
    )
    edge_target <- c(edge_target, as.integer(c(t(fallback$nn.idx))))
    edge_prototype <- c(edge_prototype, rep(0, length(missing_source) * k_fallback))
  }
  if (!length(edge_source)) {
    stop("No fine-matching candidates were available.", call. = FALSE)
  }

  squared_distance <- function(left, right) {
    rowSums((left[edge_source, , drop = FALSE] -
               right[edge_target, , drop = FALSE])^2)
  }
  latent_distance <- squared_distance(fit$scores[[source]], fit$scores[[target]])
  structure_distance <- if (!is.null(fit$preprocessing$structures)) {
    squared_distance(
      fit$preprocessing$structures[[source]],
      fit$preprocessing$structures[[target]]
    )
  } else rep(0, length(edge_source))
  position_view <- .rapid_match_view(fit, source, target, "position")
  position_distance <- if (!is.null(position_view)) {
    squared_distance(position_view$source, position_view$target)
  } else rep(0, length(edge_source))
  attribute_distance <- if (!is.null(fit$preprocessing$attributes)) {
    squared_distance(
      fit$preprocessing$attributes[[source]],
      fit$preprocessing$attributes[[target]]
    )
  } else rep(0, length(edge_source))
  cost <- view_weights[["latent"]] *
    latent_distance / .rapid_match_scale(latent_distance) +
    view_weights[["structure"]] *
    structure_distance / .rapid_match_scale(structure_distance) +
    view_weights[["position"]] *
    position_distance / .rapid_match_scale(position_distance) +
    view_weights[["attribute"]] *
    attribute_distance / .rapid_match_scale(attribute_distance) +
    view_weights[["prototype"]] * (1 - edge_prototype)
  anchor_key <- paste(anchor_pairs[, 1L], anchor_pairs[, 2L], sep = ":")
  edge_key <- paste(edge_source, edge_target, sep = ":")
  cost[edge_key %in% anchor_key] <- -Inf

  assigned_source <- rep(FALSE, n_source)
  assigned_target <- rep(FALSE, n_target)
  match_target <- rep(NA_integer_, n_source)
  match_cost <- rep(NA_real_, n_source)
  anchored <- rep(FALSE, n_source)
  if (nrow(anchor_pairs)) {
    match_target[anchor_pairs[, 1L]] <- anchor_pairs[, 2L]
    match_cost[anchor_pairs[, 1L]] <- 0
    anchored[anchor_pairs[, 1L]] <- TRUE
    assigned_source[anchor_pairs[, 1L]] <- TRUE
    assigned_target[anchor_pairs[, 2L]] <- TRUE
  }
  ordering <- order(cost, edge_source, edge_target, method = "radix")
  for (e in ordering) {
    i <- edge_source[[e]]
    j <- edge_target[[e]]
    if (assigned_source[[i]]) next
    if (identical(assignment_mode, "one_to_one") && assigned_target[[j]]) next
    if (identical(assignment_mode, "independent") && assigned_target[[j]] &&
        any(anchor_pairs[, 2L] == j)) next
    assigned_source[[i]] <- TRUE
    if (identical(assignment_mode, "one_to_one")) assigned_target[[j]] <- TRUE
    match_target[[i]] <- j
    match_cost[[i]] <- cost[[e]]
  }
  matched <- which(!is.na(match_target))
  confidence_scale <- .rapid_match_scale(match_cost[matched[!anchored[matched]]])
  confidence <- rep(0, n_source)
  confidence[matched] <- exp(-pmax(match_cost[matched], 0) / confidence_scale)
  confidence[anchored] <- 1
  assignment <- Matrix::sparseMatrix(
    i = match_target[matched],
    j = matched,
    x = rep(1, length(matched)),
    dims = c(n_target, n_source)
  )
  matches <- data.frame(
    source = seq_len(n_source),
    target = match_target,
    cost = match_cost,
    confidence = confidence,
    anchored = anchored,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      assignment = assignment,
      matches = matches,
      coverage = length(matched) / n_source,
      mean_confidence = if (length(matched)) mean(confidence[matched]) else 0,
      unmatched_source = which(is.na(match_target)),
      unmatched_target = setdiff(seq_len(n_target), match_target[matched]),
      candidate_edges = as.integer(length(edge_source)),
      candidate_cap = candidate_cap,
      candidate_views = names(candidate_state),
      view_candidate_k = view_candidate_k,
      assignment_mode = assignment_mode,
      anchors = anchor_pairs,
      from = fit$domain_names[[source]],
      to = fit$domain_names[[target]],
      weights = view_weights,
      dense_pairwise_allocated = FALSE
    ),
    class = "rapid_ma_matching"
  )
}

#' @rdname relative_transform
#' @export
relative_transform.rapid_ma_pair_fit <- function(
  fit,
  from = c("i", "j"),
  to = c("j", "i"),
  ...
) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) stop("`from` and `to` must differ.", call. = FALSE)
  if (is.null(fit$matching)) {
    stop("This pair fit was created with `fine_match = FALSE`.", call. = FALSE)
  }
  op <- fit$matching$assignment
  if (from == "j" && to == "i") op <- .normalize_rows_sparse(Matrix::t(op))
  new_align_transform(
    "perm", op, from = from, to = to,
    dim = c(fit$n1, fit$n2)
  )
}

#' @rdname pair_loss
#' @export
pair_loss.rapid_ma_pair_fit <- function(fit, X_i = NULL, X_j = NULL, ...) {
  fit$objective
}

#' @rdname latent_dim
#' @export
latent_dim.rapid_ma_pair_fit <- function(fit, ...) fit$k
