# RAPID-MA relation-aware sparse diffusion encoder --------------------------

.rapid_diffusion_control <- function(
  steps = c(0L, 1L, 2L, 4L),
  sketch_dim = 64L,
  output_dim = 32L,
  gate = c("energy", "uniform", "fixed"),
  relation_weights = NULL,
  gate_temperature = 2,
  label_gate_weight = 1,
  include_position = TRUE,
  include_attribute = TRUE,
  include_signatures = TRUE,
  direction_mode = c("bidirectional", "forward"),
  reverse_degree_cap = NULL,
  block_size = 32L,
  ridge = 1e-5,
  store_propagations = FALSE,
  seed = 1L
) {
  gate <- match.arg(gate)
  direction_mode <- match.arg(direction_mode)
  if (!is.numeric(steps) || !length(steps) || any(!is.finite(steps)) ||
      any(steps < 0) || any(abs(steps - round(steps)) > 1e-8)) {
    stop("`steps` must contain nonnegative integers.", call. = FALSE)
  }
  steps <- sort(unique(c(0L, as.integer(steps))))
  if (max(steps) > 64L) {
    stop("`steps` cannot exceed 64 in the default diffusion encoder.",
         call. = FALSE)
  }
  sketch_dim <- .rapid_int_scalar(sketch_dim, "sketch_dim")
  output_dim <- .rapid_int_scalar(output_dim, "output_dim")
  block_size <- .rapid_int_scalar(block_size, "block_size")
  seed <- .rapid_int_scalar(seed, "seed", min_value = 0L)
  gate_temperature <- .rapid_number_scalar(
    gate_temperature, "gate_temperature", lower = 0
  )
  label_gate_weight <- .rapid_number_scalar(
    label_gate_weight, "label_gate_weight", lower = 0
  )
  ridge <- .rapid_number_scalar(ridge, "ridge", lower = 0, strict = TRUE)
  if (!is.null(reverse_degree_cap)) {
    reverse_degree_cap <- .rapid_int_scalar(
      reverse_degree_cap, "reverse_degree_cap"
    )
  }
  flags <- list(
    include_position = include_position,
    include_attribute = include_attribute,
    include_signatures = include_signatures,
    store_propagations = store_propagations
  )
  for (name in names(flags)) {
    value <- flags[[name]]
    if (length(value) != 1L || !is.logical(value) || is.na(value)) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!is.null(relation_weights)) {
    if (!is.numeric(relation_weights) || any(!is.finite(relation_weights)) ||
        any(relation_weights < 0)) {
      stop("`relation_weights` must be finite and nonnegative.", call. = FALSE)
    }
    if (is.null(names(relation_weights)) || any(!nzchar(names(relation_weights))) ||
        anyDuplicated(names(relation_weights))) {
      stop("`relation_weights` must have unique non-empty names.", call. = FALSE)
    }
  }
  if (identical(gate, "fixed") && is.null(relation_weights)) {
    stop("`gate = \"fixed\"` requires `relation_weights`.", call. = FALSE)
  }

  structure(
    list(
      steps = steps,
      sketch_dim = sketch_dim,
      output_dim = output_dim,
      gate = gate,
      relation_weights = relation_weights,
      gate_temperature = gate_temperature,
      label_gate_weight = label_gate_weight,
      include_position = include_position,
      include_attribute = include_attribute,
      include_signatures = include_signatures,
      direction_mode = direction_mode,
      reverse_degree_cap = reverse_degree_cap,
      block_size = block_size,
      ridge = ridge,
      store_propagations = store_propagations,
      seed = seed
    ),
    class = "rapid_ma_diffusion_control"
  )
}

.rapid_resolve_diffusion_control <- function(control = NULL) {
  defaults <- unclass(.rapid_diffusion_control())
  if (is.null(control)) return(do.call(.rapid_diffusion_control, defaults))
  if (inherits(control, "rapid_ma_diffusion_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`control` must be NULL, a named list, or rapid_ma_diffusion_control.",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown diffusion control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  do.call(.rapid_diffusion_control, utils::modifyList(defaults, control))
}

.rapid_validate_embedding <- function(Z, n, name) {
  .rapid_validate_numeric_matrix(Z, name, n_rows = n)
  Z <- as.matrix(Z)
  if (any(!is.finite(Z))) {
    stop("`", name, "` must contain only finite values.", call. = FALSE)
  }
  Z
}

.rapid_named_view <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.list(x) && !is.null(x$view)) return(x$view)
  NULL
}

.rapid_project_view <- function(Z, width, seed, prefix) {
  Z <- as.matrix(Z)
  if (!ncol(Z)) {
    return(list(
      view = matrix(numeric(0), nrow(Z), 0L),
      projection = NULL,
      projected = FALSE,
      input_dim = 0L,
      output_dim = 0L,
      center = numeric(0),
      scale = numeric(0)
    ))
  }
  width <- min(.rapid_int_scalar(width, "width"), ncol(Z))
  projected <- ncol(Z) > width
  projection <- NULL
  if (projected) {
    keys <- paste0(prefix, ":", seq_len(ncol(Z)))
    buckets <- vapply(
      keys, .rapid_hash_bucket, integer(1), width = width, seed = seed
    )
    signs <- vapply(keys, .rapid_hash_sign, numeric(1), seed = seed)
    projection <- Matrix::sparseMatrix(
      i = seq_len(ncol(Z)), j = buckets, x = signs,
      dims = c(ncol(Z), width)
    )
    projected_view <- as.matrix(Z %*% projection)
  } else {
    projected_view <- Z
  }
  scaled <- .rapid_scale_dense(projected_view, drop_constant = TRUE)
  if (!ncol(scaled$view)) {
    scaled$view <- matrix(0, nrow(projected_view), 1L)
    colnames(scaled$view) <- paste0(prefix, "__constant")
  }
  list(
    view = scaled$view,
    projection = projection,
    projected = projected,
    input_dim = as.integer(ncol(Z)),
    output_dim = as.integer(ncol(scaled$view)),
    center = scaled$center,
    scale = scaled$scale
  )
}

.rapid_diffusion_base <- function(X,
                                  relation_state,
                                  control,
                                  precomputed = NULL,
                                  encoder = NULL) {
  n <- nrow(X)
  if (!is.null(precomputed) && !is.null(encoder)) {
    stop("Supply only one of `precomputed` and `encoder`.", call. = FALSE)
  }
  if (!is.null(encoder)) {
    if (!is.function(encoder)) stop("`encoder` must be a function.", call. = FALSE)
    encoded <- encoder(X, relation_state)
    base <- .rapid_validate_embedding(encoded, n, "encoder output")
    source <- "callback"
  } else if (!is.null(precomputed)) {
    base <- .rapid_validate_embedding(precomputed, n, "precomputed")
    source <- "precomputed"
  } else {
    feature <- .rapid_feature_view(
      X, max_dim = control$sketch_dim, seed = control$seed
    )$view
    blocks <- list(feature = feature)
    if (isTRUE(control$include_position)) {
      position <- .rapid_named_view(relation_state$views$position)
      if (!is.null(position) && ncol(position)) blocks$position <- position
    }
    if (isTRUE(control$include_attribute)) {
      attribute <- .rapid_named_view(relation_state$views$attribute)
      if (!is.null(attribute) && ncol(attribute)) blocks$attribute <- attribute
    }
    base <- do.call(cbind, blocks)
    source <- paste(names(blocks), collapse = "+")
  }

  projected <- .rapid_project_view(
    base,
    width = control$sketch_dim,
    seed = control$seed + 101L,
    prefix = "diffusion_base"
  )
  if (!ncol(projected$view)) {
    stop("The diffusion base has no informative columns.", call. = FALSE)
  }
  projected$source <- source
  projected
}

.rapid_sparse_multiply_blocked <- function(A, X, block_size = 32L) {
  X <- as.matrix(X)
  if (!ncol(X)) return(X)
  block_size <- .rapid_int_scalar(block_size, "block_size")
  out <- matrix(0, nrow(A), ncol(X))
  starts <- seq.int(1L, ncol(X), by = block_size)
  for (start in starts) {
    end <- min(ncol(X), start + block_size - 1L)
    out[, start:end] <- as.matrix(A %*% X[, start:end, drop = FALSE])
  }
  out
}

.rapid_relation_label_agreement <- function(A, labels) {
  if (is.null(labels)) return(list(agreement = NA_real_, edge_mass = 0))
  if (length(labels) != nrow(A)) {
    stop("`labels` must have one value per node.", call. = FALSE)
  }
  observed <- !is.na(labels)
  sm <- Matrix::summary(A)
  if (!nrow(sm)) return(list(agreement = NA_real_, edge_mass = 0))
  keep <- observed[sm$i] & observed[sm$j]
  if (!any(keep)) return(list(agreement = NA_real_, edge_mass = 0))
  mass <- sum(sm$x[keep])
  if (!is.finite(mass) || mass <= 0) {
    return(list(agreement = NA_real_, edge_mass = 0))
  }
  agreement <- sum(sm$x[keep] * (labels[sm$i[keep]] == labels[sm$j[keep]])) / mass
  list(agreement = as.numeric(agreement), edge_mass = as.numeric(mass))
}

.rapid_relation_quality <- function(A, base, labels = NULL) {
  row_mass <- as.numeric(Matrix::rowSums(A))
  active <- which(row_mass > 0)
  if (!length(active) || Matrix::nnzero(A) == 0L) {
    return(list(
      energy = Inf,
      connectivity = 0,
      label_agreement = NA_real_,
      labeled_edge_mass = 0,
      nnz = 0L
    ))
  }
  P <- Matrix::Diagonal(x = ifelse(row_mass > 0, 1 / row_mass, 0)) %*% A
  neighbor <- .rapid_sparse_multiply_blocked(P, base, block_size = ncol(base))
  residual <- base[active, , drop = FALSE] - neighbor[active, , drop = FALSE]
  numerator <- mean(rowSums(residual * residual))
  denominator <- mean(rowSums(base[active, , drop = FALSE]^2))
  energy <- numerator / max(denominator, 1e-12)
  label <- .rapid_relation_label_agreement(A, labels)
  list(
    energy = as.numeric(energy),
    connectivity = as.numeric(length(active) / nrow(A)),
    label_agreement = label$agreement,
    labeled_edge_mass = label$edge_mass,
    nnz = as.integer(Matrix::nnzero(A))
  )
}

.rapid_relation_weights <- function(relations, base, control, labels = NULL) {
  relation_names <- names(relations)
  if (!length(relation_names)) {
    return(list(weights = numeric(0), quality = list(), raw_scores = numeric(0)))
  }
  quality <- lapply(relations, .rapid_relation_quality, base = base, labels = labels)
  names(quality) <- relation_names
  available <- vapply(quality, function(x) x$nnz > 0L, logical(1))
  fixed <- control$relation_weights
  if (!is.null(fixed)) {
    unknown <- setdiff(names(fixed), relation_names)
    if (length(unknown)) {
      stop("Unknown relation weight name(s): ", paste(unknown, collapse = ", "),
           ".", call. = FALSE)
    }
  }

  if (identical(control$gate, "fixed")) {
    raw <- stats::setNames(rep(0, length(relations)), relation_names)
    raw[names(fixed)] <- fixed
  } else if (identical(control$gate, "uniform")) {
    raw <- stats::setNames(as.numeric(available), relation_names)
    if (!is.null(fixed)) raw[names(fixed)] <- raw[names(fixed)] * fixed
  } else {
    raw <- vapply(quality, function(x) {
      if (!is.finite(x$energy) || x$nnz == 0L) return(0)
      score <- exp(pmax(-control$gate_temperature * x$energy, -700)) *
        sqrt(max(x$connectivity, 0))
      if (is.finite(x$label_agreement)) {
        score <- score * exp(
          control$label_gate_weight * (2 * x$label_agreement - 1)
        )
      }
      score
    }, numeric(1))
    if (!is.null(fixed)) raw[names(fixed)] <- raw[names(fixed)] * fixed
  }
  raw[!available] <- 0
  total <- sum(raw)
  weights <- if (is.finite(total) && total > 0) raw / total else {
    fallback <- as.numeric(available)
    if (sum(fallback) > 0) fallback / sum(fallback) else fallback
  }
  names(weights) <- relation_names
  list(weights = weights, quality = quality, raw_scores = raw)
}

.rapid_relation_channels <- function(relations, metadata, weights, control) {
  channels <- list()
  channel_weights <- numeric(0)
  channel_relation <- character(0)
  for (name in names(relations)) {
    A <- .rapid_force_dgC(relations[[name]])
    if (Matrix::nnzero(A) == 0L || weights[[name]] <= 0) next
    is_symmetric <- isTRUE(isSymmetric(A, tol = 1e-12))
    add_reverse <- identical(control$direction_mode, "bidirectional") && !is_symmetric
    if (add_reverse) {
      cap <- control$reverse_degree_cap
      if (is.null(cap)) {
        cap <- .rapid_metadata_value(metadata[[name]], "degree_cap", 64L)
      }
      reverse <- .rapid_cap_out_degree(.rapid_force_dgC(Matrix::t(A)), cap)
      reverse <- .rapid_normalize_relation(reverse, "random_walk")$A
      channels[[paste0(name, "__forward")]] <- A
      channels[[paste0(name, "__reverse")]] <- reverse
      channel_weights <- c(channel_weights, weights[[name]] / 2, weights[[name]] / 2)
      channel_relation <- c(channel_relation, name, name)
    } else {
      channels[[name]] <- A
      channel_weights <- c(channel_weights, weights[[name]])
      channel_relation <- c(channel_relation, name)
    }
  }
  names(channel_weights) <- names(channels)
  names(channel_relation) <- names(channels)
  list(
    channels = channels,
    weights = channel_weights,
    relation = channel_relation
  )
}

.rapid_structural_signatures <- function(relations, n = NULL) {
  if (!length(relations)) {
    if (is.null(n)) n <- 0L
    return(matrix(numeric(0), nrow = n, ncol = 0L))
  }
  n <- nrow(relations[[1L]])
  blocks <- lapply(names(relations), function(name) {
    A <- relations[[name]]
    counts <- .rapid_degree_counts(A)
    cbind(
      log_out_degree = log1p(counts$out_degree),
      log_in_degree = log1p(counts$in_degree),
      isolate = as.numeric((counts$out_degree + counts$in_degree) == 0L)
    )
  })
  names(blocks) <- names(relations)
  signature <- do.call(cbind, blocks)
  colnames(signature) <- unlist(lapply(names(blocks), function(name) {
    paste0(name, c("__log_out_degree", "__log_in_degree", "__isolate"))
  }), use.names = FALSE)
  if (!nrow(signature) && n > 0L) signature <- matrix(0, n, 0L)
  signature
}

.rapid_propagate_channels <- function(channels, base, steps, block_size,
                                      store_propagations = FALSE) {
  requested <- steps[steps > 0L]
  blocks <- list()
  stored <- list()
  if (!length(requested) || !length(channels)) {
    return(list(blocks = blocks, propagations = stored, multiplies = 0L))
  }
  max_step <- max(requested)
  multiplies <- 0L
  for (name in names(channels)) {
    current <- base
    channel_store <- list()
    for (step in seq_len(max_step)) {
      current <- .rapid_sparse_multiply_blocked(
        channels[[name]], current, block_size = block_size
      )
      multiplies <- multiplies + 1L
      if (step %in% requested) {
        key <- paste0(name, "__step", step)
        blocks[[key]] <- current
        if (isTRUE(store_propagations)) {
          channel_store[[as.character(step)]] <- current
        }
      }
    }
    if (isTRUE(store_propagations)) stored[[name]] <- channel_store
  }
  list(blocks = blocks, propagations = stored, multiplies = as.integer(multiplies))
}

.rapid_compress_diffusion <- function(fused, output_dim, ridge, seed) {
  projected <- .rapid_project_view(
    fused,
    width = min(output_dim, ncol(fused)),
    seed = seed + 5003L,
    prefix = "diffusion_output"
  )
  Z <- projected$view
  denom <- max(nrow(Z) - 1L, 1L)
  covariance <- crossprod(Z) / denom
  ridge_used <- ridge
  R <- tryCatch(
    chol(covariance + diag(ridge_used, ncol(Z))),
    error = function(...) NULL
  )
  if (is.null(R)) {
    ridge_used <- max(ridge, 1e-3)
    R <- chol(covariance + diag(ridge_used, ncol(Z)))
  }
  whitening <- backsolve(R, diag(ncol(Z)))
  embedding <- Z %*% whitening
  embedding[!is.finite(embedding)] <- 0
  list(
    embedding = embedding,
    projection = projected,
    whitening = whitening,
    ridge = ridge_used
  )
}

# Encode one domain with bounded sparse relation diffusion.
.rapid_diffusion_encode <- function(
  X,
  relation_state,
  control = NULL,
  labels = NULL,
  precomputed = NULL,
  encoder = NULL
) {
  .rapid_validate_numeric_matrix(X, "X")
  n <- nrow(X)
  if (!inherits(relation_state, "rapid_ma_relations")) {
    stop("`relation_state` must come from `.rapid_prepare_relations()`.",
         call. = FALSE)
  }
  if (relation_state$n_nodes != n) {
    stop("`relation_state` and `X` have different node counts.", call. = FALSE)
  }
  if (!is.null(labels) && length(labels) != n) {
    stop("`labels` must have one value per node.", call. = FALSE)
  }
  control <- .rapid_resolve_diffusion_control(control)
  base_fit <- .rapid_diffusion_base(
    X, relation_state, control, precomputed = precomputed, encoder = encoder
  )
  base <- base_fit$view

  gated <- .rapid_relation_weights(
    relation_state$relations, base, control, labels = labels
  )
  channel_state <- .rapid_relation_channels(
    relation_state$relations,
    relation_state$metadata,
    gated$weights,
    control
  )
  propagated <- .rapid_propagate_channels(
    channel_state$channels,
    base,
    steps = control$steps,
    block_size = control$block_size,
    store_propagations = control$store_propagations
  )

  blocks <- list(base = base)
  positive_steps <- control$steps[control$steps > 0L]
  if (length(propagated$blocks)) {
    for (key in names(propagated$blocks)) {
      channel <- sub("__step[0-9]+$", "", key)
      weight <- channel_state$weights[[channel]] / max(length(positive_steps), 1L)
      scaled <- .rapid_scale_dense(
        propagated$blocks[[key]], drop_constant = TRUE
      )$view
      blocks[[key]] <- sqrt(max(weight, 0)) * scaled
    }
  }

  signatures <- .rapid_structural_signatures(relation_state$relations, n = n)
  if (isTRUE(control$include_signatures) && ncol(signatures)) {
    signature_view <- .rapid_scale_dense(signatures, drop_constant = TRUE)$view
    if (ncol(signature_view)) blocks$structural_signatures <- signature_view
  }
  fused <- do.call(cbind, blocks)
  compressed <- .rapid_compress_diffusion(
    fused,
    output_dim = control$output_dim,
    ridge = control$ridge,
    seed = control$seed
  )

  diagnostics <- list(
    n_nodes = as.integer(n),
    base_dim = as.integer(ncol(base)),
    fused_dim = as.integer(ncol(fused)),
    output_dim = as.integer(ncol(compressed$embedding)),
    relation_weights = gated$weights,
    relation_quality = gated$quality,
    channel_weights = channel_state$weights,
    channels = names(channel_state$channels),
    steps = control$steps,
    sparse_multiplies = propagated$multiplies,
    total_relation_nnz = as.integer(sum(vapply(
      relation_state$relations, Matrix::nnzero, numeric(1)
    ))),
    channel_nnz = as.integer(sum(vapply(
      channel_state$channels, Matrix::nnzero, numeric(1)
    ))),
    channel_nnz_by_name = vapply(
      channel_state$channels, Matrix::nnzero, numeric(1)
    ),
    channel_max_out_degree = vapply(
      channel_state$channels,
      function(A) {
        degree <- .rapid_degree_counts(A)$out_degree
        if (length(degree)) as.integer(max(degree)) else 0L
      },
      integer(1)
    ),
    dense_node_square_allocated = FALSE,
    seed = control$seed
  )

  structure(
    list(
      embedding = compressed$embedding,
      base = base,
      signatures = signatures,
      relation_weights = gated$weights,
      relation_quality = gated$quality,
      channel_weights = channel_state$weights,
      propagations = if (isTRUE(control$store_propagations)) {
        propagated$propagations
      } else {
        NULL
      },
      encoder_fit = list(
        base = base_fit,
        output_projection = compressed$projection,
        whitening = compressed$whitening,
        ridge = compressed$ridge,
        control = control
      ),
      diagnostics = diagnostics
    ),
    class = "rapid_ma_diffusion"
  )
}
