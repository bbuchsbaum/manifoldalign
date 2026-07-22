# RAPID-MA shared structural prototype bank -------------------------------

.rapid_prototype_control <- function(
  prototypes_per_class = 2L,
  unknown_prototypes = 2L,
  unknown_level = ".unknown",
  embedding_weight = 1,
  structure_weight = 1,
  attribute_weight = 1,
  min_mass = 1e-8,
  class_smoothing = 1e-3,
  max_recoveries = 64L,
  seed = 1L
) {
  prototypes_per_class <- .rapid_int_scalar(
    prototypes_per_class, "prototypes_per_class"
  )
  unknown_prototypes <- .rapid_int_scalar(
    unknown_prototypes, "unknown_prototypes", min_value = 0L
  )
  max_recoveries <- .rapid_int_scalar(
    max_recoveries, "max_recoveries", min_value = 0L
  )
  seed <- .rapid_int_scalar(seed, "seed", min_value = 0L)
  embedding_weight <- .rapid_number_scalar(
    embedding_weight, "embedding_weight", lower = 0
  )
  structure_weight <- .rapid_number_scalar(
    structure_weight, "structure_weight", lower = 0
  )
  attribute_weight <- .rapid_number_scalar(
    attribute_weight, "attribute_weight", lower = 0
  )
  min_mass <- .rapid_number_scalar(
    min_mass, "min_mass", lower = 0, strict = TRUE
  )
  class_smoothing <- .rapid_number_scalar(
    class_smoothing, "class_smoothing", lower = 0
  )
  if (length(unknown_level) != 1L || !is.character(unknown_level) ||
      is.na(unknown_level) || !nzchar(unknown_level)) {
    stop("`unknown_level` must be one non-empty string.", call. = FALSE)
  }
  if ((embedding_weight + structure_weight + attribute_weight) <= 0) {
    stop("At least one prototype view weight must be positive.", call. = FALSE)
  }

  structure(
    list(
      prototypes_per_class = prototypes_per_class,
      unknown_prototypes = unknown_prototypes,
      unknown_level = unknown_level,
      embedding_weight = embedding_weight,
      structure_weight = structure_weight,
      attribute_weight = attribute_weight,
      min_mass = min_mass,
      class_smoothing = class_smoothing,
      max_recoveries = max_recoveries,
      seed = seed
    ),
    class = "rapid_ma_prototype_control"
  )
}

.rapid_resolve_prototype_control <- function(control = NULL) {
  defaults <- unclass(.rapid_prototype_control())
  if (is.null(control)) return(do.call(.rapid_prototype_control, defaults))
  if (inherits(control, "rapid_ma_prototype_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`control` must be NULL, a named list, or rapid_ma_prototype_control.",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown prototype control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  do.call(.rapid_prototype_control, utils::modifyList(defaults, control))
}

.rapid_as_domain_matrices <- function(x, name, allow_null = FALSE,
                                      n_by_domain = NULL) {
  if (is.null(x)) {
    if (!allow_null) stop("`", name, "` cannot be NULL.", call. = FALSE)
    if (is.null(n_by_domain)) return(NULL)
    return(lapply(n_by_domain, function(n) matrix(numeric(0), n, 0L)))
  }
  if (is.matrix(x) || inherits(x, "Matrix")) x <- list(x)
  if (!is.list(x) || !length(x)) {
    stop("`", name, "` must be a matrix or a non-empty list of matrices.",
         call. = FALSE)
  }
  out <- lapply(seq_along(x), function(i) {
    value <- x[[i]]
    .rapid_validate_numeric_matrix(value, paste0(name, "[[", i, "]]"))
    value <- as.matrix(value)
    if (any(!is.finite(value))) {
      stop("`", name, "[[", i, "]]` must contain only finite values.",
           call. = FALSE)
    }
    value
  })
  dims <- vapply(out, ncol, integer(1))
  if (length(unique(dims)) != 1L) {
    stop("All `", name, "` matrices must have the same number of columns.",
         call. = FALSE)
  }
  if (!is.null(n_by_domain) &&
      (!identical(length(out), length(n_by_domain)) ||
       any(vapply(out, nrow, integer(1)) != n_by_domain))) {
    stop("`", name, "` must match the domain row counts.", call. = FALSE)
  }
  out
}

.rapid_as_domain_labels <- function(labels, n_by_domain) {
  if (is.null(labels)) {
    return(lapply(n_by_domain, function(n) rep(NA_character_, n)))
  }
  if (!is.list(labels)) labels <- list(labels)
  if (length(labels) != length(n_by_domain)) {
    stop("`labels` must match the number of domains.", call. = FALSE)
  }
  lapply(seq_along(labels), function(i) {
    value <- labels[[i]]
    if (length(value) != n_by_domain[[i]]) {
      stop("`labels[[", i, "]]` must have one value per node.", call. = FALSE)
    }
    as.character(value)
  })
}

.rapid_canonical_order <- function(view, labels) {
  if (!nrow(view)) return(integer(0))
  columns <- if (ncol(view)) {
    lapply(seq_len(ncol(view)), function(j) view[, j])
  } else {
    list(rep(0, nrow(view)))
  }
  label_key <- ifelse(is.na(labels), "", labels)
  do.call(order, c(list(label_key), columns, list(method = "radix", na.last = TRUE)))
}

.rapid_weighted_prototype_view <- function(embedding, structure, attribute,
                                           control) {
  blocks <- list()
  add_block <- function(x, weight) {
    if (!ncol(x) || weight <= 0) return(NULL)
    sqrt(weight) * .rapid_scale_dense(x, drop_constant = FALSE)$view
  }
  blocks$embedding <- add_block(embedding, control$embedding_weight)
  blocks$structure <- add_block(structure, control$structure_weight)
  blocks$attribute <- add_block(attribute, control$attribute_weight)
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  if (!length(blocks)) {
    stop("Prototype initialization has no non-empty weighted view.", call. = FALSE)
  }
  do.call(cbind, blocks)
}

.rapid_prototype_observations <- function(embeddings, labels = NULL,
                                          structures = NULL, attributes = NULL,
                                          control = NULL) {
  control <- .rapid_resolve_prototype_control(control)
  embedding_domains <- .rapid_as_domain_matrices(embeddings, "embeddings")
  n_by_domain <- vapply(embedding_domains, nrow, integer(1))
  if (!sum(n_by_domain)) stop("At least one observation is required.", call. = FALSE)
  structure_domains <- .rapid_as_domain_matrices(
    structures, "structures", allow_null = TRUE, n_by_domain = n_by_domain
  )
  attribute_domains <- .rapid_as_domain_matrices(
    attributes, "attributes", allow_null = TRUE, n_by_domain = n_by_domain
  )
  label_domains <- .rapid_as_domain_labels(labels, n_by_domain)

  embedding <- do.call(rbind, embedding_domains)
  structure <- do.call(rbind, structure_domains)
  attribute <- do.call(rbind, attribute_domains)
  label <- unlist(label_domains, use.names = FALSE)
  if (control$unknown_level %in% label[!is.na(label)]) {
    stop("`unknown_level` collides with an observed class label.", call. = FALSE)
  }
  selection_view <- .rapid_weighted_prototype_view(
    embedding, structure, attribute, control
  )
  canonical_order <- .rapid_canonical_order(selection_view, label)
  domain_index <- rep(seq_along(n_by_domain), n_by_domain)
  row_index <- sequence(n_by_domain)

  list(
    embedding = embedding,
    structure = structure,
    attribute = attribute,
    label = label,
    selection_view = selection_view,
    canonical_order = canonical_order,
    domain_index = domain_index,
    row_index = row_index,
    n_by_domain = n_by_domain,
    control = control
  )
}

.rapid_sqdist_to_centers <- function(X, centers) {
  if (!nrow(centers)) return(matrix(Inf, nrow(X), 0L))
  distance <- outer(rowSums(X * X), rowSums(centers * centers), "+") -
    2 * tcrossprod(X, centers)
  distance[!is.finite(distance)] <- Inf
  pmax(distance, 0)
}

.rapid_farthest_seeds <- function(view, candidates, k, canonical_rank) {
  candidates <- unique(as.integer(candidates))
  if (!length(candidates) || k <= 0L) return(integer(0))
  k <- min(as.integer(k), length(candidates))
  local <- view[candidates, , drop = FALSE]
  center <- colMeans(local)
  first_distance <- rowSums((local - rep(center, each = nrow(local)))^2)
  first_pool <- which(first_distance == min(first_distance))
  first <- first_pool[which.min(canonical_rank[candidates[first_pool]])]
  selected <- first

  while (length(selected) < k) {
    distance <- .rapid_sqdist_to_centers(local, local[selected, , drop = FALSE])
    nearest <- apply(distance, 1L, min)
    nearest[selected] <- -Inf
    pool <- which(nearest == max(nearest))
    if (!length(pool) || !is.finite(nearest[pool[[1L]]])) {
      pool <- setdiff(seq_along(candidates), selected)
    }
    next_index <- pool[which.min(canonical_rank[candidates[pool]])]
    selected <- c(selected, next_index)
  }
  candidates[selected]
}

.rapid_assign_seeded_group <- function(view, candidates, seeds) {
  distance <- .rapid_sqdist_to_centers(
    view[candidates, , drop = FALSE], view[seeds, , drop = FALSE]
  )
  assignment <- max.col(-distance, ties.method = "first")
  seed_match <- match(seeds, candidates)
  assignment[seed_match] <- seq_along(seeds)
  assignment
}

.rapid_group_centers <- function(values, candidates, assignment, k) {
  if (!ncol(values)) return(matrix(numeric(0), k, 0L))
  out <- matrix(0, k, ncol(values))
  for (j in seq_len(k)) {
    members <- candidates[assignment == j]
    out[j, ] <- colMeans(values[members, , drop = FALSE])
  }
  colnames(out) <- colnames(values)
  out
}

.rapid_balanced_prototype_prior <- function(groups) {
  group_levels <- unique(groups)
  prior <- numeric(length(groups))
  for (group in group_levels) {
    members <- which(groups == group)
    prior[members] <- 1 / length(group_levels) / length(members)
  }
  prior / sum(prior)
}

# Initialize a deterministic shared prototype bank across domains.
.rapid_initialize_prototypes <- function(
  embeddings,
  labels = NULL,
  structures = NULL,
  attributes = NULL,
  control = NULL
) {
  observations <- .rapid_prototype_observations(
    embeddings, labels, structures, attributes, control
  )
  control <- observations$control
  label <- observations$label
  known_classes <- sort(unique(label[!is.na(label)]), method = "radix")
  canonical_rank <- integer(length(label))
  canonical_rank[observations$canonical_order] <- seq_along(label)

  group_candidates <- lapply(known_classes, function(class) which(label == class))
  names(group_candidates) <- known_classes
  unknown_candidates <- which(is.na(label))
  if (length(unknown_candidates) &&
      (control$unknown_prototypes > 0L || !length(known_classes))) {
    group_candidates[[control$unknown_level]] <- unknown_candidates
  }
  if (!length(group_candidates)) {
    group_candidates[[control$unknown_level]] <- seq_along(label)
  }

  seed_groups <- list()
  for (group in names(group_candidates)) {
    requested <- if (identical(group, control$unknown_level)) {
      max(control$unknown_prototypes, 1L)
    } else {
      control$prototypes_per_class
    }
    seed_groups[[group]] <- .rapid_farthest_seeds(
      observations$selection_view,
      group_candidates[[group]],
      requested,
      canonical_rank
    )
  }
  seed_groups <- seed_groups[lengths(seed_groups) > 0L]
  if (!length(seed_groups)) stop("No prototype could be initialized.", call. = FALSE)

  K <- sum(lengths(seed_groups))
  embedding <- matrix(0, K, ncol(observations$embedding))
  structure <- matrix(0, K, ncol(observations$structure))
  attribute <- matrix(0, K, ncol(observations$attribute))
  prototype_class <- rep(names(seed_groups), lengths(seed_groups))
  seed_index <- unlist(seed_groups, use.names = FALSE)
  global_assignment <- rep(NA_integer_, length(label))
  offset <- 0L

  for (group in names(seed_groups)) {
    candidates <- group_candidates[[group]]
    seeds <- seed_groups[[group]]
    assignment <- .rapid_assign_seeded_group(
      observations$selection_view, candidates, seeds
    )
    local_k <- length(seeds)
    rows <- offset + seq_len(local_k)
    embedding[rows, ] <- .rapid_group_centers(
      observations$embedding, candidates, assignment, local_k
    )
    if (ncol(structure)) {
      structure[rows, ] <- .rapid_group_centers(
        observations$structure, candidates, assignment, local_k
      )
    }
    if (ncol(attribute)) {
      attribute[rows, ] <- .rapid_group_centers(
        observations$attribute, candidates, assignment, local_k
      )
    }
    global_assignment[candidates] <- offset + assignment
    offset <- offset + local_k
  }
  colnames(embedding) <- colnames(observations$embedding)
  colnames(structure) <- colnames(observations$structure)
  colnames(attribute) <- colnames(observations$attribute)

  class_levels <- c(known_classes, control$unknown_level)
  class_prob <- matrix(
    0, nrow = K, ncol = length(class_levels),
    dimnames = list(NULL, class_levels)
  )
  class_prob[cbind(seq_len(K), match(prototype_class, class_levels))] <- 1
  fixed_class <- prototype_class != control$unknown_level
  prior <- .rapid_balanced_prototype_prior(prototype_class)
  assignments <- split(global_assignment, observations$domain_index)
  names(assignments) <- names(if (is.list(embeddings)) embeddings else list(domain1 = embeddings))

  structure(
    list(
      embedding = embedding,
      structure = structure,
      attribute = attribute,
      class_prob = class_prob,
      class_levels = class_levels,
      prototype_class = prototype_class,
      fixed_class = fixed_class,
      prior = prior,
      mass = prior,
      active = rep(TRUE, K),
      recovered = rep(FALSE, K),
      seed_index = as.integer(seed_index),
      assignments = assignments,
      control = control,
      diagnostics = list(
        n_observations = as.integer(length(label)),
        n_domains = as.integer(length(observations$n_by_domain)),
        n_prototypes = as.integer(K),
        known_classes = known_classes,
        class_counts = table(factor(label, levels = known_classes)),
        prototype_counts = table(prototype_class),
        class_prior_mass = tapply(prior, prototype_class, sum),
        unknown_only = !length(known_classes),
        fully_supervised = !anyNA(label),
        seed = control$seed
      )
    ),
    class = "rapid_ma_prototypes"
  )
}

.rapid_as_couplings <- function(couplings, n_by_domain, K) {
  if (is.matrix(couplings) || inherits(couplings, "Matrix")) couplings <- list(couplings)
  if (!is.list(couplings) || length(couplings) != length(n_by_domain)) {
    stop("`couplings` must provide one matrix per domain.", call. = FALSE)
  }
  lapply(seq_along(couplings), function(i) {
    Q <- couplings[[i]]
    .rapid_validate_numeric_matrix(
      Q, paste0("couplings[[", i, "]]"), n_rows = n_by_domain[[i]]
    )
    if (ncol(Q) != K) {
      stop("Every coupling must have one column per prototype.", call. = FALSE)
    }
    values <- if (inherits(Q, "Matrix")) Q@x else Q
    if (any(!is.finite(values)) || any(values < 0)) {
      stop("Couplings must contain finite nonnegative mass.", call. = FALSE)
    }
    Q
  })
}

.rapid_weighted_centers <- function(Q, values, mass, fallback) {
  if (!ncol(values)) return(matrix(numeric(0), ncol(Q), 0L))
  numerator <- as.matrix(Matrix::crossprod(Q, values))
  out <- fallback
  nonempty <- which(mass > 0)
  if (length(nonempty)) {
    out[nonempty, ] <- numerator[nonempty, , drop = FALSE] / mass[nonempty]
  }
  out
}

.rapid_recovery_candidates <- function(prototype_class, observations, unknown_level) {
  if (identical(prototype_class, unknown_level)) {
    candidates <- which(is.na(observations$label))
    if (!length(candidates)) candidates <- seq_along(observations$label)
    candidates
  } else {
    which(observations$label == prototype_class)
  }
}

.rapid_recover_empty_prototypes <- function(state, observations, empty, control) {
  if (!length(empty) || control$max_recoveries == 0L) {
    return(list(state = state, recovered = integer(0), retained = empty))
  }
  recover <- head(empty, control$max_recoveries)
  retained <- setdiff(empty, recover)
  used <- integer(0)
  current_view <- .rapid_weighted_prototype_view(
    state$embedding, state$structure, state$attribute, control
  )
  observation_view <- observations$selection_view

  for (k in recover) {
    candidates <- .rapid_recovery_candidates(
      state$prototype_class[[k]], observations, control$unknown_level
    )
    if (!length(candidates)) {
      retained <- c(retained, k)
      next
    }
    unused <- setdiff(candidates, used)
    if (length(unused)) candidates <- unused
    reference <- setdiff(which(state$mass > control$min_mass), k)
    if (length(reference)) {
      distance <- .rapid_sqdist_to_centers(
        observation_view[candidates, , drop = FALSE],
        current_view[reference, , drop = FALSE]
      )
      score <- apply(distance, 1L, min)
    } else {
      score <- rowSums(observation_view[candidates, , drop = FALSE]^2)
    }
    canonical_rank <- integer(nrow(observation_view))
    canonical_rank[observations$canonical_order] <- seq_len(nrow(observation_view))
    pool <- which(score == max(score))
    chosen <- candidates[pool[which.min(canonical_rank[candidates[pool]])]]
    state$embedding[k, ] <- observations$embedding[chosen, ]
    if (ncol(state$structure)) state$structure[k, ] <- observations$structure[chosen, ]
    if (ncol(state$attribute)) state$attribute[k, ] <- observations$attribute[chosen, ]
    state$mass[[k]] <- max(control$min_mass, state$prior[[k]])
    state$active[[k]] <- TRUE
    state$recovered[[k]] <- TRUE
    used <- c(used, chosen)
    current_view[k, ] <- observation_view[chosen, ]
  }
  recovered <- setdiff(recover, retained)
  list(state = state, recovered = recovered, retained = sort(unique(retained)))
}

# Update prototype summaries from node-to-prototype coupling mass.
.rapid_update_prototypes <- function(
  state,
  couplings,
  embeddings,
  labels = NULL,
  structures = NULL,
  attributes = NULL,
  control = NULL
) {
  if (!inherits(state, "rapid_ma_prototypes")) {
    stop("`state` must come from `.rapid_initialize_prototypes()`.",
         call. = FALSE)
  }
  control <- if (is.null(control)) state$control else {
    .rapid_resolve_prototype_control(control)
  }
  observations <- .rapid_prototype_observations(
    embeddings, labels, structures, attributes, control
  )
  K <- nrow(state$embedding)
  expected_dims <- c(
    embedding = ncol(state$embedding),
    structure = ncol(state$structure),
    attribute = ncol(state$attribute)
  )
  observed_dims <- c(
    embedding = ncol(observations$embedding),
    structure = ncol(observations$structure),
    attribute = ncol(observations$attribute)
  )
  if (!identical(expected_dims, observed_dims)) {
    stop(
      "Update views must match the initialized embedding, structure, and ",
      "attribute dimensions.",
      call. = FALSE
    )
  }
  new_classes <- setdiff(
    unique(observations$label[!is.na(observations$label)]),
    state$class_levels
  )
  if (length(new_classes)) {
    stop("Update labels contain class(es) absent from the prototype bank: ",
         paste(new_classes, collapse = ", "), ".", call. = FALSE)
  }
  Q_domains <- .rapid_as_couplings(
    couplings, observations$n_by_domain, K
  )
  if (any(vapply(Q_domains, inherits, logical(1), "Matrix"))) {
    Q <- do.call(rbind, lapply(Q_domains, Matrix::Matrix, sparse = TRUE))
    raw_mass <- as.numeric(Matrix::colSums(Q))
  } else {
    Q <- do.call(rbind, Q_domains)
    raw_mass <- colSums(Q)
  }

  updated <- state
  updated$embedding <- .rapid_weighted_centers(
    Q, observations$embedding, raw_mass, state$embedding
  )
  updated$structure <- .rapid_weighted_centers(
    Q, observations$structure, raw_mass, state$structure
  )
  updated$attribute <- .rapid_weighted_centers(
    Q, observations$attribute, raw_mass, state$attribute
  )
  updated$mass <- raw_mass
  updated$active <- raw_mass > control$min_mass
  updated$recovered <- rep(FALSE, K)

  label_index <- match(observations$label, state$class_levels)
  observed <- which(!is.na(label_index))
  if (length(observed)) {
    indicator <- Matrix::sparseMatrix(
      i = observed,
      j = label_index[observed],
      x = 1,
      dims = c(nrow(Q), length(state$class_levels))
    )
    evidence <- as.matrix(Matrix::crossprod(Q, indicator))
    for (k in which(!state$fixed_class)) {
      total <- sum(evidence[k, ])
      if (total > control$min_mass) {
        probability <- evidence[k, ] + control$class_smoothing
        updated$class_prob[k, ] <- probability / sum(probability)
      }
    }
  }
  for (k in which(state$fixed_class)) {
    updated$class_prob[k, ] <- 0
    updated$class_prob[k, match(state$prototype_class[[k]], state$class_levels)] <- 1
  }

  empty <- which(raw_mass <= control$min_mass)
  recovery <- .rapid_recover_empty_prototypes(
    updated, observations, empty, control
  )
  updated <- recovery$state
  if (length(recovery$retained)) {
    updated$active[recovery$retained] <- TRUE
    updated$mass[recovery$retained] <- pmax(
      control$min_mass, updated$prior[recovery$retained]
    )
  }
  updated$control <- control
  updated$diagnostics$raw_mass <- raw_mass
  updated$diagnostics$empty_before_recovery <- as.integer(empty)
  updated$diagnostics$recovered <- as.integer(recovery$recovered)
  updated$diagnostics$retained_previous <- as.integer(recovery$retained)
  updated$diagnostics$n_recovered <- as.integer(length(recovery$recovered))
  updated$diagnostics$n_retained <- as.integer(length(recovery$retained))
  updated$diagnostics$total_coupling_mass <- as.numeric(sum(raw_mass))
  updated
}
