# RAPID-MA sparse node-to-prototype unbalanced transport ------------------

.rapid_uot_control <- function(
  q = 8L,
  epsilon = 0.1,
  rho_source = 1,
  rho_target = 1,
  latent_weight = 1,
  structure_weight = 0.5,
  attribute_weight = 0.25,
  position_weight = 0.25,
  label_weight = 2,
  label_confidence_threshold = 0.8,
  hard_labels = FALSE,
  anchor_weight = 5,
  hard_anchors = FALSE,
  block_size = 512L,
  ensure_prototype_coverage = TRUE,
  cost_cap = 1e6,
  max_iter = 1000L,
  tol = 1e-7,
  require_convergence = TRUE
) {
  q <- .rapid_int_scalar(q, "q")
  block_size <- .rapid_int_scalar(block_size, "block_size")
  max_iter <- .rapid_int_scalar(max_iter, "max_iter")
  epsilon <- .rapid_number_scalar(epsilon, "epsilon", lower = 0, strict = TRUE)
  rho_source <- .rapid_number_scalar(
    rho_source, "rho_source", lower = 0, strict = TRUE
  )
  rho_target <- .rapid_number_scalar(
    rho_target, "rho_target", lower = 0, strict = TRUE
  )
  tol <- .rapid_number_scalar(tol, "tol", lower = 0, strict = TRUE)
  cost_cap <- .rapid_number_scalar(cost_cap, "cost_cap", lower = 0, strict = TRUE)
  label_confidence_threshold <- .rapid_number_scalar(
    label_confidence_threshold,
    "label_confidence_threshold",
    lower = 0
  )
  if (label_confidence_threshold > 1) {
    stop("`label_confidence_threshold` must be <= 1.", call. = FALSE)
  }
  weights <- c(
    latent_weight = latent_weight,
    structure_weight = structure_weight,
    attribute_weight = attribute_weight,
    position_weight = position_weight,
    label_weight = label_weight,
    anchor_weight = anchor_weight
  )
  for (name in names(weights)) {
    weights[[name]] <- .rapid_number_scalar(weights[[name]], name, lower = 0)
  }
  flags <- list(
    hard_labels = hard_labels,
    hard_anchors = hard_anchors,
    ensure_prototype_coverage = ensure_prototype_coverage,
    require_convergence = require_convergence
  )
  for (name in names(flags)) {
    if (length(flags[[name]]) != 1L || !is.logical(flags[[name]]) ||
        is.na(flags[[name]])) {
      stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (weights[["latent_weight"]] + weights[["structure_weight"]] +
      weights[["attribute_weight"]] + weights[["position_weight"]] <= 0) {
    stop("At least one numeric UOT view weight must be positive.", call. = FALSE)
  }

  structure(
    c(
      list(
        q = q,
        epsilon = epsilon,
        rho_source = rho_source,
        rho_target = rho_target,
        label_confidence_threshold = label_confidence_threshold,
        hard_labels = hard_labels,
        hard_anchors = hard_anchors,
        block_size = block_size,
        ensure_prototype_coverage = ensure_prototype_coverage,
        cost_cap = cost_cap,
        max_iter = max_iter,
        tol = tol,
        require_convergence = require_convergence
      ),
      as.list(weights)
    ),
    class = "rapid_ma_uot_control"
  )
}

.rapid_resolve_uot_control <- function(control = NULL) {
  defaults <- unclass(.rapid_uot_control())
  if (is.null(control)) return(do.call(.rapid_uot_control, defaults))
  if (inherits(control, "rapid_ma_uot_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`control` must be NULL, a named list, or rapid_ma_uot_control.",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown UOT control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  do.call(.rapid_uot_control, utils::modifyList(defaults, control))
}

.rapid_uot_numeric_view <- function(node_view, prototype_view, n, K, name,
                                    weight, required = FALSE) {
  node_missing <- is.null(node_view)
  prototype_missing <- is.null(prototype_view)
  if (node_missing || prototype_missing) {
    if (isTRUE(required) || xor(node_missing, prototype_missing)) {
      stop("Both node and prototype `", name, "` views are required together.",
           call. = FALSE)
    }
    return(NULL)
  }
  .rapid_validate_numeric_matrix(node_view, paste0("node_", name), n_rows = n)
  .rapid_validate_numeric_matrix(
    prototype_view, paste0("prototype_", name), n_rows = K
  )
  node_view <- as.matrix(node_view)
  prototype_view <- as.matrix(prototype_view)
  if (ncol(node_view) != ncol(prototype_view)) {
    stop("Node and prototype `", name, "` dimensions differ.", call. = FALSE)
  }
  if (any(!is.finite(node_view)) || any(!is.finite(prototype_view))) {
    stop("`", name, "` views must contain only finite values.", call. = FALSE)
  }
  if (weight <= 0) return(NULL)
  list(node = node_view, prototype = prototype_view, weight = weight)
}

.rapid_prepare_uot_views <- function(
  embedding,
  prototypes,
  structure,
  attributes,
  positions,
  prototype_positions,
  control
) {
  if (!inherits(prototypes, "rapid_ma_prototypes")) {
    stop("`prototypes` must come from `.rapid_initialize_prototypes()`.",
         call. = FALSE)
  }
  .rapid_validate_numeric_matrix(embedding, "embedding")
  embedding <- as.matrix(embedding)
  n <- nrow(embedding)
  K <- nrow(prototypes$embedding)
  if (!K) stop("The prototype bank is empty.", call. = FALSE)

  prototype_structure <- prototypes$structure
  if (!ncol(prototype_structure)) prototype_structure <- NULL
  prototype_attribute <- prototypes$attribute
  if (!ncol(prototype_attribute)) prototype_attribute <- NULL
  views <- list(
    latent = .rapid_uot_numeric_view(
      embedding, prototypes$embedding, n, K, "latent",
      control$latent_weight, required = TRUE
    ),
    structure = .rapid_uot_numeric_view(
      structure, prototype_structure, n, K, "structure",
      control$structure_weight
    ),
    attribute = .rapid_uot_numeric_view(
      attributes, prototype_attribute, n, K, "attribute",
      control$attribute_weight
    ),
    position = .rapid_uot_numeric_view(
      positions, prototype_positions, n, K, "position",
      control$position_weight
    )
  )
  views <- views[!vapply(views, is.null, logical(1))]
  if (!length(views)) stop("No active UOT cost views remain.", call. = FALSE)
  list(views = views, n = n, K = K)
}

.rapid_prepare_uot_constraints <- function(labels, label_confidence, anchors,
                                            prototypes, n, control) {
  if (is.null(labels)) labels <- rep(NA_character_, n)
  if (length(labels) != n) stop("`labels` must have one value per node.", call. = FALSE)
  labels <- as.character(labels)
  label_index <- match(labels, prototypes$class_levels)
  unknown_labels <- unique(labels[!is.na(labels) & is.na(label_index)])
  if (length(unknown_labels)) {
    stop("Labels absent from the prototype bank: ",
         paste(unknown_labels, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(label_confidence)) {
    label_confidence <- ifelse(is.na(labels), 0, 1)
  }
  if (length(label_confidence) != n || !is.numeric(label_confidence) ||
      any(!is.finite(label_confidence)) || any(label_confidence < 0) ||
      any(label_confidence > 1)) {
    stop("`label_confidence` must contain one finite value in [0, 1] per node.",
         call. = FALSE)
  }
  confident <- !is.na(label_index) &
    label_confidence >= control$label_confidence_threshold

  if (is.null(anchors)) anchors <- rep(NA_integer_, n)
  if (length(anchors) != n || (!is.numeric(anchors) && !is.integer(anchors))) {
    stop("`anchors` must contain one prototype index or NA per node.",
         call. = FALSE)
  }
  bad_anchor <- !is.na(anchors) &
    (!is.finite(anchors) | anchors < 1 | anchors > nrow(prototypes$embedding) |
       abs(anchors - round(anchors)) > 1e-8)
  if (any(bad_anchor)) {
    stop("`anchors` must contain only NA or valid integer prototype indices.",
         call. = FALSE)
  }
  anchors <- as.integer(anchors)

  list(
    labels = labels,
    label_index = label_index,
    label_confidence = as.numeric(label_confidence),
    confident = confident,
    anchors = anchors
  )
}

.rapid_mean_squared_distance <- function(X, Y) {
  distance <- outer(rowSums(X * X), rowSums(Y * Y), "+") - 2 * tcrossprod(X, Y)
  distance <- pmax(distance, 0)
  distance / max(ncol(X), 1L)
}

.rapid_uot_cost_block <- function(rows, prototype_index, views, constraints,
                                  prototypes, control) {
  cost <- matrix(0, length(rows), length(prototype_index))
  for (view in views) {
    cost <- cost + view$weight * .rapid_mean_squared_distance(
      view$node[rows, , drop = FALSE],
      view$prototype[prototype_index, , drop = FALSE]
    )
  }

  confident <- constraints$confident[rows]
  if (any(confident) && (control$label_weight > 0 || control$hard_labels)) {
    local_rows <- which(confident)
    class_columns <- constraints$label_index[rows[local_rows]]
    probability <- matrix(0, length(local_rows), length(prototype_index))
    for (i in seq_along(local_rows)) {
      probability[i, ] <- prototypes$class_prob[
        prototype_index, class_columns[[i]]
      ]
    }
    if (control$label_weight > 0) {
      confidence <- constraints$label_confidence[rows[local_rows]]
      cost[local_rows, ] <- cost[local_rows, , drop = FALSE] +
        control$label_weight * confidence * (1 - probability)
    }
    if (isTRUE(control$hard_labels)) {
      invalid <- probability <= 0
      indexed <- which(invalid, arr.ind = TRUE)
      if (nrow(indexed)) {
        cost[cbind(local_rows[indexed[, 1L]], indexed[, 2L])] <- Inf
      }
    }
  }

  anchored <- which(!is.na(constraints$anchors[rows]))
  if (length(anchored) && (control$anchor_weight > 0 || control$hard_anchors)) {
    preferred <- constraints$anchors[rows[anchored]]
    mismatch <- outer(preferred, prototype_index, `!=`)
    if (control$anchor_weight > 0) {
      cost[anchored, ] <- cost[anchored, , drop = FALSE] +
        control$anchor_weight * mismatch
    }
    if (isTRUE(control$hard_anchors)) {
      indexed <- which(mismatch, arr.ind = TRUE)
      if (nrow(indexed)) {
        cost[cbind(anchored[indexed[, 1L]], indexed[, 2L])] <- Inf
      }
    }
  }

  finite <- is.finite(cost)
  cost[finite] <- pmin(pmax(cost[finite], 0), control$cost_cap)
  cost
}

.rapid_edges_to_csr_csc <- function(edge_i, edge_j, edge_cost, n, K) {
  if (!length(edge_i)) stop("The sparse transport support has no edges.", call. = FALSE)
  key <- as.double(edge_i) + as.double(n) * (as.double(edge_j) - 1)
  if (anyDuplicated(key)) {
    ordering <- order(key, edge_cost, method = "radix")
    keep <- !duplicated(key[ordering])
    ordering <- ordering[keep]
    edge_i <- edge_i[ordering]
    edge_j <- edge_j[ordering]
    edge_cost <- edge_cost[ordering]
  }
  if (length(edge_i) > .Machine$integer.max) {
    stop("Sparse transport support exceeds the integer edge limit.",
         call. = FALSE)
  }

  row_order <- order(edge_i, edge_j, method = "radix")
  csr_i <- edge_i[row_order]
  csr_j <- edge_j[row_order]
  csr_cost <- edge_cost[row_order]
  row_counts <- tabulate(csr_i, nbins = n)
  row_ptr <- as.integer(c(0, cumsum(row_counts)))

  col_order <- order(edge_j, edge_i, method = "radix")
  csc_i <- edge_i[col_order]
  csc_j <- edge_j[col_order]
  csc_cost <- edge_cost[col_order]
  col_counts <- tabulate(csc_j, nbins = K)
  col_ptr <- as.integer(c(0, cumsum(col_counts)))

  list(
    row_ptr = row_ptr,
    col_idx = as.integer(csr_j),
    cost = as.numeric(csr_cost),
    col_ptr = col_ptr,
    row_idx = as.integer(csc_i),
    cost_csc = as.numeric(csc_cost),
    n_rows = as.integer(n),
    n_cols = as.integer(K)
  )
}

.rapid_sparse_uot_cost <- function(views, constraints, prototypes,
                                   initial_prototypes, control) {
  n <- nrow(views[[1L]]$node)
  q <- min(control$q, length(initial_prototypes))
  max_edges <- as.double(n) * q + length(initial_prototypes)
  if (max_edges > .Machine$integer.max) {
    stop("Requested `n * q` support exceeds the R integer edge limit.",
         call. = FALSE)
  }
  edge_i <- integer(as.integer(max_edges))
  edge_j <- integer(as.integer(max_edges))
  edge_cost <- numeric(as.integer(max_edges))
  cursor <- 0L
  block_starts <- seq.int(1L, n, by = control$block_size)
  peak_block_entries <- 0

  for (start in block_starts) {
    rows <- start:min(n, start + control$block_size - 1L)
    cost <- .rapid_uot_cost_block(
      rows, initial_prototypes, views, constraints, prototypes, control
    )
    peak_block_entries <- max(peak_block_entries, length(cost))
    available <- rowSums(is.finite(cost))
    if (any(available == 0L)) {
      offending <- rows[available == 0L]
      stop("Transport constraints leave source row(s) without a candidate: ",
           paste(head(offending, 8L), collapse = ", "), ".", call. = FALSE)
    }
    for (pass in seq_len(q)) {
      chosen <- max.col(-cost, ties.method = "first")
      valid <- available >= pass
      if (!any(valid)) next
      selected_rows <- which(valid)
      count <- length(selected_rows)
      target <- cursor + seq_len(count)
      edge_i[target] <- rows[selected_rows]
      edge_j[target] <- initial_prototypes[chosen[selected_rows]]
      edge_cost[target] <- cost[cbind(selected_rows, chosen[selected_rows])]
      cost[cbind(selected_rows, chosen[selected_rows])] <- Inf
      cursor <- cursor + count
    }
  }

  if (cursor < length(edge_i)) {
    edge_i <- edge_i[seq_len(cursor)]
    edge_j <- edge_j[seq_len(cursor)]
    edge_cost <- edge_cost[seq_len(cursor)]
  }
  covered <- tabulate(edge_j, nbins = nrow(prototypes$embedding)) > 0L
  missing <- initial_prototypes[!covered[initial_prototypes]]
  unreachable <- integer(0)

  if (length(missing) && isTRUE(control$ensure_prototype_coverage)) {
    best_cost <- rep(Inf, length(missing))
    best_row <- rep(NA_integer_, length(missing))
    for (start in block_starts) {
      rows <- start:min(n, start + control$block_size - 1L)
      cost <- .rapid_uot_cost_block(
        rows, missing, views, constraints, prototypes, control
      )
      peak_block_entries <- max(peak_block_entries, length(cost))
      for (j in seq_along(missing)) {
        local <- which.min(cost[, j])
        candidate_cost <- cost[local, j]
        if (is.finite(candidate_cost) && candidate_cost < best_cost[[j]]) {
          best_cost[[j]] <- candidate_cost
          best_row[[j]] <- rows[[local]]
        }
      }
    }
    reachable <- which(is.finite(best_cost) & !is.na(best_row))
    if (length(reachable)) {
      edge_i <- c(edge_i, best_row[reachable])
      edge_j <- c(edge_j, missing[reachable])
      edge_cost <- c(edge_cost, best_cost[reachable])
    }
    unreachable <- missing[!seq_along(missing) %in% reachable]
  } else if (length(missing)) {
    unreachable <- missing
  }

  active <- setdiff(initial_prototypes, unreachable)
  keep <- edge_j %in% active
  edge_i <- edge_i[keep]
  edge_cost <- edge_cost[keep]
  global_j <- edge_j[keep]
  edge_j <- match(global_j, active)
  cost <- .rapid_edges_to_csr_csc(
    edge_i, edge_j, edge_cost, n = n, K = length(active)
  )
  local_coverage <- diff(cost$col_ptr)
  if (any(local_coverage == 0L)) {
    stop("Internal error: an active prototype lacks sparse support.",
         call. = FALSE)
  }

  list(
    cost = cost,
    prototype_index = as.integer(active),
    inactive_prototypes = as.integer(setdiff(seq_len(nrow(prototypes$embedding)), active)),
    diagnostics = list(
      mode = "sparse",
      q = as.integer(q),
      stored_edges = as.integer(length(cost$cost)),
      edge_bound = as.integer(max_edges),
      row_candidate_counts = as.integer(diff(cost$row_ptr)),
      prototype_edge_counts = as.integer(local_coverage),
      coverage_edges_added = as.integer(max(0, length(cost$cost) - cursor)),
      peak_block_entries = as.integer(peak_block_entries),
      dense_node_prototype_retained = FALSE
    )
  )
}

.rapid_dense_uot_cost <- function(views, constraints, prototypes,
                                  initial_prototypes, control) {
  n <- nrow(views[[1L]]$node)
  cost <- matrix(Inf, n, length(initial_prototypes))
  for (start in seq.int(1L, n, by = control$block_size)) {
    rows <- start:min(n, start + control$block_size - 1L)
    cost[rows, ] <- .rapid_uot_cost_block(
      rows, initial_prototypes, views, constraints, prototypes, control
    )
  }
  if (any(rowSums(is.finite(cost)) == 0L)) {
    stop("Transport constraints leave at least one source row without a candidate.",
         call. = FALSE)
  }
  reachable <- colSums(is.finite(cost)) > 0L
  active <- initial_prototypes[reachable]
  cost <- cost[, reachable, drop = FALSE]
  list(
    cost = cost,
    prototype_index = as.integer(active),
    inactive_prototypes = as.integer(setdiff(seq_len(nrow(prototypes$embedding)), active)),
    diagnostics = list(
      mode = "dense",
      q = as.integer(length(active)),
      stored_edges = as.integer(sum(is.finite(cost))),
      edge_bound = as.integer(length(cost)),
      row_candidate_counts = as.integer(rowSums(is.finite(cost))),
      prototype_edge_counts = as.integer(colSums(is.finite(cost))),
      coverage_edges_added = 0L,
      peak_block_entries = as.integer(
        min(n, control$block_size) * length(initial_prototypes)
      ),
      dense_node_prototype_retained = TRUE
    )
  )
}

# Build a composite dense or q-sparse node-to-prototype transport cost.
.rapid_uot_build_cost <- function(
  embedding,
  prototypes,
  structure = NULL,
  attributes = NULL,
  positions = NULL,
  prototype_positions = NULL,
  labels = NULL,
  label_confidence = NULL,
  anchors = NULL,
  control = NULL,
  mode = c("sparse", "dense")
) {
  mode <- match.arg(mode)
  control <- .rapid_resolve_uot_control(control)
  prepared <- .rapid_prepare_uot_views(
    embedding, prototypes, structure, attributes, positions,
    prototype_positions, control
  )
  constraints <- .rapid_prepare_uot_constraints(
    labels, label_confidence, anchors, prototypes, prepared$n, control
  )
  initial <- which(prototypes$active & prototypes$prior > 0)
  anchored <- unique(constraints$anchors[!is.na(constraints$anchors)])
  if (length(setdiff(anchored, initial))) {
    stop("An anchor refers to an inactive or zero-prior prototype.",
         call. = FALSE)
  }
  if (!length(initial)) stop("No active positive-prior prototype remains.", call. = FALSE)

  built <- if (identical(mode, "sparse")) {
    .rapid_sparse_uot_cost(
      prepared$views, constraints, prototypes, initial, control
    )
  } else {
    .rapid_dense_uot_cost(
      prepared$views, constraints, prototypes, initial, control
    )
  }
  built$diagnostics$n_nodes <- as.integer(prepared$n)
  built$diagnostics$n_prototypes <- as.integer(prepared$K)
  built$diagnostics$n_transport_prototypes <- as.integer(length(built$prototype_index))
  built$diagnostics$active_views <- names(prepared$views)
  built$diagnostics$weights <- vapply(prepared$views, `[[`, numeric(1), "weight")
  built$control <- control
  class(built) <- "rapid_ma_uot_cost"
  built
}

.rapid_uot_coupling_from_solution <- function(cost_state, solution, alpha, beta,
                                              epsilon, K) {
  cost <- cost_state$cost
  if (is.matrix(cost)) {
    n <- nrow(cost)
    row_index <- rep(seq_len(n), times = ncol(cost))
    col_index <- rep(seq_len(ncol(cost)), each = n)
    edge_cost <- as.numeric(cost)
  } else {
    n <- cost$n_rows
    row_index <- rep(seq_len(n), diff(cost$row_ptr))
    col_index <- cost$col_idx
    edge_cost <- cost$cost
  }
  log_mass <- rep(-Inf, length(edge_cost))
  positive <- alpha[row_index] > 0 & beta[col_index] > 0 & is.finite(edge_cost)
  log_mass[positive] <- log(alpha[row_index[positive]]) +
    log(beta[col_index[positive]]) +
    (solution$fbar[row_index[positive]] + solution$gbar[col_index[positive]] -
       edge_cost[positive]) / epsilon
  values <- numeric(length(log_mass))
  values[positive] <- exp(pmin(pmax(log_mass[positive], -745), 700))
  values[!is.finite(values)] <- 0
  global_col <- cost_state$prototype_index[col_index]
  coupling <- Matrix::drop0(Matrix::sparseMatrix(
    i = row_index,
    j = global_col,
    x = values,
    dims = c(n, K)
  ))
  list(
    coupling = coupling,
    transport_cost = sum(values * ifelse(is.finite(edge_cost), edge_cost, 0)),
    support_mass = values
  )
}

.rapid_same_uot_problem <- function(warm_start, cost_state, alpha, beta, control) {
  inherits(warm_start, "rapid_ma_uot") &&
    identical(warm_start$cost_state$cost, cost_state$cost) &&
    identical(warm_start$cost_state$prototype_index, cost_state$prototype_index) &&
    identical(warm_start$alpha, alpha) &&
    identical(warm_start$beta, beta) &&
    identical(
      warm_start$solver_parameters,
      unname(c(
        epsilon = control$epsilon,
        rho_source = control$rho_source,
        rho_target = control$rho_target,
        max_iter = control$max_iter,
        tol = control$tol
      ))
    )
}

# Solve sparse unbalanced transport and materialize only its bounded support.
.rapid_uot_solve <- function(
  embedding,
  prototypes,
  structure = NULL,
  attributes = NULL,
  positions = NULL,
  prototype_positions = NULL,
  labels = NULL,
  label_confidence = NULL,
  anchors = NULL,
  source_mass = NULL,
  control = NULL,
  mode = c("sparse", "dense"),
  warm_start = NULL
) {
  mode <- match.arg(mode)
  control <- .rapid_resolve_uot_control(control)
  cost_state <- .rapid_uot_build_cost(
    embedding = embedding,
    prototypes = prototypes,
    structure = structure,
    attributes = attributes,
    positions = positions,
    prototype_positions = prototype_positions,
    labels = labels,
    label_confidence = label_confidence,
    anchors = anchors,
    control = control,
    mode = mode
  )
  n <- nrow(embedding)
  if (is.null(source_mass)) source_mass <- rep(1 / n, n)
  if (length(source_mass) != n || !is.numeric(source_mass) ||
      any(!is.finite(source_mass)) || any(source_mass < 0) ||
      sum(source_mass) <= 0) {
    stop("`source_mass` must contain finite nonnegative mass with positive total.",
         call. = FALSE)
  }
  alpha <- as.numeric(source_mass)
  beta <- prototypes$prior[cost_state$prototype_index]
  beta <- as.numeric(beta / sum(beta) * sum(alpha))
  solver_parameters <- unname(c(
    epsilon = control$epsilon,
    rho_source = control$rho_source,
    rho_target = control$rho_target,
    max_iter = control$max_iter,
    tol = control$tol
  ))

  reused <- .rapid_same_uot_problem(
    warm_start, cost_state, alpha, beta, control
  )
  solution <- if (reused) {
    warm_start$solution
  } else {
    uot_ti_sinkhorn_kl(
      cost_state$cost,
      alpha = alpha,
      beta = beta,
      epsilon = control$epsilon,
      rho1 = control$rho_source,
      rho2 = control$rho_target,
      max_iter = control$max_iter,
      tol = control$tol
    )
  }
  finite_potentials <- all(is.finite(solution$fbar)) &&
    all(is.finite(solution$gbar))
  if (!finite_potentials) {
    stop(
      "TI-Sinkhorn produced non-finite potentials; increase `epsilon` or ",
      "`rho_source`/`rho_target`, or relax hard constraints.",
      call. = FALSE
    )
  }
  if (isTRUE(control$require_convergence) && !isTRUE(solution$converged)) {
    stop("TI-Sinkhorn did not converge in `max_iter` iterations (residual ",
         format(solution$residual, digits = 5), ").", call. = FALSE)
  }
  materialized <- .rapid_uot_coupling_from_solution(
    cost_state, solution, alpha, beta, control$epsilon,
    K = nrow(prototypes$embedding)
  )
  row_mass <- as.numeric(Matrix::rowSums(materialized$coupling))
  prototype_mass <- as.numeric(Matrix::colSums(materialized$coupling))
  total_mass <- sum(row_mass)
  dropped_fraction <- max(0, 1 - total_mass / sum(alpha))

  structure(
    list(
      coupling = materialized$coupling,
      cost_state = cost_state,
      solution = solution,
      alpha = alpha,
      beta = beta,
      prototype_index = cost_state$prototype_index,
      inactive_prototypes = cost_state$inactive_prototypes,
      row_mass = row_mass,
      prototype_mass = prototype_mass,
      solver_parameters = solver_parameters,
      diagnostics = utils::modifyList(
        cost_state$diagnostics,
        list(
          backend = if (!is.null(solution$backend)) solution$backend else "dense",
          converged = isTRUE(solution$converged),
          iterations = as.integer(solution$iterations),
          residual = as.numeric(solution$residual),
          total_source_mass = as.numeric(sum(alpha)),
          transported_mass = as.numeric(total_mass),
          dropped_fraction = as.numeric(dropped_fraction),
          transport_cost = as.numeric(materialized$transport_cost),
          warm_start_supplied = !is.null(warm_start),
          warm_start_reused = reused,
          coupling_nnz = as.integer(Matrix::nnzero(materialized$coupling)),
          dense_node_square_allocated = FALSE
        )
      )
    ),
    class = "rapid_ma_uot"
  )
}
