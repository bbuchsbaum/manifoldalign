# RAPID-MA prototype-level multi-relation structure matching --------------

.rapid_relmatch_control <- function(
  normalization = c("mass", "frobenius", "symmetric"),
  storage = c("auto", "sparse", "dense"),
  sparse_threshold = 0.25,
  relation_gate = c("robust", "uniform", "fixed"),
  relation_weights = NULL,
  gate_temperature = 5,
  min_mass = 1e-12,
  dense_coupling_max_entries = 1e6
) {
  normalization <- match.arg(normalization)
  storage <- match.arg(storage)
  relation_gate <- match.arg(relation_gate)
  sparse_threshold <- .rapid_number_scalar(
    sparse_threshold, "sparse_threshold", lower = 0
  )
  if (sparse_threshold > 1) {
    stop("`sparse_threshold` must be <= 1.", call. = FALSE)
  }
  gate_temperature <- .rapid_number_scalar(
    gate_temperature, "gate_temperature", lower = 0
  )
  min_mass <- .rapid_number_scalar(
    min_mass, "min_mass", lower = 0, strict = TRUE
  )
  dense_coupling_max_entries <- .rapid_number_scalar(
    dense_coupling_max_entries,
    "dense_coupling_max_entries",
    lower = 0,
    strict = TRUE
  )
  if (!is.null(relation_weights)) {
    if (!is.numeric(relation_weights) || any(!is.finite(relation_weights)) ||
        any(relation_weights < 0) || is.null(names(relation_weights)) ||
        any(!nzchar(names(relation_weights))) || anyDuplicated(names(relation_weights))) {
      stop("`relation_weights` must be finite, nonnegative, uniquely named values.",
           call. = FALSE)
    }
  }
  if (identical(relation_gate, "fixed") && is.null(relation_weights)) {
    stop("`relation_gate = \"fixed\"` requires `relation_weights`.",
         call. = FALSE)
  }

  structure(
    list(
      normalization = normalization,
      storage = storage,
      sparse_threshold = sparse_threshold,
      relation_gate = relation_gate,
      relation_weights = relation_weights,
      gate_temperature = gate_temperature,
      min_mass = min_mass,
      dense_coupling_max_entries = dense_coupling_max_entries
    ),
    class = "rapid_ma_relmatch_control"
  )
}

.rapid_resolve_relmatch_control <- function(control = NULL) {
  defaults <- unclass(.rapid_relmatch_control())
  if (is.null(control)) return(do.call(.rapid_relmatch_control, defaults))
  if (inherits(control, "rapid_ma_relmatch_control")) return(control)
  if (!is.list(control) || is.null(names(control)) || any(!nzchar(names(control)))) {
    stop("`control` must be NULL, a named list, or rapid_ma_relmatch_control.",
         call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    stop("Unknown relation-matching control field(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  do.call(.rapid_relmatch_control, utils::modifyList(defaults, control))
}

.rapid_relmatch_coupling <- function(coupling, n, control) {
  if (inherits(coupling, "rapid_ma_uot")) coupling <- coupling$coupling
  .rapid_validate_numeric_matrix(coupling, "coupling", n_rows = n)
  values <- if (inherits(coupling, "Matrix")) coupling@x else coupling
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("`coupling` must contain finite nonnegative mass.", call. = FALSE)
  }
  if (!inherits(coupling, "Matrix")) {
    entries <- as.double(nrow(coupling)) * as.double(ncol(coupling))
    if (entries > control$dense_coupling_max_entries) {
      stop(
        "Dense coupling exceeds `dense_coupling_max_entries`; pass a sparse ",
        "node-to-prototype coupling.",
        call. = FALSE
      )
    }
  }
  .rapid_force_dgC(Matrix::Matrix(coupling, sparse = TRUE))
}

.rapid_normalize_prototype_relation <- function(G, normalization, min_mass) {
  G <- .rapid_force_dgC(G)
  total <- sum(G@x)
  frobenius <- sqrt(sum(G@x * G@x))
  if (!length(G@x) || total <= min_mass || frobenius <= min_mass) {
    return(list(
      graph = .rapid_empty_relation(nrow(G)),
      scale = 0,
      total = as.numeric(total),
      frobenius = as.numeric(frobenius),
      empty = TRUE
    ))
  }
  if (identical(normalization, "mass")) {
    scale <- total
    out <- G / scale
  } else if (identical(normalization, "frobenius")) {
    scale <- frobenius
    out <- G / scale
  } else {
    out_strength <- as.numeric(Matrix::rowSums(G))
    in_strength <- as.numeric(Matrix::colSums(G))
    left <- ifelse(out_strength > min_mass, 1 / sqrt(out_strength), 0)
    right <- ifelse(in_strength > min_mass, 1 / sqrt(in_strength), 0)
    out <- Matrix::Diagonal(x = left) %*% G %*% Matrix::Diagonal(x = right)
    norm <- sqrt(sum(out@x * out@x))
    scale <- norm
    if (norm > min_mass) out <- out / norm
  }
  list(
    graph = .rapid_force_dgC(Matrix::drop0(out)),
    scale = as.numeric(scale),
    total = as.numeric(total),
    frobenius = as.numeric(frobenius),
    empty = FALSE
  )
}

.rapid_store_prototype_relation <- function(G, storage, sparse_threshold) {
  K <- nrow(G)
  density <- Matrix::nnzero(G) / max(as.double(K) * K, 1)
  use_sparse <- identical(storage, "sparse") ||
    (identical(storage, "auto") && density <= sparse_threshold)
  if (use_sparse) .rapid_force_dgC(G) else as.matrix(G)
}

# Compute normalized Q' A Q using sparse node/prototype intermediates only.
.rapid_prototype_relation <- function(A, coupling, control = NULL) {
  control <- .rapid_resolve_relmatch_control(control)
  .rapid_validate_numeric_matrix(A, "A")
  if (nrow(A) != ncol(A)) stop("`A` must be square.", call. = FALSE)
  A <- .rapid_force_dgC(A)
  if (any(!is.finite(A@x)) || any(A@x < 0)) {
    stop("`A` must contain finite nonnegative relation weights.", call. = FALSE)
  }
  Q <- .rapid_relmatch_coupling(coupling, nrow(A), control)
  K <- ncol(Q)
  row_support <- tabulate(Matrix::summary(Q)$i, nbins = nrow(Q))
  q_max <- if (length(row_support)) max(row_support) else 0L
  coupling_mass <- sum(Q@x)
  relation_nnz <- Matrix::nnzero(A)

  if (!relation_nnz) {
    graph <- .rapid_empty_relation(K)
    return(structure(
      list(
        graph = .rapid_store_prototype_relation(
          graph, control$storage, control$sparse_threshold
        ),
        raw = graph,
        diagnostics = list(
          status = "empty_relation",
          n_nodes = as.integer(nrow(A)),
          n_prototypes = as.integer(K),
          relation_nnz = 0L,
          coupling_nnz = as.integer(Matrix::nnzero(Q)),
          coupling_mass = as.numeric(coupling_mass),
          unmatched_relation_nodes = 0L,
          q_max = as.integer(q_max),
          operation_bound = 0,
          intermediate_nnz = 0L,
          graph_nnz = 0L,
          dense_node_square_allocated = FALSE
        )
      ),
      class = "rapid_ma_prototype_relation"
    ))
  }
  if (coupling_mass <= control$min_mass || Matrix::nnzero(Q) == 0L) {
    graph <- .rapid_empty_relation(K)
    return(structure(
      list(
        graph = .rapid_store_prototype_relation(
          graph, control$storage, control$sparse_threshold
        ),
        raw = graph,
        diagnostics = list(
          status = "unmatched_mass",
          n_nodes = as.integer(nrow(A)),
          n_prototypes = as.integer(K),
          relation_nnz = as.integer(relation_nnz),
          coupling_nnz = 0L,
          coupling_mass = as.numeric(coupling_mass),
          unmatched_relation_nodes = as.integer(length(unique(c(
            Matrix::summary(A)$i, Matrix::summary(A)$j
          )))),
          q_max = 0L,
          operation_bound = 0,
          intermediate_nnz = 0L,
          graph_nnz = 0L,
          dense_node_square_allocated = FALSE
        )
      ),
      class = "rapid_ma_prototype_relation"
    ))
  }

  AQ <- A %*% Q
  raw <- .rapid_force_dgC(Matrix::crossprod(Q, AQ))
  normalized <- .rapid_normalize_prototype_relation(
    raw, control$normalization, control$min_mass
  )
  incident <- unique(c(Matrix::summary(A)$i, Matrix::summary(A)$j))
  unmatched <- sum(row_support[incident] == 0L)
  status <- if (normalized$empty) "zero_projected_mass" else "ok"

  structure(
    list(
      graph = .rapid_store_prototype_relation(
        normalized$graph, control$storage, control$sparse_threshold
      ),
      raw = raw,
      diagnostics = list(
        status = status,
        n_nodes = as.integer(nrow(A)),
        n_prototypes = as.integer(K),
        relation_nnz = as.integer(relation_nnz),
        coupling_nnz = as.integer(Matrix::nnzero(Q)),
        coupling_mass = as.numeric(coupling_mass),
        unmatched_relation_nodes = as.integer(unmatched),
        q_max = as.integer(q_max),
        operation_bound = as.double(relation_nnz) * as.double(q_max)^2,
        intermediate_nnz = as.integer(Matrix::nnzero(AQ)),
        raw_graph_mass = normalized$total,
        normalization_scale = normalized$scale,
        graph_nnz = as.integer(Matrix::nnzero(normalized$graph)),
        graph_density = as.numeric(
          Matrix::nnzero(normalized$graph) / max(as.double(K) * K, 1)
        ),
        dense_node_square_allocated = FALSE
      )
    ),
    class = "rapid_ma_prototype_relation"
  )
}

.rapid_relation_list <- function(state, name) {
  if (inherits(state, "rapid_ma_relations")) return(state$relations)
  if (!is.list(state) || is.null(names(state)) || any(!nzchar(names(state)))) {
    stop("`", name, "` must be a rapid_ma_relations object or named relation list.",
         call. = FALSE)
  }
  state
}

.rapid_domain_relation_weight <- function(weights, domain, relation) {
  if (is.null(weights)) return(1)
  current <- if (is.list(weights)) weights[[domain]] else weights
  if (is.null(current) || is.null(names(current)) || is.null(current[[relation]])) {
    return(1)
  }
  value <- current[[relation]]
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < 0) {
    stop("Domain relation weights must be finite and nonnegative.", call. = FALSE)
  }
  as.numeric(value)
}

.rapid_weighted_graph_mean <- function(graphs, weights, K) {
  available <- which(weights > 0 & !vapply(graphs, is.null, logical(1)))
  if (!length(available)) return(.rapid_empty_relation(K))
  weights <- weights[available]
  weights <- weights / sum(weights)
  out <- .rapid_empty_relation(K)
  for (i in seq_along(available)) {
    out <- out + weights[[i]] * Matrix::Matrix(
      graphs[[available[[i]]]], sparse = TRUE
    )
  }
  .rapid_force_dgC(Matrix::drop0(out))
}

.rapid_graph_squared_loss <- function(A, B, min_mass) {
  difference <- Matrix::Matrix(A, sparse = TRUE) - Matrix::Matrix(B, sparse = TRUE)
  numerator <- sum(difference@x * difference@x)
  denominator <- max(sum(Matrix::Matrix(B, sparse = TRUE)@x^2), min_mass)
  as.numeric(numerator / denominator)
}

# Project each domain relation through its sparse coupling and fit barycenters.
.rapid_relmatch_fit <- function(
  relation_states,
  couplings,
  domain_weights = NULL,
  domain_relation_weights = NULL,
  control = NULL
) {
  control <- .rapid_resolve_relmatch_control(control)
  if (!is.list(relation_states) || !length(relation_states)) {
    stop("`relation_states` must be a non-empty list of domains.", call. = FALSE)
  }
  if (!is.list(couplings) || length(couplings) != length(relation_states)) {
    stop("`couplings` must provide one coupling per domain.", call. = FALSE)
  }
  M <- length(relation_states)
  domain_names <- names(relation_states)
  if (is.null(domain_names) || any(!nzchar(domain_names))) {
    domain_names <- paste0("domain", seq_len(M))
  }
  relation_lists <- lapply(seq_len(M), function(i) {
    .rapid_relation_list(relation_states[[i]], paste0("relation_states[[", i, "]]"))
  })
  names(relation_lists) <- domain_names
  relation_names <- sort(unique(unlist(lapply(relation_lists, names))), method = "radix")
  if (!length(relation_names)) stop("No named relations were supplied.", call. = FALSE)

  if (is.null(domain_weights)) domain_weights <- rep(1, M)
  if (length(domain_weights) != M || !is.numeric(domain_weights) ||
      any(!is.finite(domain_weights)) || any(domain_weights < 0) ||
      sum(domain_weights) <= 0) {
    stop("`domain_weights` must be finite nonnegative values with positive total.",
         call. = FALSE)
  }
  domain_weights <- as.numeric(domain_weights / sum(domain_weights))

  coupling_matrices <- vector("list", M)
  K <- NULL
  for (m in seq_len(M)) {
    n <- if (length(relation_lists[[m]])) nrow(relation_lists[[m]][[1L]]) else {
      if (inherits(couplings[[m]], "rapid_ma_uot")) {
        nrow(couplings[[m]]$coupling)
      } else {
        nrow(couplings[[m]])
      }
    }
    coupling_matrices[[m]] <- .rapid_relmatch_coupling(couplings[[m]], n, control)
    if (is.null(K)) K <- ncol(coupling_matrices[[m]])
    if (ncol(coupling_matrices[[m]]) != K) {
      stop("All domains must share the same prototype columns.", call. = FALSE)
    }
  }

  projected <- setNames(vector("list", M), domain_names)
  graph_lists <- setNames(vector("list", length(relation_names)), relation_names)
  contribution_lists <- setNames(vector("list", length(relation_names)), relation_names)
  for (relation in relation_names) {
    graph_lists[[relation]] <- vector("list", M)
    contribution_lists[[relation]] <- numeric(M)
  }
  for (m in seq_len(M)) {
    projected[[m]] <- setNames(vector("list", length(relation_names)), relation_names)
    for (relation in relation_names) {
      A <- relation_lists[[m]][[relation]]
      if (is.null(A)) {
        projected[[m]][[relation]] <- list(
          graph = NULL,
          raw = NULL,
          diagnostics = list(status = "missing_relation")
        )
        next
      }
      fit <- .rapid_prototype_relation(A, coupling_matrices[[m]], control)
      projected[[m]][[relation]] <- fit
      if (identical(fit$diagnostics$status, "ok")) {
        graph_lists[[relation]][[m]] <- fit$graph
        contribution_lists[[relation]][[m]] <- domain_weights[[m]] *
          .rapid_domain_relation_weight(
            domain_relation_weights, domain_names[[m]], relation
          )
      }
    }
  }

  barycenters <- setNames(vector("list", length(relation_names)), relation_names)
  losses <- matrix(
    NA_real_, M, length(relation_names),
    dimnames = list(domain_names, relation_names)
  )
  valid_relation <- logical(length(relation_names))
  names(valid_relation) <- relation_names
  for (relation in relation_names) {
    barycenter <- .rapid_weighted_graph_mean(
      graph_lists[[relation]], contribution_lists[[relation]], K
    )
    barycenters[[relation]] <- .rapid_store_prototype_relation(
      barycenter, control$storage, control$sparse_threshold
    )
    valid <- which(contribution_lists[[relation]] > 0)
    valid_relation[[relation]] <- length(valid) > 0L
    for (m in valid) {
      losses[m, relation] <- .rapid_graph_squared_loss(
        graph_lists[[relation]][[m]], barycenter, control$min_mass
      )
    }
  }

  fixed <- control$relation_weights
  if (!is.null(fixed)) {
    unknown <- setdiff(names(fixed), relation_names)
    if (length(unknown)) {
      stop("Unknown fixed relation weight(s): ", paste(unknown, collapse = ", "),
           ".", call. = FALSE)
    }
  }
  raw_gate <- setNames(numeric(length(relation_names)), relation_names)
  for (relation in relation_names) {
    observed_loss <- losses[, relation]
    observed_loss <- observed_loss[is.finite(observed_loss)]
    if (!valid_relation[[relation]] || !length(observed_loss)) next
    reliability <- if (identical(control$relation_gate, "robust")) {
      exp(pmax(-control$gate_temperature * stats::median(observed_loss), -700))
    } else {
      1
    }
    raw_gate[[relation]] <- reliability
  }
  if (identical(control$relation_gate, "fixed")) {
    raw_gate[] <- 0
    raw_gate[names(fixed)] <- fixed
  } else if (!is.null(fixed)) {
    raw_gate[names(fixed)] <- raw_gate[names(fixed)] * fixed
  }
  raw_gate[!valid_relation] <- 0
  relation_weights <- if (sum(raw_gate) > 0) raw_gate / sum(raw_gate) else raw_gate
  domain_loss <- rowSums(
    sweep(replace(losses, !is.finite(losses), 0), 2L, relation_weights, "*"),
    na.rm = TRUE
  )
  total_loss <- sum(domain_weights * domain_loss)

  structure(
    list(
      projected = projected,
      barycenters = barycenters,
      losses = losses,
      domain_loss = domain_loss,
      total_loss = as.numeric(total_loss),
      relation_weights = relation_weights,
      raw_relation_gate = raw_gate,
      domain_weights = domain_weights,
      diagnostics = list(
        n_domains = as.integer(M),
        n_prototypes = as.integer(K),
        relation_names = relation_names,
        valid_relations = names(valid_relation)[valid_relation],
        relation_status = lapply(projected, function(domain) {
          vapply(domain, function(x) x$diagnostics$status, character(1))
        }),
        total_relation_nnz = as.integer(sum(vapply(
          relation_lists,
          function(domain) sum(vapply(domain, Matrix::nnzero, numeric(1))),
          numeric(1)
        ))),
        dense_node_square_allocated = FALSE
      ),
      control = control
    ),
    class = "rapid_ma_relmatch"
  )
}
