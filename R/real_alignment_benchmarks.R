.real_alignment_protocol_version <- "1.0.0"

.real_benchmark_inst_file <- function(subdir, name) {
  path <- system.file(subdir, name, package = "manifoldalign")
  if (!nzchar(path)) {
    path <- file.path("inst", subdir, name)
  }
  if (!file.exists(path)) {
    stop("Missing real-benchmark file: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.real_sha256_file <- function(path) {
  sha256sum <- Sys.which("sha256sum")
  shasum <- Sys.which("shasum")
  if (nzchar(sha256sum)) {
    output <- system2(sha256sum, path, stdout = TRUE, stderr = TRUE)
  } else if (nzchar(shasum)) {
    output <- system2(
      shasum, c("-a", "256", path), stdout = TRUE, stderr = TRUE
    )
  } else {
    stop("A `sha256sum` or `shasum` executable is required.", call. = FALSE)
  }
  if (!length(output)) {
    stop("SHA-256 command returned no output for ", path, ".", call. = FALSE)
  }
  tolower(strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]])
}

real_alignment_benchmark_registry <- function() {
  path <- .real_benchmark_inst_file(
    "benchmarks", "real_alignment_manifest.csv"
  )
  registry <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(
    "protocol_version", "task_id", "dataset_id", "task_type",
    "source_domain", "target_domain", "stress_regime", "split_id",
    "ci_fixture", "full_runner", "source_version", "license",
    "source_url", "requires_external_data"
  )
  if (!identical(names(registry), required)) {
    stop("The real-benchmark manifest schema is invalid.", call. = FALSE)
  }
  if (anyDuplicated(registry$task_id) ||
      any(registry$protocol_version != .real_alignment_protocol_version)) {
    stop("The real-benchmark task IDs or protocol versions are invalid.",
      call. = FALSE
    )
  }
  attr(registry, "manifest_md5") <- unname(tools::md5sum(path))
  registry
}

real_alignment_artifact_registry <- function() {
  path <- .real_benchmark_inst_file(
    "benchmarks", "real_alignment_artifacts.csv"
  )
  artifacts <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(
    "dataset_id", "artifact_id", "source_revision", "relative_path",
    "bytes", "sha256", "redistribution", "license", "source_url"
  )
  if (!identical(names(artifacts), required) ||
      anyDuplicated(paste(artifacts$dataset_id, artifacts$artifact_id))) {
    stop("The real-benchmark artifact registry is invalid.", call. = FALSE)
  }
  artifacts$bytes <- as.numeric(artifacts$bytes)
  attr(artifacts, "manifest_md5") <- unname(tools::md5sum(path))
  artifacts
}

.real_mpx_parse <- function(path) {
  lines <- trimws(readLines(path, warn = FALSE, encoding = "UTF-8"))
  section <- ""
  actor_attributes <- character()
  actors <- list()
  layers <- list()
  edges <- list()

  for (line in lines) {
    if (!nzchar(line) || startsWith(line, "--")) next
    if (startsWith(line, "#")) {
      section <- toupper(sub("^#+", "", line))
      next
    }
    fields <- trimws(strsplit(line, ",", fixed = TRUE)[[1L]])
    if (identical(section, "ACTOR ATTRIBUTES")) {
      actor_attributes <- c(actor_attributes, fields[[1L]])
    } else if (identical(section, "ACTORS")) {
      values <- rep(NA_character_, length(actor_attributes))
      if (length(fields) > 1L && length(actor_attributes)) {
        ncopy <- min(length(actor_attributes), length(fields) - 1L)
        values[seq_len(ncopy)] <- fields[1L + seq_len(ncopy)]
        values[values == "NA"] <- NA_character_
      }
      actors[[length(actors) + 1L]] <- c(actor = fields[[1L]], values)
    } else if (identical(section, "LAYERS")) {
      layers[[length(layers) + 1L]] <- c(
        layer = fields[[1L]],
        directed = tolower(fields[[2L]]) == "directed"
      )
    } else if (identical(section, "EDGES")) {
      weight <- if (length(fields) >= 4L) {
        suppressWarnings(as.numeric(fields[[4L]]))
      } else {
        1
      }
      if (!is.finite(weight)) weight <- 1
      edges[[length(edges) + 1L]] <- list(
        from_actor = fields[[1L]], to_actor = fields[[2L]],
        layer = fields[[3L]], weight = weight
      )
    }
  }

  actor_names <- c("actor", actor_attributes)
  actor_frame <- if (length(actors)) {
    out <- as.data.frame(do.call(rbind, actors), stringsAsFactors = FALSE)
    names(out) <- actor_names
    out
  } else {
    data.frame(actor = character(), stringsAsFactors = FALSE)
  }
  if (anyDuplicated(actor_frame$actor)) {
    stop("MPX actor identities are duplicated in ", basename(path), ".",
      call. = FALSE
    )
  }

  edge_frame <- if (length(edges)) {
    data.frame(
      from_actor = vapply(edges, `[[`, character(1), "from_actor"),
      to_actor = vapply(edges, `[[`, character(1), "to_actor"),
      layer = vapply(edges, `[[`, character(1), "layer"),
      weight = vapply(edges, `[[`, numeric(1), "weight"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      from_actor = character(), to_actor = character(), layer = character(),
      weight = numeric(), stringsAsFactors = FALSE
    )
  }
  layer_frame <- if (length(layers)) {
    data.frame(
      layer = vapply(layers, `[[`, character(1), "layer"),
      directed = vapply(layers, function(x) identical(x[["directed"]], "TRUE"),
        logical(1)
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      layer = sort(unique(edge_frame$layer)), directed = FALSE,
      stringsAsFactors = FALSE
    )
  }
  edge_frame$directed <- layer_frame$directed[
    match(edge_frame$layer, layer_frame$layer)
  ]
  if (anyNA(edge_frame$directed)) {
    stop("MPX edge references an undeclared layer.", call. = FALSE)
  }

  # The AUCS source lists both orientations of undirected ties. Collapse those
  # records without inventing weights; directed Sampson nominations are kept.
  if (nrow(edge_frame)) {
    a <- ifelse(
      edge_frame$directed,
      edge_frame$from_actor,
      pmin(edge_frame$from_actor, edge_frame$to_actor)
    )
    b <- ifelse(
      edge_frame$directed,
      edge_frame$to_actor,
      pmax(edge_frame$from_actor, edge_frame$to_actor)
    )
    key <- paste(edge_frame$layer, a, b, sep = "\r")
    groups <- split(seq_len(nrow(edge_frame)), key)
    edge_frame <- do.call(rbind, lapply(groups, function(index) {
      first <- index[[1L]]
      data.frame(
        from_actor = a[[first]], to_actor = b[[first]],
        layer = edge_frame$layer[[first]],
        weight = max(edge_frame$weight[index]),
        directed = edge_frame$directed[[first]],
        source_records = length(index),
        stringsAsFactors = FALSE
      )
    }))
    rownames(edge_frame) <- NULL
    edge_frame <- edge_frame[order(
      edge_frame$layer, edge_frame$from_actor, edge_frame$to_actor
    ), , drop = FALSE]
  } else {
    edge_frame$source_records <- integer()
  }

  list(
    actors = actor_frame,
    layers = layer_frame,
    edges = edge_frame,
    source_path = normalizePath(path),
    source_sha256 = .real_sha256_file(path)
  )
}

.real_actor_order <- function(actor, salt) {
  if (!length(actor)) return(actor)
  key <- vapply(actor, function(id) {
    .graph_v2_hash(list(protocol = .real_alignment_protocol_version, salt, id))
  }, character(1))
  actor[order(key, actor)]
}

.real_actor_partition <- function(shared_actor, salt) {
  ordered <- .real_actor_order(shared_actor, paste0(salt, "-partition"))
  n <- length(ordered)
  if (n < 4L) {
    stop("A real network task needs at least four shared actors.", call. = FALSE)
  }
  n_train <- min(max(2L, ceiling(0.2 * n)), n - 2L)
  n_dev <- min(max(1L, floor(0.2 * n)), n - n_train - 1L)
  role <- rep("final", n)
  role[seq_len(n_train)] <- "train"
  if (n_dev) role[n_train + seq_len(n_dev)] <- "dev"
  names(role) <- ordered
  role
}

.real_actor_features <- function(actor_frame) {
  if (ncol(actor_frame) == 1L) {
    out <- matrix(0, nrow(actor_frame), 1L)
    colnames(out) <- "no_observed_attribute"
    rownames(out) <- actor_frame$actor
    return(out)
  }
  data <- actor_frame[, -1L, drop = FALSE]
  data[] <- lapply(data, function(x) {
    x <- as.character(x)
    x[is.na(x) | !nzchar(x)] <- "missing"
    factor(x, levels = sort(unique(x)))
  })
  formula <- stats::reformulate(names(data), intercept = FALSE)
  out <- stats::model.matrix(formula, data = data)
  storage.mode(out) <- "double"
  rownames(out) <- actor_frame$actor
  out
}

.real_network_fixture_file <- function(dataset_id) {
  name <- switch(
    dataset_id,
    aucs = "aucs.mpx",
    sampson_monastery = "monastery.mpx",
    stop("Unknown packaged network dataset: ", dataset_id, call. = FALSE)
  )
  .real_benchmark_inst_file(file.path("extdata", "alignment-real"), name)
}

real_alignment_network_case <- function(task_id) {
  registry <- real_alignment_benchmark_registry()
  task <- registry[registry$task_id == task_id, , drop = FALSE]
  if (nrow(task) != 1L ||
      !identical(task$task_type[[1L]], "attributed_network_correspondence")) {
    stop("Unknown attributed-network task: ", task_id, call. = FALSE)
  }
  source_path <- .real_network_fixture_file(task$dataset_id[[1L]])
  raw <- .real_mpx_parse(source_path)
  source_layer <- task$source_domain[[1L]]
  target_layer <- task$target_domain[[1L]]
  if (!all(c(source_layer, target_layer) %in% raw$layers$layer)) {
    stop("The governed network layers are missing from the source.",
      call. = FALSE
    )
  }

  layer_edges <- function(layer) {
    raw$edges[raw$edges$layer == layer, , drop = FALSE]
  }
  source_edges <- layer_edges(source_layer)
  target_edges <- layer_edges(target_layer)
  active_actor <- function(edges) {
    sort(unique(c(edges$from_actor, edges$to_actor)))
  }
  source_actor <- .real_actor_order(
    active_actor(source_edges), paste0(task_id, "-source")
  )
  target_actor <- .real_actor_order(
    active_actor(target_edges), paste0(task_id, "-target")
  )
  shared_actor <- intersect(source_actor, target_actor)
  partition <- .real_actor_partition(shared_actor, task$split_id[[1L]])
  all_features <- .real_actor_features(raw$actors)

  make_domain <- function(actor, edges, layer, domain_id) {
    index <- setNames(seq_along(actor), actor)
    relation <- data.frame(
      from = unname(index[edges$from_actor]),
      to = unname(index[edges$to_actor]),
      weight = edges$weight,
      relation = layer,
      directed = edges$directed,
      stringsAsFactors = FALSE
    )
    if (anyNA(relation$from) || anyNA(relation$to)) {
      stop("An observed edge endpoint is missing from its layer.",
        call. = FALSE
      )
    }
    role <- unname(partition[actor])
    role[is.na(role)] <- "final"
    attrs <- raw$actors[
      match(actor, raw$actors$actor), -1L, drop = FALSE
    ]
    rownames(attrs) <- NULL
    shared <- actor %in% shared_actor
    list(
      domain_id = domain_id,
      X = all_features[actor, , drop = FALSE],
      relations = relation,
      positions = matrix(numeric(), nrow = length(actor), ncol = 0L),
      attributes = attrs,
      labels = rep(NA_character_, length(actor)),
      row_id = actor,
      node_role = role,
      shared_mask = shared,
      private_mask = !shared
    )
  }

  domain1 <- make_domain(
    source_actor, source_edges, source_layer, paste0(task_id, "-source")
  )
  domain2 <- make_domain(
    target_actor, target_edges, target_layer, paste0(task_id, "-target")
  )
  map12 <- match(domain1$row_id, domain2$row_id)
  map12[domain1$private_mask] <- NA_integer_
  anchor_actor <- names(partition)[partition == "train"]
  anchors <- data.frame(
    source_row = match(anchor_actor, domain1$row_id),
    target_row = match(anchor_actor, domain2$row_id),
    actor_id = anchor_actor,
    stringsAsFactors = FALSE
  )
  anchors <- anchors[order(anchors$source_row), , drop = FALSE]
  rownames(anchors) <- NULL

  case <- list(
    protocol_version = .real_alignment_protocol_version,
    task_id = task_id,
    dataset_id = task$dataset_id[[1L]],
    task_type = task$task_type[[1L]],
    stress_regime = task$stress_regime[[1L]],
    split_id = task$split_id[[1L]],
    domains = list(domain1, domain2),
    relation_lists = list(domain1$relations, domain2$relations),
    anchors = anchors,
    truth = list(
      map12 = as.integer(map12),
      source_private_mask = domain1$private_mask,
      target_private_mask = domain2$private_mask,
      correspondence = data.frame(
        source_row = which(!is.na(map12)),
        target_row = map12[!is.na(map12)],
        actor_id = domain1$row_id[!is.na(map12)],
        stringsAsFactors = FALSE
      )
    ),
    governance = list(
      license = task$license[[1L]], source_url = task$source_url[[1L]],
      source_version = task$source_version[[1L]],
      source_sha256 = raw$source_sha256,
      identity_use = "audit_and_scoring_only",
      final_labels_available_to_fit = FALSE
    )
  )
  class(case) <- c("real_graph_alignment_case", "list")
  problems <- validate_real_alignment_case(case)
  if (nrow(problems)) {
    stop(paste(problems$message, collapse = "\n"), call. = FALSE)
  }
  case
}

validate_real_alignment_case <- function(case) {
  problems <- list()
  add <- function(code, message) {
    problems[[length(problems) + 1L]] <<- data.frame(
      code = code, message = message, stringsAsFactors = FALSE
    )
  }
  if (!inherits(case, "real_graph_alignment_case")) {
    add("class", "Object is not a real_graph_alignment_case.")
  }
  if (length(case$domains) != 2L) {
    add("domains", "A correspondence task must contain exactly two domains.")
  } else {
    for (i in seq_along(case$domains)) {
      domain <- case$domains[[i]]
      n <- nrow(domain$X)
      if (length(domain$row_id) != n || anyDuplicated(domain$row_id)) {
        add("row_identity", paste0("Domain ", i, " row IDs are invalid."))
      }
      if (length(domain$node_role) != n ||
          !all(domain$node_role %in% c("train", "dev", "final"))) {
        add("node_role", paste0("Domain ", i, " node roles are invalid."))
      }
      if (!identical(domain$shared_mask, !domain$private_mask)) {
        add("overlap_mask", paste0("Domain ", i, " overlap masks disagree."))
      }
      rel <- domain$relations
      if (nrow(rel) && (any(rel$from < 1L | rel$from > n) ||
          any(rel$to < 1L | rel$to > n) ||
          any(!is.finite(rel$weight)) || any(rel$weight < 0))) {
        add("relations", paste0("Domain ", i, " relations are invalid."))
      }
      if (any(grepl("(^|_)actor($|_)|(^|_)id($|_)", colnames(domain$X)))) {
        add("identity_leakage", "Actor identity appears in model features.")
      }
      if (any(!is.na(domain$labels[domain$node_role == "final"]))) {
        add("final_label_leakage", "Final labels are visible to fitting code.")
      }
    }
    map12 <- case$truth$map12
    if (length(map12) != nrow(case$domains[[1L]]$X)) {
      add("map_length", "Truth map length disagrees with the source domain.")
    } else {
      shared <- which(!is.na(map12))
      if (length(shared) && !identical(
        case$domains[[1L]]$row_id[shared],
        case$domains[[2L]]$row_id[map12[shared]]
      )) {
        add("map_identity", "Truth map does not preserve actor identity.")
      }
    }
    if (nrow(case$anchors)) {
      expected <- case$truth$map12[case$anchors$source_row]
      if (!identical(as.integer(case$anchors$target_row), as.integer(expected)) ||
          any(case$domains[[1L]]$node_role[case$anchors$source_row] != "train")) {
        add("anchors", "Anchors are not exact training correspondences.")
      }
    }
  }
  if (!length(problems)) {
    return(data.frame(
      code = character(), message = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, problems)
}

.real_pairwise_nn <- function(source, target) {
  source <- as.matrix(source)
  target <- as.matrix(target)
  if (ncol(source) != ncol(target)) {
    stop("Pairwise diagnostic features must have equal dimensions.",
      call. = FALSE
    )
  }
  both <- rbind(source, target)
  center <- colMeans(both)
  scale <- apply(both, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-12] <- 1
  source <- sweep(sweep(source, 2L, center), 2L, scale, "/")
  target <- sweep(sweep(target, 2L, center), 2L, scale, "/")
  d2 <- pmax(
    outer(rowSums(source^2), rowSums(target^2), "+") -
      2 * tcrossprod(source, target),
    0
  )
  max.col(-d2, ties.method = "first")
}

.real_topology_features <- function(domain) {
  n <- nrow(domain$X)
  rel <- domain$relations
  directed <- nrow(rel) && any(rel$directed)
  vertices <- data.frame(name = as.character(seq_len(n)))
  edge_frame <- data.frame(
    from = as.character(rel$from), to = as.character(rel$to),
    weight = rel$weight, stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(
    edge_frame, directed = directed, vertices = vertices
  )
  mode_out <- if (directed) "out" else "all"
  mode_in <- if (directed) "in" else "all"
  undirected <- igraph::as_undirected(graph, mode = "collapse")
  transitivity <- suppressWarnings(igraph::transitivity(
    undirected, type = "local", isolates = "zero"
  ))
  transitivity[!is.finite(transitivity)] <- 0
  out <- cbind(
    degree_out = igraph::degree(graph, mode = mode_out),
    degree_in = igraph::degree(graph, mode = mode_in),
    strength_out = igraph::strength(graph, mode = mode_out, weights = rel$weight),
    strength_in = igraph::strength(graph, mode = mode_in, weights = rel$weight),
    pagerank = igraph::page_rank(
      graph, directed = directed, weights = rel$weight
    )$vector,
    transitivity = transitivity
  )
  unname(out)
}

.real_network_control_row <- function(
  case, control_id, prediction = NULL, observed_topology_used = FALSE,
  status = "ok", status_code = "ok", status_detail = NA_character_
) {
  if (!identical(status, "ok")) {
    return(data.frame(
      protocol_version = case$protocol_version,
      task_id = case$task_id,
      dataset_id = case$dataset_id,
      control_id = control_id,
      status = status,
      status_code = status_code,
      status_detail = status_detail,
      observed_topology_used = observed_topology_used,
      n_eligible = sum(!is.na(case$truth$map12)),
      top1_accuracy = NA_real_, coverage = 0,
      anchor_top1_accuracy = NA_real_,
      nonanchor_top1_accuracy = NA_real_,
      private_rejection_accuracy = NA_real_,
      prediction_hash = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  normalized <- .graph_benchmark_normalize_prediction(case, prediction)
  if (!isTRUE(normalized$ok)) {
    return(.real_network_control_row(
      case, control_id, status = "algorithm_error",
      status_code = "malformed_output", status_detail = normalized$detail,
      observed_topology_used = observed_topology_used
    ))
  }
  data.frame(
    protocol_version = case$protocol_version,
    task_id = case$task_id,
    dataset_id = case$dataset_id,
    control_id = control_id,
    status = "ok", status_code = "ok", status_detail = NA_character_,
    observed_topology_used = observed_topology_used,
    n_eligible = normalized$n_eligible,
    top1_accuracy = normalized$top1_accuracy,
    coverage = normalized$coverage,
    anchor_top1_accuracy = normalized$anchor_top1_accuracy,
    nonanchor_top1_accuracy = normalized$nonanchor_top1_accuracy,
    private_rejection_accuracy = normalized$private_rejection_accuracy,
    prediction_hash = normalized$prediction_hash,
    stringsAsFactors = FALSE
  )
}

benchmark_real_network_controls <- function(case) {
  problems <- validate_real_alignment_case(case)
  if (nrow(problems)) {
    stop(paste(problems$message, collapse = "\n"), call. = FALSE)
  }
  d1 <- case$domains[[1L]]
  d2 <- case$domains[[2L]]
  attribute_prediction <- .real_pairwise_nn(d1$X, d2$X)
  topology_prediction <- .real_pairwise_nn(
    .real_topology_features(d1), .real_topology_features(d2)
  )
  row_prediction <- seq_len(nrow(d1$X))
  row_prediction[row_prediction > nrow(d2$X)] <- NA_integer_
  chance_prediction <- .graph_v2_with_seed(1701L, function() {
    sample.int(
      nrow(d2$X), nrow(d1$X), replace = nrow(d1$X) > nrow(d2$X)
    )
  })
  do.call(rbind, list(
    .real_network_control_row(case, "raw_attribute_knn", attribute_prediction),
    .real_network_control_row(
      case, "topology_descriptor_only", topology_prediction,
      observed_topology_used = TRUE
    ),
    .real_network_control_row(
      case, "position_only", status = "unsupported",
      status_code = "positions_not_available",
      status_detail = "This network dataset has no governed spatial positions."
    ),
    .real_network_control_row(case, "row_index", row_prediction),
    .real_network_control_row(case, "chance", chance_prediction),
    .real_network_control_row(
      case, "target_only_classifier", status = "unsupported",
      status_code = "not_a_classification_task",
      status_detail = "Target-only classification is not defined for correspondence."
    )
  ))
}

.real_nearest_centroid_fit <- function(x, y) {
  x <- as.matrix(x)
  y <- as.character(y)
  classes <- sort(unique(y[!is.na(y)]))
  centers <- do.call(rbind, lapply(classes, function(label) {
    colMeans(x[y == label, , drop = FALSE])
  }))
  rownames(centers) <- classes
  list(centers = centers)
}

.real_nearest_centroid_predict <- function(fit, x) {
  x <- as.matrix(x)
  d2 <- pmax(
    outer(rowSums(x^2), rowSums(fit$centers^2), "+") -
      2 * tcrossprod(x, fit$centers),
    0
  )
  rownames(fit$centers)[max.col(-d2, ties.method = "first")]
}

.real_classification_metrics <- function(truth, prediction) {
  truth <- as.character(truth)
  prediction <- as.character(prediction)
  valid <- !is.na(prediction)
  classes <- sort(unique(truth))
  recall <- vapply(classes, function(label) {
    index <- truth == label
    if (!any(index)) return(NA_real_)
    sum(valid[index] & prediction[index] == label) / sum(index)
  }, numeric(1))
  f1 <- vapply(classes, function(label) {
    tp <- sum(valid & truth == label & prediction == label)
    fp <- sum(valid & truth != label & prediction == label)
    fn <- sum(truth == label & (!valid | prediction != label))
    if (2 * tp + fp + fn == 0) return(NA_real_)
    2 * tp / (2 * tp + fp + fn)
  }, numeric(1))
  list(
    accuracy = sum(valid & prediction == truth) / length(truth),
    balanced_accuracy = mean(recall, na.rm = TRUE),
    macro_f1 = mean(f1, na.rm = TRUE),
    coverage = mean(valid)
  )
}

.real_so2sat_case_from_csv <- function(path, evidence_scope) {
  path <- normalizePath(path, mustWork = TRUE)
  data <- utils::read.csv(
    path, stringsAsFactors = FALSE, check.names = FALSE
  )
  required <- c(
    "split", "sample_id", "source_index", "city", "epsg", "x", "y",
    "label"
  )
  if (!all(required %in% names(data)) || anyDuplicated(data$sample_id)) {
    stop("The So2Sat feature artifact schema or identities are invalid.",
      call. = FALSE
    )
  }
  s1_columns <- grep("^s1_", names(data), value = TRUE)
  s2_columns <- grep("^s2_", names(data), value = TRUE)
  if (!length(s1_columns) || !identical(
    sub("^s1_", "", s1_columns), sub("^s2_", "", s2_columns)
  )) {
    stop("The So2Sat sensor feature contracts disagree.", call. = FALSE)
  }
  development <- data$split == "validation_west"
  final <- data$split == "testing_east"
  if (!any(development) || !any(final)) {
    stop("The So2Sat artifact lacks its frozen blocked split.", call. = FALSE)
  }
  feature_data <- data[, c(s1_columns, s2_columns), drop = FALSE]
  if (any(!is.finite(as.matrix(feature_data))) ||
      any(!is.finite(data$x)) || any(!is.finite(data$y))) {
    stop("The So2Sat feature artifact contains non-finite model inputs.",
      call. = FALSE
    )
  }
  case <- list(
    protocol_version = .real_alignment_protocol_version,
    task_id = "so2sat_lcz42_sar_optical_culture10",
    dataset_id = "so2sat_lcz42",
    fit = list(
      source_X = as.matrix(data[development, s1_columns, drop = FALSE]),
      target_X = as.matrix(data[development, s2_columns, drop = FALSE]),
      source_labels = data$label[development],
      target_labels = data$label[development],
      source_positions = as.matrix(data[development, c("x", "y")]),
      target_positions = as.matrix(data[development, c("x", "y")]),
      source_city = data$city[development],
      target_city = data$city[development],
      source_row_id = data$sample_id[development],
      target_row_id = data$sample_id[development],
      final_target_X = as.matrix(data[final, s2_columns, drop = FALSE]),
      final_target_positions = as.matrix(data[final, c("x", "y")]),
      final_target_city = data$city[final],
      final_target_row_id = data$sample_id[final]
    ),
    evaluation = list(
      final_target_labels = data$label[final],
      final_target_row_id = data$sample_id[final]
    ),
    governance = list(
      license = "CC BY 4.0",
      source_revision = "b5c817486899935e864832b93086cc87f3eee473",
      split_id = "culture10-v4.2-blocked",
      evidence_scope = evidence_scope,
      final_labels_available_to_fit = FALSE,
      feature_artifact = path,
      feature_artifact_sha256 = .real_sha256_file(path),
      n_development = sum(development),
      n_final = sum(final)
    )
  )
  if ("final_target_labels" %in% names(case$fit) ||
      !identical(case$fit$final_target_row_id,
        case$evaluation$final_target_row_id
      )) {
    stop("The So2Sat fitting view leaks or loses final identities.",
      call. = FALSE
    )
  }
  class(case) <- c("real_remote_sensing_case", "list")
  case
}

real_so2sat_ci_case <- function() {
  path <- .real_benchmark_inst_file(
    file.path("extdata", "alignment-real"),
    "so2sat_lcz42_ci_features.csv"
  )
  .real_so2sat_case_from_csv(path, evidence_scope = "ci_schema_only")
}

real_so2sat_full_case <- function(path) {
  .real_so2sat_case_from_csv(path, evidence_scope = "full_frozen_split")
}

.real_position_predict <- function(
  train_position, train_city, train_label, test_position, test_city
) {
  prediction <- rep(NA_character_, nrow(test_position))
  for (city in unique(test_city)) {
    train <- which(train_city == city)
    test <- which(test_city == city)
    if (!length(train)) next
    d2 <- pmax(
      outer(rowSums(test_position[test, , drop = FALSE]^2),
        rowSums(train_position[train, , drop = FALSE]^2), "+"
      ) - 2 * tcrossprod(
        test_position[test, , drop = FALSE],
        train_position[train, , drop = FALSE]
      ),
      0
    )
    prediction[test] <- train_label[train[max.col(-d2, ties.method = "first")]]
  }
  prediction
}

.real_class_control_row <- function(case, control_id, prediction = NULL,
                                    status = "ok", status_code = "ok",
                                    status_detail = NA_character_) {
  if (!identical(status, "ok")) {
    return(data.frame(
      protocol_version = case$protocol_version,
      task_id = case$task_id,
      dataset_id = case$dataset_id,
      control_id = control_id,
      status = status, status_code = status_code,
      status_detail = status_detail,
      n_eligible = length(case$evaluation$final_target_labels),
      accuracy = NA_real_, balanced_accuracy = NA_real_,
      macro_f1 = NA_real_, coverage = 0,
      prediction_hash = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  metrics <- .real_classification_metrics(
    case$evaluation$final_target_labels, prediction
  )
  data.frame(
    protocol_version = case$protocol_version,
    task_id = case$task_id,
    dataset_id = case$dataset_id,
    control_id = control_id,
    status = "ok", status_code = "ok", status_detail = NA_character_,
    n_eligible = length(case$evaluation$final_target_labels),
    accuracy = metrics$accuracy,
    balanced_accuracy = metrics$balanced_accuracy,
    macro_f1 = metrics$macro_f1,
    coverage = metrics$coverage,
    prediction_hash = .graph_v2_hash(list(
      row_id = case$evaluation$final_target_row_id,
      prediction = prediction
    )),
    stringsAsFactors = FALSE
  )
}

benchmark_real_so2sat_controls <- function(case = real_so2sat_ci_case()) {
  source_fit <- .real_nearest_centroid_fit(
    case$fit$source_X, case$fit$source_labels
  )
  target_fit <- .real_nearest_centroid_fit(
    case$fit$target_X, case$fit$target_labels
  )
  source_prediction <- .real_nearest_centroid_predict(
    source_fit, case$fit$final_target_X
  )
  target_prediction <- .real_nearest_centroid_predict(
    target_fit, case$fit$final_target_X
  )
  position_prediction <- .real_position_predict(
    case$fit$target_positions, case$fit$target_city,
    case$fit$target_labels, case$fit$final_target_positions,
    case$fit$final_target_city
  )
  chance_prediction <- .graph_v2_with_seed(2701L, function() {
    sample(
      sort(unique(case$fit$source_labels)),
      length(case$evaluation$final_target_labels), replace = TRUE
    )
  })
  do.call(rbind, list(
    .real_class_control_row(case, "raw_attribute_source_only", source_prediction),
    .real_class_control_row(case, "position_only", position_prediction),
    .real_class_control_row(
      case, "topology_descriptor_only", status = "unsupported",
      status_code = "observed_graph_not_available",
      status_detail = "So2Sat patch classification has no observed node graph."
    ),
    .real_class_control_row(case, "target_only_classifier", target_prediction),
    .real_class_control_row(case, "chance", chance_prediction)
  ))
}

real_alignment_status_row <- function(
  task_id, stage, status, status_code, status_detail = NA_character_,
  artifact = NA_character_
) {
  allowed <- c(
    "ok", "dataset_unavailable", "dataset_invalid", "algorithm_error"
  )
  status <- match.arg(status, allowed)
  data.frame(
    protocol_version = .real_alignment_protocol_version,
    task_id = task_id, stage = stage, artifact = artifact,
    status = status, status_code = status_code,
    status_detail = status_detail,
    failure_class = if (startsWith(status, "dataset_")) {
      "dataset"
    } else if (identical(status, "algorithm_error")) {
      "algorithm"
    } else {
      NA_character_
    },
    stringsAsFactors = FALSE
  )
}
