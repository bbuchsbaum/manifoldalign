.graph_benchmark_base_coords <- function(n, d, structure) {
  if (identical(structure, "ring")) {
    angles <- 2 * pi * seq_len(n) / n
    X <- cbind(cos(angles), sin(angles))
    if (d > 2L) {
      X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
    } else if (d == 1L) {
      X <- matrix(cos(angles), n, 1L)
    }
    return(X)
  }

  if (identical(structure, "grid")) {
    grid_size <- ceiling(sqrt(n))
    grid_coords <- expand.grid(
      x = seq_len(grid_size) / grid_size,
      y = seq_len(grid_size) / grid_size
    )[seq_len(n), ]
    X <- cbind(grid_coords$x, grid_coords$y)
    if (d > 2L) {
      X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
    } else if (d == 1L) {
      X <- matrix(grid_coords$x, n, 1L)
    }
    return(as.matrix(X))
  }

  if (identical(structure, "community")) {
    k <- 4L
    X <- matrix(0, n, max(2L, d))
    n_per <- ceiling(n / k)
    for (g in seq_len(k)) {
      idx <- seq.int((g - 1L) * n_per + 1L, min(g * n_per, n))
      center <- c(cos(2 * pi * g / k), sin(2 * pi * g / k))
      X[idx, 1:2] <- matrix(stats::rnorm(length(idx) * 2, sd = 0.2), length(idx), 2) +
        matrix(center, length(idx), 2, byrow = TRUE)
    }
    if (d > 2L) {
      X[, 3:d] <- matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L)
    }
    return(X[, seq_len(d), drop = FALSE])
  }

  matrix(stats::rnorm(n * d), n, d)
}

synthetic_graph_alignment_registry <- function(sizes = c(50L, 100L),
                                               structures = c("ring", "grid", "random", "community"),
                                               d = 3L,
                                               noise_sd = 0.05,
                                               permute_fraction = 1,
                                               n_anchors = 0L,
                                               n_reps = 3L,
                                               seed = 1L) {
  structures <- match.arg(
    structures,
    c("ring", "grid", "random", "community", "sbm", "geometric", "scale_free"),
    several.ok = TRUE
  )
  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) {
    stop("`sizes` must contain integers >= 3.", call. = FALSE)
  }

  rows <- list()
  ptr <- 1L
  for (structure in structures) {
    for (s in seq_along(sizes)) {
      n <- sizes[[s]]
      for (r in seq_len(as.integer(n_reps))) {
        generation_seed <- as.integer(seed) + match(structure, structures) * 10007L + s * 1009L + r * 17L
        rows[[ptr]] <- data.frame(
          scenario_family = "graph",
          scenario = paste(structure, n, r, sep = "_"),
          structure = structure,
          n = n,
          d = as.integer(d),
          noise_sd = as.numeric(noise_sd),
          permute_fraction = as.numeric(permute_fraction),
          n_anchors = as.integer(n_anchors),
          rep = as.integer(r),
          generation_seed = generation_seed,
          stringsAsFactors = FALSE
        )
        ptr <- ptr + 1L
      }
    }
  }

  do.call(rbind, rows)
}

synthetic_graph_alignment_case <- function(n,
                                           d = 3L,
                                           structure = c("ring", "grid", "random", "community"),
                                           noise_sd = 0.05,
                                           permute_fraction = 1,
                                           n_anchors = 0L,
                                           seed = 1L) {
  structure <- match.arg(structure)
  .local_seed(as.integer(seed))

  X1 <- .graph_benchmark_base_coords(n = as.integer(n), d = as.integer(d), structure = structure)

  n_permute <- floor(as.integer(n) * as.numeric(permute_fraction))
  perm_idx <- if (n_permute >= 2L) sort(sample.int(as.integer(n), n_permute)) else integer(0)
  perm21 <- seq_len(as.integer(n))
  if (length(perm_idx)) {
    perm21[perm_idx] <- perm_idx[sample.int(length(perm_idx))]
  }

  X2 <- X1[perm21, , drop = FALSE] +
    matrix(stats::rnorm(as.integer(n) * as.integer(d), sd = noise_sd), as.integer(n), as.integer(d))
  map12 <- match(seq_len(as.integer(n)), perm21)

  a1 <- rep(NA_integer_, as.integer(n))
  a2 <- rep(NA_integer_, as.integer(n))
  anchor_idx <- integer(0)
  if (as.integer(n_anchors) > 0L) {
    available <- setdiff(seq_len(as.integer(n)), perm_idx)
    if (length(available) < as.integer(n_anchors)) {
      available <- seq_len(as.integer(n))
    }
    anchor_idx <- sort(sample(available, min(as.integer(n_anchors), length(available))))
    a1[anchor_idx] <- seq_along(anchor_idx)
    a2[map12[anchor_idx]] <- seq_along(anchor_idx)
  }

  hd <- as_hyperdesign(
    X_list = list(domain1 = X1, domain2 = X2),
    labels = list(a1, a2),
    label_name = "anchors"
  )

  correspondences <- NULL
  if (length(anchor_idx)) {
    correspondences <- data.frame(
      domain_i = rep(1L, length(anchor_idx)),
      index_i = anchor_idx,
      domain_j = rep(2L, length(anchor_idx)),
      index_j = map12[anchor_idx],
      stringsAsFactors = FALSE
    )
  }

  list(
    X1 = X1,
    X2 = X2,
    perm21 = perm21,
    map12 = map12,
    permuted_index = perm_idx,
    anchor_idx = anchor_idx,
    anchors1 = a1,
    anchors2 = a2,
    correspondences = correspondences,
    hd = hd,
    anchors_name = "anchors",
    benchmark_meta = list(
      scenario_family = "graph",
      structure = structure,
      n = as.integer(n),
      d = as.integer(d),
      noise_sd = as.numeric(noise_sd),
      permute_fraction = as.numeric(permute_fraction),
      n_anchors = as.integer(n_anchors),
      seed = as.integer(seed)
    )
  )
}

# Frozen graph benchmark protocol v2 -------------------------------------

.graph_v2_protocol_version <- "2.0.0"
.graph_v2_preprocessing_id <- "train_numeric_center_scale_v1"

.graph_v2_with_seed <- function(seed, code) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  code()
}

.graph_v2_hash <- function(x) {
  path <- tempfile("manifoldalign-graph-v2-hash-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, version = 3, compress = FALSE)
  unname(tools::md5sum(path))
}

.graph_v2_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < lower || x > upper) {
    stop(
      "`", name, "` must be one finite number in [", lower, ", ", upper,
      "].", call. = FALSE
    )
  }
  as.numeric(x)
}

.graph_v2_scalar_integer <- function(x, name, lower = 0L) {
  value <- .graph_v2_scalar_number(x, name, lower = lower)
  if (value != floor(value)) {
    stop("`", name, "` must be integer-valued.", call. = FALSE)
  }
  as.integer(value)
}

.graph_v2_config <- function(
  n_shared,
  n_private,
  family,
  attr_dim,
  n_anchors,
  label_fraction,
  edge_delete,
  edge_add,
  edge_rewire,
  weight_noise,
  weight_scale,
  directed,
  heterophily,
  misleading_fraction,
  position_jitter,
  rigid_transform,
  missing_attribute_fraction,
  attribute_noise,
  class_imbalance,
  private_classes,
  split,
  tier,
  preprocessing_id,
  tuning_budget,
  seed,
  case_id
) {
  family <- match.arg(family, c("sbm", "geometric", "scale_free"))
  split <- match.arg(split, c("development", "final"))
  tier <- match.arg(
    tier,
    c("oracle", "small_accuracy", "medium_robustness", "large_scaling")
  )
  n_shared <- .graph_v2_scalar_integer(n_shared, "n_shared", 4L)
  if (!is.numeric(n_private) || length(n_private) != 2L || anyNA(n_private) ||
      any(!is.finite(n_private)) || any(n_private < 0) ||
      any(n_private != floor(n_private))) {
    stop("`n_private` must contain two non-negative integers.", call. = FALSE)
  }
  n_private <- as.integer(n_private)
  attr_dim <- .graph_v2_scalar_integer(attr_dim, "attr_dim", 1L)
  n_anchors <- .graph_v2_scalar_integer(n_anchors, "n_anchors", 0L)
  if (n_anchors > n_shared) {
    stop("`n_anchors` cannot exceed `n_shared`.", call. = FALSE)
  }
  tuning_budget <- .graph_v2_scalar_integer(
    tuning_budget, "tuning_budget", 0L
  )
  if (identical(split, "final") && tuning_budget != 0L) {
    stop("Final benchmark cases must have `tuning_budget = 0`.", call. = FALSE)
  }
  fraction <- function(x, name) {
    .graph_v2_scalar_number(x, name, lower = 0, upper = 1)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed) ||
      !is.logical(rigid_transform) || length(rigid_transform) != 1L ||
      is.na(rigid_transform) || !is.logical(private_classes) ||
      length(private_classes) != 1L || is.na(private_classes)) {
    stop(
      "`directed`, `rigid_transform`, and `private_classes` must be scalar logicals.",
      call. = FALSE
    )
  }
  if (!is.character(preprocessing_id) || length(preprocessing_id) != 1L ||
      is.na(preprocessing_id) || !nzchar(preprocessing_id)) {
    stop("`preprocessing_id` must be one non-empty string.", call. = FALSE)
  }
  if (!identical(preprocessing_id, .graph_v2_preprocessing_id)) {
    stop("Unknown frozen preprocessing ID: ", preprocessing_id, call. = FALSE)
  }
  if (!is.character(case_id) || length(case_id) != 1L || is.na(case_id) ||
      !nzchar(case_id)) {
    stop("`case_id` must be one non-empty string.", call. = FALSE)
  }

  list(
    protocol_version = .graph_v2_protocol_version,
    case_id = case_id,
    split = split,
    tier = tier,
    family = family,
    n_shared = n_shared,
    n_private = n_private,
    attr_dim = attr_dim,
    n_anchors = n_anchors,
    label_fraction = fraction(label_fraction, "label_fraction"),
    edge_delete = fraction(edge_delete, "edge_delete"),
    edge_add = fraction(edge_add, "edge_add"),
    edge_rewire = fraction(edge_rewire, "edge_rewire"),
    weight_noise = .graph_v2_scalar_number(
      weight_noise, "weight_noise", lower = 0
    ),
    weight_scale = .graph_v2_scalar_number(
      weight_scale, "weight_scale", lower = .Machine$double.eps
    ),
    directed = isTRUE(directed),
    heterophily = fraction(heterophily, "heterophily"),
    misleading_fraction = fraction(
      misleading_fraction, "misleading_fraction"
    ),
    position_jitter = .graph_v2_scalar_number(
      position_jitter, "position_jitter", lower = 0
    ),
    rigid_transform = isTRUE(rigid_transform),
    missing_attribute_fraction = fraction(
      missing_attribute_fraction, "missing_attribute_fraction"
    ),
    attribute_noise = .graph_v2_scalar_number(
      attribute_noise, "attribute_noise", lower = 0
    ),
    class_imbalance = fraction(class_imbalance, "class_imbalance"),
    private_classes = isTRUE(private_classes),
    preprocessing_id = preprocessing_id,
    tuning_budget = tuning_budget,
    seed = .graph_v2_scalar_integer(seed, "seed", 1L)
  )
}

.graph_v2_class_labels <- function(n, imbalance, prefix = "class") {
  k <- min(4L, n)
  major <- 1 / k + imbalance * (1 - 1 / k)
  probabilities <- c(major, rep((1 - major) / (k - 1L), k - 1L))
  counts <- floor(probabilities * n)
  counts[counts < 1L] <- 1L
  while (sum(counts) > n) {
    index <- which.max(counts[-1L]) + 1L
    counts[[index]] <- counts[[index]] - 1L
  }
  while (sum(counts) < n) {
    counts[[which.max(probabilities - counts / n)]] <-
      counts[[which.max(probabilities - counts / n)]] + 1L
  }
  sample(rep(paste0(prefix, seq_len(k)), counts), n)
}

.graph_v2_positions <- function(n, family, labels) {
  if (identical(family, "geometric")) {
    width <- ceiling(sqrt(n))
    grid <- as.matrix(expand.grid(
      x = seq_len(width) / width,
      y = seq_len(width) / width
    )[seq_len(n), , drop = FALSE])
    return(grid + matrix(stats::rnorm(n * 2L, sd = 0.015), n, 2L))
  }
  levels <- unique(labels)
  centers <- cbind(
    cos(2 * pi * seq_along(levels) / length(levels)),
    sin(2 * pi * seq_along(levels) / length(levels))
  )
  center <- centers[match(labels, levels), , drop = FALSE]
  spread <- if (identical(family, "sbm")) 0.3 else 0.75
  center + matrix(stats::rnorm(n * 2L, sd = spread), n, 2L)
}

.graph_v2_latent_features <- function(labels, positions, attr_dim) {
  levels <- unique(labels)
  centers <- matrix(
    stats::rnorm(length(levels) * attr_dim, sd = 1.5),
    length(levels), attr_dim
  )
  position_loading <- matrix(stats::rnorm(2L * attr_dim), 2L, attr_dim)
  centers[match(labels, levels), , drop = FALSE] +
    positions %*% position_loading +
    matrix(stats::rnorm(length(labels) * attr_dim, sd = 0.15),
           length(labels), attr_dim)
}

.graph_v2_orthogonal <- function(d) {
  decomposition <- qr.Q(qr(matrix(stats::rnorm(d * d), d, d)))
  if (det(decomposition) < 0) decomposition[, d] <- -decomposition[, d]
  decomposition
}

.graph_v2_canonical_edges <- function(edges, directed) {
  if (!nrow(edges)) return(edges)
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  if (!nrow(edges)) return(edges)
  if (!directed) {
    from <- pmin(edges$from, edges$to)
    to <- pmax(edges$from, edges$to)
    edges$from <- from
    edges$to <- to
  }
  key <- paste(edges$relation, edges$from, edges$to, sep = ":")
  ordering <- order(key, method = "radix")
  edges <- edges[ordering, , drop = FALSE]
  key <- key[ordering]
  first <- !duplicated(key)
  group <- cumsum(first)
  weight_sum <- rowsum(edges$weight, group, reorder = FALSE)[, 1L]
  count <- tabulate(group)
  out <- edges[first, , drop = FALSE]
  out$weight <- as.numeric(weight_sum / count)
  rownames(out) <- NULL
  out[order(out$relation, out$from, out$to), , drop = FALSE]
}

.graph_v2_base_edges <- function(positions, labels, family, directed,
                                 mean_degree = 6L) {
  n <- nrow(positions)
  mean_degree <- min(as.integer(mean_degree), n - 1L)
  if (identical(family, "geometric")) {
    neighbor <- RANN::nn2(
      positions, positions, k = min(n, mean_degree + 1L)
    )$nn.idx
    from <- rep(seq_len(n), each = ncol(neighbor) - 1L)
    to <- as.integer(c(t(neighbor[, -1L, drop = FALSE])))
  } else if (identical(family, "scale_free")) {
    graph <- igraph::sample_pa(
      n = n,
      power = 1,
      m = max(1L, floor(mean_degree / 2L)),
      directed = directed
    )
    edge_matrix <- igraph::as_edgelist(graph, names = FALSE)
    from <- as.integer(edge_matrix[, 1L])
    to <- as.integer(edge_matrix[, 2L])
  } else {
    members <- split(seq_len(n), labels)
    levels <- names(members)
    from <- rep(seq_len(n), each = mean_degree)
    to <- integer(length(from))
    cursor <- 1L
    for (i in seq_len(n)) {
      own <- members[[labels[[i]]]]
      same_count <- min(length(own) - 1L, ceiling(0.75 * mean_degree))
      same <- integer(0)
      if (same_count > 0L) {
        same <- sample(own, min(length(own), same_count + 1L))
        same <- head(same[same != i], same_count)
        if (length(same) < same_count) {
          same <- c(
            same,
            sample(
              own[!own %in% c(i, same)],
              same_count - length(same)
            )
          )
        }
      }
      other_count <- mean_degree - length(same)
      other <- integer(other_count)
      if (other_count) {
        other_levels <- setdiff(levels, labels[[i]])
        for (j in seq_len(other_count)) {
          chosen <- sample(other_levels, 1L)
          other[[j]] <- sample(members[[chosen]], 1L)
        }
      }
      chosen <- c(same, other)
      if (length(chosen) < mean_degree) {
        chosen <- c(
          chosen,
          sample(setdiff(seq_len(n), c(i, chosen)), mean_degree - length(chosen))
        )
      }
      to[cursor:(cursor + mean_degree - 1L)] <- chosen
      cursor <- cursor + mean_degree
    }
  }
  distance <- sqrt(rowSums(
    (positions[from, , drop = FALSE] - positions[to, , drop = FALSE])^2
  ))
  scale <- stats::median(distance[distance > 0])
  if (!is.finite(scale) || scale <= 0) scale <- 1
  edges <- data.frame(
    from = as.integer(from),
    to = as.integer(to),
    weight = pmax(exp(-distance / scale), 1e-6),
    relation = "structural",
    directed = isTRUE(directed),
    stringsAsFactors = FALSE
  )
  .graph_v2_canonical_edges(edges, directed)
}

.graph_v2_private_edges <- function(n_shared, n_private, directed,
                                    mean_degree = 4L) {
  n <- n_shared + n_private
  if (!n_private) {
    return(data.frame(
      from = integer(), to = integer(), weight = numeric(),
      relation = character(), directed = logical(),
      stringsAsFactors = FALSE
    ))
  }
  private <- n_shared + seq_len(n_private)
  degree <- min(mean_degree, n - 1L)
  from <- rep(private, each = degree)
  to <- unlist(lapply(private, function(i) {
    sample(setdiff(seq_len(n), i), degree)
  }), use.names = FALSE)
  edges <- data.frame(
    from = as.integer(from), to = as.integer(to),
    weight = stats::runif(length(from), 0.25, 1),
    relation = "structural", directed = isTRUE(directed),
    stringsAsFactors = FALSE
  )
  .graph_v2_canonical_edges(edges, directed)
}

.graph_v2_misleading_edges <- function(labels, fraction, directed, reference_n) {
  count <- floor(reference_n * fraction)
  if (!count) {
    return(data.frame(
      from = integer(), to = integer(), weight = numeric(),
      relation = character(), directed = logical(),
      stringsAsFactors = FALSE
    ))
  }
  n <- length(labels)
  from <- sample.int(n, count, replace = count > n)
  to <- integer(count)
  for (index in seq_len(count)) {
    candidates <- which(labels != labels[[from[[index]]]])
    if (!length(candidates)) candidates <- setdiff(seq_len(n), from[[index]])
    to[[index]] <- sample(candidates, 1L)
  }
  .graph_v2_canonical_edges(
    data.frame(
      from = from, to = to,
      weight = stats::runif(count, 0.75, 1.25),
      relation = "misleading", directed = isTRUE(directed),
      stringsAsFactors = FALSE
    ),
    directed
  )
}

.graph_v2_random_edges <- function(n, count, labels, heterophily, directed) {
  if (!count) {
    return(data.frame(
      from = integer(), to = integer(), weight = numeric(),
      relation = character(), directed = logical(),
      stringsAsFactors = FALSE
    ))
  }
  from <- sample.int(n, count, replace = count > n)
  to <- integer(count)
  for (index in seq_len(count)) {
    prefer_other <- stats::runif(1L) < heterophily
    candidates <- if (prefer_other) {
      which(labels != labels[[from[[index]]]])
    } else {
      setdiff(seq_len(n), from[[index]])
    }
    if (!length(candidates)) candidates <- setdiff(seq_len(n), from[[index]])
    to[[index]] <- sample(candidates, 1L)
  }
  data.frame(
    from = from, to = to,
    weight = stats::runif(count, 0.2, 1),
    relation = "structural", directed = isTRUE(directed),
    stringsAsFactors = FALSE
  )
}

.graph_v2_perturb_edges <- function(edges, labels, config) {
  structural <- which(edges$relation == "structural")
  base_n <- length(structural)
  delete_n <- min(base_n, floor(base_n * config$edge_delete))
  deleted <- if (delete_n) sample(structural, delete_n) else integer(0)
  kept <- edges[setdiff(seq_len(nrow(edges)), deleted), , drop = FALSE]

  structural <- which(kept$relation == "structural")
  rewire_n <- min(length(structural), floor(base_n * config$edge_rewire))
  rewired <- if (rewire_n) sample(structural, rewire_n) else integer(0)
  if (length(rewired)) {
    replacement <- .graph_v2_random_edges(
      length(labels), length(rewired), labels,
      heterophily = config$heterophily,
      directed = config$directed
    )
    kept$to[rewired] <- replacement$to
  }

  structural <- which(kept$relation == "structural")
  heterophily_n <- min(
    length(structural), floor(base_n * config$heterophily)
  )
  heterophilic <- if (heterophily_n) {
    sample(structural, heterophily_n)
  } else {
    integer(0)
  }
  for (index in heterophilic) {
    candidates <- which(labels != labels[[kept$from[[index]]]])
    if (length(candidates)) kept$to[[index]] <- sample(candidates, 1L)
  }

  add_n <- floor(base_n * config$edge_add)
  additions <- .graph_v2_random_edges(
    length(labels), add_n, labels,
    heterophily = config$heterophily,
    directed = config$directed
  )
  misleading <- .graph_v2_misleading_edges(
    labels, config$misleading_fraction, config$directed, base_n
  )
  out <- rbind(kept, additions, misleading)
  if (nrow(out)) {
    out$weight <- pmax(
      out$weight * config$weight_scale *
        exp(stats::rnorm(nrow(out), sd = config$weight_noise)),
      1e-8
    )
  }
  out <- .graph_v2_canonical_edges(out, config$directed)
  list(
    edges = out,
    counts = list(
      base = as.integer(base_n),
      deleted = as.integer(delete_n),
      rewired = as.integer(rewire_n),
      heterophilic = as.integer(heterophily_n),
      added_requested = as.integer(add_n),
      misleading = as.integer(sum(out$relation == "misleading"))
    )
  )
}

.graph_v2_mixed_attributes <- function(labels, positions, noise = 0) {
  n <- length(labels)
  data.frame(
    terrain_score = positions[, 1L] + stats::rnorm(n, sd = noise),
    elevation = positions[, 2L] + stats::rnorm(n, sd = noise),
    cover = as.character(labels),
    zone = paste0("z", 1L + (seq_len(n) %% 3L)),
    active = (seq_len(n) %% 2L) == 0L,
    stringsAsFactors = FALSE
  )
}

.graph_v2_corrupt_categories <- function(attributes, fraction) {
  if (fraction <= 0 || !nrow(attributes)) return(attributes)
  categorical <- which(vapply(attributes, is.character, logical(1)))
  for (column in categorical) {
    count <- floor(nrow(attributes) * min(fraction, 1))
    if (!count) next
    index <- sample.int(nrow(attributes), count)
    values <- unique(attributes[[column]])
    if (length(values) > 1L) {
      attributes[[column]][index] <- sample(values, length(index), replace = TRUE)
    }
  }
  attributes
}

.graph_v2_apply_missing <- function(attributes, fraction) {
  total <- nrow(attributes) * ncol(attributes)
  count <- floor(total * fraction)
  if (!count) return(attributes)
  cells <- sample.int(total, count)
  row <- ((cells - 1L) %% nrow(attributes)) + 1L
  column <- ((cells - 1L) %/% nrow(attributes)) + 1L
  for (index in seq_along(cells)) {
    attributes[[column[[index]]]][row[[index]]] <- NA
  }
  attributes
}

.graph_v2_roles <- function(n_shared, n_anchors) {
  order <- sample.int(n_shared)
  n_train <- max(n_anchors, ceiling(0.2 * n_shared))
  n_dev <- max(1L, ceiling(0.2 * n_shared))
  n_dev <- min(n_dev, n_shared - n_train)
  role <- rep("final", n_shared)
  if (n_train) role[order[seq_len(n_train)]] <- "train"
  if (n_dev) {
    index <- order[n_train + seq_len(n_dev)]
    role[index] <- "dev"
  }
  list(role = role, train_order = order[seq_len(n_train)])
}

.graph_v2_reindex_edges <- function(edges, latent_to_row) {
  out <- edges
  out$from <- latent_to_row[out$from]
  out$to <- latent_to_row[out$to]
  out[order(out$relation, out$from, out$to), , drop = FALSE]
}

.graph_v2_relation_list <- function(edges) {
  split_edges <- split(edges, edges$relation)
  lapply(split_edges, function(x) {
    x[, c("from", "to", "weight"), drop = FALSE]
  })
}

.graph_v2_edge_key <- function(edges, directed) {
  if (!nrow(edges)) return(character())
  if (!directed) {
    from <- pmin(edges$from, edges$to)
    to <- pmax(edges$from, edges$to)
  } else {
    from <- edges$from
    to <- edges$to
  }
  paste(edges$relation, from, to, sep = ":")
}

.graph_v2_case_hash_payload <- function(case) {
  list(
    protocol_version = case$protocol_version,
    case_id = case$case_id,
    split = case$split,
    tier = case$tier,
    config = case$config,
    domains = lapply(case$domains, function(domain) {
      domain[c(
        "X", "relations", "positions", "attributes", "labels",
        "row_id", "node_role", "shared_mask", "private_mask"
      )]
    }),
    truth = case$truth,
    anchors = case$anchors,
    perturbation = case$perturbation,
    preprocessing = case$preprocessing
  )
}

#' Generate a frozen-protocol synthetic attributed graph alignment case
#'
#' @keywords internal
#' @noRd
synthetic_graph_alignment_case_v2 <- function(
  n_shared,
  n_private = c(0L, 0L),
  family = c("sbm", "geometric", "scale_free"),
  attr_dim = 6L,
  n_anchors = 0L,
  label_fraction = 0.1,
  edge_delete = 0.1,
  edge_add = 0.05,
  edge_rewire = 0.05,
  weight_noise = 0.1,
  weight_scale = 1,
  directed = FALSE,
  heterophily = 0,
  misleading_fraction = 0,
  position_jitter = 0.05,
  rigid_transform = TRUE,
  missing_attribute_fraction = 0.05,
  attribute_noise = 0.05,
  class_imbalance = 0,
  private_classes = TRUE,
  split = c("development", "final"),
  tier = c("oracle", "small_accuracy", "medium_robustness", "large_scaling"),
  preprocessing_id = .graph_v2_preprocessing_id,
  tuning_budget = if (identical(split[[1L]], "final")) 0L else 12L,
  seed = 1L,
  case_id = paste0("interactive-", seed)
) {
  config <- .graph_v2_config(
    n_shared, n_private, match.arg(family), attr_dim, n_anchors,
    label_fraction, edge_delete, edge_add, edge_rewire, weight_noise,
    weight_scale, directed, heterophily, misleading_fraction,
    position_jitter, rigid_transform, missing_attribute_fraction,
    attribute_noise, class_imbalance, private_classes, match.arg(split),
    match.arg(tier), preprocessing_id, tuning_budget, seed, case_id
  )

  .graph_v2_with_seed(config$seed, function() {
    shared_labels <- .graph_v2_class_labels(
      config$n_shared, config$class_imbalance
    )
    shared_positions <- .graph_v2_positions(
      config$n_shared, config$family, shared_labels
    )
    shared_latent <- .graph_v2_latent_features(
      shared_labels, shared_positions, config$attr_dim
    )
    private_labels <- lapply(seq_len(2L), function(domain) {
      n <- config$n_private[[domain]]
      if (!n) return(character())
      if (config$private_classes) {
        rep(paste0("private_d", domain), n)
      } else {
        sample(unique(shared_labels), n, replace = TRUE)
      }
    })
    labels_latent <- list(
      c(shared_labels, private_labels[[1L]]),
      c(shared_labels, private_labels[[2L]])
    )
    private_positions <- lapply(config$n_private, function(n) {
      if (!n) return(matrix(numeric(), 0L, 2L))
      matrix(stats::runif(n * 2L, -1.5, 1.5), n, 2L)
    })
    positions_latent <- list(
      rbind(shared_positions, private_positions[[1L]]),
      rbind(shared_positions, private_positions[[2L]])
    )
    private_latent <- lapply(seq_len(2L), function(domain) {
      n <- config$n_private[[domain]]
      if (!n) return(matrix(numeric(), 0L, config$attr_dim))
      .graph_v2_latent_features(
        labels_latent[[domain]][config$n_shared + seq_len(n)],
        private_positions[[domain]], config$attr_dim
      )
    })
    latent <- list(
      rbind(shared_latent, private_latent[[1L]]),
      rbind(shared_latent, private_latent[[2L]])
    )

    rotation_feature <- .graph_v2_orthogonal(config$attr_dim)
    X_latent <- list(
      latent[[1L]] + matrix(
        stats::rnorm(length(latent[[1L]]), sd = 0.02),
        nrow(latent[[1L]]), config$attr_dim
      ),
      latent[[2L]] %*% rotation_feature + matrix(
        stats::rnorm(length(latent[[2L]]), sd = config$attribute_noise),
        nrow(latent[[2L]]), config$attr_dim
      )
    )

    theta <- if (config$rigid_transform) stats::runif(1L, -pi, pi) else 0
    rotation_position <- matrix(
      c(cos(theta), -sin(theta), sin(theta), cos(theta)),
      2L, 2L, byrow = TRUE
    )
    translation <- if (config$rigid_transform) {
      stats::runif(2L, -2, 2)
    } else {
      c(0, 0)
    }
    target_position_expected <- sweep(
      positions_latent[[2L]] %*% rotation_position,
      2L, translation, "+"
    )
    positions_latent[[2L]] <- target_position_expected + matrix(
      stats::rnorm(length(target_position_expected), sd = config$position_jitter),
      nrow(target_position_expected), 2L
    )

    attributes_latent <- list(
      .graph_v2_mixed_attributes(
        labels_latent[[1L]], positions_latent[[1L]], noise = 0.02
      ),
      .graph_v2_mixed_attributes(
        labels_latent[[2L]], positions_latent[[2L]],
        noise = config$attribute_noise
      )
    )
    attributes_latent[[2L]] <- .graph_v2_corrupt_categories(
      attributes_latent[[2L]], config$attribute_noise
    )
    attributes_latent <- lapply(
      attributes_latent,
      .graph_v2_apply_missing,
      fraction = config$missing_attribute_fraction
    )

    base_shared_edges <- .graph_v2_base_edges(
      shared_positions, shared_labels, config$family, config$directed
    )
    source_edges <- rbind(
      base_shared_edges,
      .graph_v2_private_edges(
        config$n_shared, config$n_private[[1L]], config$directed
      ),
      .graph_v2_misleading_edges(
        labels_latent[[1L]], config$misleading_fraction,
        config$directed, nrow(base_shared_edges)
      )
    )
    source_edges <- .graph_v2_canonical_edges(source_edges, config$directed)
    target_unperturbed <- rbind(
      base_shared_edges,
      .graph_v2_private_edges(
        config$n_shared, config$n_private[[2L]], config$directed
      )
    )
    target_unperturbed <- .graph_v2_canonical_edges(
      target_unperturbed, config$directed
    )
    node_roles <- .graph_v2_roles(config$n_shared, config$n_anchors)
    anchor_latent <- if (config$n_anchors) {
      sort(head(node_roles$train_order, config$n_anchors))
    } else {
      integer()
    }
    label_count <- min(
      length(node_roles$train_order),
      floor(config$n_shared * config$label_fraction)
    )
    labelled_latent <- if (label_count) {
      sort(head(node_roles$train_order, label_count))
    } else {
      integer()
    }
    sizes <- config$n_shared + config$n_private
    row_order <- lapply(sizes, sample.int)
    latent_to_row <- lapply(seq_len(2L), function(domain) {
      match(seq_len(sizes[[domain]]), row_order[[domain]])
    })
    target_state <- .graph_v2_perturb_edges(
      target_unperturbed, labels_latent[[2L]], config
    )

    shared_to_target <- latent_to_row[[2L]][seq_len(config$n_shared)]
    map12 <- rep(NA_integer_, sizes[[1L]])
    source_shared_rows <- latent_to_row[[1L]][seq_len(config$n_shared)]
    map12[source_shared_rows] <- shared_to_target

    roles_latent <- list(
      c(node_roles$role, rep("private", config$n_private[[1L]])),
      c(node_roles$role, rep("private", config$n_private[[2L]]))
    )
    observed_labels <- lapply(seq_len(2L), function(domain) {
      labels <- rep(NA_character_, sizes[[domain]])
      labels[latent_to_row[[domain]][labelled_latent]] <-
        shared_labels[labelled_latent]
      labels
    })
    anchor_vectors <- lapply(seq_len(2L), function(domain) {
      values <- rep(NA_character_, sizes[[domain]])
      values[latent_to_row[[domain]][anchor_latent]] <-
        paste0("anchor-", seq_along(anchor_latent))
      values
    })

    domains <- lapply(seq_len(2L), function(domain) {
      order <- row_order[[domain]]
      private_mask <- order > config$n_shared
      row_id <- ifelse(
        private_mask,
        paste0(
          "domain", domain, "-private-",
          order - config$n_shared
        ),
        paste0("shared-", order)
      )
      edges <- if (domain == 1L) source_edges else target_state$edges
      list(
        X = X_latent[[domain]][order, , drop = FALSE],
        relations = .graph_v2_reindex_edges(edges, latent_to_row[[domain]]),
        positions = positions_latent[[domain]][order, , drop = FALSE],
        attributes = attributes_latent[[domain]][order, , drop = FALSE],
        labels = observed_labels[[domain]],
        row_id = row_id,
        node_role = roles_latent[[domain]][order],
        shared_mask = !private_mask,
        private_mask = private_mask
      )
    })
    names(domains) <- c("domain1", "domain2")
    scaled_X <- .graph_v2_scale_pair(
      domains[[1L]]$X,
      domains[[2L]]$X,
      domains[[1L]]$node_role,
      domains[[2L]]$node_role
    )
    domains[[1L]]$X <- scaled_X$source
    domains[[2L]]$X <- scaled_X$target

    correspondence <- data.frame(
      source_row = source_shared_rows,
      target_row = shared_to_target,
      source_id = domains[[1L]]$row_id[source_shared_rows],
      target_id = domains[[2L]]$row_id[shared_to_target],
      class = shared_labels,
      node_role = node_roles$role,
      stringsAsFactors = FALSE
    )
    correspondence <- correspondence[
      order(correspondence$source_row), , drop = FALSE
    ]
    anchors <- correspondence[
      correspondence$source_row %in% latent_to_row[[1L]][anchor_latent],
      c("source_row", "target_row", "source_id", "target_id"),
      drop = FALSE
    ]
    anchor_idx <- sort(anchors$source_row)

    hd <- as_hyperdesign(
      X_list = lapply(domains, `[[`, "X"),
      labels = anchor_vectors,
      label_name = "anchors"
    )
    shared_target_edges <- target_state$edges[
      target_state$edges$from <= config$n_shared &
        target_state$edges$to <= config$n_shared &
        target_state$edges$relation == "structural",
      , drop = FALSE
    ]
    base_keys <- .graph_v2_edge_key(base_shared_edges, config$directed)
    target_keys <- .graph_v2_edge_key(shared_target_edges, config$directed)
    union_n <- length(union(base_keys, target_keys))
    jaccard <- if (union_n) {
      length(intersect(base_keys, target_keys)) / union_n
    } else {
      1
    }
    missing_count <- sum(vapply(
      attributes_latent,
      function(x) sum(is.na(x)),
      integer(1)
    ))
    attribute_cells <- sum(vapply(
      attributes_latent,
      function(x) nrow(x) * ncol(x),
      integer(1)
    ))

    case <- list(
      protocol_version = .graph_v2_protocol_version,
      case_id = config$case_id,
      split = config$split,
      tier = config$tier,
      tuning_budget = config$tuning_budget,
      config = config,
      domains = domains,
      truth = list(
        map12 = map12,
        correspondence = correspondence,
        shared_labels = shared_labels,
        domain_labels = lapply(seq_len(2L), function(domain) {
          labels_latent[[domain]][row_order[[domain]]]
        }),
        source_shared_mask = domains[[1L]]$shared_mask,
        target_shared_mask = domains[[2L]]$shared_mask,
        source_private_mask = domains[[1L]]$private_mask,
        target_private_mask = domains[[2L]]$private_mask
      ),
      anchors = anchors,
      perturbation = list(
        requested = config[c(
          "edge_delete", "edge_add", "edge_rewire", "weight_noise",
          "weight_scale", "heterophily", "misleading_fraction",
          "position_jitter", "missing_attribute_fraction",
          "attribute_noise"
        )],
        observed = c(
          target_state$counts,
          list(
            edge_delete_fraction = if (target_state$counts$base) {
              target_state$counts$deleted / target_state$counts$base
            } else 0,
            structural_edge_jaccard = jaccard,
            position_jitter_rmse = sqrt(mean(
              (positions_latent[[2L]] - target_position_expected)^2
            )),
            missing_attribute_fraction = missing_count / attribute_cells
          )
        )
      ),
      preprocessing = list(
        id = config$preprocessing_id,
        fit_on = "train_only",
        numeric = c("center", "scale"),
        categorical = "raw_with_explicit_missing",
        graph = "raw_observed_relations_no_truth_features"
      ),
      hd = hd,
      anchors_name = "anchors",
      provenance = list(
        generator = "synthetic_graph_alignment_case_v2",
        generator_version = .graph_v2_protocol_version,
        seed = config$seed,
        config_hash = .graph_v2_hash(config)
      )
    )
    class(case) <- c("graph_alignment_case_v2", "list")
    case$provenance$case_hash <- .graph_v2_hash(
      .graph_v2_case_hash_payload(case)
    )

    # Compatibility aliases used by the current benchmark dispatcher.
    case$X1 <- domains[[1L]]$X
    case$X2 <- domains[[2L]]$X
    case$map12 <- map12
    case$anchor_idx <- anchor_idx
    case$correspondences <- data.frame(
      domain_i = rep(1L, nrow(anchors)),
      index_i = anchors$source_row,
      domain_j = rep(2L, nrow(anchors)),
      index_j = anchors$target_row,
      stringsAsFactors = FALSE
    )
    case$anchors1 <- anchor_vectors[[1L]]
    case$anchors2 <- anchor_vectors[[2L]]
    case$relation_lists <- lapply(domains, function(domain) {
      .graph_v2_relation_list(domain$relations)
    })
    case$benchmark_meta <- list(
      scenario_family = "attributed_graph_v2",
      family = config$family,
      split = config$split,
      tier = config$tier,
      case_id = config$case_id,
      case_hash = case$provenance$case_hash,
      seed = config$seed
    )
    case
  })
}

.graph_v2_empty_problems <- function() {
  data.frame(
    code = character(), field = character(), message = character(),
    stringsAsFactors = FALSE
  )
}

#' Validate a graph benchmark v2 case and its generator truth
#'
#' @keywords internal
#' @noRd
validate_graph_alignment_case_v2 <- function(case, verify_hash = TRUE) {
  problems <- .graph_v2_empty_problems()
  add <- function(code, field, message) {
    problems <<- rbind(
      problems,
      data.frame(
        code = code, field = field, message = message,
        stringsAsFactors = FALSE
      )
    )
  }
  if (!inherits(case, "graph_alignment_case_v2")) {
    add("invalid_class", "case", "Case must inherit graph_alignment_case_v2.")
    return(problems)
  }
  if (!identical(case$protocol_version, .graph_v2_protocol_version)) {
    add("protocol_version", "protocol_version", "Protocol version is not v2.")
  }
  if (!case$split %in% c("development", "final")) {
    add("invalid_split", "split", "Split must be development or final.")
  }
  if (identical(case$split, "final") && case$tuning_budget != 0L) {
    add(
      "final_tuning_budget", "tuning_budget",
      "Final cases must have a zero tuning budget."
    )
  }
  if (!is.list(case$domains) || length(case$domains) != 2L) {
    add("invalid_domains", "domains", "Exactly two domains are required.")
    return(problems)
  }

  required_relation <- c("from", "to", "weight", "relation", "directed")
  for (domain_index in seq_len(2L)) {
    domain <- case$domains[[domain_index]]
    prefix <- paste0("domains[[", domain_index, "]]$")
    n <- nrow(domain$X)
    if (!is.matrix(domain$X) || !is.numeric(domain$X) || !n ||
        any(!is.finite(domain$X))) {
      add("invalid_features", paste0(prefix, "X"), "X must be a finite numeric matrix.")
    }
    if (!is.data.frame(domain$relations) ||
        !all(required_relation %in% names(domain$relations))) {
      add(
        "invalid_relations", paste0(prefix, "relations"),
        "Relations must be a weighted sparse edge-list data frame."
      )
    } else {
      relation <- domain$relations
      if (any(!is.finite(relation$from)) || any(!is.finite(relation$to)) ||
          any(relation$from < 1L | relation$from > n) ||
          any(relation$to < 1L | relation$to > n) ||
          any(relation$from == relation$to) ||
          any(!is.finite(relation$weight)) || any(relation$weight < 0) ||
          any(!nzchar(relation$relation))) {
        add(
          "invalid_relation_values", paste0(prefix, "relations"),
          "Relation endpoints, weights, or names violate the edge-list contract."
        )
      }
    }
    if (!is.matrix(domain$positions) || nrow(domain$positions) != n ||
        ncol(domain$positions) < 2L || any(!is.finite(domain$positions))) {
      add(
        "invalid_positions", paste0(prefix, "positions"),
        "Positions must be a finite matrix with one row per node."
      )
    }
    if (!is.data.frame(domain$attributes) || nrow(domain$attributes) != n ||
        !any(vapply(domain$attributes, is.numeric, logical(1))) ||
        !any(vapply(domain$attributes, is.character, logical(1)))) {
      add(
        "invalid_attributes", paste0(prefix, "attributes"),
        "Attributes must be a mixed data frame with one row per node."
      )
    }
    equal_length <- c(
      labels = length(domain$labels), row_id = length(domain$row_id),
      node_role = length(domain$node_role),
      shared_mask = length(domain$shared_mask),
      private_mask = length(domain$private_mask)
    )
    if (any(equal_length != n)) {
      add(
        "domain_length_mismatch", prefix,
        "Labels, identities, roles, and masks must match the feature rows."
      )
    }
    if (anyDuplicated(domain$row_id)) {
      add("duplicate_row_id", paste0(prefix, "row_id"), "Row identities must be unique within a domain.")
    }
    if (!identical(domain$shared_mask, !domain$private_mask)) {
      add(
        "overlap_mask_contradiction", paste0(prefix, "shared_mask"),
        "Shared and private masks must be exact complements."
      )
    }
  }

  source <- case$domains[[1L]]
  target <- case$domains[[2L]]
  map12 <- case$truth$map12
  if (!is.list(case$truth$domain_labels) ||
      length(case$truth$domain_labels) != 2L ||
      length(case$truth$domain_labels[[1L]]) != nrow(source$X) ||
      length(case$truth$domain_labels[[2L]]) != nrow(target$X)) {
    add(
      "truth_label_length", "truth$domain_labels",
      "Complete generator labels must match both domain row counts."
    )
  }
  if (length(map12) != nrow(source$X)) {
    add("truth_length", "truth$map12", "Truth map must have one entry per source row.")
  } else {
    shared <- source$shared_mask
    valid_shared <- is.finite(map12[shared]) & map12[shared] >= 1L &
      map12[shared] <= nrow(target$X)
    if (!all(valid_shared) || any(!is.na(map12[source$private_mask]))) {
      add(
        "truth_overlap", "truth$map12",
        "Shared rows must map to targets and source-private rows must map to NA."
      )
    } else if (anyDuplicated(map12[shared]) ||
               !all(target$shared_mask[map12[shared]])) {
      add(
        "truth_bijection", "truth$map12",
        "Shared correspondence must be one-to-one onto target-shared rows."
      )
    }
  }
  correspondence <- case$truth$correspondence
  if (!is.data.frame(correspondence) ||
      nrow(correspondence) != sum(source$shared_mask) ||
      any(correspondence$target_row != map12[correspondence$source_row]) ||
      any(correspondence$source_id != source$row_id[correspondence$source_row]) ||
      any(correspondence$target_id != target$row_id[correspondence$target_row])) {
    add(
      "correspondence_mismatch", "truth$correspondence",
      "Correspondence rows, identities, and map12 disagree."
    )
  }
  if (nrow(case$anchors)) {
    expected <- map12[case$anchors$source_row]
    if (any(expected != case$anchors$target_row) ||
        any(source$node_role[case$anchors$source_row] != "train")) {
      add(
        "anchor_truth_mismatch", "anchors",
        "Anchors must be true correspondences drawn only from train nodes."
      )
    }
  }
  if (!identical(case$preprocessing$id, .graph_v2_preprocessing_id) ||
      !identical(case$preprocessing$fit_on, "train_only")) {
    add(
      "preprocessing_contract", "preprocessing",
      "Frozen preprocessing must be the train-only v1 specification."
    )
  }
  if (isTRUE(verify_hash)) {
    expected_config <- .graph_v2_hash(case$config)
    expected_case <- .graph_v2_hash(.graph_v2_case_hash_payload(case))
    if (!identical(case$provenance$config_hash, expected_config)) {
      add("config_hash", "provenance$config_hash", "Configuration hash does not match.")
    }
    if (!identical(case$provenance$case_hash, expected_case)) {
      add("case_hash", "provenance$case_hash", "Case content hash does not match.")
    }
  }
  rownames(problems) <- NULL
  problems
}

#' Score a graph benchmark v2 correspondence without tuning leakage
#'
#' @keywords internal
#' @noRd
score_graph_alignment_case_v2 <- function(
  case,
  prediction,
  purpose = c("report", "tune"),
  evaluation_role = c("all_shared", "dev", "final")
) {
  purpose <- match.arg(purpose)
  evaluation_role <- match.arg(evaluation_role)
  problems <- validate_graph_alignment_case_v2(case)
  if (nrow(problems)) {
    stop(
      "Invalid graph benchmark case: ",
      paste(problems$message, collapse = "; "), call. = FALSE
    )
  }
  if (identical(case$split, "final") && identical(purpose, "tune")) {
    condition <- structure(
      list(
        message = paste0(
          "Final benchmark truth cannot be used for tuning (case ",
          case$case_id, ")."
        ),
        call = NULL,
        case_id = case$case_id
      ),
      class = c(
        "graph_benchmark_final_label_leakage", "error", "condition"
      )
    )
    stop(condition)
  }
  if (!is.numeric(prediction) || length(prediction) != nrow(case$domains[[1L]]$X)) {
    stop("`prediction` must have one target-row index per source row.", call. = FALSE)
  }
  truth <- case$truth$map12
  eligible <- case$domains[[1L]]$shared_mask
  if (!identical(evaluation_role, "all_shared")) {
    eligible <- eligible &
      case$domains[[1L]]$node_role == evaluation_role
  }
  valid <- eligible & is.finite(prediction) & prediction >= 1L &
    prediction <= nrow(case$domains[[2L]]$X)
  correct <- valid & as.integer(prediction) == truth
  data.frame(
    case_id = case$case_id,
    split = case$split,
    evaluation_role = evaluation_role,
    n_eligible = sum(eligible),
    coverage = if (any(eligible)) sum(valid) / sum(eligible) else NA_real_,
    top1_accuracy = if (any(eligible)) {
      sum(correct) / sum(eligible)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

graph_benchmark_control_registry_v2 <- function() {
  data.frame(
    control_id = c(
      "raw_attribute_knn", "position_only", "topology_descriptor_only",
      "anchor_only", "row_identity", "source_only_classifier",
      "target_only_classifier", "chance"
    ),
    task = c(
      rep("node_correspondence", 5L), rep("target_classification", 2L),
      "node_correspondence"
    ),
    observed_topology_used = c(
      FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE
    ),
    uses_target_truth_for_fit = FALSE,
    stringsAsFactors = FALSE
  )
}

.graph_v2_scale_pair <- function(source, target, source_role, target_role) {
  train <- rbind(
    source[source_role == "train", , drop = FALSE],
    target[target_role == "train", , drop = FALSE]
  )
  if (!nrow(train)) train <- rbind(source, target)
  center <- colMeans(train)
  scale <- apply(train, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= .Machine$double.eps] <- 1
  list(
    source = sweep(sweep(source, 2L, center, "-"), 2L, scale, "/"),
    target = sweep(sweep(target, 2L, center, "-"), 2L, scale, "/")
  )
}

.graph_v2_nearest <- function(source, target) {
  if (!nrow(source) || !nrow(target)) return(integer())
  as.integer(RANN::nn2(
    data = target,
    query = source,
    k = 1L
  )$nn.idx[, 1L])
}

.graph_v2_position_prediction <- function(case) {
  source <- case$domains[[1L]]$positions
  target <- case$domains[[2L]]$positions
  anchors <- case$anchors
  if (nrow(anchors) >= 2L) {
    source_anchor <- source[anchors$source_row, , drop = FALSE]
    target_anchor <- target[anchors$target_row, , drop = FALSE]
    source_center <- colMeans(source_anchor)
    target_center <- colMeans(target_anchor)
    decomposition <- svd(crossprod(
      sweep(source_anchor, 2L, source_center, "-"),
      sweep(target_anchor, 2L, target_center, "-")
    ))
    rotation <- decomposition$u %*% t(decomposition$v)
    source <- sweep(source, 2L, source_center, "-") %*% rotation
    source <- sweep(source, 2L, target_center, "+")
  }
  scaled <- .graph_v2_scale_pair(
    source, target,
    case$domains[[1L]]$node_role,
    case$domains[[2L]]$node_role
  )
  .graph_v2_nearest(scaled$source, scaled$target)
}

.graph_v2_topology_descriptors <- function(domain, relation_levels) {
  n <- nrow(domain$X)
  sum_by_node <- function(value, index) {
    out <- numeric(n)
    grouped <- rowsum(value, index, reorder = FALSE)
    out[as.integer(rownames(grouped))] <- grouped[, 1L]
    out
  }
  output <- matrix(
    0, n, 4L * length(relation_levels),
    dimnames = list(NULL, unlist(lapply(relation_levels, function(relation) {
      paste(relation, c("out_degree", "in_degree", "out_weight", "in_weight"), sep = "__")
    })))
  )
  for (relation_index in seq_along(relation_levels)) {
    edges <- domain$relations[
      domain$relations$relation == relation_levels[[relation_index]],
      , drop = FALSE
    ]
    if (!nrow(edges)) next
    offset <- 4L * (relation_index - 1L)
    output[, offset + 1L] <- tabulate(edges$from, nbins = n)
    output[, offset + 2L] <- tabulate(edges$to, nbins = n)
    output[, offset + 3L] <- sum_by_node(edges$weight, edges$from)
    output[, offset + 4L] <- sum_by_node(edges$weight, edges$to)
    output[!is.finite(output)] <- 0
    if (nrow(edges) && !any(edges$directed)) {
      output[, offset + 1L] <- output[, offset + 1L] +
        tabulate(edges$to, nbins = n)
      output[, offset + 2L] <- output[, offset + 1L]
      reverse_weight <- sum_by_node(edges$weight, edges$to)
      output[, offset + 3L] <- output[, offset + 3L] + reverse_weight
      output[, offset + 4L] <- output[, offset + 3L]
    }
  }
  output[!is.finite(output)] <- 0
  output
}

.graph_v2_target_classifier <- function(case, evaluation_role) {
  target <- case$domains[[2L]]
  truth <- case$truth$domain_labels[[2L]]
  known <- !is.na(target$labels)
  classes <- sort(unique(target$labels[known]))
  eligible <- target$node_role == evaluation_role
  if (length(classes) < 2L || !any(eligible)) {
    return(list(
      status = "unsupported", status_code = "insufficient_training_labels",
      prediction = rep(NA_character_, nrow(target$X)),
      score = NA_real_, coverage = 0, n_eligible = sum(eligible)
    ))
  }
  scaled <- .graph_v2_scale_pair(
    target$X, target$X, target$node_role, target$node_role
  )
  X <- scaled$source
  centroids <- do.call(rbind, lapply(classes, function(class) {
    colMeans(X[known & target$labels == class, , drop = FALSE])
  }))
  rownames(centroids) <- classes
  index <- .graph_v2_nearest(X, centroids)
  prediction <- classes[index]
  valid <- eligible & !is.na(prediction)
  list(
    status = "ok", status_code = "ok", prediction = prediction,
    score = if (any(valid)) mean(prediction[valid] == truth[valid]) else NA_real_,
    coverage = if (any(eligible)) sum(valid) / sum(eligible) else NA_real_,
    n_eligible = sum(eligible)
  )
}

.graph_v2_source_classifier <- function(case, evaluation_role) {
  source <- case$domains[[1L]]
  target <- case$domains[[2L]]
  truth <- case$truth$domain_labels[[2L]]
  known <- !is.na(source$labels)
  classes <- sort(unique(source$labels[known]))
  eligible <- target$node_role == evaluation_role
  if (length(classes) < 2L || !any(eligible)) {
    return(list(
      status = "unsupported", status_code = "insufficient_training_labels",
      prediction = rep(NA_character_, nrow(target$X)),
      score = NA_real_, coverage = 0, n_eligible = sum(eligible)
    ))
  }
  scaled <- .graph_v2_scale_pair(
    source$X, target$X, source$node_role, target$node_role
  )
  centroids <- do.call(rbind, lapply(classes, function(class) {
    colMeans(scaled$source[known & source$labels == class, , drop = FALSE])
  }))
  rownames(centroids) <- classes
  index <- .graph_v2_nearest(scaled$target, centroids)
  prediction <- classes[index]
  valid <- eligible & !is.na(prediction)
  list(
    status = "ok", status_code = "ok", prediction = prediction,
    score = if (any(eligible)) {
      sum(valid & prediction == truth) / sum(eligible)
    } else {
      NA_real_
    },
    coverage = if (any(eligible)) sum(valid) / sum(eligible) else NA_real_,
    n_eligible = sum(eligible)
  )
}

benchmark_graph_alignment_controls_v2 <- function(
  case,
  purpose = c("report", "tune"),
  evaluation_role = c("final", "dev", "all_shared")
) {
  purpose <- match.arg(purpose)
  evaluation_role <- match.arg(evaluation_role)
  problems <- validate_graph_alignment_case_v2(case)
  if (nrow(problems)) {
    stop(
      "Invalid graph benchmark case: ",
      paste(problems$message, collapse = "; "), call. = FALSE
    )
  }
  if (identical(case$split, "final") && identical(purpose, "tune")) {
    score_graph_alignment_case_v2(
      case, case$truth$map12,
      purpose = "tune", evaluation_role = "dev"
    )
  }
  source <- case$domains[[1L]]
  target <- case$domains[[2L]]
  scaled_features <- .graph_v2_scale_pair(
    source$X, target$X, source$node_role, target$node_role
  )
  relation_levels <- sort(unique(c(
    source$relations$relation, target$relations$relation
  )))
  source_topology <- .graph_v2_topology_descriptors(source, relation_levels)
  target_topology <- .graph_v2_topology_descriptors(target, relation_levels)
  scaled_topology <- .graph_v2_scale_pair(
    source_topology, target_topology, source$node_role, target$node_role
  )
  anchor_prediction <- rep(NA_integer_, nrow(source$X))
  anchor_prediction[case$anchors$source_row] <- case$anchors$target_row
  identity_prediction <- seq_len(nrow(source$X))
  identity_prediction[identity_prediction > nrow(target$X)] <- NA_integer_
  chance_prediction <- .graph_v2_with_seed(
    case$config$seed + 900001L,
    function() sample.int(
      nrow(target$X), nrow(source$X),
      replace = nrow(source$X) > nrow(target$X)
    )
  )
  predictions <- list(
    raw_attribute_knn = .graph_v2_nearest(
      scaled_features$source, scaled_features$target
    ),
    position_only = .graph_v2_position_prediction(case),
    topology_descriptor_only = .graph_v2_nearest(
      scaled_topology$source, scaled_topology$target
    ),
    anchor_only = anchor_prediction,
    row_identity = identity_prediction,
    chance = chance_prediction
  )
  registry <- graph_benchmark_control_registry_v2()
  rows <- lapply(names(predictions), function(control_id) {
    score <- score_graph_alignment_case_v2(
      case, predictions[[control_id]], purpose, evaluation_role
    )
    metadata <- registry[registry$control_id == control_id, , drop = FALSE]
    data.frame(
      control_id = control_id,
      task = metadata$task,
      status = "ok",
      status_code = "ok",
      observed_topology_used = metadata$observed_topology_used,
      n_eligible = score$n_eligible,
      coverage = score$coverage,
      top1_accuracy = score$top1_accuracy,
      prediction_hash = .graph_v2_hash(predictions[[control_id]]),
      stringsAsFactors = FALSE
    )
  })
  source_classifier <- .graph_v2_source_classifier(case, evaluation_role)
  rows[[length(rows) + 1L]] <- data.frame(
    control_id = "source_only_classifier",
    task = "target_classification",
    status = source_classifier$status,
    status_code = source_classifier$status_code,
    observed_topology_used = FALSE,
    n_eligible = source_classifier$n_eligible,
    coverage = source_classifier$coverage,
    top1_accuracy = source_classifier$score,
    prediction_hash = .graph_v2_hash(source_classifier$prediction),
    stringsAsFactors = FALSE
  )
  target_classifier <- .graph_v2_target_classifier(case, evaluation_role)
  rows[[length(rows) + 1L]] <- data.frame(
    control_id = "target_only_classifier",
    task = "target_classification",
    status = target_classifier$status,
    status_code = target_classifier$status_code,
    observed_topology_used = FALSE,
    n_eligible = target_classifier$n_eligible,
    coverage = target_classifier$coverage,
    top1_accuracy = target_classifier$score,
    prediction_hash = .graph_v2_hash(target_classifier$prediction),
    stringsAsFactors = FALSE
  )
  out <- do.call(rbind, rows)
  out <- out[match(registry$control_id, out$control_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.graph_v2_manifest_hashes <- c(
  development = "6030e2e2d49b0473bba78ed283ea2fe6",
  final = "2f79e4b34867edba288239fbc427a437"
)

.graph_v2_manifest_file <- function(split) {
  split <- match.arg(split, c("development", "final"))
  filename <- paste0("graph_alignment_v2_", split, "_manifest.csv")
  path <- system.file("benchmarks", filename, package = "manifoldalign")
  if (!nzchar(path)) path <- file.path("inst", "benchmarks", filename)
  if (!file.exists(path)) {
    stop("Frozen graph benchmark manifest is missing: ", filename, call. = FALSE)
  }
  path
}

validate_graph_alignment_manifest_v2 <- function(manifest, split) {
  required <- c(
    "protocol_version", "case_id", "split", "tier", "family",
    "n_shared", "n_private_1", "n_private_2", "attr_dim", "n_anchors",
    "label_fraction", "edge_delete", "edge_add", "edge_rewire",
    "weight_noise", "weight_scale", "directed", "heterophily",
    "misleading_fraction", "position_jitter", "rigid_transform",
    "missing_attribute_fraction", "attribute_noise", "class_imbalance",
    "private_classes", "preprocessing_id", "tuning_budget", "seed"
  )
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop("Manifest is missing fields: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(manifest$case_id) || anyDuplicated(manifest$seed)) {
    stop("Manifest case IDs and seeds must be unique.", call. = FALSE)
  }
  if (!all(manifest$protocol_version == .graph_v2_protocol_version) ||
      !all(manifest$split == split) ||
      !all(manifest$preprocessing_id == .graph_v2_preprocessing_id)) {
    stop("Manifest protocol, split, or preprocessing metadata is inconsistent.", call. = FALSE)
  }
  required_tiers <- c(
    "oracle", "small_accuracy", "medium_robustness", "large_scaling"
  )
  required_families <- c("sbm", "geometric", "scale_free")
  if (!setequal(unique(manifest$tier), required_tiers) ||
      !setequal(unique(manifest$family), required_families)) {
    stop("Manifest must cover every frozen tier and graph family.", call. = FALSE)
  }
  if (identical(split, "final")) {
    accuracy_cases <- manifest$tier %in% c(
      "small_accuracy", "medium_robustness"
    )
    if (nrow(manifest) != 42L || sum(accuracy_cases) != 30L ||
        any(manifest$tuning_budget != 0L)) {
      stop(
        "The frozen final manifest requires exactly 42 cases, including ",
        "30 accuracy/robustness cases, and zero tuning budget.",
        call. = FALSE
      )
    }
  } else if (any(manifest$tuning_budget < 1L)) {
    stop("Development cases require a positive tuning budget.", call. = FALSE)
  }
  for (row_index in seq_len(nrow(manifest))) {
    row <- manifest[row_index, , drop = FALSE]
    tryCatch(
      .graph_v2_config(
        row$n_shared, c(row$n_private_1, row$n_private_2), row$family,
        row$attr_dim, row$n_anchors, row$label_fraction, row$edge_delete,
        row$edge_add, row$edge_rewire, row$weight_noise, row$weight_scale,
        row$directed, row$heterophily, row$misleading_fraction,
        row$position_jitter, row$rigid_transform,
        row$missing_attribute_fraction, row$attribute_noise,
        row$class_imbalance, row$private_classes, row$split, row$tier,
        row$preprocessing_id, row$tuning_budget, row$seed, row$case_id
      ),
      error = function(error) {
        stop(
          "Invalid manifest row ", row_index, " (`", row$case_id,
          "`): ", conditionMessage(error), call. = FALSE
        )
      }
    )
  }
  invisible(TRUE)
}

graph_alignment_benchmark_manifest_v2 <- function(
  split = c("development", "final"),
  verify_hash = TRUE
) {
  split <- match.arg(split)
  path <- .graph_v2_manifest_file(split)
  observed_hash <- unname(tools::md5sum(path))
  if (isTRUE(verify_hash) &&
      !identical(observed_hash, .graph_v2_manifest_hashes[[split]])) {
    stop(
      "Frozen ", split, " manifest hash mismatch: expected ",
      .graph_v2_manifest_hashes[[split]], ", observed ", observed_hash, ".",
      call. = FALSE
    )
  }
  manifest <- utils::read.csv(path, stringsAsFactors = FALSE)
  validate_graph_alignment_manifest_v2(manifest, split)
  attr(manifest, "manifest_hash") <- observed_hash
  attr(manifest, "manifest_path") <- path
  manifest
}

graph_alignment_case_from_manifest_v2 <- function(
  case_id,
  split = c("development", "final")
) {
  split <- match.arg(split)
  manifest <- graph_alignment_benchmark_manifest_v2(split)
  row <- manifest[manifest$case_id == case_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Expected exactly one manifest row for case: ", case_id, call. = FALSE)
  }
  synthetic_graph_alignment_case_v2(
    n_shared = row$n_shared,
    n_private = c(row$n_private_1, row$n_private_2),
    family = row$family,
    attr_dim = row$attr_dim,
    n_anchors = row$n_anchors,
    label_fraction = row$label_fraction,
    edge_delete = row$edge_delete,
    edge_add = row$edge_add,
    edge_rewire = row$edge_rewire,
    weight_noise = row$weight_noise,
    weight_scale = row$weight_scale,
    directed = row$directed,
    heterophily = row$heterophily,
    misleading_fraction = row$misleading_fraction,
    position_jitter = row$position_jitter,
    rigid_transform = row$rigid_transform,
    missing_attribute_fraction = row$missing_attribute_fraction,
    attribute_noise = row$attribute_noise,
    class_imbalance = row$class_imbalance,
    private_classes = row$private_classes,
    split = row$split,
    tier = row$tier,
    preprocessing_id = row$preprocessing_id,
    tuning_budget = row$tuning_budget,
    seed = row$seed,
    case_id = row$case_id
  )
}
