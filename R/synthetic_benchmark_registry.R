.bench_or <- function(x, y) if (is.null(x)) y else x

.rand_rotation_benchmark <- function(d, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  Q
}

.synthetic_bundle_views <- function(view_mats,
                                    latent_by_view,
                                    true_labels_by_view,
                                    observed_labels_by_view = true_labels_by_view,
                                    row_ids = NULL,
                                    extras = list()) {
  if (!is.list(view_mats) || !length(view_mats)) {
    stop("`view_mats` must be a non-empty named list.", call. = FALSE)
  }
  if (is.null(names(view_mats)) || any(!nzchar(names(view_mats)))) {
    names(view_mats) <- paste0("view", seq_along(view_mats))
  }

  view_names <- names(view_mats)
  if (is.null(names(latent_by_view))) names(latent_by_view) <- view_names
  if (is.null(names(true_labels_by_view))) names(true_labels_by_view) <- view_names
  if (is.null(names(observed_labels_by_view))) names(observed_labels_by_view) <- view_names

  if (is.null(row_ids)) {
    row_ids <- stats::setNames(lapply(view_mats, function(X) seq_len(nrow(X))), view_names)
  } else if (is.null(names(row_ids))) {
    names(row_ids) <- view_names
  }

  primary_view <- view_names[[1L]]
  out <- c(
    list(
      latent = latent_by_view[[primary_view]],
      true_labels = true_labels_by_view[[primary_view]],
      observed_labels = observed_labels_by_view[[primary_view]],
      latent_by_view = latent_by_view,
      true_labels_by_view = true_labels_by_view,
      observed_labels_by_view = observed_labels_by_view,
      row_ids = row_ids
    ),
    view_mats,
    extras
  )
  out
}

.scenario_row_ids <- function(scenario, views = NULL) {
  view_names <- .bench_or(views, .bench_or(scenario$view_names, grep("^view[0-9]+$", names(scenario), value = TRUE)))
  row_ids <- scenario$row_ids
  if (is.null(row_ids)) {
    row_ids <- stats::setNames(lapply(view_names, function(view_name) seq_len(nrow(as.matrix(scenario[[view_name]])))), view_names)
  }
  row_ids[view_names]
}

.scenario_labels_by_view <- function(scenario, kind = c("observed", "true"), views = NULL) {
  kind <- match.arg(kind)
  view_names <- .bench_or(views, .bench_or(scenario$view_names, grep("^view[0-9]+$", names(scenario), value = TRUE)))
  label_obj <- if (identical(kind, "observed")) {
    .bench_or(scenario$observed_labels_by_view, .bench_or(scenario$observed_labels, .bench_or(scenario$true_labels_by_view, scenario$true_labels)))
  } else {
    .bench_or(scenario$true_labels_by_view, scenario$true_labels)
  }

  if (is.list(label_obj) && !is.null(names(label_obj))) {
    return(label_obj[view_names])
  }
  stats::setNames(lapply(view_names, function(view_name) label_obj), view_names)
}

.synthetic_linear_triplet <- function(n = 60L, noise = 0.01, seed = 1L) {
  set.seed(seed)

  n1 <- floor(n / 2)
  n2 <- n - n1

  latent_c1 <- matrix(rnorm(n1 * 2, mean = c(-0.5, -0.5), sd = 0.3), n1, 2)
  latent_c2 <- matrix(rnorm(n2 * 2, mean = c(0.5, 0.5), sd = 0.3), n2, 2)
  latent <- rbind(latent_c1, latent_c2)
  colnames(latent) <- c("u1", "u2")

  true_labels <- factor(c(rep("class1", n1), rep("class2", n2)))

  R1 <- .rand_rotation_benchmark(2, seed + 1L)
  R2 <- .rand_rotation_benchmark(2, seed + 2L)
  R3 <- .rand_rotation_benchmark(2, seed + 3L)
  b1 <- c(0.2, -0.4)
  b2 <- c(-0.7, 0.6)
  b3 <- c(0.1, 0.8)

  X1 <- latent %*% diag(c(1.0, 0.7)) %*% R1 +
    matrix(b1, n, 2, TRUE) + rnorm(n * 2, 0, noise)
  X2 <- latent %*% diag(c(0.5, 1.2)) %*% R2 +
    matrix(b2, n, 2, TRUE) + rnorm(n * 2, 0, noise)
  X3 <- latent %*% diag(c(1.5, 0.4)) %*% R3 +
    matrix(b3, n, 2, TRUE) + rnorm(n * 2, 0, noise)

  .synthetic_bundle_views(
    view_mats = list(view1 = X1, view2 = X2, view3 = X3),
    latent_by_view = list(view1 = latent, view2 = latent, view3 = latent),
    true_labels_by_view = list(view1 = true_labels, view2 = true_labels, view3 = true_labels),
    extras = list(
      class_centers = list(c(-0.5, -0.5), c(0.5, 0.5))
    )
  )
}

.synthetic_isometric_triplet <- function(n = 50L, noise = 0.02, seed = 2L) {
  set.seed(seed)
  t <- sort(runif(n, 0, 2 * pi))
  true_labels <- factor(ifelse(t < pi, "early", "late"))

  v1 <- cbind(t, 0)
  v2 <- cbind(cos(t), sin(t))
  v3 <- cbind(cos(t), sin(t), t / (2 * pi))

  v1 <- v1 + matrix(rnorm(n * 2, 0, noise), n, 2)
  v2 <- v2 + matrix(rnorm(n * 2, 0, noise), n, 2)
  v3 <- v3 + matrix(rnorm(n * 3, 0, noise), n, 3)

  .synthetic_bundle_views(
    view_mats = list(view1 = v1, view2 = v2, view3 = v3),
    latent_by_view = list(view1 = t, view2 = t, view3 = t),
    true_labels_by_view = list(view1 = true_labels, view2 = true_labels, view3 = true_labels),
    extras = list(parameter_boundary = pi)
  )
}

.synthetic_hard_triplet <- function(n_side = 10L, noise = 0.03, seed = 3L) {
  set.seed(seed)
  u <- seq(-1, 1, length.out = n_side)
  grid <- as.matrix(expand.grid(u, u))
  true_labels <- factor(ifelse(grid[, 1] < 0, "left", "right"))

  v1 <- grid + matrix(rnorm(nrow(grid) * 2, 0, noise), ncol = 2)

  theta <- (grid[, 1] + 1) * 2 * pi
  height <- grid[, 2]
  v2 <- cbind(theta * cos(theta), height, theta * sin(theta)) +
    matrix(rnorm(nrow(grid) * 3, 0, noise), ncol = 3)

  A <- matrix(c(2, 0.5, 0, 0.2), 2, 2)
  v3 <- grid %*% A
  idx_hard <- which(grid[, 1] < 0 & grid[, 2] < 0)
  v3[idx_hard, ] <- v3[idx_hard, ] +
    matrix(rnorm(length(idx_hard) * 2, 0, 0.2), ncol = 2)
  v3 <- v3 + matrix(rnorm(nrow(grid) * 2, 0, noise), ncol = 2)

  .synthetic_bundle_views(
    view_mats = list(view1 = v1, view2 = v2, view3 = v3),
    latent_by_view = list(view1 = grid, view2 = grid, view3 = grid),
    true_labels_by_view = list(view1 = true_labels, view2 = true_labels, view3 = true_labels),
    extras = list(hard_indices = idx_hard)
  )
}

.synthetic_partial_overlap_triplet <- function(n_shared = 36L,
                                               n_private = 12L,
                                               noise = 0.02,
                                               seed = 4L) {
  set.seed(seed)

  mk_latent <- function(n) matrix(rnorm(n * 2, sd = 0.6), ncol = 2)
  shared <- mk_latent(n_shared)
  private1 <- mk_latent(n_private)
  private2 <- mk_latent(n_private)
  private3 <- mk_latent(max(2L, floor(n_private / 2L)))

  make_labels <- function(latent) factor(ifelse(latent[, 1] + 0.35 * latent[, 2] > 0, "shared_right", "shared_left"))

  row_ids <- list(
    view1 = c(sprintf("shared_%03d", seq_len(n_shared)), sprintf("view1_private_%03d", seq_len(n_private))),
    view2 = c(sprintf("shared_%03d", seq_len(n_shared)), sprintf("view2_private_%03d", seq_len(n_private))),
    view3 = c(sprintf("shared_%03d", seq_len(n_shared)), sprintf("view3_private_%03d", seq_len(length(private3))))
  )

  latent1 <- rbind(shared, private1)
  latent2 <- rbind(shared, private2)
  latent3 <- rbind(shared, private3)

  R1 <- .rand_rotation_benchmark(2, seed + 11L)
  R2 <- .rand_rotation_benchmark(2, seed + 12L)
  R3 <- .rand_rotation_benchmark(2, seed + 13L)

  view1 <- latent1 %*% diag(c(1.0, 0.6)) %*% R1 + matrix(rnorm(length(latent1), sd = noise), ncol = 2)
  view2 <- latent2 %*% diag(c(0.8, 1.1)) %*% R2 + matrix(rnorm(length(latent2), sd = noise), ncol = 2)
  view3 <- latent3 %*% diag(c(1.2, 0.7)) %*% R3 + matrix(rnorm(length(latent3), sd = noise), ncol = 2)

  .synthetic_bundle_views(
    view_mats = list(view1 = view1, view2 = view2, view3 = view3),
    latent_by_view = list(view1 = latent1, view2 = latent2, view3 = latent3),
    true_labels_by_view = list(
      view1 = make_labels(latent1),
      view2 = make_labels(latent2),
      view3 = make_labels(latent3)
    ),
    row_ids = row_ids,
    extras = list(
      n_shared = n_shared,
      n_private = n_private,
      overlap_fraction = n_shared / max(length(row_ids$view1), length(row_ids$view2))
    )
  )
}

.synthetic_noisy_label_triplet <- function(n = 60L,
                                           noise = 0.015,
                                           label_noise = 0.2,
                                           seed = 5L) {
  base <- .synthetic_linear_triplet(n = n, noise = noise, seed = seed)
  view_names <- c("view1", "view2", "view3")

  corrupt <- function(labels, seed_offset) {
    set.seed(seed + seed_offset)
    labels_chr <- as.character(labels)
    lv <- sort(unique(labels_chr))
    n_flip <- max(1L, floor(length(labels_chr) * label_noise))
    flip_idx <- sample(seq_along(labels_chr), n_flip, replace = FALSE)
    labels_chr[flip_idx] <- ifelse(labels_chr[flip_idx] == lv[1L], lv[2L], lv[1L])
    factor(labels_chr, levels = lv)
  }

  observed_labels_by_view <- list(
    view1 = corrupt(base$true_labels_by_view$view1, 1L),
    view2 = corrupt(base$true_labels_by_view$view2, 2L),
    view3 = corrupt(base$true_labels_by_view$view3, 3L)
  )

  base$observed_labels <- observed_labels_by_view$view1
  base$observed_labels_by_view <- observed_labels_by_view
  base$label_noise_rate <- label_noise
  base
}

.synthetic_many_domain_suite <- function(n = 54L,
                                         noise = 0.02,
                                         n_views = 4L,
                                         seed = 6L) {
  set.seed(seed)
  latent <- matrix(rnorm(n * 2), ncol = 2)
  labels <- factor(ifelse(latent[, 1] * latent[, 2] > 0, "diagonal", "off_diagonal"))

  view_mats <- vector("list", n_views)
  names(view_mats) <- paste0("view", seq_len(n_views))
  for (i in seq_len(n_views)) {
    Ri <- .rand_rotation_benchmark(2, seed + i)
    scale_i <- diag(c(0.7 + 0.15 * i, 1.1 - 0.08 * (i - 1L)))
    offset_i <- matrix(cos(i), n, 2) * matrix(c(0.15, -0.1), n, 2, byrow = TRUE)
    warp_i <- latent
    if (i %% 2L == 0L) {
      warp_i[, 2] <- warp_i[, 2] + 0.2 * latent[, 1]^2
    }
    view_mats[[i]] <- warp_i %*% scale_i %*% Ri + offset_i + matrix(rnorm(n * 2, sd = noise), ncol = 2)
  }

  .synthetic_bundle_views(
    view_mats = view_mats,
    latent_by_view = stats::setNames(replicate(n_views, latent, simplify = FALSE), names(view_mats)),
    true_labels_by_view = stats::setNames(replicate(n_views, labels, simplify = FALSE), names(view_mats)),
    extras = list(n_domains = n_views)
  )
}

.synthetic_highdim_triplet <- function(n = 48L,
                                       dims = c(96L, 128L, 160L),
                                       noise = 0.03,
                                       seed = 7L) {
  set.seed(seed)
  latent <- matrix(rnorm(n * 3), ncol = 3)
  labels <- factor(ifelse(latent[, 1] + latent[, 2] > 0, "high", "low"))

  view_mats <- lapply(seq_along(dims), function(i) {
    Wi <- matrix(rnorm(ncol(latent) * dims[i], sd = 0.6), nrow = ncol(latent))
    Xi <- latent %*% Wi
    Xi + matrix(rnorm(n * dims[i], sd = noise), nrow = n)
  })
  names(view_mats) <- paste0("view", seq_along(dims))

  .synthetic_bundle_views(
    view_mats = view_mats,
    latent_by_view = stats::setNames(replicate(length(dims), latent, simplify = FALSE), names(view_mats)),
    true_labels_by_view = stats::setNames(replicate(length(dims), labels, simplify = FALSE), names(view_mats)),
    extras = list(feature_dims = dims)
  )
}

.synthetic_benchmark_specs <- function(profile = c("full", "fast")) {
  profile <- match.arg(profile)

  full_specs <- list(
    linear_affine = list(
      generator = .synthetic_linear_triplet,
      defaults = list(n = 60L, noise = 0.01, seed = 1L),
      family = "feature",
      latent_relation = "linear",
      supervision = "labels",
      side_information = "class_labels",
      side_information_quality = "clean",
      overlap = "full",
      n_views = 3L,
      n_samples = 60L,
      difficulty = "easy",
      notes = "Affine multi-view Gaussian blobs with stable class structure."
    ),
    isometric_curve = list(
      generator = .synthetic_isometric_triplet,
      defaults = list(n = 50L, noise = 0.02, seed = 2L),
      family = "feature",
      latent_relation = "nonlinear_isometric",
      supervision = "labels",
      side_information = "progression_labels",
      side_information_quality = "clean",
      overlap = "full",
      n_views = 3L,
      n_samples = 50L,
      difficulty = "medium",
      notes = "Shared 1D latent parameter rendered as line, circle, and helix."
    ),
    hard_nonisometric = list(
      generator = .synthetic_hard_triplet,
      defaults = list(n_side = 10L, noise = 0.03, seed = 3L),
      family = "feature",
      latent_relation = "nonlinear_nonisometric",
      supervision = "labels",
      side_information = "spatial_labels",
      side_information_quality = "clean",
      overlap = "full",
      n_views = 3L,
      n_samples = 100L,
      difficulty = "hard",
      notes = "Grid latent space with nonlinear distortion, density skew, and heteroscedastic corruption."
    ),
    partial_open_set = list(
      generator = .synthetic_partial_overlap_triplet,
      defaults = list(n_shared = 48L, n_private = 18L, noise = 0.02, seed = 4L),
      family = "feature",
      latent_relation = "linear",
      supervision = "labels",
      side_information = "class_labels",
      side_information_quality = "clean",
      overlap = "partial_open_set",
      n_views = 3L,
      n_samples = 66L,
      difficulty = "hard",
      notes = "Shared latent structure with private rows in each view and only partial row overlap."
    ),
    noisy_label_shift = list(
      generator = .synthetic_noisy_label_triplet,
      defaults = list(n = 60L, noise = 0.015, label_noise = 0.2, seed = 5L),
      family = "feature",
      latent_relation = "linear",
      supervision = "labels",
      side_information = "class_labels",
      side_information_quality = "noisy",
      overlap = "full",
      n_views = 3L,
      n_samples = 60L,
      difficulty = "medium",
      notes = "Linear shared structure with corrupted label supervision used during fitting."
    ),
    many_domain_consensus = list(
      generator = .synthetic_many_domain_suite,
      defaults = list(n = 54L, noise = 0.02, n_views = 4L, seed = 6L),
      family = "feature",
      latent_relation = "nonlinear_isometric",
      supervision = "labels",
      side_information = "class_labels",
      side_information_quality = "clean",
      overlap = "full",
      n_views = 4L,
      n_samples = 54L,
      difficulty = "medium",
      notes = "Four-domain consensus alignment problem with consistent row IDs across all views."
    ),
    highdim_stress = list(
      generator = .synthetic_highdim_triplet,
      defaults = list(n = 48L, dims = c(96L, 128L, 160L), noise = 0.03, seed = 7L),
      family = "feature",
      latent_relation = "linear",
      supervision = "labels",
      side_information = "class_labels",
      side_information_quality = "clean",
      overlap = "full",
      n_views = 3L,
      n_samples = 48L,
      difficulty = "hard",
      notes = "High-dimensional p >> n stress test with mixed feature dimensions across views."
    )
  )

  fast_specs <- list(
    linear_affine = utils::modifyList(full_specs$linear_affine, list(defaults = list(n = 20L, noise = 0.01, seed = 1L), n_samples = 20L)),
    isometric_curve = utils::modifyList(full_specs$isometric_curve, list(defaults = list(n = 25L, noise = 0.02, seed = 2L), n_samples = 25L)),
    hard_nonisometric = utils::modifyList(full_specs$hard_nonisometric, list(defaults = list(n_side = 5L, noise = 0.03, seed = 3L), n_samples = 25L)),
    partial_open_set = utils::modifyList(full_specs$partial_open_set, list(defaults = list(n_shared = 18L, n_private = 6L, noise = 0.02, seed = 4L), n_samples = 24L)),
    noisy_label_shift = utils::modifyList(full_specs$noisy_label_shift, list(defaults = list(n = 24L, noise = 0.015, label_noise = 0.2, seed = 5L), n_samples = 24L))
  )

  if (identical(profile, "full")) full_specs else fast_specs
}

synthetic_alignment_scenarios <- function(profile = c("full", "fast")) {
  profile <- match.arg(profile)
  specs <- .synthetic_benchmark_specs(profile)

  rows <- lapply(names(specs), function(name) {
    spec <- specs[[name]]
    data.frame(
      scenario = name,
      profile = profile,
      family = spec$family,
      latent_relation = spec$latent_relation,
      supervision = spec$supervision,
      side_information = spec$side_information,
      side_information_quality = spec$side_information_quality,
      overlap = spec$overlap,
      difficulty = spec$difficulty,
      n_views = as.integer(spec$n_views),
      n_samples = as.integer(spec$n_samples),
      default_seed = as.integer(spec$defaults$seed),
      notes = spec$notes,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

synthetic_alignment_scenario <- function(scenario_name, profile = c("full", "fast"), seed = NULL, ...) {
  profile <- match.arg(profile)
  specs <- .synthetic_benchmark_specs(profile)
  scenario_name <- as.character(scenario_name)[1L]
  if (!scenario_name %in% names(specs)) {
    stop(
      "Unknown synthetic benchmark scenario '", scenario_name, "'. Available scenarios: ",
      paste(names(specs), collapse = ", "),
      call. = FALSE
    )
  }
  spec <- specs[[scenario_name]]

  args <- spec$defaults
  overrides <- list(...)
  if (!is.null(seed)) {
    overrides$seed <- seed
  }
  if (length(overrides)) {
    args[names(overrides)] <- overrides
  }

  scenario <- do.call(spec$generator, args)
  scenario$view_names <- grep("^view[0-9]+$", names(scenario), value = TRUE)
  scenario$benchmark_meta <- list(
    scenario = scenario_name,
    profile = profile,
    family = spec$family,
    latent_relation = spec$latent_relation,
    supervision = spec$supervision,
    side_information = spec$side_information,
    side_information_quality = spec$side_information_quality,
    overlap = spec$overlap,
    n_views = spec$n_views,
    n_samples = spec$n_samples,
    difficulty = spec$difficulty,
    generator_args = args,
    notes = spec$notes
  )
  scenario
}

synthetic_alignment_to_hyperdesign <- function(scenario,
                                               views = NULL,
                                               label = NULL,
                                               label_col = "label",
                                               id_col = "sample_id") {
  view_names <- .bench_or(scenario$view_names, grep("^view[0-9]+$", names(scenario), value = TRUE))
  if (is.null(views)) {
    views <- view_names
  } else {
    views <- intersect(as.character(views), view_names)
  }

  if (!length(views)) {
    stop("No valid views selected for synthetic scenario.", call. = FALSE)
  }

  label_by_view <- if (is.null(label)) {
    .scenario_labels_by_view(scenario, kind = "observed", views = views)
  } else if (is.list(label) && !is.null(names(label))) {
    label[views]
  } else {
    stats::setNames(lapply(views, function(view_name) label), views)
  }
  row_ids <- .scenario_row_ids(scenario, views = views)
  domains <- lapply(views, function(view_name) {
    X <- as.matrix(scenario[[view_name]])
    design_df <- data.frame(sample_id = row_ids[[view_name]], stringsAsFactors = FALSE)
    names(design_df)[1] <- id_col
    labels_i <- label_by_view[[view_name]]
    if (!is.null(labels_i)) {
      if (length(labels_i) != nrow(X)) {
        stop("Label vector for `", view_name, "` has length ", length(labels_i), " but expected ", nrow(X), ".", call. = FALSE)
      }
      design_df[[label_col]] <- labels_i
    }

    if (requireNamespace("multidesign", quietly = TRUE)) {
      return(multidesign::multidesign(X, design_df))
    }

    list(x = X, design = design_df)
  })
  names(domains) <- views

  if (requireNamespace("multidesign", quietly = TRUE)) {
    return(multidesign::hyperdesign(domains))
  }

  class(domains) <- c("hyperdesign", "list")
  domains
}

.benchmark_metrics_to_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.data.frame(x)) {
    if (nrow(x) != 1L) {
      stop("`score_fn` must return a named list/vector or a one-row data.frame.", call. = FALSE)
    }
    return(as.list(x[1, , drop = FALSE]))
  }
  if (is.atomic(x) && !is.null(names(x))) {
    return(as.list(x))
  }
  if (is.list(x) && !is.null(names(x))) {
    return(x)
  }
  stop("`score_fn` must return a named list/vector or a one-row data.frame.", call. = FALSE)
}

.combine_benchmark_rows <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  normalized <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    if (length(missing)) {
      row[missing] <- NA
    }
    row[all_names]
  })
  do.call(rbind, lapply(normalized, function(row) {
    as.data.frame(row, stringsAsFactors = FALSE)
  }))
}

.normalize_benchmark_method_meta <- function(meta, method_name) {
  if (is.null(meta)) {
    return(list())
  }
  if (!is.list(meta) || is.null(names(meta)) || any(!nzchar(names(meta)))) {
    stop("Benchmark method metadata for `", method_name, "` must be a named list.", call. = FALSE)
  }

  bad <- names(meta)[vapply(meta, function(x) {
    is.null(x) || length(x) != 1L || (is.list(x) && !is.atomic(x))
  }, logical(1))]
  if (length(bad)) {
    stop(
      "Benchmark method metadata entries must be scalar atomic values. Invalid entries for `",
      method_name, "`: ", paste(bad, collapse = ", "),
      call. = FALSE
    )
  }

  meta
}

.normalize_benchmark_method_entry <- function(entry, method_name) {
  if (is.function(entry)) {
    return(list(fit = entry, meta = list()))
  }

  if (!is.list(entry)) {
    stop("Benchmark method `", method_name, "` must be a function or a list with `$fit`.", call. = FALSE)
  }

  fit <- .bench_or(entry$fit, .bench_or(entry$fn, entry$method))
  if (!is.function(fit)) {
    stop("Benchmark method `", method_name, "` must provide a callable `$fit`.", call. = FALSE)
  }

  meta <- .bench_or(entry$meta, .bench_or(entry$metadata, list()))
  list(fit = fit, meta = .normalize_benchmark_method_meta(meta, method_name))
}

.normalize_benchmark_methods <- function(methods) {
  if (is.function(methods)) {
    methods <- list(method = methods)
  }
  if (!is.list(methods) || !length(methods)) {
    stop("`methods` must be a non-empty named list of functions or method specs.", call. = FALSE)
  }
  if (is.null(names(methods)) || any(!nzchar(names(methods)))) {
    stop("`methods` must be named.", call. = FALSE)
  }

  stats::setNames(
    lapply(names(methods), function(method_name) {
      .normalize_benchmark_method_entry(methods[[method_name]], method_name)
    }),
    names(methods)
  )
}

.benchmark_seed_from_parts <- function(...) {
  txt <- paste(..., collapse = "::")
  ints <- utf8ToInt(txt)
  if (!length(ints)) return(1L)
  acc <- sum(as.numeric(ints) * seq_along(ints))
  as.integer(1L + (acc %% (.Machine$integer.max - 1)))
}

.as_numeric_column <- function(x) {
  if (is.factor(x)) {
    return(as.numeric(as.character(x)))
  }
  as.numeric(x)
}

.aggregate_stat_matrix <- function(x, stat_names) {
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = length(stat_names), byrow = TRUE)
  } else {
    x <- as.matrix(x)
  }
  colnames(x) <- stat_names
  x
}

.benchmark_group_key <- function(df, group_cols) {
  parts <- lapply(group_cols, function(col) {
    x <- df[[col]]
    if (is.factor(x)) {
      x <- as.character(x)
    }
    if (is.logical(x)) {
      x <- ifelse(is.na(x), "__NA__", ifelse(x, "__TRUE__", "__FALSE__"))
    } else {
      x <- ifelse(is.na(x), "__NA__", as.character(x))
    }
    x
  })
  do.call(paste, c(parts, sep = "||"))
}

synthetic_benchmark_method_catalog <- function(methods) {
  method_specs <- .normalize_benchmark_methods(methods)
  rows <- lapply(names(method_specs), function(method_name) {
    c(
      list(method = method_name),
      method_specs[[method_name]]$meta
    )
  })
  .combine_benchmark_rows(rows)
}

run_synthetic_alignment_benchmark <- function(methods,
                                              scenarios = NULL,
                                              profile = c("full", "fast"),
                                              seeds = 1L,
                                              score_fn = NULL,
                                              verbose = FALSE) {
  profile <- match.arg(profile)

  method_specs <- .normalize_benchmark_methods(methods)

  if (is.null(scenarios)) {
    scenario_names <- synthetic_alignment_scenarios(profile)$scenario
  } else if (is.data.frame(scenarios)) {
    if (!"scenario" %in% names(scenarios)) {
      stop("When `scenarios` is a data.frame it must contain a `scenario` column.", call. = FALSE)
    }
    scenario_names <- as.character(scenarios$scenario)
  } else {
    scenario_names <- as.character(scenarios)
  }

  seeds <- as.integer(seeds)
  rows <- list()
  row_idx <- 1L

  for (scenario_name in scenario_names) {
    for (seed in seeds) {
      scenario <- synthetic_alignment_scenario(scenario_name, profile = profile, seed = seed)
      scenario_meta <- .bench_or(scenario$benchmark_meta, list())
      scenario_cols <- list(
        scenario_family = .bench_or(scenario_meta$family, NA_character_),
        scenario_latent_relation = .bench_or(scenario_meta$latent_relation, NA_character_),
        scenario_supervision = .bench_or(scenario_meta$supervision, NA_character_),
        scenario_side_information = .bench_or(scenario_meta$side_information, NA_character_),
        scenario_side_information_quality = .bench_or(scenario_meta$side_information_quality, NA_character_),
        scenario_overlap = .bench_or(scenario_meta$overlap, NA_character_),
        scenario_n_views = .bench_or(scenario_meta$n_views, NA_integer_),
        scenario_n_samples = .bench_or(scenario_meta$n_samples, NA_integer_),
        scenario_difficulty = .bench_or(scenario_meta$difficulty, NA_character_)
      )

      for (method_name in names(method_specs)) {
        if (isTRUE(verbose)) {
          message(
            "run_synthetic_alignment_benchmark: scenario=", scenario_name,
            " seed=", seed,
            " method=", method_name
          )
        }

        fit <- NULL
        error_msg <- NA_character_
        method_seed <- .benchmark_seed_from_parts(profile, scenario_name, seed, method_name)
        method_spec <- method_specs[[method_name]]
        runtime_sec <- system.time({
          set.seed(method_seed)
          fit <- tryCatch(
            method_spec$fit(scenario),
            error = function(e) {
              error_msg <<- conditionMessage(e)
              NULL
            }
          )
        })[["elapsed"]]

        metrics <- list()
        if (is.na(error_msg) && !is.null(score_fn)) {
          metrics <- tryCatch(
            .benchmark_metrics_to_list(score_fn(fit, scenario, method_name)),
            error = function(e) {
              error_msg <<- paste0("score_fn failed: ", conditionMessage(e))
              list()
            }
          )
        }

        rows[[row_idx]] <- c(
          list(
            method = method_name,
            scenario = scenario_name,
            profile = profile,
            seed = seed,
            runtime_sec = as.numeric(runtime_sec),
            error = error_msg
          ),
          scenario_cols,
          method_spec$meta,
          metrics
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  .combine_benchmark_rows(rows)
}

summarize_synthetic_alignment_benchmark <- function(results,
                                                    metric_cols = NULL,
                                                    group_cols = c("method", "scenario", "profile")) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame.", call. = FALSE)
  }
  if (!all(group_cols %in% names(results))) {
    stop("Missing grouping columns in `results`.", call. = FALSE)
  }

  if (is.null(metric_cols)) {
    metric_cols <- names(results)[vapply(results, is.numeric, logical(1))]
    metric_cols <- setdiff(metric_cols, "seed")
  }
  metric_cols <- intersect(metric_cols, names(results))
  out <- unique(results[group_cols])
  rownames(out) <- NULL
  out$.group_key <- .benchmark_group_key(out, group_cols)
  result_keys <- .benchmark_group_key(results, group_cols)

  n_runs <- stats::aggregate(
    rep(1L, nrow(results)),
    by = list(.group_key = result_keys),
    FUN = sum
  )
  out$n_runs <- as.integer(.as_numeric_column(n_runs$x[match(out$.group_key, n_runs$.group_key)]))

  if ("error" %in% names(results)) {
    n_success <- stats::aggregate(
      as.integer(is.na(results$error)),
      by = list(.group_key = result_keys),
      FUN = sum
    )
    out$n_success <- as.integer(.as_numeric_column(n_success$x[match(out$.group_key, n_success$.group_key)]))
    out$n_errors <- out$n_runs - out$n_success
    out$error_rate <- ifelse(out$n_runs > 0, out$n_errors / out$n_runs, NA_real_)
  } else {
    out$n_success <- out$n_runs
    out$n_errors <- 0L
    out$error_rate <- 0
  }

  if (!length(metric_cols)) {
    out$.group_key <- NULL
    return(out)
  }

  for (metric in metric_cols) {
    agg <- stats::aggregate(
      results[[metric]],
      by = list(.group_key = result_keys),
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
    )
    stats_mat <- .aggregate_stat_matrix(agg$x, c("mean", "sd"))
    idx <- match(out$.group_key, agg$.group_key)
    out[[paste0(metric, "_mean")]] <- stats_mat[idx, "mean"]
    out[[paste0(metric, "_sd")]] <- stats_mat[idx, "sd"]
  }

  out$.group_key <- NULL
  out
}

summarize_synthetic_alignment_seed_robustness <- function(results,
                                                          metric_cols = NULL,
                                                          group_cols = c("method", "scenario", "profile")) {
  if (!is.data.frame(results)) {
    stop("`results` must be a data.frame.", call. = FALSE)
  }
  if (!all(group_cols %in% names(results))) {
    stop("Missing grouping columns in `results`.", call. = FALSE)
  }

  if (is.null(metric_cols)) {
    metric_cols <- names(results)[vapply(results, is.numeric, logical(1))]
    metric_cols <- setdiff(metric_cols, "seed")
  }
  metric_cols <- intersect(metric_cols, names(results))

  out <- summarize_synthetic_alignment_benchmark(
    results = results,
    metric_cols = character(0),
    group_cols = group_cols
  )
  out$.group_key <- .benchmark_group_key(out, group_cols)
  result_keys <- .benchmark_group_key(results, group_cols)

  for (metric in metric_cols) {
    agg <- stats::aggregate(
      results[[metric]],
      by = list(.group_key = result_keys),
      FUN = function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) {
          return(c(
            mean = NA_real_,
            sd = NA_real_,
            median = NA_real_,
            iqr = NA_real_,
            min = NA_real_,
            max = NA_real_,
            cv = NA_real_
          ))
        }

        mean_x <- mean(x)
        sd_x <- stats::sd(x)
        c(
          mean = mean_x,
          sd = sd_x,
          median = stats::median(x),
          iqr = stats::IQR(x),
          min = min(x),
          max = max(x),
          cv = if (isTRUE(all.equal(mean_x, 0))) NA_real_ else sd_x / abs(mean_x)
        )
      }
    )
    stats_mat <- .aggregate_stat_matrix(
      agg$x,
      c("mean", "sd", "median", "iqr", "min", "max", "cv")
    )
    colnames(stats_mat) <- paste0(metric, "_", colnames(stats_mat))
    idx <- match(out$.group_key, agg$.group_key)
    for (j in seq_len(ncol(stats_mat))) {
      out[[colnames(stats_mat)[j]]] <- stats_mat[idx, j]
    }
  }

  out$.group_key <- NULL
  out
}
