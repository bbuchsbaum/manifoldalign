.bench_or <- function(x, y) if (is.null(x)) y else x

.rand_rotation_benchmark <- function(d, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Q <- qr.Q(qr(matrix(rnorm(d * d), d)))
  if (det(Q) < 0) Q[, 1] <- -Q[, 1]
  Q
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

  list(
    latent = latent,
    view1 = X1,
    view2 = X2,
    view3 = X3,
    true_labels = true_labels,
    class_centers = list(c(-0.5, -0.5), c(0.5, 0.5))
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

  list(
    latent = t,
    view1 = v1,
    view2 = v2,
    view3 = v3,
    true_labels = true_labels,
    parameter_boundary = pi
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

  list(
    latent = grid,
    view1 = v1,
    view2 = v2,
    view3 = v3,
    hard_indices = idx_hard,
    true_labels = true_labels
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
      difficulty = "hard",
      notes = "Grid latent space with nonlinear distortion, density skew, and heteroscedastic corruption."
    )
  )

  fast_specs <- list(
    linear_affine = utils::modifyList(full_specs$linear_affine, list(defaults = list(n = 20L, noise = 0.01, seed = 1L))),
    isometric_curve = utils::modifyList(full_specs$isometric_curve, list(defaults = list(n = 25L, noise = 0.02, seed = 2L))),
    hard_nonisometric = utils::modifyList(full_specs$hard_nonisometric, list(defaults = list(n_side = 5L, noise = 0.03, seed = 3L)))
  )

  if (identical(profile, "full")) full_specs else fast_specs
}

synthetic_alignment_scenarios <- function(profile = c("full", "fast")) {
  profile <- match.arg(profile)
  specs <- .synthetic_benchmark_specs(profile)

  rows <- lapply(names(specs), function(name) {
    spec <- specs[[name]]
    defaults <- spec$defaults
    n_samples <- if (!is.null(defaults$n_side)) defaults$n_side^2 else defaults$n

    data.frame(
      scenario = name,
      profile = profile,
      family = spec$family,
      latent_relation = spec$latent_relation,
      supervision = spec$supervision,
      side_information = spec$side_information,
      difficulty = spec$difficulty,
      n_views = 3L,
      n_samples = as.integer(n_samples),
      default_seed = as.integer(defaults$seed),
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

  label_vec <- .bench_or(label, scenario$true_labels)
  domains <- lapply(views, function(view_name) {
    X <- as.matrix(scenario[[view_name]])
    design_df <- data.frame(sample_id = seq_len(nrow(X)), stringsAsFactors = FALSE)
    names(design_df)[1] <- id_col
    if (!is.null(label_vec)) {
      design_df[[label_col]] <- label_vec
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

run_synthetic_alignment_benchmark <- function(methods,
                                              scenarios = NULL,
                                              profile = c("full", "fast"),
                                              seeds = 1L,
                                              score_fn = NULL,
                                              verbose = FALSE) {
  profile <- match.arg(profile)

  if (is.function(methods)) {
    methods <- list(method = methods)
  }
  if (!is.list(methods) || !length(methods)) {
    stop("`methods` must be a non-empty named list of functions.", call. = FALSE)
  }
  if (is.null(names(methods)) || any(!nzchar(names(methods)))) {
    stop("`methods` must be named.", call. = FALSE)
  }
  if (!all(vapply(methods, is.function, logical(1)))) {
    stop("All entries in `methods` must be functions.", call. = FALSE)
  }

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
      for (method_name in names(methods)) {
        if (isTRUE(verbose)) {
          message(
            "run_synthetic_alignment_benchmark: scenario=", scenario_name,
            " seed=", seed,
            " method=", method_name
          )
        }

        fit <- NULL
        error_msg <- NA_character_
        runtime_sec <- system.time({
          fit <- tryCatch(
            methods[[method_name]](scenario),
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
  if (!length(metric_cols)) {
    base <- unique(results[group_cols])
    rownames(base) <- NULL
    return(base)
  }

  out <- unique(results[group_cols])
  rownames(out) <- NULL

  for (metric in metric_cols) {
    agg <- stats::aggregate(
      results[[metric]],
      by = results[group_cols],
      FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
    )
    out[[paste0(metric, "_mean")]] <- agg$x[, "mean"]
    out[[paste0(metric, "_sd")]] <- agg$x[, "sd"]
  }

  out
}
