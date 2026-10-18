# Reproducible RAPID-MA accuracy and package-comparison harness -------------

.rapid_benchmark_sample_classes <- function(n, probabilities, levels) {
  counts <- floor(n * probabilities)
  counts[[which.max(probabilities)]] <- counts[[which.max(probabilities)]] +
    n - sum(counts)
  rep(levels, counts)
}

#' @rdname benchmark_rapid_ma
#' @name benchmark_rapid_ma
#' @param n_common Number of corresponding nodes shared by the two domains.
#' @param extra_target Number of additional target-only nodes.
#' @param source_features,target_features Modality-specific feature widths.
#' @param labels_per_class Paired observed labels retained per class.
#' @export
rapid_ma_benchmark_fixture <- function(
  n_common = 120L,
  extra_target = 24L,
  source_features = 24L,
  target_features = 10L,
  labels_per_class = 4L,
  seed = 1L
) {
  values <- c(
    n_common = n_common,
    extra_target = extra_target,
    source_features = source_features,
    target_features = target_features,
    labels_per_class = labels_per_class,
    seed = seed
  )
  if (any(!is.numeric(values)) || any(!is.finite(values)) ||
      any(abs(values - round(values)) > 1e-8) ||
      n_common < 12L || extra_target < 0L || source_features < 2L ||
      target_features < 2L || labels_per_class < 1L || seed < 0L) {
    stop("Fixture sizes and seed must be valid nonnegative integers.",
         call. = FALSE)
  }
  n_common <- as.integer(n_common)
  extra_target <- as.integer(extra_target)
  source_features <- as.integer(source_features)
  target_features <- as.integer(target_features)
  labels_per_class <- as.integer(labels_per_class)
  seed <- as.integer(seed)
  .with_local_seed(seed, .rapid_benchmark_fixture_draw(
    n_common, extra_target, source_features, target_features,
    labels_per_class, seed
  ))
}

.rapid_benchmark_fixture_draw <- function(n_common, extra_target,
                                          source_features, target_features,
                                          labels_per_class, seed) {
  levels <- c("forest", "urban", "water")
  probabilities <- c(0.55, 0.30, 0.15)
  centers <- rbind(
    forest = c(-2.2, 0.2, 0.8),
    urban = c(2.0, 1.6, -0.4),
    water = c(0.3, -2.3, -1.1)
  )
  truth_common <- .rapid_benchmark_sample_classes(
    n_common, probabilities, levels
  )
  truth_extra <- if (extra_target) {
    sample(levels, extra_target, replace = TRUE, prob = probabilities)
  } else character(0)
  latent_common <- centers[truth_common, , drop = FALSE] +
    matrix(stats::rnorm(n_common * 3L, sd = 0.65), n_common)
  latent_extra <- if (extra_target) {
    centers[truth_extra, , drop = FALSE] +
      matrix(stats::rnorm(extra_target * 3L, sd = 0.75), extra_target)
  } else matrix(numeric(0), 0L, 3L)
  position_source <- latent_common[, 1:2, drop = FALSE] +
    matrix(stats::rnorm(n_common * 2L, sd = 0.12), n_common)
  position_target_unpermuted <- rbind(
    latent_common[, 1:2, drop = FALSE] +
      matrix(stats::rnorm(n_common * 2L, sd = 0.18), n_common),
    latent_extra[, 1:2, drop = FALSE] +
      matrix(stats::rnorm(extra_target * 2L, sd = 0.20), extra_target)
  )
  source_loading <- matrix(
    stats::rnorm(3L * source_features), 3L, source_features
  )
  target_loading <- matrix(
    stats::rnorm(3L * target_features), 3L, target_features
  )
  source <- latent_common %*% source_loading +
    matrix(stats::rnorm(n_common * source_features, sd = 0.28), n_common)
  target_latent <- rbind(latent_common, latent_extra)
  target <- tanh(target_latent %*% target_loading / 2) +
    matrix(
      stats::rnorm((n_common + extra_target) * target_features, sd = 0.34),
      n_common + extra_target
    )
  target_truth_unpermuted <- c(truth_common, truth_extra)
  permutation <- sample.int(n_common + extra_target)
  inverse_permutation <- order(permutation)
  target <- target[permutation, , drop = FALSE]
  position_target <- position_target_unpermuted[permutation, , drop = FALSE]
  target_truth <- target_truth_unpermuted[permutation]
  correspondence <- inverse_permutation[seq_len(n_common)]

  source_attributes <- data.frame(
    elevation = latent_common[, 3L] + stats::rnorm(n_common, sd = 0.15),
    texture = ifelse(latent_common[, 2L] > 0, "rough", "smooth"),
    stringsAsFactors = FALSE
  )
  target_attributes_unpermuted <- data.frame(
    elevation = target_latent[, 3L] +
      stats::rnorm(n_common + extra_target, sd = 0.25),
    texture = ifelse(target_latent[, 2L] > 0.2, "rough", "smooth"),
    stringsAsFactors = FALSE
  )
  target_attributes <- target_attributes_unpermuted[permutation, , drop = FALSE]
  rownames(target_attributes) <- NULL
  # Indexing the class-center matrix by truth labels names every row by its
  # class; drop those row names so held-out truth cannot travel with the data.
  rownames(source) <- rownames(target) <- NULL
  rownames(position_source) <- rownames(position_target) <- NULL
  missing_rows <- seq.int(7L, nrow(target_attributes), by = 19L)
  target_attributes$elevation[missing_rows] <- NA_real_

  labels_source <- rep(NA_character_, n_common)
  labels_target <- rep(NA_character_, n_common + extra_target)
  for (class in levels) {
    source_index <- head(which(truth_common == class), labels_per_class)
    if (length(source_index) < labels_per_class) {
      stop("`labels_per_class` exceeds a fixture class size.", call. = FALSE)
    }
    target_index <- correspondence[source_index]
    labels_source[source_index] <- class
    labels_target[target_index] <- class
  }
  list(
    data = list(source = source, target = target),
    truth = list(source = truth_common, target = target_truth),
    labels = list(source = labels_source, target = labels_target),
    positions = list(source = position_source, target = position_target),
    attributes = list(source = source_attributes, target = target_attributes),
    correspondence = correspondence,
    metadata = list(
      kind = "synthetic cross-modality land-cover fixture",
      license = "CC0 generated data; no external dataset redistributed",
      class_probabilities = stats::setNames(probabilities, levels),
      target_permutation = permutation,
      preprocessing = "method-local; observed labels fixed before fitting",
      seed = seed
    )
  )
}

.rapid_benchmark_probabilities <- function(probability, classes) {
  probability <- as.matrix(probability)
  output <- matrix(0, nrow(probability), length(classes),
                   dimnames = list(NULL, classes))
  shared <- intersect(colnames(probability), classes)
  if (length(shared)) output[, shared] <- probability[, shared, drop = FALSE]
  zero <- rowSums(output) <= 0
  if (any(zero)) output[zero, ] <- 1 / length(classes)
  output / pmax(rowSums(output), 1e-15)
}

# Held-out classification metrics shared by benchmark_rapid_ma() and the
# blocked-inductive evidence runner.
#
# Class universe for macro-F1 and mIoU: every class present in the held-out
# truth or among the held-out predictions. A prediction outside the truth
# classes (an open-set ".unknown", or a labelled class absent from the held-out
# rows) is a false positive of a class with zero F1 and IoU, and a miss for the
# true class. Missing predictions (NA) are misses only. OA and ECE (`bins`
# equal-width confidence bins on [0, 1]) use every held-out row; confidence is
# the largest row-normalized probability. `rare_class` defaults to the truth
# class with the smallest held-out support (ties broken alphabetically).
.rapid_classification_metrics <- function(truth, prediction, probability,
                                          rare_class = NULL, bins = 10L) {
  truth <- as.character(truth)
  prediction <- as.character(prediction)
  probability <- as.matrix(probability)
  if (!length(truth) || anyNA(truth) || length(prediction) != length(truth) ||
      nrow(probability) != length(truth)) {
    stop("Metrics need complete truth and one prediction and probability ",
         "row per held-out row.", call. = FALSE)
  }
  predicted <- !is.na(prediction)
  universe <- sort(unique(c(truth, prediction[predicted])), method = "radix")
  scores <- vapply(universe, function(class) {
    hit <- predicted & prediction == class
    tp <- sum(truth == class & hit)
    fp <- sum(truth != class & hit)
    fn <- sum(truth == class & !hit)
    precision <- if (tp + fp > 0) tp / (tp + fp) else 0
    recall <- if (tp + fn > 0) tp / (tp + fn) else 0
    f1 <- if (precision + recall > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
    iou <- if (tp + fp + fn > 0) tp / (tp + fp + fn) else 0
    c(recall = recall, f1 = f1, iou = iou)
  }, numeric(3))
  if (is.null(rare_class)) rare_class <- .rapid_benchmark_rare_class(truth)
  total <- rowSums(probability)
  confidence <- apply(probability, 1L, max) / ifelse(total > 0, total, 1)
  confidence[!is.finite(confidence) | total <= 0] <- 0
  correct <- as.numeric(predicted & prediction == truth)
  breaks <- seq(0, 1, length.out = bins + 1L)
  bucket <- pmin(findInterval(confidence, breaks, all.inside = TRUE), bins)
  ece <- sum(vapply(seq_len(bins), function(b) {
    index <- which(bucket == b)
    if (!length(index)) return(0)
    length(index) / length(bucket) *
      abs(mean(correct[index]) - mean(confidence[index]))
  }, numeric(1)))
  c(
    oa = mean(correct),
    macro_f1 = mean(scores["f1", ]),
    miou = mean(scores["iou", ]),
    rare_recall = if (rare_class %in% universe) {
      scores["recall", rare_class]
    } else {
      NA_real_
    },
    ece = ece
  )
}

# The evaluated class with the smallest held-out support, ties broken
# alphabetically. Computed from held-out truth at scoring time only.
.rapid_benchmark_rare_class <- function(truth) {
  support <- table(factor(truth, levels = sort(unique(truth), method = "radix")))
  names(support)[[which.min(support)]]
}

.rapid_benchmark_centroid <- function(embeddings, labels, classes) {
  train_x <- do.call(rbind, lapply(seq_along(embeddings), function(m) {
    embeddings[[m]][!is.na(labels[[m]]), , drop = FALSE]
  }))
  train_y <- unlist(lapply(labels, function(x) x[!is.na(x)]), use.names = FALSE)
  # A class without training rows has no centroid and receives zero
  # probability rather than a NaN centroid.
  present <- classes[classes %in% train_y]
  if (!length(present)) {
    stop("No class has a labelled training row.", call. = FALSE)
  }
  centroids <- do.call(rbind, lapply(present, function(class) {
    colMeans(train_x[train_y == class, , drop = FALSE])
  }))
  rownames(centroids) <- present
  predict_one <- function(X) {
    distance <- vapply(seq_len(nrow(centroids)), function(j) {
      rowSums((X - matrix(centroids[j, ], nrow(X), ncol(X), byrow = TRUE))^2)
    }, numeric(nrow(X)))
    distance <- as.matrix(distance)
    scale <- stats::median(distance[is.finite(distance) & distance > 0])
    if (!is.finite(scale) || scale <= 1e-12) scale <- 1
    probability <- exp(pmax(-distance / scale, -700))
    probability <- probability / pmax(rowSums(probability), 1e-15)
    colnames(probability) <- present
    output <- matrix(0, nrow(probability), length(classes),
                     dimnames = list(NULL, classes))
    output[, present] <- probability
    list(
      probability = output,
      prediction = present[max.col(probability, ties.method = "first")]
    )
  }
  lapply(embeddings, predict_one)
}

# Unsupervised domain-local embedding for the target-only comparator and the
# local part of the RAPID-MA teacher seed: the leading principal components
# (at most `ncomp`) of the standardized features. Every feature column
# contributes; no label is read. Constant columns standardize to zero.
.rapid_benchmark_local_embedding <- function(data, ncomp) {
  lapply(data, function(X) {
    Z <- scale(as.matrix(X), center = TRUE, scale = TRUE)
    attr(Z, "scaled:center") <- attr(Z, "scaled:scale") <- NULL
    Z[!is.finite(Z)] <- 0
    width <- max(1L, min(ncomp, ncol(Z), nrow(Z) - 1L))
    decomposition <- svd(Z, nu = 0L, nv = width)
    scores <- Z %*% decomposition$v[, seq_len(width), drop = FALSE]
    colnames(scores) <- paste0("pc_", seq_len(width))
    scores
  })
}

.rapid_benchmark_teacher_seeds <- function(
  fit,
  data,
  labels,
  classes,
  ncomp,
  local_weight
) {
  local_embedding <- .rapid_benchmark_local_embedding(data, ncomp)
  blend <- function(m, domain_labels) {
    shared <- .rapid_benchmark_centroid(
      fit$scores, domain_labels, classes
    )[[m]]$probability
    local_available <- all(classes %in% domain_labels[[m]])
    if (!local_available) return(shared)
    local <- .rapid_benchmark_centroid(
      list(local_embedding[[m]]), list(domain_labels[[m]]), classes
    )[[1L]]$probability
    local_weight * local + (1 - local_weight) * shared
  }
  initial <- lapply(seq_along(labels), blend, domain_labels = labels)
  # Leave-one-label-out folds: each fold recomputes the whole probability
  # matrix of its domain without the held-out label, so the relation-safety
  # gate propagates probabilities in which no neighbour used that label.
  # A row whose removal leaves a class without any label is not scored.
  validation <- lapply(seq_along(labels), function(m) {
    folds <- list()
    for (row in which(!is.na(labels[[m]]))) {
      held_out_labels <- labels
      held_out_labels[[m]][row] <- NA_character_
      remaining <- unlist(held_out_labels, use.names = FALSE)
      if (!all(classes %in% remaining)) next
      folds[[length(folds) + 1L]] <- list(
        rows = row,
        probabilities = blend(m, held_out_labels)
      )
    }
    if (!length(folds)) {
      stop("The relation-safety gate has no out-of-fold label row in ",
           "domain ", m, ".", call. = FALSE)
    }
    folds
  })
  list(initial = initial, validation = validation)
}

.rapid_benchmark_correspondence <- function(correspondence, n_source, n_target) {
  if (is.null(correspondence)) return(rep(NA_integer_, n_source))
  if (is.matrix(correspondence) || is.data.frame(correspondence)) {
    if (ncol(correspondence) != 2L) {
      stop("`correspondence` must have two columns.", call. = FALSE)
    }
    output <- rep(NA_integer_, n_source)
    output[as.integer(correspondence[, 1L])] <- as.integer(correspondence[, 2L])
  } else {
    if (length(correspondence) != n_source) {
      stop("Correspondence vectors need one entry per source row.", call. = FALSE)
    }
    output <- as.integer(correspondence)
  }
  bad <- !is.na(output) & (output < 1L | output > n_target)
  if (any(bad)) stop("Correspondence indices are outside target bounds.", call. = FALSE)
  output
}

.rapid_benchmark_matching <- function(source, target, correspondence,
                                      rank_cap = 100L) {
  known <- which(!is.na(correspondence))
  if (!length(known) || is.null(source) || is.null(target)) {
    return(c(hits1 = NA_real_, mrr = NA_real_, coverage = NA_real_))
  }
  rank_cap <- .rapid_int_scalar(rank_cap, "rank_cap")
  source_finite <- rowSums(!is.finite(source[known, , drop = FALSE])) == 0L
  target_finite <- rowSums(!is.finite(target)) == 0L
  valid_target <- which(target_finite)
  truth_local <- match(correspondence[known], valid_target)
  valid <- source_finite & !is.na(truth_local)
  if (!any(valid) || !length(valid_target)) {
    return(c(hits1 = NA_real_, mrr = NA_real_, coverage = mean(valid)))
  }
  query_rows <- known[valid]
  k <- min(rank_cap, length(valid_target))
  neighbors <- RANN::nn2(
    data = target[valid_target, , drop = FALSE],
    query = source[query_rows, , drop = FALSE],
    k = k,
    treetype = "kd",
    searchtype = "standard",
    eps = 0
  )$nn.idx
  if (is.null(dim(neighbors))) neighbors <- matrix(neighbors, ncol = k)
  local_truth <- truth_local[valid]
  ranks <- vapply(seq_len(nrow(neighbors)), function(row) {
    rank <- match(local_truth[[row]], neighbors[row, ])
    if (is.na(rank)) Inf else rank
  }, numeric(1))
  c(
    hits1 = mean(ranks == 1),
    mrr = mean(ifelse(is.finite(ranks), 1 / ranks, 0)),
    coverage = mean(valid)
  )
}

.rapid_benchmark_split_scores <- function(scores, sizes) {
  offsets <- c(0L, cumsum(sizes))
  lapply(seq_along(sizes), function(m) {
    scores[(offsets[[m]] + 1L):offsets[[m + 1L]], , drop = FALSE]
  })
}

.rapid_benchmark_rss <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  as.numeric(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^2
}

.rapid_benchmark_measure <- function(code) {
  gc(FALSE)
  before <- .rapid_benchmark_rss()
  start <- proc.time()[["elapsed"]]
  value <- force(code)
  elapsed <- proc.time()[["elapsed"]] - start
  after <- .rapid_benchmark_rss()
  list(
    value = value,
    runtime_sec = as.numeric(elapsed),
    retained_mb = as.numeric(utils::object.size(value)) / 1024^2,
    endpoint_rss_mb = if (all(is.na(c(before, after)))) NA_real_ else
      max(c(before, after), na.rm = TRUE),
    rss_delta_mb = if (is.finite(before) && is.finite(after)) after - before else NA_real_
  )
}

.rapid_benchmark_lema_path <- function() {
  installed <- system.file("benchmarks", "lema_reference.R", package = "manifoldalign")
  candidates <- c(installed, file.path("inst", "benchmarks", "lema_reference.R"))
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("The non-production LeMA reference file is unavailable.", call. = FALSE)
  }
  candidates[[1L]]
}

.rapid_benchmark_lema <- function(data, labels, correspondence, ncomp, seed,
                                  control = list()) {
  environment <- new.env(parent = baseenv())
  sys.source(.rapid_benchmark_lema_path(), envir = environment)
  pairs <- which(!is.na(labels[[1L]]) & !is.na(correspondence) &
                   !is.na(labels[[2L]][correspondence]))
  pairs <- pairs[labels[[1L]][pairs] == labels[[2L]][correspondence[pairs]]]
  if (!length(pairs)) {
    stop("LeMA needs at least one completed labeled cross-modality pair.",
         call. = FALSE)
  }
  paired_target <- correspondence[pairs]
  unlabeled_target <- setdiff(seq_len(nrow(data[[2L]])), paired_target)
  arguments <- utils::modifyList(
    list(
      X_high = data[[1L]][pairs, , drop = FALSE],
      X_low_labeled = data[[2L]][paired_target, , drop = FALSE],
      X_low_unlabeled = data[[2L]][unlabeled_target, , drop = FALSE],
      labels = labels[[1L]][pairs],
      ncomp = ncomp,
      seed = seed
    ),
    control
  )
  fit <- do.call(environment$lema_reference_fit, arguments)
  target_order <- c(paired_target, unlabeled_target)
  probability <- fit$probabilities[
    length(pairs) + seq_along(target_order), , drop = FALSE
  ]
  target_probability <- probability[match(seq_len(nrow(data[[2L]])), target_order),
                                    , drop = FALSE]
  target_embedding <- fit$target_embedding[
    match(seq_len(nrow(data[[2L]])), target_order), , drop = FALSE
  ]
  source_embedding <- matrix(
    NA_real_, nrow(data[[1L]]), ncol(fit$source_embedding)
  )
  source_embedding[pairs, ] <- fit$source_embedding
  embeddings <- list(source_embedding, target_embedding)
  classifier <- .rapid_benchmark_centroid(
    embeddings, labels, colnames(target_probability)
  )[[2L]]
  list(
    fit = fit,
    embeddings = embeddings,
    probability = classifier$probability,
    prediction = classifier$prediction,
    native_probability = target_probability,
    evaluation_mode = "transductive",
    classifier = "shared nearest centroid",
    pairing = "completed labeled correspondences"
  )
}

#' Reproducible RAPID-MA Accuracy and Package Comparison
#'
#' Builds a deterministic imbalanced cross-modality land-cover fixture and
#' compares RAPID-MA with LeMA and related package methods on identical
#' splits. Embedding methods use the same nearest-centroid prediction head.
#' Native UOT and sparse relation teacher rows are named explicitly so their
#' accuracy is not conflated with embedding quality.
#'
#' @details
#' Held-out truth is used only after fitting, for scoring. The class set given
#' to every prediction head is the set of classes with at least one observed
#' (training) label; held-out truth never defines it. The rare class reported
#' by `rare_recall` is the evaluated class with the smallest held-out support
#' (ties broken alphabetically) and is computed only at scoring time.
#'
#' Macro-F1 and mIoU follow one class-universe rule, shared with the
#' blocked-inductive evidence runner: they average over every class present in
#' the held-out truth or in the held-out predictions. Predictions outside the
#' truth classes, including open-set `".unknown"` predictions and labelled
#' classes absent from the held-out rows, therefore count as false positives
#' of a class whose F1 and IoU are zero, as well as misses for the true class.
#' Overall accuracy and expected calibration error (ten equal-width bins) use
#' all held-out rows.
#'
#' The `target_centroid` comparator and the domain-local part of the RAPID-MA
#' teacher seed both use a principal-component projection of the standardized
#' domain features (at most `ncomp` components, computed without labels); all
#' feature columns contribute.
#'
#' The `rapid_ma_semisup` row's relation-safety gate scores relation subsets on
#' genuinely out-of-fold teacher probabilities: for each observed label row,
#' the whole probability matrix of its domain is recomputed without that label
#' before propagation, so neighbours never carry the held-out label. Rows whose
#' removal would leave a class without any label are not scored.
#'
#' The bundled LeMA implementation is an independent dense equation-level
#' benchmark oracle and is not available as a production backend. Reported
#' endpoint RSS is a lower bound rather than a sampled peak; evidence scripts
#' should measure process peak RSS separately. Hits@1 uses exact
#' nearest-neighbour retrieval. MRR is a bounded MRR@100 lower bound; true
#' matches outside the first 100 candidates contribute zero.
#'
#' @param fixture A fixture from `rapid_ma_benchmark_fixture()`, or a
#'   compatible two-domain list containing `data`, `truth`, `labels`, optional
#'   `positions` and `attributes`, and optional source-to-target
#'   `correspondence`.
#' @param methods Character vector of supported comparison methods:
#'   `target_centroid`, `rapid_ma`, `rapid_ma_native`, `rapid_ma_semisup`,
#'   `lema`, `ssma`, `kema`, `multiscale`, `spectral_mnn`, and `token_ot`.
#'   Embedding methods use one shared nearest-centroid head; explicitly named
#'   RAPID native and semi-supervised rows retain their own prediction heads.
#' @param ncomp Shared requested embedding width.
#' @param seed Deterministic fixture or method seed. Package functions restore
#'   the caller's random-number stream on exit.
#' @param controls Optional named per-method control lists.
#' @param keep_fits Whether to retain fitted objects.
#' @return `rapid_ma_benchmark_fixture()` returns generated data, held-out
#'   truth, observed labels, positions, attributes, correspondences, and
#'   provenance. `benchmark_rapid_ma()` returns a `rapid_ma_benchmark` list with
#'   same-split accuracy, calibration, matching, runtime, memory lower bounds,
#'   failures, split provenance, and optionally fits.
#' @export
benchmark_rapid_ma <- function(
  fixture = rapid_ma_benchmark_fixture(),
  methods = c(
    "target_centroid", "rapid_ma", "rapid_ma_native",
    "rapid_ma_semisup", "lema",
    "ssma", "kema", "multiscale", "spectral_mnn", "token_ot"
  ),
  ncomp = 8L,
  seed = 1L,
  controls = list(),
  keep_fits = FALSE
) {
  supported <- c(
    "target_centroid", "rapid_ma", "rapid_ma_native",
    "rapid_ma_semisup", "lema",
    "ssma", "kema", "multiscale", "spectral_mnn", "token_ot"
  )
  unknown <- setdiff(methods, supported)
  if (length(unknown)) {
    stop("Unknown benchmark method(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  if (anyDuplicated(methods)) stop("`methods` cannot contain duplicates.", call. = FALSE)
  data <- fixture$data
  truth <- fixture$truth
  labels <- fixture$labels
  if (!is.list(data) || length(data) != 2L || !is.list(truth) ||
      length(truth) != 2L || !is.list(labels) || length(labels) != 2L) {
    stop("The benchmark currently requires two data, truth, and label domains.",
         call. = FALSE)
  }
  data <- lapply(data, as.matrix)
  sizes <- vapply(data, nrow, integer(1))
  if (any(lengths(truth) != sizes) || any(lengths(labels) != sizes)) {
    stop("Truth and label vectors must match their domain row counts.",
         call. = FALSE)
  }
  domain_names <- names(data)
  if (is.null(domain_names) || any(!nzchar(domain_names))) {
    domain_names <- c("source", "target")
    names(data) <- names(truth) <- names(labels) <- domain_names
  }
  # Prediction heads only know classes with an observed (training) label;
  # held-out truth never defines the class set.
  classes <- sort(unique(stats::na.omit(unlist(labels, use.names = FALSE))),
                  method = "radix")
  if (!length(classes)) stop("No observed labels are available.", call. = FALSE)
  target_test <- is.na(labels[[2L]])
  if (!any(target_test)) stop("The target domain has no held-out labels.", call. = FALSE)
  correspondence <- .rapid_benchmark_correspondence(
    fixture$correspondence, sizes[[1L]], sizes[[2L]]
  )
  positions <- if (is.null(fixture$positions)) list(NULL, NULL) else fixture$positions
  attributes <- if (is.null(fixture$attributes)) list(NULL, NULL) else fixture$attributes
  hd <- as_hyperdesign(data, labels = labels, label_name = "benchmark_label",
                       domain_names = domain_names)
  rank <- min(32L, min(sizes) - 1L)
  knn <- min(10L, min(sizes) - 1L)
  fits <- list()
  rows <- vector("list", length(methods))

  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    method_control <- if (!is.null(controls[[method]])) controls[[method]] else list()
    attempt <- tryCatch({
      measured <- switch(
        method,
        target_centroid = .rapid_benchmark_measure({
          embedding <- .rapid_benchmark_local_embedding(data, ncomp)
          prediction <- .rapid_benchmark_centroid(
            list(embedding[[2L]]), list(labels[[2L]]), classes
          )[[1L]]
          list(
            embeddings = embedding,
            probability = prediction$probability,
            prediction = prediction$prediction,
            fit = NULL,
            evaluation_mode = "inductive target-only",
            classifier = "target-only nearest centroid",
            matching_eligible = FALSE
          )
        }),
        rapid_ma = .rapid_benchmark_measure({
          control <- do.call(
            rapid_ma_control,
            utils::modifyList(
              list(max_iter = 4L, min_iter = 1L, seed = seed),
              method_control
            )
          )
          fit <- rapid_ma(
            data, labels = labels, positions = positions,
            attributes = attributes, ncomp = ncomp, control = control
          )
          prediction <- .rapid_benchmark_centroid(
            fit$scores, labels, classes
          )[[2L]]
          list(
            embeddings = fit$scores,
            probability = prediction$probability,
            prediction = prediction$prediction,
            fit = fit,
            evaluation_mode = "transductive",
            classifier = "shared nearest centroid",
            matching_eligible = TRUE
          )
        }),
        rapid_ma_native = .rapid_benchmark_measure({
          control <- do.call(
            rapid_ma_control,
            utils::modifyList(
              list(max_iter = 4L, min_iter = 1L, seed = seed),
              method_control
            )
          )
          fit <- rapid_ma(
            data, labels = labels, positions = positions,
            attributes = attributes, ncomp = ncomp, control = control
          )
          list(
            embeddings = fit$scores,
            probability = fit$prediction_probabilities[[2L]],
            prediction = fit$predictions[[2L]],
            fit = fit,
            evaluation_mode = "transductive open-set",
            classifier = "native sparse UOT",
            matching_eligible = TRUE
          )
        }),
        rapid_ma_semisup = .rapid_benchmark_measure({
          core_control <- if (!is.null(method_control$core)) method_control$core else list()
          teacher_control <- if (!is.null(method_control$teacher)) method_control$teacher else list()
          local_weight <- if (!is.null(method_control$seed_local_weight)) {
            method_control$seed_local_weight
          } else {
            0.75
          }
          if (!is.numeric(local_weight) || length(local_weight) != 1L ||
              !is.finite(local_weight) || local_weight < 0 || local_weight > 1) {
            stop("`seed_local_weight` must be one number in [0, 1].",
                 call. = FALSE)
          }
          fit <- rapid_ma(
            data, labels = labels, positions = positions,
            attributes = attributes, ncomp = ncomp,
            control = utils::modifyList(
              list(max_iter = 4L, min_iter = 1L, seed = seed),
              core_control
            )
          )
          seeds <- .rapid_benchmark_teacher_seeds(
            fit, data, labels, classes, ncomp, local_weight
          )
          fit <- rapid_ma_semisup(
            fit,
            initial_probabilities = seeds$initial,
            validation_probabilities = seeds$validation,
            control = utils::modifyList(
              list(
                rounds = 0L,
                propagation_steps = 2L,
                propagation_strength = 1,
                relation_validation = TRUE,
                relation_validation_margin = 1 / 6
              ),
              teacher_control
            )
          )
          list(
            embeddings = fit$scores,
            probability = fit$prediction_probabilities[[2L]],
            prediction = fit$predictions[[2L]],
            fit = fit,
            evaluation_mode = "transductive",
            classifier = paste0(
              "blended centroid (", format(local_weight, trim = TRUE),
              " local) + validated sparse relation teacher"
            ),
            matching_eligible = TRUE
          )
        }),
        lema = .rapid_benchmark_measure(
          .rapid_benchmark_lema(
            data, labels, correspondence, ncomp, seed, method_control
          )
        ),
        ssma = .rapid_benchmark_measure({
          fit <- ssma_align(
            hd, y = benchmark_label, ncomp = ncomp,
            control = do.call(
              ssma_align_control,
              utils::modifyList(
                list(knn = knn, rank_per_domain = rank,
                     solver = "reduced", seed = seed, verbose = FALSE),
                method_control
              )
            )
          )
          embedding <- .rapid_benchmark_split_scores(fit$s, sizes)
          prediction <- .rapid_benchmark_centroid(embedding, labels, classes)[[2L]]
          list(embeddings = embedding, probability = prediction$probability,
               prediction = prediction$prediction, fit = fit,
               evaluation_mode = "transductive",
               classifier = "shared nearest centroid",
               matching_eligible = TRUE)
        }),
        kema = .rapid_benchmark_measure({
          arguments <- utils::modifyList(
            list(data = hd, ncomp = ncomp, knn = min(5L, knn),
                 solver = "regression", sample_frac = 1,
                 lambda = 0.01),
            method_control
          )
          fit <- rlang::inject(kema(!!!arguments, y = benchmark_label))
          embedding <- .rapid_benchmark_split_scores(as.matrix(fit$s), sizes)
          prediction <- .rapid_benchmark_centroid(embedding, labels, classes)[[2L]]
          list(embeddings = embedding, probability = prediction$probability,
               prediction = prediction$prediction, fit = fit,
               evaluation_mode = "transductive",
               classifier = "shared nearest centroid",
               matching_eligible = TRUE)
        }),
        multiscale = .rapid_benchmark_measure({
          fit <- multiscale_manifold_align(
            hd, y = benchmark_label, ncomp = ncomp, verbose = FALSE,
            control = do.call(
              multiscale_manifold_align_control,
              utils::modifyList(
                list(knn = knn, rank_per_domain = rank, max_levels = 4L,
                     tune = FALSE, seed = seed, verbose = FALSE),
                method_control
              )
            )
          )
          embedding <- .rapid_benchmark_split_scores(as.matrix(fit$s), sizes)
          prediction <- .rapid_benchmark_centroid(embedding, labels, classes)[[2L]]
          list(embeddings = embedding, probability = prediction$probability,
               prediction = prediction$prediction, fit = fit,
               evaluation_mode = "transductive",
               classifier = "shared nearest centroid",
               matching_eligible = TRUE)
        }),
        spectral_mnn = .rapid_benchmark_measure({
          fit <- spectral_mnn_align(
            hd, y = benchmark_label, ncomp = ncomp, verbose = FALSE,
            control = do.call(
              spectral_mnn_align_control,
              utils::modifyList(
                list(knn = knn, k_embed = max(2L * ncomp, ncomp + 4L),
                     seed = seed, verbose = FALSE),
                method_control
              )
            )
          )
          embedding <- .rapid_benchmark_split_scores(as.matrix(fit$s), sizes)
          prediction <- .rapid_benchmark_centroid(embedding, labels, classes)[[2L]]
          list(embeddings = embedding, probability = prediction$probability,
               prediction = prediction$prediction, fit = fit,
               evaluation_mode = "transductive",
               classifier = "shared nearest centroid",
               matching_eligible = TRUE)
        }),
        token_ot = .rapid_benchmark_measure({
          fit <- token_ot_graph_align(
            hd, ncomp = ncomp, verbose = FALSE,
            control = do.call(
              token_ot_graph_align_control,
              utils::modifyList(
                list(graph_knn = knn, candidate_k = min(32L, sizes[[2L]]),
                     views = c("eigenmap", "hks"), token_mode = "view_only",
                     seed = seed, verbose = FALSE),
                method_control
              )
            )
          )
          embedding <- .rapid_benchmark_split_scores(as.matrix(fit$s), sizes)
          prediction <- .rapid_benchmark_centroid(embedding, labels, classes)[[2L]]
          list(embeddings = embedding, probability = prediction$probability,
               prediction = prediction$prediction, fit = fit,
               evaluation_mode = "transductive",
               classifier = "shared nearest centroid",
               matching_eligible = TRUE)
        })
      )
      value <- measured$value
      classification <- .rapid_classification_metrics(
        truth[[2L]][target_test], value$prediction[target_test],
        value$probability[target_test, , drop = FALSE]
      )
      matching <- if (isFALSE(value$matching_eligible)) {
        c(hits1 = NA_real_, mrr = NA_real_, coverage = NA_real_)
      } else {
        .rapid_benchmark_matching(
          value$embeddings[[1L]], value$embeddings[[2L]], correspondence
        )
      }
      rows[[method_index]] <- data.frame(
        method = method,
        status = "ok",
        evaluation_mode = value$evaluation_mode,
        classifier = value$classifier,
        runtime_sec = measured$runtime_sec,
        retained_mb = measured$retained_mb,
        endpoint_rss_mb = measured$endpoint_rss_mb,
        rss_delta_mb = measured$rss_delta_mb,
        oa = classification[["oa"]],
        macro_f1 = classification[["macro_f1"]],
        miou = classification[["miou"]],
        rare_recall = classification[["rare_recall"]],
        ece = classification[["ece"]],
        hits1 = matching[["hits1"]],
        mrr = matching[["mrr"]],
        coverage = matching[["coverage"]],
        error = NA_character_,
        stringsAsFactors = FALSE
      )
      fits[[method]] <- value$fit
      TRUE
    }, error = function(error) {
      rows[[method_index]] <<- data.frame(
        method = method,
        status = "failed",
        evaluation_mode = NA_character_,
        classifier = NA_character_,
        runtime_sec = NA_real_, retained_mb = NA_real_,
        endpoint_rss_mb = NA_real_, rss_delta_mb = NA_real_,
        oa = NA_real_, macro_f1 = NA_real_, miou = NA_real_,
        rare_recall = NA_real_, ece = NA_real_, hits1 = NA_real_,
        mrr = NA_real_, coverage = NA_real_,
        error = conditionMessage(error), stringsAsFactors = FALSE
      )
      FALSE
    })
    invisible(attempt)
  }
  results <- do.call(rbind, rows)
  rownames(results) <- NULL
  rare_class <- .rapid_benchmark_rare_class(truth[[2L]][target_test])
  structure(
    list(
      results = results,
      fits = if (isTRUE(keep_fits)) fits else NULL,
      split = list(
        target_test_rows = which(target_test),
        target_train_rows = which(!target_test),
        fit_label_counts = vapply(
          labels, function(x) sum(!is.na(x)), integer(1)
        ),
        correspondence = correspondence,
        rare_class = rare_class,
        seed = as.integer(seed)
      ),
      provenance = list(
        fixture = fixture$metadata,
        class_levels = classes,
        ncomp = as.integer(ncomp),
        requested_methods = methods,
        memory = paste(
          "retained_mb is serialized R object size; endpoint_rss_mb is an",
          "endpoint RSS lower bound, not a sampled peak"
        ),
        matching = paste(
          "Hits@1 is exact nearest-neighbour retrieval; mrr is a bounded",
          "MRR@100 lower bound with ranks beyond 100 contributing zero"
        ),
        lema_reference = "Hong et al. 2019 equation-level dense oracle",
        class_levels_source = "observed (training) labels",
        metrics = paste(
          "macro-F1 and mIoU average over classes present in held-out truth",
          "or held-out predictions; out-of-truth predictions are false",
          "positives; rare class = smallest held-out support"
        ),
        local_embedding = paste(
          "leading principal components (at most ncomp) of the standardized",
          "domain features; no labels"
        ),
        rapid_teacher_seed = paste(
          "fixed 0.75 domain-local and 0.25 shared-embedding centroid blend;",
          "two sparse relation-propagation hops with a leave-one-label-out",
          "safety gate (each fold recomputes the domain's probabilities)"
        )
      )
    ),
    class = "rapid_ma_benchmark"
  )
}
