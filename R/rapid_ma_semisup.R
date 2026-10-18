# RAPID-MA class-balanced semi-supervised teacher refinement --------------

#' Control RAPID-MA Semi-Supervised Teacher Refinement
#'
#' Controls deterministic relation propagation, class-balanced guarded
#' pseudo-label selection, EMA teacher updates, and rollback.
#'
#' @param rounds Maximum pseudo-label refresh rounds.
#' @param ema_decay Teacher retention in `[0, 1)` for the exponential moving
#'   average (EMA) teacher.
#' @param confidence_threshold Base class confidence threshold.
#' @param minimum_threshold Lower bound for minority-class thresholds.
#' @param balance_strength Amount by which rare classes receive lower thresholds.
#' @param max_pseudo_fraction Hard maximum fraction of eligible unlabelled
#'   nodes per domain that can be pseudo-labelled in one round. The budget is
#'   split evenly across classes (at least one candidate per class while budget
#'   remains); when the per-class floor would exceed the total budget, the most
#'   confident candidates are kept.
#' @param min_margin Minimum top-one versus top-two probability margin.
#' @param uncertainty_power Exponent applied to confidence-above-threshold weights.
#' @param min_weight Minimum accepted pseudo-label weight.
#' @param relation_consistency Whether candidates must agree with trusted
#'   relation neighbours.
#' @param relation_validation Whether explicitly supplied out-of-fold
#'   probabilities guard against misleading relation channels.
#' @param relation_validation_margin Minimum balanced-accuracy improvement over
#'   the full relation set before selecting a safer subset.
#' @param relation_validation_max_relations Maximum highest-weight relations
#'   used to construct bounded single and pair candidates.
#' @param propagation_steps Fixed number of sparse relation-propagation steps.
#' @param propagation_strength Convex weight assigned to relation neighbours.
#' @param min_relation_agreement Minimum neighbour probability for the candidate class.
#' @param relation_gate_min Ignore relations below this learned gate weight.
#' @param label_weight_multiplier Multiplier for pseudo-label transport costs.
#' @param consistency_weight,relation_weight Nonnegative refinement score weights.
#' @param rollback_tolerance Allowed score increase before rollback.
#' @return A `rapid_ma_semisup_control` object.
#' @export
rapid_ma_semisup_control <- function(
  rounds = 3L,
  ema_decay = 0.8,
  confidence_threshold = 0.8,
  minimum_threshold = 0.6,
  balance_strength = 0.08,
  max_pseudo_fraction = 0.5,
  min_margin = 0.15,
  uncertainty_power = 1,
  min_weight = 0.05,
  relation_consistency = TRUE,
  relation_validation = FALSE,
  relation_validation_margin = 1 / 6,
  relation_validation_max_relations = 6L,
  propagation_steps = 2L,
  propagation_strength = 0.35,
  min_relation_agreement = 0.5,
  relation_gate_min = 0.05,
  label_weight_multiplier = 1.5,
  consistency_weight = 0.5,
  relation_weight = 0.1,
  rollback_tolerance = 1e-6
) {
  rounds <- .rapid_int_scalar(rounds, "rounds", min_value = 0L)
  propagation_steps <- .rapid_int_scalar(
    propagation_steps, "propagation_steps", min_value = 0L
  )
  relation_validation_max_relations <- .rapid_int_scalar(
    relation_validation_max_relations,
    "relation_validation_max_relations"
  )
  values_01 <- c(
    ema_decay = ema_decay,
    confidence_threshold = confidence_threshold,
    minimum_threshold = minimum_threshold,
    max_pseudo_fraction = max_pseudo_fraction,
    min_margin = min_margin,
    min_weight = min_weight,
    propagation_strength = propagation_strength,
    min_relation_agreement = min_relation_agreement,
    relation_gate_min = relation_gate_min,
    relation_validation_margin = relation_validation_margin
  )
  for (name in names(values_01)) {
    values_01[[name]] <- .rapid_number_scalar(values_01[[name]], name, lower = 0)
    if (values_01[[name]] > 1) {
      stop("`", name, "` must be <= 1.", call. = FALSE)
    }
  }
  if (values_01[["ema_decay"]] >= 1) {
    stop("`ema_decay` must be < 1.", call. = FALSE)
  }
  if (values_01[["minimum_threshold"]] >
      values_01[["confidence_threshold"]]) {
    stop("`minimum_threshold` cannot exceed `confidence_threshold`.",
         call. = FALSE)
  }
  balance_strength <- .rapid_number_scalar(
    balance_strength, "balance_strength", lower = 0
  )
  uncertainty_power <- .rapid_number_scalar(
    uncertainty_power, "uncertainty_power", lower = 0
  )
  label_weight_multiplier <- .rapid_number_scalar(
    label_weight_multiplier, "label_weight_multiplier", lower = 0,
    strict = TRUE
  )
  consistency_weight <- .rapid_number_scalar(
    consistency_weight, "consistency_weight", lower = 0
  )
  relation_weight <- .rapid_number_scalar(
    relation_weight, "relation_weight", lower = 0
  )
  rollback_tolerance <- .rapid_number_scalar(
    rollback_tolerance, "rollback_tolerance", lower = 0
  )
  if (length(relation_consistency) != 1L || !is.logical(relation_consistency) ||
      is.na(relation_consistency)) {
    stop("`relation_consistency` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(relation_validation) != 1L || !is.logical(relation_validation) ||
      is.na(relation_validation)) {
    stop("`relation_validation` must be TRUE or FALSE.", call. = FALSE)
  }

  structure(
    list(
      rounds = rounds,
      ema_decay = values_01[["ema_decay"]],
      confidence_threshold = values_01[["confidence_threshold"]],
      minimum_threshold = values_01[["minimum_threshold"]],
      balance_strength = balance_strength,
      max_pseudo_fraction = values_01[["max_pseudo_fraction"]],
      min_margin = values_01[["min_margin"]],
      uncertainty_power = uncertainty_power,
      min_weight = values_01[["min_weight"]],
      relation_consistency = relation_consistency,
      relation_validation = relation_validation,
      relation_validation_margin =
        values_01[["relation_validation_margin"]],
      relation_validation_max_relations = relation_validation_max_relations,
      propagation_steps = propagation_steps,
      propagation_strength = values_01[["propagation_strength"]],
      min_relation_agreement = values_01[["min_relation_agreement"]],
      relation_gate_min = values_01[["relation_gate_min"]],
      label_weight_multiplier = label_weight_multiplier,
      consistency_weight = consistency_weight,
      relation_weight = relation_weight,
      rollback_tolerance = rollback_tolerance
    ),
    class = "rapid_ma_semisup_control"
  )
}

# Weight of a named relation, defaulting to 1 when the relation has no entry.
# `[[` on a named numeric vector errors for an absent name instead of
# returning NULL, so the name is matched explicitly.
.rapid_named_weight <- function(weights, name, default = 1) {
  if (is.null(weights) || is.null(names(weights)) ||
      !name %in% names(weights)) {
    return(default)
  }
  weights[[name]]
}

.rapid_propagate_probabilities <- function(
  probabilities,
  labels,
  relation_states,
  diffusion_fits,
  control
) {
  output <- .rapid_fix_true_probabilities(
    probabilities, labels, colnames(probabilities[[1L]])
  )
  if (control$propagation_steps == 0L || control$propagation_strength == 0) {
    return(output)
  }
  base <- output
  for (step in seq_len(control$propagation_steps)) {
    output <- lapply(seq_along(output), function(m) {
      probability <- output[[m]]
      numerator <- matrix(0, nrow(probability), ncol(probability))
      denominator <- numeric(nrow(probability))
      relations <- relation_states[[m]]$relations
      weights <- diffusion_fits[[m]]$relation_weights
      for (name in names(relations)) {
        weight <- .rapid_named_weight(weights, name)
        if (!is.finite(weight) || weight < control$relation_gate_min) next
        A <- relations[[name]]
        mass <- as.numeric(Matrix::rowSums(A))
        numerator <- numerator + weight * as.matrix(A %*% probability)
        denominator <- denominator + weight * mass
      }
      informed <- denominator > 0
      if (!any(informed)) return(probability)
      neighbour <- probability
      neighbour[informed, ] <-
        numerator[informed, , drop = FALSE] / denominator[informed]
      value <- (1 - control$propagation_strength) * base[[m]] +
        control$propagation_strength * neighbour
      value <- value / pmax(rowSums(value), 1e-15)
      value
    })
    output <- .rapid_fix_true_probabilities(
      output, labels, colnames(output[[1L]])
    )
  }
  names(output) <- names(probabilities)
  output
}

.rapid_resolve_semisup_control <- function(control = NULL) {
  .rapid_resolve_settings(
    control, rapid_ma_semisup_control, "rapid_ma_semisup_control", "semi-supervised",
    "`control` must be NULL, a named list, or rapid_ma_semisup_control()."
  )
}

.rapid_semisup_masks <- function(fit, labels, train_mask, mode) {
  sizes <- fit$domain_sizes
  domain_names <- fit$domain_names
  if (is.null(labels)) labels <- fit$labels
  labels <- .rapid_validate_domain_labels(labels, sizes, domain_names)
  if (is.null(train_mask)) {
    if (identical(mode, "inductive")) {
      stop("`train_mask` is required in inductive mode.", call. = FALSE)
    }
    train_mask <- lapply(sizes, function(n) rep(TRUE, n))
  }
  if (!is.list(train_mask) || length(train_mask) != length(sizes)) {
    stop("`train_mask` must provide one logical vector per domain.",
         call. = FALSE)
  }
  train_mask <- lapply(seq_along(train_mask), function(m) {
    value <- train_mask[[m]]
    if (length(value) != sizes[[m]] || !is.logical(value) || anyNA(value)) {
      stop("Domain `", domain_names[[m]],
           "`: train_mask must be a complete logical vector.", call. = FALSE)
    }
    value
  })

  if (identical(mode, "inductive")) {
    leaked_in_fit <- sum(vapply(seq_along(fit$labels), function(m) {
      sum(!train_mask[[m]] & !is.na(fit$labels[[m]]))
    }, integer(1)))
    if (leaked_in_fit > 0L) {
      stop(
        "Inductive refinement requires a core fit created without labels ",
        "outside `train_mask`.",
        call. = FALSE
      )
    }
  }
  ignored <- integer(length(labels))
  if (identical(mode, "inductive")) {
    for (m in seq_along(labels)) {
      ignored[[m]] <- sum(!train_mask[[m]] & !is.na(labels[[m]]))
      labels[[m]][!train_mask[[m]]] <- NA_character_
    }
  }
  eligible <- if (identical(mode, "transductive")) {
    lapply(sizes, function(n) rep(TRUE, n))
  } else {
    train_mask
  }
  list(
    labels = labels,
    train_mask = train_mask,
    eligible = eligible,
    ignored_labels = ignored
  )
}

.rapid_fix_true_probabilities <- function(probabilities, labels, class_levels) {
  lapply(seq_along(probabilities), function(m) {
    probability <- probabilities[[m]]
    observed <- which(!is.na(labels[[m]]))
    if (length(observed)) {
      class_index <- match(labels[[m]][observed], class_levels)
      valid <- !is.na(class_index)
      observed <- observed[valid]
      class_index <- class_index[valid]
      probability[observed, ] <- 0
      probability[cbind(observed, class_index)] <- 1
    }
    probability
  })
}

# Validate one class-probability matrix for domain `m` and expand it to the
# fit's class levels (absent classes get zero mass); rows are renormalized.
.rapid_semisup_probability_matrix <- function(fit, m, probability) {
  class_levels <- fit$prototypes$class_levels
  if (!is.numeric(probability) || length(dim(probability)) != 2L ||
      nrow(probability) != fit$domain_sizes[[m]]) {
    stop(
      "Domain `", fit$domain_names[[m]],
      "`: initial probabilities must be a numeric matrix with one row per node.",
      call. = FALSE
    )
  }
  probability <- as.matrix(probability)
  probability_classes <- colnames(probability)
  if (is.null(probability_classes) || any(!nzchar(probability_classes)) ||
      anyDuplicated(probability_classes)) {
    stop("Initial probability matrices need unique class column names.",
         call. = FALSE)
  }
  unknown <- setdiff(probability_classes, class_levels)
  if (length(unknown)) {
    stop(
      "Initial probabilities contain class(es) absent from the fit: ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (any(!is.finite(probability)) || any(probability < 0)) {
    stop("Initial probabilities must be finite and nonnegative.",
         call. = FALSE)
  }
  output <- matrix(
    0,
    nrow(probability),
    length(class_levels),
    dimnames = list(NULL, class_levels)
  )
  output[, probability_classes] <- probability
  mass <- rowSums(output)
  if (any(mass <= 0)) {
    stop("Every initial probability row must have positive mass.",
         call. = FALSE)
  }
  output / mass
}

.rapid_semisup_initial_probabilities <- function(fit, initial_probabilities) {
  if (is.null(initial_probabilities)) {
    return(list(
      probabilities = fit$prediction_probabilities,
      source = "native_transport"
    ))
  }
  if (!is.list(initial_probabilities) ||
      length(initial_probabilities) != length(fit$domain_sizes)) {
    stop(
      "`initial_probabilities` must provide one matrix per domain.",
      call. = FALSE
    )
  }
  probabilities <- lapply(seq_along(initial_probabilities), function(m) {
    .rapid_semisup_probability_matrix(fit, m, initial_probabilities[[m]])
  })
  names(probabilities) <- fit$domain_names
  list(probabilities = probabilities, source = "external")
}

# Normalize `validation_probabilities` to out-of-fold folds per domain. Each
# domain entry is either one probability matrix (a single fold that scores
# every observed label row; the caller guarantees that no row of it used any
# scored label) or a list of folds `list(rows = <row indices>, probabilities =
# <matrix>)`, where each fold's full matrix was computed without the labels
# of its `rows`. Only fold rows are scored, after propagating that fold's own
# matrix, so relation neighbours never carry a scored row's label.
.rapid_semisup_validation_folds <- function(fit, validation_probabilities) {
  if (!is.list(validation_probabilities) ||
      length(validation_probabilities) != length(fit$domain_sizes)) {
    stop(
      "`validation_probabilities` must provide one matrix or fold list per ",
      "domain.",
      call. = FALSE
    )
  }
  lapply(seq_along(validation_probabilities), function(m) {
    value <- validation_probabilities[[m]]
    n <- fit$domain_sizes[[m]]
    if (length(dim(value)) == 2L) {
      return(list(list(
        rows = seq_len(n),
        probabilities = .rapid_semisup_probability_matrix(fit, m, value)
      )))
    }
    if (!is.list(value) || !length(value)) {
      stop(
        "Domain `", fit$domain_names[[m]], "`: validation probabilities must ",
        "be a matrix or a non-empty list of folds.",
        call. = FALSE
      )
    }
    folds <- lapply(value, function(fold) {
      if (!is.list(fold) || is.null(fold$rows) || is.null(fold$probabilities)) {
        stop(
          "Domain `", fit$domain_names[[m]], "`: every validation fold needs ",
          "`rows` and `probabilities`.",
          call. = FALSE
        )
      }
      rows <- fold$rows
      if (!is.numeric(rows) || !length(rows) || anyNA(rows) ||
          any(abs(rows - round(rows)) > 1e-8) || any(rows < 1 | rows > n) ||
          anyDuplicated(rows)) {
        stop(
          "Domain `", fit$domain_names[[m]], "`: fold rows must be unique ",
          "row indices.",
          call. = FALSE
        )
      }
      list(
        rows = as.integer(rows),
        probabilities = .rapid_semisup_probability_matrix(
          fit, m, fold$probabilities
        )
      )
    })
    scored <- unlist(lapply(folds, `[[`, "rows"), use.names = FALSE)
    if (anyDuplicated(scored)) {
      stop(
        "Domain `", fit$domain_names[[m]], "`: validation folds must have ",
        "disjoint rows.",
        call. = FALSE
      )
    }
    folds
  })
}

.rapid_semisup_store_predictions <- function(fit, probabilities) {
  fit$prediction_probabilities <- probabilities
  fit$predictions <- lapply(probabilities, function(probability) {
    stats::setNames(
      colnames(probability)[max.col(probability, ties.method = "first")],
      rownames(probability)
    )
  })
  fit$prediction_confidence <- lapply(probabilities, function(probability) {
    apply(probability, 1L, max)
  })
  names(fit$prediction_probabilities) <- names(fit$predictions) <-
    names(fit$prediction_confidence) <- fit$domain_names
  fit
}

.rapid_semisup_relation_candidates <- function(names, weights, max_relations) {
  if (!length(names)) return(list(none = character(0), full = character(0)))
  weights <- weights[names]
  weights[!is.finite(weights)] <- 0
  ranked <- names[order(-weights, names, method = "radix")]
  ranked <- head(ranked, max_relations)
  candidates <- list(none = character(0), full = names)
  for (name in ranked) candidates[[paste0("only_", name)]] <- name
  if (length(ranked) >= 2L) {
    pairs <- utils::combn(ranked, 2L, simplify = FALSE)
    for (pair in pairs) {
      candidates[[paste0("pair_", paste(pair, collapse = "+"))]] <- pair
    }
  }
  key <- vapply(candidates, function(value) {
    paste(sort(unique(value), method = "radix"), collapse = "\r")
  }, character(1))
  candidates[!duplicated(key)]
}

.rapid_semisup_validation_score <- function(probability, labels, classes) {
  observed <- which(!is.na(labels) & labels %in% classes)
  if (!length(observed)) return(c(accuracy = NA_real_, log_loss = NA_real_))
  known <- probability[, classes, drop = FALSE]
  prediction <- classes[max.col(known, ties.method = "first")]
  recalls <- vapply(classes, function(class) {
    rows <- observed[labels[observed] == class]
    if (!length(rows)) return(NA_real_)
    mean(prediction[rows] == class)
  }, numeric(1))
  class_index <- match(labels[observed], classes)
  counts <- table(factor(labels[observed], levels = classes))
  weight <- 1 / pmax(as.numeric(counts[class_index]), 1)
  weight <- weight / mean(weight)
  c(
    accuracy = mean(recalls, na.rm = TRUE),
    log_loss = mean(
      weight * -log(pmax(known[cbind(observed, class_index)], 1e-15))
    )
  )
}

.rapid_semisup_validate_relations <- function(
  fit,
  validation_probabilities,
  labels,
  control
) {
  diagnostics <- vector("list", length(fit$domain_names))
  names(diagnostics) <- fit$domain_names
  if (!isTRUE(control$relation_validation)) {
    return(list(fit = fit, diagnostics = diagnostics, status = "disabled"))
  }
  if (is.null(validation_probabilities)) {
    stop(
      "`relation_validation = TRUE` requires out-of-fold ",
      "`validation_probabilities`.",
      call. = FALSE
    )
  }
  folds <- .rapid_semisup_validation_folds(fit, validation_probabilities)
  class_levels <- fit$prototypes$class_levels
  unknown <- fit$prototypes$control$unknown_level
  classes <- setdiff(class_levels, unknown)

  for (m in seq_along(fit$domain_names)) {
    relation_names <- names(fit$relations[[m]]$relations)
    weights <- fit$diffusion[[m]]$relation_weights
    candidates <- .rapid_semisup_relation_candidates(
      relation_names,
      weights,
      control$relation_validation_max_relations
    )
    score <- matrix(
      NA_real_, length(candidates), 2L,
      dimnames = list(names(candidates), c("accuracy", "log_loss"))
    )
    n <- fit$domain_sizes[[m]]
    scored <- rep(FALSE, n)
    for (fold in folds[[m]]) scored[fold$rows] <- TRUE
    scored_labels <- labels[[m]]
    scored_labels[!scored] <- NA_character_
    for (candidate in names(candidates)) {
      keep <- intersect(candidates[[candidate]], relation_names)
      relation_state <- fit$relations[[m]]
      diffusion_fit <- fit$diffusion[[m]]
      relation_state$relations <- relation_state$relations[keep]
      diffusion_fit$relation_weights <- diffusion_fit$relation_weights[keep]
      assembled <- matrix(
        0, n, length(class_levels), dimnames = list(NULL, class_levels)
      )
      for (fold in folds[[m]]) {
        propagated <- .rapid_propagate_probabilities(
          list(fold$probabilities),
          list(rep(NA_character_, n)),
          list(relation_state),
          list(diffusion_fit),
          control
        )[[1L]]
        assembled[fold$rows, ] <- propagated[fold$rows, , drop = FALSE]
      }
      score[candidate, ] <- .rapid_semisup_validation_score(
        assembled, scored_labels, classes
      )
    }
    full_index <- match("full", rownames(score))
    valid <- which(is.finite(score[, "accuracy"]) &
                     is.finite(score[, "log_loss"]))
    selected <- "full"
    switched <- FALSE
    if (length(valid) && !is.na(full_index) &&
        is.finite(score[full_index, "accuracy"])) {
      ordering <- valid[order(
        -score[valid, "accuracy"],
        score[valid, "log_loss"],
        rownames(score)[valid],
        method = "radix"
      )]
      best <- rownames(score)[ordering[[1L]]]
      improvement <- score[best, "accuracy"] - score[full_index, "accuracy"]
      if (!identical(best, "full") &&
          improvement + 1e-12 >= control$relation_validation_margin) {
        selected <- best
        switched <- TRUE
      }
    }
    keep <- intersect(candidates[[selected]], relation_names)
    fit$relations[[m]]$relations <- fit$relations[[m]]$relations[keep]
    fit$diffusion[[m]]$relation_weights <-
      fit$diffusion[[m]]$relation_weights[keep]
    diagnostics[[m]] <- list(
      selected = selected,
      selected_relations = keep,
      switched = switched,
      margin = control$relation_validation_margin,
      scores = as.data.frame(score),
      candidate_relations = candidates,
      validation_rows = sum(!is.na(scored_labels) & scored_labels %in% classes),
      validation_folds = length(folds[[m]])
    )
  }
  list(fit = fit, diagnostics = diagnostics, status = "validated")
}

.rapid_relation_agreement <- function(probability, choice, relation_state,
                                      diffusion_fit, control) {
  n <- nrow(probability)
  numerator <- denominator <- numeric(n)
  relations <- relation_state$relations
  weights <- diffusion_fit$relation_weights
  for (name in names(relations)) {
    weight <- .rapid_named_weight(weights, name)
    if (!is.finite(weight) || weight < control$relation_gate_min) next
    A <- relations[[name]]
    mass <- as.numeric(Matrix::rowSums(A))
    neighbour <- as.matrix(A %*% probability)
    numerator <- numerator + weight * neighbour[cbind(seq_len(n), choice)]
    denominator <- denominator + weight * mass
  }
  agreement <- rep(1, n)
  informed <- denominator > 0
  agreement[informed] <- numerator[informed] / denominator[informed]
  pmin(pmax(agreement, 0), 1)
}

.rapid_class_thresholds <- function(labels, classes, control) {
  observed <- unlist(labels, use.names = FALSE)
  counts <- table(factor(observed[!is.na(observed)], levels = classes))
  maximum <- max(as.numeric(counts), 1)
  scarcity <- log(pmax(as.numeric(counts), 1) / maximum)
  threshold <- control$confidence_threshold + control$balance_strength * scarcity
  threshold <- pmin(control$confidence_threshold,
                    pmax(control$minimum_threshold, threshold))
  stats::setNames(as.numeric(threshold), classes)
}

.rapid_select_pseudolabels <- function(probabilities, true_labels, eligible,
                                       relation_states, diffusion_fits,
                                       prototypes, control) {
  unknown <- prototypes$control$unknown_level
  classes <- setdiff(prototypes$class_levels, unknown)
  thresholds <- .rapid_class_thresholds(true_labels, classes, control)
  pseudo_labels <- lapply(true_labels, function(x) rep(NA_character_, length(x)))
  pseudo_weight <- lapply(true_labels, function(x) numeric(length(x)))
  selected <- relation_agreement <- vector("list", length(probabilities))

  for (m in seq_along(probabilities)) {
    known_probability <- probabilities[[m]][, classes, drop = FALSE]
    choice_local <- max.col(known_probability, ties.method = "first")
    choice <- match(classes[choice_local], prototypes$class_levels)
    confidence <- known_probability[cbind(seq_len(nrow(known_probability)), choice_local)]
    if (ncol(known_probability) > 1L) {
      second <- apply(known_probability, 1L, function(x) sort(x, decreasing = TRUE)[[2L]])
    } else {
      second <- rep(0, nrow(known_probability))
    }
    margin <- confidence - second
    entropy <- -rowSums(probabilities[[m]] * log(pmax(probabilities[[m]], 1e-15)))
    uncertainty <- 1 - entropy / max(log(ncol(probabilities[[m]])), 1)
    uncertainty <- pmin(pmax(uncertainty, 0), 1)
    agreement <- if (isTRUE(control$relation_consistency)) {
      .rapid_relation_agreement(
        probabilities[[m]], choice, relation_states[[m]],
        diffusion_fits[[m]], control
      )
    } else {
      rep(1, nrow(probabilities[[m]]))
    }
    relation_agreement[[m]] <- agreement
    chosen_class <- classes[choice_local]
    row_threshold <- thresholds[chosen_class]
    candidate <- eligible[[m]] & is.na(true_labels[[m]]) &
      confidence >= row_threshold & margin >= control$min_margin &
      agreement >= control$min_relation_agreement
    raw_weight <- ((confidence - row_threshold) /
      pmax(1 - row_threshold, 1e-12))^control$uncertainty_power
    raw_weight <- pmin(pmax(raw_weight, 0), 1) * uncertainty * agreement
    candidate <- candidate & raw_weight >= control$min_weight

    # `max_pseudo_fraction` is a hard cap on the rows selected per domain and
    # round. The budget is split evenly across classes with a floor of one
    # candidate per class; if that floor overshoots the budget, the most
    # confident candidates across classes are kept.
    budget <- floor(
      control$max_pseudo_fraction * sum(eligible[[m]] & is.na(true_labels[[m]]))
    )
    limit <- if (budget <= 0) {
      0L
    } else {
      max(1L, floor(budget / max(length(classes), 1L)))
    }
    keep <- integer(0)
    for (class in classes) {
      index <- which(candidate & chosen_class == class)
      if (!length(index)) next
      ordering <- order(-confidence[index], -margin[index], index, method = "radix")
      keep <- c(keep, head(index[ordering], limit))
    }
    if (length(keep) > budget) {
      ordering <- order(-confidence[keep], -margin[keep], keep, method = "radix")
      keep <- keep[ordering][seq_len(budget)]
    }
    selected[[m]] <- sort(keep)
    if (length(keep)) {
      pseudo_labels[[m]][keep] <- chosen_class[keep]
      pseudo_weight[[m]][keep] <- raw_weight[keep]
    }
  }
  list(
    labels = pseudo_labels,
    weight = pseudo_weight,
    selected = selected,
    thresholds = thresholds,
    relation_agreement = relation_agreement,
    counts = table(factor(
      unlist(pseudo_labels, use.names = FALSE), levels = classes
    ))
  )
}

.rapid_semisup_score <- function(probabilities, true_labels, selection, control) {
  true_loss <- numeric(0)
  pseudo_loss <- pseudo_weights <- relation_penalty <- numeric(0)
  for (m in seq_along(probabilities)) {
    observed <- which(!is.na(true_labels[[m]]))
    if (length(observed)) {
      index <- match(true_labels[[m]][observed], colnames(probabilities[[m]]))
      counts <- table(true_labels[[m]][observed])
      weight <- 1 / as.numeric(counts[true_labels[[m]][observed]])
      weight <- weight / mean(weight)
      true_loss <- c(true_loss, weight * -log(pmax(
        probabilities[[m]][cbind(observed, index)], 1e-15
      )))
    }
    selected <- selection$selected[[m]]
    if (length(selected)) {
      index <- match(selection$labels[[m]][selected], colnames(probabilities[[m]]))
      weight <- selection$weight[[m]][selected]
      pseudo_loss <- c(pseudo_loss, -log(pmax(
        probabilities[[m]][cbind(selected, index)], 1e-15
      )))
      pseudo_weights <- c(pseudo_weights, weight)
      relation_penalty <- c(
        relation_penalty, 1 - selection$relation_agreement[[m]][selected]
      )
    }
  }
  labeled <- if (length(true_loss)) mean(true_loss) else 0
  consistency <- if (length(pseudo_loss)) {
    sum(pseudo_weights * pseudo_loss) / sum(pseudo_weights)
  } else {
    0
  }
  relation <- if (length(relation_penalty)) mean(relation_penalty) else 0
  components <- c(labeled = labeled, consistency = consistency, relation = relation)
  list(
    score = labeled + control$consistency_weight * consistency +
      control$relation_weight * relation,
    components = components
  )
}

.rapid_semisup_transport <- function(fit, merged_labels, label_confidence, control) {
  uot_control <- unclass(fit$control$uot)
  uot_control$label_weight <- uot_control$label_weight *
    control$label_weight_multiplier
  transport <- vector("list", length(fit$scores))
  for (m in seq_along(fit$scores)) {
    transport[[m]] <- tryCatch(
      .rapid_uot_solve(
        fit$scores[[m]],
        fit$prototypes,
        structure = if (is.null(fit$preprocessing$structures)) {
          NULL
        } else {
          fit$preprocessing$structures[[m]]
        },
        attributes = if (is.null(fit$preprocessing$attributes)) {
          NULL
        } else {
          fit$preprocessing$attributes[[m]]
        },
        labels = merged_labels[[m]],
        label_confidence = label_confidence[[m]],
        control = uot_control,
        warm_start = fit$transport[[m]]
      ),
      error = function(e) {
        stop("Domain `", fit$domain_names[[m]], "` semi-supervised transport: ",
             conditionMessage(e), call. = FALSE)
      }
    )
  }
  names(transport) <- fit$domain_names
  transport
}

#' Refine a RAPID-MA Fit with Guarded Class-Balanced Pseudo-Labels
#'
#' Combines learned sparse relation propagation with class-balanced confidence,
#' margin, uncertainty, and neighbour-agreement guards. Accepted pseudo-labels
#' update sparse node-to-prototype transport; any score increase is rolled
#' back.
#'
#' @details
#' Observed labels are fixed throughout refinement. In inductive mode, held-out
#' label values are removed before thresholds, selection, or diagnostics are
#' computed, and pseudo-labels are restricted to training rows.
#'
#' With `control$relation_validation = TRUE`, every candidate relation subset
#' (none, all, each single relation, and each pair of the highest-weight
#' relations) is scored by propagating `validation_probabilities` and measuring
#' class-balanced accuracy on the observed labels. The full relation set is
#' replaced only when a subset improves on it by at least
#' `control$relation_validation_margin`. Because propagation mixes each row with
#' its relation neighbours, replacing only a labelled row's own probability by a
#' leave-one-out value is not out-of-fold: its neighbours' probabilities would
#' still have used its label. Folds avoid this: each fold's whole matrix is
#' propagated and only that fold's rows are scored, so per-row (leave-one-out)
#' or grouped refits give genuinely out-of-fold scores. Rows outside every fold
#' are not scored.
#'
#' @param fit A fitted `rapid_ma` object.
#' @param labels Optional training-label vectors. In inductive mode, labels
#'   outside `train_mask` are discarded before any statistic is computed.
#' @param train_mask Optional logical vectors identifying training rows.
#' @param mode Either `"transductive"` or `"inductive"`.
#' @param initial_probabilities Optional named-class probability matrices, one
#'   per domain. A supervised downstream classifier can seed the sparse relation
#'   teacher without providing additional labels.
#' @param validation_probabilities Optional out-of-fold probabilities, used
#'   only when relation validation is enabled
#'   (`control$relation_validation = TRUE`). One entry per domain: either a
#'   list of folds, each `list(rows = <row indices>, probabilities = <matrix>)`
#'   whose full probability matrix was computed without the labels of its
#'   `rows`, or a single probability matrix, which must then be out-of-fold for
#'   every observed label row and all of its relation neighbours (see
#'   Details).
#' @param control Settings returned by [rapid_ma_semisup_control()].
#' @return The updated `rapid_ma` object with refined probabilities and a
#'   `semisup` provenance record.
#' @export
rapid_ma_semisup <- function(
  fit,
  labels = NULL,
  train_mask = NULL,
  mode = c("transductive", "inductive"),
  initial_probabilities = NULL,
  validation_probabilities = NULL,
  control = rapid_ma_semisup_control()
) {
  if (!inherits(fit, "rapid_ma")) {
    stop("`fit` must be a fitted rapid_ma object.", call. = FALSE)
  }
  mode <- match.arg(mode)
  control <- .rapid_resolve_semisup_control(control)
  masks <- .rapid_semisup_masks(fit, labels, train_mask, mode)
  true_labels <- masks$labels
  unknown <- fit$prototypes$control$unknown_level
  known_classes <- setdiff(fit$prototypes$class_levels, unknown)
  if (!length(known_classes) || !any(vapply(
    true_labels, function(x) any(!is.na(x)), logical(1)
  ))) {
    fit$semisup <- list(
      status = "structural_only",
      mode = mode,
      history = data.frame(),
      pseudo_labels = lapply(true_labels, function(x) rep(NA_character_, length(x))),
      provenance = list(
        ignored_labels = masks$ignored_labels
      ),
      control = control
    )
    return(fit)
  }

  initialization <- .rapid_semisup_initial_probabilities(
    fit, initial_probabilities
  )
  relation_validation <- .rapid_semisup_validate_relations(
    fit, validation_probabilities, true_labels, control
  )
  fit <- relation_validation$fit
  probabilities <- .rapid_fix_true_probabilities(
    initialization$probabilities, true_labels, fit$prototypes$class_levels
  )
  probabilities <- .rapid_propagate_probabilities(
    probabilities, true_labels, fit$relations, fit$diffusion, control
  )
  fit <- .rapid_semisup_store_predictions(fit, probabilities)
  teacher <- probabilities
  base_objective <- fit$objective
  pseudo_labels <- lapply(true_labels, function(x) rep(NA_character_, length(x)))
  pseudo_weight <- lapply(true_labels, function(x) numeric(length(x)))
  history <- vector("list", control$rounds)
  accepted_rounds <- 0L
  status <- if (control$rounds == 0L) {
    if (control$propagation_steps > 0L && control$propagation_strength > 0) {
      "propagation_only"
    } else {
      "disabled"
    }
  } else {
    "max_rounds"
  }
  last_score <- NA_real_

  if (control$rounds > 0L) {
    for (round in seq_len(control$rounds)) {
      selection <- .rapid_select_pseudolabels(
        teacher, true_labels, masks$eligible, fit$relations, fit$diffusion,
        fit$prototypes, control
      )
      n_selected <- sum(lengths(selection$selected))
      if (!n_selected) {
        status <- "no_confident_pseudolabels"
        history[[round]] <- data.frame(
          round = round, selected = 0L, accepted = FALSE,
          score_before = NA_real_, score_after = NA_real_,
          labeled = NA_real_, consistency = NA_real_, relation = NA_real_,
          status = status, stringsAsFactors = FALSE
        )
        history <- history[seq_len(round)]
        break
      }
      merged <- confidence <- vector("list", length(true_labels))
      for (m in seq_along(true_labels)) {
        merged[[m]] <- true_labels[[m]]
        confidence[[m]] <- as.numeric(!is.na(true_labels[[m]]))
        selected <- selection$selected[[m]]
        if (length(selected)) {
          merged[[m]][selected] <- selection$labels[[m]][selected]
          confidence[[m]][selected] <- selection$weight[[m]][selected]
        }
      }
      before <- .rapid_semisup_score(probabilities, true_labels, selection, control)
      candidate_transport <- .rapid_semisup_transport(
        fit, merged, confidence, control
      )
      candidate_prediction <- .rapid_prediction_state(
        candidate_transport, fit$prototypes, merged
      )
      candidate_probabilities <- .rapid_fix_true_probabilities(
        candidate_prediction$probabilities,
        true_labels,
        fit$prototypes$class_levels
      )
      candidate_probabilities <- .rapid_propagate_probabilities(
        candidate_probabilities, true_labels, fit$relations, fit$diffusion,
        control
      )
      after <- .rapid_semisup_score(
        candidate_probabilities, true_labels, selection, control
      )
      accepted <- after$score <= before$score + control$rollback_tolerance
      history[[round]] <- data.frame(
        round = round,
        selected = as.integer(n_selected),
        accepted = accepted,
        score_before = before$score,
        score_after = after$score,
        labeled = after$components[["labeled"]],
        consistency = after$components[["consistency"]],
        relation = after$components[["relation"]],
        status = if (accepted) "accepted" else "rolled_back",
        stringsAsFactors = FALSE
      )
      if (!accepted) {
        status <- "rolled_back"
        history <- history[seq_len(round)]
        break
      }

      accepted_rounds <- accepted_rounds + 1L
      last_score <- after$score
      fit$transport <- candidate_transport
      fit$couplings <- lapply(candidate_transport, `[[`, "coupling")
      probabilities <- candidate_probabilities
      teacher <- Map(function(old, new) {
        value <- control$ema_decay * old + (1 - control$ema_decay) * new
        value / pmax(rowSums(value), 1e-15)
      }, teacher, candidate_probabilities)
      pseudo_labels <- selection$labels
      pseudo_weight <- selection$weight
      fit$dropped_mass <- vapply(
        candidate_transport,
        function(x) x$diagnostics$dropped_fraction,
        numeric(1)
      )
      fit <- .rapid_semisup_store_predictions(fit, probabilities)
      has_relations <- any(vapply(
        fit$relations, function(x) length(x$relations) > 0L, logical(1)
      ))
      if (has_relations) {
        fit$relation_match <- .rapid_relmatch_fit(
          fit$relations,
          fit$transport,
          domain_relation_weights = stats::setNames(
            lapply(fit$diffusion, `[[`, "relation_weights"),
            fit$domain_names
          ),
          control = fit$control$relmatch
        )
        fit$relation_barycenters <- fit$relation_match$barycenters
      }
      names(fit$couplings) <- names(fit$prediction_probabilities) <-
        names(fit$predictions) <- names(fit$prediction_confidence) <-
        fit$domain_names
    }
  }

  history <- if (length(history) && !all(vapply(history, is.null, logical(1)))) {
    do.call(rbind, history[!vapply(history, is.null, logical(1))])
  } else {
    data.frame()
  }
  fit$semisup <- list(
    status = status,
    mode = mode,
    accepted_rounds = accepted_rounds,
    base_objective = base_objective,
    score = last_score,
    history = history,
    pseudo_labels = pseudo_labels,
    pseudo_weight = pseudo_weight,
    teacher_probabilities = teacher,
    class_thresholds = .rapid_class_thresholds(true_labels, known_classes, control),
    provenance = list(
      ignored_labels = masks$ignored_labels,
      train_rows = vapply(masks$train_mask, sum, integer(1)),
      eligible_rows = vapply(masks$eligible, sum, integer(1)),
      true_label_counts = lapply(true_labels, function(x) table(x, useNA = "no")),
      probability_initialization = initialization$source,
      relation_validation = relation_validation$diagnostics
    ),
    control = control
  )
  fit
}
