make_rapid_semisup_domain <- function(n, p, seed, noise = 0.9) {
  set.seed(seed)
  truth <- rep(c("a", "b", "c"), length.out = n)
  centers <- rbind(
    a = c(-1.5, 0, 0),
    b = c(1.5, 0, 0),
    c = c(0, 1.8, 0)
  )
  latent <- centers[truth, , drop = FALSE] +
    matrix(stats::rnorm(n * 3L, sd = noise), n)
  features <- latent %*% matrix(stats::rnorm(3L * p), 3L, p) +
    matrix(stats::rnorm(n * p, sd = 0.8), n)
  positions <- latent[, 1:2, drop = FALSE] +
    matrix(stats::rnorm(n * 2L, sd = 0.25), n)
  attributes <- data.frame(
    u = latent[, 1L] + stats::rnorm(n, sd = 0.4),
    v = latent[, 2L] + stats::rnorm(n, sd = 0.4)
  )
  list(X = features, truth = truth, positions = positions, attributes = attributes)
}

sparse_rapid_semisup_labels <- function(truth, per_class = 2L) {
  labels <- rep(NA_character_, length(truth))
  for (class in unique(truth)) {
    index <- which(truth == class)
    labels[sample(index, per_class)] <- class
  }
  labels
}

fit_rapid_semisup_fixture <- function(seed = 8L, n = c(120L, 135L)) {
  d1 <- make_rapid_semisup_domain(n[[1L]], 18L, 1000L + seed)
  d2 <- make_rapid_semisup_domain(n[[2L]], 12L, 2000L + seed)
  set.seed(3000L + seed)
  labels <- list(
    a = sparse_rapid_semisup_labels(d1$truth),
    b = sparse_rapid_semisup_labels(d2$truth)
  )
  fit <- rapid_ma(
    list(a = d1$X, b = d2$X),
    labels = labels,
    positions = list(d1$positions, d2$positions),
    attributes = list(d1$attributes, d2$attributes),
    ncomp = 8L,
    control = list(
      max_iter = 2L,
      min_iter = 1L,
      uot = list(q = 8L, epsilon = 0.18, label_weight = 3)
    )
  )
  list(fit = fit, domains = list(a = d1, b = d2), labels = labels)
}

test_that("semi-supervised control validates bounded refinement settings", {
  control <- rapid_ma_semisup_control(
    rounds = 4L,
    propagation_steps = 3L,
    propagation_strength = 0.4,
    max_pseudo_fraction = 0.25
  )

  expect_s3_class(control, "rapid_ma_semisup_control")
  expect_identical(control$rounds, 4L)
  expect_identical(control$propagation_steps, 3L)
  expect_equal(control$propagation_strength, 0.4)
  expect_error(rapid_ma_semisup_control(ema_decay = 1), "must be < 1")
  expect_error(
    rapid_ma_semisup_control(minimum_threshold = 0.9, confidence_threshold = 0.8),
    "cannot exceed"
  )
  expect_error(
    rapid_ma_semisup_control(propagation_strength = 1.1),
    "must be <= 1"
  )
  expect_error(
    manifoldalign:::.rapid_resolve_semisup_control(list(not_a_field = 1)),
    "Unknown semi-supervised"
  )
  expect_error(
    rapid_ma_semisup_control(relation_validation_margin = 1.1),
    "must be <= 1"
  )
})

test_that("out-of-fold safety validation removes a misleading relation", {
  labels <- list(c("a", "a", "a", "b", "b", "b"))
  probability <- matrix(
    c(
      0.4, 0.6, 0, 0.9, 0.1, 0, 0.9, 0.1, 0,
      0.6, 0.4, 0, 0.1, 0.9, 0, 0.1, 0.9, 0
    ),
    6L,
    3L,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", ".unknown"))
  )
  good <- Matrix::sparseMatrix(
    i = 1:6,
    j = c(2L, 3L, 2L, 5L, 6L, 5L),
    x = 1,
    dims = c(6L, 6L)
  )
  bad <- Matrix::sparseMatrix(
    i = 1:6,
    j = c(4L, 5L, 6L, 1L, 2L, 3L),
    x = 1,
    dims = c(6L, 6L)
  )
  fit <- list(
    domain_names = "d1",
    domain_sizes = c(d1 = 6L),
    prototypes = list(
      class_levels = colnames(probability),
      control = list(unknown_level = ".unknown")
    ),
    relations = list(list(relations = list(good = good, bad = bad))),
    diffusion = list(list(relation_weights = c(good = 0.5, bad = 0.5))),
    prediction_probabilities = list(probability)
  )
  control <- rapid_ma_semisup_control(
    rounds = 0L,
    propagation_steps = 1L,
    propagation_strength = 1,
    relation_gate_min = 0,
    relation_validation = TRUE,
    relation_validation_margin = 0.2
  )

  selected <- manifoldalign:::.rapid_semisup_validate_relations(
    fit, list(probability), labels, control
  )

  expect_identical(
    selected$diagnostics$d1$selected_relations,
    "good"
  )
  expect_true(selected$diagnostics$d1$switched)
  expect_identical(names(selected$fit$relations[[1L]]$relations), "good")
  expect_error(
    manifoldalign:::.rapid_semisup_validate_relations(
      fit, NULL, labels, control
    ),
    "requires out-of-fold"
  )
})

test_that("sparse relation propagation corrects locally inconsistent probabilities", {
  adjacency <- Matrix::bdiag(
    matrix(1, 4L, 4L) - diag(4L),
    matrix(1, 4L, 4L) - diag(4L)
  )
  probability <- rbind(
    c(0.95, 0.04, 0.01), c(0.85, 0.10, 0.05),
    c(0.35, 0.60, 0.05), c(0.80, 0.15, 0.05),
    c(0.04, 0.95, 0.01), c(0.55, 0.40, 0.05),
    c(0.10, 0.85, 0.05), c(0.15, 0.80, 0.05)
  )
  colnames(probability) <- c("a", "b", ".unknown")
  labels <- list(c("a", rep(NA_character_, 3L), "b", rep(NA_character_, 3L)))
  relation_states <- list(list(relations = list(spatial = adjacency)))
  diffusion_fits <- list(list(relation_weights = c(spatial = 1)))
  control <- rapid_ma_semisup_control(
    propagation_steps = 3L,
    propagation_strength = 0.7,
    relation_gate_min = 0
  )

  propagated <- manifoldalign:::.rapid_propagate_probabilities(
    list(probability), labels, relation_states, diffusion_fits, control
  )[[1L]]

  expect_gt(propagated[3L, "a"], propagated[3L, "b"])
  expect_gt(propagated[6L, "b"], propagated[6L, "a"])
  expect_identical(unname(propagated[1L, ]), c(1, 0, 0))
  expect_identical(unname(propagated[5L, ]), c(0, 1, 0))
  expect_equal(rowSums(propagated), rep(1, 8L), tolerance = 1e-12)
})

test_that("external supervised probabilities can seed propagation without new labels", {
  fixture <- fit_rapid_semisup_fixture(seed = 10L, n = c(42L, 45L))
  fit <- fixture$fit
  classes <- setdiff(
    fit$prototypes$class_levels,
    fit$prototypes$control$unknown_level
  )
  initial <- lapply(seq_along(fit$domain_sizes), function(m) {
    probability <- matrix(
      1 / length(classes),
      fit$domain_sizes[[m]],
      length(classes),
      dimnames = list(NULL, classes)
    )
    observed <- which(!is.na(fixture$labels[[m]]))
    probability[observed, ] <- 0
    probability[cbind(
      observed,
      match(fixture$labels[[m]][observed], classes)
    )] <- 1
    probability
  })

  propagated <- rapid_ma_semisup(
    fit,
    initial_probabilities = initial,
    control = list(
      rounds = 0L,
      propagation_steps = 2L,
      propagation_strength = 1
    )
  )

  expect_identical(propagated$semisup$status, "propagation_only")
  expect_identical(
    propagated$semisup$provenance$probability_initialization,
    "external"
  )
  expect_false(propagated$semisup$provenance$test_labels_read)
  expect_equal(
    unlist(lapply(propagated$prediction_probabilities, rowSums),
           use.names = FALSE),
    rep(1, sum(fit$domain_sizes)),
    tolerance = 1e-12
  )
  expect_error(
    rapid_ma_semisup(fit, initial_probabilities = list(initial[[1L]])),
    "one matrix per domain"
  )
})

test_that("pseudo-label selection is guarded, balanced, and can be disabled", {
  adjacency <- Matrix::bdiag(
    matrix(1, 4L, 4L) - diag(4L),
    matrix(1, 4L, 4L) - diag(4L)
  )
  probability <- rbind(
    c(1, 0, 0), c(0.90, 0.08, 0.02), c(0.85, 0.10, 0.05), c(0.6, 0.35, 0.05),
    c(0, 1, 0), c(0.08, 0.90, 0.02), c(0.10, 0.85, 0.05), c(0.35, 0.6, 0.05)
  )
  colnames(probability) <- c("a", "b", ".unknown")
  labels <- list(c("a", rep(NA_character_, 3L), "b", rep(NA_character_, 3L)))
  relation_states <- list(list(relations = list(spatial = adjacency)))
  diffusion_fits <- list(list(relation_weights = c(spatial = 1)))
  prototypes <- list(
    class_levels = colnames(probability),
    control = list(unknown_level = ".unknown")
  )
  control <- rapid_ma_semisup_control(
    confidence_threshold = 0.65,
    minimum_threshold = 0.55,
    max_pseudo_fraction = 0.5,
    min_margin = 0.1,
    min_weight = 0,
    min_relation_agreement = 0.5,
    relation_gate_min = 0
  )
  selection <- manifoldalign:::.rapid_select_pseudolabels(
    list(probability), labels, list(rep(TRUE, 8L)), relation_states,
    diffusion_fits, prototypes, control
  )

  expect_lte(length(selection$selected[[1L]]), 2L)
  expect_false(any(selection$selected[[1L]] %in% c(1L, 5L)))
  expect_setequal(stats::na.omit(selection$labels[[1L]]), c("a", "b"))

  disabled <- manifoldalign:::.rapid_select_pseudolabels(
    list(probability), labels, list(rep(TRUE, 8L)), relation_states,
    diffusion_fits, prototypes,
    rapid_ma_semisup_control(max_pseudo_fraction = 0, min_weight = 0)
  )
  expect_length(disabled$selected[[1L]], 0L)

  thresholds <- manifoldalign:::.rapid_class_thresholds(
    list(c(rep("a", 8L), "b")), c("a", "b"), control
  )
  expect_lt(thresholds[["b"]], thresholds[["a"]])
})

test_that("guarded refinement improves a hard held-out structured fixture", {
  fixture <- fit_rapid_semisup_fixture(seed = 8L)
  fit <- fixture$fit
  truth <- unlist(lapply(fixture$domains, `[[`, "truth"), use.names = FALSE)
  observed <- !is.na(unlist(fixture$labels, use.names = FALSE))
  accuracy_before <- mean(unlist(fit$predictions, use.names = FALSE)[!observed] ==
                            truth[!observed])

  refined <- rapid_ma_semisup(
    fit,
    control = list(
      rounds = 3L,
      confidence_threshold = 0.7,
      minimum_threshold = 0.55,
      min_margin = 0.05,
      min_relation_agreement = 0.2,
      min_weight = 0.005,
      label_weight_multiplier = 3,
      propagation_strength = 0.45,
      propagation_steps = 3L
    )
  )
  accuracy_after <- mean(
    unlist(refined$predictions, use.names = FALSE)[!observed] == truth[!observed]
  )

  expect_gt(accuracy_after, accuracy_before)
  expect_gte(accuracy_after - accuracy_before, 0.03)
  expect_gt(refined$semisup$accepted_rounds, 0L)
  expect_true(all(
    refined$semisup$history$score_after[refined$semisup$history$accepted] <=
      refined$semisup$history$score_before[refined$semisup$history$accepted] +
      refined$semisup$control$rollback_tolerance
  ))
  for (m in seq_along(fixture$labels)) {
    index <- which(!is.na(fixture$labels[[m]]))
    expect_identical(refined$predictions[[m]][index], fixture$labels[[m]][index])
  }

  repeated <- rapid_ma_semisup(fit, control = refined$semisup$control)
  expect_identical(repeated$predictions, refined$predictions)
  expect_equal(repeated$semisup$history, refined$semisup$history, tolerance = 1e-12)
})

test_that("inductive masks make held-out label values observationally irrelevant", {
  fixture <- fit_rapid_semisup_fixture(seed = 4L, n = c(60L, 66L))
  fit <- fixture$fit
  train_mask <- lapply(seq_along(fixture$labels), function(m) {
    mask <- seq_along(fixture$labels[[m]]) <= floor(0.65 * length(fixture$labels[[m]]))
    mask[!is.na(fixture$labels[[m]])] <- TRUE
    mask
  })
  labels_one <- labels_two <- fixture$labels
  for (m in seq_along(labels_one)) {
    held_out <- !train_mask[[m]]
    labels_one[[m]][held_out] <- fixture$domains[[m]]$truth[held_out]
    labels_two[[m]][held_out] <- rev(fixture$domains[[m]]$truth[held_out])
  }
  control <- list(
    rounds = 2L,
    confidence_threshold = 0.65,
    minimum_threshold = 0.55,
    min_margin = 0.05,
    min_relation_agreement = 0.2,
    min_weight = 0.001
  )

  first <- rapid_ma_semisup(
    fit, labels = labels_one, train_mask = train_mask,
    mode = "inductive", control = control
  )
  second <- rapid_ma_semisup(
    fit, labels = labels_two, train_mask = train_mask,
    mode = "inductive", control = control
  )

  expect_false(first$semisup$provenance$test_labels_read)
  expect_identical(
    unname(first$semisup$provenance$ignored_labels),
    vapply(train_mask, function(x) sum(!x), integer(1))
  )
  expect_identical(second$predictions, first$predictions)
  expect_equal(second$prediction_probabilities, first$prediction_probabilities,
               tolerance = 1e-12)
  for (m in seq_along(train_mask)) {
    expect_true(all(is.na(first$semisup$pseudo_labels[[m]][!train_mask[[m]]])))
  }

  leaked <- fit
  leaked$labels[[1L]][which(!train_mask[[1L]])[[1L]]] <- "a"
  expect_error(
    rapid_ma_semisup(leaked, train_mask = train_mask, mode = "inductive"),
    "created without labels outside"
  )
})

test_that("fully unlabeled fits remain explicitly structural-only", {
  d1 <- make_rapid_semisup_domain(30L, 6L, 5001L)
  d2 <- make_rapid_semisup_domain(33L, 7L, 5002L)
  fit <- rapid_ma(
    list(a = d1$X, b = d2$X),
    positions = list(d1$positions, d2$positions),
    ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L)
  )
  refined <- rapid_ma_semisup(fit)

  expect_identical(refined$semisup$status, "structural_only")
  expect_identical(refined$predictions, fit$predictions)
  expect_false(refined$semisup$provenance$test_labels_read)
})
