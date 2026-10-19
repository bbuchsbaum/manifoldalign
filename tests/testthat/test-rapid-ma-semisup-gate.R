# Relation-safety gate: no-transfer reference, comparator, and leakage.

gate_stub_fit <- function(relations, weights, n) {
  list(
    domain_names = "d1",
    domain_sizes = c(d1 = n),
    prototypes = list(
      class_levels = c("a", "b", ".unknown"),
      control = list(unknown_level = ".unknown")
    ),
    relations = list(list(relations = relations)),
    diffusion = list(list(relation_weights = weights))
  )
}

gate_stub_probability <- function(wrong_rows, n_per_class = 6L) {
  n <- 2L * n_per_class
  probability <- matrix(0, n, 3L, dimnames = list(NULL, c("a", "b", ".unknown")))
  probability[seq_len(n_per_class), c("a", "b")] <- rep(c(0.9, 0.1), each = n_per_class)
  probability[n_per_class + seq_len(n_per_class), c("a", "b")] <-
    rep(c(0.1, 0.9), each = n_per_class)
  for (row in wrong_rows) {
    probability[row, c("a", "b")] <- rev(probability[row, c("a", "b")]) * 0.8 + 0.1
  }
  probability
}

gate_stub_control <- function(margin = 1 / 6) {
  rapid_ma_semisup_control(
    rounds = 0L,
    propagation_steps = 1L,
    propagation_strength = 1,
    relation_gate_min = 0,
    relation_validation = TRUE,
    relation_validation_margin = margin
  )
}

# Six labels per class; a misleading relation sends row 2 (class a) to row 7
# (class b). Every other row is its own neighbour.
gate_stub_misleading <- function() {
  labels <- list(rep(c("a", "b"), each = 6L))
  bad <- Matrix::sparseMatrix(
    i = 1:12, j = c(1L, 7L, 3:12), x = 1, dims = c(12L, 12L)
  )
  list(
    labels = labels,
    fit = gate_stub_fit(list(bad = bad), c(bad = 1), 12L)
  )
}

misleading_benchmark_fixture <- function(seed) {
  fixture <- rapid_ma_benchmark_fixture(
    n_common = 48L, extra_target = 12L, labels_per_class = 3L, seed = seed
  )
  noise <- manifoldalign:::.with_local_seed(
    900000L + seed, stats::rnorm(length(fixture$data[[2L]]))
  )
  fixture$data[[2L]] <- fixture$data[[2L]] +
    matrix(0.75 * noise, nrow(fixture$data[[2L]]))
  permutation <- manifoldalign:::.with_local_seed(
    800000L + seed, sample.int(nrow(fixture$data[[2L]]))
  )
  fixture$positions[[2L]] <- fixture$positions[[2L]][permutation, , drop = FALSE]
  fixture$attributes[[2L]] <- fixture$attributes[[2L]][permutation, , drop = FALSE]
  fixture
}

test_that("harmful structure is compared with the no-relation reference", {
  # Regression: the full relation set used to be the comparator, so a
  # harmful full set was kept unless a subset beat it by the margin (1/6).
  # Here dropping the relation is better by one held-out row (1/12) only.
  stub <- gate_stub_misleading()
  probability <- gate_stub_probability(wrong_rows = 1L)
  gate <- manifoldalign:::.rapid_semisup_validate_relations(
    stub$fit, list(probability), stub$labels, gate_stub_control()
  )
  diagnostics <- gate$diagnostics$d1
  expect_equal(diagnostics$scores["none", "accuracy"], 11 / 12)
  expect_equal(diagnostics$scores["full", "accuracy"], 10 / 12)
  expect_identical(diagnostics$reference, "none")
  expect_identical(diagnostics$comparator, "none")
  expect_identical(diagnostics$selected, "none")
  expect_identical(diagnostics$decision, "rejected_transfer")
  expect_true(diagnostics$switched)
  expect_length(gate$fit$relations[[1L]]$relations, 0L)
  expect_length(gate$fit$diffusion[[1L]]$relation_weights, 0L)
})

test_that("the gate selects the no-transfer fallback when structure misleads", {
  stub <- gate_stub_misleading()
  # The transferred seed misses rows 1 and 3; the domain-local fallback
  # misses only row 1.
  seed <- gate_stub_probability(wrong_rows = c(1L, 3L))
  fallback <- gate_stub_probability(wrong_rows = 1L)
  folds <- list(list(rows = 1:12, probabilities = seed, fallback = fallback))
  gate <- manifoldalign:::.rapid_semisup_validate_relations(
    stub$fit, list(folds), stub$labels, gate_stub_control(),
    fallback_probabilities = list(fallback)
  )
  diagnostics <- gate$diagnostics$d1
  scores <- diagnostics$scores
  expect_identical(diagnostics$reference, "fallback")
  expect_identical(diagnostics$selected, "fallback")
  expect_identical(diagnostics$teacher_seed, "fallback")
  expect_true(gate$fallback_selected[["d1"]])
  expect_length(diagnostics$selected_relations, 0L)
  # Out of fold, the selected teacher is never below the fallback.
  expect_gte(
    scores[diagnostics$selected, "accuracy"], scores["fallback", "accuracy"]
  )
  expect_lt(scores["full", "accuracy"], scores["fallback", "accuracy"])
  expect_lt(scores["none", "accuracy"], scores["fallback", "accuracy"])

  # Informative structure: a relation that repairs row 1 from its class
  # neighbour makes the full set better than the fallback, so it is kept.
  good <- Matrix::sparseMatrix(
    i = 1:12, j = c(4L, 2:12), x = 1, dims = c(12L, 12L)
  )
  informative <- gate_stub_fit(list(good = good), c(good = 1), 12L)
  kept <- manifoldalign:::.rapid_semisup_validate_relations(
    informative, list(folds), stub$labels, gate_stub_control(),
    fallback_probabilities = list(fallback)
  )
  expect_identical(kept$diagnostics$d1$selected, "full")
  expect_identical(kept$diagnostics$d1$decision, "kept_full")
  expect_false(kept$fallback_selected[["d1"]])
  expect_identical(names(kept$fit$relations[[1L]]$relations), "good")
})

test_that("the candidate selector never ends below the no-transfer reference", {
  select <- manifoldalign:::.rapid_semisup_select_candidate
  score <- function(accuracy, log_loss, names) {
    matrix(c(accuracy, log_loss), ncol = 2L,
           dimnames = list(names, c("accuracy", "log_loss")))
  }
  names <- c("none", "full", "only_x", "fallback")
  # Full not below the reference: kept unless an alternative beats it by
  # the margin without a worse log-loss.
  choice <- select(score(c(0.5, 0.75, 0.8, 0.75), c(1, 1, 0.9, 1), names),
                   "fallback", 1 / 6)
  expect_identical(choice$selected, "full")
  expect_identical(choice$decision, "kept_full")
  choice <- select(score(c(0.5, 0.75, 0.95, 0.75), c(1, 1, 0.9, 1), names),
                   "fallback", 1 / 6)
  expect_identical(choice$selected, "only_x")
  # A larger accuracy bought with a worse log-loss is not selected.
  choice <- select(score(c(0.5, 0.75, 0.95, 0.75), c(1, 1, 1.2, 1), names),
                   "fallback", 1 / 6)
  expect_identical(choice$selected, "full")
  # An unscoreable reference cannot be validated: no transfer.
  choice <- select(score(c(NA, NA, NA, NA), c(NA, NA, NA, NA), names),
                   "fallback", 1 / 6)
  expect_identical(choice$selected, "fallback")
  expect_identical(choice$decision, "unvalidated")
  # An exact tie with the comparator is not a switch, even with margin 0.
  choice <- select(score(c(0.5, 0.5, 0.75, 0.75), c(1, 1, 0.9, 0.9), names),
                   "fallback", 0)
  expect_identical(choice$selected, "fallback")

  # Property: over random score tables, the selected candidate's balanced
  # accuracy is at least the reference's, and a switch away from the
  # comparator always clears the margin without a worse log-loss.
  manifoldalign:::.with_local_seed(20261019L, {
    for (i in seq_len(300L)) {
      accuracy <- sample((0:12) / 12, 4L, replace = TRUE)
      log_loss <- round(stats::runif(4L, 0.2, 1.5), 2)
      accuracy[stats::runif(4L) < 0.1] <- NA
      margin <- sample(c(0, 1 / 12, 1 / 6), 1L)
      table <- score(accuracy, log_loss, names)
      choice <- select(table, "fallback", margin)
      if (!is.finite(table["fallback", "accuracy"])) {
        expect_identical(choice$selected, "fallback")
        next
      }
      expect_gte(table[choice$selected, "accuracy"] + 1e-12,
                 table["fallback", "accuracy"])
      if (!identical(choice$selected, choice$comparator)) {
        expect_gte(
          table[choice$selected, "accuracy"] -
            table[choice$comparator, "accuracy"] + 1e-12,
          margin
        )
        expect_lte(table[choice$selected, "log_loss"],
                   table[choice$comparator, "log_loss"] + 1e-12)
      }
    }
  })
})

test_that("fallback inputs are validated and never silently ignored", {
  stub <- gate_stub_misleading()
  seed <- gate_stub_probability(wrong_rows = c(1L, 3L))
  fallback <- gate_stub_probability(wrong_rows = 1L)
  control <- gate_stub_control()
  validate <- manifoldalign:::.rapid_semisup_validate_relations
  expect_error(
    validate(stub$fit, list(seed), stub$labels, control,
             fallback_probabilities = list(fallback)),
    "not a single probability matrix"
  )
  expect_error(
    validate(stub$fit, list(list(list(rows = 1:12, probabilities = seed))),
             stub$labels, control, fallback_probabilities = list(fallback)),
    "out-of-fold `fallback`"
  )
  expect_error(
    validate(
      stub$fit,
      list(list(list(rows = 1:12, probabilities = seed, fallback = fallback))),
      stub$labels, control
    ),
    "provides no fallback"
  )
  expect_error(
    validate(stub$fit, list(seed), stub$labels,
             rapid_ma_semisup_control(rounds = 0L),
             fallback_probabilities = list(fallback)),
    "used only by the relation-safety gate"
  )

  fixture <- misleading_benchmark_fixture(2002L)
  fit <- rapid_ma(
    fixture$data, labels = fixture$labels, ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L, seed = 2002L)
  )
  probabilities <- fit$prediction_probabilities
  expect_error(
    rapid_ma_semisup(
      fit, control = list(rounds = 0L),
      fallback_probabilities = probabilities
    ),
    "used only by the relation-safety gate"
  )
  expect_error(
    rapid_ma_semisup(
      fit, validation_probabilities = probabilities,
      control = list(rounds = 0L, relation_validation = TRUE),
      fallback_probabilities = probabilities[1L]
    ),
    "one matrix \\(or NULL\\) per domain"
  )
})

test_that("under misleading target structure the teacher reduces to the fallback", {
  fixture <- misleading_benchmark_fixture(2002L)
  result <- suppressWarnings(benchmark_rapid_ma(
    fixture,
    methods = c("target_centroid", "rapid_ma_semisup"),
    ncomp = 8L,
    seed = 2002L,
    keep_fits = TRUE
  ))
  expect_identical(result$results$status, c("ok", "ok"))
  fit <- result$fits$rapid_ma_semisup
  gate <- fit$semisup$provenance$relation_validation
  for (domain in gate) {
    expect_identical(domain$reference, "fallback")
    expect_gte(
      domain$scores[domain$selected, "accuracy"] + 1e-12,
      domain$scores["fallback", "accuracy"]
    )
  }
  expect_identical(gate$target$selected, "fallback")
  expect_true(fit$semisup$provenance$fallback_selected[["target"]])

  # Exact oracle: with the fallback selected, the target teacher is the
  # target-only nearest centroid on every held-out row.
  test_rows <- result$split$target_test_rows
  classes <- result$provenance$class_levels
  local <- manifoldalign:::.rapid_benchmark_local_embedding(
    lapply(fixture$data, as.matrix), 8L
  )
  oracle <- manifoldalign:::.rapid_benchmark_centroid(
    list(local[[2L]]), list(fixture$labels[[2L]]), classes
  )[[1L]]
  expect_equal(
    unname(fit$prediction_probabilities$target[test_rows, classes]),
    unname(oracle$probability[test_rows, classes]),
    tolerance = 1e-12
  )
  expect_identical(
    unname(fit$predictions$target[test_rows]), oracle$prediction[test_rows]
  )
  rows <- split(result$results, result$results$method)
  expect_equal(rows$rapid_ma_semisup$oa, rows$target_centroid$oa)
  expect_equal(rows$rapid_ma_semisup$macro_f1, rows$target_centroid$macro_f1)
})

test_that("fallback folds are out-of-fold and held-out labels are never read", {
  fixture <- misleading_benchmark_fixture(2003L)
  data <- lapply(fixture$data, as.matrix)
  fit <- rapid_ma(
    data, labels = fixture$labels, positions = fixture$positions,
    attributes = fixture$attributes, ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L, seed = 2003L)
  )
  classes <- sort(unique(stats::na.omit(unlist(fixture$labels))),
                  method = "radix")
  seeds <- manifoldalign:::.rapid_benchmark_teacher_seeds(
    fit, data, fixture$labels, classes, 4L, 0.75
  )
  # Each fold's seed and fallback matrices are unchanged when the fold's own
  # label changes, so neither the row nor its neighbours used it.
  for (m in seq_along(fixture$labels)) {
    expect_false(is.null(seeds$fallback[[m]]))
    for (fold in head(seeds$validation[[m]], 3L)) {
      altered <- fixture$labels
      altered[[m]][fold$rows] <- setdiff(classes, altered[[m]][fold$rows])[[1L]]
      reseeded <- manifoldalign:::.rapid_benchmark_teacher_seeds(
        fit, data, altered, classes, 4L, 0.75
      )
      match_fold <- Filter(
        function(x) identical(x$rows, fold$rows), reseeded$validation[[m]]
      )
      expect_length(match_fold, 1L)
      expect_identical(match_fold[[1L]]$probabilities, fold$probabilities)
      expect_identical(match_fold[[1L]]$fallback, fold$fallback)
    }
  }

  # Inductive refinement: label values outside `train_mask` are discarded
  # before the gate, so changing them changes nothing.
  train_mask <- lapply(fixture$labels, function(x) !is.na(x))
  with_truth <- fixture$truth
  scrambled <- lapply(fixture$truth, function(x) {
    rotated <- c(classes[-1L], classes[[1L]])
    rotated[match(x, classes)]
  })
  for (m in seq_along(scrambled)) {
    scrambled[[m]][train_mask[[m]]] <- fixture$labels[[m]][train_mask[[m]]]
    with_truth[[m]][train_mask[[m]]] <- fixture$labels[[m]][train_mask[[m]]]
  }
  refine <- function(labels) {
    rapid_ma_semisup(
      fit, labels = labels, train_mask = train_mask, mode = "inductive",
      initial_probabilities = seeds$initial,
      validation_probabilities = seeds$validation,
      fallback_probabilities = seeds$fallback,
      control = list(
        rounds = 0L, propagation_steps = 2L, propagation_strength = 1,
        relation_validation = TRUE
      )
    )
  }
  first <- refine(with_truth)
  second <- refine(scrambled)
  expect_gt(sum(first$semisup$provenance$ignored_labels), 0L)
  expect_identical(first$prediction_probabilities, second$prediction_probabilities)
  expect_identical(
    first$semisup$provenance$relation_validation,
    second$semisup$provenance$relation_validation
  )

  # Benchmark: held-out truth is only used for scoring.
  original <- suppressWarnings(benchmark_rapid_ma(
    fixture, methods = "rapid_ma_semisup", ncomp = 4L, seed = 2003L,
    keep_fits = TRUE
  ))
  relabelled <- fixture
  test_rows <- original$split$target_test_rows
  relabelled$truth[[2L]][test_rows] <- scrambled[[2L]][test_rows]
  shuffled <- suppressWarnings(benchmark_rapid_ma(
    relabelled, methods = "rapid_ma_semisup", ncomp = 4L, seed = 2003L,
    keep_fits = TRUE
  ))
  expect_identical(
    original$fits$rapid_ma_semisup$prediction_probabilities,
    shuffled$fits$rapid_ma_semisup$prediction_probabilities
  )
  expect_identical(
    original$fits$rapid_ma_semisup$semisup$provenance$relation_validation,
    shuffled$fits$rapid_ma_semisup$semisup$provenance$relation_validation
  )
})
