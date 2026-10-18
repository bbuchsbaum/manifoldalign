test_that("benchmark fixture is deterministic, structured, and provenance-rich", {
  first <- make_tiny_rapid_benchmark_fixture(seed = 19L)
  second <- make_tiny_rapid_benchmark_fixture(seed = 19L)

  expect_identical(first, second)
  expect_named(first, c(
    "data", "truth", "labels", "positions", "attributes",
    "correspondence", "metadata"
  ))
  expect_equal(length(first$correspondence), nrow(first$data[[1L]]))
  expect_true(all(first$truth[[1L]] ==
                    first$truth[[2L]][first$correspondence]))
  expect_true(anyNA(first$attributes[[2L]]$elevation))
  expect_match(first$metadata$license, "CC0")
})

test_that("dense LeMA reference obeys its equation-level constraints", {
  fixture <- make_tiny_rapid_benchmark_fixture(seed = 21L)
  fit <- manifoldalign:::.rapid_benchmark_lema(
    fixture$data,
    fixture$labels,
    fixture$correspondence,
    ncomp = 5L,
    seed = 21L,
    control = list(max_iter = 4L, admm_max_iter = 40L)
  )$fit

  expect_s3_class(fit, "lema_reference_fit")
  expect_false(fit$provenance$production_backend)
  expect_true(fit$provenance$dense_joint_graph)
  expect_lt(fit$constraints$orthogonality_error, 1e-7)
  expect_lt(fit$constraints$graph_symmetry_error, 1e-12)
  expect_gte(fit$constraints$graph_minimum, 0)
  expect_lt(fit$constraints$laplacian_row_sum_error, 1e-10)
  expect_equal(
    fit$history$total,
    fit$history$fidelity + fit$history$ridge + fit$history$graph,
    tolerance = 1e-9
  )
  expect_error(
    rapid_ma(
      fixture$data,
      labels = fixture$labels,
      control = list(backend = "lema_reference")
    ),
    "Unknown RAPID-MA control"
  )
})

test_that("benchmark uses explicit heads and keeps method failures visible", {
  fixture <- make_tiny_rapid_benchmark_fixture(seed = 23L)
  comparison <- benchmark_rapid_ma(
    fixture,
    methods = c(
      "target_centroid", "rapid_ma", "rapid_ma_native",
      "rapid_ma_semisup", "lema", "ssma"
    ),
    ncomp = 5L,
    seed = 23L
  )

  expect_s3_class(comparison, "rapid_ma_benchmark")
  expect_true(all(comparison$results$status == "ok"))
  expect_setequal(comparison$results$method, c(
    "target_centroid", "rapid_ma", "rapid_ma_native",
    "rapid_ma_semisup", "lema", "ssma"
  ))
  expect_identical(
    comparison$results$classifier[comparison$results$method == "rapid_ma"],
    "shared nearest centroid"
  )
  expect_match(
    comparison$results$classifier[
      comparison$results$method == "rapid_ma_semisup"
    ],
    "relation teacher"
  )
  expect_true(is.na(
    comparison$results$hits1[
      comparison$results$method == "target_centroid"
    ]
  ))
  expect_true(all(comparison$results$oa >= 0 & comparison$results$oa <= 1))
  expect_identical(
    comparison$provenance$class_levels,
    sort(unique(stats::na.omit(unlist(fixture$labels))), method = "radix")
  )
  expect_identical(
    unname(comparison$split$fit_label_counts),
    unname(vapply(fixture$labels, function(x) sum(!is.na(x)), integer(1)))
  )

  failed <- benchmark_rapid_ma(
    fixture,
    methods = "ssma",
    controls = list(ssma = list(not_a_control = TRUE)),
    ncomp = 4L
  )$results
  expect_identical(failed$status, "failed")
  expect_match(failed$error, "unused argument|Unknown")
})

test_that("held-out truth values cannot affect RAPID fitting", {
  original <- make_tiny_rapid_benchmark_fixture(seed = 29L)
  altered <- original
  held_out <- rapid_benchmark_held_out(original)
  altered$truth[[2L]][held_out] <- rev(altered$truth[[2L]][held_out])

  first <- benchmark_rapid_ma(
    original,
    methods = "rapid_ma_semisup",
    ncomp = 5L,
    seed = 29L,
    keep_fits = TRUE
  )
  second <- benchmark_rapid_ma(
    altered,
    methods = "rapid_ma_semisup",
    ncomp = 5L,
    seed = 29L,
    keep_fits = TRUE
  )

  expect_equal(
    first$fits$rapid_ma_semisup$scores,
    second$fits$rapid_ma_semisup$scores,
    tolerance = 1e-12
  )
  expect_equal(
    first$fits$rapid_ma_semisup$prediction_probabilities,
    second$fits$rapid_ma_semisup$prediction_probabilities,
    tolerance = 1e-12
  )
})

test_that("bounded matching is exact within its declared retrieval cap", {
  source <- rbind(c(0, 0), c(2, 0), c(5, 0))
  target <- rbind(c(0.1, 0), c(1.9, 0), c(4, 0), c(5.1, 0))
  correspondence <- c(1L, 2L, 4L)

  full <- manifoldalign:::.rapid_benchmark_matching(
    source, target, correspondence, rank_cap = 4L
  )
  capped <- manifoldalign:::.rapid_benchmark_matching(
    source, target, correspondence, rank_cap = 1L
  )

  expect_equal(full[["hits1"]], 1)
  expect_equal(full[["mrr"]], 1)
  expect_equal(full[["coverage"]], 1)
  expect_equal(capped[["hits1"]], 1)
  expect_equal(capped[["mrr"]], 1)
})

test_that("held-out metrics follow one shared class-universe rule", {
  truth <- c("a", "a", "b", "b", "c")
  prediction <- c("a", ".unknown", "b", "a", NA)
  probability <- rbind(
    c(0.85, 0.15), c(0.65, 0.35), c(0.95, 0.05), c(0.55, 0.45), c(0.75, 0.25)
  )
  metrics <- manifoldalign:::.rapid_classification_metrics(
    truth, prediction, probability
  )
  # Universe {.unknown, a, b, c}. a: P = R = 1/2, IoU = 1/3; b: P = 1,
  # R = 1/2, IoU = 1/2; c (NA prediction) and .unknown (pure false positive):
  # F1 = IoU = 0.
  expect_equal(metrics[["oa"]], 2 / 5)
  expect_equal(metrics[["macro_f1"]], (1 / 2 + 2 / 3) / 4)
  expect_equal(metrics[["miou"]], (1 / 3 + 1 / 2) / 4)
  expect_equal(metrics[["rare_recall"]], 0)
  expect_equal(metrics[["ece"]], (0.15 + 0.65 + 0.05 + 0.55 + 0.75) / 5)

  # A labelled class predicted but absent from the held-out truth is charged.
  closed <- manifoldalign:::.rapid_classification_metrics(
    c("a", "b"), c("a", "b"), diag(2L)
  )
  charged <- manifoldalign:::.rapid_classification_metrics(
    c("a", "b"), c("a", "z"), diag(2L)
  )
  expect_equal(closed[["macro_f1"]], 1)
  expect_equal(charged[["macro_f1"]], 1 / 3)
  expect_equal(charged[["miou"]], 1 / 3)
})

test_that("held-out-only classes never enter the prediction heads", {
  fixture <- make_tiny_rapid_benchmark_fixture(seed = 31L)
  altered <- fixture
  held_out <- which(rapid_benchmark_held_out(fixture))
  altered$truth[[2L]][held_out[1:3]] <- "swamp"
  methods <- c("target_centroid", "rapid_ma", "rapid_ma_semisup")
  first <- benchmark_rapid_ma(
    fixture, methods = methods, ncomp = 4L, seed = 31L, keep_fits = TRUE
  )
  second <- benchmark_rapid_ma(
    altered, methods = methods, ncomp = 4L, seed = 31L, keep_fits = TRUE
  )
  expect_true(all(second$results$status == "ok"))
  expect_identical(second$provenance$class_levels, c("forest", "urban", "water"))
  expect_equal(second$fits$rapid_ma$s, first$fits$rapid_ma$s, tolerance = 1e-12)
  expect_equal(
    second$fits$rapid_ma_semisup$prediction_probabilities,
    first$fits$rapid_ma_semisup$prediction_probabilities,
    tolerance = 1e-12
  )
  # The unseen class is scored (never predicted), lowering macro-F1.
  expect_true(all(second$results$macro_f1 < first$results$macro_f1))
  expect_null(rownames(fixture$data[[1L]]))
  expect_null(rownames(fixture$data[[2L]]))
})

test_that("the target-only comparator uses every feature column", {
  fixture <- make_tiny_rapid_benchmark_fixture(seed = 33L)
  permuted <- fixture
  set.seed(330)
  order <- sample.int(ncol(fixture$data[[2L]]))
  permuted$data[[2L]] <- fixture$data[[2L]][, order, drop = FALSE]
  first <- benchmark_rapid_ma(
    fixture, methods = "target_centroid", ncomp = 3L, seed = 33L
  )$results
  second <- benchmark_rapid_ma(
    permuted, methods = "target_centroid", ncomp = 3L, seed = 33L
  )$results
  for (metric in c("oa", "macro_f1", "miou", "rare_recall", "ece")) {
    expect_equal(second[[metric]], first[[metric]], tolerance = 1e-8)
  }
  embedding <- manifoldalign:::.rapid_benchmark_local_embedding(
    list(fixture$data[[2L]]), 3L
  )[[1L]]
  standardized <- scale(fixture$data[[2L]])
  expect_identical(ncol(embedding), 3L)
  expect_equal(
    sum(embedding^2),
    sum(svd(standardized)$d[1:3]^2),
    tolerance = 1e-10
  )
})

test_that("teacher validation folds are genuinely out-of-fold", {
  fixture <- make_tiny_rapid_benchmark_fixture(seed = 37L)
  fit <- rapid_ma(
    fixture$data, labels = fixture$labels, ncomp = 4L,
    control = list(max_iter = 1L, min_iter = 1L, seed = 37L)
  )
  classes <- sort(unique(stats::na.omit(unlist(fixture$labels))),
                  method = "radix")
  seeds <- manifoldalign:::.rapid_benchmark_teacher_seeds(
    fit, fixture$data, fixture$labels, classes, 4L, 0.75
  )
  for (m in seq_along(fixture$labels)) {
    folds <- seeds$validation[[m]]
    expect_gt(length(folds), 0L)
    for (fold in head(folds, 3L)) {
      row <- fold$rows
      expect_length(row, 1L)
      altered <- fixture$labels
      altered[[m]][row] <- setdiff(classes, altered[[m]][row])[[1L]]
      reseeded <- manifoldalign:::.rapid_benchmark_teacher_seeds(
        fit, fixture$data, altered, classes, 4L, 0.75
      )
      match_fold <- Filter(
        function(x) identical(x$rows, row), reseeded$validation[[m]]
      )
      expect_length(match_fold, 1L)
      # Every row of the fold's matrix, including relation neighbours, is
      # unchanged when the held-out label changes.
      expect_identical(match_fold[[1L]]$probabilities, fold$probabilities)
    }
  }
})

test_that("benchmark fixtures and runs restore the caller's random stream", {
  set.seed(4040)
  expected <- stats::runif(2L)
  set.seed(4040)
  make_tiny_rapid_benchmark_fixture(seed = 41L)
  expect_identical(stats::runif(2L), expected)
})
