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
  expect_false(comparison$split$test_labels_read_during_fit)
  expect_true(comparison$provenance$no_test_label_tuning)

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
