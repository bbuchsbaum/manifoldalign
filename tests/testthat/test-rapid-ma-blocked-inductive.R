test_that("spatial blocks are deterministic and geometry-equivariant", {
  position <- cbind(
    x = c(0, 1, 2, 8, 9, 10, 0, 5, 10),
    y = c(0, 0, 0, 5, 5, 5, 10, 10, 10)
  )

  block <- manifoldalign:::.rapid_spatial_blocks(position, grid = c(3L, 2L))
  permutation <- c(9, 2, 5, 1, 7, 3, 8, 6, 4)
  transformed <- sweep(position, 2L, c(4, 7), "*")
  transformed <- sweep(transformed, 2L, c(11, 30), "+")

  expect_length(block, nrow(position))
  expect_false(anyNA(block))
  expect_lte(length(unique(block)), 6L)
  expect_identical(
    manifoldalign:::.rapid_spatial_blocks(
      position[permutation, , drop = FALSE], grid = c(3L, 2L)
    )[order(permutation)],
    block
  )
  expect_identical(
    manifoldalign:::.rapid_spatial_blocks(transformed, grid = c(3L, 2L)),
    block
  )
})

test_that("blocked split is strictly disjoint and cannot read test labels", {
  position <- as.matrix(expand.grid(x = seq_len(12), y = seq_len(8)))
  labels <- rep(c("forest", "urban", "water"), length.out = nrow(position))
  block <- manifoldalign:::.rapid_spatial_blocks(position, grid = c(4L, 2L))
  test_block <- max(block)
  calibration_block <- sort(unique(block))[[length(unique(block)) - 1L]]

  first <- manifoldalign:::.rapid_blocked_inductive_split(
    labels, position,
    test_block = test_block,
    calibration_block = calibration_block,
    grid = c(4L, 2L),
    labels_per_class = 3L,
    calibration_per_class = 2L,
    seed = 73L
  )
  altered <- labels
  altered[first$test_rows] <- rev(altered[first$test_rows])
  second <- manifoldalign:::.rapid_blocked_inductive_split(
    altered, position,
    test_block = test_block,
    calibration_block = calibration_block,
    grid = c(4L, 2L),
    labels_per_class = 3L,
    calibration_per_class = 2L,
    seed = 73L
  )

  expect_identical(first$fit_rows, second$fit_rows)
  expect_identical(first$labeled_rows, second$labeled_rows)
  expect_identical(first$calibration_rows, second$calibration_rows)
  expect_identical(first$test_rows, second$test_rows)
  expect_length(intersect(first$fit_rows, first$calibration_rows), 0L)
  expect_length(intersect(first$fit_rows, first$test_rows), 0L)
  expect_length(intersect(first$calibration_rows, first$test_rows), 0L)
  expect_true(all(first$labeled_rows %in% first$fit_rows))
  expect_true(all(table(labels[first$labeled_rows]) <= 3L))
  expect_true(all(table(labels[first$calibration_rows]) <= 2L))
  expect_false(first$provenance$test_labels_read_during_split)
})

test_that("blocked split caps unlabeled fit state without dropping labels", {
  position <- as.matrix(expand.grid(x = seq_len(100), y = seq_len(40)))
  labels <- rep(c("a", "b", "c", "d"), length.out = nrow(position))
  split <- manifoldalign:::.rapid_blocked_inductive_split(
    labels, position,
    test_block = 20L,
    calibration_block = 19L,
    grid = c(5L, 4L),
    labels_per_class = 4L,
    calibration_per_class = 3L,
    max_fit_rows = 250L,
    seed = 81L
  )

  expect_lte(length(split$fit_rows), 250L)
  expect_true(all(split$labeled_rows %in% split$fit_rows))
  expect_equal(length(split$block_id), nrow(position))
  expect_lt(as.numeric(object.size(split)), 20 * as.numeric(object.size(position)))
})

test_that("temperature scaling is normalized and cannot worsen calibration NLL", {
  probability <- rbind(
    c(0.99, 0.01), c(0.98, 0.02), c(0.97, 0.03),
    c(0.03, 0.97), c(0.02, 0.98), c(0.01, 0.99)
  )
  colnames(probability) <- c("a", "b")
  labels <- c("b", "a", "a", "a", "b", "b")

  fit <- manifoldalign:::.rapid_temperature_fit(probability, labels)
  calibrated <- manifoldalign:::.rapid_temperature_apply(
    probability, fit$temperature
  )
  identity <- manifoldalign:::.rapid_temperature_apply(probability, 1)

  expect_equal(identity, probability, tolerance = 1e-12)
  expect_true(is.finite(fit$temperature))
  expect_gte(fit$temperature, 1)
  expect_lte(fit$nll_after, fit$nll_before + 1e-12)
  expect_equal(rowSums(calibrated), rep(1, nrow(calibrated)), tolerance = 1e-12)
  expect_identical(max.col(calibrated), max.col(probability))
})

test_that("inductive interpolation is exact, bounded, and permutation stable", {
  training <- rbind(c(0, 0), c(1, 0), c(0, 1), c(2, 2))
  values <- rbind(c(1, 0), c(0.8, 0.2), c(0.7, 0.3), c(0, 1))
  query <- rbind(training[2, ], c(1.8, 1.9))

  first <- manifoldalign:::.rapid_inductive_interpolate(
    training, values, query, k = 2L
  )
  permutation <- c(4, 2, 1, 3)
  second <- manifoldalign:::.rapid_inductive_interpolate(
    training[permutation, , drop = FALSE],
    values[permutation, , drop = FALSE],
    query,
    k = 2L
  )

  expect_equal(first$value[1, ], values[2, ], tolerance = 1e-12)
  expect_equal(first$value, second$value, tolerance = 1e-12)
  expect_equal(rowSums(first$value), rep(1, nrow(query)), tolerance = 1e-12)
  expect_equal(ncol(first$neighbours), 2L)
})

test_that("blocked validation rejects ambiguous or invalid geometry", {
  position <- cbind(x = 1:8, y = rep(1:2, each = 4))
  labels <- rep(c("a", "b"), 4)

  expect_error(
    manifoldalign:::.rapid_spatial_blocks(
      rbind(position, c(NA_real_, 1)), grid = c(2L, 2L)
    ),
    "finite"
  )
  expect_error(
    manifoldalign:::.rapid_blocked_inductive_split(
      labels, position,
      test_block = 1L,
      calibration_block = 1L,
      grid = c(2L, 2L)
    ),
    "different"
  )
})
