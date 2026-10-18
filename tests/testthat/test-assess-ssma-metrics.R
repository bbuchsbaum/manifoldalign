# assess_ssma() metric contract: denominators, unavailable decoders,
# Procrustes with translation, sampling of singleton candidate sets.

test_that("top-1 accuracy counts missing and out-of-range predictions as wrong", {
  res <- manifoldalign:::.assess_ssma_top1(c(1L, NA, 3L, 99L), 1:4)
  expect_equal(res$accuracy, 0.5)
  expect_equal(res$coverage, 0.5)
  all_missing <- manifoldalign:::.assess_ssma_top1(rep(NA_integer_, 3), 1:3)
  expect_equal(all_missing$accuracy, 0)
  expect_equal(all_missing$coverage, 0)
  sub <- manifoldalign:::.assess_ssma_top1(c(1L, NA, 3L, 4L), 1:4, idx = c(2L, 4L))
  expect_equal(sub$accuracy, 0.5)
})

test_that("an unavailable assignment solver is reported as unavailable", {
  E <- matrix(stats::rnorm(12), 6, 2)
  lap <- manifoldalign:::.assess_ssma_lap_assignment(E, E, clue_available = FALSE)
  expect_identical(lap$status, "unavailable")
  ok <- manifoldalign:::.assess_ssma_lap_assignment(E, E, clue_available = TRUE)
  expect_identical(ok$status, "ok")
  expect_identical(ok$prediction, 1:6)
})

test_that("Procrustes decoding applies the translation it estimates", {
  set.seed(3)
  E1 <- matrix(stats::rnorm(40), 20, 2)
  theta <- 0.7
  Q0 <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  shift <- c(5, -3)
  # E2 is a rotated and shifted copy of E1: E1 = E2 Q0' - shift Q0'.
  E2 <- sweep(E1 %*% Q0, 2L, shift, "+")
  anchors <- c(2L, 5L, 9L, 14L)          # a subset whose mean is not zero
  map <- manifoldalign:::.assess_ssma_procrustes(E1, E2, anchors, anchors)
  E2p <- manifoldalign:::.assess_ssma_apply_procrustes(map, E2)
  expect_lt(max(abs(E2p - E1)), 1e-10)
  # Rotation alone (the old behaviour) leaves a residual offset.
  expect_gt(max(abs(E2 %*% map$Q - E1)), 1)
})

test_that("anchor sampling never reinterprets a single candidate as 1:x", {
  for (sd in 1:30) {
    set.seed(sd)
    expect_identical(manifoldalign:::.assess_ssma_sample(37L, 1L), 37L)
  }
  set.seed(1)
  draw <- manifoldalign:::.assess_ssma_sample(c(4L, 9L, 11L), 2L)
  expect_true(all(draw %in% c(4L, 9L, 11L)))
})

test_that("assess_ssma reports decoder status and leaves the caller RNG alone", {
  skip_if_not_installed("clue")
  set.seed(8)
  before <- .Random.seed
  res <- assess_ssma(
    sizes = 12L, d = 3L, noise_sd = 0.02, structure = "ring",
    permute_fraction = 0.9, n_anchors = 1L, holdout_fraction = 0,
    n_reps = 1L, solvers = "reduced", serial = FALSE, ssma_knn = 5L,
    ssma_rank_per_domain = 3L, ncomp = 2L, decode_k = 2L
  )
  expect_identical(.Random.seed, before)
  expect_true(all(c("lap_status", "top1_procrustes_lap_coverage") %in% names(res)))
  expect_equal(res$n_train_anchors, 1L)
  expect_true(res$lap_status %in% c("ok", "not_run"))
  res2 <- assess_ssma(
    sizes = 12L, d = 3L, noise_sd = 0.02, structure = "ring",
    permute_fraction = 0.9, n_anchors = 1L, holdout_fraction = 0,
    n_reps = 1L, solvers = "reduced", serial = FALSE, ssma_knn = 5L,
    ssma_rank_per_domain = 3L, ncomp = 2L, decode_k = 2L
  )
  expect_equal(res$top1_raw_nn, res2$top1_raw_nn)
})
