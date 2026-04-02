library(testthat)
library(manifoldalign)

test_that("graph benchmark registry enumerates scenario metadata cleanly", {
  reg <- manifoldalign:::synthetic_graph_alignment_registry(
    sizes = c(12L, 18L),
    structures = c("ring", "grid"),
    d = 4L,
    noise_sd = 0.02,
    permute_fraction = 0.5,
    n_anchors = 3L,
    n_reps = 2L,
    seed = 11L
  )

  expect_true(is.data.frame(reg))
  expect_equal(nrow(reg), 8)
  expect_true(all(c(
    "scenario_family", "scenario", "structure", "n", "d", "noise_sd",
    "permute_fraction", "n_anchors", "rep", "generation_seed"
  ) %in% names(reg)))
  expect_true(all(reg$scenario_family == "graph"))
  expect_true(all(reg$structure %in% c("ring", "grid")))
})

test_that("graph benchmark case returns consistent maps and anchor correspondences", {
  case <- manifoldalign:::synthetic_graph_alignment_case(
    n = 20L,
    d = 3L,
    structure = "community",
    noise_sd = 0.01,
    permute_fraction = 0.6,
    n_anchors = 4L,
    seed = 22L
  )

  expect_equal(dim(case$X1), c(20, 3))
  expect_equal(dim(case$X2), c(20, 3))
  expect_equal(length(case$map12), 20)
  expect_equal(sort(case$map12), seq_len(20))
  expect_equal(length(case$anchor_idx), 4)
  expect_true(is.data.frame(case$correspondences))
  expect_equal(nrow(case$correspondences), 4)
  expect_s3_class(case$hd, "hyperdesign")
  expect_identical(case$benchmark_meta$structure, "community")
})

test_that("graph benchmark output includes shared scenario metadata columns", {
  res <- benchmark_graph_alignment_methods(
    sizes = 18L,
    d = 3L,
    noise_sd = 0.01,
    structure = "grid",
    permute_fraction = 1,
    n_anchors = 3L,
    methods = "token_ot_graph",
    candidate_k = 8L,
    n_levels = 1L,
    n_reps = 1L,
    seed = 7L,
    verbose = FALSE
  )

  expect_true(is.data.frame(res))
  expect_true(all(c(
    "method", "scenario_family", "scenario", "structure", "n", "d",
    "noise_sd", "permute_fraction", "n_anchors", "rep",
    "generation_seed", "decode_mode", "runtime_sec", "top1_accuracy",
    "coverage", "error"
  ) %in% names(res)))
  expect_true(all(res$scenario_family == "graph"))
  expect_true(all(res$structure == "grid"))
  expect_true(all(res$decode_mode == "native"))
})

test_that("graph benchmark can sweep multiple structures in one call", {
  res <- benchmark_graph_alignment_methods(
    sizes = 16L,
    d = 3L,
    noise_sd = 0.01,
    structure = c("ring", "community"),
    permute_fraction = 1,
    n_anchors = 3L,
    methods = "token_ot_graph",
    candidate_k = 8L,
    n_levels = 1L,
    n_reps = 1L,
    seed = 17L,
    verbose = FALSE
  )

  expect_true(is.data.frame(res))
  expect_setequal(unique(res$structure), c("ring", "community"))
  expect_equal(nrow(res), 2L)
  expect_true(all(grepl("^(ring|community)_16_1$", res$scenario)))
})
