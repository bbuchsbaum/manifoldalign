library(testthat)
library(manifoldalign)

test_that("synthetic benchmark registry exposes stable scenario metadata", {
  reg_full <- manifoldalign:::synthetic_alignment_scenarios(profile = "full")
  reg_fast <- manifoldalign:::synthetic_alignment_scenarios(profile = "fast")

  expect_equal(
    reg_full$scenario,
    c("linear_affine", "isometric_curve", "hard_nonisometric")
  )
  expect_equal(reg_fast$scenario, reg_full$scenario)
  expect_true(all(reg_full$family == "feature"))
  expect_true(all(reg_fast$n_samples < reg_full$n_samples))
})

test_that("synthetic benchmark scenarios support parameter overrides and metadata", {
  scen <- manifoldalign:::synthetic_alignment_scenario(
    "linear_affine",
    profile = "full",
    seed = 99,
    n = 18,
    noise = 0.05
  )

  expect_equal(nrow(scen$latent), 18)
  expect_equal(length(scen$view_names), 3)
  expect_identical(scen$benchmark_meta$scenario, "linear_affine")
  expect_identical(scen$benchmark_meta$profile, "full")
  expect_equal(scen$benchmark_meta$generator_args$seed, 99)
  expect_equal(scen$benchmark_meta$generator_args$n, 18)
})

test_that("synthetic benchmark hyperdesign conversion preserves rows and labels", {
  skip_if_not_installed("multidesign")

  scen <- manifoldalign:::synthetic_alignment_scenario("isometric_curve", profile = "fast", seed = 7)
  hd <- manifoldalign:::synthetic_alignment_to_hyperdesign(
    scen,
    views = c("view1", "view2"),
    label_col = "condition",
    id_col = "row_id"
  )

  expect_s3_class(hd, "hyperdesign")
  expect_equal(length(hd), 2)
  expect_equal(nrow(hd[[1]]$x), nrow(scen$view1))
  expect_equal(nrow(hd[[2]]$x), nrow(scen$view2))
  expect_true(all(c("row_id", "condition") %in% names(hd[[1]]$design)))
  expect_equal(as.character(hd[[1]]$design$condition), as.character(scen$true_labels))
})

test_that("synthetic benchmark runner and summary return a unified result schema", {
  methods <- list(
    ok = function(scenario) scenario$view1,
    fail = function(scenario) stop("intentional failure")
  )

  res <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = methods,
    scenarios = c("linear_affine", "isometric_curve"),
    profile = "fast",
    seeds = 1:2,
    score_fn = function(fit, scenario, method) {
      c(
        n_rows = nrow(scenario$view1),
        n_cols = ncol(as.matrix(fit)),
        latent_dim = if (is.matrix(scenario$latent)) ncol(scenario$latent) else 1L
      )
    }
  )

  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 8)
  expect_true(all(c("method", "scenario", "profile", "seed", "runtime_sec", "error") %in% names(res)))
  expect_true(all(res$method %in% c("ok", "fail")))
  expect_true(any(!is.na(res$error)))
  expect_true(any(is.na(res$error)))

  ok_rows <- res[res$method == "ok", , drop = FALSE]
  expect_true(all(is.finite(ok_rows$n_rows)))
  expect_true(all(is.finite(ok_rows$n_cols)))

  summary_df <- manifoldalign:::summarize_synthetic_alignment_benchmark(
    ok_rows,
    metric_cols = c("runtime_sec", "n_rows", "n_cols", "latent_dim")
  )

  expect_true(is.data.frame(summary_df))
  expect_equal(nrow(summary_df), 2)
  expect_true(all(c("runtime_sec_mean", "n_rows_mean", "n_cols_mean", "latent_dim_mean") %in% names(summary_df)))
})
