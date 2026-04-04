library(testthat)
library(manifoldalign)

test_that("synthetic benchmark registry exposes stable scenario metadata", {
  reg_full <- manifoldalign:::synthetic_alignment_scenarios(profile = "full")
  reg_fast <- manifoldalign:::synthetic_alignment_scenarios(profile = "fast")

  expect_equal(
    reg_full$scenario,
    c(
      "linear_affine", "isometric_curve", "hard_nonisometric",
      "partial_open_set", "noisy_label_shift", "many_domain_consensus",
      "highdim_stress"
    )
  )
  expect_equal(
    reg_fast$scenario,
    c(
      "linear_affine", "isometric_curve", "hard_nonisometric",
      "partial_open_set", "noisy_label_shift"
    )
  )
  expect_true(all(reg_full$family == "feature"))
  common <- merge(
    reg_fast[, c("scenario", "n_samples")],
    reg_full[, c("scenario", "n_samples")],
    by = "scenario",
    suffixes = c("_fast", "_full"),
    sort = FALSE
  )
  expect_true(all(common$n_samples_fast < common$n_samples_full))
  expect_true(all(c("overlap", "side_information_quality", "n_views") %in% names(reg_full)))
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

test_that("partial-overlap scenarios preserve per-view row IDs and observed labels", {
  skip_if_not_installed("multidesign")

  scen <- manifoldalign:::synthetic_alignment_scenario("partial_open_set", profile = "fast", seed = 41)
  expect_true(length(intersect(scen$row_ids$view1, scen$row_ids$view2)) < length(scen$row_ids$view1))
  expect_true(any(!scen$row_ids$view1 %in% scen$row_ids$view2))

  hd <- manifoldalign:::synthetic_alignment_to_hyperdesign(
    scen,
    views = c("view1", "view2"),
    label_col = "condition",
    id_col = "sample_id"
  )

  expect_equal(as.character(hd[[1]]$design$sample_id), as.character(scen$row_ids$view1))
  expect_equal(as.character(hd[[2]]$design$sample_id), as.character(scen$row_ids$view2))
  expect_equal(length(intersect(hd[[1]]$design$sample_id, hd[[2]]$design$sample_id)), scen$n_shared)
})

test_that("noisy-label scenarios separate observed and true labels", {
  scen <- manifoldalign:::synthetic_alignment_scenario("noisy_label_shift", profile = "fast", seed = 17)

  expect_true(any(as.character(scen$observed_labels) != as.character(scen$true_labels)))
  expect_equal(length(scen$observed_labels_by_view), 3L)
  expect_true(all(vapply(scen$observed_labels_by_view, length, integer(1)) == nrow(scen$view1)))
  expect_equal(levels(scen$observed_labels), levels(scen$true_labels))
})

test_that("many-domain and high-dimensional scenarios expose expected structure", {
  many <- manifoldalign:::synthetic_alignment_scenario("many_domain_consensus", profile = "full", seed = 11)
  highdim <- manifoldalign:::synthetic_alignment_scenario("highdim_stress", profile = "full", seed = 12)

  expect_equal(length(many$view_names), 4L)
  expect_true(all(vapply(many$view_names, function(view) ncol(many[[view]]), integer(1)) == 2L))

  hd_dims <- vapply(highdim$view_names, function(view) ncol(highdim[[view]]), integer(1))
  hd_rows <- vapply(highdim$view_names, function(view) nrow(highdim[[view]]), integer(1))
  expect_true(all(hd_dims > hd_rows))
  expect_true(all(c("feature_dims", "row_ids", "true_labels_by_view") %in% names(highdim)))
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
  expect_true(all(c("n_runs", "n_success", "n_errors", "error_rate") %in% names(summary_df)))
})

test_that("synthetic benchmark runner is invariant to method order for stochastic methods", {
  stochastic_methods_1 <- list(
    alpha = function(scenario) matrix(stats::rnorm(6), 3, 2),
    beta = function(scenario) matrix(stats::rnorm(6), 3, 2)
  )
  stochastic_methods_2 <- rev(stochastic_methods_1)

  score_fn <- function(fit, scenario, method) {
    c(sum_fit = sum(fit), mean_fit = mean(fit))
  }

  res1 <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = stochastic_methods_1,
    scenarios = "linear_affine",
    profile = "fast",
    seeds = 7L,
    score_fn = score_fn
  )
  res2 <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = stochastic_methods_2,
    scenarios = "linear_affine",
    profile = "fast",
    seeds = 7L,
    score_fn = score_fn
  )

  res1 <- res1[order(res1$method), c("method", "sum_fit", "mean_fit")]
  res2 <- res2[order(res2$method), c("method", "sum_fit", "mean_fit")]
  rownames(res1) <- NULL
  rownames(res2) <- NULL

  expect_equal(res1, res2)
})

test_that("synthetic benchmark runner propagates method and scenario metadata", {
  methods <- list(
    tagged = list(
      fit = function(scenario) matrix(1, 3, 2),
      meta = list(
        method_family = "toy_family",
        supervision_regime = "toy_supervision",
        tuned = FALSE
      )
    )
  )

  res <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = methods,
    scenarios = "linear_affine",
    profile = "fast",
    seeds = 9L,
    score_fn = function(fit, scenario, method) c(n_rows = nrow(scenario$view1))
  )

  expect_equal(res$method_family, "toy_family")
  expect_equal(res$supervision_regime, "toy_supervision")
  expect_false(res$tuned)
  expect_equal(res$scenario_family, "feature")
  expect_equal(res$scenario_latent_relation, "linear")
  expect_equal(res$scenario_difficulty, "easy")

  catalog <- manifoldalign:::synthetic_benchmark_method_catalog(methods)
  expect_equal(catalog$method, "tagged")
  expect_equal(catalog$method_family, "toy_family")
})

test_that("synthetic benchmark summary aligns aggregates by group keys", {
  results <- data.frame(
    method = c("beta", "alpha", "beta", "alpha"),
    scenario = c("s2", "s1", "s2", "s1"),
    profile = "fast",
    seed = c(1L, 1L, 2L, 2L),
    runtime_sec = c(1, 10, 3, 14),
    error = NA_character_,
    stringsAsFactors = FALSE
  )

  summary_df <- manifoldalign:::summarize_synthetic_alignment_benchmark(
    results,
    metric_cols = "runtime_sec"
  )

  expect_equal(summary_df$method, c("beta", "alpha"))
  expect_equal(summary_df$scenario, c("s2", "s1"))
  expect_equal(summary_df$runtime_sec_mean, c(2, 12))
  expect_equal(round(summary_df$runtime_sec_sd, 8), c(round(sqrt(2), 8), round(sqrt(8), 8)))
  expect_equal(summary_df$n_runs, c(2L, 2L))
  expect_equal(summary_df$n_success, c(2L, 2L))
  expect_equal(summary_df$error_rate, c(0, 0))
})

test_that("synthetic benchmark seed robustness summary reports spread and failures", {
  results <- data.frame(
    method = c("alpha", "alpha", "alpha", "beta"),
    scenario = c("s1", "s1", "s1", "s1"),
    profile = "fast",
    seed = 1:4,
    runtime_sec = c(1, 2, 3, 4),
    centroid_acc = c(0.8, 0.9, NA, 0.5),
    error = c(NA, NA, "fit failed", NA),
    stringsAsFactors = FALSE
  )

  robust <- manifoldalign:::summarize_synthetic_alignment_seed_robustness(
    results,
    metric_cols = c("runtime_sec", "centroid_acc")
  )

  alpha <- robust[robust$method == "alpha", , drop = FALSE]
  expect_equal(alpha$n_runs, 3L)
  expect_equal(alpha$n_success, 2L)
  expect_equal(alpha$n_errors, 1L)
  expect_equal(alpha$error_rate, 1 / 3)
  expect_equal(alpha$centroid_acc_mean, 0.85)
  expect_equal(alpha$centroid_acc_median, 0.85)
  expect_equal(alpha$centroid_acc_iqr, 0.05)
  expect_true(is.finite(alpha$centroid_acc_cv))
})
