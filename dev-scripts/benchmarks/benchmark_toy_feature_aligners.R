#!/usr/bin/env Rscript

# benchmark_toy_feature_aligners.R
# ------------------------------------------------------------
# Compare feature-domain alignment methods on the shared synthetic
# benchmark registry used by the toy benchmark tests.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
})

`%||%` <- function(a, b) if (is.null(a) || identical(a, "")) b else a
resolve_script_path <- function(default_path) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(sub("^--file=", "", file_arg[[1L]]))
  }

  frames <- sys.frames()
  if (length(frames) && !is.null(frames[[1L]]$ofile)) {
    return(frames[[1L]]$ofile)
  }

  default_path
}

script_path <- resolve_script_path("dev-scripts/benchmarks/benchmark_toy_feature_aligners.R")
here <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(here, "toy_feature_benchmark_common.R"))

main <- function() {
  out_results <- file.path(here, "toy_feature_aligners_results.csv")
  out_summary <- file.path(here, "toy_feature_aligners_summary.csv")
  out_regime_ranks <- file.path(here, "toy_feature_aligners_regime_ranks.csv")
  out_multiscale <- file.path(here, "toy_feature_aligners_multiscale_ablation.csv")
  out_method_catalog <- file.path(here, "toy_feature_aligners_method_catalog.csv")

  profile <- match.arg(Sys.getenv("MANIFOLDALIGN_TOY_PROFILE", "fast"), c("full", "fast"))
  seeds <- as.integer(strsplit(Sys.getenv("MANIFOLDALIGN_TOY_SEEDS", "1,2,3"), ",")[[1]])
  seeds <- seeds[is.finite(seeds)]
  scenario_txt <- Sys.getenv("MANIFOLDALIGN_TOY_SCENARIOS", "")
  scenarios <- manifoldalign:::synthetic_alignment_scenarios(profile = profile)$scenario
  if (nzchar(scenario_txt)) {
    wanted <- trimws(strsplit(scenario_txt, ",")[[1]])
    wanted <- wanted[nzchar(wanted)]
    missing <- setdiff(wanted, scenarios)
    if (length(missing)) {
      stop("Unknown benchmark scenario(s): ", paste(missing, collapse = ", "), call. = FALSE)
    }
    scenarios <- wanted
  }
  methods <- build_toy_feature_methods(include_multiscale_ablations = TRUE)

  if (!length(methods)) {
    stop("No benchmark methods available. Install multidesign/multivarious/kernlab as needed.", call. = FALSE)
  }

  cat("Toy feature aligner benchmark\n")
  cat("----------------------------\n")
  cat("profile: ", profile, "\n", sep = "")
  cat("seeds: ", paste(seeds, collapse = ","), "\n", sep = "")
  cat("scenarios: ", paste(scenarios, collapse = ", "), "\n", sep = "")
  cat("methods: ", paste(names(methods), collapse = ", "), "\n\n", sep = "")

  results <- manifoldalign:::run_synthetic_alignment_benchmark(
    methods = methods,
    scenarios = scenarios,
    profile = profile,
    seeds = seeds,
    score_fn = score_toy_feature_fit,
    verbose = TRUE
  )

  method_catalog <- manifoldalign:::synthetic_benchmark_method_catalog(methods)

  summary_group_cols <- intersect(
    c(
      "method", "scenario", "profile",
      "scenario_family", "scenario_latent_relation", "scenario_supervision",
      "scenario_side_information", "scenario_side_information_quality",
      "scenario_overlap", "scenario_n_views", "scenario_n_samples",
      "scenario_difficulty",
      "method_family", "supervision_regime", "side_information",
      "dimensionality_constraint", "variant_family", "variant", "backend",
      "tuned", "kernel", "landmarking", "scale_selection"
    ),
    names(results)
  )

  summary_df <- manifoldalign:::summarize_synthetic_alignment_benchmark(
    results[is.na(results$error), , drop = FALSE],
    metric_cols = c("runtime_sec", "max_abs_latent_cor", "centroid_acc"),
    group_cols = summary_group_cols
  )

  regime_rank_df <- summarize_toy_feature_regime_ranks(summary_df)
  multiscale_ablation_df <- regime_rank_df[regime_rank_df$method_family == "multiscale", , drop = FALSE]

  write_csv_safe(results, out_results)
  write_csv_safe(summary_df, out_summary)
  write_csv_safe(regime_rank_df, out_regime_ranks)
  write_csv_safe(multiscale_ablation_df, out_multiscale)
  write_csv_safe(method_catalog, out_method_catalog)

  cat("Wrote:\n")
  cat("  - ", out_results, "\n", sep = "")
  cat("  - ", out_summary, "\n", sep = "")
  cat("  - ", out_regime_ranks, "\n", sep = "")
  cat("  - ", out_multiscale, "\n", sep = "")
  cat("  - ", out_method_catalog, "\n\n", sep = "")

  cat("Regime leaders:\n")
  leaders <- regime_rank_df[
    !duplicated(regime_rank_df[c("scenario", "profile", "supervision_regime", "dimensionality_constraint")]),
    c("scenario", "supervision_regime", "dimensionality_constraint", "method", "overall_rank",
      "runtime_sec_mean", "max_abs_latent_cor_mean", "centroid_acc_mean"),
    drop = FALSE
  ]
  print(leaders)

  invisible(list(
    results = results,
    summary = summary_df,
    regime_ranks = regime_rank_df,
    multiscale_ablation = multiscale_ablation_df,
    method_catalog = method_catalog
  ))
}

if (sys.nframe() == 0) {
  main()
}
