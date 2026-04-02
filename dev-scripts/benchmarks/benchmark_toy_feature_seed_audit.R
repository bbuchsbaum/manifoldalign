#!/usr/bin/env Rscript

# benchmark_toy_feature_seed_audit.R
# ------------------------------------------------------------
# Audit per-seed stability for the synthetic toy feature benchmark.
# Intended to distinguish genuine seed sensitivity from reporting bugs.
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

script_path <- resolve_script_path("dev-scripts/benchmarks/benchmark_toy_feature_seed_audit.R")
here <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(here, "toy_feature_benchmark_common.R"))

parse_method_filter <- function(methods, requested = NULL) {
  all_names <- names(methods)
  if (is.null(requested) || !nzchar(requested)) {
    requested_names <- unique(c(
      "ssma_exact",
      grep("^gpca_", all_names, value = TRUE),
      grep("^kema_", all_names, value = TRUE),
      grep("^multiscale_", all_names, value = TRUE)
    ))
    return(methods[requested_names])
  }

  requested_names <- trimws(strsplit(requested, ",")[[1]])
  requested_names <- requested_names[nzchar(requested_names)]
  missing <- setdiff(requested_names, all_names)
  if (length(missing)) {
    stop("Unknown seed-audit method(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  methods[requested_names]
}

main <- function() {
  out_results <- file.path(here, "toy_feature_seed_audit_results.csv")
  out_summary <- file.path(here, "toy_feature_seed_audit_summary.csv")
  out_regime_ranks <- file.path(here, "toy_feature_seed_audit_regime_ranks.csv")

  profile <- match.arg(Sys.getenv("MANIFOLDALIGN_TOY_PROFILE", "fast"), c("full", "fast"))
  seed_txt <- Sys.getenv("MANIFOLDALIGN_SEED_AUDIT_SEEDS", paste(1:12, collapse = ","))
  seeds <- as.integer(strsplit(seed_txt, ",")[[1]])
  seeds <- seeds[is.finite(seeds)]

  methods <- build_toy_feature_methods(include_multiscale_ablations = TRUE)
  methods <- parse_method_filter(methods, Sys.getenv("MANIFOLDALIGN_SEED_AUDIT_METHODS", ""))
  scenarios <- manifoldalign:::synthetic_alignment_scenarios(profile = profile)$scenario

  if (!length(methods)) {
    stop("No seed-audit methods available after filtering.", call. = FALSE)
  }

  cat("Toy feature seed audit\n")
  cat("----------------------\n")
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

  group_cols <- intersect(
    c(
      "method", "scenario", "profile",
      "scenario_family", "scenario_latent_relation", "scenario_supervision",
      "scenario_side_information", "scenario_difficulty",
      "method_family", "supervision_regime", "side_information",
      "dimensionality_constraint", "variant_family", "variant", "backend",
      "tuned", "kernel", "landmarking", "scale_selection"
    ),
    names(results)
  )

  summary_df <- manifoldalign:::summarize_synthetic_alignment_seed_robustness(
    results = results,
    metric_cols = c("runtime_sec", "max_abs_latent_cor", "centroid_acc"),
    group_cols = group_cols
  )
  regime_rank_df <- summarize_toy_feature_regime_ranks(summary_df)

  write_csv_safe(results, out_results)
  write_csv_safe(summary_df, out_summary)
  write_csv_safe(regime_rank_df, out_regime_ranks)

  cat("Wrote:\n")
  cat("  - ", out_results, "\n", sep = "")
  cat("  - ", out_summary, "\n", sep = "")
  cat("  - ", out_regime_ranks, "\n\n", sep = "")

  cat("Highest centroid-accuracy variability:\n")
  variability <- summary_df[order(-summary_df$centroid_acc_cv, summary_df$method), c(
    "scenario", "method", "supervision_regime", "centroid_acc_mean",
    "centroid_acc_sd", "centroid_acc_iqr", "centroid_acc_cv",
    "max_abs_latent_cor_mean", "max_abs_latent_cor_sd"
  ), drop = FALSE]
  print(utils::head(variability, 10L))

  invisible(list(
    results = results,
    summary = summary_df,
    regime_ranks = regime_rank_df
  ))
}

if (sys.nframe() == 0) {
  main()
}
