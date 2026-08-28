#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  package_root <- Sys.getenv("MANIFOLDALIGN_PACKAGE_ROOT", ".")
  helper_root <- Sys.getenv("MANIFOLDALIGN_HELPER_ROOT", ".")
  devtools::load_all(package_root, quiet = TRUE)
})

source(file.path(helper_root, "dev-scripts/benchmarks/toy_feature_benchmark_common.R"))

args <- commandArgs(trailingOnly = TRUE)
profile <- if (length(args) >= 1L) args[[1L]] else "full"
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path("artifacts", "kema-eigencore-benchmark-2026-08-27", paste0("toy-feature-", profile))
seeds <- if (length(args) >= 3L) as.integer(strsplit(args[[3L]], ",", fixed = TRUE)[[1L]]) else 1:3

profile <- match.arg(profile, c("fast", "full"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_methods <- build_toy_feature_methods(include_multiscale_ablations = FALSE)
method_names <- c(
  "kema_linear_label",
  "kema_rbf_full_label",
  "kema_rbf_reduced_label"
)
methods <- all_methods[method_names]
if (!identical(names(methods), method_names)) {
  stop("The expected KEMA benchmark variants are not all available.", call. = FALSE)
}

score_with_fidelity <- function(fit, scenario, method) {
  base <- score_toy_feature_fit(fit, scenario, method)
  c(
    as.list(base),
    list(
      fidelity_passed = isTRUE(fit$fidelity$passed),
      max_rel_residual = fit$fidelity$max_rel_residual,
      max_b_orth_offdiag = fit$fidelity$max_b_orth_offdiag,
      eigensolver = fit$fidelity$eigensolver
    )
  )
}

scenarios <- manifoldalign:::synthetic_alignment_scenarios(profile = profile)$scenario
results <- manifoldalign:::run_synthetic_alignment_benchmark(
  methods = methods,
  scenarios = scenarios,
  profile = profile,
  seeds = seeds,
  score_fn = score_with_fidelity,
  verbose = TRUE
)

group_cols <- intersect(
  c(
    "method", "scenario", "profile",
    "scenario_family", "scenario_latent_relation", "scenario_supervision",
    "scenario_side_information", "scenario_side_information_quality",
    "scenario_overlap", "scenario_n_views", "scenario_n_samples",
    "scenario_difficulty", "method_family", "supervision_regime",
    "side_information", "dimensionality_constraint", "variant_family",
    "variant", "backend", "tuned", "kernel", "landmarking",
    "scale_selection"
  ),
  names(results)
)
ok <- is.na(results$error)
summary <- manifoldalign:::summarize_synthetic_alignment_benchmark(
  results[ok, , drop = FALSE],
  metric_cols = c(
    "runtime_sec", "max_abs_latent_cor", "centroid_acc",
    "max_rel_residual", "max_b_orth_offdiag"
  ),
  group_cols = group_cols
)

utils::write.csv(results, file.path(out_dir, "raw_results.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "summary.csv"), row.names = FALSE)

receipt <- c(
  paste0("generated: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
  paste0("package_root: ", normalizePath(package_root, mustWork = TRUE)),
  paste0("helper_root: ", normalizePath(helper_root, mustWork = TRUE)),
  paste0("profile: ", profile),
  paste0("methods: ", paste(method_names, collapse = ",")),
  paste0("scenarios: ", paste(scenarios, collapse = ",")),
  paste0("seeds: ", paste(seeds, collapse = ",")),
  paste0("eigencore: ", as.character(utils::packageVersion("eigencore"))),
  paste0("successful_fits: ", sum(ok)),
  paste0("failed_fits: ", sum(!ok)),
  paste0("successful_fidelity_passes: ", sum(results$fidelity_passed %in% TRUE, na.rm = TRUE))
)
writeLines(receipt, file.path(out_dir, "receipt.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

print(summary[, intersect(
  c(
    "method", "scenario", "runtime_sec_mean", "max_abs_latent_cor_mean",
    "centroid_acc_mean", "max_rel_residual_mean", "max_b_orth_offdiag_mean"
  ),
  names(summary)
), drop = FALSE])
