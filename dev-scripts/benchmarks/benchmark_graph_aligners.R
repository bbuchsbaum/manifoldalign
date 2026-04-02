#!/usr/bin/env Rscript

# benchmark_graph_aligners.R
# ------------------------------------------------------------
# Compare graph alignment methods available in manifoldalign:
#   - token_ot_graph_align (Token-OT Graph Align)
#   - fpgw (Fused-Partial Gromov-Wasserstein)
#   - parrot
#   - cone_align
#   - grasp
#
# Writes:
#   - graph_aligners_benchmark_results.csv
#   - graph_aligners_benchmark_summary.csv
# ------------------------------------------------------------

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION") &&
      grepl("^Package:\\s*manifoldalign", paste(readLines("DESCRIPTION", n = 5L), collapse = "\n"))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(manifoldalign)
})

write_csv_safe <- function(df, path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  utils::write.csv(df, file = path, row.names = FALSE)
  invisible(path)
}

summarize_results <- function(res) {
  stopifnot(is.data.frame(res))
  res <- res[is.na(res$error), , drop = FALSE]
  if (!nrow(res)) return(res)

  agg <- aggregate(
    cbind(top1_accuracy, runtime_sec) ~ method + structure + n,
    data = res,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))
  )
  data.frame(
    method = agg$method,
    structure = agg$structure,
    n = agg$n,
    top1_mean = agg$top1_accuracy[, "mean"],
    top1_sd = agg$top1_accuracy[, "sd"],
    runtime_mean_sec = agg$runtime_sec[, "mean"],
    runtime_sd_sec = agg$runtime_sec[, "sd"],
    stringsAsFactors = FALSE
  )
}

main <- function() {
  here <- dirname(normalizePath(sys.frame(1)$ofile %||% "dev-scripts/benchmarks/benchmark_graph_aligners.R", mustWork = FALSE))
  out_results <- file.path(here, "graph_aligners_benchmark_results.csv")
  out_summary <- file.path(here, "graph_aligners_benchmark_summary.csv")

  # Tuned to be representative across canonical graph families without turning
  # this script into a full nightly benchmark sweep.
  sizes <- c(50L, 100L, 200L)
  structures <- c("ring", "grid", "community")
  n_reps <- 3L
  n_anchors <- 10L

  ctrl <- list(
    candidate_k = 80L,
    n_levels = 2L,
    prior_mode = "hard",
    coarsen_method = "louvain",
    views = "raw",
    token_mode = "view_only"
  )

  cat("Graph aligner benchmark\n")
  cat("-----------------------\n")
  cat("sizes: ", paste(sizes, collapse = ","), "\n", sep = "")
  cat("structures: ", paste(structures, collapse = ","), "\n", sep = "")
  cat("n_reps: ", n_reps, "\n", sep = "")
  cat("n_anchors: ", n_anchors, "\n", sep = "")
  cat("token_ot settings: ", paste(names(ctrl), unlist(ctrl), sep = "=", collapse = ", "), "\n\n", sep = "")

  res <- benchmark_graph_alignment_methods(
    sizes = sizes,
    d = 3L,
    noise_sd = 0.01,
    structure = structures,
    permute_fraction = 1,
    n_anchors = n_anchors,
    methods = c("token_ot_graph", "fpgw", "parrot", "cone_align", "grasp"),
    candidate_k = ctrl$candidate_k,
    n_levels = ctrl$n_levels,
    prior_mode = ctrl$prior_mode,
    coarsen_method = ctrl$coarsen_method,
    views = ctrl$views,
    token_mode = ctrl$token_mode,
    fpgw_omega1 = 0.5,
    fpgw_epsilon = 0.01,
    fpgw_max_iter = 50L,
    fpgw_inner_max_iter = 20L,
    fpgw_tol = 1e-6,
    n_reps = n_reps,
    seed = 1L,
    verbose = FALSE
  )

  sum_df <- summarize_results(res)

  write_csv_safe(res, out_results)
  write_csv_safe(sum_df, out_summary)

  cat("Wrote:\n")
  cat("  - ", out_results, "\n", sep = "")
  cat("  - ", out_summary, "\n\n", sep = "")

  print(sum_df)

  invisible(list(results = res, summary = sum_df))
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || identical(a, "")) b else a

if (sys.nframe() == 0) {
  main()
}
