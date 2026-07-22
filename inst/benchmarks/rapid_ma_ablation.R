#!/usr/bin/env Rscript

# Frozen structure and prediction-head ablation for the degraded-target fixture.

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args)) args[[1L]] else {
  file.path("docs", "rapid-ma-ablation")
}
seeds <- 4001:4020
if (nzchar(Sys.getenv("RAPID_MA_SEEDS"))) {
  seeds <- as.integer(strsplit(Sys.getenv("RAPID_MA_SEEDS"), ",", fixed = TRUE)[[1L]])
}

if (!"manifoldalign" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(manifoldalign)
  }
}

degraded_fixture <- function(seed) {
  fixture <- rapid_ma_benchmark_fixture(seed = seed)
  set.seed(900000L + seed)
  fixture$data[[2L]] <- fixture$data[[2L]] + matrix(
    stats::rnorm(length(fixture$data[[2L]]), sd = 0.75),
    nrow(fixture$data[[2L]])
  )
  fixture$metadata$perturbation <-
    "target feature Gaussian noise, sd=0.75"
  fixture
}

run_variant <- function(fixture, variant, seed) {
  method <- "rapid_ma_semisup"
  controls <- list()
  if (identical(variant, "target_centroid")) method <- "target_centroid"
  if (identical(variant, "dense_lema")) method <- "lema"
  if (identical(variant, "no_position")) fixture$positions <- NULL
  if (identical(variant, "no_attribute")) fixture$attributes <- NULL
  if (identical(variant, "feature_relations_only")) {
    fixture$positions <- NULL
    fixture$attributes <- NULL
  }
  if (identical(variant, "misleading_target_structure")) {
    set.seed(800000L + seed)
    permutation <- sample.int(nrow(fixture$data[[2L]]))
    fixture$positions[[2L]] <- fixture$positions[[2L]][permutation, , drop = FALSE]
    fixture$attributes[[2L]] <- fixture$attributes[[2L]][permutation, , drop = FALSE]
  }
  if (identical(variant, "no_relation_propagation")) {
    controls <- list(rapid_ma_semisup = list(
      teacher = list(propagation_steps = 0L, propagation_strength = 0)
    ))
  }
  if (identical(variant, "domain_local_seed_only")) {
    controls <- list(rapid_ma_semisup = list(seed_local_weight = 1))
  }
  if (identical(variant, "aligned_seed_only")) {
    controls <- list(rapid_ma_semisup = list(seed_local_weight = 0))
  }
  result <- benchmark_rapid_ma(
    fixture,
    methods = method,
    ncomp = 8L,
    seed = seed,
    controls = controls,
    keep_fits = identical(method, "rapid_ma_semisup")
  )
  row <- result$results
  row$selected_relation_candidate <- NA_character_
  row$selected_target_relations <- NA_character_
  row$relation_safety_switched <- NA
  if (identical(method, "rapid_ma_semisup")) {
    validation <- result$fits[[method]]$semisup$provenance$
      relation_validation[[2L]]
    row$selected_relation_candidate <- validation$selected
    row$selected_target_relations <- paste(
      validation$selected_relations, collapse = "+"
    )
    row$relation_safety_switched <- validation$switched
  }
  row$variant <- variant
  row$replicate_seed <- seed
  row
}

variants <- c(
  "rapid_full", "target_centroid", "dense_lema", "no_position",
  "no_attribute", "feature_relations_only", "misleading_target_structure",
  "no_relation_propagation", "domain_local_seed_only", "aligned_seed_only"
)
rows <- list()
index <- 0L
for (seed in seeds) {
  fixture <- degraded_fixture(seed)
  for (variant in variants) {
    index <- index + 1L
    rows[[index]] <- suppressWarnings(run_variant(fixture, variant, seed))
  }
  message("completed seed=", seed)
}
raw <- do.call(rbind, rows)
rownames(raw) <- NULL

metrics <- c(
  "runtime_sec", "retained_mb", "endpoint_rss_mb", "rss_delta_mb",
  "oa", "macro_f1", "miou", "rare_recall", "ece",
  "hits1", "mrr", "coverage"
)
groups <- split(seq_len(nrow(raw)), raw$variant)
summary <- do.call(rbind, lapply(groups, function(rows) {
  block <- raw[rows, , drop = FALSE]
  output <- data.frame(
    variant = block$variant[[1L]],
    method = block$method[[1L]],
    classifier = block$classifier[[which(!is.na(block$classifier))[[1L]]]],
    n_replicates = nrow(block),
    n_ok = sum(block$status == "ok"),
    stringsAsFactors = FALSE
  )
  for (metric in metrics) {
    values <- block[[metric]][block$status == "ok"]
    output[[paste0(metric, "_mean")]] <- if (any(is.finite(values))) {
      mean(values, na.rm = TRUE)
    } else {
      NA_real_
    }
    output[[paste0(metric, "_sd")]] <- if (sum(is.finite(values)) > 1L) {
      stats::sd(values, na.rm = TRUE)
    } else {
      NA_real_
    }
  }
  output
}))
rownames(summary) <- NULL
summary <- summary[order(-summary$macro_f1_mean), ]

raw_path <- paste0(output_prefix, "-raw.csv")
summary_path <- paste0(output_prefix, "-summary.csv")
metadata_path <- paste0(output_prefix, "-metadata.txt")
dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(raw, raw_path, row.names = FALSE, na = "")
utils::write.csv(summary, summary_path, row.names = FALSE, na = "")

source_paths <- c(
  "R/rapid_ma_benchmark.R", "R/rapid_ma.R", "R/rapid_ma_semisup.R",
  "inst/benchmarks/lema_reference.R", "inst/benchmarks/rapid_ma_ablation.R"
)
source_paths <- source_paths[file.exists(source_paths)]
metadata <- c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE), sep = "="),
  paste("seeds", paste(seeds, collapse = ","), sep = "="),
  "development_seeds=201:205,2001:2020",
  "target_feature_noise_sd=0.75",
  "test_label_tuning=false",
  paste("R.version", R.version.string, sep = "="),
  paste("platform", R.version$platform, sep = "="),
  paste(names(tools::md5sum(source_paths)), tools::md5sum(source_paths), sep = " md5=")
)
writeLines(metadata, metadata_path)

print(summary[, c(
  "variant", "n_ok", "oa_mean", "macro_f1_mean", "miou_mean",
  "rare_recall_mean", "ece_mean", "runtime_sec_mean"
)], row.names = FALSE)
