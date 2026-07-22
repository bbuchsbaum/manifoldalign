#!/usr/bin/env Rscript

# Frozen same-split accuracy benchmark for RAPID-MA and package baselines.
# The default seeds (3001:3020) are disjoint from the development seeds used
# while constructing the fixture and benchmark harness.

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args)) args[[1L]] else {
  file.path("docs", "rapid-ma-validation")
}
profile <- Sys.getenv("RAPID_MA_PROFILE", "full")
seeds <- if (identical(profile, "quick")) 3001:3003 else 3001:3020
if (nzchar(Sys.getenv("RAPID_MA_SEEDS"))) {
  seeds <- as.integer(strsplit(Sys.getenv("RAPID_MA_SEEDS"), ",", fixed = TRUE)[[1L]])
}
if (any(!is.finite(seeds)) || !length(seeds)) {
  stop("RAPID_MA_SEEDS must be a comma-separated integer sequence.")
}

if (!"manifoldalign" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(manifoldalign)
  }
}

methods <- c(
  "target_centroid", "rapid_ma", "rapid_ma_native", "rapid_ma_semisup",
  "lema", "ssma", "kema", "multiscale", "spectral_mnn", "token_ot"
)
method_timeout_sec <- as.numeric(Sys.getenv("RAPID_MA_METHOD_TIMEOUT", "5"))
if (!is.finite(method_timeout_sec) || method_timeout_sec <= 0) {
  stop("RAPID_MA_METHOD_TIMEOUT must be a positive number of seconds.")
}

make_scenario <- function(seed, scenario) {
  fixture <- rapid_ma_benchmark_fixture(seed = seed)
  if (identical(scenario, "degraded_target")) {
    set.seed(900000L + seed)
    fixture$data[[2L]] <- fixture$data[[2L]] + matrix(
      stats::rnorm(length(fixture$data[[2L]]), sd = 0.75),
      nrow(fixture$data[[2L]])
    )
    fixture$metadata$perturbation <-
      "target feature Gaussian noise, sd=0.75; positions and attributes unchanged"
  }
  fixture
}

failed_method_row <- function(method, elapsed, peak, peak_delta, message) {
  data.frame(
    method = method,
    status = "failed",
    evaluation_mode = NA_character_,
    classifier = NA_character_,
    runtime_sec = elapsed,
    retained_mb = NA_real_,
    endpoint_rss_mb = NA_real_,
    rss_delta_mb = NA_real_,
    oa = NA_real_,
    macro_f1 = NA_real_,
    miou = NA_real_,
    rare_recall = NA_real_,
    ece = NA_real_,
    hits1 = NA_real_,
    mrr = NA_real_,
    coverage = NA_real_,
    error = message,
    sampled_peak_rss_mb = peak,
    sampled_peak_delta_mb = peak_delta,
    stringsAsFactors = FALSE
  )
}

sample_process_rss <- function(pid) {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle(pid))[["rss"]]) / 1024^2,
    error = function(...) NA_real_
  )
}

run_isolated_method <- function(fixture, method, seed, timeout_sec) {
  start <- proc.time()[["elapsed"]]
  job <- parallel::mcparallel(
    benchmark_rapid_ma(
      fixture,
      methods = method,
      ncomp = 8L,
      seed = seed
    ),
    name = paste(method, seed, sep = "-"),
    mc.set.seed = FALSE,
    silent = TRUE
  )
  initial <- sample_process_rss(job$pid)
  peak <- initial
  value <- NULL
  repeat {
    collected <- parallel::mccollect(job, wait = FALSE)
    current <- sample_process_rss(job$pid)
    if (is.finite(current)) peak <- max(peak, current, na.rm = TRUE)
    if (!is.null(collected)) {
      value <- collected[[1L]]
      break
    }
    elapsed <- proc.time()[["elapsed"]] - start
    if (elapsed >= timeout_sec) {
      parallel:::mckill(job, signal = 15L)
      parallel::mccollect(job, wait = TRUE, timeout = 1)
      delta <- if (is.finite(initial) && is.finite(peak)) peak - initial else NA_real_
      return(failed_method_row(
        method, elapsed, peak, delta,
        paste0("method exceeded isolated ", timeout_sec, " second timeout")
      ))
    }
    Sys.sleep(0.02)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  delta <- if (is.finite(initial) && is.finite(peak)) peak - initial else NA_real_
  if (inherits(value, "try-error") ||
      !is.list(value) || is.null(value$results)) {
    return(failed_method_row(
      method, elapsed, peak, delta,
      if (inherits(value, "try-error")) as.character(value) else
        "isolated method returned no benchmark result"
    ))
  }
  row <- value$results
  row$sampled_peak_rss_mb <- peak
  row$sampled_peak_delta_mb <- delta
  row
}

rows <- list()
index <- 0L
for (scenario in c("structured", "degraded_target")) {
  for (seed in seeds) {
    fixture <- make_scenario(seed, scenario)
    method_rows <- lapply(methods, function(method) {
      suppressWarnings(run_isolated_method(
        fixture, method, seed, method_timeout_sec
      ))
    })
    block <- do.call(rbind, method_rows)
    rownames(block) <- NULL
    block$scenario <- scenario
    block$replicate_seed <- seed
    block$n_source <- nrow(fixture$data[[1L]])
    block$n_target <- nrow(fixture$data[[2L]])
    block$labels_per_class <- 4L
    block$test_labels_read_during_fit <- FALSE
    index <- index + 1L
    rows[[index]] <- block
    message(
      "completed scenario=", scenario,
      " seed=", seed,
      " failures=", sum(block$status != "ok")
    )
  }
}
raw <- do.call(rbind, rows)
rownames(raw) <- NULL

metric_names <- c(
  "runtime_sec", "retained_mb", "endpoint_rss_mb", "rss_delta_mb",
  "sampled_peak_rss_mb", "sampled_peak_delta_mb",
  "oa", "macro_f1", "miou", "rare_recall", "ece",
  "hits1", "mrr", "coverage"
)
groups <- split(seq_len(nrow(raw)), interaction(
  raw$scenario, raw$method, drop = TRUE, lex.order = TRUE
))
summary_rows <- lapply(groups, function(rows) {
  block <- raw[rows, , drop = FALSE]
  classifier_index <- which(!is.na(block$classifier))
  output <- data.frame(
    scenario = block$scenario[[1L]],
    method = block$method[[1L]],
    classifier = if (length(classifier_index)) {
      block$classifier[[classifier_index[[1L]]]]
    } else {
      NA_character_
    },
    n_replicates = nrow(block),
    n_ok = sum(block$status == "ok"),
    n_failed = sum(block$status != "ok"),
    stringsAsFactors = FALSE
  )
  for (metric in metric_names) {
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
})
summary <- do.call(rbind, summary_rows)
rownames(summary) <- NULL
summary <- summary[order(summary$scenario, -summary$macro_f1_mean), ]

raw_path <- paste0(output_prefix, "-raw.csv")
summary_path <- paste0(output_prefix, "-summary.csv")
metadata_path <- paste0(output_prefix, "-metadata.txt")
dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(raw, raw_path, row.names = FALSE, na = "")
utils::write.csv(summary, summary_path, row.names = FALSE, na = "")

source_paths <- c(
  "R/rapid_ma_benchmark.R", "R/rapid_ma.R", "R/rapid_ma_semisup.R",
  "inst/benchmarks/lema_reference.R", "inst/benchmarks/rapid_ma_accuracy.R"
)
source_paths <- source_paths[file.exists(source_paths)]
metadata <- c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE), sep = "="),
  paste("profile", profile, sep = "="),
  paste("seeds", paste(seeds, collapse = ","), sep = "="),
  paste("method_timeout_sec", method_timeout_sec, sep = "="),
  "development_seeds=1:5,9,17,19,21,23,29,101:120,201:205,1001:1020",
  "test_label_tuning=false",
  "shared_embedding_classifier=nearest_centroid",
  paste(
    "rapid_teacher_preset=seed_local_weight:0.75,rounds:0,",
    "propagation_steps:2,propagation_strength:1,relation_validation:true,",
    "relation_validation_margin:0.1666667",
    sep = ""
  ),
  paste(
    "memory_note=sampled peak RSS uses isolated forked methods at 20 ms intervals;",
    "absolute RSS includes inherited package state"
  ),
  paste("R.version", R.version.string, sep = "="),
  paste("platform", R.version$platform, sep = "="),
  paste(names(tools::md5sum(source_paths)), tools::md5sum(source_paths), sep = " md5=")
)
writeLines(metadata, metadata_path)

print(summary[, c(
  "scenario", "method", "n_ok", "oa_mean", "macro_f1_mean", "miou_mean",
  "rare_recall_mean", "ece_mean", "hits1_mean", "mrr_mean", "runtime_sec_mean"
)], row.names = FALSE)
