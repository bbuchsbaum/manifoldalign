#!/usr/bin/env Rscript

# Frozen exact-node retrieval benchmark for RAPID-MA's bounded matcher.
# Correspondence truth is used only after every fit and match has completed.

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- if (length(args)) args[[1L]] else {
  file.path("docs", "rapid-ma-matching")
}
profile <- Sys.getenv("RAPID_MA_PROFILE", "full")
if (!profile %in% c("quick", "full")) {
  stop("RAPID_MA_PROFILE must be quick or full.")
}
default_seed_start <- if (profile == "quick") "7101" else "7201"
seed_start <- as.integer(Sys.getenv("RAPID_MA_SEED_START", default_seed_start))
seed_count <- as.integer(Sys.getenv(
  "RAPID_MA_SEED_COUNT", if (profile == "quick") "3" else "20"
))
n_common <- as.integer(Sys.getenv("RAPID_MA_MATCH_N_COMMON", "300"))
if (any(!is.finite(c(seed_start, seed_count, n_common))) ||
    seed_start < 0L || seed_count < 1L || n_common < 12L) {
  stop("Matching seed and size controls must be valid positive integers.")
}
seeds <- seq.int(seed_start, length.out = seed_count)

if (!"manifoldalign" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(manifoldalign)
  }
}

score_matches <- function(target, truth) {
  matched <- !is.na(target)
  correct <- matched & target == truth
  c(
    hits1 = mean(correct),
    matched_hits1 = if (any(matched)) mean(correct[matched]) else NA_real_,
    mrr = mean(correct),
    coverage = mean(matched)
  )
}

measure_match <- function(fit, arguments) {
  gc(FALSE)
  start <- proc.time()[["elapsed"]]
  value <- do.call(rapid_ma_match, c(list(fit = fit), arguments))
  elapsed <- proc.time()[["elapsed"]] - start
  list(
    value = value,
    runtime_sec = elapsed,
    retained_mb = as.numeric(utils::object.size(value)) / 1024^2
  )
}

failure_row <- function(seed, method, message) {
  data.frame(
    seed = seed,
    method = method,
    status = "failed",
    assignment_mode = NA_character_,
    uses_alignment = NA,
    uses_position = NA,
    fit_runtime_sec = NA_real_,
    match_runtime_sec = NA_real_,
    matcher_retained_mb = NA_real_,
    hits1 = NA_real_, matched_hits1 = NA_real_, mrr = NA_real_,
    coverage = NA_real_, candidate_edges = NA_real_, candidate_cap = NA_real_,
    dense_pairwise_allocated = NA,
    correspondence_read_during_fit_or_match = FALSE,
    error = message,
    stringsAsFactors = FALSE
  )
}

rows <- list()
row_index <- 0L
for (seed in seeds) {
  fixture <- rapid_ma_benchmark_fixture(
    n_common = n_common,
    extra_target = as.integer(ceiling(0.2 * n_common)),
    labels_per_class = 4L,
    seed = seed
  )
  comparison <- benchmark_rapid_ma(
    fixture,
    methods = c("rapid_ma_semisup", "lema", "ssma", "spectral_mnn", "token_ot"),
    ncomp = 8L,
    seed = seed,
    keep_fits = TRUE
  )
  rapid_row <- comparison$results[
    comparison$results$method == "rapid_ma_semisup", , drop = FALSE
  ]
  fit <- comparison$fits$rapid_ma_semisup

  configurations <- list(
    rapid_prototype_one_to_one = list(
      candidate_views = character(0),
      assignment = "one_to_one",
      position_weight = 0
    ),
    rapid_multiview_one_to_one = list(assignment = "one_to_one"),
    rapid_multiview_independent = list(assignment = "independent")
  )
  generated <- list()
  if (!is.null(fit) && identical(rapid_row$status, "ok")) {
    for (name in names(configurations)) {
      generated[[name]] <- tryCatch(
        measure_match(fit, configurations[[name]]),
        error = identity
      )
    }
  }
  position_start <- proc.time()[["elapsed"]]
  position_state <- if (is.null(fit)) NULL else list(
    source = fit$relations$source$views$position$view,
    target = fit$relations$target$views$position$view
  )
  position_metrics <- if (is.null(position_state$source) ||
                          is.null(position_state$target)) {
    NULL
  } else {
    manifoldalign:::.rapid_benchmark_matching(
      position_state$source,
      position_state$target,
      fixture$correspondence,
      rank_cap = 100L
    )
  }
  position_runtime <- proc.time()[["elapsed"]] - position_start

  # Correspondence truth is first accessed for custom matcher scoring here,
  # after all RAPID fits and candidate generation have completed.
  truth <- fixture$correspondence
  for (name in names(configurations)) {
    row_index <- row_index + 1L
    measured <- generated[[name]]
    if (is.null(measured) || inherits(measured, "error")) {
      rows[[row_index]] <- failure_row(
        seed, name,
        if (inherits(measured, "error")) conditionMessage(measured) else
          "RAPID fit unavailable"
      )
      next
    }
    matching <- measured$value
    metric <- score_matches(matching$matches$target, truth)
    rows[[row_index]] <- data.frame(
      seed = seed,
      method = name,
      status = "ok",
      assignment_mode = matching$assignment_mode,
      uses_alignment = TRUE,
      uses_position = "position" %in% matching$candidate_views ||
        matching$weights[["position"]] > 0,
      fit_runtime_sec = rapid_row$runtime_sec,
      match_runtime_sec = measured$runtime_sec,
      matcher_retained_mb = measured$retained_mb,
      hits1 = metric[["hits1"]],
      matched_hits1 = metric[["matched_hits1"]],
      mrr = metric[["mrr"]],
      coverage = metric[["coverage"]],
      candidate_edges = matching$candidate_edges,
      candidate_cap = matching$candidate_cap,
      dense_pairwise_allocated = matching$dense_pairwise_allocated,
      correspondence_read_during_fit_or_match = FALSE,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  for (method in c(
    "rapid_ma_semisup", "lema", "ssma", "spectral_mnn", "token_ot"
  )) {
    row_index <- row_index + 1L
    value <- comparison$results[
      comparison$results$method == method, , drop = FALSE
    ]
    if (!nrow(value) || !identical(value$status, "ok")) {
      rows[[row_index]] <- failure_row(
        seed, paste0(method, "_embedding"),
        if (nrow(value)) value$error else "benchmark result unavailable"
      )
      next
    }
    rows[[row_index]] <- data.frame(
      seed = seed,
      method = paste0(method, "_embedding"),
      status = "ok",
      assignment_mode = "independent ranked retrieval",
      uses_alignment = TRUE,
      uses_position = FALSE,
      fit_runtime_sec = value$runtime_sec,
      match_runtime_sec = NA_real_,
      matcher_retained_mb = NA_real_,
      hits1 = value$hits1,
      matched_hits1 = value$hits1,
      mrr = value$mrr,
      coverage = value$coverage,
      candidate_edges = n_common * min(100L, nrow(fixture$data$target)),
      candidate_cap = min(100L, nrow(fixture$data$target)),
      dense_pairwise_allocated = FALSE,
      correspondence_read_during_fit_or_match = identical(method, "lema"),
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  row_index <- row_index + 1L
  if (is.null(position_metrics)) {
    rows[[row_index]] <- failure_row(
      seed, "position_knn_diagnostic", "position view unavailable"
    )
  } else {
    rows[[row_index]] <- data.frame(
      seed = seed,
      method = "position_knn_diagnostic",
      status = "ok",
      assignment_mode = "independent ranked retrieval",
      uses_alignment = FALSE,
      uses_position = TRUE,
      fit_runtime_sec = 0,
      match_runtime_sec = position_runtime,
      matcher_retained_mb = NA_real_,
      hits1 = position_metrics[["hits1"]],
      matched_hits1 = position_metrics[["hits1"]],
      mrr = position_metrics[["mrr"]],
      coverage = position_metrics[["coverage"]],
      candidate_edges = n_common * min(100L, nrow(fixture$data$target)),
      candidate_cap = min(100L, nrow(fixture$data$target)),
      dense_pairwise_allocated = FALSE,
      correspondence_read_during_fit_or_match = FALSE,
      error = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  message("completed matching seed=", seed)
}

raw <- do.call(rbind, rows)
rownames(raw) <- NULL
groups <- split(seq_len(nrow(raw)), raw$method)
summary <- do.call(rbind, lapply(groups, function(index) {
  block <- raw[index, , drop = FALSE]
  ok <- block$status == "ok"
  output <- data.frame(
    method = block$method[[1L]],
    n_runs = nrow(block),
    n_ok = sum(ok),
    n_failed = sum(!ok),
    stringsAsFactors = FALSE
  )
  for (metric in c(
    "fit_runtime_sec", "match_runtime_sec", "matcher_retained_mb",
    "hits1", "matched_hits1", "mrr", "coverage", "candidate_edges"
  )) {
    value <- block[[metric]][ok]
    finite <- value[is.finite(value)]
    output[[paste0(metric, "_mean")]] <- if (length(finite)) {
      mean(finite)
    } else NA_real_
    output[[paste0(metric, "_sd")]] <- if (length(finite) > 1L) {
      stats::sd(finite)
    } else NA_real_
  }
  output
}))
rownames(summary) <- NULL
summary <- summary[order(-summary$hits1_mean, summary$method), ]

dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
raw_path <- paste0(output_prefix, "-raw.csv")
summary_path <- paste0(output_prefix, "-summary.csv")
metadata_path <- paste0(output_prefix, "-metadata.txt")
utils::write.csv(raw, raw_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)

source_files <- c(
  "R/adapters_rapid_ma.R", "R/rapid_ma_benchmark.R",
  "inst/benchmarks/rapid_ma_matching.R"
)
source_hash <- tools::md5sum(source_files[file.exists(source_files)])
git_commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(...) NA_character_
)
if (length(git_commit) != 1L || !nzchar(git_commit)) git_commit <- NA_character_
metadata <- c(
  "benchmark=RAPID-MA bounded exact-node matching",
  paste0(
    "evaluation_role=",
    if (profile == "quick") "development" else "frozen held-out seeds"
  ),
  "development_seeds=7101,7102,7103",
  paste0("implementation_commit=", git_commit),
  paste0("profile=", profile),
  paste0("seeds=", paste(seeds, collapse = ",")),
  paste0("n_common=", n_common),
  paste0("n_target=", n_common + ceiling(0.2 * n_common)),
  "labels_per_class=4",
  "rapid_matching_truth_read_during_fit_or_match=false",
  "lema_receives_paired_landmark_correspondence=true",
  "rapid_candidate_cap=32",
  "rapid_view_candidate_k=12",
  "dense_cross_domain_pair_matrix=false",
  paste0("R.version=", R.version.string),
  paste0("platform=", R.version$platform),
  paste0(names(source_hash), " md5=", source_hash)
)
writeLines(metadata, metadata_path)

print(summary[, c(
  "method", "n_ok", "hits1_mean", "matched_hits1_mean", "mrr_mean",
  "coverage_mean", "fit_runtime_sec_mean", "match_runtime_sec_mean"
)], row.names = FALSE)
