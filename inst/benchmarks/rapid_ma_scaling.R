#!/usr/bin/env Rscript

# Isolated complete-pipeline scaling comparison. Each task constructs the same
# deterministic fixture inside a forked process, runs one method, and is sampled
# for RSS every 20 ms. Timeout failures are data, not dropped observations.

args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(
  VECLIB_MAXIMUM_THREADS = "1",
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1"
)
output_prefix <- if (length(args)) args[[1L]] else {
  file.path("docs", "rapid-ma-scaling")
}
profile <- Sys.getenv("RAPID_MA_PROFILE", "full")
timeout_sec <- as.numeric(Sys.getenv("RAPID_MA_METHOD_TIMEOUT", "10"))
large_timeout_sec <- as.numeric(Sys.getenv("RAPID_MA_LARGE_TIMEOUT", "240"))
memory_limit_gb <- as.numeric(Sys.getenv("RAPID_MA_MEMORY_LIMIT_GB", "Inf"))
matching_eval_max_common <- as.integer(Sys.getenv(
  "RAPID_MA_MATCHING_EVAL_MAX_COMMON", "100000"
))
if (any(!is.finite(c(timeout_sec, large_timeout_sec))) ||
    any(c(timeout_sec, large_timeout_sec) <= 0)) {
  stop("Scaling timeouts must be positive numbers of seconds.")
}
if (is.na(memory_limit_gb) || memory_limit_gb <= 0 ||
    is.na(matching_eval_max_common) || matching_eval_max_common < 0L) {
  stop("Memory and matching-evaluation limits must be positive.")
}

if (!"manifoldalign" %in% loadedNamespaces()) {
  if (requireNamespace("devtools", quietly = TRUE) && file.exists("DESCRIPTION")) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(manifoldalign)
  }
}

all_methods <- c(
  "rapid_ma_semisup", "lema", "ssma", "kema", "multiscale",
  "spectral_mnn", "token_ot"
)
if (identical(profile, "quick")) {
  common_sizes <- c(250L, 500L)
  rapid_sizes <- c(1000L)
  repetitions <- 1L
} else {
  common_sizes <- c(250L, 500L, 1000L, 2500L, 5000L, 10000L)
  rapid_sizes <- c(50000L, 100000L)
  repetitions <- 1:3
}
common_size_control <- Sys.getenv("RAPID_MA_COMMON_SIZES")
if (nzchar(common_size_control)) {
  common_sizes <- if (identical(tolower(common_size_control), "none")) {
    integer(0)
  } else {
    as.integer(strsplit(common_size_control, ",", fixed = TRUE)[[1L]])
  }
}
if (nzchar(Sys.getenv("RAPID_MA_RAPID_SIZES"))) {
  rapid_sizes <- as.integer(strsplit(
    Sys.getenv("RAPID_MA_RAPID_SIZES"), ",", fixed = TRUE
  )[[1L]])
}
if (nzchar(Sys.getenv("RAPID_MA_REPETITIONS"))) {
  repetitions <- seq_len(as.integer(Sys.getenv("RAPID_MA_REPETITIONS")))
}

common_tasks <- lapply(common_sizes, function(size) {
  size_repetitions <- if (size <= 1000L) repetitions else repetitions[[1L]]
  expand.grid(
    method = all_methods,
    n_common = size,
    repetition = size_repetitions,
    stringsAsFactors = FALSE
  )
})
tasks <- if (length(common_tasks)) do.call(rbind, common_tasks) else NULL
tasks <- rbind(
  tasks,
  expand.grid(
    method = "rapid_ma_semisup",
    n_common = rapid_sizes,
    repetition = 1L,
    stringsAsFactors = FALSE
  )
)
tasks$seed <- 6000L + tasks$repetition
tasks <- tasks[order(tasks$n_common, tasks$repetition, tasks$method), ]
rownames(tasks) <- NULL

sample_process_rss <- function(pid) {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  tryCatch(
    as.numeric(ps::ps_memory_info(ps::ps_handle(pid))[["rss"]]) / 1024^2,
    error = function(...) NA_real_
  )
}

failure_row <- function(task, elapsed, peak, error, status = "failed") {
  data.frame(
    method = task$method,
    status = status,
    evaluation_mode = NA_character_,
    classifier = NA_character_,
    runtime_sec = elapsed,
    retained_mb = NA_real_,
    endpoint_rss_mb = NA_real_,
    rss_delta_mb = NA_real_,
    oa = NA_real_, macro_f1 = NA_real_, miou = NA_real_,
    rare_recall = NA_real_, ece = NA_real_,
    hits1 = NA_real_, mrr = NA_real_, coverage = NA_real_,
    error = error,
    sampled_peak_rss_mb = peak,
    matching_evaluation_skipped = task$n_common > matching_eval_max_common,
    memory_limit_gb = memory_limit_gb,
    n_common = task$n_common,
    n_source = task$n_common,
    n_target = task$n_common + as.integer(ceiling(0.2 * task$n_common)),
    total_nodes = 2L * task$n_common + as.integer(ceiling(0.2 * task$n_common)),
    repetition = task$repetition,
    seed = task$seed,
    stringsAsFactors = FALSE
  )
}

run_task <- function(task) {
  task <- as.list(task)
  limit <- if (task$n_common >= 50000L) large_timeout_sec else timeout_sec
  start <- proc.time()[["elapsed"]]
  job <- parallel::mcparallel({
    fixture <- rapid_ma_benchmark_fixture(
      n_common = task$n_common,
      extra_target = as.integer(ceiling(0.2 * task$n_common)),
      labels_per_class = 4L,
      seed = task$seed
    )
    matching_evaluation_skipped <- task$n_common > matching_eval_max_common
    if (matching_evaluation_skipped) fixture$correspondence <- NULL
    benchmark_rapid_ma(
      fixture,
      methods = task$method,
      ncomp = 8L,
      seed = task$seed
    )
  }, mc.set.seed = FALSE, silent = TRUE)
  peak <- sample_process_rss(job$pid)
  value <- NULL
  last_progress <- 0
  repeat {
    collected <- parallel::mccollect(job, wait = FALSE)
    current <- sample_process_rss(job$pid)
    if (is.finite(current)) peak <- max(peak, current, na.rm = TRUE)
    if (!is.null(collected)) {
      value <- collected[[1L]]
      break
    }
    elapsed <- proc.time()[["elapsed"]] - start
    if (elapsed - last_progress >= 30) {
      message(
        "running method=", task$method,
        " n=", task$n_common,
        " elapsed_sec=", round(elapsed, 1),
        " peak_rss_mb=", round(peak, 1)
      )
      last_progress <- elapsed
    }
    if (is.finite(memory_limit_gb) && is.finite(peak) &&
        peak > memory_limit_gb * 1024) {
      parallel:::mckill(job, signal = 15L)
      parallel::mccollect(job, wait = TRUE, timeout = 1)
      return(failure_row(
        task, elapsed, peak,
        paste0("method exceeded isolated ", memory_limit_gb, " GB RSS ceiling"),
        status = "memory_limit"
      ))
    }
    if (elapsed >= limit) {
      parallel:::mckill(job, signal = 15L)
      parallel::mccollect(job, wait = TRUE, timeout = 1)
      return(failure_row(
        task, elapsed, peak,
        paste0("method exceeded isolated ", limit, " second timeout"),
        status = "timeout"
      ))
    }
    Sys.sleep(0.02)
  }
  elapsed <- proc.time()[["elapsed"]] - start
  if (is.finite(memory_limit_gb) && is.finite(peak) &&
      peak > memory_limit_gb * 1024) {
    return(failure_row(
      task, elapsed, peak,
      paste0("method exceeded isolated ", memory_limit_gb, " GB RSS ceiling"),
      status = "memory_limit"
    ))
  }
  if (inherits(value, "try-error") || !is.list(value) || is.null(value$results)) {
    return(failure_row(
      task, elapsed, peak,
      if (inherits(value, "try-error")) as.character(value) else
        "isolated method returned no benchmark result"
    ))
  }
  row <- value$results
  row$sampled_peak_rss_mb <- peak
  row$matching_evaluation_skipped <- task$n_common > matching_eval_max_common
  row$memory_limit_gb <- memory_limit_gb
  row$n_common <- task$n_common
  row$n_source <- task$n_common
  row$n_target <- task$n_common + as.integer(ceiling(0.2 * task$n_common))
  row$total_nodes <- row$n_source + row$n_target
  row$repetition <- task$repetition
  row$seed <- task$seed
  row
}

raw_path <- paste0(output_prefix, "-raw.csv")
summary_path <- paste0(output_prefix, "-summary.csv")
metadata_path <- paste0(output_prefix, "-metadata.txt")
dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)

rows <- vector("list", nrow(tasks))
for (index in seq_len(nrow(tasks))) {
  rows[[index]] <- suppressWarnings(run_task(tasks[index, , drop = FALSE]))
  raw_checkpoint <- do.call(rbind, rows[seq_len(index)])
  utils::write.csv(raw_checkpoint, raw_path, row.names = FALSE, na = "")
  message(
    "completed ", index, "/", nrow(tasks),
    " method=", tasks$method[[index]],
    " n=", tasks$n_common[[index]],
    " status=", rows[[index]]$status[[1L]]
  )
}
raw <- do.call(rbind, rows)
rownames(raw) <- NULL

groups <- split(
  seq_len(nrow(raw)),
  interaction(raw$method, raw$n_common, drop = TRUE, lex.order = TRUE)
)
summary <- do.call(rbind, lapply(groups, function(rows) {
  block <- raw[rows, , drop = FALSE]
  ok <- block$status == "ok"
  data.frame(
    method = block$method[[1L]],
    n_common = block$n_common[[1L]],
    total_nodes = block$total_nodes[[1L]],
    n_replicates = nrow(block),
    n_ok = sum(ok),
    n_timeout = sum(block$status == "timeout"),
    n_memory_limit = sum(block$status == "memory_limit"),
    n_failed = sum(!ok & !block$status %in% c("timeout", "memory_limit")),
    runtime_sec_mean = if (any(ok)) mean(block$runtime_sec[ok]) else NA_real_,
    runtime_sec_sd = if (sum(ok) > 1L) stats::sd(block$runtime_sec[ok]) else NA_real_,
    sampled_peak_rss_mb_mean = if (any(ok & is.finite(block$sampled_peak_rss_mb))) {
      mean(block$sampled_peak_rss_mb[ok], na.rm = TRUE)
    } else NA_real_,
    endpoint_rss_delta_mb_mean = if (any(ok & is.finite(block$rss_delta_mb))) {
      mean(block$rss_delta_mb[ok], na.rm = TRUE)
    } else NA_real_,
    retained_mb_mean = if (any(ok)) mean(block$retained_mb[ok]) else NA_real_,
    oa_mean = if (any(ok)) mean(block$oa[ok]) else NA_real_,
    macro_f1_mean = if (any(ok)) mean(block$macro_f1[ok]) else NA_real_,
    stringsAsFactors = FALSE
  )
}))
rownames(summary) <- NULL
summary <- summary[order(summary$n_common, summary$method), ]
utils::write.csv(summary, summary_path, row.names = FALSE, na = "")

source_paths <- c(
  "R/rapid_ma_benchmark.R", "R/rapid_ma.R", "R/rapid_ma_semisup.R",
  "inst/benchmarks/lema_reference.R", "inst/benchmarks/rapid_ma_scaling.R"
)
source_paths <- source_paths[file.exists(source_paths)]
memory_gb <- if (requireNamespace("ps", quietly = TRUE) &&
                 "ps_system_memory" %in% getNamespaceExports("ps")) {
  tryCatch(ps::ps_system_memory()[["total"]] / 1024^3, error = function(...) NA_real_)
} else NA_real_
metadata <- c(
  paste("generated_utc", format(Sys.time(), tz = "UTC", usetz = TRUE), sep = "="),
  paste("profile", profile, sep = "="),
  paste("common_sizes", paste(common_sizes, collapse = ","), sep = "="),
  paste("rapid_sizes", paste(rapid_sizes, collapse = ","), sep = "="),
  paste("repetitions", paste(repetitions, collapse = ","), sep = "="),
  paste("method_timeout_sec", timeout_sec, sep = "="),
  paste("large_timeout_sec", large_timeout_sec, sep = "="),
  paste("memory_limit_gb", memory_limit_gb, sep = "="),
  paste("matching_eval_max_common", matching_eval_max_common, sep = "="),
  "large_matching_note=correspondence evaluator skipped above matching_eval_max_common",
  "rss_sampling_interval_ms=20",
  "memory_note=absolute child RSS includes inherited R and package state",
  paste("machine", Sys.info()[["machine"]], sep = "="),
  paste("sysname", Sys.info()[["sysname"]], sep = "="),
  paste("release", Sys.info()[["release"]], sep = "="),
  paste("logical_cores", parallel::detectCores(logical = TRUE), sep = "="),
  paste("system_memory_gb", format(memory_gb, digits = 6), sep = "="),
  paste("R.version", R.version.string, sep = "="),
  paste("platform", R.version$platform, sep = "="),
  paste(names(tools::md5sum(source_paths)), tools::md5sum(source_paths), sep = " md5=")
)
writeLines(metadata, metadata_path)

print(summary, row.names = FALSE)
