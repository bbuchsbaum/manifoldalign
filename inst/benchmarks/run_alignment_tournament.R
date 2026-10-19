#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
option_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- arguments[startsWith(arguments, prefix)]
  if (!length(hit)) return(default)
  sub(prefix, "", hit[[length(hit)]], fixed = TRUE)
}

script_argument <- commandArgs()[startsWith(commandArgs(), "--file=")]
script_path <- if (length(script_argument)) {
  normalizePath(sub("--file=", "", script_argument[[1L]], fixed = TRUE))
} else {
  NA_character_
}
repo_root <- if (file.exists("DESCRIPTION")) {
  normalizePath(".")
} else if (!is.na(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  stop("Run this script from the manifoldalign repository root.", call. = FALSE)
}
setwd(repo_root)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("The tournament runner requires the suggested package `pkgload`.", call. = FALSE)
}
if (!requireNamespace("ps", quietly = TRUE)) {
  stop("The tournament runner requires the suggested package `ps` for peak RSS.", call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE)

profile <- match.arg(
  option_value("profile", "smoke"), c("smoke", "oracle", "final")
)
default_output <- if (identical(profile, "final")) {
  file.path(repo_root, "docs", "evidence", "alignment-tournament-v2")
} else {
  file.path(tempdir(), paste0("alignment-tournament-v2-", profile))
}
output_dir <- normalizePath(
  option_value("output", default_output), mustWork = FALSE
)
resume <- !tolower(option_value("resume", "true")) %in% c("false", "0", "no")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(output_dir, "run.log")

log_event <- function(...) {
  message <- paste0(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", paste0(..., collapse = ""))
  cat(message, "\n", file = log_file, append = TRUE)
  cat(message, "\n")
  flush.console()
}

manifest <- manifoldalign:::graph_alignment_benchmark_manifest_v2("final")
case_override <- option_value("cases")
if (!is.null(case_override)) {
  case_ids <- strsplit(case_override, ",", fixed = TRUE)[[1L]]
} else if (identical(profile, "smoke")) {
  case_ids <- c(
    "final-oracle-sbm-01",
    "final-small-accuracy-sbm-01"
  )
} else if (identical(profile, "oracle")) {
  case_ids <- manifest$case_id[manifest$tier == "oracle"]
} else {
  case_ids <- manifest$case_id
}
unknown_cases <- setdiff(case_ids, manifest$case_id)
if (length(unknown_cases)) {
  stop("Unknown final case(s): ", paste(unknown_cases, collapse = ", "), call. = FALSE)
}
manifest <- manifest[match(case_ids, manifest$case_id), , drop = FALSE]
attr(manifest, "manifest_hash") <- manifoldalign:::.graph_v2_manifest_hashes[["final"]]

method_override <- option_value("methods")
methods <- if (is.null(method_override)) {
  manifoldalign:::.alignment_tournament_methods_v2()
} else {
  strsplit(method_override, ",", fixed = TRUE)[[1L]]
}
unknown_methods <- setdiff(methods, manifoldalign:::.alignment_tournament_methods_v2())
if (length(unknown_methods)) {
  stop("Unknown tournament method(s): ", paste(unknown_methods, collapse = ", "), call. = FALSE)
}

sample_worker_rss_mb <- function(pid) {
  root <- tryCatch(ps::ps_handle(pid), error = function(error) NULL)
  if (is.null(root)) return(NA_real_)
  handles <- c(
    list(root),
    tryCatch(ps::ps_children(root, recursive = TRUE), error = function(error) list())
  )
  rss <- vapply(handles, function(handle) {
    info <- tryCatch(ps::ps_memory_info(handle), error = function(error) NULL)
    if (is.null(info)) return(NA_real_)
    as.numeric(info[["rss"]]) / 1024^2
  }, numeric(1))
  if (all(!is.finite(rss))) NA_real_ else sum(rss[is.finite(rss)])
}

run_external <- function(case, method) {
  applicability <- manifoldalign:::.alignment_tournament_applicability_v2(
    method, case
  )
  if (!isTRUE(applicability$intended_to_run)) {
    return(manifoldalign:::.alignment_tournament_base_row_v2(
      case, method, applicability
    ))
  }
  if (.Platform$OS.type == "windows") {
    row <- manifoldalign:::.alignment_tournament_run_case_v2(case, method)
    row$process_elapsed_sec <- row$runtime_sec
    row$memory_measurement <- "unavailable_on_windows_fallback"
    row$peak_rss_mb <- NA_real_
    return(row)
  }

  started <- Sys.time()
  job <- parallel::mcparallel(
    manifoldalign:::.alignment_tournament_run_case_v2(case, method),
    mc.set.seed = FALSE,
    silent = TRUE
  )
  pid <- as.integer(job$pid)
  budget <- as.numeric(applicability$contract$wall_clock_budget_sec)
  peak <- NA_real_
  collected <- NULL
  timed_out <- FALSE
  repeat {
    rss <- sample_worker_rss_mb(pid)
    if (is.finite(rss)) peak <- if (is.finite(peak)) max(peak, rss) else rss
    collected <- suppressWarnings(parallel::mccollect(job, wait = FALSE))
    if (!is.null(collected)) break
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (is.finite(budget) && elapsed >= budget) {
      timed_out <- TRUE
      try(tools::pskill(pid, signal = 15L), silent = TRUE)
      Sys.sleep(0.2)
      try(tools::pskill(pid, signal = 9L), silent = TRUE)
      suppressWarnings(try(parallel::mccollect(job, wait = FALSE), silent = TRUE))
      break
    }
    Sys.sleep(0.05)
  }
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  if (timed_out) {
    row <- manifoldalign:::.alignment_tournament_failure_row_v2(
      case, method, "timeout", "wall_clock_budget_exceeded",
      paste0("External worker exceeded the frozen ", budget, " second budget."),
      runtime_sec = elapsed, peak_rss_mb = peak,
      memory_measurement = "external_process_peak_rss"
    )
  } else if (!length(collected) || inherits(collected[[1L]], "try-error")) {
    detail <- if (length(collected)) {
      as.character(collected[[1L]])
    } else {
      "External worker exited without returning a result."
    }
    row <- manifoldalign:::.alignment_tournament_failure_row_v2(
      case, method, "crash", "external_worker_crash", detail,
      runtime_sec = elapsed, peak_rss_mb = peak,
      memory_measurement = "external_process_peak_rss"
    )
  } else {
    row <- collected[[1L]]
    row$process_elapsed_sec <- elapsed
    row$worker_pid <- pid
    row$peak_rss_mb <- peak
    row$memory_measurement <- "external_process_peak_rss"
  }
  row$worker_pid <- pid
  row$process_elapsed_sec <- elapsed
  row
}

checkpoint_path <- function(case_id, method) {
  file.path(checkpoint_dir, paste0(case_id, "__", method, ".rds"))
}
checkpoint_is_current <- function(row, case, method) {
  if (!is.data.frame(row) || nrow(row) != 1L) return(FALSE)
  expected <- manifoldalign:::.alignment_tournament_base_row_v2(case, method)
  identical(row$tournament_version, expected$tournament_version) &&
    identical(row$manifest_hash, expected$manifest_hash) &&
    identical(row$case_hash, expected$case_hash) &&
    identical(row$parameter_set_hash, expected$parameter_set_hash)
}

log_event(
  "START profile=", profile, " cases=", length(case_ids),
  " methods=", length(methods), " resume=", resume,
  " output=", output_dir
)

rows <- list()
controls <- list()
pointer <- 0L
control_pointer <- 0L
for (case_index in seq_along(case_ids)) {
  case_id <- case_ids[[case_index]]
  log_event("CASE ", case_index, "/", length(case_ids), " ", case_id)
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(case_id, "final")
  for (method_index in seq_along(methods)) {
    method <- methods[[method_index]]
    path <- checkpoint_path(case_id, method)
    row <- NULL
    if (resume && file.exists(path)) {
      candidate <- tryCatch(readRDS(path), error = function(error) NULL)
      if (checkpoint_is_current(candidate, case, method)) row <- candidate
    }
    if (is.null(row)) {
      log_event(
        "RUN ", case_index, "/", length(case_ids), " ",
        method_index, "/", length(methods), " ", case_id, " ", method
      )
      row <- run_external(case, method)
      saveRDS(row, path, version = 3)
    } else {
      log_event("RESUME ", case_id, " ", method, " status=", row$status)
    }
    pointer <- pointer + 1L
    rows[[pointer]] <- row
    log_event(
      "RESULT ", case_id, " ", method, " status=", row$status,
      " top1=", ifelse(is.finite(row$top1_accuracy), format(row$top1_accuracy, digits = 4), "NA"),
      " elapsed=", format(row$process_elapsed_sec, digits = 4),
      " peak_rss_mb=", ifelse(is.finite(row$peak_rss_mb), format(row$peak_rss_mb, digits = 5), "NA")
    )
  }

  control_started <- proc.time()[[3L]]
  control <- manifoldalign:::benchmark_graph_alignment_controls_v2(
    case, purpose = "report", evaluation_role = "final"
  )
  control_elapsed <- proc.time()[[3L]] - control_started
  control$case_id <- case$case_id
  control$case_hash <- case$provenance$case_hash
  control$tier <- case$tier
  control$family <- case$config$family
  control$overlap_regime <- if (
    any(case$truth$source_private_mask) || any(case$truth$target_private_mask)
  ) "partial_overlap" else "closed_world"
  control$generation_seed <- case$config$seed
  control$batch_runtime_sec <- control_elapsed
  control_pointer <- control_pointer + 1L
  controls[[control_pointer]] <- control
}

results <- do.call(rbind, rows)
rownames(results) <- NULL
control_results <- do.call(rbind, controls)
rownames(control_results) <- NULL

summaries <- manifoldalign:::.alignment_tournament_summarize_v2(results)
problems <- manifoldalign:::.alignment_tournament_validate_results_v2(
  results, manifest, methods, require_external_memory = TRUE
)
ledger <- manifoldalign:::.alignment_tournament_parameter_ledger_v2(
  manifest, methods
)

utils::write.csv(results, file.path(output_dir, "raw_results.csv"), row.names = FALSE, na = "")
utils::write.csv(control_results, file.path(output_dir, "controls.csv"), row.names = FALSE, na = "")
utils::write.csv(manifest, file.path(output_dir, "manifest.csv"), row.names = FALSE, na = "")
utils::write.csv(ledger, file.path(output_dir, "tuning_ledger.csv"), row.names = FALSE, na = "")
utils::write.csv(summaries$by_regime, file.path(output_dir, "summary_by_regime.csv"), row.names = FALSE, na = "")
utils::write.csv(summaries$paired_vs_rapid_ma, file.path(output_dir, "paired_vs_rapid_ma.csv"), row.names = FALSE, na = "")
utils::write.csv(summaries$pareto_by_regime, file.path(output_dir, "pareto_by_regime.csv"), row.names = FALSE, na = "")
utils::write.csv(summaries$scaling, file.path(output_dir, "scaling_exponents.csv"), row.names = FALSE, na = "")
utils::write.csv(problems, file.path(output_dir, "validation_problems.csv"), row.names = FALSE, na = "")

session_lines <- capture.output(sessionInfo())
writeLines(session_lines, file.path(output_dir, "sessionInfo.txt"))
git_commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE),
  error = function(error) NA_character_
)
git_status <- tryCatch(
  system2("git", c("status", "--short"), stdout = TRUE, stderr = TRUE),
  error = function(error) NA_character_
)
writeLines(
  c(
    paste0("tournament_version: ", manifoldalign:::.alignment_tournament_version_v2),
    paste0("profile: ", profile),
    paste0("started_or_resumed: ", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
    paste0("manifest_hash: ", manifoldalign:::.graph_v2_manifest_hashes[["final"]]),
    paste0("git_commit: ", paste(git_commit, collapse = " ")),
    paste0("case_count: ", length(case_ids)),
    paste0("method_count: ", length(methods)),
    paste0("result_count: ", nrow(results)),
    paste0("validation_problem_count: ", nrow(problems)),
    "regeneration_command:",
    paste0(
      "Rscript inst/benchmarks/run_alignment_tournament.R --profile=", profile,
      " --output=", output_dir, " --resume=true"
    ),
    "git_status:",
    git_status
  ),
  file.path(output_dir, "run_receipt.txt")
)

plot_summary <- summaries$by_regime[
  summaries$by_regime$n_intended > 0L &
    is.finite(summaries$by_regime$intent_top1) &
    is.finite(summaries$by_regime$median_runtime_sec),
  , drop = FALSE
]
if (nrow(plot_summary)) {
  grDevices::png(
    file.path(output_dir, "accuracy_runtime.png"),
    width = 1400, height = 900, res = 144
  )
  graphics::plot(
    pmax(plot_summary$median_runtime_sec, 1e-4),
    plot_summary$intent_top1,
    log = "x", pch = 19,
    xlab = "Median external runtime (seconds, log scale)",
    ylab = "Intent-to-run Hits@1",
    main = paste0("Alignment tournament: ", profile)
  )
  graphics::text(
    pmax(plot_summary$median_runtime_sec, 1e-4),
    plot_summary$intent_top1,
    labels = plot_summary$method_id, pos = 3, cex = 0.65
  )
  grDevices::dev.off()
}

status_table <- as.data.frame(table(results$status), stringsAsFactors = FALSE)
names(status_table) <- c("status", "count")
readme <- c(
  "# Alignment tournament v2 evidence",
  "",
  paste0("Profile: `", profile, "`; final manifest MD5: `", manifoldalign:::.graph_v2_manifest_hashes[["final"]], "`."),
  "",
  paste0("This receipt contains ", nrow(results), " typed method-case rows across ", length(case_ids), " cases and ", length(methods), " methods."),
  "Failures and timeouts remain in intent-to-run denominators; unsupported tasks remain visible but are not scored as algorithm failures.",
  "Matching uses native top-1 followed by a declared transport or common-embedding ranking. Classification is a frozen source-label downstream probe, not the method's private classifier.",
  "The checked interpretation and lifecycle implications are in `INTERPRETATION.md`.",
  "Checksums cover published receipt files; resumable binary checkpoints and the debug log remain local operational state.",
  "",
  "## Status counts",
  "",
  "| Status | Count |",
  "|---|---:|",
  paste0("| ", status_table$status, " | ", status_table$count, " |"),
  "",
  "## Regeneration",
  "",
  "```sh",
  paste0(
    "Rscript inst/benchmarks/run_alignment_tournament.R --profile=", profile,
    " --output=", output_dir, " --resume=true"
  ),
  "```",
  "",
  "See `validation_problems.csv`; an official receipt is valid only when it has zero rows."
)
writeLines(readme, file.path(output_dir, "README.md"))

log_event(
  "DONE rows=", nrow(results), " validation_problems=", nrow(problems),
  " output=", output_dir
)
artifact_files <- list.files(
  output_dir, full.names = TRUE, recursive = TRUE
)
artifact_files <- artifact_files[!file.info(artifact_files)$isdir]
artifact_relative <- substring(artifact_files, nchar(output_dir) + 2L)
publishable <- basename(artifact_files) != "checksums.csv" &
  artifact_relative != "run.log" &
  !startsWith(artifact_relative, "checkpoints/")
artifact_files <- artifact_files[publishable]
artifact_relative <- artifact_relative[publishable]
checksums <- data.frame(
  file = artifact_relative,
  bytes = as.numeric(file.info(artifact_files)$size),
  md5 = unname(tools::md5sum(artifact_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(checksums, file.path(output_dir, "checksums.csv"), row.names = FALSE)

if (nrow(problems)) quit(status = 1L)
