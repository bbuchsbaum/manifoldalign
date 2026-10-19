#!/usr/bin/env Rscript

# Run governed real-data controls and emit auditable evidence receipts. The
# packaged So2Sat fixture checks schemas and leakage boundaries only. Passing
# --full-so2sat derives and scores the complete frozen validation/testing split.

find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(
          path, "inst", "benchmarks", "real_alignment_manifest.csv"
        ))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the manifoldalign package root.", call. = FALSE)
    }
    path <- parent
  }
}

parse_args <- function(args, root) {
  output_dir <- file.path(root, "artifacts", "alignment-real")
  data_dir <- Sys.getenv("MANIFOLDALIGN_REAL_DATA_DIR", unset = "")
  full_so2sat <- FALSE
  rebuild <- FALSE
  strict_data <- FALSE
  for (arg in args) {
    if (startsWith(arg, "--output-dir=")) {
      output_dir <- sub("^--output-dir=", "", arg)
    } else if (startsWith(arg, "--data-dir=")) {
      data_dir <- sub("^--data-dir=", "", arg)
    } else if (identical(arg, "--full-so2sat")) {
      full_so2sat <- TRUE
    } else if (identical(arg, "--rebuild-so2sat-features")) {
      rebuild <- TRUE
    } else if (identical(arg, "--strict-data")) {
      strict_data <- TRUE
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", output_dir)) {
    output_dir <- file.path(root, output_dir)
  }
  list(
    output_dir = normalizePath(output_dir, mustWork = FALSE),
    data_dir = if (nzchar(data_dir)) {
      normalizePath(data_dir, mustWork = FALSE)
    } else {
      ""
    },
    full_so2sat = full_so2sat,
    rebuild = rebuild,
    strict_data = strict_data
  )
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

bind_rows <- function(rows) {
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(NULL)
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (name in missing) row[[name]] <- NA
    row[, columns, drop = FALSE]
  })
  do.call(rbind, rows)
}

network_metadata <- function(case) {
  d1 <- case$domains[[1L]]
  d2 <- case$domains[[2L]]
  data.frame(
    protocol_version = case$protocol_version,
    task_id = case$task_id,
    dataset_id = case$dataset_id,
    evidence_scope = "full_packaged_network",
    split_id = case$split_id,
    n_source = nrow(d1$X), n_target = nrow(d2$X),
    n_shared = sum(!is.na(case$truth$map12)),
    n_source_private = sum(d1$private_mask),
    n_target_private = sum(d2$private_mask),
    n_anchors = nrow(case$anchors),
    n_source_edges = nrow(d1$relations),
    n_target_edges = nrow(d2$relations),
    source_sha256 = case$governance$source_sha256,
    stringsAsFactors = FALSE
  )
}

so2sat_metadata <- function(case) {
  data.frame(
    protocol_version = case$protocol_version,
    task_id = case$task_id,
    dataset_id = case$dataset_id,
    evidence_scope = case$governance$evidence_scope,
    split_id = case$governance$split_id,
    n_development = case$governance$n_development,
    n_final = case$governance$n_final,
    n_features_per_sensor = ncol(case$fit$source_X),
    feature_artifact_sha256 = case$governance$feature_artifact_sha256,
    final_labels_available_to_fit =
      case$governance$final_labels_available_to_fit,
    stringsAsFactors = FALSE
  )
}

verify_so2sat_inputs <- function(data_dir, artifacts, sha256_file) {
  required <- c("validation", "validation_geo", "testing", "testing_geo")
  selected <- artifacts[
    artifacts$dataset_id == "so2sat_lcz42" &
      artifacts$artifact_id %in% required,
    , drop = FALSE
  ]
  selected <- selected[match(required, selected$artifact_id), , drop = FALSE]
  rows <- vector("list", nrow(selected))
  for (i in seq_len(nrow(selected))) {
    path <- file.path(data_dir, selected$relative_path[[i]])
    status <- "ok"
    code <- "verified"
    detail <- NA_character_
    observed_sha <- NA_character_
    observed_bytes <- if (file.exists(path)) file.info(path)$size else NA_real_
    if (!file.exists(path)) {
      status <- "dataset_unavailable"
      code <- "artifact_missing"
      detail <- paste0("Missing ", selected$relative_path[[i]], ".")
    } else if (!identical(
      as.numeric(observed_bytes), as.numeric(selected$bytes[[i]])
    )) {
      status <- "dataset_invalid"
      code <- "byte_count_mismatch"
      detail <- paste0(
        "Expected ", selected$bytes[[i]], " bytes; observed ",
        observed_bytes, "."
      )
    } else {
      observed_sha <- sha256_file(path)
      if (!identical(observed_sha, selected$sha256[[i]])) {
        status <- "dataset_invalid"
        code <- "checksum_mismatch"
        detail <- paste0(
          "Expected SHA-256 ", selected$sha256[[i]],
          "; observed ", observed_sha, "."
        )
      }
    }
    rows[[i]] <- data.frame(
      protocol_version = "1.0.0",
      task_id = "so2sat_lcz42_sar_optical_culture10",
      stage = "input_validation", artifact = selected$artifact_id[[i]],
      status = status, status_code = code, status_detail = detail,
      failure_class = if (identical(status, "ok")) NA_character_ else "dataset",
      expected_bytes = selected$bytes[[i]], observed_bytes = observed_bytes,
      expected_sha256 = selected$sha256[[i]], observed_sha256 = observed_sha,
      path = path, stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

receipt <- function(paths, sha256_file) {
  paths <- paths[file.exists(paths)]
  data.frame(
    file = basename(paths),
    bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, sha256_file, character(1)),
    stringsAsFactors = FALSE
  )
}

main <- function() {
  root <- find_package_root()
  args <- parse_args(commandArgs(trailingOnly = TRUE), root)
  dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The real benchmark runner requires pkgload.", call. = FALSE)
  }
  pkgload::load_all(root, quiet = TRUE)
  ns <- asNamespace("manifoldalign")
  fn <- function(name) get(name, envir = ns, inherits = FALSE)
  registry <- fn("real_alignment_benchmark_registry")()
  artifacts <- fn("real_alignment_artifact_registry")()
  sha256_file <- fn(".real_sha256_file")
  status_row <- fn("real_alignment_status_row")

  registry_path <- file.path(args$output_dir, "real-benchmark-registry.csv")
  artifact_path <- file.path(args$output_dir, "real-artifact-registry.csv")
  write_csv(registry, registry_path)
  write_csv(artifacts, artifact_path)

  statuses <- list()
  metadata <- list()
  network_controls <- list()
  network_tasks <- registry$task_id[
    registry$task_type == "attributed_network_correspondence"
  ]
  for (task_id in network_tasks) {
    case <- tryCatch(
      fn("real_alignment_network_case")(task_id),
      error = function(e) e
    )
    if (inherits(case, "error")) {
      statuses[[length(statuses) + 1L]] <- status_row(
        task_id, "case_construction", "dataset_invalid", "fixture_invalid",
        conditionMessage(case)
      )
      next
    }
    metadata[[length(metadata) + 1L]] <- network_metadata(case)
    controls <- tryCatch(
      fn("benchmark_real_network_controls")(case),
      error = function(e) e
    )
    if (inherits(controls, "error")) {
      statuses[[length(statuses) + 1L]] <- status_row(
        task_id, "controls", "algorithm_error", "control_failure",
        conditionMessage(controls)
      )
    } else {
      controls$evidence_scope <- "full_packaged_network"
      network_controls[[length(network_controls) + 1L]] <- controls
      statuses[[length(statuses) + 1L]] <- status_row(
        task_id, "controls", "ok", "complete"
      )
    }
  }

  ci_case <- tryCatch(fn("real_so2sat_ci_case")(), error = function(e) e)
  so2sat_controls <- list()
  if (inherits(ci_case, "error")) {
    statuses[[length(statuses) + 1L]] <- status_row(
      "so2sat_lcz42_sar_optical_culture10", "ci_fixture",
      "dataset_invalid", "fixture_invalid", conditionMessage(ci_case)
    )
  } else {
    metadata[[length(metadata) + 1L]] <- so2sat_metadata(ci_case)
    controls <- tryCatch(
      fn("benchmark_real_so2sat_controls")(ci_case),
      error = function(e) e
    )
    if (inherits(controls, "error")) {
      statuses[[length(statuses) + 1L]] <- status_row(
        ci_case$task_id, "ci_controls", "algorithm_error",
        "control_failure", conditionMessage(controls)
      )
    } else {
      controls$evidence_scope <- "ci_schema_only"
      so2sat_controls[[length(so2sat_controls) + 1L]] <- controls
      statuses[[length(statuses) + 1L]] <- status_row(
        ci_case$task_id, "ci_controls", "ok", "schema_only_complete"
      )
    }
  }

  full_input_status <- NULL
  if (isTRUE(args$full_so2sat)) {
    if (!nzchar(args$data_dir) || !dir.exists(args$data_dir)) {
      statuses[[length(statuses) + 1L]] <- status_row(
        "so2sat_lcz42_sar_optical_culture10", "full_inputs",
        "dataset_unavailable", "data_directory_missing",
        "Provide --data-dir or MANIFOLDALIGN_REAL_DATA_DIR."
      )
    } else {
      full_input_status <- verify_so2sat_inputs(
        args$data_dir, artifacts, sha256_file
      )
      invalid <- full_input_status$status != "ok"
      if (any(invalid)) {
        statuses[[length(statuses) + 1L]] <- status_row(
          "so2sat_lcz42_sar_optical_culture10", "full_inputs",
          if (any(full_input_status$status == "dataset_invalid")) {
            "dataset_invalid"
          } else {
            "dataset_unavailable"
          },
          "artifact_validation_failed",
          paste(full_input_status$status_detail[invalid], collapse = " ")
        )
      } else if (!requireNamespace("hdf5r", quietly = TRUE)) {
        statuses[[length(statuses) + 1L]] <- status_row(
          "so2sat_lcz42_sar_optical_culture10", "feature_extraction",
          "dataset_unavailable", "hdf5r_missing",
          "Install hdf5r to materialize full So2Sat features."
        )
      } else {
        feature_path <- file.path(
          args$output_dir, "so2sat_lcz42_full_features.csv"
        )
        extraction_error <- NULL
        if (isTRUE(args$rebuild) || !file.exists(feature_path)) {
          rscript <- file.path(R.home("bin"), "Rscript")
          output <- tryCatch(
            system2(rscript, c(
              file.path(root, "inst", "benchmarks",
                "extract_so2sat_ci_fixture.R"
              ),
              paste0("--data-dir=", args$data_dir),
              paste0("--output=", feature_path),
              "--per-class=all"
            ), stdout = TRUE, stderr = TRUE),
            error = function(e) e
          )
          if (inherits(output, "error") ||
              !identical(attr(output, "status") %||% 0L, 0L)) {
            extraction_error <- if (inherits(output, "error")) {
              conditionMessage(output)
            } else {
              paste(output, collapse = "\n")
            }
          }
          if (!inherits(output, "error")) {
            writeLines(output, file.path(
              args$output_dir, "so2sat-feature-extraction.log"
            ))
          }
        }
        if (!is.null(extraction_error)) {
          statuses[[length(statuses) + 1L]] <- status_row(
            "so2sat_lcz42_sar_optical_culture10", "feature_extraction",
            "dataset_invalid", "feature_extraction_failed", extraction_error
          )
        } else {
          expected_feature <- artifacts[
            artifacts$dataset_id == "so2sat_lcz42" &
              artifacts$artifact_id == "full_feature_receipt",
            , drop = FALSE
          ]
          observed_bytes <- if (file.exists(feature_path)) {
            as.numeric(file.info(feature_path)$size)
          } else {
            NA_real_
          }
          observed_sha <- if (file.exists(feature_path)) {
            sha256_file(feature_path)
          } else {
            NA_character_
          }
          receipt_valid <- nrow(expected_feature) == 1L &&
            identical(observed_bytes, as.numeric(expected_feature$bytes[[1L]])) &&
            identical(observed_sha, expected_feature$sha256[[1L]])
          if (!receipt_valid) {
            statuses[[length(statuses) + 1L]] <- status_row(
              "so2sat_lcz42_sar_optical_culture10", "feature_extraction",
              "dataset_invalid", "feature_receipt_mismatch",
              paste0(
                "Generated feature receipt differs: ", observed_bytes,
                " bytes / ", observed_sha, "."
              ),
              basename(feature_path)
            )
          } else {
            full_case <- tryCatch(
              fn("real_so2sat_full_case")(feature_path),
              error = function(e) e
            )
            if (inherits(full_case, "error")) {
              statuses[[length(statuses) + 1L]] <- status_row(
                "so2sat_lcz42_sar_optical_culture10", "full_case",
                "dataset_invalid", "feature_artifact_invalid",
                conditionMessage(full_case), basename(feature_path)
              )
            } else {
              metadata[[length(metadata) + 1L]] <- so2sat_metadata(full_case)
              controls <- tryCatch(
                fn("benchmark_real_so2sat_controls")(full_case),
                error = function(e) e
              )
              if (inherits(controls, "error")) {
                statuses[[length(statuses) + 1L]] <- status_row(
                  full_case$task_id, "full_controls", "algorithm_error",
                  "control_failure", conditionMessage(controls)
                )
              } else {
                controls$evidence_scope <- "full_frozen_split"
                so2sat_controls[[length(so2sat_controls) + 1L]] <- controls
                statuses[[length(statuses) + 1L]] <- status_row(
                  full_case$task_id, "full_controls", "ok", "complete",
                  artifact = basename(feature_path)
                )
              }
            }
          }
        }
      }
    }
  } else {
    statuses[[length(statuses) + 1L]] <- status_row(
      "so2sat_lcz42_sar_optical_culture10", "full_inputs",
      "dataset_unavailable", "full_run_not_requested",
      "Pass --full-so2sat for release-grade accuracy evidence."
    )
  }

  output_files <- c(registry_path, artifact_path)
  tables <- list(
    "real-dataset-status.csv" = bind_rows(statuses),
    "real-case-metadata.csv" = bind_rows(metadata),
    "real-network-controls.csv" = bind_rows(network_controls),
    "real-so2sat-controls.csv" = bind_rows(so2sat_controls),
    "real-so2sat-input-validation.csv" = full_input_status
  )
  for (name in names(tables)) {
    if (!is.null(tables[[name]])) {
      path <- file.path(args$output_dir, name)
      write_csv(tables[[name]], path)
      output_files <- c(output_files, path)
    }
  }
  session_path <- file.path(args$output_dir, "session-info.txt")
  writeLines(c(
    paste0("run_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("package_root: ", root),
    paste0("full_so2sat_requested: ", args$full_so2sat),
    "",
    capture.output(sessionInfo())
  ), session_path)
  output_files <- c(output_files, session_path)
  write_csv(
    receipt(output_files, sha256_file),
    file.path(args$output_dir, "run-receipt.csv")
  )
  message("Real alignment evidence written to ", args$output_dir)
  status_table <- bind_rows(statuses)
  if (isTRUE(args$strict_data) &&
      any(status_table$status %in% c("dataset_unavailable", "dataset_invalid"))) {
    return(2L)
  }
  if (any(status_table$status == "algorithm_error")) 1L else 0L
}

`%||%` <- function(x, y) if (is.null(x)) y else x

status <- main()
quit(save = "no", status = status, runLast = FALSE)
