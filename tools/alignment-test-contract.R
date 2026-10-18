# Base-R helpers for the enforced alignment test contract -----------------

.alignment_ci_empty_problems <- function() {
  data.frame(
    severity = character(),
    code = character(),
    method = character(),
    path = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
}

.alignment_ci_problem <- function(
  code,
  message,
  method = NA_character_,
  path = NA_character_,
  severity = "error"
) {
  data.frame(
    severity = severity,
    code = code,
    method = method,
    path = path,
    message = message,
    stringsAsFactors = FALSE
  )
}

.alignment_ci_unquote <- function(x) {
  x <- trimws(x)
  if (
    nchar(x) >= 2L &&
      substr(x, 1L, 1L) %in% c("\"", "'") &&
      substr(x, nchar(x), nchar(x)) == substr(x, 1L, 1L)
  ) {
    x <- substr(x, 2L, nchar(x) - 1L)
  }
  x
}

.alignment_ci_inline_value <- function(x) {
  x <- trimws(x)
  if (identical(x, "[]")) {
    return(character())
  }
  if (grepl("^\\[.*\\]$", x)) {
    x <- substr(x, 2L, nchar(x) - 1L)
    if (!nzchar(trimws(x))) {
      return(character())
    }
    return(vapply(
      strsplit(x, ",", fixed = TRUE)[[1L]],
      .alignment_ci_unquote,
      character(1)
    ))
  }
  if (grepl("^[0-9]+$", x)) {
    return(as.integer(x))
  }
  .alignment_ci_unquote(x)
}

#' Read the alignment CI contract from ALIGNMENT_TEST_MATRIX.yml
#'
#' This intentionally parses only the small, flat `ci_contract` schema. It
#' avoids making the test contract itself dependent on an optional YAML
#' package. The broader fixture inventory remains ordinary YAML for people and
#' external tooling.
#'
#' @param path Path to `ALIGNMENT_TEST_MATRIX.yml`.
#' @return A list containing contract metadata and method records.
#' @keywords internal
read_alignment_test_contract <- function(path) {
  lines <- readLines(path, warn = FALSE)
  start <- match("ci_contract:", lines)
  if (is.na(start)) {
    stop("ALIGNMENT_TEST_MATRIX.yml has no ci_contract section.", call. = FALSE)
  }

  header <- list()
  records <- list()
  current <- NULL
  in_records <- FALSE

  flush_record <- function() {
    if (!is.null(current)) {
      records[[length(records) + 1L]] <<- current
    }
  }

  for (line in lines[seq.int(start + 1L, length(lines))]) {
    if (grepl("^[^[:space:]#]", line)) {
      break
    }
    if (grepl("^  method_records:[[:space:]]*$", line)) {
      in_records <- TRUE
      next
    }
    if (in_records && grepl("^    - method:[[:space:]]*", line)) {
      flush_record()
      value <- sub("^    - method:[[:space:]]*", "", line)
      current <- list(method = .alignment_ci_inline_value(value))
      next
    }
    if (in_records && grepl("^      [a-z][a-z0-9_]*:[[:space:]]*", line)) {
      key <- sub("^      ([a-z][a-z0-9_]*):.*$", "\\1", line)
      value <- sub("^      [a-z][a-z0-9_]*:[[:space:]]*", "", line)
      current[[key]] <- .alignment_ci_inline_value(value)
      next
    }
    if (!in_records && grepl("^  [a-z][a-z0-9_]*:[[:space:]]*", line)) {
      key <- sub("^  ([a-z][a-z0-9_]*):.*$", "\\1", line)
      value <- sub("^  [a-z][a-z0-9_]*:[[:space:]]*", "", line)
      header[[key]] <- .alignment_ci_inline_value(value)
    }
  }
  flush_record()

  names(records) <- vapply(records, function(x) x$method, character(1))
  structure(
    c(header, list(method_records = records, source = normalizePath(path))),
    class = "alignment_test_contract"
  )
}

.alignment_ci_test_titles <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hits <- grep("test_that\\([[:space:]]*\"[^\"]+\"", lines, value = TRUE)
  sub('^.*test_that\\([[:space:]]*\"([^\"]+)\".*$', "\\1", hits)
}

.alignment_ci_mandatory_parts <- function(x) {
  split <- strsplit(x, "::", fixed = TRUE)[[1L]]
  if (length(split) != 2L || !all(nzchar(split))) {
    return(list(path = NA_character_, test = NA_character_))
  }
  list(path = split[[1L]], test = split[[2L]])
}

.alignment_ci_suggests <- function(root) {
  desc <- read.dcf(file.path(root, "DESCRIPTION"), fields = "Suggests")
  value <- unname(desc[1L, 1L])
  if (is.na(value) || !nzchar(value)) {
    return(character())
  }
  trimws(sub("[[:space:]]*\\(.*$", "", strsplit(value, ",")[[1L]]))
}

.alignment_ci_mandatory_block <- function(path, title) {
  lines <- readLines(path, warn = FALSE)
  escaped <- gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", title)
  starts <- grep(
    paste0('test_that\\([[:space:]]*\"', escaped, '\"'),
    lines
  )
  if (length(starts) != 1L) {
    return(character())
  }
  later <- grep("test_that\\([[:space:]]*\"", lines)
  later <- later[later > starts]
  end <- if (length(later)) later[[1L]] - 1L else length(lines)
  lines[seq.int(starts, end)]
}

.alignment_ci_tautologies <- function(root, pattern_ids) {
  known <- c(
    "nonnegative_correct_count",
    "nonnegative_length",
    "lower_bound_only_accuracy"
  )
  unknown <- setdiff(pattern_ids, known)
  if (length(unknown)) {
    stop(
      "Unknown forbidden_tautologies IDs: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  files <- list.files(
    file.path(root, "tests", "testthat"),
    pattern = "^test.*\\.[rR]$",
    full.names = TRUE
  )
  out <- list()
  for (path in files) {
    lines <- readLines(path, warn = FALSE)
    for (i in seq_along(lines)) {
      line <- lines[[i]]
      hit <- character()
      if (
        "nonnegative_correct_count" %in% pattern_ids &&
          grepl(
            "expect_(true|gte)\\([^#]*(n_correct|correct_count)[^#]*(>=[[:space:]]*0|,[[:space:]]*0[[:space:]]*\\))",
            line
          )
      ) {
        hit <- c(hit, "nonnegative_correct_count")
      }
      if (
        "nonnegative_length" %in% pattern_ids &&
          grepl(
            "expect_(true|gte)\\([^#]*length\\([^#]*\\)[^#]*(>=[[:space:]]*0|,[[:space:]]*0[[:space:]]*\\))",
            line
          )
      ) {
        hit <- c(hit, "nonnegative_length")
      }
      if (
        "lower_bound_only_accuracy" %in% pattern_ids &&
          grepl(
            "expect_true\\([^#]*(accuracy|n_correct|correct_count)[^#]*>=[[:space:]]*0",
            line
          ) &&
          !grepl("<=[[:space:]]*1", line)
      ) {
        hit <- c(hit, "lower_bound_only_accuracy")
      }
      if (length(hit)) {
        out[[length(out) + 1L]] <- data.frame(
          pattern = paste(unique(hit), collapse = ","),
          path = sub(paste0("^", root, "/?"), "", path),
          line = i,
          text = trimws(line),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(out)) {
    return(data.frame(
      pattern = character(), path = character(), line = integer(),
      text = character(), stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, out)
}

#' Validate the alignment CI contract against the method registry
#'
#' @param contract Output from [read_alignment_test_contract()].
#' @param registry Output from `alignment_method_registry()`.
#' @param root Package root.
#' @return A data frame of contract problems; zero rows means valid.
#' @keywords internal
validate_alignment_test_contract <- function(contract, registry, root = ".") {
  root <- normalizePath(root, mustWork = TRUE)
  problems <- list()
  add <- function(...) {
    problems[[length(problems) + 1L]] <<- .alignment_ci_problem(...)
  }

  required_fields <- c(
    "families", "pr_files", "nightly_files", "mandatory_tests",
    "release_evidence", "known_gaps"
  )
  allowed <- contract$allowed_families
  stable_required <- contract$stable_required_families
  records <- contract$method_records

  if (!identical(contract$schema_version, 1L)) {
    add("schema_version", "ci_contract schema_version must be 1.")
  }
  if (!length(allowed) || length(setdiff(stable_required, allowed))) {
    add(
      "family_schema",
      "stable_required_families must be a non-empty subset of allowed_families."
    )
  }
  if (anyDuplicated(names(records))) {
    add("duplicate_method", "Method contract IDs must be unique.")
  }
  missing_records <- setdiff(names(registry), names(records))
  extra_records <- setdiff(names(records), names(registry))
  if (length(missing_records)) {
    add(
      "missing_method",
      paste("Registry methods missing from CI contract:", paste(missing_records, collapse = ", "))
    )
  }
  if (length(extra_records)) {
    add(
      "unknown_method",
      paste("CI contract methods absent from registry:", paste(extra_records, collapse = ", "))
    )
  }

  suggests <- .alignment_ci_suggests(root)
  for (method in intersect(names(records), names(registry))) {
    record <- records[[method]]
    entry <- registry[[method]]
    missing_fields <- setdiff(required_fields, names(record))
    if (length(missing_fields)) {
      add(
        "missing_field",
        paste("Missing fields:", paste(missing_fields, collapse = ", ")),
        method = method
      )
      next
    }
    unknown_families <- setdiff(record$families, allowed)
    if (length(unknown_families)) {
      add(
        "unknown_family",
        paste("Unknown test families:", paste(unknown_families, collapse = ", ")),
        method = method
      )
    }
    if (identical(entry$lifecycle$state, "experimental") && !length(record$known_gaps)) {
      add(
        "unacknowledged_gap",
        "Experimental methods must name their remaining evidence gaps.",
        method = method
      )
    }

    paths <- unique(c(record$pr_files, record$nightly_files, record$release_evidence))
    for (path in paths) {
      if (!file.exists(file.path(root, path))) {
        add(
          "missing_artifact",
          "Declared test or evidence artifact does not exist.",
          method = method,
          path = path
        )
      }
    }

    for (mandatory in record$mandatory_tests) {
      parts <- .alignment_ci_mandatory_parts(mandatory)
      if (is.na(parts$path)) {
        add(
          "malformed_mandatory_test",
          "Mandatory tests must use path::exact test title.",
          method = method,
          path = mandatory
        )
        next
      }
      full_path <- file.path(root, parts$path)
      if (!file.exists(full_path)) {
        add(
          "missing_mandatory_file",
          "Mandatory test file does not exist.",
          method = method,
          path = parts$path
        )
        next
      }
      if (!parts$path %in% record$pr_files) {
        add(
          "mandatory_not_in_pr",
          "Every mandatory test must run in the PR tier.",
          method = method,
          path = parts$path
        )
      }
      if (!parts$test %in% .alignment_ci_test_titles(full_path)) {
        add(
          "missing_mandatory_test",
          paste("No test_that block has title:", parts$test),
          method = method,
          path = parts$path
        )
        next
      }

      block <- .alignment_ci_mandatory_block(full_path, parts$test)
      optional_skip <- vapply(suggests, function(pkg) {
        any(grepl(
          paste0("skip_if_not_installed\\([[:space:]]*[\"']", pkg, "[\"']"),
          block
        ))
      }, logical(1))
      if (any(optional_skip)) {
        add(
          "optional_skip_in_mandatory_test",
          paste(
            "Mandatory test directly skips optional dependencies:",
            paste(suggests[optional_skip], collapse = ", ")
          ),
          method = method,
          path = parts$path
        )
      }
    }

    correctness <- entry$lifecycle$evidence$correctness
    if (identical(correctness, "oracle")) {
      if (!"oracle" %in% record$families || !length(record$mandatory_tests)) {
        add(
          "unbacked_oracle_claim",
          "Registry oracle claims require an oracle family and a mandatory PR test.",
          method = method
        )
      }
      mandatory_paths <- unique(vapply(
        record$mandatory_tests,
        function(x) .alignment_ci_mandatory_parts(x)$path,
        character(1)
      ))
      unregistered <- setdiff(
        mandatory_paths,
        entry$lifecycle$evidence$oracle_ids
      )
      if (length(unregistered)) {
        add(
          "oracle_registry_path_missing",
          paste(
            "Mandatory oracle paths absent from registry evidence:",
            paste(unregistered, collapse = ", ")
          ),
          method = method
        )
      }
    }

    if (identical(entry$lifecycle$state, "stable")) {
      absent <- setdiff(stable_required, record$families)
      if (length(absent)) {
        add(
          "stable_family_gap",
          paste("Stable method lacks required families:", paste(absent, collapse = ", ")),
          method = method
        )
      }
      if (!length(record$mandatory_tests)) {
        add(
          "stable_oracle_gap",
          "Stable methods require at least one mandatory PR oracle.",
          method = method
        )
      }
      if (!length(record$release_evidence)) {
        add(
          "stable_release_evidence_gap",
          "Stable methods require a release evidence artifact.",
          method = method
        )
      }
    }
  }

  tautologies <- tryCatch(
    .alignment_ci_tautologies(root, contract$forbidden_tautologies),
    error = function(e) e
  )
  if (inherits(tautologies, "error")) {
    add("tautology_schema", conditionMessage(tautologies))
  } else if (nrow(tautologies)) {
    for (i in seq_len(nrow(tautologies))) {
      add(
        "tautological_expectation",
        paste0(
          "Forbidden ", tautologies$pattern[[i]], " expectation at line ",
          tautologies$line[[i]], ": ", tautologies$text[[i]]
        ),
        path = tautologies$path[[i]]
      )
    }
  }

  if (!length(problems)) {
    return(.alignment_ci_empty_problems())
  }
  out <- do.call(rbind, problems)
  rownames(out) <- NULL
  out
}

alignment_ci_file_methods <- function(contract, tier = "pr") {
  tier <- match.arg(tier, c("pr", "nightly", "release"))
  records <- contract$method_records
  out <- list()
  for (method in names(records)) {
    fields <- switch(
      tier,
      pr = "pr_files",
      nightly = c("pr_files", "nightly_files"),
      release = c("pr_files", "nightly_files")
    )
    paths <- unique(unlist(records[[method]][fields], use.names = FALSE))
    for (path in paths) {
      out[[path]] <- unique(c(out[[path]], method))
    }
  }
  out
}

alignment_ci_mandatory_index <- function(contract) {
  out <- list()
  for (method in names(contract$method_records)) {
    mandatory <- contract$method_records[[method]]$mandatory_tests
    for (item in mandatory) {
      parts <- .alignment_ci_mandatory_parts(item)
      key <- paste(parts$path, parts$test, sep = "::")
      out[[key]] <- unique(c(out[[key]], method))
    }
  }
  out
}

.alignment_ci_count <- function(x) {
  if (is.logical(x)) {
    return(as.integer(isTRUE(x)))
  }
  value <- suppressWarnings(as.integer(x))
  if (length(value) != 1L || is.na(value)) 0L else value
}

.alignment_ci_result_status <- function(row) {
  if (.alignment_ci_count(row$error) > 0L) {
    "error"
  } else if (.alignment_ci_count(row$failed) > 0L) {
    "failed"
  } else if (.alignment_ci_count(row$skipped) > 0L) {
    "skipped"
  } else if (.alignment_ci_count(row$warning) > 0L) {
    "warning"
  } else {
    "passed"
  }
}

.alignment_ci_severity <- function(status, mandatory, tier) {
  if (status %in% c("error", "failed") ||
      (identical(status, "skipped") && isTRUE(mandatory))) {
    return("critical")
  }
  if (identical(status, "warning") ||
      (identical(status, "skipped") && tier != "pr")) {
    return("warning")
  }
  if (identical(status, "skipped")) "info" else "ok"
}

#' Convert testthat results into method-aware CI evidence
#'
#' @param results A `testthat_results` object or equivalent data frame.
#' @param contract Alignment test contract.
#' @param tier One of `pr`, `nightly`, or `release`.
#' @return Detailed result rows with method, status, severity, and mandatory
#'   flags. Shared test files produce one row per declared method.
#' @keywords internal
alignment_ci_results_detail <- function(results, contract, tier = "pr") {
  tier <- match.arg(tier, c("pr", "nightly", "release"))
  raw <- as.data.frame(results)
  file_methods <- alignment_ci_file_methods(contract, tier)
  mandatory_index <- alignment_ci_mandatory_index(contract)
  rows <- list()

  for (i in seq_len(nrow(raw))) {
    file <- as.character(raw$file[[i]])
    path <- if (grepl("tests/testthat/", file, fixed = TRUE)) {
      sub("^.*tests/testthat/", "tests/testthat/", file)
    } else {
      file.path("tests", "testthat", file)
    }
    methods <- file_methods[[path]]
    if (is.null(methods) || !length(methods)) {
      methods <- "<package>"
    }
    title <- as.character(raw$test[[i]])
    mandatory_methods <- mandatory_index[[paste(path, title, sep = "::")]]
    status <- .alignment_ci_result_status(raw[i, , drop = FALSE])

    for (method in methods) {
      mandatory <- method %in% mandatory_methods
      rows[[length(rows) + 1L]] <- data.frame(
        tier = tier,
        method = method,
        file = path,
        test = title,
        status = status,
        severity = .alignment_ci_severity(status, mandatory, tier),
        mandatory = mandatory,
        failed = .alignment_ci_count(raw$failed[[i]]),
        errors = .alignment_ci_count(raw$error[[i]]),
        warnings = .alignment_ci_count(raw$warning[[i]]),
        skipped = .alignment_ci_count(raw$skipped[[i]]),
        passed = .alignment_ci_count(raw$passed[[i]]),
        duration_sec = as.numeric(raw$real[[i]]),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(rows)) {
    return(data.frame(
      tier = character(), method = character(), file = character(),
      test = character(), status = character(), severity = character(),
      mandatory = logical(), failed = integer(), errors = integer(),
      warnings = integer(), skipped = integer(), passed = integer(),
      duration_sec = numeric(), stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

alignment_ci_results_summary <- function(detail) {
  if (!nrow(detail)) {
    return(data.frame(
      tier = character(), method = character(), severity = character(),
      status = character(), count = integer(), stringsAsFactors = FALSE
    ))
  }
  out <- stats::aggregate(
    rep.int(1L, nrow(detail)),
    detail[c("tier", "method", "severity", "status")],
    sum
  )
  names(out)[[ncol(out)]] <- "count"
  out[order(out$method, out$severity, out$status), , drop = FALSE]
}

alignment_ci_contract_status <- function(contract, tier = "pr") {
  tier <- match.arg(tier, c("pr", "nightly", "release"))
  rows <- list()
  for (method in names(contract$method_records)) {
    record <- contract$method_records[[method]]
    for (gap in record$known_gaps) {
      rows[[length(rows) + 1L]] <- data.frame(
        tier = tier,
        method = method,
        kind = "known_gap",
        item = gap,
        status = "open",
        severity = if (tier == "pr") "info" else "warning",
        stringsAsFactors = FALSE
      )
    }
    if (tier == "release") {
      if (length(record$release_evidence)) {
        for (path in record$release_evidence) {
          rows[[length(rows) + 1L]] <- data.frame(
            tier = tier,
            method = method,
            kind = "release_evidence",
            item = path,
            status = "available",
            severity = "ok",
            stringsAsFactors = FALSE
          )
        }
      } else {
        rows[[length(rows) + 1L]] <- data.frame(
          tier = tier,
          method = method,
          kind = "release_evidence",
          item = "none_declared",
          status = "not_available",
          severity = "warning",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(data.frame(
      tier = character(), method = character(), kind = character(),
      item = character(), status = character(), severity = character(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
