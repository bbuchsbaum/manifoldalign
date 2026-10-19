#!/usr/bin/env Rscript

find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (
      file.exists(file.path(path, "DESCRIPTION")) &&
        file.exists(file.path(
          path, "tests", "testthat", "ALIGNMENT_TEST_MATRIX.yml"
        ))
    ) {
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
  tier <- "pr"
  output_dir <- Sys.getenv("ALIGNMENT_EVIDENCE_DIR", unset = "")
  for (arg in args) {
    if (arg %in% c("pr", "nightly", "release")) {
      tier <- arg
    } else if (startsWith(arg, "--tier=")) {
      tier <- sub("^--tier=", "", arg)
    } else if (startsWith(arg, "--output-dir=")) {
      output_dir <- sub("^--output-dir=", "", arg)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  if (!tier %in% c("pr", "nightly", "release")) {
    stop("Tier must be one of pr, nightly, or release.", call. = FALSE)
  }
  if (!nzchar(output_dir)) {
    output_dir <- file.path(root, "artifacts", "alignment-ci", tier)
  } else if (!grepl("^(/|[A-Za-z]:[/\\\\])", output_dir)) {
    output_dir <- file.path(root, output_dir)
  }
  list(tier = tier, output_dir = output_dir)
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "")
}

missing_mandatory_rows <- function(detail, contract, tier) {
  rows <- list()
  for (method in names(contract$method_records)) {
    for (item in contract$method_records[[method]]$mandatory_tests) {
      parts <- .alignment_ci_mandatory_parts(item)
      observed <- detail$method == method &
        detail$file == parts$path &
        detail$test == parts$test
      if (!any(observed)) {
        rows[[length(rows) + 1L]] <- data.frame(
          tier = tier,
          method = method,
          file = parts$path,
          test = parts$test,
          status = "missing",
          severity = "critical",
          mandatory = TRUE,
          failed = 0L,
          errors = 1L,
          warnings = 0L,
          skipped = 0L,
          passed = 0L,
          duration_sec = 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(rows)) {
    return(detail[0, , drop = FALSE])
  }
  do.call(rbind, rows)
}

write_markdown_summary <- function(summary, contract_status, path, tier) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    paste0("# Alignment evidence: ", tier),
    "",
    "## Test outcomes",
    "",
    "| Method | Severity | Status | Count |",
    "|---|---|---|---:|"
  ), con)
  if (nrow(summary)) {
    apply(summary, 1L, function(row) {
      writeLines(paste0(
        "| ", row[["method"]], " | ", row[["severity"]], " | ",
        " ", row[["status"]], " | ", row[["count"]], " |"
      ), con)
    })
  } else {
    writeLines("| <package> | critical | no_results | 1 |", con)
  }
  writeLines(c("", "## Open evidence gaps", ""), con)
  gaps <- contract_status[contract_status$kind == "known_gap", , drop = FALSE]
  if (nrow(gaps)) {
    writeLines(paste0("- `", gaps$method, "`: `", gaps$item, "`"), con)
  } else {
    writeLines("- None declared.", con)
  }
}

main <- function() {
  root <- find_package_root()
  args <- parse_args(commandArgs(trailingOnly = TRUE), root)
  tier <- args$tier
  output_dir <- args$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  source(file.path(root, "tools", "alignment-test-contract.R"))
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The alignment CI runner requires the testthat dependency pkgload.", call. = FALSE)
  }
  pkgload::load_all(root, quiet = TRUE)

  contract <- read_alignment_test_contract(file.path(
    root, "tests", "testthat", "ALIGNMENT_TEST_MATRIX.yml"
  ))
  registry <- getExportedValue("manifoldalign", "alignment_method_registry")()
  contract_problems <- validate_alignment_test_contract(
    contract, registry, root = root
  )
  write_csv(
    contract_problems,
    file.path(output_dir, "alignment-contract-problems.csv")
  )
  if (nrow(contract_problems)) {
    print(contract_problems, row.names = FALSE)
    message("Alignment CI contract is invalid; tests were not run.")
    return(2L)
  }

  options(
    manifoldalign.run_benchmarks = tier != "pr",
    manifoldalign.ci_tier = tier
  )
  Sys.setenv(MANIFOLDALIGN_CI_TIER = tier)
  if (tier == "release") {
    Sys.setenv(
      MANIFOLDALIGN_RUN_STRESS = "1",
      WITH_PARROT_PERF_TESTS = "1"
    )
  }

  message("Running alignment evidence tier: ", tier)
  results <- tryCatch(
    testthat::test_local(
      path = root,
      reporter = "silent",
      stop_on_failure = FALSE,
      stop_on_warning = FALSE,
      load_package = "source",
      shuffle = FALSE
    ),
    error = function(e) e
  )

  if (inherits(results, "error")) {
    detail <- data.frame(
      tier = tier,
      method = "<package>",
      file = NA_character_,
      test = conditionMessage(results),
      status = "runner_error",
      severity = "critical",
      mandatory = FALSE,
      failed = 0L,
      errors = 1L,
      warnings = 0L,
      skipped = 0L,
      passed = 0L,
      duration_sec = 0,
      stringsAsFactors = FALSE
    )
  } else {
    detail <- alignment_ci_results_detail(results, contract, tier = tier)
    absent <- missing_mandatory_rows(detail, contract, tier)
    if (nrow(absent)) {
      detail <- rbind(detail, absent)
    }
  }

  # Source-layout and evidence audits cannot resolve checkout paths from an
  # installed namespace. Run them explicitly in every tier; their failures
  # are critical and remain in the same evidence tables as runtime tests.
  source_results <- testthat::test_dir(
    file.path(root, "tools", "tests"), reporter = "silent",
    stop_on_failure = FALSE, stop_on_warning = FALSE, shuffle = FALSE
  )
  source_detail <- alignment_ci_results_detail(source_results, contract, tier)
  source_detail$file <- sub("^tests/testthat/", "tools/tests/", source_detail$file)
  detail <- rbind(detail, source_detail)

  summary <- alignment_ci_results_summary(detail)
  contract_status <- alignment_ci_contract_status(contract, tier = tier)
  write_csv(detail, file.path(output_dir, "alignment-test-results.csv"))
  write_csv(summary, file.path(output_dir, "alignment-test-summary.csv"))
  write_csv(
    contract_status,
    file.path(output_dir, "alignment-contract-status.csv")
  )
  markdown <- file.path(output_dir, "alignment-test-summary.md")
  write_markdown_summary(summary, contract_status, markdown, tier)

  print(summary, row.names = FALSE)
  message("Evidence written to ", output_dir)
  github_summary <- Sys.getenv("GITHUB_STEP_SUMMARY", unset = "")
  if (nzchar(github_summary)) {
    cat(readLines(markdown, warn = FALSE), sep = "\n", file = github_summary, append = TRUE)
    cat("\n", file = github_summary, append = TRUE)
  }

  if (any(detail$severity == "critical")) 1L else 0L
}

status <- main()
quit(save = "no", status = status, runLast = FALSE)
