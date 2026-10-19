library(testthat)
library(manifoldalign)

test_that("generated comparison table cannot drift from registry", {
  root <- test_path("..", "..")
  # The tracked docs source is authoritative; no untracked root copy.
  paths <- file.path(root, "docs", "METHOD_COMPARISON_TABLE.md")
  copies <- lapply(paths, readLines, warn = FALSE)
  for (lines in copies) {
    begin <- match("<!-- BEGIN GENERATED METHOD TABLE -->", lines)
    end <- match("<!-- END GENERATED METHOD TABLE -->", lines)
    expect_false(is.na(begin))
    expect_false(is.na(end))
    expect_gt(end, begin)
    expect_identical(
      lines[seq.int(begin + 1L, end - 1L)],
      manifoldalign:::.alignment_registry_markdown()
    )
  }
})


test_that("real-data scripts pin resume, license, and evidence-scope contracts", {
  root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
  prepare <- readLines(file.path(
    root, "inst", "benchmarks", "prepare_real_alignment_data.R"
  ), warn = FALSE)
  extract <- readLines(file.path(
    root, "inst", "benchmarks", "extract_so2sat_ci_fixture.R"
  ), warn = FALSE)
  runner <- readLines(file.path(
    root, "inst", "benchmarks", "run_real_alignment_benchmarks.R"
  ), warn = FALSE)

  expect_true(any(grepl("--continue-at", prepare, fixed = TRUE)))
  expect_true(any(grepl("MANIFOLDALIGN_ACCEPT_SO2SAT_LICENSE", prepare,
    fixed = TRUE
  )))
  expect_true(any(grepl("observed_sha", prepare, fixed = TRUE)))
  expect_true(any(grepl("per-class=all", runner, fixed = TRUE)))
  expect_true(any(grepl("full_frozen_split", runner, fixed = TRUE)))
  expect_true(any(grepl("select_balanced", extract, fixed = TRUE)))
})


test_that("vignette metadata is parsed with a clean UTF-8 encoding", {
  root <- test_path("..", "..")
  meta <- tools::pkgVignettes(dir = root)
  files <- list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$")
  expect_length(meta$docs, length(files))
  expect_true(all(meta$encodings == "UTF-8"))
  expect_true(all(meta$engines == "knitr::rmarkdown"))
})


test_that("bundled vignette benchmarks retain the tracked evidence bytes", {
  root <- test_path("..", "..")
  files <- c("toy_feature_aligners_regime_ranks.csv",
             "toy_feature_seed_audit_summary.csv")
  for (file in files) {
    original <- file.path(root, "dev-scripts", "benchmarks", file)
    bundled <- file.path(root, "inst", "extdata", "benchmark-guidance", file)
    expect_identical(unname(tools::md5sum(original)),
                     unname(tools::md5sum(bundled)))
  }
})
