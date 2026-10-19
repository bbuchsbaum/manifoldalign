library(testthat)
library(manifoldalign)

alignment_ci_root <- normalizePath(test_path("..", ".."), mustWork = TRUE)
# The contract tooling lives in the source tree and is not part of the built
# package; the PR evidence tier runs these tests from the source checkout.
skip_if_not(
  file.exists(file.path(alignment_ci_root, "tools", "alignment-test-contract.R")),
  "alignment CI contract tooling is only available in the source checkout"
)
sys.source(
  file.path(alignment_ci_root, "tools", "alignment-test-contract.R"),
  envir = environment()
)

read_ci_contract <- function() {
  read_alignment_test_contract(
    file.path(alignment_ci_root, "tests", "testthat", "ALIGNMENT_TEST_MATRIX.yml")
  )
}

test_that("alignment test matrix is an enforced registry-complete contract", {
  contract <- read_ci_contract()
  registry <- alignment_method_registry()
  problems <- validate_alignment_test_contract(
    contract, registry, root = alignment_ci_root
  )

  expect_s3_class(contract, "alignment_test_contract")
  expect_setequal(names(contract$method_records), names(registry))
  expect_equal(
    nrow(problems),
    0L,
    info = paste(problems$message, collapse = "\n")
  )
})

test_that("stable promotion fails when any computational family disappears", {
  contract <- read_ci_contract()
  registry <- alignment_method_registry()
  registry$rapid_ma$lifecycle$state <- "stable"
  contract$method_records$rapid_ma$families <- setdiff(
    contract$method_records$rapid_ma$families,
    "metamorphic"
  )

  problems <- validate_alignment_test_contract(
    contract, registry, root = alignment_ci_root
  )
  rapid <- problems[problems$method == "rapid_ma", , drop = FALSE]

  expect_true("stable_family_gap" %in% rapid$code)
  expect_match(
    rapid$message[rapid$code == "stable_family_gap"],
    "metamorphic"
  )
})

test_that("oracle evidence claims require a mandatory PR path", {
  contract <- read_ci_contract()
  registry <- alignment_method_registry()
  contract$method_records$global_geo_align$mandatory_tests <- character()

  problems <- validate_alignment_test_contract(
    contract, registry, root = alignment_ci_root
  )

  expect_true(any(
    problems$method == "global_geo_align" &
      problems$code == "unbacked_oracle_claim"
  ))
})

test_that("mandatory tests cannot directly skip a suggested dependency", {
  contract <- read_ci_contract()
  registry <- alignment_method_registry()
  contract$method_records$lowrank_align$mandatory_tests <- paste0(
    "tests/testthat/test-lowrank_align.R::",
    "lowrank_align agrees across solver backends"
  )

  problems <- validate_alignment_test_contract(
    contract, registry, root = alignment_ci_root
  )

  expect_true(any(
    problems$method == "lowrank_align" &
      problems$code == "optional_skip_in_mandatory_test"
  ))
})

test_that("tautological correctness counts are detected", {
  root <- tempfile("alignment-ci-tautology-")
  dir.create(file.path(root, "tests", "testthat"), recursive = TRUE)
  bad_line <- paste(
    c("expect_gte(", "n_correct", ", 0", ")"),
    collapse = ""
  )
  writeLines(
    c('test_that("bad evidence", {', bad_line, "})"),
    file.path(root, "tests", "testthat", "test-bad.R")
  )

  hits <- .alignment_ci_tautologies(
    root,
    c("nonnegative_correct_count", "nonnegative_length")
  )

  expect_equal(nrow(hits), 1L)
  expect_identical(hits$pattern, "nonnegative_correct_count")
})

test_that("mandatory skips become method-level critical evidence", {
  contract <- read_ci_contract()
  title <- paste0(
    "full-support sparse and dense UOT solutions agree"
  )
  fake <- data.frame(
    file = "test-rapid-ma-uot.R",
    context = "",
    test = title,
    nb = 1L,
    failed = 0L,
    skipped = TRUE,
    error = FALSE,
    warning = 0L,
    user = 0,
    system = 0,
    real = 0,
    passed = 0L,
    stringsAsFactors = FALSE
  )

  detail <- alignment_ci_results_detail(fake, contract, tier = "pr")
  rapid <- detail[detail$method == "rapid_ma", , drop = FALSE]

  expect_equal(nrow(rapid), 1L)
  expect_true(rapid$mandatory)
  expect_identical(rapid$status, "skipped")
  expect_identical(rapid$severity, "critical")
})

test_that("release status exposes checked evidence and remaining named gaps", {
  status <- alignment_ci_contract_status(read_ci_contract(), tier = "release")

  expect_false(any(
    status$method == "fpgw" &
      status$kind == "known_gap" &
      status$item == "frozen_v2_tournament_pending"
  ))
  expect_true(any(
    status$method == "fpgw" &
      status$kind == "release_evidence" &
      status$item == "docs/evidence/alignment-tournament-v2/raw_results.csv" &
      status$status == "available"
  ))
  expect_true(any(
    status$method == "rapid_ma" &
      status$kind == "release_evidence" &
      status$status == "available"
  ))
  expect_true(any(
    status$method == "rapid_ma" &
      status$kind == "known_gap" &
      status$item == "position_gated_fusion_missing"
  ))
})
