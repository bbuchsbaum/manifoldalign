library(testthat)
library(manifoldalign)

stable_lifecycle_record <- function(method = "stable_method") {
  list(
    method = method,
    state = "stable",
    recommended = TRUE,
    quarantined = FALSE,
    api_enabled = TRUE,
    findings = character(),
    evidence = list(
      correctness = "oracle",
      accuracy = "comparative",
      scalability = "comparative",
      oracle_ids = "tests/oracle.R",
      benchmark_ids = "benchmarks/comparative-v1.csv",
      reviewed_by = "independent-reviewer",
      reviewed_on = "2026-07-22"
    )
  )
}

experimental_lifecycle_record <- function(method = "experimental_method") {
  record <- stable_lifecycle_record(method)
  record$state <- "experimental"
  record$recommended <- FALSE
  record$evidence$correctness <- "smoke"
  record$evidence$accuracy <- "none"
  record$evidence$scalability <- "none"
  record$evidence$oracle_ids <- character()
  record$evidence$benchmark_ids <- character()
  record$evidence$reviewed_by <- NULL
  record$evidence$reviewed_on <- NULL
  record
}

deprecated_lifecycle_record <- function(method = "old_method") {
  record <- experimental_lifecycle_record(method)
  record$state <- "deprecated"
  record$api_enabled <- TRUE
  record$findings <- "low_accuracy"
  record$deprecation <- list(
    reason = "dominance",
    replacement = "better_method",
    warning_since = "0.2.0",
    warning_through = "0.3.0",
    removal_at = "1.0.0",
    migration = "Replace old_method() with better_method().",
    archive = "inst/legacy/old_method",
    dominance = list(
      applicable_regimes = c("small_dense", "large_sparse", "partial_overlap"),
      unique_supported_capability = FALSE,
      evidence_ids = "benchmarks/dominance-v1.csv"
    )
  )
  record
}

problem_codes <- function(record, ...) {
  validate_alignment_lifecycle(record, ...)$code
}

test_that("policy defines states, evidence dimensions, and trigger semantics", {
  policy <- alignment_lifecycle_policy()

  expect_identical(
    policy$states$state,
    c("stable", "experimental", "legacy_reference", "deprecated", "removed")
  )
  expect_identical(
    names(policy$evidence_levels),
    c("correctness", "accuracy", "scalability")
  )
  expect_true("incoherent_implementation" %in% policy$deprecation_reasons)
  expect_true(all(c(
    "fabricated_validation", "failed_numerical_invariant",
    "silent_unsupported_input", "wrong_result_risk"
  ) %in% policy$quarantine_triggers$finding[
    policy$quarantine_triggers$quarantine_required
  ]))
  expect_false(policy$quarantine_triggers$quarantine_required[
    policy$quarantine_triggers$finding == "low_accuracy"
  ])
})

test_that("a stable record needs oracle and comparative evidence", {
  record <- stable_lifecycle_record()

  expect_s3_class(validate_alignment_lifecycle(record),
                  "alignment_lifecycle_problems")
  expect_equal(nrow(validate_alignment_lifecycle(record)), 0L)
  expect_true(assert_alignment_release_ready(record))

  record$evidence$correctness <- "invariant"
  record$evidence$benchmark_ids <- character()
  expect_true(all(c("insufficient_evidence", "stable_benchmark_missing") %in%
                    problem_codes(record)))
})

test_that("low accuracy is not treated as a correctness failure", {
  record <- experimental_lifecycle_record()
  record$findings <- "low_accuracy"

  problems <- validate_alignment_lifecycle(record)
  expect_equal(nrow(problems), 0L)
  expect_false("quarantine_required" %in% problems$code)
})

test_that("each correctness trigger forces immediate quarantine", {
  policy <- alignment_lifecycle_policy()
  triggers <- policy$quarantine_triggers$finding[
    policy$quarantine_triggers$quarantine_required
  ]

  for (trigger in triggers) {
    record <- experimental_lifecycle_record(trigger)
    record$findings <- trigger
    expect_true(
      "quarantine_required" %in% problem_codes(record),
      info = trigger
    )
  }
})

test_that("quarantine requires a deadline, plan, and safe API state", {
  record <- experimental_lifecycle_record("fabricated_validator")
  record$findings <- "fabricated_validation"
  record$quarantined <- TRUE
  record$api_enabled <- TRUE
  record$quarantine <- list(
    triggers = "fabricated_validation",
    reason = "The validator returned fixture constants.",
    repair_by = "2026-08-22",
    repair_plan = "Replace it with a dependency-independent exact oracle."
  )

  expect_true("unsafe_api_enabled" %in% problem_codes(record))
  record$api_enabled <- FALSE
  expect_equal(nrow(validate_alignment_lifecycle(record)), 0L)

  record$quarantine$repair_by <- "22-08-2026"
  expect_true("invalid_repair_deadline" %in% problem_codes(record))
})

test_that("release gate rejects recommended deprecated or quarantined methods", {
  deprecated <- deprecated_lifecycle_record()
  deprecated$recommended <- TRUE
  expect_true("state_not_recommendable" %in% problem_codes(deprecated))
  expect_error(
    assert_alignment_release_ready(deprecated),
    class = "alignment_release_error"
  )

  quarantined <- experimental_lifecycle_record("risky_method")
  quarantined$findings <- "silent_unsupported_input"
  quarantined$quarantined <- TRUE
  quarantined$recommended <- TRUE
  quarantined$quarantine <- list(
    triggers = "silent_unsupported_input",
    reason = "An advertised edge-weight input was ignored.",
    repair_by = "2026-08-22",
    repair_plan = "Reject the input until its semantics are implemented."
  )
  expect_true("quarantined_recommended" %in% problem_codes(quarantined))
  expect_error(
    assert_alignment_release_ready(quarantined),
    class = "alignment_release_error"
  )
})

test_that("deprecation cadence and warning text identify replacement and removal", {
  record <- deprecated_lifecycle_record()
  expect_equal(nrow(validate_alignment_lifecycle(record)), 0L)

  message <- manifoldalign:::.alignment_deprecation_message(record)
  expect_match(message, "better_method", fixed = TRUE)
  expect_match(message, "1.0.0", fixed = TRUE)

  record$deprecation$warning_through <- "0.2.1"
  expect_true("warning_window_too_short" %in% problem_codes(record))
  record$deprecation$warning_through <- "0.3.0"
  record$deprecation$removal_at <- "2.0.0"
  expect_true("invalid_removal_version" %in% problem_codes(record))
})

test_that("performance deprecation requires three regimes and no unique niche", {
  record <- deprecated_lifecycle_record()
  record$deprecation$dominance$applicable_regimes <- c("small", "large")
  record$deprecation$dominance$unique_supported_capability <- TRUE

  codes <- problem_codes(record)
  expect_true("dominance_regimes_insufficient" %in% codes)
  expect_true("unique_capability_blocks_deprecation" %in% codes)

  record$deprecation$dominance$applicable_regimes <- c(
    "small", "large", "partial"
  )
  record$deprecation$dominance$unique_supported_capability <- FALSE
  expect_equal(nrow(validate_alignment_lifecycle(record)), 0L)
})

test_that("state transitions prevent bypassing deprecation cadence", {
  experimental <- experimental_lifecycle_record()
  expect_equal(
    nrow(validate_alignment_lifecycle(experimental, previous_state = "stable")),
    0L
  )

  removed <- deprecated_lifecycle_record()
  removed$state <- "removed"
  removed$api_enabled <- FALSE
  expect_true(
    "invalid_transition" %in%
      problem_codes(removed, previous_state = "stable")
  )
  expect_equal(
    nrow(validate_alignment_lifecycle(removed, previous_state = "deprecated")),
    0L
  )
})

test_that("stable restoration from quarantine requires reviewed fresh evidence", {
  record <- stable_lifecycle_record("repaired_method")

  expect_true(
    "restoration_evidence_missing" %in%
      problem_codes(record, previous_quarantined = TRUE)
  )
  record$restoration <- list(
    resolved_triggers = "failed_numerical_invariant",
    oracle_ids = "tests/repaired-oracle.R",
    benchmark_ids = "benchmarks/repaired-v2.csv",
    reviewed_by = "independent-reviewer",
    reviewed_on = "2026-07-22"
  )
  expect_equal(
    nrow(validate_alignment_lifecycle(record, previous_quarantined = TRUE)),
    0L
  )
})

test_that("release gate rejects duplicate lifecycle records", {
  record <- stable_lifecycle_record("duplicate")
  condition <- expect_error(
    assert_alignment_release_ready(list(record, record)),
    class = "alignment_release_error"
  )
  expect_true("duplicate_method" %in% condition$problems$code)
})
