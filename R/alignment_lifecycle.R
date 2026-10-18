# Alignment-method lifecycle policy ---------------------------------------

.alignment_lifecycle_problem <- function(method, code, field, message) {
  method_label <- if (is.character(method) && length(method) == 1L &&
                      !is.na(method) && nzchar(method)) {
    method
  } else {
    "<unknown>"
  }
  data.frame(
    method = method_label,
    code = as.character(code),
    severity = "error",
    field = as.character(field),
    message = as.character(message),
    stringsAsFactors = FALSE
  )
}

.alignment_empty_problems <- function() {
  out <- data.frame(
    method = character(),
    code = character(),
    severity = character(),
    field = character(),
    message = character(),
    stringsAsFactors = FALSE
  )
  class(out) <- c("alignment_lifecycle_problems", "data.frame")
  out
}

.alignment_scalar_character <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
}

.alignment_scalar_logical <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

.alignment_nonempty_character <- function(x) {
  is.character(x) && length(x) > 0L && all(!is.na(x) & nzchar(x))
}

.alignment_valid_date <- function(x) {
  if (!.alignment_scalar_character(x)) return(FALSE)
  parsed <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
  !is.na(parsed) && identical(format(parsed, "%Y-%m-%d"), x)
}

.alignment_version_parts <- function(x) {
  if (!.alignment_scalar_character(x) ||
      !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", x)) {
    return(NULL)
  }
  as.integer(strsplit(x, ".", fixed = TRUE)[[1L]])
}

.alignment_level_rank <- function(policy, dimension, level) {
  levels <- policy$evidence_levels[[dimension]]
  if (is.null(levels) || !.alignment_scalar_character(level)) return(NA_integer_)
  index <- match(level, levels$level)
  if (is.na(index)) NA_integer_ else levels$rank[[index]]
}

.alignment_record_field <- function(record, name, default = NULL) {
  if (is.list(record) && name %in% names(record)) record[[name]] else default
}

.alignment_or <- function(x, default) {
  if (is.null(x)) default else x
}

#' Alignment-method lifecycle policy
#'
#' Return and enforce the package-wide policy used to distinguish numerical
#' correctness from empirical performance. The policy defines five lifecycle
#' states, independent evidence levels for correctness, accuracy, and
#' scalability, immediate quarantine triggers, permitted state transitions,
#' and deprecation cadence.
#'
#' `validate_alignment_lifecycle()` validates one method record. A record is a
#' named list with `method`, `state`, `recommended`, `quarantined`,
#' `api_enabled`, `findings`, and `evidence` fields. Stable records require
#' oracle-level correctness plus comparative accuracy and scalability evidence.
#' Deprecated and removed records additionally require a `deprecation` list;
#' see `docs/ALIGNMENT_LIFECYCLE.md` for the complete schema.
#'
#' `assert_alignment_release_ready()` validates a list of records and raises an
#' `alignment_release_error` when any lifecycle invariant is violated. In
#' particular, deprecated, removed, experimental, legacy/reference, or
#' quarantined methods cannot be presented as recommended.
#'
#' @param record A named list containing lifecycle metadata for one method.
#' @param previous_state Optional prior lifecycle state. When supplied, the
#'   transition to `record$state` is validated.
#' @param previous_quarantined Whether the method is being restored from an
#'   active quarantine. Restoration to stable requires independent review and
#'   explicit oracle and benchmark evidence.
#' @param policy A policy object returned by `alignment_lifecycle_policy()`.
#' @param records One lifecycle record or a list of lifecycle records.
#'
#' @return `alignment_lifecycle_policy()` returns a named policy list.
#'   `validate_alignment_lifecycle()` returns an
#'   `alignment_lifecycle_problems` data frame (zero rows means valid).
#'   `assert_alignment_release_ready()` invisibly returns `TRUE`, or raises an
#'   error containing the combined problems in its `problems` field.
#'
#' @examples
#' record <- list(
#'   method = "example_aligner",
#'   state = "stable",
#'   recommended = TRUE,
#'   quarantined = FALSE,
#'   api_enabled = TRUE,
#'   findings = character(),
#'   evidence = list(
#'     correctness = "oracle",
#'     accuracy = "comparative",
#'     scalability = "comparative",
#'     oracle_ids = "tests/testthat/test-example-oracle.R",
#'     benchmark_ids = "benchmarks/example-v1.csv",
#'     reviewed_by = "maintainer",
#'     reviewed_on = "2026-07-22"
#'   )
#' )
#' validate_alignment_lifecycle(record)
#' assert_alignment_release_ready(record)
#'
#' @name alignment_lifecycle_policy
NULL

#' @rdname alignment_lifecycle_policy
#' @export
alignment_lifecycle_policy <- function() {
  states <- data.frame(
    state = c(
      "stable", "experimental", "legacy_reference", "deprecated", "removed"
    ),
    label = c(
      "stable", "experimental", "legacy/reference", "deprecated", "removed"
    ),
    recommendation_allowed = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    public_api_allowed = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  evidence_levels <- list(
    correctness = data.frame(
      level = c("none", "smoke", "invariant", "oracle"),
      rank = 0:3,
      stringsAsFactors = FALSE
    ),
    accuracy = data.frame(
      level = c("none", "single_fixture", "comparative", "release"),
      rank = 0:3,
      stringsAsFactors = FALSE
    ),
    scalability = data.frame(
      level = c("none", "smoke", "comparative", "release"),
      rank = 0:3,
      stringsAsFactors = FALSE
    )
  )

  state_requirements <- data.frame(
    state = states$state,
    correctness = c("oracle", "smoke", "none", "none", "none"),
    accuracy = c("comparative", "none", "none", "none", "none"),
    scalability = c("comparative", "none", "none", "none", "none"),
    stringsAsFactors = FALSE
  )

  transitions <- data.frame(
    from = c(
      "stable", "stable", "stable",
      "experimental", "experimental", "experimental",
      "legacy_reference", "legacy_reference",
      "deprecated", "deprecated"
    ),
    to = c(
      "experimental", "legacy_reference", "deprecated",
      "stable", "legacy_reference", "deprecated",
      "experimental", "deprecated",
      "experimental", "removed"
    ),
    stringsAsFactors = FALSE
  )

  quarantine_triggers <- data.frame(
    finding = c(
      "fabricated_validation",
      "failed_numerical_invariant",
      "silent_unsupported_input",
      "wrong_result_risk",
      "low_accuracy"
    ),
    quarantine_required = c(TRUE, TRUE, TRUE, TRUE, FALSE),
    disable_api = c(TRUE, TRUE, FALSE, TRUE, FALSE),
    rationale = c(
      "Validation that does not compute the claimed quantity is not evidence.",
      "A materially failed invariant makes numerical results untrustworthy.",
      "Unsupported inputs must be rejected rather than silently ignored.",
      "A credible wrong-result risk requires containment before release.",
      "Accuracy is comparative evidence, not by itself a correctness verdict."
    ),
    stringsAsFactors = FALSE
  )

  list(
    schema_version = "1.0.0",
    states = states,
    evidence_levels = evidence_levels,
    state_requirements = state_requirements,
    transitions = transitions,
    quarantine_triggers = quarantine_triggers,
    deprecation_reasons = c(
      "correctness", "fabricated_validation", "unsupported_contract",
      "wrong_result_risk", "incoherent_implementation", "dominance"
    ),
    deprecation_required_fields = c(
      "reason", "replacement", "warning_since", "warning_through",
      "removal_at", "migration", "archive"
    ),
    quarantine_required_fields = c(
      "triggers", "reason", "repair_by", "repair_plan"
    ),
    restoration_required_fields = c(
      "resolved_triggers", "oracle_ids", "benchmark_ids", "reviewed_by",
      "reviewed_on"
    ),
    dominance_minimum_regimes = 3L
  )
}

.alignment_check_deprecation <- function(record, method, policy, add_problem) {
  deprecation <- .alignment_record_field(record, "deprecation")
  if (!is.list(deprecation)) {
    add_problem(
      "deprecation_metadata_missing", "deprecation",
      "Deprecated and removed methods require a `deprecation` metadata list."
    )
    return(invisible(NULL))
  }

  for (field in policy$deprecation_required_fields) {
    if (!(field %in% names(deprecation)) ||
        !.alignment_scalar_character(deprecation[[field]])) {
      add_problem(
        "deprecation_metadata_missing", paste0("deprecation$", field),
        paste0("Deprecation metadata requires a non-empty `", field, "` field.")
      )
    }
  }

  reason <- deprecation$reason
  if (.alignment_scalar_character(reason) &&
      !(reason %in% policy$deprecation_reasons)) {
    add_problem(
      "invalid_deprecation_reason", "deprecation$reason",
      paste0(
        "`deprecation$reason` must be one of: ",
        paste(policy$deprecation_reasons, collapse = ", "), "."
      )
    )
  }

  if (.alignment_scalar_character(deprecation$replacement) &&
      identical(deprecation$replacement, method)) {
    add_problem(
      "invalid_replacement", "deprecation$replacement",
      "A deprecated method must name a distinct replacement."
    )
  }

  warning_since <- .alignment_version_parts(deprecation$warning_since)
  warning_through <- .alignment_version_parts(deprecation$warning_through)
  removal_at <- .alignment_version_parts(deprecation$removal_at)
  if (is.null(warning_since)) {
    add_problem(
      "invalid_warning_version", "deprecation$warning_since",
      "`warning_since` must be a three-part package version such as 0.2.0."
    )
  }
  if (is.null(warning_through)) {
    add_problem(
      "invalid_warning_version", "deprecation$warning_through",
      "`warning_through` must be a three-part package version such as 0.3.0."
    )
  }
  if (is.null(removal_at)) {
    add_problem(
      "invalid_removal_version", "deprecation$removal_at",
      "`removal_at` must be a three-part package version such as 1.0.0."
    )
  }

  if (!is.null(warning_since) && !is.null(warning_through)) {
    one_minor <- warning_since
    one_minor[[2L]] <- one_minor[[2L]] + 1L
    one_minor[[3L]] <- 0L
    if (warning_through[[1L]] != warning_since[[1L]] ||
        utils::compareVersion(
          paste(warning_through, collapse = "."),
          paste(one_minor, collapse = ".")
        ) < 0L) {
      add_problem(
        "warning_window_too_short", "deprecation$warning_through",
        "The warning and migration path must span at least one minor release."
      )
    }
  }

  if (!is.null(warning_since) && !is.null(removal_at) &&
      removal_at[[1L]] != warning_since[[1L]] + 1L) {
    add_problem(
      "invalid_removal_version", "deprecation$removal_at",
      "Removal must target the next major package version."
    )
  }

  if (.alignment_scalar_character(reason)) {
    findings <- .alignment_record_field(record, "findings", character())
    finding_for_reason <- c(
      correctness = "failed_numerical_invariant",
      fabricated_validation = "fabricated_validation",
      unsupported_contract = "silent_unsupported_input",
      wrong_result_risk = "wrong_result_risk"
    )
    if (reason %in% names(finding_for_reason) &&
        !(unname(finding_for_reason[[reason]]) %in% findings)) {
      add_problem(
        "deprecation_reason_unsubstantiated", "deprecation$reason",
        paste0(
          "Deprecation reason `", reason, "` requires finding `",
          unname(finding_for_reason[[reason]]), "`."
        )
      )
    }
  }

  if (identical(reason, "dominance")) {
    dominance <- deprecation$dominance
    if (!is.list(dominance)) {
      add_problem(
        "dominance_evidence_missing", "deprecation$dominance",
        "Performance-based deprecation requires a `dominance` evidence list."
      )
      return(invisible(NULL))
    }
    regimes <- dominance$applicable_regimes
    if (!.alignment_nonempty_character(regimes) ||
        length(unique(regimes)) < policy$dominance_minimum_regimes) {
      add_problem(
        "dominance_regimes_insufficient", "deprecation$dominance$applicable_regimes",
        paste0(
          "Dominance must be demonstrated in at least ",
          policy$dominance_minimum_regimes, " applicable regimes."
        )
      )
    }
    unique_capability <- dominance$unique_supported_capability
    if (!identical(unique_capability, FALSE)) {
      add_problem(
        "unique_capability_blocks_deprecation",
        "deprecation$dominance$unique_supported_capability",
        "A method with a unique supported capability cannot be deprecated solely for dominance."
      )
    }
    if (!.alignment_nonempty_character(dominance$evidence_ids)) {
      add_problem(
        "dominance_evidence_missing", "deprecation$dominance$evidence_ids",
        "Dominance claims require one or more reproducible benchmark artifacts."
      )
    }
  }
  invisible(NULL)
}

#' @rdname alignment_lifecycle_policy
#' @export
validate_alignment_lifecycle <- function(
  record,
  previous_state = NULL,
  previous_quarantined = FALSE,
  policy = alignment_lifecycle_policy()
) {
  problems <- .alignment_empty_problems()
  method <- .alignment_record_field(record, "method", "<unknown>")
  add_problem <- function(code, field, message) {
    problems <<- rbind(
      problems,
      .alignment_lifecycle_problem(method, code, field, message)
    )
  }

  if (!is.list(record) || is.data.frame(record)) {
    add_problem(
      "invalid_record", "record",
      "A lifecycle record must be a named list for exactly one method."
    )
    class(problems) <- c("alignment_lifecycle_problems", "data.frame")
    return(problems)
  }

  required <- c(
    "method", "state", "recommended", "quarantined", "api_enabled",
    "findings", "evidence"
  )
  for (field in required) {
    if (!(field %in% names(record))) {
      add_problem(
        "missing_field", field,
        paste0("Lifecycle metadata requires a `", field, "` field.")
      )
    }
  }

  if (!.alignment_scalar_character(method)) {
    add_problem(
      "invalid_method", "method", "`method` must be one non-empty string."
    )
  }

  state <- .alignment_record_field(record, "state")
  valid_state <- .alignment_scalar_character(state) &&
    state %in% policy$states$state
  if (!valid_state) {
    add_problem(
      "invalid_state", "state",
      paste0(
        "`state` must be one of: ",
        paste(policy$states$state, collapse = ", "), "."
      )
    )
  }

  for (field in c("recommended", "quarantined", "api_enabled")) {
    if (!.alignment_scalar_logical(.alignment_record_field(record, field))) {
      add_problem(
        "invalid_flag", field,
        paste0("`", field, "` must be TRUE or FALSE.")
      )
    }
  }

  findings <- .alignment_record_field(record, "findings", character())
  if (!is.character(findings) || anyNA(findings)) {
    add_problem(
      "invalid_findings", "findings", "`findings` must be a character vector."
    )
    findings <- character()
  }

  evidence <- .alignment_record_field(record, "evidence")
  if (!is.list(evidence)) {
    add_problem(
      "invalid_evidence", "evidence", "`evidence` must be a named list."
    )
  } else {
    for (dimension in names(policy$evidence_levels)) {
      level <- evidence[[dimension]]
      if (is.na(.alignment_level_rank(policy, dimension, level))) {
        add_problem(
          "invalid_evidence_level", paste0("evidence$", dimension),
          paste0(
            "`evidence$", dimension, "` must be one of: ",
            paste(policy$evidence_levels[[dimension]]$level, collapse = ", "),
            "."
          )
        )
      }
    }

    if (valid_state) {
      requirement <- policy$state_requirements[
        policy$state_requirements$state == state, , drop = FALSE
      ]
      for (dimension in names(policy$evidence_levels)) {
        observed <- .alignment_level_rank(policy, dimension, evidence[[dimension]])
        minimum <- .alignment_level_rank(
          policy, dimension, requirement[[dimension]][[1L]]
        )
        if (!is.na(observed) && observed < minimum) {
          add_problem(
            "insufficient_evidence", paste0("evidence$", dimension),
            paste0(
              "State `", state, "` requires ", dimension, " evidence at level `",
              requirement[[dimension]][[1L]], "` or higher."
            )
          )
        }
      }
    }

    if (identical(state, "stable")) {
      if (!.alignment_nonempty_character(evidence$oracle_ids)) {
        add_problem(
          "stable_oracle_missing", "evidence$oracle_ids",
          "Stable methods require reproducible exact or differential oracle artifacts."
        )
      }
      if (!.alignment_nonempty_character(evidence$benchmark_ids)) {
        add_problem(
          "stable_benchmark_missing", "evidence$benchmark_ids",
          "Stable methods require reproducible comparative accuracy and scale artifacts."
        )
      }
      if (!.alignment_scalar_character(evidence$reviewed_by)) {
        add_problem(
          "stable_review_missing", "evidence$reviewed_by",
          "Stable evidence requires a named reviewer."
        )
      }
      if (!.alignment_valid_date(evidence$reviewed_on)) {
        add_problem(
          "stable_review_missing", "evidence$reviewed_on",
          "Stable evidence requires an ISO review date (YYYY-MM-DD)."
        )
      }
    }
  }

  recommended <- .alignment_record_field(record, "recommended")
  quarantined <- .alignment_record_field(record, "quarantined")
  api_enabled <- .alignment_record_field(record, "api_enabled")
  if (valid_state && identical(recommended, TRUE)) {
    allowed <- policy$states$recommendation_allowed[
      match(state, policy$states$state)
    ]
    if (!isTRUE(allowed)) {
      add_problem(
        "state_not_recommendable", "recommended",
        paste0("A method in state `", state, "` cannot be presented as recommended.")
      )
    }
  }
  if (identical(quarantined, TRUE) && identical(recommended, TRUE)) {
    add_problem(
      "quarantined_recommended", "recommended",
      "A quarantined method cannot be presented as recommended."
    )
  }
  if (identical(state, "stable") && identical(quarantined, TRUE)) {
    add_problem(
      "stable_quarantined", "state",
      "An actively quarantined method cannot have stable status."
    )
  }
  if (identical(state, "removed") && identical(api_enabled, TRUE)) {
    add_problem(
      "removed_api_enabled", "api_enabled",
      "A removed method cannot expose an enabled public API."
    )
  }

  trigger_rows <- policy$quarantine_triggers[
    policy$quarantine_triggers$quarantine_required, , drop = FALSE
  ]
  triggered <- intersect(findings, trigger_rows$finding)
  if (length(triggered) > 0L && identical(quarantined, FALSE)) {
    add_problem(
      "quarantine_required", "quarantined",
      paste0(
        "Finding(s) require immediate quarantine: ",
        paste(triggered, collapse = ", "), "."
      )
    )
  }

  if (identical(quarantined, TRUE)) {
    quarantine <- .alignment_record_field(record, "quarantine")
    if (!is.list(quarantine)) {
      add_problem(
        "quarantine_metadata_missing", "quarantine",
        "An active quarantine requires a `quarantine` metadata list."
      )
    } else {
      for (field in policy$quarantine_required_fields) {
        value <- quarantine[[field]]
        valid <- if (identical(field, "triggers")) {
          .alignment_nonempty_character(value)
        } else {
          .alignment_scalar_character(value)
        }
        if (!valid) {
          add_problem(
            "quarantine_metadata_missing", paste0("quarantine$", field),
            paste0("Active quarantine requires a non-empty `", field, "` field.")
          )
        }
      }
      if (.alignment_nonempty_character(quarantine$triggers)) {
        unknown <- setdiff(
          quarantine$triggers,
          trigger_rows$finding
        )
        if (length(unknown) > 0L) {
          add_problem(
            "invalid_quarantine_trigger", "quarantine$triggers",
            paste0("Unknown quarantine trigger(s): ", paste(unknown, collapse = ", "), ".")
          )
        }
        if (length(setdiff(triggered, quarantine$triggers)) > 0L) {
          add_problem(
            "quarantine_trigger_unrecorded", "quarantine$triggers",
            "All immediate findings must be recorded in `quarantine$triggers`."
          )
        }
        if (length(setdiff(quarantine$triggers, findings)) > 0L) {
          add_problem(
            "quarantine_trigger_unsubstantiated", "quarantine$triggers",
            "Each quarantine trigger must also appear in the method's `findings`."
          )
        }
        disabling <- trigger_rows$finding[trigger_rows$disable_api]
        if (any(quarantine$triggers %in% disabling) && identical(api_enabled, TRUE)) {
          add_problem(
            "unsafe_api_enabled", "api_enabled",
            "This quarantine trigger requires the affected public API to be disabled."
          )
        }
      }
      if (!is.null(quarantine$repair_by) &&
          !.alignment_valid_date(quarantine$repair_by)) {
        add_problem(
          "invalid_repair_deadline", "quarantine$repair_by",
          "`quarantine$repair_by` must be an ISO date (YYYY-MM-DD)."
        )
      }
    }
  }

  if (identical(state, "deprecated") || identical(state, "removed")) {
    .alignment_check_deprecation(record, method, policy, add_problem)
  }

  if (!is.null(previous_state)) {
    if (!.alignment_scalar_character(previous_state) ||
        !(previous_state %in% policy$states$state)) {
      add_problem(
        "invalid_previous_state", "previous_state",
        "`previous_state` must be a recognized lifecycle state."
      )
    } else if (valid_state && !identical(previous_state, state)) {
      allowed <- any(
        policy$transitions$from == previous_state &
          policy$transitions$to == state
      )
      if (!allowed) {
        add_problem(
          "invalid_transition", "state",
          paste0(
            "Lifecycle transition `", previous_state, "` -> `", state,
            "` is not permitted."
          )
        )
      }
    }
  }

  if (isTRUE(previous_quarantined) && identical(state, "stable") &&
      identical(quarantined, FALSE)) {
    restoration <- .alignment_record_field(record, "restoration")
    if (!is.list(restoration)) {
      add_problem(
        "restoration_evidence_missing", "restoration",
        "Restoration from quarantine to stable requires a `restoration` evidence list."
      )
    } else {
      for (field in policy$restoration_required_fields) {
        value <- restoration[[field]]
        valid <- if (field %in% c(
          "resolved_triggers", "oracle_ids", "benchmark_ids"
        )) {
          .alignment_nonempty_character(value)
        } else if (identical(field, "reviewed_on")) {
          .alignment_valid_date(value)
        } else {
          .alignment_scalar_character(value)
        }
        if (!valid) {
          add_problem(
            "restoration_evidence_missing", paste0("restoration$", field),
            paste0("Stable restoration requires a valid `", field, "` field.")
          )
        }
      }
    }
  }

  class(problems) <- c("alignment_lifecycle_problems", "data.frame")
  rownames(problems) <- NULL
  problems
}

.alignment_records <- function(records) {
  if (is.list(records) && all(c("method", "state") %in% names(records))) {
    return(list(records))
  }
  if (!is.list(records)) {
    stop("`records` must be one lifecycle record or a list of records.", call. = FALSE)
  }
  records
}

.alignment_release_condition <- function(problems) {
  detail <- paste0(
    "- ", problems$method, " [", problems$code, "] ", problems$message,
    collapse = "\n"
  )
  structure(
    list(
      message = paste0("Alignment lifecycle release gate failed:\n", detail),
      call = NULL,
      problems = problems
    ),
    class = c("alignment_release_error", "error", "condition")
  )
}

#' @rdname alignment_lifecycle_policy
#' @export
assert_alignment_release_ready <- function(
  records,
  policy = alignment_lifecycle_policy()
) {
  records <- .alignment_records(records)
  problems <- .alignment_empty_problems()
  methods <- character(length(records))
  for (index in seq_along(records)) {
    problems <- rbind(
      problems,
      validate_alignment_lifecycle(records[[index]], policy = policy)
    )
    record_method <- .alignment_record_field(
      records[[index]], "method", "<unknown>"
    )
    methods[[index]] <- if (.alignment_scalar_character(record_method)) {
      record_method
    } else {
      "<unknown>"
    }
  }
  duplicated_methods <- unique(methods[duplicated(methods)])
  for (method in duplicated_methods) {
    problems <- rbind(
      problems,
      .alignment_lifecycle_problem(
        method, "duplicate_method", "method",
        "Lifecycle release metadata must contain exactly one record per method."
      )
    )
  }
  class(problems) <- c("alignment_lifecycle_problems", "data.frame")
  rownames(problems) <- NULL
  if (nrow(problems) > 0L) stop(.alignment_release_condition(problems))
  invisible(TRUE)
}

.alignment_deprecation_message <- function(record, package = "manifoldalign") {
  method <- .alignment_record_field(record, "method", "This method")
  deprecation <- .alignment_record_field(record, "deprecation", list())
  replacement <- .alignment_or(
    deprecation$replacement, "a supported replacement"
  )
  removal_at <- .alignment_or(
    deprecation$removal_at, "the documented removal release"
  )
  migration <- .alignment_or(
    deprecation$migration,
    "See the lifecycle documentation for migration guidance."
  )
  status <- if (identical(.alignment_record_field(record, "api_enabled"), FALSE)) {
    " is disabled and deprecated"
  } else {
    " is deprecated"
  }
  paste0(
    "`", method, "()`", status, "; use `", replacement, "()` instead. ",
    "Removal is scheduled for ", package, " ", removal_at, ". ", migration
  )
}
