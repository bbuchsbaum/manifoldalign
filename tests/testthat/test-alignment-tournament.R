make_tournament_perfect_ranking <- function(case, prediction) {
  n_target <- nrow(case$domains[[2L]]$X)
  ranking <- matrix(NA_integer_, length(prediction), n_target)
  for (row in which(!is.na(prediction))) {
    ranking[row, ] <- c(
      prediction[[row]], setdiff(seq_len(n_target), prediction[[row]])
    )
  }
  list(index = ranking, complete = TRUE)
}

permute_tournament_target <- function(case, permutation) {
  inverse <- match(seq_along(permutation), permutation)
  target <- case$domains[[2L]]
  row_fields <- c(
    "X", "positions", "attributes", "labels", "row_id", "node_role",
    "shared_mask", "private_mask"
  )
  for (field in row_fields) {
    value <- target[[field]]
    target[[field]] <- if (is.data.frame(value) || is.matrix(value)) {
      value[permutation, , drop = FALSE]
    } else {
      value[permutation]
    }
  }
  target$relations$from <- inverse[target$relations$from]
  target$relations$to <- inverse[target$relations$to]
  case$domains[[2L]] <- target
  mapped <- case$truth$map12
  mapped[!is.na(mapped)] <- inverse[mapped[!is.na(mapped)]]
  case$truth$map12 <- mapped
  case$truth$target_shared_mask <- target$shared_mask
  case$truth$target_private_mask <- target$private_mask
  case$truth$domain_labels[[2L]] <- case$truth$domain_labels[[2L]][permutation]
  case$anchors$target_row <- inverse[case$anchors$target_row]
  case
}

test_that("tournament roster and resource envelopes are frozen and complete", {
  methods <- manifoldalign:::.alignment_tournament_methods_v2()
  envelope <- manifoldalign:::.alignment_tournament_resource_ceiling_v2()

  expect_length(methods, 14L)
  expect_identical(anyDuplicated(methods), 0L)
  expect_setequal(envelope$method_id, methods)
  expect_identical(anyDuplicated(envelope$method_id), 0L)
  expect_true(all(envelope$max_total_nodes >= 0L))
  expect_identical(
    envelope$max_total_nodes[envelope$method_id == "linear_sim_embed"],
    0L
  )
})

test_that("applicability separates semantic and resource refusals", {
  partial <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-small-accuracy-sbm-01", "final"
  )
  closed_medium <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-medium-robustness-sbm-closed-01", "final"
  )

  expect_true(
    manifoldalign:::.alignment_tournament_applicability_v2(
      "rapid_ma", partial
    )$intended_to_run
  )
  token <- manifoldalign:::.alignment_tournament_applicability_v2(
    "token_ot_graph_align", partial
  )
  expect_false(token$intended_to_run)
  expect_identical(token$status_code, "unsupported_partial_overlap")

  fpgw <- manifoldalign:::.alignment_tournament_applicability_v2(
    "fpgw", closed_medium
  )
  expect_false(fpgw$intended_to_run)
  expect_identical(fpgw$status_code, "unsupported_resource_envelope")

  lse <- manifoldalign:::.alignment_tournament_applicability_v2(
    "linear_sim_embed", partial
  )
  expect_false(lse$intended_to_run)
  expect_identical(lse$status, "not_applicable")
})

test_that("perfect correspondence scores every matching partition exactly", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-small-accuracy-sbm-01", "final"
  )
  prediction <- case$truth$map12
  score <- manifoldalign:::.alignment_tournament_score_matching_v2(
    case,
    prediction,
    make_tournament_perfect_ranking(case, prediction),
    "independent_oracle"
  )

  expect_equal(score$matching_n_eligible, case$config$n_shared)
  expect_equal(score$coverage, 1)
  expect_equal(score$top1_accuracy, 1)
  expect_equal(score$hits_at_1, 1)
  expect_equal(score$hits_at_5, 1)
  expect_equal(score$hits_at_10, 1)
  expect_equal(score$mrr, 1)
  expect_equal(score$mrr_at_10, 1)
  expect_equal(score$anchor_top1_accuracy, 1)
  expect_equal(score$nonanchor_top1_accuracy, 1)
  expect_equal(score$private_rejection_accuracy, 1)
  expect_true(score$one_to_one_feasible)
  expect_equal(score$one_to_one_rate, 1)
  expect_gte(score$weighted_structure_conservation, 0)
  expect_lte(score$weighted_structure_conservation, 1)
})

test_that("missing and duplicate predictions cannot inflate accuracy", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-small-accuracy-sbm-01", "final"
  )
  prediction <- case$truth$map12
  shared <- which(!is.na(prediction))
  prediction[shared[[1L]]] <- NA_integer_
  prediction[shared[[2L]]] <- prediction[shared[[3L]]]
  score <- manifoldalign:::.alignment_tournament_score_matching_v2(
    case, prediction
  )

  expect_equal(score$matching_n_eligible, case$config$n_shared)
  expect_equal(score$coverage, (case$config$n_shared - 1L) / case$config$n_shared)
  expect_equal(score$top1_accuracy, (case$config$n_shared - 2L) / case$config$n_shared)
  expect_false(score$one_to_one_feasible)
  expect_lt(score$one_to_one_rate, 1)
  expect_error(
    manifoldalign:::.alignment_tournament_score_matching_v2(
      case, prediction[-1L]
    ),
    "expected"
  )
})

test_that("target row permutation leaves oracle metrics invariant", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-small-accuracy-geometric-closed-01", "final"
  )
  prediction <- case$truth$map12
  baseline <- manifoldalign:::.alignment_tournament_score_matching_v2(
    case, prediction, make_tournament_perfect_ranking(case, prediction)
  )
  set.seed(812L)
  permutation <- sample.int(nrow(case$domains[[2L]]$X))
  permuted <- permute_tournament_target(case, permutation)
  changed_prediction <- permuted$truth$map12
  changed <- manifoldalign:::.alignment_tournament_score_matching_v2(
    permuted,
    changed_prediction,
    make_tournament_perfect_ranking(permuted, changed_prediction)
  )

  fields <- c(
    "coverage", "top1_accuracy", "hits_at_1", "hits_at_5",
    "hits_at_10", "mrr", "one_to_one_rate",
    "weighted_structure_conservation"
  )
  expect_equal(baseline[fields], changed[fields], tolerance = 1e-12)
})

test_that("classification probe reports calibrated multi-class metrics", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-small-accuracy-sbm-closed-01", "final"
  )
  embedding <- list(
    source = case$domains[[1L]]$X,
    target = case$domains[[2L]]$X
  )
  probabilities <- manifoldalign:::.alignment_tournament_centroid_probabilities_v2(
    case, embedding
  )
  score <- manifoldalign:::.alignment_tournament_score_classification_v2(
    case, probabilities
  )

  expect_identical(score$classification_status, "ok")
  bounded <- unlist(score[c(
    "classification_coverage", "classification_oa",
    "classification_balanced_accuracy", "classification_macro_f1",
    "classification_miou", "classification_ece"
  )])
  expect_true(all(is.finite(bounded)))
  expect_true(all(bounded >= 0 & bounded <= 1))
  expect_true(is.finite(score$classification_nll))
  expect_gte(score$classification_nll, 0)
})

test_that("every applicable package adapter completes the tiny oracle contract", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-sbm-01", "final"
  )
  rows <- lapply(
    manifoldalign:::.alignment_tournament_methods_v2(),
    function(method) manifoldalign:::.alignment_tournament_run_case_v2(
      case, method
    )
  )
  results <- do.call(rbind, rows)
  applicable <- results$intended_to_run

  expect_true(all(results$status[applicable] == "ok"))
  expect_true(all(is.finite(results$top1_accuracy[applicable])))
  expect_true(all(results$matching_status[applicable] == "ok"))
  expect_true(all(results$parameter_trials_used[applicable] == 1L))
  expect_identical(
    results$status[results$method_id == "linear_sim_embed"],
    "not_applicable"
  )
})

test_that("intent-to-run summary keeps failures in its denominator", {
  case_a <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-sbm-01", "final"
  )
  case_b <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-sbm-02", "final"
  )
  success <- manifoldalign:::.alignment_tournament_base_row_v2(
    case_a, "rapid_ma"
  )
  success$status <- "ok"
  success$matching_status <- "ok"
  success$top1_accuracy <- 1
  success$intended_to_run <- TRUE
  success$runtime_sec <- 1
  success$peak_rss_mb <- 10
  failure <- manifoldalign:::.alignment_tournament_failure_row_v2(
    case_b, "rapid_ma", "timeout", "wall_clock_budget_exceeded",
    "synthetic timeout", runtime_sec = 30, peak_rss_mb = 12
  )
  results <- rbind(success, failure)
  summary <- manifoldalign:::.alignment_tournament_summary_by_regime_v2(results)

  expect_equal(summary$n_intended, 2L)
  expect_equal(summary$n_ok, 1L)
  expect_equal(summary$conditional_top1, 1)
  expect_equal(summary$intent_top1, 0.5)
  expect_equal(summary$completion_rate, 0.5)
})

test_that("grid validator catches omissions and accepts complete typed rows", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-sbm-01", "final"
  )
  methods <- c("rapid_ma", "linear_sim_embed")
  results <- do.call(rbind, lapply(methods, function(method) {
    manifoldalign:::.alignment_tournament_run_case_v2(case, method)
  }))
  manifest <- manifoldalign:::graph_alignment_benchmark_manifest_v2("final")
  manifest <- manifest[manifest$case_id == case$case_id, , drop = FALSE]

  problems <- manifoldalign:::.alignment_tournament_validate_results_v2(
    results, manifest, methods, require_external_memory = FALSE
  )
  expect_equal(nrow(problems), 0L)

  incomplete <- manifoldalign:::.alignment_tournament_validate_results_v2(
    results[-1L, , drop = FALSE], manifest, methods,
    require_external_memory = FALSE
  )
  expect_true("missing_rows" %in% incomplete$code)
})
