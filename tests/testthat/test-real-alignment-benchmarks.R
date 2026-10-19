library(testthat)
library(manifoldalign)

test_that("real benchmark registry freezes three complementary tasks", {
  registry <- manifoldalign:::real_alignment_benchmark_registry()

  expect_identical(registry$protocol_version, rep("1.0.0", 3L))
  expect_setequal(registry$task_id, c(
    "so2sat_lcz42_sar_optical_culture10",
    "aucs_coauthor_lunch",
    "monastery_like_dislike"
  ))
  expect_setequal(registry$stress_regime, c(
    "modality_conflict", "attributes_informative", "topology_informative"
  ))
  expect_true(all(nzchar(registry$source_version)))
  expect_true(all(nzchar(registry$license)))
  expect_true(all(grepl("^https://", registry$source_url)))
  expect_identical(
    registry$split_id[registry$dataset_id == "so2sat_lcz42"],
    "culture10-v4.2-blocked"
  )
})

test_that("real artifact registry carries exact provenance", {
  artifacts <- manifoldalign:::real_alignment_artifact_registry()

  expect_identical(
    anyDuplicated(paste(artifacts$dataset_id, artifacts$artifact_id)), 0L
  )
  expect_true(all(is.finite(artifacts$bytes) & artifacts$bytes > 0))
  expect_true(all(grepl("^[0-9a-f]{64}$", artifacts$sha256)))
  expect_true(all(nzchar(artifacts$source_revision)))
  expect_true(all(nzchar(artifacts$license)))

  source <- artifacts[artifacts$artifact_id == "multinet_source", ]
  expect_equal(source$bytes, 2251834)
  expect_identical(
    source$sha256,
    "732128032fb9d2ef6b04d9577553eb08b5dbe7078a3681b1c4a53f0edb07039a"
  )

  for (dataset_id in c("aucs", "sampson_monastery")) {
    row <- artifacts[
      artifacts$dataset_id == dataset_id &
        artifacts$artifact_id == "packaged_fixture",
      , drop = FALSE
    ]
    expect_equal(nrow(row), 1L)
    path <- manifoldalign:::.real_network_fixture_file(dataset_id)
    expect_equal(as.numeric(file.info(path)$size), row$bytes)
    expect_identical(manifoldalign:::.real_sha256_file(path), row$sha256)
  }

  so2sat <- artifacts[
    artifacts$dataset_id == "so2sat_lcz42" &
      artifacts$artifact_id == "ci_feature_fixture",
    , drop = FALSE
  ]
  fixture <- manifoldalign:::.real_benchmark_inst_file(
    file.path("extdata", "alignment-real"),
    "so2sat_lcz42_ci_features.csv"
  )
  expect_equal(nrow(so2sat), 1L)
  expect_equal(as.numeric(file.info(fixture)$size), so2sat$bytes)
  expect_identical(manifoldalign:::.real_sha256_file(fixture), so2sat$sha256)

  full <- artifacts[
    artifacts$dataset_id == "so2sat_lcz42" &
      artifacts$artifact_id == "full_feature_receipt",
    , drop = FALSE
  ]
  expect_equal(nrow(full), 1L)
  expect_equal(full$bytes, 19531044)
  expect_identical(
    full$sha256,
    "9eeb19fce6bb45b1fd4d3c685d03b6f4b8704dd6ec1cd6de059c7aa699fdce78"
  )
})

test_that("real network cases preserve observed structure and exact identities", {
  expected <- list(
    aucs_coauthor_lunch = c(25L, 60L, 25L, 21L, 193L),
    monastery_like_dislike = c(18L, 17L, 17L, 55L, 47L)
  )

  for (task_id in names(expected)) {
    case <- manifoldalign:::real_alignment_network_case(task_id)
    d1 <- case$domains[[1L]]
    d2 <- case$domains[[2L]]
    observed <- c(
      nrow(d1$X), nrow(d2$X), sum(!is.na(case$truth$map12)),
      nrow(d1$relations), nrow(d2$relations)
    )

    expect_s3_class(case, "real_graph_alignment_case")
    expect_identical(observed, expected[[task_id]])
    expect_gt(sum(d1$private_mask) + sum(d2$private_mask), 0L)
    expect_false(identical(d1$row_id, d2$row_id))
    expect_identical(
      d1$row_id[case$truth$correspondence$source_row],
      d2$row_id[case$truth$correspondence$target_row]
    )
    expect_false(any(grepl("actor|(^|_)id($|_)", colnames(d1$X))))
    expect_false(any(grepl("actor|(^|_)id($|_)", colnames(d2$X))))
    expect_true(all(is.finite(d1$relations$weight)))
    expect_true(all(is.finite(d2$relations$weight)))
    expect_equal(nrow(manifoldalign:::validate_real_alignment_case(case)), 0L)
  }

  aucs <- manifoldalign:::real_alignment_network_case("aucs_coauthor_lunch")
  monastery <- manifoldalign:::real_alignment_network_case(
    "monastery_like_dislike"
  )
  expect_false(any(aucs$domains[[1L]]$relations$directed))
  expect_false(any(aucs$domains[[2L]]$relations$directed))
  expect_true(all(monastery$domains[[1L]]$relations$directed))
  expect_true(all(monastery$domains[[2L]]$relations$directed))
  expect_setequal(monastery$domains[[1L]]$relations$weight, 1:3)
  expect_setequal(monastery$domains[[2L]]$relations$weight, 1:3)
})

test_that("real network validator catches identity, label, and edge leakage", {
  case <- manifoldalign:::real_alignment_network_case("aucs_coauthor_lunch")

  wrong_map <- case
  shared <- which(!is.na(wrong_map$truth$map12))
  wrong_map$truth$map12[shared[[1L]]] <- wrong_map$truth$map12[shared[[2L]]]
  expect_true("map_identity" %in%
    manifoldalign:::validate_real_alignment_case(wrong_map)$code)

  leaked_label <- case
  final <- which(leaked_label$domains[[1L]]$node_role == "final")[[1L]]
  leaked_label$domains[[1L]]$labels[[final]] <- "leaked"
  expect_true("final_label_leakage" %in%
    manifoldalign:::validate_real_alignment_case(leaked_label)$code)

  leaked_id <- case
  leaked_id$domains[[1L]]$X <- cbind(
    leaked_id$domains[[1L]]$X,
    actor_id = seq_len(nrow(leaked_id$domains[[1L]]$X))
  )
  expect_true("identity_leakage" %in%
    manifoldalign:::validate_real_alignment_case(leaked_id)$code)

  bad_edge <- case
  bad_edge$domains[[1L]]$relations$from[[1L]] <- 0L
  expect_true("relations" %in%
    manifoldalign:::validate_real_alignment_case(bad_edge)$code)
})

test_that("real network controls are complete deterministic diagnostics", {
  expected_controls <- c(
    "raw_attribute_knn", "topology_descriptor_only", "position_only",
    "row_index", "chance", "target_only_classifier"
  )
  for (task_id in c("aucs_coauthor_lunch", "monastery_like_dislike")) {
    case <- manifoldalign:::real_alignment_network_case(task_id)
    first <- manifoldalign:::benchmark_real_network_controls(case)
    second <- manifoldalign:::benchmark_real_network_controls(case)

    expect_setequal(first$control_id, expected_controls)
    expect_identical(first, second)
    expect_true(all(first$status[first$control_id %in% c(
      "raw_attribute_knn", "topology_descriptor_only", "row_index", "chance"
    )] == "ok"))
    expect_true(first$observed_topology_used[
      first$control_id == "topology_descriptor_only"
    ])
    expect_identical(
      first$status_code[first$control_id == "position_only"],
      "positions_not_available"
    )
    expect_true(all(is.finite(first$top1_accuracy[first$status == "ok"])))
    expect_true(all(nzchar(first$prediction_hash[first$status == "ok"])))
  }
})

test_that("So2Sat CI fixture enforces blocked inductive leakage boundaries", {
  case <- manifoldalign:::real_so2sat_ci_case()

  expect_s3_class(case, "real_remote_sensing_case")
  expect_identical(case$governance$evidence_scope, "ci_schema_only")
  expect_identical(case$governance$split_id, "culture10-v4.2-blocked")
  expect_false("final_target_labels" %in% names(case$fit))
  expect_identical(
    case$fit$final_target_row_id,
    case$evaluation$final_target_row_id
  )
  expect_length(intersect(
    case$fit$source_row_id, case$evaluation$final_target_row_id
  ), 0L)
  expect_identical(dim(case$fit$source_X), dim(case$fit$target_X))
  expect_equal(ncol(case$fit$source_X), ncol(case$fit$final_target_X))
  expect_true(all(is.finite(case$fit$source_X)))
  expect_true(all(is.finite(case$fit$target_X)))
  expect_true(all(is.finite(case$fit$final_target_X)))
  expect_equal(case$governance$n_development, 68L)
  expect_equal(case$governance$n_final, 68L)

  controls <- manifoldalign:::benchmark_real_so2sat_controls(case)
  expect_setequal(controls$control_id, c(
    "raw_attribute_source_only", "position_only",
    "topology_descriptor_only", "target_only_classifier", "chance"
  ))
  expect_identical(
    controls$status_code[controls$control_id == "topology_descriptor_only"],
    "observed_graph_not_available"
  )
  expect_true(all(is.finite(controls$accuracy[controls$status == "ok"])))
  expect_true(all(nzchar(controls$prediction_hash[controls$status == "ok"])))
  expect_identical(
    controls,
    manifoldalign:::benchmark_real_so2sat_controls(case)
  )
})

test_that("full So2Sat loader cannot silently masquerade as CI evidence", {
  ci <- manifoldalign:::real_so2sat_ci_case()
  copy <- tempfile(fileext = ".csv")
  file.copy(ci$governance$feature_artifact, copy)

  full <- manifoldalign:::real_so2sat_full_case(copy)
  expect_identical(full$governance$evidence_scope, "full_frozen_split")
  expect_identical(full$governance$feature_artifact_sha256,
    ci$governance$feature_artifact_sha256
  )
  expect_false("final_target_labels" %in% names(full$fit))
})

test_that("dataset and algorithm failures remain distinguishable", {
  unavailable <- manifoldalign:::real_alignment_status_row(
    "task", "input", "dataset_unavailable", "missing"
  )
  invalid <- manifoldalign:::real_alignment_status_row(
    "task", "input", "dataset_invalid", "checksum"
  )
  algorithm <- manifoldalign:::real_alignment_status_row(
    "task", "fit", "algorithm_error", "solver"
  )
  ok <- manifoldalign:::real_alignment_status_row(
    "task", "score", "ok", "complete"
  )

  expect_identical(unavailable$failure_class, "dataset")
  expect_identical(invalid$failure_class, "dataset")
  expect_identical(algorithm$failure_class, "algorithm")
  expect_true(is.na(ok$failure_class))
})
