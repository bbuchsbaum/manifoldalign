library(testthat)
library(manifoldalign)

test_that("graph benchmark registry enumerates scenario metadata cleanly", {
  reg <- manifoldalign:::synthetic_graph_alignment_registry(
    sizes = c(12L, 18L),
    structures = c("ring", "grid"),
    d = 4L,
    noise_sd = 0.02,
    permute_fraction = 0.5,
    n_anchors = 3L,
    n_reps = 2L,
    seed = 11L
  )

  expect_true(is.data.frame(reg))
  expect_equal(nrow(reg), 8)
  expect_true(all(c(
    "scenario_family", "scenario", "structure", "n", "d", "noise_sd",
    "permute_fraction", "n_anchors", "rep", "generation_seed"
  ) %in% names(reg)))
  expect_true(all(reg$scenario_family == "graph"))
  expect_true(all(reg$structure %in% c("ring", "grid")))
})

test_that("graph benchmark case returns consistent maps and anchor correspondences", {
  case <- manifoldalign:::synthetic_graph_alignment_case(
    n = 20L,
    d = 3L,
    structure = "community",
    noise_sd = 0.01,
    permute_fraction = 0.6,
    n_anchors = 4L,
    seed = 22L
  )

  expect_equal(dim(case$X1), c(20, 3))
  expect_equal(dim(case$X2), c(20, 3))
  expect_equal(length(case$map12), 20)
  expect_equal(sort(case$map12), seq_len(20))
  expect_equal(length(case$anchor_idx), 4)
  expect_true(is.data.frame(case$correspondences))
  expect_equal(nrow(case$correspondences), 4)
  expect_s3_class(case$hd, "hyperdesign")
  expect_identical(case$benchmark_meta$structure, "community")
})

test_that("graph benchmark output includes shared scenario metadata columns", {
  res <- benchmark_graph_alignment_methods(
    sizes = 18L,
    d = 3L,
    noise_sd = 0.01,
    structure = "grid",
    permute_fraction = 1,
    n_anchors = 3L,
    methods = "token_ot_graph",
    candidate_k = 8L,
    n_levels = 1L,
    n_reps = 1L,
    seed = 7L,
    verbose = FALSE
  )

  expect_true(is.data.frame(res))
  expect_true(all(c(
    "method", "method_id", "input_class", "observed_topology_used",
    "adapter_id", "topology_disposition", "observed_edge_count",
    "adapter_input_hash", "row_identity_hash", "row_identity_preserved",
    "prediction_hash", "search_space_id", "search_space_hash",
    "parameter_search_budget", "parameter_trials_used",
    "wall_clock_budget_sec", "peak_rss_mb", "memory_measurement",
    "protocol_version", "benchmark_case_id",
    "benchmark_case_hash", "benchmark_split", "benchmark_tier",
    "benchmark_family", "preprocessing_id", "scenario_family", "scenario",
    "structure", "n", "d", "noise_sd", "permute_fraction", "n_anchors", "rep",
    "generation_seed", "decode_mode", "runtime_sec", "top1_accuracy",
    "coverage", "n_eligible", "anchor_top1_accuracy",
    "nonanchor_top1_accuracy", "private_rejection_accuracy", "converged",
    "status", "status_code", "status_detail", "error"
  ) %in% names(res)))
  expect_true(all(res$scenario_family == "attributed_graph_v2"))
  expect_true(all(res$structure == "grid"))
  expect_true(all(res$decode_mode == "native"))
  expect_true(all(res$method_id == "token_ot_graph_align"))
  expect_true(all(res$input_class == "feature_derived_graph"))
  expect_true(all(res$topology_disposition == "feature_derived_only"))
  expect_false(any(res$observed_topology_used))
  expect_true(all(res$observed_edge_count == 0L))
  expect_true(all(res$protocol_version == "2.0.0"))
  expect_true(all(res$benchmark_family == "geometric"))
  expect_true(all(res$preprocessing_id == "train_numeric_center_scale_v1"))
  expect_true(all(nzchar(res$benchmark_case_hash)))
  expect_true(all(res$status == "ok"))
  expect_true(all(res$status_code == "ok"))
  expect_true(all(is.na(res$status_detail)))
  expect_true(all(res$row_identity_preserved))
  expect_true(all(nzchar(res$adapter_input_hash)))
  expect_true(all(nzchar(res$row_identity_hash)))
  expect_true(all(nzchar(res$prediction_hash)))
  expect_identical(res$parameter_trials_used, 1L)
  expect_identical(res$memory_measurement, "external_process_required")
})

test_that("graph benchmark can sweep multiple structures in one call", {
  res <- benchmark_graph_alignment_methods(
    sizes = 16L,
    d = 3L,
    noise_sd = 0.01,
    structure = c("ring", "community"),
    permute_fraction = 1,
    n_anchors = 3L,
    methods = "token_ot_graph",
    candidate_k = 8L,
    n_levels = 1L,
    n_reps = 1L,
    seed = 17L,
    verbose = FALSE
  )

  expect_true(is.data.frame(res))
  expect_setequal(unique(res$structure), c("ring", "community"))
  expect_equal(nrow(res), 2L)
  expect_true(all(grepl("^(ring|community)_16_1$", res$scenario)))
})

test_that("RAPID-MA benchmark rows prove that supplied v2 edges were consumed", {
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.01,
    structure = "sbm",
    n_anchors = 3L,
    methods = "rapid_ma",
    candidate_k = 8L,
    n_reps = 1L,
    seed = 19L,
    verbose = FALSE
  )

  expect_identical(res$method_id, "rapid_ma")
  expect_identical(res$input_class, "supplied_graph")
  expect_identical(res$topology_disposition, "consumed")
  expect_true(res$observed_topology_used)
  expect_gt(res$observed_edge_count, 0L)
  expect_identical(res$scenario_family, "attributed_graph_v2")
  expect_identical(res$protocol_version, "2.0.0")
  expect_identical(res$benchmark_family, "sbm")
  expect_identical(res$status, "ok")
  expect_true(is.finite(res$top1_accuracy))
  expect_true(is.na(res$error))
})

test_that("interactive dispatcher refuses non-frozen partial permutation", {
  expect_error(
    benchmark_graph_alignment_methods(
      sizes = 8L,
      structure = "sbm",
      permute_fraction = 0.5,
      methods = "rapid_ma",
      n_reps = 1L
    ),
    "protocol v2 requires `permute_fraction = 1`",
    fixed = TRUE
  )
})

test_that("methods requiring anchors report unsupported without running solvers", {
  res <- benchmark_graph_alignment_methods(
    sizes = 8L,
    d = 2L,
    structure = "ring",
    n_anchors = 0L,
    methods = c("parrot", "lra", "gpca"),
    n_reps = 1L,
    seed = 31L,
    verbose = FALSE
  )

  expect_equal(nrow(res), 3L)
  expect_setequal(
    res$method_id,
    c("parrot", "lowrank_align", "gpca_align")
  )
  expect_true(all(res$status == "unsupported"))
  expect_true(all(res$status_code == "requires_anchors"))
  expect_true(all(is.na(res$runtime_sec)))
  expect_true(all(is.na(res$top1_accuracy)))
  expect_true(all(res$coverage == 0))
  expect_true(all(is.na(res$error)))
  expect_true(all(nzchar(res$status_detail)))
})

test_that("FPGW uses the compiled transport backend without an lpSolve gate", {
  res <- benchmark_graph_alignment_methods(
    sizes = 8L,
    d = 2L,
    structure = "random",
    methods = "fpgw",
    n_reps = 1L,
    seed = 32L,
    verbose = FALSE
  )

  expect_identical(res$method_id, "fpgw")
  expect_identical(res$input_class, "metric_measure")
  expect_identical(res$topology_disposition, "not_applicable")
  expect_false(res$observed_topology_used)
  expect_identical(res$status, "ok")
  expect_identical(res$status_code, "ok")
  expect_true(is.finite(res$runtime_sec))
  expect_true(is.finite(res$top1_accuracy))
  expect_true(is.na(res$error))
  expect_true(res$row_identity_preserved)
  expect_false(is.na(res$converged))
})

test_that("supplied-edge mutations change only native topology adapter inputs", {
  baseline <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 24L,
    family = "sbm",
    attr_dim = 4L,
    n_anchors = 4L,
    directed = TRUE,
    seed = 3201L,
    case_id = "adapter-edge-baseline"
  )
  changed <- baseline
  changed$domains[[2L]]$relations$weight[[1L]] <-
    1.5 * changed$domains[[2L]]$relations$weight[[1L]]

  methods <- c(
    "rapid_ma", "token_ot_graph", "fpgw", "gw", "spectral_mnn",
    "ssma", "kema", "global_geo", "lra", "gpca", "cone_align",
    "grasp", "parrot", "lse"
  )
  before <- vapply(methods, function(method) {
    manifoldalign:::.graph_benchmark_adapter_contract_v2(
      method, baseline
    )$adapter_input_hash
  }, character(1))
  after <- vapply(methods, function(method) {
    manifoldalign:::.graph_benchmark_adapter_contract_v2(
      method, changed
    )$adapter_input_hash
  }, character(1))

  expect_identical(names(before), methods)
  expect_identical(unname(before != after), methods == "rapid_ma")
})

test_that("topology preflight rejects unsupported signed directed and weighted inputs", {
  info <- manifoldalign:::alignment_method_info("rapid_ma")

  signed <- manifoldalign:::.graph_benchmark_validate_topology_contract(
    info,
    list(signed = TRUE, weighted = TRUE, directed = FALSE, multirelation = FALSE)
  )
  expect_false(signed$supported)
  expect_identical(signed$status_code, "unsupported_signed_weights")

  no_direction <- info
  no_direction$topology$directed <- "unsupported"
  directed <- manifoldalign:::.graph_benchmark_validate_topology_contract(
    no_direction,
    list(signed = FALSE, weighted = TRUE, directed = TRUE, multirelation = FALSE)
  )
  expect_false(directed$supported)
  expect_identical(directed$status_code, "unsupported_directed_topology")

  no_weights <- info
  no_weights$topology$nonnegative_weights <- "unsupported"
  weighted <- manifoldalign:::.graph_benchmark_validate_topology_contract(
    no_weights,
    list(signed = FALSE, weighted = TRUE, directed = FALSE, multirelation = FALSE)
  )
  expect_false(weighted$supported)
  expect_identical(weighted$status_code, "unsupported_weighted_topology")
})

test_that("normalized adapter predictions preserve row identities and partitions", {
  case <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 20L,
    n_private = c(4L, 3L),
    family = "geometric",
    n_anchors = 5L,
    seed = 3202L,
    case_id = "adapter-identity-contract"
  )
  prediction <- case$truth$map12
  normalized <- manifoldalign:::.graph_benchmark_normalize_prediction(
    case, prediction
  )

  expect_true(normalized$ok)
  expect_true(normalized$row_identity_preserved)
  expect_equal(normalized$n_eligible, 20L)
  expect_equal(normalized$top1_accuracy, 1)
  expect_equal(normalized$coverage, 1)
  expect_equal(normalized$anchor_top1_accuracy, 1)
  expect_equal(normalized$nonanchor_top1_accuracy, 1)
  expect_equal(normalized$private_rejection_accuracy, 1)
  expect_identical(
    normalized$records$source_row_id,
    case$domains[[1L]]$row_id
  )
  expect_identical(
    normalized$records$source_is_anchor,
    seq_len(nrow(case$domains[[1L]]$X)) %in% case$anchors$source_row
  )
  expect_identical(
    normalized$records$source_is_private,
    case$truth$source_private_mask
  )
  expect_true(nzchar(normalized$prediction_hash))
})

test_that("adapter contracts freeze method search spaces and equal tier budgets", {
  case <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 20L,
    family = "sbm",
    n_anchors = 4L,
    tier = "small_accuracy",
    tuning_budget = 12L,
    seed = 3203L,
    case_id = "adapter-budget-contract"
  )
  methods <- c(
    "rapid_ma", "token_ot_graph", "fpgw", "gw", "spectral_mnn",
    "ssma", "kema", "global_geo", "lra", "gpca", "cone_align",
    "grasp", "parrot", "lse"
  )
  contracts <- lapply(methods, function(method) {
    manifoldalign:::.graph_benchmark_adapter_contract_v2(method, case)
  })

  expect_true(all(vapply(
    contracts, function(x) x$parameter_search_budget == 12L, logical(1)
  )))
  expect_true(all(vapply(
    contracts, function(x) x$wall_clock_budget_sec == 120, logical(1)
  )))
  expect_true(all(nzchar(vapply(
    contracts, `[[`, character(1), "search_space_id"
  ))))
  expect_true(all(nzchar(vapply(
    contracts, `[[`, character(1), "search_space_hash"
  ))))
})

test_that("extended package comparators use truthful cross-domain contracts", {
  res <- benchmark_graph_alignment_methods(
    sizes = 12L,
    d = 3L,
    noise_sd = 0.01,
    structure = "sbm",
    n_anchors = 4L,
    methods = c("gw", "spectral_mnn", "kema", "global_geo", "lse"),
    n_reps = 1L,
    seed = 41L,
    verbose = FALSE
  )

  expect_setequal(res$method_id, c(
    "gromov_wasserstein", "spectral_mnn_align", "kema",
    "global_geo_align", "linear_sim_embed"
  ))
  aligned <- res$method_id != "linear_sim_embed"
  expect_true(all(res$status[aligned] == "ok"))
  expect_true(all(is.finite(res$top1_accuracy[aligned])))
  expect_true(all(res$row_identity_preserved[aligned]))
  expect_false(any(res$observed_topology_used))
  expect_identical(
    res$topology_disposition[res$method_id == "gromov_wasserstein"],
    "not_applicable"
  )
  expect_true(all(
    res$topology_disposition[res$method_id %in% c(
      "spectral_mnn_align", "kema", "global_geo_align"
    )] == "feature_derived_only"
  ))
  lse <- res[res$method_id == "linear_sim_embed", ]
  expect_identical(lse$status, "unsupported")
  expect_identical(lse$status_code, "unsupported_domain_count")
  expect_true(is.na(lse$top1_accuracy))
})

test_that("v2 case carries observed graphs, mixed node data, and overlap truth", {
  case <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 30L,
    n_private = c(4L, 7L),
    family = "sbm",
    attr_dim = 5L,
    n_anchors = 4L,
    directed = TRUE,
    misleading_fraction = 0.1,
    seed = 101L,
    case_id = "contract-smoke"
  )

  expect_s3_class(case, "graph_alignment_case_v2")
  expect_identical(case$protocol_version, "2.0.0")
  expect_equal(nrow(manifoldalign:::validate_graph_alignment_case_v2(case)), 0L)
  expect_equal(
    vapply(case$domains, function(x) nrow(x$X), integer(1)),
    c(domain1 = 34L, domain2 = 37L)
  )
  for (domain in case$domains) {
    expect_true(all(c(
      "from", "to", "weight", "relation", "directed"
    ) %in% names(domain$relations)))
    expect_gt(nrow(domain$relations), 0L)
    expect_true(all(is.finite(domain$relations$weight)))
    expect_true(all(domain$relations$weight >= 0))
    expect_true(all(domain$relations$directed))
    expect_true(any(vapply(domain$attributes, is.numeric, logical(1))))
    expect_true(any(vapply(domain$attributes, is.character, logical(1))))
    expect_identical(domain$shared_mask, !domain$private_mask)
    expect_identical(length(domain$row_id), nrow(domain$X))
  }
  expect_true("misleading" %in% case$domains[[1L]]$relations$relation)
  expect_equal(sum(case$truth$source_private_mask), 4L)
  expect_equal(sum(case$truth$target_private_mask), 7L)
  expect_true(all(is.na(case$truth$map12[case$truth$source_private_mask])))
})

test_that("v2 correspondence and anchors are exact after independent row permutations", {
  case <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 36L,
    n_private = c(5L, 3L),
    family = "geometric",
    n_anchors = 6L,
    seed = 102L,
    case_id = "truth-permutation"
  )
  source <- case$domains[[1L]]
  target <- case$domains[[2L]]
  shared <- which(source$shared_mask)
  mapped <- case$truth$map12[shared]

  expect_equal(source$row_id[shared], target$row_id[mapped])
  expect_equal(sort(mapped), which(target$shared_mask))
  expect_equal(
    case$anchors$target_row,
    case$truth$map12[case$anchors$source_row]
  )
  expect_true(all(source$node_role[case$anchors$source_row] == "train"))
  expect_equal(nrow(case$anchors), 6L)
})

test_that("v2 hashes are deterministic and detect content mutation", {
  args <- list(
    n_shared = 28L,
    n_private = c(2L, 4L),
    family = "scale_free",
    attr_dim = 4L,
    n_anchors = 3L,
    directed = TRUE,
    seed = 103L,
    case_id = "hash-repeat"
  )
  first <- do.call(manifoldalign:::synthetic_graph_alignment_case_v2, args)
  second <- do.call(manifoldalign:::synthetic_graph_alignment_case_v2, args)

  expect_identical(first$provenance$config_hash, second$provenance$config_hash)
  expect_identical(first$provenance$case_hash, second$provenance$case_hash)
  expect_identical(first$domains, second$domains)

  changed <- first
  changed$domains[[1L]]$relations$weight[[1L]] <- -1
  codes <- manifoldalign:::validate_graph_alignment_case_v2(changed)$code
  expect_true("invalid_relation_values" %in% codes)
  expect_true("case_hash" %in% codes)
})

test_that("edge severity changes observed graph while holding feature clouds fixed", {
  make_case <- function(edge_delete) {
    manifoldalign:::synthetic_graph_alignment_case_v2(
      n_shared = 80L,
      family = "sbm",
      attr_dim = 6L,
      edge_delete = edge_delete,
      edge_add = 0,
      edge_rewire = 0,
      weight_noise = 0,
      heterophily = 0,
      misleading_fraction = 0,
      position_jitter = 0,
      missing_attribute_fraction = 0,
      attribute_noise = 0,
      seed = 104L,
      case_id = paste0("delete-", edge_delete)
    )
  }
  mild <- make_case(0.1)
  severe <- make_case(0.6)

  expect_identical(mild$domains[[1L]]$X, severe$domains[[1L]]$X)
  expect_identical(mild$domains[[2L]]$X, severe$domains[[2L]]$X)
  expect_gt(
    severe$perturbation$observed$edge_delete_fraction,
    mild$perturbation$observed$edge_delete_fraction
  )
  expect_lt(
    severe$perturbation$observed$structural_edge_jaccard,
    mild$perturbation$observed$structural_edge_jaccard
  )
  expect_false(identical(
    mild$domains[[2L]]$relations,
    severe$domains[[2L]]$relations
  ))
})

test_that("weight rescaling is measured independently of topology", {
  make_case <- function(weight_scale) {
    manifoldalign:::synthetic_graph_alignment_case_v2(
      n_shared = 42L,
      family = "geometric",
      edge_delete = 0,
      edge_add = 0,
      edge_rewire = 0,
      weight_noise = 0,
      weight_scale = weight_scale,
      position_jitter = 0,
      missing_attribute_fraction = 0,
      attribute_noise = 0,
      seed = 105L,
      case_id = paste0("weight-", weight_scale)
    )
  }
  baseline <- make_case(1)
  doubled <- make_case(2)
  edge_columns <- c("from", "to", "relation")

  expect_identical(
    baseline$domains[[2L]]$relations[, edge_columns],
    doubled$domains[[2L]]$relations[, edge_columns]
  )
  expect_equal(
    doubled$domains[[2L]]$relations$weight,
    2 * baseline$domains[[2L]]$relations$weight,
    tolerance = 1e-12
  )
})

test_that("frozen manifests have immutable hashes and cover every v2 tier/family", {
  development <- manifoldalign:::graph_alignment_benchmark_manifest_v2(
    "development"
  )
  final <- manifoldalign:::graph_alignment_benchmark_manifest_v2("final")

  expect_equal(nrow(development), 12L)
  expect_equal(nrow(final), 42L)
  expect_equal(
    sum(final$tier %in% c("small_accuracy", "medium_robustness")),
    30L
  )
  overlap <- with(final, n_private_1 + n_private_2 > 0L)
  expect_true(any(overlap))
  expect_true(any(!overlap & final$tier != "oracle"))
  expect_true(all(development$tuning_budget > 0L))
  expect_true(all(final$tuning_budget == 0L))
  expect_length(intersect(development$seed, final$seed), 0L)
  expect_setequal(
    unique(final$tier),
    c("oracle", "small_accuracy", "medium_robustness", "large_scaling")
  )
  expect_setequal(unique(final$family), c("sbm", "geometric", "scale_free"))
  expect_identical(
    attr(development, "manifest_hash"),
    manifoldalign:::.graph_v2_manifest_hashes[["development"]]
  )
  expect_identical(
    attr(final, "manifest_hash"),
    manifoldalign:::.graph_v2_manifest_hashes[["final"]]
  )
})

test_that("final truth cannot enter a tuning score", {
  final <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-sbm-01", "final"
  )
  development <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "dev-oracle-sbm", "development"
  )

  expect_error(
    manifoldalign:::score_graph_alignment_case_v2(
      final, final$truth$map12, purpose = "tune"
    ),
    class = "graph_benchmark_final_label_leakage"
  )
  final_score <- manifoldalign:::score_graph_alignment_case_v2(
    final, final$truth$map12, purpose = "report"
  )
  dev_score <- manifoldalign:::score_graph_alignment_case_v2(
    development, development$truth$map12,
    purpose = "tune", evaluation_role = "dev"
  )
  expect_equal(final_score$top1_accuracy, 1)
  expect_equal(dev_score$top1_accuracy, 1)
  expect_gt(dev_score$n_eligible, 0L)

  missing <- final$truth$map12
  eligible <- which(!is.na(missing))
  missing[eligible[[1L]]] <- NA_integer_
  missing_score <- manifoldalign:::score_graph_alignment_case_v2(
    final, missing, purpose = "report"
  )
  expect_equal(
    missing_score$top1_accuracy,
    (length(eligible) - 1L) / length(eligible)
  )
  expect_equal(
    missing_score$coverage,
    (length(eligible) - 1L) / length(eligible)
  )

  absent <- rep(NA_real_, length(final$truth$map12))
  absent_score <- manifoldalign:::score_graph_alignment_case_v2(
    final, absent, purpose = "report"
  )
  expect_equal(absent_score$top1_accuracy, 0)
  expect_equal(absent_score$coverage, 0)
})

test_that("v2 control baselines are complete, typed, and deterministic", {
  case <- manifoldalign:::synthetic_graph_alignment_case_v2(
    n_shared = 60L,
    n_private = c(4L, 6L),
    family = "sbm",
    attr_dim = 6L,
    n_anchors = 6L,
    label_fraction = 0.25,
    edge_delete = 0.15,
    edge_add = 0.1,
    seed = 106L,
    case_id = "control-baselines"
  )
  registry <- manifoldalign:::graph_benchmark_control_registry_v2()
  first <- manifoldalign:::benchmark_graph_alignment_controls_v2(
    case, evaluation_role = "final"
  )
  second <- manifoldalign:::benchmark_graph_alignment_controls_v2(
    case, evaluation_role = "final"
  )

  expect_setequal(first$control_id, c(
    "raw_attribute_knn", "position_only", "topology_descriptor_only",
    "anchor_only", "row_identity", "source_only_classifier",
    "target_only_classifier", "chance"
  ))
  expect_identical(first$control_id, registry$control_id)
  expect_identical(first$prediction_hash, second$prediction_hash)
  expect_true(all(first$status %in% c("ok", "unsupported")))
  expect_true(all(first$coverage >= 0 & first$coverage <= 1))
  expect_identical(
    first$observed_topology_used,
    first$control_id == "topology_descriptor_only"
  )
  expect_true(all(
    first$top1_accuracy[is.finite(first$top1_accuracy)] >= 0 &
      first$top1_accuracy[is.finite(first$top1_accuracy)] <= 1
  ))
  target_only <- first[first$control_id == "target_only_classifier", ]
  expect_identical(target_only$task, "target_classification")
  expect_identical(target_only$status, "ok")
  source_only <- first[first$control_id == "source_only_classifier", ]
  expect_identical(source_only$task, "target_classification")
  expect_identical(source_only$status, "ok")
  expect_true(is.finite(source_only$top1_accuracy))
})

test_that("control runner inherits the final-label tuning boundary", {
  case <- manifoldalign:::graph_alignment_case_from_manifest_v2(
    "final-oracle-geometric-01", "final"
  )
  expect_error(
    manifoldalign:::benchmark_graph_alignment_controls_v2(
      case, purpose = "tune", evaluation_role = "dev"
    ),
    class = "graph_benchmark_final_label_leakage"
  )
})
