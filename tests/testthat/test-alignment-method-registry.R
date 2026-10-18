library(testthat)
library(manifoldalign)

matrix_method_ids <- function(path = test_path("ALIGNMENT_TEST_MATRIX.yml")) {
  lines <- readLines(path, warn = FALSE)
  start <- match("methods:", lines)
  end <- match("utilities:", lines)
  stopifnot(!is.na(start), !is.na(end), end > start)
  block <- lines[seq.int(start + 1L, end - 1L)]
  keys <- grepl("^  [a-z][a-z0-9_]*:$", block)
  sub(":$", "", trimws(block[keys]))
}

test_that("registry is schema-valid and matches the enforced method matrix", {
  registry <- alignment_method_registry()
  problems <- validate_alignment_method_registry(registry)

  expect_s3_class(registry, "alignment_method_registry")
  expect_s3_class(problems, "alignment_registry_problems")
  expect_equal(nrow(problems), 0L, info = paste(problems$message, collapse = "\n"))
  expect_setequal(names(registry), matrix_method_ids())
  expect_equal(length(registry), length(unique(names(registry))))
})

test_that("every exported concrete aligner constructor maps exactly once", {
  registry <- alignment_method_registry()
  exports <- getNamespaceExports("manifoldalign")
  exported <- setdiff(grep("_aligner$", exports, value = TRUE), "new_aligner")
  declared <- vapply(registry, function(entry) {
    entry$adapter$constructor
  }, character(1))
  declared <- declared[!is.na(declared)]

  expect_setequal(declared, exported)
  expect_identical(anyDuplicated(declared), 0L)

  for (constructor in declared) {
    entry <- alignment_method_info(constructor, registry)
    algo <- getExportedValue("manifoldalign", constructor)()
    caps <- aligner_capabilities(algo)
    expect_identical(caps$group, entry$adapter$group, info = constructor)
    expect_identical(
      caps$supports_multi,
      entry$adapter$supports_multi,
      info = constructor
    )
  }
})

test_that("public entry points are exported and map to one canonical method", {
  registry <- alignment_method_registry()
  exports <- getNamespaceExports("manifoldalign")
  entry_points <- unlist(lapply(registry, `[[`, "entry_points"), use.names = FALSE)

  expect_true(all(entry_points %in% exports))
  expect_identical(anyDuplicated(entry_points), 0L)
  for (entry_point in entry_points) {
    expect_identical(
      alignment_method_info(entry_point, registry)$id,
      unname(names(registry)[vapply(
        registry,
        function(entry) entry_point %in% entry$entry_points,
        logical(1)
      )]),
      info = entry_point
    )
  }
})

test_that("graph input contracts distinguish observed, derived, and metric geometry", {
  rapid <- alignment_method_info("rapid_ma")
  expect_identical(rapid$input$primary_class, "supplied_graph")
  expect_identical(rapid$topology$supplied, "native")
  expect_identical(rapid$topology$nonnegative_weights, "native")
  expect_identical(rapid$topology$directed, "native")
  expect_identical(rapid$topology$multirelation, "native")
  expect_identical(rapid$topology$signed_weights, "unsupported")
  expect_identical(rapid$node_data$categorical_attributes, "native")
  expect_identical(rapid$node_data$missing_attributes, "native")
  expect_identical(rapid$node_data$spatial_positions, "native")

  derived <- c(
    "grasp", "cone_align", "parrot", "token_ot_graph_align",
    "spectral_mnn_align", "ssma_align"
  )
  expect_true(all(vapply(derived, function(method) {
    entry <- alignment_method_info(method)
    identical(entry$input$primary_class, "feature_derived_graph") &&
      identical(entry$topology$supplied, "unsupported")
  }, logical(1))))

  for (method in c("gromov_wasserstein", "fpgw")) {
    entry <- alignment_method_info(method)
    expect_identical(entry$input$primary_class, "metric_measure")
    expect_identical(entry$topology$supplied, "not_applicable")
  }
})

test_that("lifecycle gate accepts the registry but recommends no unproven method", {
  registry <- alignment_method_registry()
  lifecycles <- lapply(registry, `[[`, "lifecycle")

  expect_true(assert_alignment_release_ready(lifecycles))
  expect_false(any(vapply(lifecycles, `[[`, logical(1), "recommended")))
  expect_identical(
    alignment_method_info("multiscale")$lifecycle$state,
    "deprecated"
  )
  expect_true(
    "kema_orig" %in% names(alignment_method_info("kema")$deprecated_entry_points)
  )
  expect_identical(
    alignment_method_info("global_geo_align")$lifecycle$state,
    "legacy_reference"
  )
  expect_identical(
    alignment_method_info("global_geo_align")$supervision$partial_overlap,
    "unsupported"
  )
  expect_true(all(vapply(c("cone_align", "cone_align_multiple"), function(id) {
    identical(alignment_method_info(id)$lifecycle$state, "legacy_reference")
  }, logical(1))))
})

test_that("scale ceilings expose current, stale, failed, and unvalidated evidence", {
  rapid <- alignment_method_info("rapid_ma")$validated_scale_ceiling
  multiscale <- alignment_method_info(
    "multiscale_manifold_align"
  )$validated_scale_ceiling

  expect_identical(rapid$status, "validated")
  expect_equal(rapid$n_total, 22600L)
  expect_equal(rapid$attempted_n_total, 22600L)
  expect_identical(multiscale$status, "failed")
  expect_equal(multiscale$attempted_n_total, 550L)
  expect_true(any(alignment_method_table()$scale_status == "unvalidated"))
})

test_that("tournament methods cite the frozen raw receipt and limitations", {
  tournament_methods <- c(
    "rapid_ma", "token_ot_graph_align", "fpgw", "gromov_wasserstein",
    "spectral_mnn_align", "ssma_align", "kema", "global_geo_align",
    "lowrank_align", "gpca_align", "cone_align", "grasp", "parrot",
    "linear_sim_embed"
  )
  evidence_id <- "docs/evidence/alignment-tournament-v2/raw_results.csv"
  expect_true(all(vapply(tournament_methods, function(id) {
    evidence_id %in% alignment_method_info(id)$lifecycle$evidence$benchmark_ids
  }, logical(1))))
  expect_setequal(
    alignment_method_info("rapid_ma")$lifecycle$findings,
    c(
      "position_control_outperforms_final_matching",
      "blocked_inductive_improvement_gate_failed",
      "million_node_attempt_timed_out"
    )
  )
})

test_that("registry validator catches contract contradictions", {
  registry <- alignment_method_registry()
  registry$rapid_ma$topology$supplied <- "unsupported"
  registry$fpgw$topology$supplied <- "native"
  registry$ssma_align$adapter$supports_multi <- FALSE

  codes <- validate_alignment_method_registry(registry)$code
  expect_true("input_topology_contradiction" %in% codes)
  expect_true("metric_graph_contradiction" %in% codes)
  expect_true("adapter_multi_contradiction" %in% codes)
})

test_that("flattened method table preserves required comparison fields", {
  table <- alignment_method_table()
  expect_equal(nrow(table), length(alignment_method_registry()))
  expect_true(all(c(
    "lifecycle", "correctness", "accuracy", "scalability",
    "primary_input", "supplied_topology", "nonnegative_weights",
    "signed_weights", "directed", "multirelation", "numeric_attributes",
    "categorical_attributes", "missing_attributes", "spatial_positions",
    "labels", "anchors", "partial_overlap", "unequal_sizes",
    "private_classes", "domains", "inductive_oos", "outputs", "runtime",
    "storage", "optional_dependencies", "scale_ceiling", "scale_evidence"
  ) %in% names(table)))
})

test_that("generated comparison table cannot drift from registry", {
  root <- test_path("..", "..")
  paths <- c(
    file.path(root, "METHOD_COMPARISON_TABLE.md"),
    file.path(root, "docs", "METHOD_COMPARISON_TABLE.md")
  )
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
  expect_identical(copies[[1L]], copies[[2L]])
})
