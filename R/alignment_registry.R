# Canonical alignment method registry -------------------------------------

.alignment_registry_topology <- function(
  supplied = "unsupported",
  nonnegative_weights = "unsupported",
  signed_weights = "unsupported",
  directed = "unsupported",
  multirelation = "unsupported"
) {
  list(
    supplied = supplied,
    nonnegative_weights = nonnegative_weights,
    signed_weights = signed_weights,
    directed = directed,
    multirelation = multirelation
  )
}

.alignment_registry_node_data <- function(
  numeric_attributes = "feature_matrix",
  categorical_attributes = "unsupported",
  missing_attributes = "unsupported",
  spatial_positions = "unsupported"
) {
  list(
    numeric_attributes = numeric_attributes,
    categorical_attributes = categorical_attributes,
    missing_attributes = missing_attributes,
    spatial_positions = spatial_positions
  )
}

.alignment_registry_supervision <- function(
  labels = "unsupported",
  anchors = "unsupported",
  requirement = "none",
  partial_overlap = "unsupported",
  unequal_sizes = "unsupported",
  private_classes = "not_applicable"
) {
  list(
    labels = labels,
    anchors = anchors,
    requirement = requirement,
    partial_overlap = partial_overlap,
    unequal_sizes = unequal_sizes,
    private_classes = private_classes
  )
}

.alignment_registry_domains <- function(min = 2L, max = 2L, native_multi = FALSE) {
  list(min = as.integer(min), max = max, native_multi = isTRUE(native_multi))
}

.alignment_registry_adapter <- function(
  constructor = NA_character_,
  group = NA_character_,
  supports_multi = FALSE
) {
  list(
    constructor = constructor,
    group = group,
    supports_multi = isTRUE(supports_multi)
  )
}

.alignment_registry_evidence <- function(
  correctness = "smoke",
  accuracy = "none",
  scalability = "none",
  oracle_ids = character(),
  benchmark_ids = character()
) {
  list(
    correctness = correctness,
    accuracy = accuracy,
    scalability = scalability,
    oracle_ids = oracle_ids,
    benchmark_ids = benchmark_ids,
    reviewed_by = "alignment-evidence-audit",
    reviewed_on = "2026-07-22"
  )
}

.alignment_registry_lifecycle <- function(
  id,
  state = "experimental",
  correctness = "smoke",
  accuracy = "none",
  scalability = "none",
  oracle_ids = character(),
  benchmark_ids = character(),
  findings = character()
) {
  list(
    method = id,
    state = state,
    recommended = FALSE,
    quarantined = FALSE,
    api_enabled = TRUE,
    findings = findings,
    evidence = .alignment_registry_evidence(
      correctness = correctness,
      accuracy = accuracy,
      scalability = scalability,
      oracle_ids = oracle_ids,
      benchmark_ids = benchmark_ids
    )
  )
}

.alignment_registry_scale <- function(
  status = "unvalidated",
  n_total = NA_integer_,
  attempted_n_total = NA_integer_,
  evidence_id = NA_character_,
  note = "No provenance-complete scale ceiling has been established."
) {
  list(
    status = status,
    n_total = as.integer(n_total),
    attempted_n_total = as.integer(attempted_n_total),
    unit = "total_nodes",
    evidence_id = evidence_id,
    note = note
  )
}

.alignment_registry_complexity <- function(runtime, storage) {
  list(runtime = runtime, storage = storage)
}

.alignment_registry_entry <- function(
  id,
  display_name,
  entry_points,
  aliases = character(),
  benchmark_aliases = character(),
  deprecated_entry_points = list(),
  adapter = .alignment_registry_adapter(),
  lifecycle,
  primary_input_class,
  accepted_input_classes = primary_input_class,
  topology,
  node_data = .alignment_registry_node_data(),
  supervision = .alignment_registry_supervision(),
  domains = .alignment_registry_domains(),
  inductive_oos = "unsupported",
  outputs,
  complexity,
  optional_dependencies = character(),
  validated_scale_ceiling = .alignment_registry_scale(),
  summary
) {
  list(
    id = id,
    display_name = display_name,
    entry_points = entry_points,
    aliases = aliases,
    benchmark_aliases = benchmark_aliases,
    deprecated_entry_points = deprecated_entry_points,
    adapter = adapter,
    lifecycle = lifecycle,
    input = list(
      primary_class = primary_input_class,
      accepted_classes = accepted_input_classes
    ),
    topology = topology,
    node_data = node_data,
    supervision = supervision,
    domains = domains,
    inductive_oos = inductive_oos,
    outputs = outputs,
    complexity = complexity,
    optional_dependencies = optional_dependencies,
    validated_scale_ceiling = validated_scale_ceiling,
    summary = summary
  )
}

#' Canonical alignment method registry
#'
#' Returns the machine-readable source of truth for public alignment methods.
#' Registry records distinguish supplied graphs from feature-derived topology
#' and metric-measure geometry; describe topology, node-data, supervision,
#' output, complexity, dependency, and scale contracts; and embed lifecycle
#' metadata validated by [validate_alignment_lifecycle()].
#'
#' `alignment_method_info()` resolves a canonical method ID, exported entry
#' point, adapter name, or benchmark alias to exactly one record.
#' `alignment_method_table()` returns a flattened data frame suitable for
#' documentation and audits. `validate_alignment_method_registry()` returns
#' typed schema and contradiction problems.
#'
#' @param method One canonical ID, exported entry point, adapter name, or
#'   benchmark alias.
#' @param registry A registry returned by `alignment_method_registry()`.
#'
#' @return `alignment_method_registry()` returns a named list with class
#'   `alignment_method_registry`. `alignment_method_info()` returns one record.
#'   `alignment_method_table()` returns a data frame.
#'   `validate_alignment_method_registry()` returns an
#'   `alignment_registry_problems` data frame; zero rows means valid.
#'
#' @examples
#' registry <- alignment_method_registry()
#' rapid <- alignment_method_info("rapid_ma", registry)
#' rapid$topology$supplied
#' alignment_method_table(registry)[, c("id", "lifecycle", "primary_input")]
#'
#' @name alignment_method_registry
NULL

#' @rdname alignment_method_registry
#' @export
alignment_method_registry <- function() {
  feature_topology <- .alignment_registry_topology(
    supplied = "unsupported",
    nonnegative_weights = "derived_only",
    signed_weights = "unsupported",
    directed = "unsupported",
    multirelation = "unsupported"
  )
  metric_topology <- .alignment_registry_topology(
    supplied = "not_applicable",
    nonnegative_weights = "not_applicable",
    signed_weights = "not_applicable",
    directed = "not_applicable",
    multirelation = "not_applicable"
  )
  no_topology <- .alignment_registry_topology(
    supplied = "unsupported",
    nonnegative_weights = "unsupported",
    signed_weights = "unsupported",
    directed = "unsupported",
    multirelation = "unsupported"
  )
  feature_nodes <- .alignment_registry_node_data()
  no_supervision <- .alignment_registry_supervision(
    unequal_sizes = "native"
  )

  registry <- list(
    rapid_ma = .alignment_registry_entry(
      id = "rapid_ma",
      display_name = "RAPID-MA",
      entry_points = c(
        "rapid_ma", "rapid_ma_semisup", "rapid_ma_match", "rapid_ma_aligner"
      ),
      aliases = "rapid",
      adapter = .alignment_registry_adapter("rapid_ma_aligner", "perm", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "rapid_ma", correctness = "invariant", accuracy = "comparative",
        scalability = "comparative",
        oracle_ids = c(
          "tests/testthat/test-rapid-ma.R",
          "tests/testthat/test-rapid-ma-uot.R",
          "tests/testthat/test-rapid-ma-relations.R"
        ),
        benchmark_ids = c(
          "docs/BENCHMARK_PROTOCOL.md",
          "inst/benchmarks/graph_alignment_v2_final_manifest.csv",
          "docs/rapid-ma-validation-summary.csv",
          "docs/rapid-ma-scaling-summary.csv",
          "docs/rapid-ma-matching-summary.csv"
        ),
        findings = c(
          "position_control_outperforms_final_matching",
          "blocked_inductive_improvement_gate_failed",
          "million_node_attempt_timed_out"
        )
      ),
      primary_input_class = "supplied_graph",
      accepted_input_classes = c("supplied_graph", "feature_manifold"),
      topology = .alignment_registry_topology(
        supplied = "native", nonnegative_weights = "native",
        signed_weights = "unsupported", directed = "native",
        multirelation = "native"
      ),
      node_data = .alignment_registry_node_data(
        numeric_attributes = "native", categorical_attributes = "native",
        missing_attributes = "native", spatial_positions = "native"
      ),
      supervision = .alignment_registry_supervision(
        labels = "optional", anchors = "optional",
        partial_overlap = "native", unequal_sizes = "native",
        private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "couplings", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(E s t + N q I) for bounded relation width, sketch width s, prototypes q",
        "O(E + N(s + q)) additional sparse storage"
      ),
      validated_scale_ceiling = .alignment_registry_scale(
        status = "validated", n_total = 22600L,
        attempted_n_total = 22600L,
        evidence_id = "docs/evidence/alignment-tournament-v2/raw_results.csv",
        note = paste0(
          "The current full-scoring external-process tournament completed ",
          "22,600 total nodes. The regenerated scaling runner completed ",
          "110,000 total nodes on a 4-core machine (220,000 exceeded its ",
          "four-minute cap there; an older run on a faster machine completed ",
          "it), and the 1,000,002-node attempt timed out."
        )
      ),
      summary = paste0(
        "Sparse structure-aware prototype alignment; the only current method ",
        "with native supplied directed nonnegative multirelations plus positions ",
        "and mixed attributes. Experimental, not a universal accuracy winner."
      )
    ),
    token_ot_graph_align = .alignment_registry_entry(
      id = "token_ot_graph_align",
      display_name = "Token-OT Graph Align",
      entry_points = c("token_ot_graph_align", "token_ot_graph_aligner"),
      aliases = "token_ot_graph",
      benchmark_aliases = "token_ot_graph",
      adapter = .alignment_registry_adapter(
        "token_ot_graph_aligner", "perm", FALSE
      ),
      lifecycle = .alignment_registry_lifecycle(
        "token_ot_graph_align", correctness = "invariant",
        accuracy = "comparative", scalability = "comparative",
        oracle_ids = "tests/testthat/test-token-ot-graph-align-contracts.R",
        benchmark_ids = c(
          "docs/rapid-ma-scaling-summary.csv",
          "dev-scripts/benchmarks/graph_aligners_benchmark_results.csv"
        )
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        anchors = "optional", partial_overlap = "unsupported",
        unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, 2L, FALSE),
      inductive_oos = "unsupported",
      outputs = c("embeddings", "couplings", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(E k + I C) with bounded candidate support C",
        "O(E + C + N k) for sparse candidates and embeddings"
      ),
      optional_dependencies = c("RcppAnnoy", "Rnanoflann"),
      validated_scale_ceiling = .alignment_registry_scale(
        "validated", 5500L, 11000L, "docs/rapid-ma-scaling-summary.csv",
        "Regenerated 2026-10-19 (4-core x86_64, ten-second cap): 5,500 total nodes completed; 11,000 timed out."
      ),
      summary = paste0(
        "Sparse token/candidate OT over graphs derived from the feature matrix; ",
        "it does not consume an observed adjacency contract."
      )
    ),
    spectral_mnn_align = .alignment_registry_entry(
      id = "spectral_mnn_align",
      display_name = "Spectral MNN Align",
      entry_points = c("spectral_mnn_align", "spectral_mnn_aligner"),
      aliases = "spectral_mnn",
      adapter = .alignment_registry_adapter(
        "spectral_mnn_aligner", "O", TRUE
      ),
      lifecycle = .alignment_registry_lifecycle(
        "spectral_mnn_align", correctness = "invariant",
        accuracy = "comparative", scalability = "comparative",
        oracle_ids = "tests/testthat/test-spectral_mnn_align.R",
        benchmark_ids = c(
          "docs/rapid-ma-scaling-summary.csv",
          "dev-scripts/benchmarks/global_geo_vs_baselines_summary.csv"
        )
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "optional", anchors = "optional",
        partial_overlap = "unsupported", unequal_sizes = "native",
        private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(E k + N log N + N k^2)",
        "O(E + N k)"
      ),
      validated_scale_ceiling = .alignment_registry_scale(
        "validated", 11000L, 22000L, "docs/rapid-ma-scaling-summary.csv",
        "Regenerated 2026-10-19 (4-core x86_64, ten-second cap): 11,000 total nodes completed; 22,000 timed out."
      ),
      summary = "Fast HKS/MNN spectral alignment on feature-derived kNN graphs."
    ),
    ssma_align = .alignment_registry_entry(
      id = "ssma_align",
      display_name = "SSMA",
      entry_points = c("ssma_align", "ssma_aligner"),
      aliases = "ssma",
      benchmark_aliases = "ssma",
      adapter = .alignment_registry_adapter("ssma_aligner", "GL", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "ssma_align", correctness = "invariant", accuracy = "comparative",
        scalability = "comparative",
        oracle_ids = "tests/testthat/test-ssma-align.R",
        benchmark_ids = c(
          "docs/rapid-ma-scaling-summary.csv",
          "dev-scripts/benchmarks/toy_feature_seed_audit_results.csv"
        )
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "optional", anchors = "optional",
        partial_overlap = "unsupported", unequal_sizes = "native",
        private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(E k + N r^2) for reduced rank r",
        "O(E + N r)"
      ),
      validated_scale_ceiling = .alignment_registry_scale(
        "validated", 11000L, 22000L, "docs/rapid-ma-scaling-summary.csv",
        "Regenerated 2026-10-19 (4-core x86_64, ten-second cap): 11,000 total nodes completed; 22,000 timed out."
      ),
      summary = "Semi-supervised spectral alignment using feature-derived graphs and correspondences."
    ),
    global_geo_align = .alignment_registry_entry(
      id = "global_geo_align",
      display_name = "GlobalGeo",
      entry_points = c("global_geo_align", "global_geo_aligner"),
      aliases = "global_geo",
      adapter = .alignment_registry_adapter(
        "global_geo_aligner", "GL", TRUE
      ),
      lifecycle = .alignment_registry_lifecycle(
        "global_geo_align", state = "legacy_reference",
        correctness = "oracle", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = "tests/testthat/test-global_geo_align.R",
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv",
        findings = c(
          "final_tournament_completion_8_of_21",
          "disconnected_landmark_graph_errors",
          "partial_overlap_contract_withdrawn"
        )
      ),
      primary_input_class = "feature_derived_graph",
      accepted_input_classes = c("feature_derived_graph", "feature_manifold"),
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "optional", anchors = "optional",
        requirement = "labels_or_anchors", partial_overlap = "unsupported",
        unequal_sizes = "native", private_classes = "unsupported"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(E k + N r^2 + r^3)",
        "O(E + N r + r^2)"
      ),
      summary = paste0(
        "Legacy/reference correspondence-guided alignment for connected ",
        "feature-derived global geometry; partial-overlap claims are withdrawn."
      )
    ),
    kema = .alignment_registry_entry(
      id = "kema",
      display_name = "KEMA / REKEMA",
      entry_points = c("kema", "kema_orig", "kema_aligner"),
      aliases = "rekema",
      deprecated_entry_points = list(
        kema_orig = list(
          replacement = "kema", warning_since = "0.1.0",
          warning_through = "0.2.0", removal_at = "1.0.0"
        )
      ),
      adapter = .alignment_registry_adapter("kema_aligner", "kernel", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "kema", correctness = "oracle", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-kema-numerical-validation.R",
          "tests/testthat/test-kema-scalable-smoke.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_kema_lra_results.csv"
      ),
      primary_input_class = "feature_derived_graph",
      accepted_input_classes = c("feature_derived_graph", "feature_manifold"),
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "required", partial_overlap = "unsupported",
        unequal_sizes = "native", private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(N^3) full; reduced/operator paths depend on landmark rank",
        "O(N^2) full kernel/graph storage; reduced kernels remain stored"
      ),
      validated_scale_ceiling = .alignment_registry_scale(
        "failed", NA_integer_, 550L, "docs/rapid-ma-scaling-summary.csv",
        "Regenerated 2026-10-19: certified eigensolves rejected all three 550-node runs; 1,100 to 5,500 nodes timed out and larger tasks exceed max_dense_bytes."
      ),
      summary = paste0(
        "Semi-supervised kernel manifold alignment with repaired generalized-",
        "eigen fidelity; kema_orig() is a deprecated duplicate entry point."
      )
    ),
    gpca_align = .alignment_registry_entry(
      id = "gpca_align",
      display_name = "GPCA Align",
      entry_points = "gpca_align",
      aliases = "gpca",
      benchmark_aliases = "gpca",
      lifecycle = .alignment_registry_lifecycle(
        "gpca_align", correctness = "invariant", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = "tests/testthat/test-gpca_align.R",
        benchmark_ids = c(
          "dev-scripts/benchmarks/toy_feature_seed_audit_results.csv",
          "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv"
        )
      ),
      primary_input_class = "feature_manifold",
      topology = no_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "required", unequal_sizes = "native",
        private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(N r^2 + r^3) after sparse similarity assembly",
        "O(nnz(S) + N r + p r)"
      ),
      summary = "Supervised generalized-PCA alignment for numeric feature manifolds."
    ),
    lowrank_align = .alignment_registry_entry(
      id = "lowrank_align",
      display_name = "Low-rank Alignment",
      entry_points = "lowrank_align",
      aliases = "lra",
      benchmark_aliases = "lra",
      lifecycle = .alignment_registry_lifecycle(
        "lowrank_align", correctness = "invariant", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = "tests/testthat/test-lowrank_align.R",
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_kema_lra_results.csv"
      ),
      primary_input_class = "feature_manifold",
      topology = no_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "required", unequal_sizes = "native",
        private_classes = "experimental"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(N r^2 + r^3) operator or explicit reduced solve",
        "O(nnz(S) + N r + p r)"
      ),
      summary = "Semi-supervised low-rank linear alignment driven by a supplied label similarity."
    ),
    linear_sim_embed = .alignment_registry_entry(
      id = "linear_sim_embed",
      display_name = "Linear Similarity Embedding",
      entry_points = "linear_sim_embed",
      aliases = "lse",
      lifecycle = .alignment_registry_lifecycle(
        "linear_sim_embed", correctness = "oracle", accuracy = "single_fixture",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test_linear_sim_embed_consistency.R",
          "tests/testthat/test-lse-recovery.R",
          "tests/testthat/test-lse-orthogonality.R"
        ),
        benchmark_ids = "tests/testthat/test-toy-manifold-benchmarks-fast.R"
      ),
      primary_input_class = "feature_manifold",
      topology = no_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "unsupported", anchors = "unsupported",
        partial_overlap = "not_applicable", unequal_sizes = "not_applicable"
      ),
      domains = .alignment_registry_domains(1L, 1L, FALSE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(I N^2 k) for a dense target-similarity objective",
        "O(N^2 + p k)"
      ),
      summary = "Single-matrix linear embedding that fits a supplied or derived target similarity."
    ),
    generalized_procrustes = .alignment_registry_entry(
      id = "generalized_procrustes",
      display_name = "Generalized Procrustes",
      entry_points = c("generalized_procrustes", "gpa_aligner"),
      aliases = "gpa",
      adapter = .alignment_registry_adapter("gpa_aligner", "O", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "generalized_procrustes", correctness = "oracle",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = "tests/testthat/test-generalized_procrustes.R",
        benchmark_ids = "dev-scripts/benchmarks/toy_feature_aligners_results.csv"
      ),
      primary_input_class = "feature_manifold",
      topology = no_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "required", partial_overlap = "native",
        unequal_sizes = "native", private_classes = "unsupported"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(I sum_d p k_d^2) plus SVD initialization",
        "O(sum_d p k_d)"
      ),
      summary = "Orthogonal task alignment with native partial task observation."
    ),
    grasp = .alignment_registry_entry(
      id = "grasp",
      display_name = "GRASP",
      entry_points = c("grasp", "grasp_aligner"),
      benchmark_aliases = "grasp",
      adapter = .alignment_registry_adapter("grasp_aligner", "O", FALSE),
      lifecycle = .alignment_registry_lifecycle(
        "grasp", correctness = "invariant", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-grasp.R",
          "tests/testthat/test-grasp-regressions.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/graph_aligners_benchmark_results.csv"
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = no_supervision,
      domains = .alignment_registry_domains(2L, 2L, FALSE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(E k + N k^2 + n_max^3) including dense assignment",
        "O(E + N k + n_1 n_2)"
      ),
      summary = "Functional-map graph matching on kNN graphs derived from node features."
    ),
    grasp_multiset = .alignment_registry_entry(
      id = "grasp_multiset",
      display_name = "GRASP Multiset",
      entry_points = "grasp_multiset",
      lifecycle = .alignment_registry_lifecycle(
        "grasp_multiset", correctness = "invariant", accuracy = "single_fixture",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-grasp_multiset.R",
          "tests/testthat/test-grasp-regressions.R"
        )
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = no_supervision,
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(sum_d E_d k + I D k^3 + assignment)",
        "O(sum_d(E_d + n_d k) + pairwise assignment costs)"
      ),
      summary = "Multi-domain extension of feature-derived GRASP functional maps."
    ),
    cone_align = .alignment_registry_entry(
      id = "cone_align",
      display_name = "CONE-Align",
      entry_points = "cone_align",
      benchmark_aliases = "cone_align",
      lifecycle = .alignment_registry_lifecycle(
        "cone_align", state = "legacy_reference",
        correctness = "invariant", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = "tests/testthat/test-cone_align.R",
        benchmark_ids = "dev-scripts/benchmarks/graph_aligners_benchmark_results.csv",
        findings = c(
          "final_tournament_completion_10_of_12",
          "latent_dimension_procrustes_errors"
        )
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        anchors = "optional", partial_overlap = "unsupported",
        unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, 2L, FALSE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(E k + I(n_1 n_2 k + assignment))",
        "O(E + N k + n_1 n_2)"
      ),
      summary = paste0(
        "Legacy/reference alternating Procrustes/assignment alignment; restore ",
        "to experimental after latent dimensions are reconciled on frozen cases."
      )
    ),
    cone_align_multiple = .alignment_registry_entry(
      id = "cone_align_multiple",
      display_name = "CONE-Align Multiple",
      entry_points = c("cone_align_multiple", "cone_multi_aligner"),
      aliases = "cone_multi",
      adapter = .alignment_registry_adapter("cone_multi_aligner", "O", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "cone_align_multiple", state = "legacy_reference",
        correctness = "invariant",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = "tests/testthat/test-cone_align_multiple.R",
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv",
        findings = "inherits_pairwise_latent_dimension_repair_gate"
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = no_supervision,
      domains = .alignment_registry_domains(3L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(sum_d E_d k + I sum_d n_ref n_d k)",
        "O(sum_d(E_d + n_d k + n_ref n_d))"
      ),
      summary = "Reference-based multi-domain CONE on feature-derived graphs."
    ),
    parrot = .alignment_registry_entry(
      id = "parrot",
      display_name = "PARROT",
      entry_points = c("parrot", "parrot_aligner"),
      benchmark_aliases = "parrot",
      adapter = .alignment_registry_adapter("parrot_aligner", "perm", FALSE),
      lifecycle = .alignment_registry_lifecycle(
        "parrot", correctness = "invariant", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-parrot-computational-contracts.R",
          "tests/testthat/test-parrot-equivalence.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/graph_aligners_benchmark_results.csv"
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        anchors = "required", partial_overlap = "unsupported",
        unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, 2L, FALSE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "couplings", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(I n_1 n_2 + random-walk feature construction)",
        "O(n_1 n_2 + n_1^2 + n_2^2)"
      ),
      optional_dependencies = "Rcpp",
      summary = paste0(
        "Anchor-guided dense OT with random-walk position features built from ",
        "feature-derived kNN networks; no supplied adjacency input."
      )
    ),
    gromov_wasserstein = .alignment_registry_entry(
      id = "gromov_wasserstein",
      display_name = "Gromov-Wasserstein",
      entry_points = "gromov_wasserstein",
      aliases = "gw",
      lifecycle = .alignment_registry_lifecycle(
        "gromov_wasserstein", correctness = "invariant",
        accuracy = "comparative", scalability = "smoke",
        oracle_ids = "tests/testthat/test-gromov_wasserstein.R"
      ),
      primary_input_class = "metric_measure",
      topology = metric_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        partial_overlap = "unsupported", unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = "couplings",
      complexity = .alignment_registry_complexity(
        "O(I(n_1^2 n_2 + n_1 n_2^2))",
        "O(n_1^2 + n_2^2 + n_1 n_2)"
      ),
      summary = "Balanced metric-measure OT built from within-domain distances, not observed graph edges."
    ),
    fpgw = .alignment_registry_entry(
      id = "fpgw",
      display_name = "Fused Partial Gromov-Wasserstein",
      entry_points = "fpgw",
      benchmark_aliases = "fpgw",
      lifecycle = .alignment_registry_lifecycle(
        "fpgw", correctness = "oracle", accuracy = "comparative",
        scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-fpgw-oracle.R",
          "tests/testthat/test-fpgw-validation.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/graph_aligners_benchmark_results.csv",
        findings = "penalized_lambda_withdrawn_exact_mass_oracle_ready"
      ),
      primary_input_class = "metric_measure",
      topology = metric_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        partial_overlap = "native", unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = "couplings",
      complexity = .alignment_registry_complexity(
        "O(I(n_1^2 n_2 + n_1 n_2^2) + I LP(n_1 n_2))",
        "O(n_1^2 + n_2^2 + n_1 n_2)"
      ),
      optional_dependencies = "lpSolve",
      summary = paste0(
        "Feature-fused metric-measure OT with oracle-verified classical and ",
        "fixed-rho partial modes. It does not consume observed topology; the ",
        "incorrect lambda-penalized route is withdrawn."
      )
    ),
    coupled_diagonalization = .alignment_registry_entry(
      id = "coupled_diagonalization",
      display_name = "Coupled Diagonalization",
      entry_points = "coupled_diagonalization",
      lifecycle = .alignment_registry_lifecycle(
        "coupled_diagonalization", correctness = "oracle",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-coupled-diagonalization-oracle.R",
          "tests/testthat/test-coupled_diagonalization.R",
          "tests/testthat/test-coupled_diag_mathematical.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/coupled_diagonalization_refactor_results.csv"
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        anchors = "optional", partial_overlap = "experimental",
        unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(sum_d E_d k + I sum_d n_d k^2)",
        "O(sum_d(E_d + n_d k + p_d k))"
      ),
      summary = "Joint spectral-basis optimization with optional cross-domain correspondence matrices."
    ),
    mma_align_multiple = .alignment_registry_entry(
      id = "mma_align_multiple",
      display_name = "Multiset Manifold Alignment",
      entry_points = "mma_align_multiple",
      aliases = "mma",
      lifecycle = .alignment_registry_lifecycle(
        "mma_align_multiple", correctness = "invariant",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-mma-align-multiple.R",
          "tests/testthat/test-mma-eigengate.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/toy_feature_aligners_results.csv"
      ),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = no_supervision,
      domains = .alignment_registry_domains(3L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms", "node_correspondences"),
      complexity = .alignment_registry_complexity(
        "O(sum_d E_d k + I sum_d n_d K^2 + assignment)",
        "O(sum_d(E_d + n_d k + n_d K))"
      ),
      summary = "Unsupervised multiset spectral embedding with probabilistic registration."
    ),
    multiscale_manifold_align = .alignment_registry_entry(
      id = "multiscale_manifold_align",
      display_name = "Multiscale Manifold Align",
      entry_points = "multiscale_manifold_align",
      aliases = "multiscale",
      lifecycle = .multiscale_lifecycle_record(),
      primary_input_class = "feature_derived_graph",
      topology = feature_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        labels = "optional", anchors = "optional",
        partial_overlap = "experimental", unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "experimental",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "Declared sparse spectral hierarchy; current solver fails mandatory scale cases",
        "Declared O(E + N r); not validated"
      ),
      validated_scale_ceiling = .alignment_registry_scale(
        "failed", NA_integer_, 550L, "docs/rapid-ma-scaling-raw.csv",
        "The first 550-total-node scaling cases fail with a malformed solver formula."
      ),
      summary = "Deprecated incoherent experimental hierarchy; retained only for reproducibility until 1.0.0."
    ),
    multiset_uot_align = .alignment_registry_entry(
      id = "multiset_uot_align",
      display_name = "Multiset UOT",
      entry_points = c("multiset_uot_align", "uot_fit_pair"),
      aliases = c("uot", "uot_alignment"),
      lifecycle = .alignment_registry_lifecycle(
        "multiset_uot_align", correctness = "oracle",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = c(
          "tests/testthat/test-uot-ti-sinkhorn.R",
          "tests/testthat/test-uot-synthetic.R",
          "tests/testthat/test-uot-stress.R"
        ),
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_baselines_results.csv"
      ),
      primary_input_class = "metric_measure",
      accepted_input_classes = c("metric_measure", "feature_manifold"),
      topology = metric_topology,
      node_data = .alignment_registry_node_data(
        numeric_attributes = "native", spatial_positions = "native"
      ),
      supervision = .alignment_registry_supervision(
        anchors = "optional", partial_overlap = "native",
        unequal_sizes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "unsupported",
      outputs = c("couplings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(I nnz(C)) sparse or O(I n m) dense",
        "O(nnz(C)) sparse or O(n m) dense"
      ),
      summary = "Sparse/dense unbalanced transport to pairwise or shared template supports."
    ),
    sumcor_align = .alignment_registry_entry(
      id = "sumcor_align",
      display_name = "Event-keyed SUMCOR",
      entry_points = c("sumcor_align", "sumcor_aligner", "calibrate_domain"),
      aliases = "sumcor",
      adapter = .alignment_registry_adapter("sumcor_aligner", "encoder", TRUE),
      lifecycle = .alignment_registry_lifecycle(
        "sumcor_align", correctness = "oracle", accuracy = "single_fixture",
        scalability = "none",
        oracle_ids = "tests/testthat/test-sumcor-align.R"
      ),
      primary_input_class = "feature_manifold",
      topology = no_topology,
      node_data = feature_nodes,
      supervision = .alignment_registry_supervision(
        anchors = "required", requirement = "event_identifiers",
        partial_overlap = "native", unequal_sizes = "native",
        private_classes = "native"
      ),
      domains = .alignment_registry_domains(2L, Inf, TRUE),
      inductive_oos = "native",
      outputs = c("embeddings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(sum_s n_s p_s^2 + I m^2 E k p) for m domains, E events and basis width p",
        "O(sum_s n_s p_s + m E p)"
      ),
      summary = paste0(
        "Linear per-domain encoders maximising cross-domain agreement of ",
        "event-level means keyed by event identifiers, with exact ",
        "out-of-sample encoding and frozen-reference calibration of new domains."
      )
    ),
    ot_procrustes = .alignment_registry_entry(
      id = "ot_procrustes",
      display_name = "OT-Procrustes",
      entry_points = "ot_procrustes_aligner",
      adapter = .alignment_registry_adapter(
        "ot_procrustes_aligner", "O", FALSE
      ),
      lifecycle = .alignment_registry_lifecycle(
        "ot_procrustes", correctness = "invariant",
        accuracy = "single_fixture", scalability = "smoke",
        oracle_ids = "tests/testthat/test-ot-procrustes-aligner.R",
        benchmark_ids = "dev-scripts/benchmarks/global_geo_vs_kema_lra_results.csv"
      ),
      primary_input_class = "feature_manifold",
      accepted_input_classes = c("feature_manifold", "metric_measure"),
      topology = metric_topology,
      node_data = feature_nodes,
      supervision = no_supervision,
      domains = .alignment_registry_domains(2L, 2L, FALSE),
      inductive_oos = "native",
      outputs = c("couplings", "transforms"),
      complexity = .alignment_registry_complexity(
        "O(I(n_1 n_2 d + d^3))",
        "O(n_1 n_2 + d^2)"
      ),
      summary = "Alternating feature-space OT and orthogonal Procrustes alignment."
    )
  )

  tournament_methods <- c(
    "rapid_ma", "token_ot_graph_align", "fpgw", "gromov_wasserstein",
    "spectral_mnn_align", "ssma_align", "kema", "global_geo_align",
    "lowrank_align", "gpca_align", "cone_align", "grasp", "parrot",
    "linear_sim_embed"
  )
  tournament_evidence <- c(
    "docs/evidence/alignment-tournament-v2/raw_results.csv",
    "docs/evidence/alignment-tournament-v2/INTERPRETATION.md"
  )
  for (method in tournament_methods) {
    registry[[method]]$lifecycle$evidence$benchmark_ids <- unique(c(
      registry[[method]]$lifecycle$evidence$benchmark_ids,
      tournament_evidence
    ))
  }

  class(registry) <- c("alignment_method_registry", "list")
  registry
}

.alignment_registry_lookup_names <- function(entry) {
  unique(c(
    entry$id, entry$entry_points, entry$aliases, entry$benchmark_aliases
  ))
}

#' @rdname alignment_method_registry
#' @export
alignment_method_info <- function(
  method,
  registry = alignment_method_registry()
) {
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("`method` must be one non-empty string.", call. = FALSE)
  }
  matched <- vapply(
    registry,
    function(entry) method %in% .alignment_registry_lookup_names(entry),
    logical(1)
  )
  if (!any(matched)) {
    stop("Unknown alignment method or alias: `", method, "`.", call. = FALSE)
  }
  if (sum(matched) != 1L) {
    stop("Ambiguous alignment method alias: `", method, "`.", call. = FALSE)
  }
  unclass(registry)[[which(matched)]]
}

.alignment_registry_compact <- function(x) {
  if (!length(x)) "none" else paste(x, collapse = ", ")
}

.alignment_registry_scale_label <- function(scale) {
  if (identical(scale$status, "validated")) {
    paste0(format(scale$n_total, big.mark = ","), " ", scale$unit)
  } else if (identical(scale$status, "stale") && is.finite(scale$n_total)) {
    paste0("stale: ", format(scale$n_total, big.mark = ","), " ", scale$unit)
  } else if (identical(scale$status, "failed")) {
    paste0("failed by ", format(scale$attempted_n_total, big.mark = ","), " ", scale$unit)
  } else {
    "unvalidated"
  }
}

#' @rdname alignment_method_registry
#' @export
alignment_method_table <- function(registry = alignment_method_registry()) {
  rows <- lapply(registry, function(entry) {
    max_domains <- if (is.infinite(entry$domains$max)) "Inf" else entry$domains$max
    data.frame(
      id = entry$id,
      method = entry$display_name,
      lifecycle = entry$lifecycle$state,
      correctness = entry$lifecycle$evidence$correctness,
      accuracy = entry$lifecycle$evidence$accuracy,
      scalability = entry$lifecycle$evidence$scalability,
      primary_input = entry$input$primary_class,
      supplied_topology = entry$topology$supplied,
      nonnegative_weights = entry$topology$nonnegative_weights,
      signed_weights = entry$topology$signed_weights,
      directed = entry$topology$directed,
      multirelation = entry$topology$multirelation,
      numeric_attributes = entry$node_data$numeric_attributes,
      categorical_attributes = entry$node_data$categorical_attributes,
      missing_attributes = entry$node_data$missing_attributes,
      spatial_positions = entry$node_data$spatial_positions,
      labels = entry$supervision$labels,
      anchors = entry$supervision$anchors,
      partial_overlap = entry$supervision$partial_overlap,
      unequal_sizes = entry$supervision$unequal_sizes,
      private_classes = entry$supervision$private_classes,
      domains = paste0(entry$domains$min, "..", max_domains),
      inductive_oos = entry$inductive_oos,
      outputs = .alignment_registry_compact(entry$outputs),
      runtime = entry$complexity$runtime,
      storage = entry$complexity$storage,
      optional_dependencies = .alignment_registry_compact(
        entry$optional_dependencies
      ),
      scale_ceiling = .alignment_registry_scale_label(
        entry$validated_scale_ceiling
      ),
      scale_status = entry$validated_scale_ceiling$status,
      scale_evidence = entry$validated_scale_ceiling$evidence_id,
      summary = entry$summary,
      stringsAsFactors = FALSE
    )
  })
  rownames_out <- names(registry)
  out <- do.call(rbind, rows)
  rownames(out) <- rownames_out
  out
}

.alignment_registry_empty_problems <- function() {
  out <- data.frame(
    method = character(), code = character(), field = character(),
    message = character(), stringsAsFactors = FALSE
  )
  class(out) <- c("alignment_registry_problems", "data.frame")
  out
}

.alignment_registry_add_problem <- function(problems, method, code, field, message) {
  rbind(
    problems,
    data.frame(
      method = method, code = code, field = field, message = message,
      stringsAsFactors = FALSE
    )
  )
}

#' @rdname alignment_method_registry
#' @export
validate_alignment_method_registry <- function(
  registry = alignment_method_registry()
) {
  problems <- .alignment_registry_empty_problems()
  add <- function(method, code, field, message) {
    problems <<- .alignment_registry_add_problem(
      problems, method, code, field, message
    )
  }
  required <- c(
    "id", "display_name", "entry_points", "aliases", "benchmark_aliases",
    "deprecated_entry_points", "adapter", "lifecycle", "input", "topology",
    "node_data", "supervision", "domains", "inductive_oos", "outputs",
    "complexity", "optional_dependencies", "validated_scale_ceiling",
    "summary"
  )
  input_classes <- c(
    "feature_manifold", "supplied_graph", "feature_derived_graph",
    "metric_measure"
  )
  topology_levels <- c(
    "native", "derived_only", "experimental", "unsupported", "not_applicable"
  )
  node_levels <- c(
    "native", "feature_matrix", "experimental", "unsupported",
    "not_applicable"
  )
  supervision_levels <- c("required", "optional", "unsupported", "not_applicable")
  support_levels <- c("native", "experimental", "unsupported", "not_applicable")
  output_levels <- c("embeddings", "transforms", "couplings", "node_correspondences")
  oos_levels <- c("native", "experimental", "unsupported")
  scale_levels <- c("validated", "stale", "unvalidated", "failed")

  if (!is.list(registry) || !length(registry)) {
    add("<registry>", "invalid_registry", "registry", "Registry must be a non-empty named list.")
    class(problems) <- c("alignment_registry_problems", "data.frame")
    return(problems)
  }

  all_entry_points <- character()
  all_lookup_names <- character()
  declared_adapters <- character()
  for (index in seq_along(registry)) {
    entry <- registry[[index]]
    method <- if (is.list(entry) && .alignment_scalar_character(entry$id)) {
      entry$id
    } else {
      paste0("<entry-", index, ">")
    }
    if (!is.list(entry)) {
      add(method, "invalid_entry", "entry", "Each registry entry must be a named list.")
      next
    }
    missing <- setdiff(required, names(entry))
    for (field in missing) {
      add(method, "missing_field", field, paste0("Registry entry is missing `", field, "`."))
    }
    if (length(missing)) next
    if (!identical(names(registry)[[index]], entry$id)) {
      add(method, "name_id_mismatch", "id", "Registry list name must equal entry$id.")
    }
    if (!.alignment_nonempty_character(entry$entry_points)) {
      add(method, "invalid_entry_points", "entry_points", "At least one exported entry point is required.")
    }
    all_entry_points <- c(all_entry_points, entry$entry_points)
    all_lookup_names <- c(all_lookup_names, .alignment_registry_lookup_names(entry))

    lifecycle_problems <- validate_alignment_lifecycle(entry$lifecycle)
    if (nrow(lifecycle_problems)) {
      for (row in seq_len(nrow(lifecycle_problems))) {
        add(
          method, paste0("lifecycle_", lifecycle_problems$code[[row]]),
          paste0("lifecycle$", lifecycle_problems$field[[row]]),
          lifecycle_problems$message[[row]]
        )
      }
    }
    if (!identical(entry$lifecycle$method, entry$id)) {
      add(method, "lifecycle_id_mismatch", "lifecycle$method", "Lifecycle method must equal registry ID.")
    }

    if (!entry$input$primary_class %in% input_classes ||
        any(!entry$input$accepted_classes %in% input_classes) ||
        !entry$input$primary_class %in% entry$input$accepted_classes) {
      add(method, "invalid_input_class", "input", "Input classes are invalid or omit the primary class.")
    }
    if (!identical(sort(names(entry$topology)), sort(c(
      "supplied", "nonnegative_weights", "signed_weights", "directed", "multirelation"
    ))) || any(!unlist(entry$topology, use.names = FALSE) %in% topology_levels)) {
      add(method, "invalid_topology", "topology", "Topology fields are incomplete or use an invalid support level.")
    }
    if (!identical(sort(names(entry$node_data)), sort(c(
      "numeric_attributes", "categorical_attributes", "missing_attributes", "spatial_positions"
    ))) || any(!unlist(entry$node_data, use.names = FALSE) %in% node_levels)) {
      add(method, "invalid_node_data", "node_data", "Node-data fields are incomplete or invalid.")
    }
    if (!entry$supervision$labels %in% supervision_levels ||
        !entry$supervision$anchors %in% supervision_levels ||
        any(!unlist(entry$supervision[c(
          "partial_overlap", "unequal_sizes", "private_classes"
        )], use.names = FALSE) %in% support_levels)) {
      add(method, "invalid_supervision", "supervision", "Supervision/sample support fields are invalid.")
    }
    if (!entry$inductive_oos %in% oos_levels) {
      add(method, "invalid_oos", "inductive_oos", "Invalid inductive/OOS support level.")
    }
    if (!.alignment_nonempty_character(entry$outputs) ||
        any(!entry$outputs %in% output_levels)) {
      add(method, "invalid_outputs", "outputs", "Outputs must use recognized semantics.")
    }
    if (!.alignment_scalar_character(entry$complexity$runtime) ||
        !.alignment_scalar_character(entry$complexity$storage)) {
      add(method, "invalid_complexity", "complexity", "Runtime and storage classes must be declared.")
    }
    scale <- entry$validated_scale_ceiling
    if (!is.list(scale) || !scale$status %in% scale_levels ||
        !.alignment_scalar_character(scale$unit) ||
        !.alignment_scalar_character(scale$note)) {
      add(method, "invalid_scale_ceiling", "validated_scale_ceiling", "Scale ceiling metadata is incomplete.")
    } else if (identical(scale$status, "validated") &&
               (!is.finite(scale$n_total) || scale$n_total < 1L ||
                !.alignment_scalar_character(scale$evidence_id))) {
      add(method, "unproven_scale_ceiling", "validated_scale_ceiling", "Validated ceilings require a positive size and evidence ID.")
    }

    if (identical(entry$input$primary_class, "supplied_graph") &&
        !identical(entry$topology$supplied, "native")) {
      add(method, "input_topology_contradiction", "topology$supplied", "A supplied-graph primary input must consume supplied topology natively.")
    }
    if (identical(entry$input$primary_class, "feature_derived_graph") &&
        identical(entry$topology$supplied, "native")) {
      add(method, "input_topology_contradiction", "topology$supplied", "A feature-derived graph method cannot also claim native supplied topology.")
    }
    if (identical(entry$input$primary_class, "metric_measure") &&
        !identical(entry$topology$supplied, "not_applicable")) {
      add(method, "metric_graph_contradiction", "topology$supplied", "Metric-measure methods must not claim observed topology without an explicit graph API.")
    }
    if (identical(entry$topology$signed_weights, "native") &&
        !identical(entry$topology$nonnegative_weights, "native")) {
      add(method, "weight_support_contradiction", "topology$signed_weights", "Native signed weights imply native nonnegative weights.")
    }

    adapter <- entry$adapter
    if (.alignment_scalar_character(adapter$constructor)) {
      declared_adapters <- c(declared_adapters, adapter$constructor)
      if (!adapter$group %in% c("O", "GL", "perm", "kernel", "encoder")) {
        add(method, "invalid_adapter_group", "adapter$group", "Adapter group is invalid.")
      }
      if (!identical(adapter$supports_multi, entry$domains$native_multi)) {
        add(method, "adapter_multi_contradiction", "adapter$supports_multi", "Adapter and registry native-multi claims disagree.")
      }
    }

    if (length(entry$deprecated_entry_points)) {
      for (alias in names(entry$deprecated_entry_points)) {
        metadata <- entry$deprecated_entry_points[[alias]]
        fields <- c("replacement", "warning_since", "warning_through", "removal_at")
        if (!all(fields %in% names(metadata)) ||
            any(!vapply(metadata[fields], .alignment_scalar_character, logical(1)))) {
          add(method, "invalid_deprecated_entry_point", paste0("deprecated_entry_points$", alias), "Deprecated entry points require replacement and warning/removal versions.")
        }
      }
    }
  }

  duplicated_entry_points <- unique(all_entry_points[duplicated(all_entry_points)])
  for (entry_point in duplicated_entry_points) {
    add(entry_point, "duplicate_entry_point", "entry_points", "An exported entry point maps to more than one method.")
  }
  duplicated_lookup <- unique(all_lookup_names[duplicated(all_lookup_names)])
  for (lookup in duplicated_lookup) {
    add(lookup, "duplicate_lookup_name", "aliases", "A method lookup name is ambiguous.")
  }

  exports <- getNamespaceExports("manifoldalign")
  missing_exports <- setdiff(unique(all_entry_points), exports)
  for (entry_point in missing_exports) {
    add(entry_point, "stale_entry_point", "entry_points", "Registry entry point is not exported by the package.")
  }
  exported_adapters <- setdiff(
    grep("_aligner$", exports, value = TRUE),
    "new_aligner"
  )
  for (constructor in setdiff(exported_adapters, declared_adapters)) {
    add(constructor, "missing_adapter", "adapter$constructor", "Exported aligner constructor has no registry entry.")
  }
  for (constructor in setdiff(declared_adapters, exported_adapters)) {
    add(constructor, "stale_adapter", "adapter$constructor", "Declared adapter constructor is not exported.")
  }
  if (anyDuplicated(declared_adapters)) {
    for (constructor in unique(declared_adapters[duplicated(declared_adapters)])) {
      add(constructor, "duplicate_adapter", "adapter$constructor", "Adapter constructor maps to multiple registry methods.")
    }
  }

  class(problems) <- c("alignment_registry_problems", "data.frame")
  rownames(problems) <- NULL
  problems
}

.alignment_registry_markdown <- function(registry = alignment_method_registry()) {
  table <- alignment_method_table(registry)
  clean <- function(x) gsub("\\|", "/", x)
  lines <- c(
    "| Method | Lifecycle/evidence | Input contract | Structure and node data | Supervision/shape | Output, cost, scale |",
    "|---|---|---|---|---|---|"
  )
  for (row in seq_len(nrow(table))) {
    evidence <- paste0(
      table$lifecycle[[row]], "; C/A/S=", table$correctness[[row]], "/",
      table$accuracy[[row]], "/", table$scalability[[row]]
    )
    structure_data <- paste0(
      "topology=", table$supplied_topology[[row]], "; weights=",
      table$nonnegative_weights[[row]], "; directed=", table$directed[[row]],
      "; multirelation=", table$multirelation[[row]], "; attributes=",
      table$numeric_attributes[[row]], "/", table$categorical_attributes[[row]],
      "; positions=", table$spatial_positions[[row]]
    )
    supervision <- paste0(
      "labels=", table$labels[[row]], "; anchors=", table$anchors[[row]],
      "; overlap=", table$partial_overlap[[row]], "; unequal=",
      table$unequal_sizes[[row]], "; domains=", table$domains[[row]],
      "; OOS=", table$inductive_oos[[row]]
    )
    output <- paste0(
      table$outputs[[row]], "; ", table$runtime[[row]], "; storage ",
      table$storage[[row]], "; scale ", table$scale_ceiling[[row]]
    )
    lines <- c(
      lines,
      paste0(
        "| ", clean(table$method[[row]]), " (`", table$id[[row]], "`) | ",
        clean(evidence), " | ", clean(table$primary_input[[row]]), " | ",
        clean(structure_data), " | ", clean(supervision), " | ",
        clean(output), " |"
      )
    )
  }
  lines
}
