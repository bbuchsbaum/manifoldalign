.graph_benchmark_wall_clock_budgets_v2 <- c(
  oracle = 30,
  small_accuracy = 120,
  medium_robustness = 600,
  large_scaling = 1800
)

.graph_benchmark_search_spaces_v2 <- list(
  rapid_ma = list(
    ncomp = c(8L, 16L), candidate_k = c(32L, 64L, 128L),
    max_iter = c(2L, 4L)
  ),
  token_ot_graph_align = list(
    ncomp = c(8L, 16L), candidate_k = c(32L, 64L),
    token_mode = c("view_only", "view_plus_neighbors"), n_levels = c(1L, 2L)
  ),
  fpgw = list(
    omega1 = c(0.1, 0.5, 0.9), epsilon = c(0, 0.01),
    max_iter = c(50L, 100L)
  ),
  gromov_wasserstein = list(
    epsilon = c(0.01, 0.05, 0.1), max_iter = c(50L, 100L),
    scale = c("per_domain", "global")
  ),
  spectral_mnn_align = list(
    ncomp = c(8L, 16L), knn = c(5L, 10L, 20L),
    q_descriptors = c(12L, 24L)
  ),
  kema = list(
    ncomp = c(8L, 16L), knn = c(5L, 10L),
    u = c(0.25, 0.5, 0.75)
  ),
  global_geo_align = list(
    ncomp = c(8L, 16L), knn = c(5L, 10L, 20L),
    k_embed = c(12L, 24L)
  ),
  linear_sim_embed = list(
    status = "not_applicable_to_cross_domain_alignment"
  ),
  ssma_align = list(
    ncomp = c(8L, 16L), knn = c(5L, 10L, 20L),
    rank_per_domain = c(32L, 64L)
  ),
  lowrank_align = list(
    ncomp = c(8L, 16L), mu = c(0.25, 0.5, 0.75),
    lambda = c(0.001, 0.01)
  ),
  gpca_align = list(
    ncomp = c(8L, 16L), u = c(0.25, 0.5, 0.75),
    lambda = c(0.01, 0.1)
  ),
  cone_align = list(
    ncomp = c(8L, 16L), sigma = c(0.5, 1), lambda = c(0.05, 0.1)
  ),
  grasp = list(
    ncomp = c(8L, 16L), q_descriptors = c(20L, 30L, 40L)
  ),
  parrot = list(
    ncomp = c(8L, 16L), sigma = c(0.1, 0.2),
    lambda = c(0.05, 0.1, 0.2)
  )
)

.graph_benchmark_relation_profile <- function(case) {
  relations <- lapply(case$domains, `[[`, "relations")
  weights <- unlist(lapply(relations, `[[`, "weight"), use.names = FALSE)
  directed <- unlist(lapply(relations, `[[`, "directed"), use.names = FALSE)
  relation_names <- unlist(
    lapply(relations, `[[`, "relation"), use.names = FALSE
  )
  list(
    edge_count = sum(vapply(relations, nrow, integer(1))),
    signed = any(weights < 0),
    weighted = any(abs(weights - 1) > 1e-12),
    directed = any(directed),
    multirelation = length(unique(relation_names)) > 1L
  )
}

.graph_benchmark_node_profile <- function(case) {
  attributes <- unlist(
    lapply(case$domains, function(domain) as.list(domain$attributes)),
    recursive = FALSE
  )
  list(
    numeric_attributes = any(vapply(attributes, is.numeric, logical(1))),
    categorical_attributes = any(vapply(
      attributes,
      function(x) is.character(x) || is.factor(x),
      logical(1)
    )),
    missing_attributes = any(vapply(attributes, anyNA, logical(1))),
    spatial_positions = any(vapply(
      case$domains,
      function(domain) !is.null(domain$positions) && ncol(domain$positions) > 0L,
      logical(1)
    ))
  )
}

.graph_benchmark_reject <- function(code, detail) {
  list(supported = FALSE, status_code = code, status_detail = detail)
}

.graph_benchmark_validate_topology_contract <- function(info, relation_profile) {
  if (!identical(info$input$primary_class, "supplied_graph")) {
    return(list(supported = TRUE, status_code = "ok", status_detail = NA_character_))
  }
  if (!identical(info$topology$supplied, "native")) {
    return(.graph_benchmark_reject(
      "unsupported_supplied_topology",
      "The method does not accept caller-supplied graph topology."
    ))
  }
  if (isTRUE(relation_profile$signed) &&
      !identical(info$topology$signed_weights, "native")) {
    return(.graph_benchmark_reject(
      "unsupported_signed_weights",
      "The supplied graph has signed weights, which this adapter does not support."
    ))
  }
  if (isTRUE(relation_profile$weighted) &&
      !identical(info$topology$nonnegative_weights, "native")) {
    return(.graph_benchmark_reject(
      "unsupported_weighted_topology",
      "The supplied graph has non-unit weights, which this adapter does not support."
    ))
  }
  if (isTRUE(relation_profile$directed) &&
      !identical(info$topology$directed, "native")) {
    return(.graph_benchmark_reject(
      "unsupported_directed_topology",
      "The supplied graph is directed, which this adapter does not support."
    ))
  }
  if (isTRUE(relation_profile$multirelation) &&
      !identical(info$topology$multirelation, "native")) {
    return(.graph_benchmark_reject(
      "unsupported_multirelation_topology",
      "The supplied graph has multiple relation types, which this adapter does not support."
    ))
  }
  list(supported = TRUE, status_code = "ok", status_detail = NA_character_)
}

.graph_benchmark_adapter_preflight <- function(method, case) {
  info <- alignment_method_info(method)
  n_domains <- length(case$domains)
  if (n_domains < info$domains$min || n_domains > info$domains$max) {
    return(.graph_benchmark_reject(
      "unsupported_domain_count",
      paste0(
        info$display_name, " supports ", info$domains$min, "..",
        info$domains$max, " domains, but this case has ", n_domains, "."
      )
    ))
  }
  relation_profile <- .graph_benchmark_relation_profile(case)
  topology <- .graph_benchmark_validate_topology_contract(info, relation_profile)
  if (!isTRUE(topology$supported)) return(topology)

  if (identical(info$input$primary_class, "supplied_graph")) {
    node_profile <- .graph_benchmark_node_profile(case)
    for (capability in names(node_profile)) {
      if (isTRUE(node_profile[[capability]]) &&
          !identical(info$node_data[[capability]], "native")) {
        return(.graph_benchmark_reject(
          paste0("unsupported_", capability),
          paste0(
            "The supplied case contains ", gsub("_", " ", capability),
            ", which this adapter does not support."
          )
        ))
      }
    }
  }

  unequal_sizes <- nrow(case$domains[[1L]]$X) != nrow(case$domains[[2L]]$X)
  partial_overlap <- any(case$truth$source_private_mask) ||
    any(case$truth$target_private_mask)
  if (unequal_sizes && identical(info$supervision$unequal_sizes, "unsupported")) {
    return(.graph_benchmark_reject(
      "unsupported_unequal_sizes",
      "The method does not support unequal domain sizes."
    ))
  }
  if (partial_overlap &&
      identical(info$supervision$partial_overlap, "unsupported")) {
    return(.graph_benchmark_reject(
      "unsupported_partial_overlap",
      "The method does not support domain-private nodes or partial overlap."
    ))
  }

  requires_anchor_labels <- info$id %in% c(
    "parrot", "lowrank_align", "gpca_align", "kema", "global_geo_align"
  )
  if (requires_anchor_labels && !nrow(case$anchors)) {
    return(.graph_benchmark_reject(
      "requires_anchors",
      paste0(info$display_name, " requires observed anchor supervision in this benchmark.")
    ))
  }
  list(supported = TRUE, status_code = "ok", status_detail = NA_character_)
}

.graph_benchmark_adapter_payload <- function(method, case) {
  info <- alignment_method_info(method)
  features <- lapply(case$domains, `[[`, "X")
  anchors <- case$anchors[, c("source_row", "target_row"), drop = FALSE]
  model_input <- switch(
    info$id,
    rapid_ma = list(
      features = features,
      labels = lapply(case$domains, `[[`, "labels"),
      relations = lapply(case$domains, `[[`, "relations"),
      positions = lapply(case$domains, `[[`, "positions"),
      attributes = lapply(case$domains, `[[`, "attributes"),
      anchors = anchors
    ),
    token_ot_graph_align = list(features = features, anchors = anchors),
    spectral_mnn_align = list(features = features, correspondences = anchors),
    ssma_align = list(features = features, correspondences = anchors),
    kema = list(
      features = features,
      anchor_labels = list(case$anchors1, case$anchors2)
    ),
    global_geo_align = list(features = features, correspondences = anchors),
    lowrank_align = list(
      features = features,
      anchor_labels = list(case$anchors1, case$anchors2)
    ),
    gpca_align = list(
      features = features,
      anchor_labels = list(case$anchors1, case$anchors2)
    ),
    cone_align = list(features = features, anchors = anchors),
    parrot = list(features = features, anchors = anchors),
    grasp = list(features = features),
    fpgw = list(features = features),
    gromov_wasserstein = list(features = features),
    linear_sim_embed = list(features = features),
    stop("No graph benchmark payload contract for method `", info$id, "`.", call. = FALSE)
  )
  identity_contract <- list(
    domains = lapply(case$domains, function(domain) {
      list(
        row_id = domain$row_id,
        node_role = domain$node_role,
        shared_mask = domain$shared_mask,
        private_mask = domain$private_mask
      )
    }),
    anchors = anchors
  )
  list(model_input = model_input, identity_contract = identity_contract)
}

.graph_benchmark_adapter_contract_v2 <- function(method, case) {
  info <- alignment_method_info(method)
  profile <- .graph_benchmark_relation_profile(case)
  preflight <- .graph_benchmark_adapter_preflight(info$id, case)
  payload <- .graph_benchmark_adapter_payload(info$id, case)
  supplied <- identical(info$input$primary_class, "supplied_graph")
  topology_disposition <- if (supplied) {
    "consumed"
  } else if (identical(info$input$primary_class, "feature_derived_graph")) {
    "feature_derived_only"
  } else {
    "not_applicable"
  }
  search_space <- .graph_benchmark_search_spaces_v2[[info$id]]
  if (is.null(search_space)) {
    stop("Missing frozen search space for method `", info$id, "`.", call. = FALSE)
  }
  list(
    adapter_id = paste0("graph-v2-", info$id, "-adapter-v1"),
    method_id = info$id,
    input_class = info$input$primary_class,
    topology_disposition = topology_disposition,
    observed_topology_used = supplied && isTRUE(preflight$supported),
    observed_edge_count = if (supplied && isTRUE(preflight$supported)) {
      profile$edge_count
    } else {
      0L
    },
    adapter_input_hash = .graph_v2_hash(payload$model_input),
    row_identity_hash = .graph_v2_hash(payload$identity_contract),
    search_space_id = paste0("graph-v2-", info$id, "-search-v1"),
    search_space_hash = .graph_v2_hash(search_space),
    parameter_search_budget = as.integer(case$tuning_budget),
    wall_clock_budget_sec = unname(
      .graph_benchmark_wall_clock_budgets_v2[[case$tier]]
    ),
    preflight = preflight
  )
}

.graph_benchmark_normalize_prediction <- function(case, prediction) {
  source <- case$domains[[1L]]
  target <- case$domains[[2L]]
  n_source <- nrow(source$X)
  n_target <- nrow(target$X)
  if (length(prediction) != n_source) {
    return(list(
      ok = FALSE,
      detail = paste0(
        "Malformed node mapping: expected ", n_source,
        " predictions, received ", length(prediction), "."
      )
    ))
  }
  if (!is.numeric(prediction) && !is.integer(prediction)) {
    return(list(ok = FALSE, detail = "Node mapping must be numeric or integer."))
  }
  finite <- is.finite(prediction)
  if (any(finite & prediction != floor(prediction))) {
    return(list(ok = FALSE, detail = "Node mapping contains non-integer target indices."))
  }

  prediction <- as.integer(prediction)
  valid <- !is.na(prediction) & prediction >= 1L & prediction <= n_target
  truth <- as.integer(case$truth$map12)
  eligible <- !is.na(truth)
  anchor <- seq_len(n_source) %in% case$anchors$source_row
  correct <- eligible & valid & prediction == truth
  predicted_row_id <- rep(NA_character_, n_source)
  predicted_row_id[valid] <- target$row_id[prediction[valid]]
  predicted_private <- rep(NA, n_source)
  predicted_private[valid] <- target$private_mask[prediction[valid]]
  truth_row_id <- rep(NA_character_, n_source)
  truth_row_id[eligible] <- target$row_id[truth[eligible]]

  records <- data.frame(
    source_index = seq_len(n_source),
    source_row_id = source$row_id,
    source_role = source$node_role,
    source_is_anchor = anchor,
    source_is_private = source$private_mask,
    truth_target_index = truth,
    truth_target_row_id = truth_row_id,
    predicted_target_index = prediction,
    predicted_target_row_id = predicted_row_id,
    predicted_target_is_private = predicted_private,
    correct = correct,
    stringsAsFactors = FALSE
  )
  subset_accuracy <- function(index) {
    if (!any(index)) return(NA_real_)
    sum(correct[index]) / sum(index)
  }
  private <- source$private_mask
  identity_preserved <- !anyDuplicated(records$source_row_id) &&
    identical(records$source_row_id, source$row_id) &&
    all(records$predicted_target_row_id[valid] == target$row_id[prediction[valid]])

  list(
    ok = TRUE,
    records = records,
    prediction_hash = .graph_v2_hash(records),
    row_identity_preserved = identity_preserved,
    n_eligible = sum(eligible),
    top1_accuracy = subset_accuracy(eligible),
    coverage = if (any(eligible)) sum(valid & eligible) / sum(eligible) else 0,
    anchor_top1_accuracy = subset_accuracy(eligible & anchor),
    nonanchor_top1_accuracy = subset_accuracy(eligible & !anchor),
    private_rejection_accuracy = if (any(private)) {
      mean(!valid[private])
    } else {
      NA_real_
    }
  )
}

#' Benchmark alignment methods on attributed graph protocol v2 cases
#'
#' Generates paired attributed graph domains with observed weighted relations,
#' positions, mixed attributes, row identities, overlap truth, and controlled
#' perturbations, then compares alignment methods available in this package.
#'
#' Interactive calls use development cases from frozen protocol v2 semantics.
#' Final benchmark cases and seeds are held in immutable manifests and are not
#' available for tuning.
#'
#' @param sizes Integer vector of node counts to benchmark.
#' @param d Integer feature dimension for synthetic node features.
#' @param noise_sd Numeric noise standard deviation added to target features.
#' @param structure Character vector of synthetic geometries for base node
#'   coordinates or v2 graph families: `"ring"`, `"grid"`, `"random"`,
#'   `"community"`, `"sbm"`, `"geometric"`, or `"scale_free"`. Legacy names
#'   map to v2 geometric, scale-free, or SBM cases.
#' @param permute_fraction Legacy row-permutation fraction. Protocol v2 freezes
#'   full row permutation, so this must currently be `1`.
#' @param n_anchors Integer number of anchor correspondences to add (0 allowed).
#' @param methods Character vector of methods to run. Options include:
#'   `"rapid_ma"`, `"token_ot_graph"`, `"cone_align"`, `"grasp"`,
#'   `"parrot"`, `"fpgw"`, `"gw"`, `"spectral_mnn"`, `"ssma"`,
#'   `"kema"`, `"global_geo"`, `"lra"`, `"gpca"`, and `"lse"`.
#'   LSE is recorded as unsupported because its registered contract is a
#'   single-domain target-similarity problem, not cross-domain alignment.
#' @param decode_mode Decoding mode for mapping source nodes to targets:
#'   `"native"` uses each method's native assignment/transport decoding;
#'   `"common_nn"` uses nearest-neighbor decoding on aligned scores for all
#'   methods; `"common_procrustes_lap"` uses anchor-Procrustes (when anchors are
#'   available) followed by LAP/Hungarian decoding on aligned scores.
#' @param ssma_procrustes Logical; when `FALSE` (default), SSMA ignores
#'   `"common_procrustes_lap"` and decodes in its own common space via
#'   nearest-neighbor (`"common_nn"`). Set `TRUE` to explicitly benchmark an
#'   `SSMA+Procrustes` decode variant.
#' @param candidate_k Integer number of candidates per source node for
#'   `token_ot_graph_align()`.
#' @param coarsen_method Coarsening method for multilevel Token-OT alignment:
#'   `"kmeans"` or `"louvain"`.
#' @param token_mode Tokenization mode passed to [token_ot_graph_align_control()].
#' @param views Views passed to [token_ot_graph_align_control()].
#' @param n_levels Integer number of multilevel hierarchy levels for
#'   `token_ot_graph_align()`.
#' @param prior_mode Prior lifting mode for multilevel (`"none"`, `"hard"`, `"soft"`).
#' @param fpgw_omega1 Weight of the feature term for `fpgw()` (0–1).
#' @param fpgw_epsilon Entropic warm-start regularization for `fpgw()`.
#' @param fpgw_max_iter Maximum Frank-Wolfe iterations for `fpgw()`.
#' @param fpgw_inner_max_iter Inner iterations for `fpgw()` updates.
#' @param fpgw_tol Convergence tolerance for `fpgw()`.
#' @param ssma_solver Solver passed to [ssma_align_control()] for `ssma_align()`.
#' @param ssma_knn Integer k for within-domain graph construction in SSMA.
#' @param ssma_rank_per_domain Integer per-domain rank cap in SSMA reduction.
#' @param ssma_use_serial Logical; if `TRUE`, enable serial decontamination in SSMA.
#' @param ssma_lag_exclusion Non-negative lag exclusion radius for SSMA serial
#'   hard masking when `ssma_use_serial = TRUE`.
#' @param lra_mu Mixing parameter for `lowrank_align()` (`mu` argument).
#' @param lra_lambda Ridge parameter passed to `lowrank_align()`.
#' @param lra_sv_thresh Singular-value threshold passed to `lowrank_align()`.
#' @param lra_solver Solver backend passed to `lowrank_align()`.
#' @param gpca_u Within/between trade-off parameter passed to `gpca_align()`.
#' @param gpca_lambda Ridge regularization passed to `gpca_align()`.
#' @param n_reps Integer number of replications per size.
#' @param seed Integer seed for reproducibility.
#' @param verbose Logical; print per-run progress.
#'
#' @return A data.frame with per-run runtime and accuracy metrics plus typed
#'   outcome fields. `status` is one of `"ok"`, `"unsupported"`,
#'   `"not_available"`, or `"error"`; `status_code` is machine-readable and
#'   `status_detail` is explanatory text. Accuracy is meaningful only when
#'   `status == "ok"`. `method_id`, `input_class`, and
#'   `topology_disposition` link each row to the method registry and distinguish
#'   supplied-graph, feature-derived, and metric-measure comparisons.
#'   `adapter_input_hash` makes edge-sensitivity auditable. Identity-bearing
#'   prediction hashes, anchor/non-anchor accuracy, private-node rejection,
#'   convergence, frozen search-space IDs, equal tier wall-clock budgets, and
#'   external-memory-measurement status normalize result semantics. Protocol
#'   version, case ID/hash, split, tier, family, preprocessing ID, and observed
#'   edge count preserve benchmark provenance. The legacy
#'   `error` column is retained and populated only for execution/output errors.
#' @export
benchmark_graph_alignment_methods <- function(
  sizes = c(50L, 100L),
  d = 3L,
  noise_sd = 0.05,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 1,
  n_anchors = 0L,
  methods = c("rapid_ma", "token_ot_graph", "fpgw", "cone_align", "grasp", "parrot", "ssma", "lra", "gpca"),
  decode_mode = c("native", "common_nn", "common_procrustes_lap"),
  ssma_procrustes = FALSE,
  candidate_k = 80L,
  coarsen_method = c("kmeans", "louvain"),
  token_mode = c("view_only", "view_plus_neighbors"),
  views = "raw",
  n_levels = 1L,
  prior_mode = c("none", "hard", "soft"),
  fpgw_omega1 = 0.5,
  fpgw_epsilon = 0.01,
  fpgw_max_iter = 50L,
  fpgw_inner_max_iter = 20L,
  fpgw_tol = 1e-6,
  ssma_solver = c("reduced", "operator"),
  ssma_knn = 10L,
  ssma_rank_per_domain = 64L,
  ssma_use_serial = FALSE,
  ssma_lag_exclusion = 1L,
  lra_mu = 0.5,
  lra_lambda = 0.01,
  lra_sv_thresh = 1,
  lra_solver = c("operator", "explicit"),
  gpca_u = 0.5,
  gpca_lambda = 0.1,
  n_reps = 3L,
  seed = 1L,
  verbose = FALSE
) {
  structure <- match.arg(
    structure,
    c("ring", "grid", "random", "community", "sbm", "geometric", "scale_free"),
    several.ok = TRUE
  )
  decode_mode <- match.arg(decode_mode)
  coarsen_method <- match.arg(coarsen_method)
  token_mode <- match.arg(token_mode)
  prior_mode <- match.arg(prior_mode)
  ssma_solver <- match.arg(ssma_solver)
  lra_solver <- match.arg(lra_solver)

  chk::chk_number(d)
  chk::chk_true(d >= 1)
  chk::chk_number(noise_sd)
  chk::chk_true(noise_sd >= 0)
  chk::chk_number(permute_fraction)
  chk::chk_true(is.finite(permute_fraction) && permute_fraction > 0 && permute_fraction <= 1)
  if (!identical(as.numeric(permute_fraction), 1)) {
    stop(
      "Graph benchmark protocol v2 requires `permute_fraction = 1`; partial ",
      "row permutation is not a frozen v2 regime.",
      call. = FALSE
    )
  }
  chk::chk_number(n_anchors)
  chk::chk_true(n_anchors >= 0)
  chk::chk_character(methods)
  chk::chk_number(candidate_k)
  chk::chk_true(candidate_k >= 1)
  chk::chk_logical(ssma_procrustes)
  chk::chk_number(n_levels)
  chk::chk_true(n_levels >= 1)
  chk::chk_number(fpgw_omega1)
  chk::chk_true(is.finite(fpgw_omega1) && fpgw_omega1 >= 0 && fpgw_omega1 <= 1)
  chk::chk_number(fpgw_epsilon)
  chk::chk_true(is.finite(fpgw_epsilon) && fpgw_epsilon >= 0)
  chk::chk_number(fpgw_max_iter)
  chk::chk_true(fpgw_max_iter >= 1)
  chk::chk_number(fpgw_inner_max_iter)
  chk::chk_true(fpgw_inner_max_iter >= 1)
  chk::chk_number(fpgw_tol)
  chk::chk_true(is.finite(fpgw_tol) && fpgw_tol > 0)
  chk::chk_number(ssma_knn)
  chk::chk_true(ssma_knn >= 1)
  chk::chk_number(ssma_rank_per_domain)
  chk::chk_true(ssma_rank_per_domain >= 1)
  chk::chk_logical(ssma_use_serial)
  chk::chk_number(ssma_lag_exclusion)
  chk::chk_true(is.finite(ssma_lag_exclusion) && ssma_lag_exclusion >= 0)
  chk::chk_number(lra_mu)
  chk::chk_true(is.finite(lra_mu) && lra_mu >= 0 && lra_mu <= 1)
  chk::chk_number(lra_lambda)
  chk::chk_true(is.finite(lra_lambda) && lra_lambda >= 0)
  chk::chk_number(lra_sv_thresh)
  chk::chk_true(is.finite(lra_sv_thresh))
  chk::chk_number(gpca_u)
  chk::chk_true(is.finite(gpca_u) && gpca_u >= 0 && gpca_u <= 1)
  chk::chk_number(gpca_lambda)
  chk::chk_true(is.finite(gpca_lambda) && gpca_lambda >= 0)
  chk::chk_number(n_reps)
  chk::chk_true(n_reps >= 1)
  chk::chk_number(seed)
  chk::chk_logical(verbose)

  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) stop("`sizes` must contain integers >= 3.", call. = FALSE)

  allowed_methods <- c(
    "rapid_ma", "token_ot_graph", "fpgw", "cone_align", "grasp",
    "parrot", "ssma", "lra", "gpca", "gw", "spectral_mnn",
    "kema", "global_geo", "lse"
  )
  bad <- setdiff(methods, allowed_methods)
  if (length(bad)) stop("Unknown method(s): ", paste(bad, collapse = ", "), call. = FALSE)

  # Helper to extract top-1 predictions from a transport plan
  top1_from_plan <- function(W) {
    if (inherits(W, "Matrix")) {
      sm <- Matrix::summary(W)
      if (!nrow(sm)) return(rep(NA_integer_, nrow(W)))
      ord <- order(sm$i, -sm$x)
      sm <- sm[ord, ]
      keep <- !duplicated(sm$i)
      out <- rep(NA_integer_, nrow(W))
      out[sm$i[keep]] <- sm$j[keep]
      return(out)
    }
    apply(as.matrix(W), 1, which.max)
  }

  # Helper to extract top-1 predictions from aligned embedding blocks
  top1_from_scores <- function(E1, E2) {
    E1 <- as.matrix(E1)
    E2 <- as.matrix(E2)
    if (!nrow(E1) || !nrow(E2)) return(integer(0))
    if (requireNamespace("RANN", quietly = TRUE)) {
      return(as.integer(RANN::nn2(data = E2, query = E1, k = 1L)$nn.idx[, 1L]))
    }

    # Fallback without RANN: exact pairwise squared distances.
    s1 <- rowSums(E1 * E1)
    s2 <- rowSums(E2 * E2)
    D <- outer(s1, s2, "+") - 2 * (E1 %*% t(E2))
    max.col(-D, ties.method = "first")
  }

  # Helper to compute an orthogonal map from anchor correspondences.
  solve_procrustes_rotation <- function(E1, E2, idx1, idx2) {
    idx1 <- as.integer(idx1)
    idx2 <- as.integer(idx2)
    ok <- is.finite(idx1) & is.finite(idx2) &
      idx1 >= 1L & idx1 <= nrow(E1) &
      idx2 >= 1L & idx2 <= nrow(E2)
    idx1 <- idx1[ok]
    idx2 <- idx2[ok]
    if (length(idx1) < 2L) return(NULL)

    A <- scale(E1[idx1, , drop = FALSE], center = TRUE, scale = FALSE)
    B <- scale(E2[idx2, , drop = FALSE], center = TRUE, scale = FALSE)
    sv <- tryCatch(svd(crossprod(B, A)), error = function(e) NULL)
    if (is.null(sv) || is.null(sv$u) || is.null(sv$v)) return(NULL)

    U <- sv$u
    V <- sv$v
    Q <- U %*% t(V)
    if (ncol(Q) > 1L && det(Q) < 0) {
      U[, ncol(U)] <- -U[, ncol(U)]
      Q <- U %*% t(V)
    }
    Q
  }

  # Helper for LAP/Hungarian decoding from embedding distances.
  lap_from_scores <- function(E1, E2) {
    E1 <- as.matrix(E1)
    E2 <- as.matrix(E2)
    if (!nrow(E1) || !nrow(E2)) return(integer(0))

    if (!requireNamespace("clue", quietly = TRUE)) {
      return(top1_from_scores(E1, E2))
    }

    s1 <- rowSums(E1 * E1)
    s2 <- rowSums(E2 * E2)
    D <- pmax(outer(s1, s2, "+") - 2 * (E1 %*% t(E2)), 0)
    nr <- nrow(D)
    nc <- ncol(D)

    if (nr == nc) {
      return(as.integer(clue::solve_LSAP(D)))
    }

    pad <- max(D) + 1
    if (!is.finite(pad)) pad <- 1
    nmax <- max(nr, nc)
    Dpad <- matrix(pad, nmax, nmax)
    Dpad[seq_len(nr), seq_len(nc)] <- D
    sol <- as.integer(clue::solve_LSAP(Dpad))
    out <- rep(NA_integer_, nr)
    valid <- sol[seq_len(nr)] <= nc
    out[valid] <- sol[seq_len(nr)][valid]
    out
  }

  # Common-score decoder used for apples-to-apples comparisons.
  decode_from_scores_mode <- function(
    scores, n, mode, anchor_idx, anchor_target_idx
  ) {
    if (is.null(scores)) return(NULL)
    scores <- as.matrix(scores)
    if (nrow(scores) < 2L * n) return(NULL)

    E1 <- scores[seq_len(n), , drop = FALSE]
    E2 <- scores[n + seq_len(n), , drop = FALSE]

    if (identical(mode, "common_nn")) {
      return(top1_from_scores(E1, E2))
    }

    if (identical(mode, "common_procrustes_lap")) {
      E2d <- E2
      if (length(anchor_idx) >= 2L) {
        Q <- solve_procrustes_rotation(
          E1, E2, anchor_idx, anchor_target_idx
        )
        if (!is.null(Q)) E2d <- E2 %*% Q
      }
      return(lap_from_scores(E1, E2d))
    }

    top1_from_scores(E1, E2)
  }

  # Label-similarity helper for semi-supervised embedding methods (LRA/GPCA):
  # shared non-missing labels are connected; unlabeled rows remain isolated.
  label_similarity <- function(labels) {
    labels <- as.character(labels)
    n <- length(labels)
    ok <- !is.na(labels)
    S <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    if (sum(ok) >= 2L) {
      M <- outer(labels[ok], labels[ok], FUN = "==")
      S[which(ok), which(ok)] <- Matrix::Matrix(M * 1, sparse = TRUE)
    }
    Matrix::diag(S) <- 0
    Matrix::forceSymmetric(S, uplo = "U")
  }

  rows <- list()
  ptr <- 0L
  specs <- synthetic_graph_alignment_registry(
    sizes = sizes,
    structures = structure,
    d = as.integer(d),
    noise_sd = noise_sd,
    permute_fraction = permute_fraction,
    n_anchors = as.integer(n_anchors),
    n_reps = as.integer(n_reps),
    seed = as.integer(seed)
  )

  make_result_row <- function(
    method,
    spec,
    runtime_sec,
    top1_accuracy = NA_real_,
    coverage = 0,
    error = NA_character_,
    status = if (is.na(error)) "ok" else "error",
    status_code = if (identical(status, "ok")) "ok" else "execution_error",
    status_detail = if (is.na(error)) NA_character_ else error,
    observed_topology_used = FALSE,
    observed_edge_count = 0L,
    prediction_hash = NA_character_,
    row_identity_preserved = NA,
    n_eligible = NA_integer_,
    anchor_top1_accuracy = NA_real_,
    nonanchor_top1_accuracy = NA_real_,
    private_rejection_accuracy = NA_real_,
    converged = NA
  ) {
    status <- match.arg(
      status,
      c("ok", "unsupported", "not_available", "error")
    )
    if (identical(status, "error") &&
        identical(status_code, "execution_error") &&
        !is.na(error) && grepl("elapsed time limit", error, fixed = TRUE)) {
      status_code <- "wall_clock_budget_exceeded"
    }
    info <- alignment_method_info(method)
    contract <- .graph_benchmark_adapter_contract_v2(info$id, case)
    if (!identical(
      isTRUE(observed_topology_used),
      isTRUE(contract$observed_topology_used)
    ) || !identical(
      as.integer(observed_edge_count),
      as.integer(contract$observed_edge_count)
    )) {
      stop(
        "Result topology metadata disagrees with adapter contract for `",
        info$id, "`.", call. = FALSE
      )
    }
    if (!identical(status, "error")) {
      error <- NA_character_
    }
    data.frame(
      method = method,
      method_id = info$id,
      input_class = info$input$primary_class,
      adapter_id = contract$adapter_id,
      topology_disposition = contract$topology_disposition,
      observed_topology_used = contract$observed_topology_used,
      observed_edge_count = as.integer(contract$observed_edge_count),
      adapter_input_hash = contract$adapter_input_hash,
      row_identity_hash = contract$row_identity_hash,
      row_identity_preserved = as.logical(row_identity_preserved),
      prediction_hash = prediction_hash,
      search_space_id = contract$search_space_id,
      search_space_hash = contract$search_space_hash,
      parameter_search_budget = contract$parameter_search_budget,
      parameter_trials_used = if (is.finite(runtime_sec)) 1L else 0L,
      wall_clock_budget_sec = contract$wall_clock_budget_sec,
      peak_rss_mb = NA_real_,
      memory_measurement = "external_process_required",
      protocol_version = if ("protocol_version" %in% names(spec)) {
        spec$protocol_version[[1L]]
      } else {
        NA_character_
      },
      benchmark_case_id = if ("benchmark_case_id" %in% names(spec)) {
        spec$benchmark_case_id[[1L]]
      } else {
        NA_character_
      },
      benchmark_case_hash = if ("benchmark_case_hash" %in% names(spec)) {
        spec$benchmark_case_hash[[1L]]
      } else {
        NA_character_
      },
      benchmark_split = if ("benchmark_split" %in% names(spec)) {
        spec$benchmark_split[[1L]]
      } else {
        NA_character_
      },
      benchmark_tier = if ("benchmark_tier" %in% names(spec)) {
        spec$benchmark_tier[[1L]]
      } else {
        NA_character_
      },
      benchmark_family = if ("benchmark_family" %in% names(spec)) {
        spec$benchmark_family[[1L]]
      } else {
        NA_character_
      },
      preprocessing_id = if ("preprocessing_id" %in% names(spec)) {
        spec$preprocessing_id[[1L]]
      } else {
        NA_character_
      },
      scenario_family = spec$scenario_family[[1L]],
      scenario = spec$scenario[[1L]],
      structure = spec$structure[[1L]],
      n = spec$n[[1L]],
      d = spec$d[[1L]],
      noise_sd = spec$noise_sd[[1L]],
      permute_fraction = spec$permute_fraction[[1L]],
      n_anchors = spec$n_anchors[[1L]],
      rep = spec$rep[[1L]],
      generation_seed = spec$generation_seed[[1L]],
      decode_mode = decode_mode,
      runtime_sec = runtime_sec,
      top1_accuracy = top1_accuracy,
      coverage = coverage,
      n_eligible = as.integer(n_eligible),
      anchor_top1_accuracy = anchor_top1_accuracy,
      nonanchor_top1_accuracy = nonanchor_top1_accuracy,
      private_rejection_accuracy = private_rejection_accuracy,
      converged = as.logical(converged),
      status = status,
      status_code = status_code,
      status_detail = status_detail,
      error = error,
      stringsAsFactors = FALSE
    )
  }

  make_prediction_row <- function(
    method, spec, runtime_sec, pred, truth,
    observed_topology_used = FALSE, observed_edge_count = 0L,
    converged = NA
  ) {
    if (!identical(as.integer(truth), as.integer(case$truth$map12))) {
      detail <- "Adapter truth vector disagrees with the identity-bearing v2 case."
      return(make_result_row(
        method, spec, runtime_sec,
        status = "error", status_code = "malformed_output",
        status_detail = detail, error = detail,
        observed_topology_used = observed_topology_used,
        observed_edge_count = observed_edge_count
      ))
    }
    normalized <- .graph_benchmark_normalize_prediction(case, pred)
    if (!isTRUE(normalized$ok)) {
      return(make_result_row(
        method, spec, runtime_sec,
        status = "error", status_code = "malformed_output",
        status_detail = normalized$detail, error = normalized$detail,
        observed_topology_used = observed_topology_used,
        observed_edge_count = observed_edge_count
      ))
    }
    make_result_row(
      method, spec, runtime_sec,
      top1_accuracy = normalized$top1_accuracy,
      coverage = normalized$coverage,
      observed_topology_used = observed_topology_used,
      observed_edge_count = observed_edge_count,
      prediction_hash = normalized$prediction_hash,
      row_identity_preserved = normalized$row_identity_preserved,
      n_eligible = normalized$n_eligible,
      anchor_top1_accuracy = normalized$anchor_top1_accuracy,
      nonanchor_top1_accuracy = normalized$nonanchor_top1_accuracy,
      private_rejection_accuracy = normalized$private_rejection_accuracy,
      converged = converged
    )
  }

  adapter_is_ready <- function(request_name, method = request_name) {
    if (!request_name %in% methods) return(FALSE)
    contract <- .graph_benchmark_adapter_contract_v2(method, case)
    if (isTRUE(contract$preflight$supported)) return(TRUE)
    ptr <<- ptr + 1L
    rows[[ptr]] <<- make_result_row(
      method, spec, NA_real_,
      status = "unsupported",
      status_code = contract$preflight$status_code,
      status_detail = contract$preflight$status_detail,
      observed_topology_used = contract$observed_topology_used,
      observed_edge_count = contract$observed_edge_count
    )
    FALSE
  }

  run_with_budget <- function(code) {
    budget <- unname(.graph_benchmark_wall_clock_budgets_v2[[case$tier]])
    if (!is.finite(budget)) return(force(code))
    setTimeLimit(cpu = Inf, elapsed = budget, transient = TRUE)
    on.exit(
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE),
      add = TRUE
    )
    force(code)
  }

  for (row_id in seq_len(nrow(specs))) {
    spec <- specs[row_id, , drop = FALSE]
    n <- spec$n[[1L]]
    r <- spec$rep[[1L]]

    family <- switch(
      spec$structure[[1L]],
      ring = "geometric",
      grid = "geometric",
      random = "scale_free",
      community = "sbm",
      sbm = "sbm",
      geometric = "geometric",
      scale_free = "scale_free"
    )
    tier <- if (n <= 30L) {
      "oracle"
    } else if (n <= 250L) {
      "small_accuracy"
    } else if (n <= 5000L) {
      "medium_robustness"
    } else {
      "large_scaling"
    }
    case <- synthetic_graph_alignment_case_v2(
      n_shared = n,
      n_private = c(0L, 0L),
      family = family,
      attr_dim = spec$d[[1L]],
      n_anchors = spec$n_anchors[[1L]],
      label_fraction = min(
        0.25,
        max(0, spec$n_anchors[[1L]] / spec$n[[1L]])
      ),
      edge_delete = min(0.35, spec$noise_sd[[1L]]),
      edge_add = min(0.2, spec$noise_sd[[1L]] / 2),
      edge_rewire = min(0.2, spec$noise_sd[[1L]] / 2),
      weight_noise = spec$noise_sd[[1L]],
      directed = identical(family, "scale_free"),
      position_jitter = spec$noise_sd[[1L]],
      missing_attribute_fraction = min(0.25, spec$noise_sd[[1L]]),
      attribute_noise = spec$noise_sd[[1L]],
      class_imbalance = if (identical(family, "sbm")) 0.2 else 0,
      private_classes = FALSE,
      split = "development",
      tier = tier,
      tuning_budget = 12L,
      seed = spec$generation_seed[[1L]],
      case_id = paste0("interactive-", spec$scenario[[1L]])
    )

    spec$scenario_family <- case$benchmark_meta$scenario_family
    spec$protocol_version <- case$protocol_version
    spec$benchmark_family <- family
    spec$benchmark_case_id <- case$case_id
    spec$benchmark_case_hash <- case$provenance$case_hash
    spec$benchmark_split <- case$split
    spec$benchmark_tier <- case$tier
    spec$preprocessing_id <- case$preprocessing$id

    X1 <- case$X1
    X2 <- case$X2
    map12 <- case$map12
    anchor_idx <- case$anchor_idx
    anchor_target_idx <- case$anchors$target_row
    hd <- case$hd
    anchors <- case$anchors_name

    # Embedding dimensions for graph-based aligners should depend on graph size,
    # not the raw feature dimension used to generate coordinates.
    ncomp_align <- min(10L, as.integer(n) - 2L)

    if (isTRUE(verbose)) {
      message("benchmark_graph_alignment_methods: structure=", spec$structure[[1L]], " n=", n, " rep=", r)
    }

      # RAPID-MA: the supplied-graph comparator in protocol v2.
      if (adapter_is_ready("rapid_ma")) {
        observed_edge_count <- sum(vapply(
          case$domains,
          function(domain) nrow(domain$relations),
          integer(1)
        ))
        t0 <- proc.time()[[3]]
        res <- tryCatch(run_with_budget({
          fit <- suppressMessages(rapid_ma(
            lapply(case$domains, `[[`, "X"),
            labels = lapply(case$domains, `[[`, "labels"),
            relations = case$relation_lists,
            positions = lapply(case$domains, `[[`, "positions"),
            attributes = lapply(case$domains, `[[`, "attributes"),
            ncomp = ncomp_align,
            control = list(
              max_iter = 2L,
              min_iter = 1L,
              seed = as.integer(spec$generation_seed[[1L]]),
              relation = list(
                build_feature = FALSE,
                custom_symmetrize = if (isTRUE(case$config$directed)) {
                  "preserve"
                } else {
                  "mutual"
                }
              )
            )
          ))
          anchor_pairs <- if (nrow(case$anchors)) {
            as.matrix(case$anchors[, c("source_row", "target_row")])
          } else {
            NULL
          }
          matching <- rapid_ma_match(
            fit,
            anchors = anchor_pairs,
            candidate_cap = min(n, max(8L, as.integer(candidate_k))),
            prototype_bucket_cap = min(n, max(16L, 2L * as.integer(candidate_k))),
            candidate_views = c(
              "position", "attribute", "structure", "latent"
            ),
            assignment = "one_to_one"
          )
          list(fit = fit, matching = matching)
        }), error = function(e) e)
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "rapid_ma", spec, runtime,
            status = "error", status_code = "execution_error",
            status_detail = conditionMessage(res),
            error = conditionMessage(res),
            observed_topology_used = TRUE,
            observed_edge_count = observed_edge_count
          )
        } else {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_prediction_row(
            "rapid_ma", spec, runtime,
            res$matching$matches$target, map12,
            observed_topology_used = TRUE,
            observed_edge_count = observed_edge_count
          )
        }
      }

      # Token-OT Graph Align
      if (adapter_is_ready("token_ot_graph", "token_ot_graph_align")) {
        ctrl <- token_ot_graph_align_control(
          n_levels = n_levels,
          prior_mode = prior_mode,
          views = views,
          candidate_k = candidate_k,
          coarsen_method = coarsen_method,
          token_mode = token_mode,
          verbose = FALSE
        )

        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(token_ot_graph_align(
            hd, anchors = anchors, ncomp = ncomp_align, control = ctrl
          )),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("token_ot_graph_align", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          pred_native <- as.integer(res$assignment)
          if (all(is.na(pred_native)) && !is.null(res$transport_plan)) {
            pred_native <- top1_from_plan(res$transport_plan)
          }
          pred <- pred_native
          if (!identical(decode_mode, "native")) {
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ptr <- ptr + 1L
          rows[[ptr]] <- make_prediction_row(
            "token_ot_graph_align", spec, runtime, pred, map12
          )
        }
      }

      # FPGW (fused partial Gromov-Wasserstein)
      if (adapter_is_ready("fpgw")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(
            fpgw(
              hd,
              omega1 = fpgw_omega1,
              epsilon = fpgw_epsilon,
              max_iter = as.integer(fpgw_max_iter),
              inner_max_iter = as.integer(fpgw_inner_max_iter),
              tol = fpgw_tol,
              verbose = FALSE
            )
          )),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("fpgw", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          plan <- NULL
          if (!is.null(res$transport_plans) && length(res$transport_plans) >= 1L) {
            plan <- res$transport_plans[[1L]]
          }
          if (is.null(plan)) {
            detail <- "Missing transport plan from fpgw() result"
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row(
              "fpgw", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            pred_native <- top1_from_plan(plan)
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ptr <- ptr + 1L
            rows[[ptr]] <- make_prediction_row(
              "fpgw", spec, runtime, pred, map12,
              converged = all(res$converged)
            )
          }
        }
      }

      # Balanced Gromov-Wasserstein metric-measure comparator.
      if (adapter_is_ready("gw", "gromov_wasserstein")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(gromov_wasserstein(
            hd,
            epsilon = 0.05,
            max_iter = as.integer(fpgw_max_iter),
            inner_max_iter = as.integer(fpgw_inner_max_iter),
            tol = fpgw_tol,
            verbose = FALSE
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "gromov_wasserstein", spec, runtime, coverage = 0,
            error = conditionMessage(res)
          )
        } else if (is.null(res$transport_plans) ||
                   !length(res$transport_plans)) {
          detail <- "Missing transport plan from gromov_wasserstein() result"
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "gromov_wasserstein", spec, runtime,
            status = "error", status_code = "malformed_output",
            status_detail = detail, error = detail
          )
        } else {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_prediction_row(
            "gromov_wasserstein", spec, runtime,
            top1_from_plan(res$transport_plans[[1L]]), map12,
            converged = all(res$converged)
          )
        }
      }

      # SSMA (semi-supervised manifold alignment baseline)
      if (adapter_is_ready("ssma", "ssma_align")) {
        ssma_corr <- NULL
        if (length(anchor_idx)) {
          ssma_corr <- case$correspondences
          ssma_corr$weight <- rep(1, nrow(ssma_corr))
          ssma_corr$source <- rep("anchor", nrow(ssma_corr))
        }

        ssma_serial_index <- if (isTRUE(ssma_use_serial)) list(seq_len(n), seq_len(n)) else NULL
        ssma_ctrl <- ssma_align_control(
          knn = as.integer(ssma_knn),
          rank_per_domain = as.integer(ssma_rank_per_domain),
          solver = ssma_solver,
          verbose = FALSE,
          serial_control = ssma_serial_control(
            enabled = isTRUE(ssma_use_serial),
            row_whiten = if (isTRUE(ssma_use_serial)) "ar1" else "none",
            lag_mode = if (isTRUE(ssma_use_serial)) "hard" else "none",
            lag_exclusion = as.integer(ssma_lag_exclusion),
            preserve_degree = TRUE
          )
        )

        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(
            ssma_align(
              hd,
              correspondences = ssma_corr,
              serial_index = ssma_serial_index,
              ncomp = ncomp_align,
              control = ssma_ctrl
            )
          )),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("ssma_align", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          pred_native <- decode_from_scores_mode(res$s, n = n, mode = "common_nn", anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
          if (is.null(pred_native)) {
            detail <- "Missing or malformed score matrix from ssma_align()"
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row(
              "ssma_align", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              ssma_mode <- decode_mode
              if (identical(ssma_mode, "common_procrustes_lap") && !isTRUE(ssma_procrustes)) {
                ssma_mode <- "common_nn"
              }
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = ssma_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ptr <- ptr + 1L
            rows[[ptr]] <- make_prediction_row(
              "ssma_align", spec, runtime, pred, map12
            )
          }
        }
      }

      # Spectral MNN feature-derived graph comparator.
      if (adapter_is_ready("spectral_mnn", "spectral_mnn_align")) {
        spectral_corr <- if (nrow(case$anchors)) case$correspondences else NULL
        spectral_ncomp <- min(4L, ncomp_align)
        spectral_ctrl <- spectral_mnn_align_control(
          knn = min(as.integer(ssma_knn), n - 1L),
          q_descriptors = 12L,
          mnn_k = min(8L, n - 1L),
          max_anchors = max(32L, 4L * n),
          refine_rounds = 1L,
          seed = as.integer(spec$generation_seed[[1L]]),
          verbose = FALSE
        )
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(spectral_mnn_align(
            hd,
            correspondences = spectral_corr,
            ncomp = spectral_ncomp,
            control = spectral_ctrl
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "spectral_mnn_align", spec, runtime, coverage = 0,
            error = conditionMessage(res)
          )
        } else {
          pred <- decode_from_scores_mode(
            res$s, n = n, mode = if (identical(decode_mode, "native")) {
              "common_nn"
            } else {
              decode_mode
            },
            anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx
          )
          ptr <- ptr + 1L
          if (is.null(pred)) {
            detail <- "Missing or malformed score matrix from spectral_mnn_align()"
            rows[[ptr]] <- make_result_row(
              "spectral_mnn_align", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            rows[[ptr]] <- make_prediction_row(
              "spectral_mnn_align", spec, runtime, pred, map12
            )
          }
        }
      }

      # Paper-faithful KEMA, supervised by the frozen anchor labels.
      if (adapter_is_ready("kema")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(kema(
            hd,
            y = anchors,
            ncomp = min(4L, ncomp_align),
            knn = min(5L, n - 1L),
            sample_frac = 1,
            solver = "exact",
            backend = "auto"
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "kema", spec, runtime, coverage = 0,
            error = conditionMessage(res)
          )
        } else {
          pred <- decode_from_scores_mode(
            res$s, n = n, mode = if (identical(decode_mode, "native")) {
              "common_nn"
            } else {
              decode_mode
            },
            anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx
          )
          ptr <- ptr + 1L
          if (is.null(pred)) {
            detail <- "Missing or malformed score matrix from kema()"
            rows[[ptr]] <- make_result_row(
              "kema", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            rows[[ptr]] <- make_prediction_row(
              "kema", spec, runtime, pred, map12,
              converged = isTRUE(res$fidelity$passed)
            )
          }
        }
      }

      # GlobalGeo with explicit frozen anchor correspondences.
      if (adapter_is_ready("global_geo", "global_geo_align")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(global_geo_align(
            hd,
            correspondences = case$correspondences,
            preproc = NULL,
            ncomp = min(4L, ncomp_align),
            control = global_geo_align_control(
              knn = min(8L, n - 1L),
              n_landmarks_total = min(2L * n, max(8L, n)),
              k_embed = min(8L, n - 1L),
              seed = as.integer(spec$generation_seed[[1L]]),
              verbose = FALSE
            )
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "global_geo_align", spec, runtime, coverage = 0,
            error = conditionMessage(res)
          )
        } else {
          pred <- decode_from_scores_mode(
            res$s, n = n, mode = if (identical(decode_mode, "native")) {
              "common_nn"
            } else {
              decode_mode
            },
            anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx
          )
          ptr <- ptr + 1L
          if (is.null(pred)) {
            detail <- "Missing or malformed score matrix from global_geo_align()"
            rows[[ptr]] <- make_result_row(
              "global_geo_align", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            rows[[ptr]] <- make_prediction_row(
              "global_geo_align", spec, runtime, pred, map12
            )
          }
        }
      }

      # LRA (low-rank alignment; semi-supervised labels)
      if (adapter_is_ready("lra", "lowrank_align")) {
        if (as.integer(n_anchors) <= 0L) {
          detail <- "LRA benchmark requires n_anchors > 0"
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "lowrank_align", spec, NA_real_,
            status = "unsupported", status_code = "requires_anchors",
            status_detail = detail
          )
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            run_with_budget(suppressWarnings(
              suppressMessages(
                lowrank_align(
                  hd,
                  y = anchors,
                  ncomp = ncomp_align,
                  simfun = label_similarity,
                  mu = lra_mu,
                  lambda = lra_lambda,
                  sv_thresh = lra_sv_thresh,
                  solver = lra_solver
                )
              )
            )),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("lowrank_align", spec, runtime, coverage = 0, error = as.character(res$message))
          } else {
            scores <- if (!is.null(res$s)) as.matrix(res$s) else NULL
            if (is.null(scores) || nrow(scores) < 2L * n) {
              detail <- "Missing or malformed score matrix from lowrank_align()"
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row(
                "lowrank_align", spec, runtime,
                status = "error", status_code = "malformed_output",
                status_detail = detail, error = detail
              )
            } else {
              pred_native <- decode_from_scores_mode(scores, n = n, mode = "common_nn", anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
              pred <- pred_native
              if (!identical(decode_mode, "native")) {
                pred_common <- decode_from_scores_mode(scores, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
                if (!is.null(pred_common)) pred <- pred_common
              }
              ptr <- ptr + 1L
              rows[[ptr]] <- make_prediction_row(
                "lowrank_align", spec, runtime, pred, map12
              )
            }
          }
        }
      }

      # GPCA alignment (semi-supervised labels)
      if (adapter_is_ready("gpca", "gpca_align")) {
        if (as.integer(n_anchors) <= 0L) {
          detail <- "GPCA benchmark requires n_anchors > 0"
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "gpca_align", spec, NA_real_,
            status = "unsupported", status_code = "requires_anchors",
            status_detail = detail
          )
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            run_with_budget(suppressWarnings(
              suppressMessages(
                gpca_align(
                  hd,
                  y = anchors,
                  ncomp = ncomp_align,
                  simfun = label_similarity,
                  u = gpca_u,
                  lambda = gpca_lambda,
                  control = gpca_align_control(verbose = FALSE)
                )
              )
            )),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("gpca_align", spec, runtime, coverage = 0, error = as.character(res$message))
          } else {
            scores <- if (!is.null(res$s)) as.matrix(res$s) else NULL
            if (is.null(scores) || nrow(scores) < 2L * n) {
              detail <- "Missing or malformed score matrix from gpca_align()"
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row(
                "gpca_align", spec, runtime,
                status = "error", status_code = "malformed_output",
                status_detail = detail, error = detail
              )
            } else {
              pred_native <- decode_from_scores_mode(scores, n = n, mode = "common_nn", anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
              pred <- pred_native
              if (!identical(decode_mode, "native")) {
                pred_common <- decode_from_scores_mode(scores, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
                if (!is.null(pred_common)) pred <- pred_common
              }
              ptr <- ptr + 1L
              rows[[ptr]] <- make_prediction_row(
                "gpca_align", spec, runtime, pred, map12
              )
            }
          }
        }
      }

      # CONE-Align
      if (adapter_is_ready("cone_align")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(cone_align(
            hd, anchors = anchors, ncomp = ncomp_align
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("cone_align", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          pred_native <- as.integer(res$assignment)
          pred <- pred_native
          if (!identical(decode_mode, "native")) {
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ptr <- ptr + 1L
          rows[[ptr]] <- make_prediction_row(
            "cone_align", spec, runtime, pred, map12
          )
        }
      }

      # GRASP
      if (adapter_is_ready("grasp")) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          run_with_budget(suppressMessages(grasp(
            hd, ncomp = ncomp_align, q_descriptors = 30
          ))),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("grasp", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          pred_native <- as.integer(res$assignment)
          pred <- pred_native
          if (!identical(decode_mode, "native")) {
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ptr <- ptr + 1L
          rows[[ptr]] <- make_prediction_row(
            "grasp", spec, runtime, pred, map12
          )
        }
      }

      # PARROT (requires anchors)
      if (adapter_is_ready("parrot")) {
        if (as.integer(n_anchors) <= 0L) {
          detail <- "PARROT requires n_anchors > 0"
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row(
            "parrot", spec, NA_real_,
            status = "unsupported", status_code = "requires_anchors",
            status_detail = detail
          )
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            run_with_budget(suppressMessages(parrot(
              hd, anchors = anchors, ncomp = ncomp_align, max_iter = 30
            ))),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("parrot", spec, runtime, coverage = 0, error = as.character(res$message))
          } else if (is.null(res$transport_plan)) {
            detail <- "Missing transport plan from parrot() result"
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row(
              "parrot", spec, runtime,
              status = "error", status_code = "malformed_output",
              status_detail = detail, error = detail
            )
          } else {
            pred_native <- top1_from_plan(res$transport_plan)
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, anchor_target_idx = anchor_target_idx)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ptr <- ptr + 1L
            rows[[ptr]] <- make_prediction_row(
              "parrot", spec, runtime, pred, map12
            )
          }
        }
      }

      # LSE is a single-domain target-similarity estimator. Calling preflight
      # records an explicit unsupported row instead of pretending it is a
      # two-domain graph aligner.
      invisible(adapter_is_ready("lse", "linear_sim_embed"))
  }

  do.call(rbind, rows[seq_len(ptr)])
}
