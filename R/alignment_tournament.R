.alignment_tournament_version_v2 <- "2.0.0-tournament.1"

.alignment_tournament_methods_v2 <- function() {
  c(
    "rapid_ma", "token_ot_graph_align", "fpgw",
    "gromov_wasserstein", "spectral_mnn_align", "ssma_align",
    "kema", "global_geo_align", "lowrank_align", "gpca_align",
    "cone_align", "grasp", "parrot", "linear_sim_embed"
  )
}

.alignment_tournament_resource_ceiling_v2 <- function() {
  data.frame(
    method_id = .alignment_tournament_methods_v2(),
    max_total_nodes = c(
      25000L, 2000L, 300L, 300L, 2000L, 2000L, 300L,
      300L, 300L, 300L, 300L, 300L, 300L, 0L
    ),
    envelope_id = c(
      "rapid-sparse-25k-v1", "token-sparse-2k-v1",
      "fpgw-dense-300-v1", "gw-dense-300-v1",
      "spectral-sparse-2k-v1", "ssma-reduced-2k-v1",
      "kema-dense-300-v1", "globalgeo-small-300-v1",
      "lra-small-300-v1", "gpca-small-300-v1",
      "cone-dense-assignment-300-v1", "grasp-dense-assignment-300-v1",
      "parrot-dense-transport-300-v1", "lse-not-cross-domain-v1"
    ),
    stringsAsFactors = FALSE
  )
}

.alignment_tournament_parameters_v2 <- function(method, case) {
  method <- alignment_method_info(method)$id
  n1 <- nrow(case$domains[[1L]]$X)
  n2 <- nrow(case$domains[[2L]]$X)
  ncomp <- max(2L, min(8L, n1 - 2L, n2 - 2L))
  candidate_k <- min(64L, n2)
  switch(
    method,
    rapid_ma = list(
      ncomp = ncomp, candidate_k = candidate_k,
      max_iter = if (case$tier %in% c("small_accuracy", "medium_robustness")) 4L else 2L,
      candidate_views = c("position", "attribute", "structure", "latent"),
      assignment = "one_to_one"
    ),
    token_ot_graph_align = list(
      ncomp = ncomp, candidate_k = candidate_k, n_levels = 1L,
      prior_mode = "none", views = "raw", token_mode = "view_only",
      coarsen_method = "kmeans"
    ),
    fpgw = list(
      omega1 = 0.5, epsilon = 0.01, max_iter = 50L,
      inner_max_iter = 20L, tol = 1e-6
    ),
    gromov_wasserstein = list(
      epsilon = 0.05, max_iter = 50L, inner_max_iter = 20L,
      tol = 1e-6
    ),
    spectral_mnn_align = list(
      ncomp = min(4L, ncomp), knn = min(10L, n1 - 1L, n2 - 1L),
      q_descriptors = 12L, mnn_k = min(8L, n1 - 1L, n2 - 1L),
      refine_rounds = 1L
    ),
    ssma_align = list(
      ncomp = ncomp, knn = min(10L, n1 - 1L, n2 - 1L),
      rank_per_domain = min(64L, n1 - 1L, n2 - 1L), solver = "reduced"
    ),
    kema = list(
      ncomp = min(4L, ncomp), knn = min(5L, n1 - 1L, n2 - 1L),
      sample_frac = 1, solver = "exact", backend = "auto"
    ),
    global_geo_align = list(
      ncomp = min(4L, ncomp), knn = min(8L, n1 - 1L, n2 - 1L),
      n_landmarks_total = min(n1 + n2, max(8L, n1)),
      k_embed = min(8L, n1 - 1L, n2 - 1L)
    ),
    lowrank_align = list(
      ncomp = ncomp, mu = 0.5, lambda = 0.01,
      sv_thresh = 1, solver = "operator"
    ),
    gpca_align = list(
      ncomp = ncomp, u = 0.5, lambda = 0.1
    ),
    cone_align = list(ncomp = ncomp),
    grasp = list(ncomp = ncomp, q_descriptors = min(30L, n1 - 1L, n2 - 1L)),
    parrot = list(ncomp = ncomp, max_iter = 30L),
    linear_sim_embed = list(status = "not_applicable_to_cross_domain_alignment")
  )
}

.alignment_tournament_applicability_v2 <- function(method, case) {
  info <- alignment_method_info(method)
  contract <- .graph_benchmark_adapter_contract_v2(info$id, case)
  if (!isTRUE(contract$preflight$supported)) {
    not_applicable <- identical(info$id, "linear_sim_embed") ||
      grepl("not_applicable", contract$preflight$status_code, fixed = TRUE)
    return(list(
      intended_to_run = FALSE,
      status = if (not_applicable) "not_applicable" else "unsupported",
      status_code = contract$preflight$status_code,
      status_detail = contract$preflight$status_detail,
      contract = contract
    ))
  }
  envelope <- .alignment_tournament_resource_ceiling_v2()
  envelope <- envelope[envelope$method_id == info$id, , drop = FALSE]
  n_total <- sum(vapply(case$domains, function(x) nrow(x$X), integer(1)))
  if (nrow(envelope) != 1L || n_total > envelope$max_total_nodes) {
    return(list(
      intended_to_run = FALSE,
      status = "unsupported",
      status_code = "unsupported_resource_envelope",
      status_detail = paste0(
        "The frozen ", if (nrow(envelope)) envelope$envelope_id else "missing",
        " resource envelope does not admit ", n_total, " total nodes."
      ),
      contract = contract
    ))
  }
  list(
    intended_to_run = TRUE, status = "pending", status_code = "pending",
    status_detail = NA_character_, contract = contract
  )
}

.alignment_tournament_base_row_v2 <- function(case, method, applicability = NULL) {
  info <- alignment_method_info(method)
  if (is.null(applicability)) {
    applicability <- .alignment_tournament_applicability_v2(info$id, case)
  }
  contract <- applicability$contract
  parameters <- .alignment_tournament_parameters_v2(info$id, case)
  n1 <- nrow(case$domains[[1L]]$X)
  n2 <- nrow(case$domains[[2L]]$X)
  overlap <- if (any(case$truth$source_private_mask) ||
                 any(case$truth$target_private_mask)) {
    "partial_overlap"
  } else {
    "closed_world"
  }
  data.frame(
    tournament_version = .alignment_tournament_version_v2,
    protocol_version = case$protocol_version,
    manifest_hash = unname(.graph_v2_manifest_hashes[[case$split]]),
    case_id = case$case_id,
    case_hash = case$provenance$case_hash,
    split = case$split,
    tier = case$tier,
    family = case$config$family,
    overlap_regime = overlap,
    generation_seed = as.integer(case$config$seed),
    n_source = n1,
    n_target = n2,
    n_total = n1 + n2,
    n_shared = as.integer(case$config$n_shared),
    n_private_source = as.integer(case$config$n_private[[1L]]),
    n_private_target = as.integer(case$config$n_private[[2L]]),
    method_id = info$id,
    method = info$display_name,
    input_class = info$input$primary_class,
    adapter_id = contract$adapter_id,
    adapter_input_hash = contract$adapter_input_hash,
    row_identity_hash = contract$row_identity_hash,
    topology_disposition = contract$topology_disposition,
    observed_topology_used = contract$observed_topology_used,
    observed_edge_count = as.integer(contract$observed_edge_count),
    intended_to_run = isTRUE(applicability$intended_to_run),
    applicability_status = applicability$status,
    applicability_code = applicability$status_code,
    applicability_detail = applicability$status_detail,
    status = applicability$status,
    status_code = applicability$status_code,
    status_detail = applicability$status_detail,
    error_class = NA_character_,
    error_message = NA_character_,
    warning_count = 0L,
    warning_messages = NA_character_,
    parameter_set_id = paste0("graph-v2-", info$id, "-fixed-v1"),
    parameter_set_hash = .graph_v2_hash(parameters),
    parameter_provenance = "frozen_before_final_execution",
    parameter_search_budget = as.integer(contract$parameter_search_budget),
    parameter_trials_used = 0L,
    wall_clock_budget_sec = as.numeric(contract$wall_clock_budget_sec),
    runtime_sec = 0,
    process_elapsed_sec = 0,
    worker_pid = NA_integer_,
    peak_rss_mb = 0,
    memory_measurement = if (isTRUE(applicability$intended_to_run)) {
      "pending_external_process"
    } else {
      "not_run"
    },
    retained_state_mb = 0,
    candidate_edges = 0,
    dense_pairwise_allocated = NA,
    matching_status = if (isTRUE(applicability$intended_to_run)) "pending" else applicability$status,
    matching_n_eligible = sum(!is.na(case$truth$map12)),
    coverage = 0,
    top1_accuracy = NA_real_,
    hits_at_1 = NA_real_,
    hits_at_5 = NA_real_,
    hits_at_10 = NA_real_,
    mrr = NA_real_,
    mrr_at_10 = NA_real_,
    ranking_depth = 0L,
    ranking_complete = FALSE,
    ranking_decoder = NA_character_,
    one_to_one_feasible = NA,
    one_to_one_rate = NA_real_,
    anchor_top1_accuracy = NA_real_,
    nonanchor_top1_accuracy = NA_real_,
    private_rejection_accuracy = NA_real_,
    weighted_structure_conservation = NA_real_,
    prediction_hash = NA_character_,
    classification_status = if (isTRUE(applicability$intended_to_run)) "pending" else applicability$status,
    classification_n_eligible = sum(
      case$domains[[2L]]$node_role %in% c("final", "private")
    ),
    classification_coverage = 0,
    classification_oa = NA_real_,
    classification_balanced_accuracy = NA_real_,
    classification_macro_f1 = NA_real_,
    classification_miou = NA_real_,
    classification_rare_class_recall = NA_real_,
    classification_seen_accuracy = NA_real_,
    classification_unseen_accuracy = NA_real_,
    classification_ece = NA_real_,
    classification_nll = NA_real_,
    classification_decoder = NA_character_,
    converged = NA,
    solver_residual = NA_real_,
    row_marginal_error = NA_real_,
    column_marginal_error = NA_real_,
    objective_start = NA_real_,
    objective_end = NA_real_,
    objective_nonincreasing = NA,
    stringsAsFactors = FALSE
  )
}

.alignment_tournament_top1_from_plan_v2 <- function(plan) {
  if (inherits(plan, "Matrix")) {
    summary <- Matrix::summary(plan)
    out <- rep(NA_integer_, nrow(plan))
    if (!nrow(summary)) return(out)
    summary <- summary[order(summary$i, -summary$x, summary$j), , drop = FALSE]
    best <- summary[!duplicated(summary$i), , drop = FALSE]
    out[best$i] <- best$j
    return(out)
  }
  plan <- as.matrix(plan)
  if (!nrow(plan) || !ncol(plan)) return(integer())
  out <- max.col(plan, ties.method = "first")
  out[rowSums(is.finite(plan)) == 0L] <- NA_integer_
  as.integer(out)
}

.alignment_tournament_embedding_blocks_v2 <- function(scores, n1, n2) {
  if (is.list(scores) && length(scores) >= 2L) {
    left <- as.matrix(scores[[1L]])
    right <- as.matrix(scores[[2L]])
  } else {
    scores <- as.matrix(scores)
    if (nrow(scores) != n1 + n2) return(NULL)
    left <- scores[seq_len(n1), , drop = FALSE]
    right <- scores[n1 + seq_len(n2), , drop = FALSE]
  }
  if (nrow(left) != n1 || nrow(right) != n2 ||
      !ncol(left) || ncol(left) != ncol(right) ||
      any(!is.finite(left)) || any(!is.finite(right))) {
    return(NULL)
  }
  list(source = left, target = right)
}

.alignment_tournament_distance_ranking_v2 <- function(embedding, max_depth = 10L) {
  left <- embedding$source
  right <- embedding$target
  n1 <- nrow(left)
  n2 <- nrow(right)
  complete <- n2 <= 2000L && as.double(n1) * as.double(n2) <= 2e6
  depth <- if (complete) n2 else min(as.integer(max_depth), n2)
  if (!depth) {
    return(list(index = matrix(integer(), n1, 0L), complete = complete))
  }
  if (!complete && requireNamespace("RANN", quietly = TRUE)) {
    index <- RANN::nn2(data = right, query = left, k = depth)$nn.idx
    if (is.null(dim(index))) index <- matrix(index, ncol = depth)
    return(list(index = matrix(as.integer(index), nrow = n1), complete = FALSE))
  }
  left_norm <- rowSums(left * left)
  right_norm <- rowSums(right * right)
  distance <- pmax(outer(left_norm, right_norm, "+") - 2 * tcrossprod(left, right), 0)
  index <- t(apply(distance, 1L, order, method = "radix"))
  if (!is.matrix(index)) index <- matrix(index, nrow = n1)
  list(index = matrix(as.integer(index[, seq_len(depth), drop = FALSE]), nrow = n1), complete = complete)
}

.alignment_tournament_plan_ranking_v2 <- function(plan, max_depth = 10L) {
  n1 <- nrow(plan)
  n2 <- ncol(plan)
  complete <- n2 <= 2000L && as.double(n1) * as.double(n2) <= 2e6
  depth <- if (complete) n2 else min(as.integer(max_depth), n2)
  dense <- as.matrix(plan)
  index <- t(apply(dense, 1L, function(x) {
    head(order(-x, seq_along(x), method = "radix"), depth)
  }))
  if (!is.matrix(index)) index <- matrix(index, nrow = n1)
  list(index = matrix(as.integer(index), nrow = n1), complete = complete)
}

.alignment_tournament_native_first_v2 <- function(prediction, ranking, n_target) {
  prediction <- as.integer(prediction)
  candidates <- ranking$index
  out <- matrix(NA_integer_, nrow(candidates), ncol(candidates))
  for (i in seq_len(nrow(out))) {
    p <- prediction[[i]]
    if (is.na(p) || p < 1L || p > n_target) next
    values <- c(p, candidates[i, candidates[i, ] != p])
    out[i, ] <- head(values, ncol(out))
  }
  list(index = out, complete = ranking$complete)
}

.alignment_tournament_weighted_structure_v2 <- function(case, prediction) {
  source <- case$domains[[1L]]$relations
  target <- case$domains[[2L]]$relations
  n_target <- nrow(case$domains[[2L]]$X)
  prediction <- as.integer(prediction)
  valid_node <- !is.na(prediction) & prediction >= 1L & prediction <= n_target
  keep <- valid_node[source$from] & valid_node[source$to]
  if (!any(keep)) return(NA_real_)
  edge_key <- function(relation, from, to, directed) {
    lo <- ifelse(directed, from, pmin(from, to))
    hi <- ifelse(directed, to, pmax(from, to))
    paste(relation, as.integer(directed), lo, hi, sep = "|")
  }
  target_key <- edge_key(
    target$relation, target$from, target$to, target$directed
  )
  target_weight <- tapply(abs(target$weight), target_key, max)
  mapped_key <- edge_key(
    source$relation[keep], prediction[source$from[keep]],
    prediction[source$to[keep]], source$directed[keep]
  )
  observed <- unname(target_weight[mapped_key])
  observed[is.na(observed)] <- 0
  expected <- abs(source$weight[keep])
  denominator <- sum(expected)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  sum(pmin(expected, observed)) / denominator
}

.alignment_tournament_score_matching_v2 <- function(
  case, prediction, ranking = NULL, ranking_decoder = NA_character_
) {
  normalized <- .graph_benchmark_normalize_prediction(case, prediction)
  if (!isTRUE(normalized$ok)) {
    stop(normalized$detail, call. = FALSE)
  }
  truth <- as.integer(case$truth$map12)
  eligible <- !is.na(truth)
  n_target <- nrow(case$domains[[2L]]$X)
  valid <- !is.na(prediction) & prediction >= 1L & prediction <= n_target
  correct <- eligible & valid & as.integer(prediction) == truth
  anchor <- seq_along(truth) %in% case$anchors$source_row
  subset_accuracy <- function(index) {
    if (!any(index)) return(NA_real_)
    sum(correct[index]) / sum(index)
  }
  valid_targets <- as.integer(prediction[valid])
  one_to_one_rate <- if (length(valid_targets)) {
    length(unique(valid_targets)) / length(valid_targets)
  } else {
    NA_real_
  }
  truth_rank <- rep(Inf, length(truth))
  ranking_depth <- 0L
  ranking_complete <- FALSE
  if (!is.null(ranking) && is.matrix(ranking$index) &&
      nrow(ranking$index) == length(truth)) {
    ranking_depth <- ncol(ranking$index)
    ranking_complete <- isTRUE(ranking$complete)
    for (i in which(eligible & valid)) {
      rank <- match(truth[[i]], ranking$index[i, ])
      if (!is.na(rank)) truth_rank[[i]] <- rank
    }
  }
  reciprocal <- ifelse(is.finite(truth_rank), 1 / truth_rank, 0)
  private <- case$domains[[1L]]$private_mask
  list(
    matching_status = "ok",
    matching_n_eligible = sum(eligible),
    coverage = if (any(eligible)) sum(valid & eligible) / sum(eligible) else 0,
    top1_accuracy = if (any(eligible)) sum(correct) / sum(eligible) else NA_real_,
    hits_at_1 = if (ranking_depth >= 1L && any(eligible)) mean(truth_rank[eligible] <= 1L) else NA_real_,
    hits_at_5 = if (ranking_depth >= 5L && any(eligible)) mean(truth_rank[eligible] <= 5L) else NA_real_,
    hits_at_10 = if (ranking_depth >= 10L && any(eligible)) mean(truth_rank[eligible] <= 10L) else NA_real_,
    mrr = if (ranking_complete && any(eligible)) mean(reciprocal[eligible]) else NA_real_,
    mrr_at_10 = if (ranking_depth >= 10L && any(eligible)) {
      mean(ifelse(truth_rank[eligible] <= 10L, reciprocal[eligible], 0))
    } else {
      NA_real_
    },
    ranking_depth = ranking_depth,
    ranking_complete = ranking_complete,
    ranking_decoder = ranking_decoder,
    one_to_one_feasible = if (length(valid_targets)) !anyDuplicated(valid_targets) else NA,
    one_to_one_rate = one_to_one_rate,
    anchor_top1_accuracy = subset_accuracy(eligible & anchor),
    nonanchor_top1_accuracy = subset_accuracy(eligible & !anchor),
    private_rejection_accuracy = if (any(private)) mean(!valid[private]) else NA_real_,
    weighted_structure_conservation = .alignment_tournament_weighted_structure_v2(case, prediction),
    prediction_hash = normalized$prediction_hash
  )
}

.alignment_tournament_centroid_probabilities_v2 <- function(case, embedding) {
  labels <- as.character(case$domains[[1L]]$labels)
  known <- !is.na(labels) & nzchar(labels)
  classes <- sort(unique(labels[known]))
  if (length(classes) < 2L || sum(known) < 2L) return(NULL)
  centroids <- do.call(rbind, lapply(classes, function(class) {
    colMeans(embedding$source[known & labels == class, , drop = FALSE])
  }))
  left <- rowSums(embedding$target * embedding$target)
  right <- rowSums(centroids * centroids)
  distance <- pmax(outer(left, right, "+") - 2 * tcrossprod(embedding$target, centroids), 0)
  positive <- distance[is.finite(distance) & distance > 0]
  scale <- if (length(positive)) stats::median(positive) else 1
  if (!is.finite(scale) || scale <= .Machine$double.eps) scale <- 1
  logits <- -distance / scale
  logits <- sweep(logits, 1L, apply(logits, 1L, max), "-")
  probabilities <- exp(pmax(logits, -700))
  probabilities <- probabilities / rowSums(probabilities)
  colnames(probabilities) <- classes
  probabilities
}

.alignment_tournament_transport_probabilities_v2 <- function(case, plan) {
  labels <- as.character(case$domains[[1L]]$labels)
  known <- !is.na(labels) & nzchar(labels)
  classes <- sort(unique(labels[known]))
  if (length(classes) < 2L || nrow(plan) != length(labels)) return(NULL)
  indicator <- matrix(0, nrow(plan), length(classes), dimnames = list(NULL, classes))
  indicator[cbind(which(known), match(labels[known], classes))] <- 1
  score <- t(as.matrix(crossprod(indicator, plan)))
  prior <- colSums(indicator) / sum(indicator)
  total <- rowSums(score)
  empty <- !is.finite(total) | total <= 0
  if (any(!empty)) score[!empty, ] <- score[!empty, , drop = FALSE] / total[!empty]
  if (any(empty)) score[empty, ] <- matrix(prior, sum(empty), length(prior), byrow = TRUE)
  colnames(score) <- classes
  score
}

.alignment_tournament_score_classification_v2 <- function(case, probabilities) {
  target <- case$domains[[2L]]
  truth <- as.character(case$truth$domain_labels[[2L]])
  eligible <- target$node_role %in% c("final", "private")
  if (is.null(probabilities) || !is.matrix(probabilities) ||
      nrow(probabilities) != length(truth) || !ncol(probabilities) ||
      is.null(colnames(probabilities))) {
    return(list(
      classification_status = "not_available",
      classification_n_eligible = sum(eligible),
      classification_coverage = 0
    ))
  }
  probabilities[!is.finite(probabilities) | probabilities < 0] <- 0
  total <- rowSums(probabilities)
  valid <- is.finite(total) & total > 0
  probabilities[valid, ] <- probabilities[valid, , drop = FALSE] / total[valid]
  prediction <- rep(NA_character_, length(truth))
  confidence <- numeric(length(truth))
  if (any(valid)) {
    best <- max.col(probabilities[valid, , drop = FALSE], ties.method = "first")
    prediction[valid] <- colnames(probabilities)[best]
    confidence[valid] <- probabilities[valid, , drop = FALSE][cbind(seq_along(best), best)]
  }
  correct <- eligible & valid & prediction == truth
  classes <- sort(unique(c(truth[eligible], prediction[eligible & valid])))
  recall <- f1 <- iou <- rep(NA_real_, length(classes))
  for (j in seq_along(classes)) {
    class <- classes[[j]]
    tp <- sum(eligible & prediction == class & truth == class, na.rm = TRUE)
    fp <- sum(eligible & prediction == class & truth != class, na.rm = TRUE)
    fn <- sum(eligible & truth == class & (is.na(prediction) | prediction != class))
    recall[[j]] <- if (tp + fn) tp / (tp + fn) else NA_real_
    f1[[j]] <- if (2 * tp + fp + fn) 2 * tp / (2 * tp + fp + fn) else NA_real_
    iou[[j]] <- if (tp + fp + fn) tp / (tp + fp + fn) else NA_real_
  }
  source_labels <- as.character(case$domains[[1L]]$labels)
  source_labels <- source_labels[!is.na(source_labels)]
  seen_classes <- unique(source_labels)
  source_counts <- table(source_labels)
  rare_classes <- names(source_counts)[source_counts <= stats::median(source_counts)]
  rare <- eligible & truth %in% rare_classes
  seen <- eligible & truth %in% seen_classes
  unseen <- eligible & !truth %in% seen_classes
  subset_accuracy <- function(index) {
    if (!any(index)) return(NA_real_)
    sum(correct[index]) / sum(index)
  }
  bins <- pmin(10L, floor(pmin(confidence, 1 - .Machine$double.eps) * 10) + 1L)
  ece <- 0
  if (any(eligible)) {
    for (bin in seq_len(10L)) {
      index <- eligible & bins == bin
      if (any(index)) {
        ece <- ece + sum(index) / sum(eligible) *
          abs(mean(correct[index]) - mean(confidence[index]))
      }
    }
  } else {
    ece <- NA_real_
  }
  truth_probability <- rep(.Machine$double.eps, length(truth))
  truth_column <- match(truth, colnames(probabilities))
  has_truth <- valid & !is.na(truth_column)
  truth_probability[has_truth] <- probabilities[cbind(which(has_truth), truth_column[has_truth])]
  list(
    classification_status = "ok",
    classification_n_eligible = sum(eligible),
    classification_coverage = if (any(eligible)) sum(valid & eligible) / sum(eligible) else 0,
    classification_oa = subset_accuracy(eligible),
    classification_balanced_accuracy = mean(recall, na.rm = TRUE),
    classification_macro_f1 = mean(f1, na.rm = TRUE),
    classification_miou = mean(iou, na.rm = TRUE),
    classification_rare_class_recall = subset_accuracy(rare),
    classification_seen_accuracy = subset_accuracy(seen),
    classification_unseen_accuracy = subset_accuracy(unseen),
    classification_ece = ece,
    classification_nll = if (any(eligible)) {
      mean(-log(pmax(truth_probability[eligible], .Machine$double.eps)))
    } else {
      NA_real_
    }
  )
}

.alignment_tournament_label_similarity_v2 <- function(labels) {
  labels <- as.character(labels)
  ok <- !is.na(labels)
  out <- Matrix::sparseMatrix(
    i = integer(), j = integer(), x = numeric(),
    dims = c(length(labels), length(labels))
  )
  if (sum(ok) >= 2L) {
    out[which(ok), which(ok)] <- Matrix::Matrix(
      outer(labels[ok], labels[ok], "==") * 1, sparse = TRUE
    )
  }
  Matrix::diag(out) <- 0
  Matrix::forceSymmetric(out, uplo = "U")
}

.alignment_tournament_nearest_v2 <- function(embedding) {
  as.integer(RANN::nn2(
    data = embedding$target, query = embedding$source, k = 1L
  )$nn.idx[, 1L])
}

.alignment_tournament_plan_diagnostics_v2 <- function(plan) {
  if (is.null(plan)) {
    return(list(row_marginal_error = NA_real_, column_marginal_error = NA_real_))
  }
  row_mass <- as.numeric(Matrix::rowSums(plan))
  column_mass <- as.numeric(Matrix::colSums(plan))
  total <- sum(row_mass)
  if (!is.finite(total) || total <= 0) {
    return(list(row_marginal_error = Inf, column_marginal_error = Inf))
  }
  list(
    row_marginal_error = max(abs(row_mass / total - 1 / length(row_mass))),
    column_marginal_error = max(abs(column_mass / total - 1 / length(column_mass)))
  )
}

.alignment_tournament_objective_trace_v2 <- function(object) {
  candidates <- list(
    object$objective_history,
    object$loss_history,
    object$objectives,
    object$losses,
    object$objective
  )
  for (candidate in candidates) {
    if (is.null(candidate)) next
    if (is.data.frame(candidate)) {
      column <- intersect(c("objective", "loss", "value"), names(candidate))
      if (!length(column)) next
      candidate <- candidate[[column[[1L]]]]
    }
    candidate <- as.numeric(candidate)
    candidate <- candidate[is.finite(candidate)]
    if (length(candidate)) return(candidate)
  }
  numeric()
}

.alignment_tournament_convergence_v2 <- function(object) {
  if (!is.null(object$convergence$status)) {
    return(object$convergence$status %in% c("converged", "tolerance"))
  }
  value <- object$converged
  if (is.null(value)) return(NA)
  if (is.matrix(value) && nrow(value) == ncol(value)) {
    value <- value[row(value) != col(value)]
  }
  value <- as.logical(value)
  value <- value[!is.na(value)]
  if (!length(value)) NA else all(value)
}

.alignment_tournament_adapter_v2 <- function(case, method, parameters) {
  method <- alignment_method_info(method)$id
  hd <- case$hd
  anchors <- case$anchors_name
  n1 <- nrow(case$domains[[1L]]$X)
  n2 <- nrow(case$domains[[2L]]$X)
  result <- switch(
    method,
    rapid_ma = {
      fit <- suppressMessages(rapid_ma(
        lapply(case$domains, `[[`, "X"),
        labels = lapply(case$domains, `[[`, "labels"),
        relations = case$relation_lists,
        positions = lapply(case$domains, `[[`, "positions"),
        attributes = lapply(case$domains, `[[`, "attributes"),
        ncomp = parameters$ncomp,
        control = list(
          max_iter = parameters$max_iter,
          min_iter = 1L,
          seed = as.integer(case$config$seed),
          relation = list(
            build_feature = FALSE,
            custom_symmetrize = if (isTRUE(case$config$directed)) "preserve" else "mutual"
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
        candidate_cap = parameters$candidate_k,
        prototype_bucket_cap = min(
          n2, max(16L, 2L * parameters$candidate_k)
        ),
        candidate_views = parameters$candidate_views,
        view_candidate_k = min(12L, n2),
        assignment = parameters$assignment
      )
      list(fit = fit, matching = matching)
    },
    token_ot_graph_align = {
      control <- token_ot_graph_align_control(
        n_levels = parameters$n_levels,
        prior_mode = parameters$prior_mode,
        views = parameters$views,
        candidate_k = parameters$candidate_k,
        coarsen_method = parameters$coarsen_method,
        token_mode = parameters$token_mode,
        verbose = FALSE
      )
      suppressMessages(token_ot_graph_align(
        hd, anchors = anchors, ncomp = parameters$ncomp, control = control
      ))
    },
    fpgw = suppressMessages(fpgw(
      hd,
      omega1 = parameters$omega1,
      epsilon = parameters$epsilon,
      max_iter = parameters$max_iter,
      inner_max_iter = parameters$inner_max_iter,
      tol = parameters$tol,
      verbose = FALSE
    )),
    gromov_wasserstein = suppressMessages(gromov_wasserstein(
      hd,
      epsilon = parameters$epsilon,
      max_iter = parameters$max_iter,
      inner_max_iter = parameters$inner_max_iter,
      tol = parameters$tol,
      verbose = FALSE
    )),
    spectral_mnn_align = {
      correspondences <- if (nrow(case$anchors)) case$correspondences else NULL
      control <- spectral_mnn_align_control(
        knn = parameters$knn,
        q_descriptors = parameters$q_descriptors,
        mnn_k = parameters$mnn_k,
        max_anchors = max(32L, 4L * max(n1, n2)),
        refine_rounds = parameters$refine_rounds,
        seed = as.integer(case$config$seed),
        verbose = FALSE
      )
      suppressMessages(spectral_mnn_align(
        hd, correspondences = correspondences,
        ncomp = parameters$ncomp, control = control
      ))
    },
    ssma_align = {
      correspondences <- NULL
      if (nrow(case$anchors)) {
        correspondences <- case$correspondences
        correspondences$weight <- 1
        correspondences$source <- "anchor"
      }
      control <- ssma_align_control(
        knn = parameters$knn,
        rank_per_domain = parameters$rank_per_domain,
        solver = parameters$solver,
        verbose = FALSE,
        serial_control = ssma_serial_control(enabled = FALSE)
      )
      suppressMessages(ssma_align(
        hd, correspondences = correspondences,
        ncomp = parameters$ncomp, control = control
      ))
    },
    kema = suppressMessages(kema(
      hd,
      y = anchors,
      ncomp = parameters$ncomp,
      knn = parameters$knn,
      sample_frac = parameters$sample_frac,
      solver = parameters$solver,
      backend = parameters$backend
    )),
    global_geo_align = suppressMessages(global_geo_align(
      hd,
      correspondences = case$correspondences,
      preproc = NULL,
      ncomp = parameters$ncomp,
      control = global_geo_align_control(
        knn = parameters$knn,
        n_landmarks_total = parameters$n_landmarks_total,
        k_embed = parameters$k_embed,
        seed = as.integer(case$config$seed),
        verbose = FALSE
      )
    )),
    lowrank_align = suppressWarnings(suppressMessages(lowrank_align(
      hd,
      y = anchors,
      ncomp = parameters$ncomp,
      simfun = .alignment_tournament_label_similarity_v2,
      mu = parameters$mu,
      lambda = parameters$lambda,
      sv_thresh = parameters$sv_thresh,
      solver = parameters$solver
    ))),
    gpca_align = suppressWarnings(suppressMessages(gpca_align(
      hd,
      y = anchors,
      ncomp = parameters$ncomp,
      simfun = .alignment_tournament_label_similarity_v2,
      u = parameters$u,
      lambda = parameters$lambda,
      control = gpca_align_control(verbose = FALSE)
    ))),
    cone_align = suppressMessages(cone_align(
      hd, anchors = anchors, ncomp = parameters$ncomp
    )),
    grasp = suppressMessages(grasp(
      hd, ncomp = parameters$ncomp,
      q_descriptors = parameters$q_descriptors
    )),
    parrot = suppressMessages(parrot(
      hd, anchors = anchors, ncomp = parameters$ncomp,
      max_iter = parameters$max_iter
    )),
    stop("No runnable tournament adapter for method `", method, "`.", call. = FALSE)
  )

  primary <- if (identical(method, "rapid_ma")) result$fit else result
  plan <- primary[["transport_plan", exact = TRUE]]
  plans <- primary[["transport_plans", exact = TRUE]]
  if (is.null(plan) && !is.null(plans) && length(plans)) {
    plan <- plans[[1L]]
  }
  embedding <- if (identical(method, "rapid_ma")) {
    .alignment_tournament_embedding_blocks_v2(primary$scores, n1, n2)
  } else if (!is.null(primary$s)) {
    .alignment_tournament_embedding_blocks_v2(primary$s, n1, n2)
  } else {
    NULL
  }
  prediction <- if (identical(method, "rapid_ma")) {
    as.integer(result$matching$matches$target)
  } else if (!is.null(plan) && nrow(plan) == n1 && ncol(plan) == n2) {
    .alignment_tournament_top1_from_plan_v2(plan)
  } else if (!is.null(primary$assignment) &&
             length(primary$assignment) == n1) {
    as.integer(primary$assignment)
  } else if (!is.null(embedding)) {
    .alignment_tournament_nearest_v2(embedding)
  } else {
    stop("The adapter returned neither a valid mapping, plan, nor embedding.", call. = FALSE)
  }
  if (length(prediction) != n1) {
    stop(
      "The adapter returned ", length(prediction),
      " predictions for ", n1, " source rows.", call. = FALSE
    )
  }
  objective <- .alignment_tournament_objective_trace_v2(primary)
  list(
    object = result,
    primary = primary,
    prediction = prediction,
    plan = plan,
    embedding = embedding,
    candidate_edges = if (identical(method, "rapid_ma")) {
      result$matching$candidate_edges
    } else if (!is.null(primary$candidate_edges)) {
      as.integer(primary$candidate_edges)
    } else if (!is.null(plan) && inherits(plan, "sparseMatrix")) {
      as.integer(Matrix::nnzero(plan))
    } else if (!is.null(plan)) {
      as.integer(length(plan))
    } else {
      0L
    },
    dense_pairwise_allocated = if (identical(method, "rapid_ma")) {
      isTRUE(result$matching$dense_pairwise_allocated)
    } else if (!is.null(plan)) {
      !inherits(plan, "sparseMatrix")
    } else {
      NA
    },
    converged = .alignment_tournament_convergence_v2(primary),
    objective_start = if (length(objective)) objective[[1L]] else NA_real_,
    objective_end = if (length(objective)) objective[[length(objective)]] else NA_real_,
    objective_nonincreasing = if (length(objective) >= 2L) {
      all(diff(objective) <= 1e-8 * pmax(1, abs(objective[-length(objective)])))
    } else {
      NA
    }
  )
}

.alignment_tournament_fill_row_v2 <- function(row, values) {
  for (name in intersect(names(values), names(row))) {
    row[[name]] <- values[[name]]
  }
  row
}

.alignment_tournament_run_case_v2 <- function(case, method) {
  applicability <- .alignment_tournament_applicability_v2(method, case)
  row <- .alignment_tournament_base_row_v2(case, method, applicability)
  if (!isTRUE(applicability$intended_to_run)) return(row)
  parameters <- .alignment_tournament_parameters_v2(method, case)
  warnings <- character()
  started <- proc.time()[[3L]]
  adapter <- tryCatch(
    withCallingHandlers(
      .alignment_tournament_adapter_v2(case, method, parameters),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(error) error
  )
  runtime <- proc.time()[[3L]] - started
  row$runtime_sec <- as.numeric(runtime)
  row$parameter_trials_used <- 1L
  row$warning_count <- length(warnings)
  row$warning_messages <- if (length(warnings)) {
    paste(unique(warnings), collapse = " | ")
  } else {
    NA_character_
  }
  if (inherits(adapter, "error")) {
    row$status <- "error"
    row$status_code <- "execution_error"
    row$status_detail <- conditionMessage(adapter)
    row$error_class <- class(adapter)[[1L]]
    row$error_message <- conditionMessage(adapter)
    row$matching_status <- "error"
    row$classification_status <- "not_run"
    return(row)
  }

  n_target <- nrow(case$domains[[2L]]$X)
  ranking <- NULL
  ranking_decoder <- NA_character_
  if (!is.null(adapter$plan) &&
      nrow(adapter$plan) == length(adapter$prediction) &&
      ncol(adapter$plan) == n_target) {
    ranking <- .alignment_tournament_plan_ranking_v2(adapter$plan)
    ranking_decoder <- "native_first_transport"
  } else if (!is.null(adapter$embedding)) {
    ranking <- .alignment_tournament_distance_ranking_v2(adapter$embedding)
    ranking_decoder <- "native_first_common_embedding"
  }
  if (!is.null(ranking)) {
    ranking <- .alignment_tournament_native_first_v2(
      adapter$prediction, ranking, n_target
    )
  }
  matching <- tryCatch(
    .alignment_tournament_score_matching_v2(
      case, adapter$prediction, ranking, ranking_decoder
    ),
    error = function(error) error
  )
  if (inherits(matching, "error")) {
    row$status <- "error"
    row$status_code <- "malformed_output"
    row$status_detail <- conditionMessage(matching)
    row$error_class <- class(matching)[[1L]]
    row$error_message <- conditionMessage(matching)
    row$matching_status <- "error"
    row$classification_status <- "not_run"
    return(row)
  }

  probabilities <- NULL
  classification_decoder <- NA_character_
  if (!is.null(adapter$embedding)) {
    probabilities <- .alignment_tournament_centroid_probabilities_v2(
      case, adapter$embedding
    )
    classification_decoder <- "source_label_common_centroid_v1"
  } else if (!is.null(adapter$plan) &&
             nrow(adapter$plan) == nrow(case$domains[[1L]]$X) &&
             ncol(adapter$plan) == n_target) {
    probabilities <- .alignment_tournament_transport_probabilities_v2(
      case, adapter$plan
    )
    classification_decoder <- "source_label_transport_v1"
  }
  classification <- .alignment_tournament_score_classification_v2(
    case, probabilities
  )
  classification$classification_decoder <- if (is.null(probabilities)) {
    NA_character_
  } else {
    classification_decoder
  }
  marginal <- .alignment_tournament_plan_diagnostics_v2(adapter$plan)
  row <- .alignment_tournament_fill_row_v2(row, matching)
  row <- .alignment_tournament_fill_row_v2(row, classification)
  row <- .alignment_tournament_fill_row_v2(row, marginal)
  row$status <- "ok"
  row$status_code <- "ok"
  row$status_detail <- NA_character_
  row$applicability_status <- "applicable"
  row$applicability_code <- "ok"
  row$memory_measurement <- "external_process_required"
  row$retained_state_mb <- as.numeric(utils::object.size(adapter$object)) / 1024^2
  row$candidate_edges <- as.integer(adapter$candidate_edges)
  row$dense_pairwise_allocated <- as.logical(adapter$dense_pairwise_allocated)
  row$converged <- as.logical(adapter$converged)
  row$objective_start <- adapter$objective_start
  row$objective_end <- adapter$objective_end
  row$objective_nonincreasing <- as.logical(adapter$objective_nonincreasing)
  row
}

.alignment_tournament_failure_row_v2 <- function(
  case, method, status, status_code, status_detail,
  runtime_sec, peak_rss_mb = NA_real_, memory_measurement = "external_process"
) {
  applicability <- .alignment_tournament_applicability_v2(method, case)
  row <- .alignment_tournament_base_row_v2(case, method, applicability)
  row$parameter_trials_used <- if (isTRUE(applicability$intended_to_run)) 1L else 0L
  row$status <- status
  row$status_code <- status_code
  row$status_detail <- status_detail
  row$error_class <- status
  row$error_message <- status_detail
  row$matching_status <- status
  row$classification_status <- "not_run"
  row$runtime_sec <- as.numeric(runtime_sec)
  row$process_elapsed_sec <- as.numeric(runtime_sec)
  row$peak_rss_mb <- as.numeric(peak_rss_mb)
  row$memory_measurement <- memory_measurement
  row
}

.alignment_tournament_required_columns_v2 <- function() {
  names(.alignment_tournament_base_row_v2(
    synthetic_graph_alignment_case_v2(
      n_shared = 8L, family = "sbm", attr_dim = 2L,
      n_anchors = 2L, split = "development", tier = "oracle",
      tuning_budget = 1L, seed = 9182L,
      case_id = "tournament-schema-probe"
    ),
    "rapid_ma"
  ))
}

.alignment_tournament_validate_results_v2 <- function(
  results,
  manifest,
  methods = .alignment_tournament_methods_v2(),
  require_external_memory = TRUE
) {
  problems <- data.frame(
    code = character(), detail = character(), stringsAsFactors = FALSE
  )
  add <- function(code, detail) {
    problems <<- rbind(
      problems,
      data.frame(code = code, detail = detail, stringsAsFactors = FALSE)
    )
  }
  required <- .alignment_tournament_required_columns_v2()
  missing <- setdiff(required, names(results))
  if (length(missing)) {
    add("missing_columns", paste(missing, collapse = ", "))
    return(problems)
  }
  expected <- expand.grid(
    case_id = manifest$case_id,
    method_id = methods,
    stringsAsFactors = FALSE
  )
  expected_key <- paste(expected$case_id, expected$method_id, sep = "\r")
  observed_key <- paste(results$case_id, results$method_id, sep = "\r")
  if (anyDuplicated(observed_key)) add("duplicate_rows", "Method-case keys are not unique.")
  missing_key <- setdiff(expected_key, observed_key)
  extra_key <- setdiff(observed_key, expected_key)
  if (length(missing_key)) {
    add("missing_rows", paste0(length(missing_key), " intended grid rows are absent."))
  }
  if (length(extra_key)) {
    add("extra_rows", paste0(length(extra_key), " rows are outside the intended grid."))
  }
  allowed_status <- c(
    "ok", "unsupported", "not_applicable", "not_available",
    "error", "timeout", "crash"
  )
  if (any(!results$status %in% allowed_status)) {
    add("invalid_status", "Rows contain pending or unknown status values.")
  }
  intended <- results$intended_to_run
  if (any(intended & results$status %in% c("unsupported", "not_applicable"))) {
    add("intent_contradiction", "Runnable rows cannot finish as semantic non-applicability.")
  }
  if (any(!is.finite(results$runtime_sec) | results$runtime_sec < 0)) {
    add("runtime_missing", "Every row must carry a finite non-negative runtime.")
  }
  if (any(!nzchar(results$parameter_set_hash)) ||
      any(!nzchar(results$adapter_input_hash)) ||
      any(!nzchar(results$case_hash))) {
    add("provenance_missing", "Parameter, adapter, and case hashes are mandatory.")
  }
  if (any(results$split == "final" & results$parameter_search_budget != 0L)) {
    add("final_tuning", "Final rows must have zero parameter-search budget.")
  }
  if (any(intended & results$parameter_trials_used != 1L)) {
    add("trial_count", "Every intended row must record exactly one frozen trial.")
  }
  if (any(!intended & results$parameter_trials_used != 0L)) {
    add("unsupported_trial", "Non-runnable rows must not claim a solver trial.")
  }
  ok <- results$status == "ok"
  bounded <- c(
    "coverage", "top1_accuracy", "hits_at_1", "hits_at_5",
    "hits_at_10", "mrr", "mrr_at_10", "one_to_one_rate",
    "anchor_top1_accuracy", "nonanchor_top1_accuracy",
    "private_rejection_accuracy", "weighted_structure_conservation",
    "classification_coverage", "classification_oa",
    "classification_balanced_accuracy", "classification_macro_f1",
    "classification_miou", "classification_rare_class_recall",
    "classification_seen_accuracy", "classification_unseen_accuracy",
    "classification_ece"
  )
  for (column in bounded) {
    value <- results[[column]]
    finite <- is.finite(value)
    if (any(finite & (value < -1e-12 | value > 1 + 1e-12))) {
      add("metric_bounds", paste0(column, " leaves [0, 1]."))
    }
  }
  if (any(ok & results$matching_status != "ok")) {
    add("successful_matching", "Successful rows require a scored matching result.")
  }
  if (any(ok & !is.finite(results$top1_accuracy))) {
    add("successful_accuracy", "Successful rows require finite top-1 accuracy.")
  }
  failed <- intended & !ok
  if (any(failed & is.finite(results$top1_accuracy))) {
    add("failed_accuracy", "Failed/timeout/crash rows must not retain accuracy values.")
  }
  if (isTRUE(require_external_memory)) {
    if (any(intended & (!is.finite(results$peak_rss_mb) |
                        results$peak_rss_mb <= 0))) {
      add("peak_rss_missing", "Every intended row needs positive external peak RSS.")
    }
    if (any(intended & results$memory_measurement != "external_process_peak_rss")) {
      add("memory_provenance", "Intended rows require external-process RSS provenance.")
    }
  }
  rownames(problems) <- NULL
  problems
}

.alignment_tournament_mean_ci_v2 <- function(x, level = 0.95) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA, low = NA, high = NA, sd = NA))
  average <- mean(x)
  deviation <- if (length(x) >= 2L) stats::sd(x) else NA_real_
  if (length(x) < 2L || !is.finite(deviation)) {
    return(c(mean = average, low = NA, high = NA, sd = deviation))
  }
  critical <- stats::qt(1 - (1 - level) / 2, df = length(x) - 1L)
  error <- critical * deviation / sqrt(length(x))
  c(mean = average, low = average - error, high = average + error, sd = deviation)
}

.alignment_tournament_summary_by_regime_v2 <- function(results) {
  key_columns <- c("tier", "family", "overlap_regime", "method_id", "method")
  groups <- split(results, interaction(results[key_columns], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(group) {
    intended <- group$intended_to_run
    ok <- intended & group$status == "ok"
    intent_top1 <- ifelse(ok, group$top1_accuracy, 0)
    intent_classification <- ifelse(
      ok & group$classification_status == "ok", group$classification_oa, 0
    )
    top1_ci <- .alignment_tournament_mean_ci_v2(intent_top1[intended])
    data.frame(
      tier = group$tier[[1L]],
      family = group$family[[1L]],
      overlap_regime = group$overlap_regime[[1L]],
      method_id = group$method_id[[1L]],
      method = group$method[[1L]],
      n_rows = nrow(group),
      n_intended = sum(intended),
      n_ok = sum(ok),
      n_error = sum(intended & group$status == "error"),
      n_timeout = sum(intended & group$status == "timeout"),
      n_crash = sum(intended & group$status == "crash"),
      completion_rate = if (any(intended)) sum(ok) / sum(intended) else NA_real_,
      conditional_top1 = if (any(ok)) mean(group$top1_accuracy[ok]) else NA_real_,
      intent_top1 = unname(top1_ci[["mean"]]),
      intent_top1_ci_low = unname(top1_ci[["low"]]),
      intent_top1_ci_high = unname(top1_ci[["high"]]),
      intent_top1_sd = unname(top1_ci[["sd"]]),
      conditional_mrr = if (any(ok & is.finite(group$mrr))) {
        mean(group$mrr[ok & is.finite(group$mrr)])
      } else {
        NA_real_
      },
      conditional_coverage = if (any(ok)) mean(group$coverage[ok]) else NA_real_,
      conditional_classification_oa = if (any(ok & group$classification_status == "ok")) {
        mean(group$classification_oa[ok & group$classification_status == "ok"])
      } else {
        NA_real_
      },
      intent_classification_oa = if (any(intended)) {
        mean(intent_classification[intended])
      } else {
        NA_real_
      },
      median_runtime_sec = if (any(intended)) stats::median(group$runtime_sec[intended]) else NA_real_,
      median_peak_rss_mb = if (any(intended)) stats::median(group$peak_rss_mb[intended]) else NA_real_,
      max_peak_rss_mb = if (any(intended)) max(group$peak_rss_mb[intended]) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$tier, out$overlap_regime, out$family, out$method_id), , drop = FALSE]
}

.alignment_tournament_paired_vs_rapid_v2 <- function(results) {
  regime <- interaction(
    results$tier, results$family, results$overlap_regime,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(results, regime)
  rows <- list()
  pointer <- 0L
  for (group in groups) {
    rapid <- group[group$method_id == "rapid_ma", , drop = FALSE]
    if (!nrow(rapid)) next
    for (method in setdiff(unique(group$method_id), "rapid_ma")) {
      comparator <- group[group$method_id == method, , drop = FALSE]
      pair <- merge(
        rapid, comparator, by = "case_id", suffixes = c("_rapid", "_method")
      )
      pair <- pair[pair$intended_to_run_rapid & pair$intended_to_run_method, , drop = FALSE]
      if (!nrow(pair)) next
      rapid_score <- ifelse(pair$status_rapid == "ok", pair$top1_accuracy_rapid, 0)
      method_score <- ifelse(pair$status_method == "ok", pair$top1_accuracy_method, 0)
      difference <- method_score - rapid_score
      ci <- .alignment_tournament_mean_ci_v2(difference)
      deviation <- stats::sd(difference)
      pointer <- pointer + 1L
      rows[[pointer]] <- data.frame(
        tier = pair$tier_rapid[[1L]],
        family = pair$family_rapid[[1L]],
        overlap_regime = pair$overlap_regime_rapid[[1L]],
        method_id = method,
        reference_method_id = "rapid_ma",
        n_paired = nrow(pair),
        mean_difference = unname(ci[["mean"]]),
        difference_ci_low = unname(ci[["low"]]),
        difference_ci_high = unname(ci[["high"]]),
        cohen_dz = if (is.finite(deviation) && deviation > 0) {
          mean(difference) / deviation
        } else {
          NA_real_
        },
        wins = sum(difference > 0),
        ties = sum(difference == 0),
        losses = sum(difference < 0),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(
      tier = character(), family = character(), overlap_regime = character(),
      method_id = character(), reference_method_id = character(),
      n_paired = integer(), mean_difference = numeric(),
      difference_ci_low = numeric(), difference_ci_high = numeric(),
      cohen_dz = numeric(), wins = integer(), ties = integer(), losses = integer(),
      stringsAsFactors = FALSE
    ))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.alignment_tournament_pareto_v2 <- function(summary) {
  regime <- interaction(
    summary$tier, summary$family, summary$overlap_regime,
    drop = TRUE, lex.order = TRUE
  )
  groups <- split(summary, regime)
  rows <- lapply(groups, function(group) {
    group <- group[group$n_intended > 0L, , drop = FALSE]
    if (!nrow(group)) return(NULL)
    group$pareto <- TRUE
    for (i in seq_len(nrow(group))) {
      for (j in seq_len(nrow(group))) {
        if (i == j) next
        no_worse <-
          group$intent_top1[[j]] >= group$intent_top1[[i]] &&
          group$completion_rate[[j]] >= group$completion_rate[[i]] &&
          group$median_runtime_sec[[j]] <= group$median_runtime_sec[[i]] &&
          group$median_peak_rss_mb[[j]] <= group$median_peak_rss_mb[[i]]
        strictly_better <-
          group$intent_top1[[j]] > group$intent_top1[[i]] ||
          group$completion_rate[[j]] > group$completion_rate[[i]] ||
          group$median_runtime_sec[[j]] < group$median_runtime_sec[[i]] ||
          group$median_peak_rss_mb[[j]] < group$median_peak_rss_mb[[i]]
        if (isTRUE(no_worse && strictly_better)) {
          group$pareto[[i]] <- FALSE
          break
        }
      }
    }
    group[, c(
      "tier", "family", "overlap_regime", "method_id", "n_intended",
      "completion_rate", "intent_top1", "median_runtime_sec",
      "median_peak_rss_mb", "pareto"
    ), drop = FALSE]
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.alignment_tournament_scaling_v2 <- function(results) {
  rows <- lapply(split(results, results$method_id), function(group) {
    group <- group[
      group$intended_to_run & group$status == "ok" &
        is.finite(group$runtime_sec) & group$runtime_sec > 0,
      , drop = FALSE
    ]
    if (nrow(group) < 3L || length(unique(group$n_total)) < 2L) return(NULL)
    fit <- stats::lm(log(runtime_sec) ~ log(n_total), data = group)
    data.frame(
      method_id = group$method_id[[1L]],
      n_cases = nrow(group),
      min_total_nodes = min(group$n_total),
      max_total_nodes = max(group$n_total),
      runtime_scaling_exponent = unname(stats::coef(fit)[[2L]]),
      adjusted_r_squared = summary(fit)$adj.r.squared,
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.alignment_tournament_summarize_v2 <- function(results) {
  summary <- .alignment_tournament_summary_by_regime_v2(results)
  list(
    by_regime = summary,
    paired_vs_rapid_ma = .alignment_tournament_paired_vs_rapid_v2(results),
    pareto_by_regime = .alignment_tournament_pareto_v2(summary),
    scaling = .alignment_tournament_scaling_v2(results)
  )
}

.alignment_tournament_parameter_ledger_v2 <- function(
  manifest,
  methods = .alignment_tournament_methods_v2()
) {
  rows <- list()
  pointer <- 0L
  for (case_id in manifest$case_id) {
    case <- graph_alignment_case_from_manifest_v2(case_id, manifest$split[[1L]])
    for (method in methods) {
      parameters <- .alignment_tournament_parameters_v2(method, case)
      applicability <- .alignment_tournament_applicability_v2(method, case)
      pointer <- pointer + 1L
      rows[[pointer]] <- data.frame(
        case_id = case_id,
        method_id = method,
        parameter_set_id = paste0("graph-v2-", method, "-fixed-v1"),
        parameter_set_hash = .graph_v2_hash(parameters),
        parameters = paste(capture.output(dput(parameters)), collapse = ""),
        intended_to_run = applicability$intended_to_run,
        applicability_code = applicability$status_code,
        wall_clock_budget_sec = applicability$contract$wall_clock_budget_sec,
        parameter_search_budget = applicability$contract$parameter_search_budget,
        provenance = "frozen_before_final_execution",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
