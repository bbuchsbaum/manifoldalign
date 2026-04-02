#' Cross-Validated Row-Wise Alignment Scoring
#'
#' Runs explicit row-index cross-validation for multi-domain alignment models by
#' delegating fold construction and execution to `multidesign::cv_rows()` and
#' `multidesign::cross_validate()`. For each fold, held-out rows are scored by
#' how well their cross-block neighbours agree in an external feature space,
#' using either latent [oos_predict()] projection or method-specific
#' `predict(..., type = "weights")` support.
#'
#' This is designed for SSMA-style methods that return a
#' `multiblock_biprojector`, but it also works for transport-style fit objects
#' that expose weight-based prediction against training targets.
#'
#' @param data A `hyperdesign`-compatible object. If `data` is not already a
#'   `multidesign` hyperdesign, it is coerced to one internally.
#' @param rows Explicit held-out row specification forwarded to
#'   `multidesign::cv_rows()`. For hyperdesigns, each fold should be a named list
#'   mapping block names (or positions) to held-out row indices.
#' @param fit_fn Function taking the analysis split for a fold and returning a
#'   fitted alignment object.
#' @param features External feature matrices aligned to the original rows of
#'   each block. Supply either:
#'   1) a named list of matrices, one per block, or
#'   2) a hyperdesign-compatible object whose `$x` matrices are treated as
#'   features.
#' @param k Positive integer number of latent nearest neighbours used when
#'   scoring held-out rows.
#' @param feature_similarity Similarity used for the external feature space.
#'   Either `"cosine"` (default) or `"correlation"`.
#' @param target_pool Which rows are available as retrieval targets for each
#'   held-out query block: `"analysis"` (default) uses the training rows from
#'   the fitted fold, `"assessment"` uses only other held-out rows, and `"both"`
#'   concatenates the two pools when the fit supports latent out-of-sample
#'   projection.
#' @param prediction_mode How to obtain cross-block neighbour rankings. `"auto"`
#'   (default) uses latent `oos_predict()` projection for embedding-style fits
#'   and prefers `predict(..., type = "weights")` for transport-style fits when
#'   `target_pool = "analysis"`. `"embedding"` forces latent neighbour search
#'   and `"weights"` forces weight-based ranking against training targets.
#'
#' @return A `multidesign::cv_result` with one row per fold and the following
#'   metrics:
#'   \itemize{
#'     \item `mean_top1_similarity`: mean feature similarity of the closest
#'       latent neighbour across directed held-out block pairs.
#'     \item `mean_topk_similarity`: mean feature similarity of the top-`k`
#'       latent neighbours across directed held-out block pairs.
#'     \item `oracle_top1_similarity`: best achievable top-1 feature similarity
#'       for the same held-out rows when latent geometry is ignored.
#'     \item `oracle_topk_similarity`: best achievable top-`k` feature
#'       similarity for the same held-out rows when latent geometry is ignored.
#'     \item `top1_gap` / `topk_gap`: oracle minus observed similarity.
#'     \item `n_queries`: total number of held-out query rows scored.
#'     \item `n_pairs`: number of directed block pairs scored in the fold.
#'   }
#'
#' @details
#' The external feature matrices are never used for fitting. They are only used
#' for held-out evaluation, which keeps the score honest as long as the
#' matching held-out rows are excluded from the training split. If your fitting
#' procedure uses explicit correspondence tables keyed by original rows, include
#' a stable row identifier in each block's design before calling this function.
#'
#' @export
cv_alignment_rows <- function(data,
                              rows,
                              fit_fn,
                              features,
                              k = 5L,
                              feature_similarity = c("cosine", "correlation"),
                              target_pool = c("analysis", "assessment", "both"),
                              prediction_mode = c("auto", "embedding", "weights")) {
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop("cv_alignment_rows() requires the 'multidesign' package.", call. = FALSE)
  }
  if (!is.function(fit_fn)) {
    stop("`fit_fn` must be a function.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }
  feature_similarity <- match.arg(feature_similarity)
  target_pool <- match.arg(target_pool)
  prediction_mode <- match.arg(prediction_mode)
  if (identical(prediction_mode, "weights") && !identical(target_pool, "analysis")) {
    stop("`prediction_mode = \"weights\"` requires `target_pool = \"analysis\"`.", call. = FALSE)
  }

  hd <- .cv_alignment_as_multidesign_hyperdesign(data)
  feature_blocks <- .cv_alignment_resolve_feature_blocks(features, reference = hd)

  folds <- multidesign::cv_rows(hd, rows = rows)

  multidesign::cross_validate(
    folds,
    fit_fn = fit_fn,
    score_fn = function(model, assessment, fold, fold_id) {
      .cv_alignment_score_fold(
        fit = model,
        fold = fold,
        feature_blocks = feature_blocks,
        k = as.integer(k),
        feature_similarity = feature_similarity,
        target_pool = target_pool,
        prediction_mode = prediction_mode,
        fold_id = fold_id
      )
    }
  )
}

.cv_alignment_as_multidesign_hyperdesign <- function(data) {
  if (!inherits(data, "hyperdesign")) {
    stop("`data` must be a hyperdesign-compatible object.", call. = FALSE)
  }

  if (!is.null(data$blocks) && all(vapply(data$blocks, inherits, logical(1), what = "multidesign"))) {
    return(data)
  }

  resolved <- resolve_hyperdesign(data)
  domains_raw <- if (!is.null(data$blocks)) data$blocks else unclass(data)

  blocks <- lapply(seq_along(resolved$domains), function(i) {
    dom_raw <- domains_raw[[i]]
    dom <- resolved$domains[[i]]

    if (inherits(dom_raw, "multidesign")) {
      return(dom_raw)
    }

    design <- dom$design
    if (is.null(design)) {
      design <- data.frame(row_id = seq_len(nrow(dom$x)))
    } else {
      design <- as.data.frame(design)
    }

    column_design <- if (is.list(dom_raw) && !is.null(dom_raw$column_design)) {
      dom_raw$column_design
    } else {
      NULL
    }

    multidesign::multidesign(as.matrix(dom$x), design, column_design = column_design)
  })

  multidesign::hyperdesign(blocks, block_names = resolved$domain_names)
}

.cv_alignment_resolve_feature_blocks <- function(features, reference) {
  ref <- resolve_hyperdesign(reference)
  domain_names <- ref$domain_names
  domain_sizes <- vapply(ref$domains, function(dom) nrow(dom$x), integer(1))

  feature_domains <- NULL
  if (inherits(features, "hyperdesign")) {
    feature_domains <- resolve_hyperdesign(features)$domains
  } else if (is.list(features) && length(features) > 0 &&
             all(vapply(features, function(x) {
               inherits(x, "multidesign") || (is.list(x) && !is.null(x$x))
             }, logical(1)))) {
    feature_domains <- lapply(features, function(block) {
      if (inherits(block, "multidesign")) {
        list(x = block$x)
      } else {
        list(x = block$x)
      }
    })
  }

  if (!is.null(feature_domains)) {
    features <- lapply(feature_domains, function(dom) as.matrix(dom$x))
  }

  if (!is.list(features) || length(features) != length(domain_names)) {
    stop(
      "`features` must be a list with one feature matrix per block, or a hyperdesign-compatible object of the same length as `data`.",
      call. = FALSE
    )
  }

  out <- vector("list", length(domain_names))
  names(out) <- domain_names

  feature_names <- names(features)
  if (!is.null(feature_names) && all(nzchar(feature_names))) {
    missing_blocks <- setdiff(domain_names, feature_names)
    if (length(missing_blocks)) {
      stop("`features` is missing block(s): ", paste(missing_blocks, collapse = ", "), ".", call. = FALSE)
    }
    for (nm in domain_names) {
      out[[nm]] <- as.matrix(features[[nm]])
    }
  } else {
    for (i in seq_along(domain_names)) {
      out[[i]] <- as.matrix(features[[i]])
    }
  }

  for (i in seq_along(out)) {
    if (!is.numeric(out[[i]])) {
      stop("All feature blocks must be numeric matrices.", call. = FALSE)
    }
    if (nrow(out[[i]]) != domain_sizes[[i]]) {
      stop(
        "Feature block `", domain_names[[i]], "` has ", nrow(out[[i]]),
        " rows but the corresponding data block has ", domain_sizes[[i]], ".",
        call. = FALSE
      )
    }
  }

  out
}

.cv_alignment_score_fold <- function(fit,
                                     fold,
                                     feature_blocks,
                                     k,
                                     feature_similarity,
                                     target_pool,
                                     prediction_mode,
                                     fold_id = NULL) {
  if (!is.list(fold) || is.null(fold$analysis) || is.null(fold$assessment) || is.null(fold$held_out)) {
    stop("Fold is missing analysis/assessment/held_out components from multidesign::cv_rows().", call. = FALSE)
  }

  all_block_names <- names(feature_blocks)
  if (length(all_block_names) < 2L) {
    stop("cv_alignment_rows() requires at least two blocks for cross-block scoring.", call. = FALSE)
  }

  assessment_resolved <- resolve_hyperdesign(fold$assessment)
  analysis_resolved <- resolve_hyperdesign(fold$analysis)
  query_block_names <- assessment_resolved$domain_names

  if (!length(query_block_names)) {
    stop("No held-out query rows were available in this fold.", call. = FALSE)
  }

  row_sets <- .cv_alignment_row_sets(feature_blocks, fold$held_out)
  assessment_projection_cache <- new.env(parent = emptyenv())
  analysis_score_cache <- NULL
  mode_order <- .cv_alignment_mode_order(fit, target_pool, prediction_mode)

  pair_metrics <- list()
  pair_idx <- 0L

  for (query_block in query_block_names) {
    Xi <- .cv_alignment_block_matrix(assessment_resolved, query_block)
    Fi <- feature_blocks[[query_block]][row_sets[[query_block]]$assessment, , drop = FALSE]

    if (!nrow(Xi) || !nrow(Fi)) {
      next
    }

    target_blocks <- setdiff(all_block_names, query_block)
    for (target_block in target_blocks) {
      sim_mat <- NULL
      nn_idx <- NULL
      errors <- character(0)

      for (mode in mode_order) {
        if (identical(mode, "weights")) {
          weight_res <- tryCatch(
            .cv_alignment_predict_target_weights(
              fit = fit,
              newdata = Xi,
              from = query_block,
              to = target_block,
              block_names = all_block_names,
              fold_id = fold_id
            ),
            error = function(e) e
          )

          if (inherits(weight_res, "error")) {
            errors <- c(errors, paste0("weights: ", conditionMessage(weight_res)))
            if (identical(prediction_mode, "weights")) {
              stop(conditionMessage(weight_res), call. = FALSE)
            }
            next
          }

          target_features <- feature_blocks[[target_block]][row_sets[[target_block]]$analysis, , drop = FALSE]
          if (!nrow(target_features)) {
            next
          }
          if (ncol(weight_res) != nrow(target_features)) {
            stop(
              "Weight prediction for block `", query_block, "` -> `", target_block,
              "` returned ", ncol(weight_res), " target columns but analysis target pool has ",
              nrow(target_features), " rows.",
              call. = FALSE
            )
          }

          sim_mat <- .cv_alignment_similarity_matrix(Fi, target_features, method = feature_similarity)
          k_use <- min(as.integer(k), ncol(weight_res))
          nn_idx <- .cv_alignment_topk_weight_indices(weight_res, k = k_use)
          break
        }

        embedding_res <- tryCatch({
          if (is.null(analysis_score_cache)) {
            analysis_score_cache <- .cv_alignment_extract_analysis_scores(
              fit = fit,
              analysis = analysis_resolved
            )
          }
          query_scores <- .cv_alignment_get_assessment_projection(
            fit = fit,
            assessment = assessment_resolved,
            block_name = query_block,
            block_names = all_block_names,
            cache = assessment_projection_cache,
            fold_id = fold_id
          )
          target_data <- .cv_alignment_target_pool_embeddings(
            fit = fit,
            analysis = analysis_resolved,
            assessment = assessment_resolved,
            feature_blocks = feature_blocks,
            row_sets = row_sets,
            target_block = target_block,
            target_pool = target_pool,
            analysis_scores = analysis_score_cache,
            projection_cache = assessment_projection_cache,
            block_names = all_block_names,
            fold_id = fold_id
          )
          list(
            query_scores = query_scores,
            target_scores = target_data$scores,
            target_features = target_data$features
          )
        }, error = function(e) e)

        if (inherits(embedding_res, "error")) {
          errors <- c(errors, paste0("embedding: ", conditionMessage(embedding_res)))
          if (identical(prediction_mode, "embedding")) {
            stop(conditionMessage(embedding_res), call. = FALSE)
          }
          next
        }

        if (!nrow(embedding_res$query_scores) || !nrow(embedding_res$target_scores)) {
          next
        }

        sim_mat <- .cv_alignment_similarity_matrix(
          Fi,
          embedding_res$target_features,
          method = feature_similarity
        )
        k_use <- min(as.integer(k), nrow(embedding_res$target_scores))
        nn_idx <- .cv_alignment_knn_indices(
          embedding_res$query_scores,
          embedding_res$target_scores,
          k = k_use
        )
        break
      }

      if (is.null(sim_mat) || is.null(nn_idx)) {
        if (length(errors)) {
          stop(
            "Failed to score block `", query_block, "` against `", target_block,
            "` in fold ", if (is.null(fold_id)) "<unknown>" else fold_id,
            ". Attempted modes: ", paste(errors, collapse = " ; "),
            call. = FALSE
          )
        }
        next
      }

      if (is.null(dim(nn_idx))) {
        nn_idx <- matrix(nn_idx, ncol = 1L)
      }

      top1 <- sim_mat[cbind(seq_len(nrow(sim_mat)), nn_idx[, 1L])]
      topk <- vapply(seq_len(nrow(sim_mat)), function(r) {
        mean(sim_mat[r, nn_idx[r, ], drop = TRUE])
      }, numeric(1))
      oracle_top1 <- apply(sim_mat, 1L, max)
      oracle_topk <- apply(sim_mat, 1L, function(x) {
        mean(sort(x, decreasing = TRUE)[seq_len(ncol(nn_idx))])
      })

      pair_idx <- pair_idx + 1L
      pair_metrics[[pair_idx]] <- c(
        top1_sum = sum(top1),
        topk_sum = sum(topk),
        oracle_top1_sum = sum(oracle_top1),
        oracle_topk_sum = sum(oracle_topk),
        n_queries = nrow(sim_mat),
        n_pairs = 1
      )
    }
  }

  pair_metrics <- Filter(Negate(is.null), pair_metrics)
  if (!length(pair_metrics)) {
    stop("No valid cross-block held-out comparisons were available for scoring.", call. = FALSE)
  }

  totals <- Reduce(`+`, pair_metrics)
  n_queries <- unname(totals[["n_queries"]])
  n_pairs <- unname(totals[["n_pairs"]])

  mean_top1 <- unname(totals[["top1_sum"]] / n_queries)
  mean_topk <- unname(totals[["topk_sum"]] / n_queries)
  oracle_top1 <- unname(totals[["oracle_top1_sum"]] / n_queries)
  oracle_topk <- unname(totals[["oracle_topk_sum"]] / n_queries)

  c(
    mean_top1_similarity = mean_top1,
    mean_topk_similarity = mean_topk,
    oracle_top1_similarity = oracle_top1,
    oracle_topk_similarity = oracle_topk,
    top1_gap = oracle_top1 - mean_top1,
    topk_gap = oracle_topk - mean_topk,
    n_queries = n_queries,
    n_pairs = n_pairs
  )
}

.cv_alignment_mode_order <- function(fit, target_pool, prediction_mode) {
  if (identical(prediction_mode, "weights")) {
    return("weights")
  }
  if (identical(prediction_mode, "embedding")) {
    return("embedding")
  }

  if (identical(target_pool, "analysis") &&
      inherits(fit, c("fpgw", "gromov_wasserstein"))) {
    return(c("weights", "embedding"))
  }

  if (identical(target_pool, "analysis")) {
    return(c("embedding", "weights"))
  }

  "embedding"
}

.cv_alignment_row_sets <- function(feature_blocks, held_out) {
  out <- vector("list", length(feature_blocks))
  names(out) <- names(feature_blocks)

  for (block_name in names(feature_blocks)) {
    n <- nrow(feature_blocks[[block_name]])
    idx <- held_out[[block_name]]
    if (is.null(idx)) {
      idx <- integer(0)
    } else {
      idx <- as.integer(idx)
      idx <- idx[!duplicated(idx)]
    }
    out[[block_name]] <- list(
      analysis = setdiff(seq_len(n), idx),
      assessment = idx
    )
  }

  out
}

.cv_alignment_block_matrix <- function(resolved, block_name) {
  idx <- match(block_name, resolved$domain_names)
  if (is.na(idx)) {
    stop("Block `", block_name, "` is not present in the requested split.", call. = FALSE)
  }
  as.matrix(resolved$domains[[idx]]$x)
}

.cv_alignment_extract_analysis_scores <- function(fit, analysis) {
  scores <- fit$s
  if (is.null(scores)) {
    stop("Fit object does not contain training scores in `s`.", call. = FALSE)
  }
  scores <- as.matrix(scores)

  block_names <- analysis$domain_names
  block_sizes <- vapply(analysis$domains, function(dom) nrow(dom$x), integer(1))
  total_rows <- sum(block_sizes)

  if (nrow(scores) != total_rows) {
    stop(
      "Fit stores ", nrow(scores), " training score rows but the analysis split has ",
      total_rows, ".",
      call. = FALSE
    )
  }

  starts <- cumsum(c(1L, head(block_sizes, -1L)))
  out <- vector("list", length(block_names))
  names(out) <- block_names

  for (i in seq_along(block_names)) {
    rows <- seq.int(starts[[i]], length.out = block_sizes[[i]])
    out[[block_names[[i]]]] <- scores[rows, , drop = FALSE]
  }

  out
}

.cv_alignment_get_assessment_projection <- function(fit,
                                                    assessment,
                                                    block_name,
                                                    block_names,
                                                    cache,
                                                    fold_id = NULL) {
  if (exists(block_name, envir = cache, inherits = FALSE)) {
    return(get(block_name, envir = cache, inherits = FALSE))
  }

  Xi <- .cv_alignment_block_matrix(assessment, block_name)
  side_idx <- match(block_name, block_names)
  if (is.na(side_idx)) {
    stop("Unknown block `", block_name, "`.", call. = FALSE)
  }

  proj <- tryCatch(
    as.matrix(oos_predict(fit, Xi, side = block_name)),
    error = function(e_name) {
      tryCatch(
        as.matrix(oos_predict(fit, Xi, side = side_idx)),
        error = function(e_idx) {
          stop(
            "Failed to project held-out rows for block `", block_name, "` in fold ",
            if (is.null(fold_id)) "<unknown>" else fold_id,
            ". Ensure `fit_fn` returns an object with oos_predict() support. Original errors: ",
            e_name$message, " / ", e_idx$message,
            call. = FALSE
          )
        }
      )
    }
  )

  assign(block_name, proj, envir = cache)
  proj
}

.cv_alignment_target_pool_embeddings <- function(fit,
                                                 analysis,
                                                 assessment,
                                                 feature_blocks,
                                                 row_sets,
                                                 target_block,
                                                 target_pool,
                                                 analysis_scores,
                                                 projection_cache,
                                                 block_names,
                                                 fold_id = NULL) {
  score_parts <- list()
  feature_parts <- list()

  if (target_pool %in% c("analysis", "both")) {
    target_scores <- analysis_scores[[target_block]]
    target_features <- feature_blocks[[target_block]][row_sets[[target_block]]$analysis, , drop = FALSE]

    if (nrow(target_scores) != nrow(target_features)) {
      stop(
        "Training score rows for block `", target_block, "` do not match the analysis target feature rows.",
        call. = FALSE
      )
    }

    if (nrow(target_scores)) {
      score_parts[[length(score_parts) + 1L]] <- target_scores
      feature_parts[[length(feature_parts) + 1L]] <- target_features
    }
  }

  if (target_pool %in% c("assessment", "both") &&
      target_block %in% assessment$domain_names &&
      length(row_sets[[target_block]]$assessment)) {
    target_scores <- .cv_alignment_get_assessment_projection(
      fit = fit,
      assessment = assessment,
      block_name = target_block,
      block_names = block_names,
      cache = projection_cache,
      fold_id = fold_id
    )
    target_features <- feature_blocks[[target_block]][row_sets[[target_block]]$assessment, , drop = FALSE]

    if (nrow(target_scores) != nrow(target_features)) {
      stop(
        "Held-out target rows for block `", target_block, "` do not match projected assessment rows.",
        call. = FALSE
      )
    }

    if (nrow(target_scores)) {
      score_parts[[length(score_parts) + 1L]] <- target_scores
      feature_parts[[length(feature_parts) + 1L]] <- target_features
    }
  }

  if (!length(score_parts)) {
    stop(
      "No target rows were available for block `", target_block,
      "` under `target_pool = \"", target_pool, "\"`.",
      call. = FALSE
    )
  }

  list(
    scores = do.call(rbind, score_parts),
    features = do.call(rbind, feature_parts)
  )
}

.cv_alignment_predict_target_weights <- function(fit,
                                                 newdata,
                                                 from,
                                                 to,
                                                 block_names,
                                                 fold_id = NULL) {
  from_idx <- match(from, block_names)
  to_idx <- match(to, block_names)
  if (is.na(from_idx) || is.na(to_idx)) {
    stop("Unknown block name in weight prediction.", call. = FALSE)
  }

  pred <- tryCatch(
    stats::predict(fit, newdata = newdata, from = from, to = to, type = "weights"),
    error = function(e_name) {
      tryCatch(
        stats::predict(fit, newdata = newdata, from = from_idx, to = to_idx, type = "weights"),
        error = function(e_idx) {
          stop(
            "Failed to obtain weight predictions for block `", from, "` -> `", to,
            "` in fold ", if (is.null(fold_id)) "<unknown>" else fold_id,
            ". Original errors: ", e_name$message, " / ", e_idx$message,
            call. = FALSE
          )
        }
      )
    }
  )

  pred <- as.matrix(pred)
  if (nrow(pred) != nrow(newdata)) {
    stop(
      "Weight prediction for block `", from, "` -> `", to,
      "` returned ", nrow(pred), " rows for ", nrow(newdata), " queries.",
      call. = FALSE
    )
  }

  pred
}

.cv_alignment_topk_weight_indices <- function(weights, k) {
  weights <- as.matrix(weights)

  if (k == 1L) {
    return(matrix(max.col(weights, ties.method = "first"), ncol = 1L))
  }

  t(vapply(seq_len(nrow(weights)), function(i) {
    utils::head(order(weights[i, ], decreasing = TRUE), k)
  }, integer(k)))
}

.cv_alignment_knn_indices <- function(query, target, k) {
  query <- as.matrix(query)
  target <- as.matrix(target)

  if (requireNamespace("RANN", quietly = TRUE)) {
    return(RANN::nn2(data = target, query = query, k = k)$nn.idx)
  }

  s1 <- rowSums(query * query)
  s2 <- rowSums(target * target)
  D <- outer(s1, s2, "+") - 2 * (query %*% t(target))

  if (k == 1L) {
    return(matrix(max.col(-D, ties.method = "first"), ncol = 1L))
  }

  t(vapply(seq_len(nrow(D)), function(i) {
    utils::head(order(D[i, ], decreasing = FALSE), k)
  }, integer(k)))
}

.cv_alignment_similarity_matrix <- function(x, y, method = c("cosine", "correlation")) {
  method <- match.arg(method)
  x <- as.matrix(x)
  y <- as.matrix(y)

  if (ncol(x) != ncol(y)) {
    stop("Feature blocks must share the same number of columns for similarity scoring.", call. = FALSE)
  }

  if (method == "correlation") {
    x <- sweep(x, 1L, rowMeans(x), FUN = "-")
    y <- sweep(y, 1L, rowMeans(y), FUN = "-")
  }

  x_norm <- sqrt(rowSums(x * x))
  y_norm <- sqrt(rowSums(y * y))
  x_norm[x_norm <= 0] <- 1
  y_norm[y_norm <= 0] <- 1

  (x / x_norm) %*% t(y / y_norm)
}
