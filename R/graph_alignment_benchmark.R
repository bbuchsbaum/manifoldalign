#' Benchmark graph alignment methods on synthetic permuted graphs
#'
#' Generates paired domains where the target is a permuted/noisy copy of the
#' source and compares multiple alignment methods available in this package.
#'
#' Benchmarks are intended for interactive profiling and method comparison
#' rather than automated CI.
#'
#' @param sizes Integer vector of node counts to benchmark.
#' @param d Integer feature dimension for synthetic node features.
#' @param noise_sd Numeric noise standard deviation added to target features.
#' @param structure Character vector of synthetic geometries for base node
#'   coordinates: `"ring"`, `"grid"`, `"random"`, or `"community"`.
#' @param permute_fraction Fraction of nodes to permute (in `(0,1]`).
#' @param n_anchors Integer number of anchor correspondences to add (0 allowed).
#' @param methods Character vector of methods to run. Options include:
#'   `"token_ot_graph"`, `"cone_align"`, `"grasp"`, `"parrot"`, `"fpgw"`,
#'   `"ssma"`, `"lra"`, `"gpca"`.
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
#' @return A data.frame with per-run runtime and accuracy metrics.
#' @export
benchmark_graph_alignment_methods <- function(
  sizes = c(50L, 100L),
  d = 3L,
  noise_sd = 0.05,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 1,
  n_anchors = 0L,
  methods = c("token_ot_graph", "fpgw", "cone_align", "grasp", "parrot", "ssma", "lra", "gpca"),
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
  structure <- match.arg(structure, c("ring", "grid", "random", "community"), several.ok = TRUE)
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

  allowed_methods <- c("token_ot_graph", "fpgw", "cone_align", "grasp", "parrot", "ssma", "lra", "gpca")
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

  # Helper to evaluate a predicted mapping
  eval_top1 <- function(pred, truth) {
    pred <- as.integer(pred)
    truth <- as.integer(truth)
    ok <- !is.na(pred) & pred >= 1L & pred <= length(truth)
    if (!any(ok)) return(list(top1 = NA_real_, coverage = 0))
    list(
      top1 = mean(pred[ok] == truth[ok]),
      coverage = mean(ok)
    )
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
  decode_from_scores_mode <- function(scores, n, mode, anchor_idx, map12) {
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
        Q <- solve_procrustes_rotation(E1, E2, anchor_idx, map12[anchor_idx])
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

  make_result_row <- function(method, spec, runtime_sec, top1_accuracy = NA_real_, coverage = 0, error = NA_character_) {
    data.frame(
      method = method,
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
      error = error,
      stringsAsFactors = FALSE
    )
  }

  for (row_id in seq_len(nrow(specs))) {
    spec <- specs[row_id, , drop = FALSE]
    n <- spec$n[[1L]]
    r <- spec$rep[[1L]]

    case <- synthetic_graph_alignment_case(
      n = n,
      d = spec$d[[1L]],
      structure = spec$structure[[1L]],
      noise_sd = spec$noise_sd[[1L]],
      permute_fraction = spec$permute_fraction[[1L]],
      n_anchors = spec$n_anchors[[1L]],
      seed = spec$generation_seed[[1L]]
    )

    X1 <- case$X1
    X2 <- case$X2
    map12 <- case$map12
    anchor_idx <- case$anchor_idx
    hd <- case$hd
    anchors <- case$anchors_name

    # Embedding dimensions for graph-based aligners should depend on graph size,
    # not the raw feature dimension used to generate coordinates.
    ncomp_align <- min(10L, as.integer(n) - 2L)

    if (isTRUE(verbose)) {
      message("benchmark_graph_alignment_methods: structure=", spec$structure[[1L]], " n=", n, " rep=", r)
    }

      # Token-OT Graph Align
      if ("token_ot_graph" %in% methods) {
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
          token_ot_graph_align(hd, anchors = anchors, ncomp = ncomp_align, control = ctrl),
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
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("token_ot_graph_align", spec, runtime, ev$top1, ev$coverage)
        }
      }

      # FPGW (fused partial Gromov-Wasserstein)
      if ("fpgw" %in% methods) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          suppressMessages(
            fpgw(
              hd,
              omega1 = fpgw_omega1,
              epsilon = fpgw_epsilon,
              max_iter = as.integer(fpgw_max_iter),
              inner_max_iter = as.integer(fpgw_inner_max_iter),
              tol = fpgw_tol,
              verbose = FALSE
            )
          ),
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
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("fpgw", spec, runtime, coverage = 0, error = "Missing transport plan from fpgw() result")
          } else {
            pred_native <- top1_from_plan(plan)
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ev <- eval_top1(pred, map12)
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("fpgw", spec, runtime, ev$top1, ev$coverage)
          }
        }
      }

      # SSMA (semi-supervised manifold alignment baseline)
      if ("ssma" %in% methods) {
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
          suppressMessages(
            ssma_align(
              hd,
              correspondences = ssma_corr,
              serial_index = ssma_serial_index,
              ncomp = ncomp_align,
              control = ssma_ctrl
            )
          ),
          error = function(e) e
        )
        runtime <- proc.time()[[3]] - t0

        if (inherits(res, "error")) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("ssma_align", spec, runtime, coverage = 0, error = as.character(res$message))
        } else {
          pred_native <- decode_from_scores_mode(res$s, n = n, mode = "common_nn", anchor_idx = anchor_idx, map12 = map12)
          if (is.null(pred_native)) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("ssma_align", spec, runtime, coverage = 0, error = "Missing or malformed score matrix from ssma_align()")
          } else {
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              ssma_mode <- decode_mode
              if (identical(ssma_mode, "common_procrustes_lap") && !isTRUE(ssma_procrustes)) {
                ssma_mode <- "common_nn"
              }
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = ssma_mode, anchor_idx = anchor_idx, map12 = map12)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ev <- eval_top1(pred, map12)
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("ssma_align", spec, runtime, ev$top1, ev$coverage)
          }
        }
      }

      # LRA (low-rank alignment; semi-supervised labels)
      if ("lra" %in% methods) {
        if (as.integer(n_anchors) <= 0L) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("lowrank_align", spec, NA_real_, coverage = 0, error = "LRA benchmark requires n_anchors > 0")
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            suppressWarnings(
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
            ),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("lowrank_align", spec, runtime, coverage = 0, error = as.character(res$message))
          } else {
            scores <- if (!is.null(res$s)) as.matrix(res$s) else NULL
            if (is.null(scores) || nrow(scores) < 2L * n) {
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row("lowrank_align", spec, runtime, coverage = 0, error = "Missing or malformed score matrix from lowrank_align()")
            } else {
              pred_native <- decode_from_scores_mode(scores, n = n, mode = "common_nn", anchor_idx = anchor_idx, map12 = map12)
              pred <- pred_native
              if (!identical(decode_mode, "native")) {
                pred_common <- decode_from_scores_mode(scores, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
                if (!is.null(pred_common)) pred <- pred_common
              }
              ev <- eval_top1(pred, map12)
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row("lowrank_align", spec, runtime, ev$top1, ev$coverage)
            }
          }
        }
      }

      # GPCA alignment (semi-supervised labels)
      if ("gpca" %in% methods) {
        if (as.integer(n_anchors) <= 0L) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("gpca_align", spec, NA_real_, coverage = 0, error = "GPCA benchmark requires n_anchors > 0")
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            suppressWarnings(
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
            ),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("gpca_align", spec, runtime, coverage = 0, error = as.character(res$message))
          } else {
            scores <- if (!is.null(res$s)) as.matrix(res$s) else NULL
            if (is.null(scores) || nrow(scores) < 2L * n) {
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row("gpca_align", spec, runtime, coverage = 0, error = "Missing or malformed score matrix from gpca_align()")
            } else {
              pred_native <- decode_from_scores_mode(scores, n = n, mode = "common_nn", anchor_idx = anchor_idx, map12 = map12)
              pred <- pred_native
              if (!identical(decode_mode, "native")) {
                pred_common <- decode_from_scores_mode(scores, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
                if (!is.null(pred_common)) pred <- pred_common
              }
              ev <- eval_top1(pred, map12)
              ptr <- ptr + 1L
              rows[[ptr]] <- make_result_row("gpca_align", spec, runtime, ev$top1, ev$coverage)
            }
          }
        }
      }

      # CONE-Align
      if ("cone_align" %in% methods) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          suppressMessages(cone_align(hd, anchors = anchors, ncomp = ncomp_align)),
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
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("cone_align", spec, runtime, ev$top1, ev$coverage)
        }
      }

      # GRASP
      if ("grasp" %in% methods) {
        t0 <- proc.time()[[3]]
        res <- tryCatch(
          suppressMessages(grasp(hd, ncomp = ncomp_align, q_descriptors = 30)),
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
            pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
            if (!is.null(pred_common)) pred <- pred_common
          }
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("grasp", spec, runtime, ev$top1, ev$coverage)
        }
      }

      # PARROT (requires anchors)
      if ("parrot" %in% methods) {
        if (as.integer(n_anchors) <= 0L) {
          ptr <- ptr + 1L
          rows[[ptr]] <- make_result_row("parrot", spec, NA_real_, coverage = 0, error = "PARROT requires n_anchors > 0")
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            suppressMessages(parrot(hd, anchors = anchors, ncomp = ncomp_align, max_iter = 30)),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("parrot", spec, runtime, coverage = 0, error = as.character(res$message))
          } else {
            pred_native <- top1_from_plan(res$transport_plan)
            pred <- pred_native
            if (!identical(decode_mode, "native")) {
              pred_common <- decode_from_scores_mode(res$s, n = n, mode = decode_mode, anchor_idx = anchor_idx, map12 = map12)
              if (!is.null(pred_common)) pred <- pred_common
            }
            ev <- eval_top1(pred, map12)
            ptr <- ptr + 1L
            rows[[ptr]] <- make_result_row("parrot", spec, runtime, ev$top1, ev$coverage)
          }
        }
      }
  }

  do.call(rbind, rows[seq_len(ptr)])
}
