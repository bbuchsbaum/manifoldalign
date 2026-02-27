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
#' @param structure Synthetic geometry for base node coordinates:
#'   `"ring"`, `"grid"`, `"random"`, or `"community"`.
#' @param permute_fraction Fraction of nodes to permute (in `(0,1]`).
#' @param n_anchors Integer number of anchor correspondences to add (0 allowed).
#' @param methods Character vector of methods to run. Options include:
#'   `"token_ot_graph"`, `"cone_align"`, `"grasp"`, `"parrot"`, `"fpgw"`.
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
  methods = c("token_ot_graph", "fpgw", "cone_align", "grasp", "parrot"),
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
  n_reps = 3L,
  seed = 1L,
  verbose = FALSE
) {
  structure <- match.arg(structure)
  coarsen_method <- match.arg(coarsen_method)
  token_mode <- match.arg(token_mode)
  prior_mode <- match.arg(prior_mode)

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
  chk::chk_number(n_reps)
  chk::chk_true(n_reps >= 1)
  chk::chk_number(seed)
  chk::chk_logical(verbose)

  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) stop("`sizes` must contain integers >= 3.", call. = FALSE)

  allowed_methods <- c("token_ot_graph", "fpgw", "cone_align", "grasp", "parrot")
  bad <- setdiff(methods, allowed_methods)
  if (length(bad)) stop("Unknown method(s): ", paste(bad, collapse = ", "), call. = FALSE)

  # Helper to generate base coordinates
  make_base <- function(n, d, structure) {
    if (structure == "ring") {
      angles <- 2 * pi * (seq_len(n)) / n
      X <- cbind(cos(angles), sin(angles))
      if (d > 2L) {
        X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
      } else if (d == 1L) {
        X <- matrix(cos(angles), n, 1L)
      }
      return(X)
    }
    if (structure == "grid") {
      grid_size <- ceiling(sqrt(n))
      grid_coords <- expand.grid(
        x = seq_len(grid_size) / grid_size,
        y = seq_len(grid_size) / grid_size
      )[seq_len(n), ]
      X <- cbind(grid_coords$x, grid_coords$y)
      if (d > 2L) {
        X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
      } else if (d == 1L) {
        X <- matrix(grid_coords$x, n, 1L)
      }
      return(as.matrix(X))
    }
    if (structure == "community") {
      k <- 4L
      X <- matrix(0, n, max(2L, d))
      n_per <- ceiling(n / k)
      for (g in seq_len(k)) {
        idx <- seq.int((g - 1L) * n_per + 1L, min(g * n_per, n))
        center <- c(cos(2 * pi * g / k), sin(2 * pi * g / k))
        X[idx, 1:2] <- matrix(stats::rnorm(length(idx) * 2, sd = 0.2), length(idx), 2) +
          matrix(center, length(idx), 2, byrow = TRUE)
      }
      if (d > 2L) {
        X[, 3:d] <- matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L)
      }
      return(X[, seq_len(d), drop = FALSE])
    }
    # random
    matrix(stats::rnorm(n * d), n, d)
  }

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

  rows <- list()
  ptr <- 0L

  for (s in seq_along(sizes)) {
    n <- sizes[[s]]
    for (r in seq_len(as.integer(n_reps))) {
      set.seed(as.integer(seed) + s * 1009L + r * 17L)

      X1 <- make_base(n, as.integer(d), structure = structure)

      # Permute a subset of nodes; default permutes all nodes.
      n_permute <- floor(n * as.numeric(permute_fraction))
      perm_idx <- if (n_permute >= 2L) sort(sample.int(n, n_permute)) else integer(0)
      perm21 <- seq_len(n)
      if (length(perm_idx)) {
        perm21[perm_idx] <- perm_idx[sample.int(length(perm_idx))]
      }

      X2 <- X1[perm21, , drop = FALSE] + matrix(stats::rnorm(n * as.integer(d), sd = noise_sd), n, as.integer(d))
      map12 <- match(seq_len(n), perm21)

      # Embedding dimensions for graph-based aligners should depend on graph size,
      # not the raw feature dimension used to generate coordinates.
      ncomp_align <- min(10L, as.integer(n) - 2L)

      # Anchors (optional)
      a1 <- rep(NA_integer_, n)
      a2 <- rep(NA_integer_, n)
      if (as.integer(n_anchors) > 0L) {
        # Prefer anchors outside permuted set (if possible)
        available <- setdiff(seq_len(n), perm_idx)
        if (length(available) < as.integer(n_anchors)) available <- seq_len(n)
        anchor_idx <- sort(sample(available, min(as.integer(n_anchors), length(available))))
        a1[anchor_idx] <- seq_along(anchor_idx)
        a2[map12[anchor_idx]] <- seq_along(anchor_idx)
      }

      hd <- as_hyperdesign(
        X_list = list(domain1 = X1, domain2 = X2),
        labels = list(a1, a2),
        label_name = "anchors"
      )

      if (isTRUE(verbose)) {
        message("benchmark_graph_alignment_methods: n=", n, " rep=", r)
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
          rows[[ptr]] <- data.frame(
            method = "token_ot_graph_align",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = NA_real_,
            coverage = 0,
            error = as.character(res$message),
            stringsAsFactors = FALSE
          )
        } else {
          pred <- as.integer(res$assignment)
          if (all(is.na(pred)) && !is.null(res$transport_plan)) {
            pred <- top1_from_plan(res$transport_plan)
          }
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- data.frame(
            method = "token_ot_graph_align",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = ev$top1,
            coverage = ev$coverage,
            error = NA_character_,
            stringsAsFactors = FALSE
          )
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
          rows[[ptr]] <- data.frame(
            method = "fpgw",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = NA_real_,
            coverage = 0,
            error = as.character(res$message),
            stringsAsFactors = FALSE
          )
        } else {
          plan <- NULL
          if (!is.null(res$transport_plans) && length(res$transport_plans) >= 1L) {
            plan <- res$transport_plans[[1L]]
          }
          if (is.null(plan)) {
            ptr <- ptr + 1L
            rows[[ptr]] <- data.frame(
              method = "fpgw",
              n = n,
              rep = r,
              runtime_sec = runtime,
              top1_accuracy = NA_real_,
              coverage = 0,
              error = "Missing transport plan from fpgw() result",
              stringsAsFactors = FALSE
            )
          } else {
            pred <- top1_from_plan(plan)
            ev <- eval_top1(pred, map12)
            ptr <- ptr + 1L
            rows[[ptr]] <- data.frame(
              method = "fpgw",
              n = n,
              rep = r,
              runtime_sec = runtime,
              top1_accuracy = ev$top1,
              coverage = ev$coverage,
              error = NA_character_,
              stringsAsFactors = FALSE
            )
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
          rows[[ptr]] <- data.frame(
            method = "cone_align",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = NA_real_,
            coverage = 0,
            error = as.character(res$message),
            stringsAsFactors = FALSE
          )
        } else {
          pred <- as.integer(res$assignment)
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- data.frame(
            method = "cone_align",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = ev$top1,
            coverage = ev$coverage,
            error = NA_character_,
            stringsAsFactors = FALSE
          )
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
          rows[[ptr]] <- data.frame(
            method = "grasp",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = NA_real_,
            coverage = 0,
            error = as.character(res$message),
            stringsAsFactors = FALSE
          )
        } else {
          pred <- as.integer(res$assignment)
          ev <- eval_top1(pred, map12)
          ptr <- ptr + 1L
          rows[[ptr]] <- data.frame(
            method = "grasp",
            n = n,
            rep = r,
            runtime_sec = runtime,
            top1_accuracy = ev$top1,
            coverage = ev$coverage,
            error = NA_character_,
            stringsAsFactors = FALSE
          )
        }
      }

      # PARROT (requires anchors)
      if ("parrot" %in% methods) {
        if (as.integer(n_anchors) <= 0L) {
          ptr <- ptr + 1L
          rows[[ptr]] <- data.frame(
            method = "parrot",
            n = n,
            rep = r,
            runtime_sec = NA_real_,
            top1_accuracy = NA_real_,
            coverage = 0,
            error = "PARROT requires n_anchors > 0",
            stringsAsFactors = FALSE
          )
        } else {
          t0 <- proc.time()[[3]]
          res <- tryCatch(
            suppressMessages(parrot(hd, anchors = anchors, ncomp = ncomp_align, max_iter = 30)),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          if (inherits(res, "error")) {
            ptr <- ptr + 1L
            rows[[ptr]] <- data.frame(
              method = "parrot",
              n = n,
              rep = r,
              runtime_sec = runtime,
              top1_accuracy = NA_real_,
              coverage = 0,
              error = as.character(res$message),
              stringsAsFactors = FALSE
            )
          } else {
            pred <- top1_from_plan(res$transport_plan)
            ev <- eval_top1(pred, map12)
            ptr <- ptr + 1L
            rows[[ptr]] <- data.frame(
              method = "parrot",
              n = n,
              rep = r,
              runtime_sec = runtime,
              top1_accuracy = ev$top1,
              coverage = ev$coverage,
              error = NA_character_,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  do.call(rbind, rows[seq_len(ptr)])
}
