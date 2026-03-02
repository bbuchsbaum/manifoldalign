#' Assess SSMA Alignment Quality on Synthetic Paired Domains
#'
#' Runs `ssma_align()` on synthetic paired-domain problems and reports diagnostics
#' that separate embedding quality from decoding quality.
#'
#' @param sizes Integer vector of node counts to assess.
#' @param d Integer feature dimension for synthetic node features.
#' @param noise_sd Numeric noise standard deviation added to target features.
#' @param structure Synthetic geometry for base coordinates:
#'   `"ring"`, `"grid"`, `"random"`, or `"community"`.
#' @param permute_fraction Fraction of nodes to permute (in `(0,1]`).
#' @param n_anchors Number of synthetic anchor pairs to sample.
#' @param holdout_fraction Fraction of anchors reserved for holdout evaluation.
#' @param n_reps Integer number of replications per size.
#' @param solvers Character vector of SSMA solvers to assess:
#'   `"reduced"` and/or `"operator"`.
#' @param serial Logical vector indicating whether to run serial decontamination
#'   variants (`FALSE`, `TRUE`, or both).
#' @param ssma_knn Integer k used for SSMA within-domain graph construction.
#' @param ssma_rank_per_domain Integer rank cap for SSMA per-domain reduction.
#' @param ssma_lag_exclusion Integer lag exclusion for hard serial masking when
#'   `serial = TRUE`.
#' @param ncomp Optional integer embedding dimension. If `NULL`, uses
#'   `min(10, n - 2)`.
#' @param decode_k Integer `k` used for top-k retrieval metrics.
#' @param seed Integer base random seed.
#' @param verbose Logical; print per-run progress.
#'
#' @return A data.frame with per-run SSMA diagnostics including:
#'   runtime, embedding finite-ness, raw/procrustes decoding metrics, holdout
#'   anchor recovery, eigen diagnostics, serial diagnostics, and cross-solver
#'   subspace distance when both solvers are run.
#' @export
assess_ssma <- function(
  sizes = c(50L, 100L),
  d = 3L,
  noise_sd = 0.05,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 1,
  n_anchors = 8L,
  holdout_fraction = 0.25,
  n_reps = 3L,
  solvers = c("reduced", "operator"),
  serial = c(FALSE, TRUE),
  ssma_knn = 10L,
  ssma_rank_per_domain = 64L,
  ssma_lag_exclusion = 1L,
  ncomp = NULL,
  decode_k = 5L,
  seed = 1L,
  verbose = FALSE
) {
  structure <- match.arg(structure)

  chk::chk_number(d)
  chk::chk_true(d >= 1)
  chk::chk_number(noise_sd)
  chk::chk_true(noise_sd >= 0)
  chk::chk_number(permute_fraction)
  chk::chk_true(is.finite(permute_fraction) && permute_fraction > 0 && permute_fraction <= 1)
  chk::chk_number(n_anchors)
  chk::chk_true(n_anchors >= 0)
  chk::chk_number(holdout_fraction)
  chk::chk_true(is.finite(holdout_fraction) && holdout_fraction >= 0 && holdout_fraction < 1)
  chk::chk_number(n_reps)
  chk::chk_true(n_reps >= 1)
  chk::chk_number(ssma_knn)
  chk::chk_true(ssma_knn >= 1)
  chk::chk_number(ssma_rank_per_domain)
  chk::chk_true(ssma_rank_per_domain >= 1)
  chk::chk_number(ssma_lag_exclusion)
  chk::chk_true(ssma_lag_exclusion >= 0)
  if (!is.null(ncomp)) {
    chk::chk_number(ncomp)
    chk::chk_true(ncomp >= 1)
  }
  chk::chk_number(decode_k)
  chk::chk_true(decode_k >= 1)
  chk::chk_number(seed)
  chk::chk_logical(verbose)

  sizes <- as.integer(sizes)
  sizes <- sizes[is.finite(sizes) & sizes >= 3L]
  if (!length(sizes)) stop("`sizes` must contain integers >= 3.", call. = FALSE)

  if (!is.character(solvers) || !length(solvers)) {
    stop("`solvers` must be a non-empty character vector.", call. = FALSE)
  }
  allowed_solvers <- c("reduced", "operator")
  bad_solvers <- setdiff(solvers, allowed_solvers)
  if (length(bad_solvers)) {
    stop("Unknown solver(s): ", paste(bad_solvers, collapse = ", "), call. = FALSE)
  }
  solvers <- unique(solvers)

  if (!is.logical(serial) || !length(serial) || any(is.na(serial))) {
    stop("`serial` must be a non-empty logical vector without NA.", call. = FALSE)
  }
  serial <- unique(as.logical(serial))

  make_base <- function(n, d, structure) {
    if (structure == "ring") {
      angles <- 2 * pi * seq_len(n) / n
      X <- cbind(cos(angles), sin(angles))
      if (d > 2L) {
        X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
      } else if (d == 1L) {
        X <- matrix(cos(angles), n, 1L)
      }
      return(X)
    }
    if (structure == "grid") {
      g <- ceiling(sqrt(n))
      coords <- expand.grid(
        x = seq_len(g) / g,
        y = seq_len(g) / g
      )[seq_len(n), ]
      X <- cbind(coords$x, coords$y)
      if (d > 2L) {
        X <- cbind(X, matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L))
      } else if (d == 1L) {
        X <- matrix(coords$x, n, 1L)
      }
      return(as.matrix(X))
    }
    if (structure == "community") {
      k <- 4L
      X <- matrix(0, n, max(2L, d))
      n_per <- ceiling(n / k)
      for (g in seq_len(k)) {
        idx <- seq.int((g - 1L) * n_per + 1L, min(g * n_per, n))
        ctr <- c(cos(2 * pi * g / k), sin(2 * pi * g / k))
        X[idx, 1:2] <- matrix(stats::rnorm(length(idx) * 2, sd = 0.2), length(idx), 2) +
          matrix(ctr, length(idx), 2, byrow = TRUE)
      }
      if (d > 2L) {
        X[, 3:d] <- matrix(stats::runif(n * (d - 2L), -0.1, 0.1), n, d - 2L)
      }
      return(X[, seq_len(d), drop = FALSE])
    }
    matrix(stats::rnorm(n * d), n, d)
  }

  top1_from_scores <- function(E1, E2) {
    E1 <- as.matrix(E1)
    E2 <- as.matrix(E2)
    if (!nrow(E1) || !nrow(E2)) return(integer(0))

    if (requireNamespace("RANN", quietly = TRUE)) {
      return(as.integer(RANN::nn2(data = E2, query = E1, k = 1L)$nn.idx[, 1L]))
    }

    s1 <- rowSums(E1 * E1)
    s2 <- rowSums(E2 * E2)
    D <- outer(s1, s2, "+") - 2 * (E1 %*% t(E2))
    max.col(-D, ties.method = "first")
  }

  recall_at_k <- function(E1, E2, truth, k) {
    truth <- as.integer(truth)
    ok <- is.finite(truth) & truth >= 1L & truth <= nrow(E2)
    if (!any(ok)) return(NA_real_)
    k_use <- min(as.integer(k), nrow(E2))

    if (requireNamespace("RANN", quietly = TRUE)) {
      nn <- RANN::nn2(data = E2, query = E1, k = k_use)$nn.idx
      if (k_use == 1L) {
        hit <- nn[, 1L] == truth
      } else {
        hit <- rowSums(nn == truth) > 0
      }
      return(mean(hit[ok]))
    }

    s1 <- rowSums(E1 * E1)
    s2 <- rowSums(E2 * E2)
    D <- outer(s1, s2, "+") - 2 * (E1 %*% t(E2))
    hit <- logical(nrow(D))
    for (i in seq_len(nrow(D))) {
      idx <- utils::head(order(D[i, ], decreasing = FALSE), k_use)
      hit[i] <- truth[i] %in% idx
    }
    mean(hit[ok])
  }

  solve_procrustes_rotation <- function(E1, E2, idx1, idx2) {
    idx1 <- as.integer(idx1)
    idx2 <- as.integer(idx2)
    keep <- is.finite(idx1) & is.finite(idx2) &
      idx1 >= 1L & idx1 <= nrow(E1) &
      idx2 >= 1L & idx2 <= nrow(E2)
    idx1 <- idx1[keep]
    idx2 <- idx2[keep]
    if (length(idx1) < 2L) return(NULL)

    A <- scale(E1[idx1, , drop = FALSE], center = TRUE, scale = FALSE)
    B <- scale(E2[idx2, , drop = FALSE], center = TRUE, scale = FALSE)

    C <- crossprod(B, A)
    sv <- tryCatch(svd(C), error = function(e) NULL)
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

  lap_assignment <- function(E1, E2) {
    if (!requireNamespace("clue", quietly = TRUE)) {
      return(rep(NA_integer_, nrow(E1)))
    }
    s1 <- rowSums(E1 * E1)
    s2 <- rowSums(E2 * E2)
    D <- outer(s1, s2, "+") - 2 * (E1 %*% t(E2))
    D <- pmax(D, 0)
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

  subspace_distance <- function(S1, S2) {
    S1 <- as.matrix(S1)
    S2 <- as.matrix(S2)
    if (!ncol(S1) || !ncol(S2)) return(NA_real_)
    Q1 <- qr.Q(qr(S1))
    Q2 <- qr.Q(qr(S2))
    k <- min(ncol(Q1), ncol(Q2))
    if (k < 1L) return(NA_real_)
    sv <- svd(crossprod(Q1[, seq_len(k), drop = FALSE], Q2[, seq_len(k), drop = FALSE]), nu = 0, nv = 0)$d
    sqrt(sum(1 - pmin(sv, 1)^2))
  }

  evaluate_top1 <- function(pred, truth, idx = NULL) {
    pred <- as.integer(pred)
    truth <- as.integer(truth)
    if (!is.null(idx)) {
      pred <- pred[idx]
      truth <- truth[idx]
    }
    ok <- !is.na(pred) & pred >= 1L & pred <= length(truth)
    if (!any(ok)) return(NA_real_)
    mean(pred[ok] == truth[ok])
  }

  rows <- list()
  ptr <- 0L
  score_cache <- list()

  for (s in seq_along(sizes)) {
    n <- sizes[[s]]
    ncomp_use <- if (is.null(ncomp)) min(10L, as.integer(n) - 2L) else as.integer(ncomp)

    for (r in seq_len(as.integer(n_reps))) {
      set.seed(as.integer(seed) + s * 1009L + r * 17L)

      X1 <- make_base(n, as.integer(d), structure)

      n_permute <- floor(n * as.numeric(permute_fraction))
      perm_idx <- if (n_permute >= 2L) sort(sample.int(n, n_permute)) else integer(0)
      perm21 <- seq_len(n)
      if (length(perm_idx)) {
        perm21[perm_idx] <- perm_idx[sample.int(length(perm_idx))]
      }
      X2 <- X1[perm21, , drop = FALSE] +
        matrix(stats::rnorm(n * as.integer(d), sd = noise_sd), n, as.integer(d))
      map12 <- match(seq_len(n), perm21)

      # Anchor split for train/holdout diagnostics.
      train_idx <- integer(0)
      holdout_idx <- integer(0)
      if (as.integer(n_anchors) > 0L) {
        available <- setdiff(seq_len(n), perm_idx)
        if (length(available) < as.integer(n_anchors)) available <- seq_len(n)
        all_anchor_idx <- sort(sample(available, min(as.integer(n_anchors), length(available))))

        n_hold <- floor(length(all_anchor_idx) * holdout_fraction)
        if (length(all_anchor_idx) >= 2L && holdout_fraction > 0 && n_hold < 1L) n_hold <- 1L
        if (n_hold >= length(all_anchor_idx)) n_hold <- max(0L, length(all_anchor_idx) - 1L)

        if (n_hold > 0L) {
          holdout_idx <- sort(sample(all_anchor_idx, n_hold))
        }
        train_idx <- setdiff(all_anchor_idx, holdout_idx)
      }

      train_corr <- if (length(train_idx)) {
        data.frame(
          domain_i = rep(1L, length(train_idx)),
          index_i = train_idx,
          domain_j = rep(2L, length(train_idx)),
          index_j = map12[train_idx],
          weight = rep(1, length(train_idx)),
          source = rep("anchor_train", length(train_idx)),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }

      hd <- as_hyperdesign(list(domain1 = X1, domain2 = X2))

      for (serial_mode in serial) {
        serial_index <- if (isTRUE(serial_mode)) list(seq_len(n), seq_len(n)) else NULL

        for (solver in solvers) {
          if (isTRUE(verbose)) {
            message("assess_ssma: n=", n, " rep=", r, " solver=", solver, " serial=", serial_mode)
          }

          ctrl <- ssma_align_control(
            knn = as.integer(ssma_knn),
            rank_per_domain = as.integer(ssma_rank_per_domain),
            solver = solver,
            verbose = FALSE,
            serial_control = ssma_serial_control(
              enabled = isTRUE(serial_mode),
              row_whiten = if (isTRUE(serial_mode)) "ar1" else "none",
              lag_mode = if (isTRUE(serial_mode)) "hard" else "none",
              lag_exclusion = as.integer(ssma_lag_exclusion),
              preserve_degree = TRUE
            )
          )

          t0 <- proc.time()[[3]]
          fit <- tryCatch(
            suppressMessages(
              ssma_align(
                hd,
                correspondences = train_corr,
                serial_index = serial_index,
                ncomp = ncomp_use,
                control = ctrl
              )
            ),
            error = function(e) e
          )
          runtime <- proc.time()[[3]] - t0

          ptr <- ptr + 1L

          if (inherits(fit, "error")) {
            rows[[ptr]] <- data.frame(
              method = "ssma_align",
              n = n,
              rep = r,
              solver = solver,
              serial_enabled = isTRUE(serial_mode),
              n_train_anchors = length(train_idx),
              n_holdout_anchors = length(holdout_idx),
              runtime_sec = runtime,
              is_multiblock = FALSE,
              scores_finite = FALSE,
              top1_raw_nn = NA_real_,
              topk_raw_nn = NA_real_,
              top1_procrustes_nn = NA_real_,
              topk_procrustes_nn = NA_real_,
              top1_procrustes_lap = NA_real_,
              top1_holdout_raw_nn = NA_real_,
              top1_holdout_procrustes_nn = NA_real_,
              top1_holdout_procrustes_lap = NA_real_,
              paired_mse_raw = NA_real_,
              paired_mse_procrustes = NA_real_,
              eig_n_positive = NA_integer_,
              eig_min_positive = NA_real_,
              eig_max = NA_real_,
              serial_masked_fraction_domain1 = NA_real_,
              serial_masked_fraction_domain2 = NA_real_,
              serial_ar1_phi_domain1 = NA_real_,
              serial_ar1_phi_domain2 = NA_real_,
              subspace_dist_to_reduced = NA_real_,
              error = as.character(fit$message),
              stringsAsFactors = FALSE
            )
            next
          }

          scores <- if (!is.null(fit$s)) as.matrix(fit$s) else NULL
          ok_scores <- !is.null(scores) && nrow(scores) >= 2L * n && all(is.finite(scores))

          if (!ok_scores) {
            rows[[ptr]] <- data.frame(
              method = "ssma_align",
              n = n,
              rep = r,
              solver = solver,
              serial_enabled = isTRUE(serial_mode),
              n_train_anchors = length(train_idx),
              n_holdout_anchors = length(holdout_idx),
              runtime_sec = runtime,
              is_multiblock = inherits(fit, "multiblock_biprojector"),
              scores_finite = FALSE,
              top1_raw_nn = NA_real_,
              topk_raw_nn = NA_real_,
              top1_procrustes_nn = NA_real_,
              topk_procrustes_nn = NA_real_,
              top1_procrustes_lap = NA_real_,
              top1_holdout_raw_nn = NA_real_,
              top1_holdout_procrustes_nn = NA_real_,
              top1_holdout_procrustes_lap = NA_real_,
              paired_mse_raw = NA_real_,
              paired_mse_procrustes = NA_real_,
              eig_n_positive = NA_integer_,
              eig_min_positive = NA_real_,
              eig_max = NA_real_,
              serial_masked_fraction_domain1 = NA_real_,
              serial_masked_fraction_domain2 = NA_real_,
              serial_ar1_phi_domain1 = NA_real_,
              serial_ar1_phi_domain2 = NA_real_,
              subspace_dist_to_reduced = NA_real_,
              error = "Missing or malformed score matrix from ssma_align()",
              stringsAsFactors = FALSE
            )
            next
          }

          E1 <- scores[seq_len(n), , drop = FALSE]
          E2 <- scores[n + seq_len(n), , drop = FALSE]

          pred_raw <- top1_from_scores(E1, E2)
          top1_raw <- evaluate_top1(pred_raw, map12)
          topk_raw <- recall_at_k(E1, E2, map12, k = as.integer(decode_k))
          mse_raw <- mean(rowSums((E1 - E2[map12, , drop = FALSE])^2))

          # Rotation from train anchors, else fallback to discovered correspondences.
          rot_idx1 <- integer(0)
          rot_idx2 <- integer(0)
          if (length(train_idx)) {
            rot_idx1 <- train_idx
            rot_idx2 <- map12[train_idx]
          } else if (is.data.frame(fit$correspondences) && nrow(fit$correspondences)) {
            cc <- fit$correspondences
            keep12 <- cc$domain_i == 1L & cc$domain_j == 2L
            keep21 <- cc$domain_i == 2L & cc$domain_j == 1L
            idx1 <- c(cc$index_i[keep12], cc$index_j[keep21])
            idx2 <- c(cc$index_j[keep12], cc$index_i[keep21])
            good <- is.finite(idx1) & is.finite(idx2) &
              idx1 >= 1L & idx1 <= n &
              idx2 >= 1L & idx2 <= n
            rot_idx1 <- as.integer(idx1[good])
            rot_idx2 <- as.integer(idx2[good])
          }

          Q <- solve_procrustes_rotation(E1, E2, rot_idx1, rot_idx2)

          top1_pnn <- NA_real_
          topk_pnn <- NA_real_
          top1_lap <- NA_real_
          mse_p <- NA_real_
          hold_raw <- NA_real_
          hold_pnn <- NA_real_
          hold_lap <- NA_real_

          if (!is.null(Q)) {
            E2p <- E2 %*% Q
            pred_pnn <- top1_from_scores(E1, E2p)
            pred_lap <- lap_assignment(E1, E2p)

            top1_pnn <- evaluate_top1(pred_pnn, map12)
            topk_pnn <- recall_at_k(E1, E2p, map12, k = as.integer(decode_k))
            top1_lap <- evaluate_top1(pred_lap, map12)
            mse_p <- mean(rowSums((E1 - E2p[map12, , drop = FALSE])^2))

            if (length(holdout_idx)) {
              hold_raw <- evaluate_top1(pred_raw, map12, idx = holdout_idx)
              hold_pnn <- evaluate_top1(pred_pnn, map12, idx = holdout_idx)
              hold_lap <- evaluate_top1(pred_lap, map12, idx = holdout_idx)
            }
          } else if (length(holdout_idx)) {
            hold_raw <- evaluate_top1(pred_raw, map12, idx = holdout_idx)
          }

          eigvals <- fit$generalized_eigenvalues
          eigvals <- if (is.null(eigvals)) numeric(0) else as.numeric(eigvals)
          eig_pos <- eigvals[is.finite(eigvals) & eigvals > 0]

          serial1 <- if (is.list(fit$serial) && length(fit$serial) >= 1L) fit$serial[[1L]] else NULL
          serial2 <- if (is.list(fit$serial) && length(fit$serial) >= 2L) fit$serial[[2L]] else NULL

          rows[[ptr]] <- data.frame(
            method = "ssma_align",
            n = n,
            rep = r,
            solver = solver,
            serial_enabled = isTRUE(serial_mode),
            n_train_anchors = length(train_idx),
            n_holdout_anchors = length(holdout_idx),
            runtime_sec = runtime,
            is_multiblock = inherits(fit, "multiblock_biprojector"),
            scores_finite = TRUE,
            top1_raw_nn = top1_raw,
            topk_raw_nn = topk_raw,
            top1_procrustes_nn = top1_pnn,
            topk_procrustes_nn = topk_pnn,
            top1_procrustes_lap = top1_lap,
            top1_holdout_raw_nn = hold_raw,
            top1_holdout_procrustes_nn = hold_pnn,
            top1_holdout_procrustes_lap = hold_lap,
            paired_mse_raw = mse_raw,
            paired_mse_procrustes = mse_p,
            eig_n_positive = length(eig_pos),
            eig_min_positive = if (length(eig_pos)) min(eig_pos) else NA_real_,
            eig_max = if (length(eigvals)) max(eigvals[is.finite(eigvals)]) else NA_real_,
            serial_masked_fraction_domain1 = if (!is.null(serial1$masked_fraction)) as.numeric(serial1$masked_fraction) else NA_real_,
            serial_masked_fraction_domain2 = if (!is.null(serial2$masked_fraction)) as.numeric(serial2$masked_fraction) else NA_real_,
            serial_ar1_phi_domain1 = if (!is.null(serial1$ar1_phi)) as.numeric(serial1$ar1_phi) else NA_real_,
            serial_ar1_phi_domain2 = if (!is.null(serial2$ar1_phi)) as.numeric(serial2$ar1_phi) else NA_real_,
            subspace_dist_to_reduced = NA_real_,
            error = NA_character_,
            stringsAsFactors = FALSE
          )

          key <- paste(n, r, as.integer(serial_mode), sep = "_")
          if (is.null(score_cache[[key]])) score_cache[[key]] <- list()
          score_cache[[key]][[solver]] <- list(row = ptr, scores = scores)
        }
      }
    }
  }

  if (ptr == 0L) {
    return(data.frame())
  }

  out <- do.call(rbind, rows[seq_len(ptr)])

  for (key in names(score_cache)) {
    item <- score_cache[[key]]
    if (!all(c("reduced", "operator") %in% names(item))) next
    dist <- subspace_distance(item$reduced$scores, item$operator$scores)
    out$subspace_dist_to_reduced[item$reduced$row] <- dist
    out$subspace_dist_to_reduced[item$operator$row] <- dist
  }

  out
}
