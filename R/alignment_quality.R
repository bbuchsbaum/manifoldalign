#' Alignment Quality Metrics
#'
#' Compute compact, label-free quality diagnostics for multiset alignment
#' results. Supports mma_align_multiple in both reference and consensus modes.
#'
#' Metrics include per-domain mean squared alignment error under final
#' assignments (if available) and a soft counterpart using posteriors. For
#' reference mode, errors are measured to the reference domain; for consensus
#' mode, errors are measured to the consensus centers.
#'
#' @param object An alignment result, e.g. from \code{mma_align_multiple}
#' @return A list with fields:
#'   - mode: "reference" or "consensus"
#'   - per_domain: data.frame with per-domain metrics (n, mse, mse_soft,
#'                 matched, coverage, mean_confidence, converged, iterations)
#'   - global: aggregate metrics (mse_weighted, mse_soft_weighted, K, ref_idx,
#'             mean_confidence, all_converged)
#' @examples
#' \donttest{
#' # Simulate three domains with shared structure
#' set.seed(42)
#' X1 <- matrix(rnorm(100), 50, 2)
#' X2 <- matrix(rnorm(100), 50, 2)
#' X3 <- matrix(rnorm(100), 50, 2)
#'
#' # Create hyperdesign object
#' hd <- as_hyperdesign(list(X1, X2, X3))
#'
#' # Perform alignment (requires mma_align_multiple)
#' # result <- mma_align_multiple(hd, ncomp = 2)
#' # quality <- alignment_quality(result)
#' }
#' @export
alignment_quality <- function(object) {
  UseMethod("alignment_quality")
}

#' @export
alignment_quality.default <- function(object) {
  stop("alignment_quality: unsupported object type.")
}

#' @export
alignment_quality.mma_align_multiple <- function(object) {
  s <- object$s
  K <- object$K %||% ncol(s)
  post <- object$posteriors
  stopifnot(is.list(post), length(post) >= 2)

  mode <- if (!is.null(object$consensus)) "consensus" else "reference"
  ref_idx <- object$ref_idx %||% 1L
  m <- length(post)

  # -- Determine per-domain sample counts and row splits in scores --
  n_vec <- integer(m)
  if (mode == "reference") {
    for (g in seq_len(m)) {
      if (g == ref_idx) {
        n_vec[g] <- nrow(post[[g]])
      } else {
        n_vec[g] <- ncol(post[[g]])
      }
    }
  } else {
    for (g in seq_len(m)) n_vec[g] <- ncol(post[[g]])
  }
  if (sum(n_vec) != nrow(s)) {
    # Fallback: evenly partition (should not happen)
    warning("alignment_quality: could not infer row splits from posteriors; using fallback.")
    avg <- floor(nrow(s) / m)
    n_vec[] <- avg
    n_vec[1] <- n_vec[1] + (nrow(s) - sum(n_vec))
  }
  starts <- c(1L, cumsum(n_vec) + 1L); starts <- starts[-length(starts)]
  ends <- cumsum(n_vec)
  rng <- Map(seq.int, starts, ends)

  # -- Prepare confidence per domain from signature_align if present --
  aligns <- object$signature_align %||% object$hist_align
  conf <- rep(NA_real_, m)
  if (is.list(aligns) && length(aligns) == m) {
    for (g in seq_len(m)) {
      pc <- tryCatch(aligns[[g]]$pair_confidence, error = function(e) NULL)
      if (!is.null(pc)) conf[g] <- mean(pc, na.rm = TRUE)
    }
  }

  # -- EM stats per domain (reference mode has per-domain, consensus has global) --
  ems <- object$em_stats
  converged <- rep(NA, m); iters <- rep(NA_integer_, m)
  if (mode == "reference" && is.list(ems) && length(ems) == m) {
    for (g in seq_len(m)) {
      converged[g] <- isTRUE(ems[[g]]$converged)
      iters[g] <- as.integer(ems[[g]]$iterations %||% NA_integer_)
    }
  } else if (mode == "consensus" && is.list(ems)) {
    converged[] <- isTRUE(ems$converged)
    iters[] <- as.integer(ems$iterations %||% NA_integer_)
  }

  # -- Compute per-domain errors --
  mse <- rep(NA_real_, m)
  mse_soft <- rep(NA_real_, m)
  matched <- rep(NA_integer_, m)
  coverage <- rep(NA_real_, m)

  if (mode == "reference") {
    Eref <- s[rng[[ref_idx]], , drop = FALSE]
    rr <- rowSums(Eref^2)
    for (g in seq_len(m)) {
      Eg <- s[rng[[g]], , drop = FALSE]
      if (g == ref_idx) next
      # Soft MSE via posteriors
      A <- post[[g]] # nr x ng
      yg <- rowSums(Eg^2)
      D2 <- outer(rr, yg, "+") - 2 * (Eref %*% t(Eg))
      denom <- sum(A)
      if (is.finite(denom) && denom > 0) {
        mse_soft[g] <- sum(A * D2) / (K * denom)
      }
      # Hard MSE via final assignment if available
      perm <- NULL
      if (!is.null(object$assignment) && length(object$assignment) >= g) {
        perm <- object$assignment[[g]]
      }
      if (!is.null(perm) && length(perm) == nrow(Eg)) {
        ok <- which(!is.na(perm) & perm >= 1 & perm <= nrow(Eref))
        matched[g] <- length(ok)
        if (matched[g] > 0) {
          coverage[g] <- matched[g] / length(perm)
          dif <- Eref[perm[ok], , drop = FALSE] - Eg[ok, , drop = FALSE]
          mse[g] <- mean(dif^2)
        }
      }
    }
  } else {
    C <- object$consensus
    cr <- rowSums(C^2)
    for (g in seq_len(m)) {
      Eg <- s[rng[[g]], , drop = FALSE]
      yg <- rowSums(Eg^2)
      D2 <- outer(cr, yg, "+") - 2 * (C %*% t(Eg))
      # Soft MSE via posteriors
      A <- post[[g]] # n_centers x n_g
      denom <- sum(A)
      if (is.finite(denom) && denom > 0) {
        mse_soft[g] <- sum(A * D2) / (K * denom)
      }
      # Hard MSE via assignment if available
      perm <- NULL
      if (!is.null(object$assignment) && length(object$assignment) >= g) {
        perm <- object$assignment[[g]]
      }
      if (!is.null(perm) && length(perm) == nrow(Eg)) {
        ok <- which(!is.na(perm) & perm >= 1 & perm <= nrow(C))
        matched[g] <- length(ok)
        if (matched[g] > 0) {
          coverage[g] <- matched[g] / length(perm)
          dif <- C[perm[ok], , drop = FALSE] - Eg[ok, , drop = FALSE]
          mse[g] <- mean(dif^2)
        }
      }
    }
  }

  # -- Aggregate --
  w_m <- ifelse(is.na(matched) | matched <= 0, 0, matched)
  w_s <- n_vec
  mse_weighted <- if (sum(w_m) > 0) sum(w_m * replace_na(mse, 0)) / sum(w_m) else NA_real_
  mse_soft_weighted <- if (sum(w_s) > 0) sum(w_s * replace_na(mse_soft, 0)) / sum(w_s) else NA_real_
  mean_conf <- if (all(is.na(conf))) NA_real_ else mean(conf, na.rm = TRUE)
  all_conv <- if (all(is.na(converged))) NA else all(converged, na.rm = TRUE)

  per_domain <- data.frame(
    domain = seq_len(m),
    n = n_vec,
    mse = mse,
    mse_soft = mse_soft,
    matched = matched,
    coverage = coverage,
    mean_confidence = conf,
    converged = converged,
    iterations = iters,
    stringsAsFactors = FALSE
  )

  list(
    mode = mode,
    per_domain = per_domain,
    global = list(
      mse_weighted = mse_weighted,
      mse_soft_weighted = mse_soft_weighted,
      K = K,
      ref_idx = ref_idx,
      mean_confidence = mean_conf,
      all_converged = all_conv
    )
  )
}

`%||%` <- function(a, b) if (is.null(a)) b else a

replace_na <- function(x, value) { x[is.na(x)] <- value; x }

