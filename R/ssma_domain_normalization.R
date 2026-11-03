# Product-Stiefel alternative to SSMA's globally normalized eigenproblem.
# This is a different, nonconvex objective constraint, not an eigenvector rescale.
.ssma_retract <- function(x) {
  # Positive-diagonal, unpivoted QR is a smooth Stiefel retraction.
  # Scaling prevents large rejected trial steps overflowing QR.
  scale <- max(abs(x))
  if (!is.finite(scale) || scale == 0) stop("Nonfinite or zero Stiefel trial.", call. = FALSE)
  z <- qr(x / scale, tol = 0, LAPACK = FALSE)
  q <- qr.Q(z, complete = FALSE)
  signs <- sign(diag(qr.R(z))); signs[signs == 0] <- 1
  q <- sweep(q, 2, signs, "*")
  if (any(!is.finite(q)) || max(abs(crossprod(q) - diag(ncol(x)))) > 1e-8) {
    stop("Stiefel QR retraction failed.", call. = FALSE)
  }
  q
}

.ssma_domain_gradient <- function(Y, G, ii) {
  do.call(rbind, lapply(ii, function(j) {
    y <- Y[j, , drop = FALSE]; g <- G[j, , drop = FALSE]
    c <- crossprod(y, g)
    g - y %*% ((c + t(c)) / 2)
  }))
}

# Square Stiefel blocks have two disconnected orthogonal components. Their
# diagonal quadratic is constant, so a polar update solves a block exactly,
# allowing reflections that tangent steps cannot reach.
.ssma_update_square_blocks <- function(Y, H, ii) {
  updates <- 0L
  for (j in ii) {
    if (length(j) != ncol(Y)) next
    other <- setdiff(seq_len(nrow(Y)), j)
    C <- H[j, other, drop = FALSE] %*% Y[other, , drop = FALSE]
    sv <- svd(-C)
    candidate <- sv$u %*% t(sv$v)
    old <- sum(Y[j, , drop = FALSE] * C)
    new <- sum(candidate * C)
    if (old - new > 32 * .Machine$double.eps * max(1, abs(old), abs(new))) {
      Y[j, ] <- candidate
      updates <- updates + 1L
    }
  }
  list(Y = Y, updates = updates)
}

.ssma_solve_domain_scores <- function(reduced, ncomp, tol, maxiter) {
  q <- as.integer(ncomp)
  metric <- lapply(reduced$basis, function(b) {
    e <- eigen(crossprod(b$T), symmetric = TRUE)
    keep <- which(e$values > max(e$values) * 1e-10 & e$values > 0)
    if (length(keep) < q) {
      stop("domain_scores normalization needs at least ncomp independent score directions in every domain.", call. = FALSE)
    }
    sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "/")
  })
  # Rectangular maps discard only numerical score-null coordinates.
  Q <- as.matrix(Matrix::bdiag(metric))
  H <- crossprod(Q, reduced$A_r %*% Q)
  H <- (H + t(H)) / 2
  ranks <- vapply(metric, ncol, integer(1))
  offsets <- c(0L, cumsum(ranks))
  ii <- lapply(seq_along(ranks), function(i) seq.int(offsets[i] + 1L, offsets[i + 1L]))

  # Remove feasible-set constants before scaling or measuring stationarity.
  # Their contribution remains in all reported objective values.
  constant <- 0
  for (j in ii) {
    block <- H[j, j, drop = FALSE]
    if (length(j) == q) {
      constant <- constant + sum(diag(block))
      H[j, j] <- 0
    } else {
      shift <- mean(diag(block))
      constant <- constant + q * shift
      diag(block) <- diag(block) - shift
      H[j, j] <- block
    }
  }
  scale <- max(rowSums(abs(H)))
  if (!is.finite(scale)) stop("Nonfinite domain_scores objective.", call. = FALSE)
  if (scale == 0) scale <- 1
  H <- H / scale
  retract <- function(Y) do.call(rbind, lapply(ii, function(j) .ssma_retract(Y[j, , drop = FALSE])))
  # Spectral directions supply orientation; retract before evaluating cost.
  e <- eigen(H, symmetric = TRUE)
  seed <- e$vectors[, order(e$values)[seq_len(q)], drop = FALSE]
  for (j in ii) {
    if (qr(seed[j, , drop = FALSE], tol = 1e-10)$rank < q) {
      local <- eigen(H[j, j, drop = FALSE], symmetric = TRUE)
      seed[j, ] <- local$vectors[, order(local$values)[seq_len(q)], drop = FALSE]
    }
  }
  Y <- retract(seed)
  HY <- H %*% Y
  f <- sum(Y * HY)
  initial <- f
  history <- f
  step <- 1
  converged <- FALSE
  square_updates <- 0L
  for (iter in seq_len(maxiter)) {
    G <- .ssma_domain_gradient(Y, 2 * HY, ii)
    gnorm <- sqrt(sum(G^2))
    relative_gradient <- gnorm / max(1, sqrt(sum((2 * HY)^2)))
    if (relative_gradient <= tol) {
      candidate <- .ssma_update_square_blocks(Y, H, ii)
      if (candidate$updates == 0L) { converged <- TRUE; break }
      Y <- candidate$Y; HY <- H %*% Y; f <- sum(Y * HY)
      square_updates <- square_updates + candidate$updates
      history <- c(history, f)
      step <- 1
      next
    }
    accepted <- FALSE
    for (ls in seq_len(60L)) {
      next_Y <- tryCatch(retract(Y - step * G), error = function(e) NULL)
      if (is.null(next_Y)) { step <- step / 2; next }
      next_HY <- H %*% next_Y
      next_f <- sum(next_Y * next_HY)
      if (is.finite(next_f) && next_f <= f - 1e-4 * step * gnorm^2) {
        accepted <- TRUE; break
      }
      step <- step / 2
    }
    if (!accepted) break
    next_G <- .ssma_domain_gradient(next_Y, 2 * next_HY, ii)
    delta <- next_Y - Y; dg <- next_G - G
    sy <- sum(delta * dg)
    # Safeguarded Barzilai-Borwein trial; Armijo enforces descent.
    step <- if (sy > 0) min(1e8, max(1e-8, sum(delta^2) / sy)) else 1
    Y <- next_Y; HY <- next_HY; f <- next_f
    history <- c(history, f)
  }
  G <- .ssma_domain_gradient(Y, 2 * HY, ii)
  relative_gradient <- sqrt(sum(G^2)) / max(1, sqrt(sum((2 * HY)^2)))
  if (!converged) {
    stop(sprintf("domain_scores SSMA did not converge (relative tangent gradient %.3g after %d iterations).", relative_gradient, iter), call. = FALSE)
  }
  alpha <- Q %*% Y
  blocks <- .ssma_split_blocks(alpha, reduced$ranks)
  error <- max(vapply(seq_along(blocks), function(i) {
    max(abs(crossprod(reduced$basis[[i]]$T %*% blocks[[i]]) - diag(q)))
  }, numeric(1)))
  if (!is.finite(error) || error > 1e-7) stop("domain_scores score normalization failed.", call. = FALSE)
  list(alpha = alpha, eigenvalues = NULL, eigenvalues_all = NULL, ridge = 0,
       domain_diagnostics = list(converged = converged, iterations = iter,
         relative_tangent_gradient = relative_gradient, score_gram_max_error = error,
         objective = f * scale + constant, initial_objective = initial * scale + constant,
         objective_history = history * scale + constant, retained_score_ranks = ranks,
         objective_scale = scale, objective_constant = constant,
         square_block_updates = square_updates,
         interpretation = "Stationary point on a product of score-normalized Stiefel manifolds with full-rank block orientation checks; no global-optimum guarantee."))
}
