# Dense small-data reference for Hong et al. (2019) LeMA -------------------
#
# This is an independent implementation of Eqs. (2)-(11) and Algorithms 1-4
# in Hong et al., ISPRS JPRS 147 (2019), 193-205,
# doi:10.1016/j.isprsjprs.2018.10.006. It is intentionally stored under
# inst/benchmarks: it is an accuracy oracle, not a production backend. The
# dense joint graph and dense feature-space solves make its memory quadratic.

lema_reference_validate_matrix <- function(X, name) {
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("`", name, "` must be a numeric matrix.", call. = FALSE)
  }
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 1L || ncol(X) < 1L || any(!is.finite(X))) {
    stop("`", name, "` must be non-empty, numeric, and finite.", call. = FALSE)
  }
  X
}

lema_reference_scale_fit <- function(X) {
  center <- colMeans(X)
  scale <- apply(X, 2L, stats::sd)
  scale[!is.finite(scale) | scale < 1e-12] <- 1
  list(
    x = sweep(sweep(X, 2L, center, "-"), 2L, scale, "/"),
    center = center,
    scale = scale
  )
}

lema_reference_one_hot <- function(labels, levels) {
  index <- match(labels, levels)
  output <- matrix(0, length(levels), length(labels),
                   dimnames = list(levels, NULL))
  output[cbind(index, seq_along(index))] <- 1
  output
}

lema_reference_pairwise_sqdist <- function(H) {
  norm <- colSums(H * H)
  distance <- outer(norm, norm, "+") - 2 * crossprod(H)
  distance[distance < 0 & distance > -1e-10] <- 0
  distance
}

lema_reference_capped_allocation <- function(cost, cap, mass) {
  if (!length(cost)) return(numeric(0))
  cap <- rep_len(as.numeric(cap), length(cost))
  cap[!is.finite(cap) | cap < 0] <- 0
  mass <- min(max(as.numeric(mass), 0), sum(cap))
  output <- numeric(length(cost))
  if (mass <= 0) return(output)
  ordering <- order(cost, seq_along(cost), method = "radix")
  remaining <- mass
  for (index in ordering) {
    amount <- min(cap[[index]], remaining)
    output[[index]] <- amount
    remaining <- remaining - amount
    if (remaining <= 1e-12) break
  }
  output
}

lema_reference_supervised_block <- function(labels) {
  counts <- table(labels)
  outer(labels, labels, function(a, b) as.numeric(a == b)) /
    matrix(as.numeric(counts[labels]), length(labels), length(labels))
}

lema_reference_update_graph <- function(H, labels, graph_mass = NULL) {
  n <- length(labels)
  total <- ncol(H)
  nu <- total - 2L * n
  if (nu < 0L) stop("LeMA graph dimensions are inconsistent.", call. = FALSE)
  class_count <- table(labels)
  row_cap <- 1 / as.numeric(class_count[labels])
  if (is.null(graph_mass)) graph_mass <- max(1, min(n, max(nu, 1L)))
  graph_mass <- as.numeric(graph_mass)
  Z <- lema_reference_pairwise_sqdist(H)
  W <- matrix(0, total, total)
  h <- seq_len(n)
  m <- n + seq_len(n)
  fixed <- lema_reference_supervised_block(labels)
  W[h, h] <- fixed
  W[h, m] <- fixed
  W[m, h] <- fixed
  W[m, m] <- fixed

  learned_mass <- c(high_unlabeled = 0, low_unlabeled = 0, unlabeled = 0)
  if (nu > 0L) {
    u <- 2L * n + seq_len(nu)
    cap_cross <- matrix(row_cap, n, nu)
    high <- matrix(
      lema_reference_capped_allocation(
        Z[h, u, drop = FALSE], cap_cross, graph_mass
      ),
      n,
      nu
    )
    low <- matrix(
      lema_reference_capped_allocation(
        Z[m, u, drop = FALSE], cap_cross, graph_mass
      ),
      n,
      nu
    )
    # Algorithm 1 aligns W_HU and W_MU by their elementwise maximum.
    aligned <- pmax(high, low)
    W[h, u] <- aligned
    W[m, u] <- aligned
    W[u, h] <- t(aligned)
    W[u, m] <- t(aligned)
    learned_mass[["high_unlabeled"]] <- sum(aligned)
    learned_mass[["low_unlabeled"]] <- sum(aligned)

    if (nu > 1L) {
      upper <- which(upper.tri(matrix(FALSE, nu, nu)), arr.ind = TRUE)
      upper_cost <- Z[u, u, drop = FALSE][upper]
      cap <- 1 / min(as.numeric(class_count))
      value <- lema_reference_capped_allocation(
        upper_cost, cap, graph_mass / 2
      )
      Wuu <- matrix(0, nu, nu)
      Wuu[upper] <- value
      Wuu <- Wuu + t(Wuu)
      W[u, u] <- Wuu
      learned_mass[["unlabeled"]] <- sum(Wuu)
    }
  }
  W <- (W + t(W)) / 2
  L <- diag(rowSums(W), nrow(W)) - W
  list(
    W = W,
    L = L,
    distance = Z,
    learned_mass = learned_mass,
    fixed_block = fixed
  )
}

lema_reference_update_p <- function(theta, X_label, Y, alpha) {
  E <- theta %*% X_label
  system <- E %*% t(E) + alpha * diag(nrow(E))
  Y %*% t(E) %*% solve(system)
}

lema_reference_polar_rows <- function(A) {
  decomposition <- svd(A, nu = nrow(A), nv = nrow(A))
  decomposition$u %*% t(decomposition$v)
}

lema_reference_update_theta <- function(
  theta,
  P,
  Y,
  X_label,
  X_all,
  L,
  beta,
  max_iter,
  tolerance
) {
  d <- nrow(theta)
  D <- ncol(theta)
  J <- theta %*% X_label
  G <- theta
  lambda1 <- matrix(0, d, ncol(X_label))
  lambda2 <- matrix(0, d, D)
  mu <- 1e-3
  residual <- Inf
  iteration <- 0L
  graph_covariance <- X_all %*% L %*% t(X_all)
  graph_covariance <- (graph_covariance + t(graph_covariance)) / 2

  for (iteration in seq_len(max_iter)) {
    J <- solve(
      crossprod(P) + mu * diag(d),
      t(P) %*% Y + mu * theta %*% X_label - lambda1
    )
    numerator <- mu * J %*% t(X_label) +
      lambda1 %*% t(X_label) + mu * G + lambda2
    denominator <- mu * X_label %*% t(X_label) +
      mu * diag(D) + beta * graph_covariance
    denominator <- (denominator + t(denominator)) / 2
    theta <- t(solve(denominator + 1e-10 * diag(D), t(numerator)))
    G <- lema_reference_polar_rows(theta - lambda2 / mu)
    residual_j <- sqrt(sum((J - theta %*% X_label)^2))
    residual_g <- sqrt(sum((G - theta)^2))
    residual <- max(residual_j, residual_g)
    lambda1 <- lambda1 + mu * (J - theta %*% X_label)
    lambda2 <- lambda2 + mu * (G - theta)
    mu <- min(1.5 * mu, 1e6)
    if (residual <= tolerance) break
  }
  # Return the orthogonal auxiliary variable, which satisfies Eq. (2).
  list(theta = G, iterations = as.integer(iteration), residual = residual)
}

lema_reference_objective <- function(
  P, theta, X_label, X_all, Y, W, alpha, beta
) {
  fidelity <- 0.5 * sum((Y - P %*% theta %*% X_label)^2)
  ridge <- alpha / 2 * sum(P * P)
  distance <- lema_reference_pairwise_sqdist(theta %*% X_all)
  graph <- beta / 4 * sum(W * distance)
  c(total = fidelity + ridge + graph,
    fidelity = fidelity, ridge = ridge, graph = graph)
}

lema_reference_softmax <- function(score) {
  score <- sweep(score, 1L, apply(score, 1L, max), "-")
  probability <- exp(score)
  probability / pmax(rowSums(probability), 1e-15)
}

# Small-data LeMA reference. Inputs follow the paper's completed paired-label
# setting: X_high and X_low_labeled contain N corresponding labeled samples;
# X_low_unlabeled supplies the target-domain unlabeled landmarks/test rows.
lema_reference_fit <- function(
  X_high,
  X_low_labeled,
  X_low_unlabeled,
  labels,
  ncomp = 10L,
  alpha = 0.01,
  beta = 0.01,
  graph_mass = NULL,
  max_iter = 20L,
  tolerance = 1e-5,
  admm_max_iter = 100L,
  admm_tolerance = 1e-6,
  seed = 1L,
  verbose = FALSE
) {
  X_high <- lema_reference_validate_matrix(X_high, "X_high")
  X_low_labeled <- lema_reference_validate_matrix(
    X_low_labeled, "X_low_labeled"
  )
  if (is.null(X_low_unlabeled)) {
    X_low_unlabeled <- matrix(numeric(0), 0L, ncol(X_low_labeled))
  } else if (!is.matrix(X_low_unlabeled) || nrow(X_low_unlabeled) == 0L) {
    X_low_unlabeled <- matrix(numeric(0), 0L, ncol(X_low_labeled))
  } else {
    X_low_unlabeled <- lema_reference_validate_matrix(
      X_low_unlabeled, "X_low_unlabeled"
    )
  }
  n <- nrow(X_high)
  if (nrow(X_low_labeled) != n || length(labels) != n) {
    stop("High- and low-modality labeled rows and `labels` must have equal length.",
         call. = FALSE)
  }
  if (ncol(X_low_unlabeled) != ncol(X_low_labeled)) {
    stop("Labeled and unlabeled low-modality rows need equal feature width.",
         call. = FALSE)
  }
  labels <- as.character(labels)
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("The paired LeMA reference requires complete non-empty labels.",
         call. = FALSE)
  }
  scalar_positive <- function(value, name, zero = FALSE) {
    valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
      if (zero) value >= 0 else value > 0
    if (!valid) stop("`", name, "` must be a finite ",
                     if (zero) "nonnegative" else "positive", " scalar.",
                     call. = FALSE)
    as.numeric(value)
  }
  alpha <- scalar_positive(alpha, "alpha")
  beta <- scalar_positive(beta, "beta", zero = TRUE)
  tolerance <- scalar_positive(tolerance, "tolerance")
  admm_tolerance <- scalar_positive(admm_tolerance, "admm_tolerance")
  ncomp <- as.integer(ncomp)
  max_iter <- as.integer(max_iter)
  admm_max_iter <- as.integer(admm_max_iter)
  seed <- as.integer(seed)
  if (!is.finite(ncomp) || ncomp < 1L || !is.finite(max_iter) || max_iter < 1L ||
      !is.finite(admm_max_iter) || admm_max_iter < 1L) {
    stop("`ncomp`, `max_iter`, and `admm_max_iter` must be positive integers.",
         call. = FALSE)
  }

  high_scale <- lema_reference_scale_fit(X_high)
  low_all <- rbind(X_low_labeled, X_low_unlabeled)
  low_scale <- lema_reference_scale_fit(low_all)
  high <- high_scale$x
  low <- low_scale$x
  nu <- nrow(X_low_unlabeled)
  D <- ncol(high) + ncol(low)
  ncomp <- min(ncomp, D, 2L * n + nu)
  X_label <- rbind(
    cbind(t(high), matrix(0, ncol(high), n)),
    cbind(matrix(0, ncol(low), n), t(low[seq_len(n), , drop = FALSE]))
  )
  X_all <- rbind(
    cbind(t(high), matrix(0, ncol(high), n + nu)),
    cbind(matrix(0, ncol(low), n), t(low))
  )
  class_levels <- sort(unique(labels), method = "radix")
  Y_single <- lema_reference_one_hot(labels, class_levels)
  Y <- cbind(Y_single, Y_single)

  set.seed(seed)
  initial <- svd(X_all, nu = ncomp, nv = 0L)$u[, seq_len(ncomp), drop = FALSE]
  theta <- t(initial)
  graph <- lema_reference_update_graph(
    theta %*% X_all, labels, graph_mass = graph_mass
  )
  P <- lema_reference_update_p(theta, X_label, Y, alpha)
  history <- vector("list", max_iter)
  status <- "max_iter"
  previous <- Inf
  theta_diagnostics <- NULL

  for (iteration in seq_len(max_iter)) {
    P <- lema_reference_update_p(theta, X_label, Y, alpha)
    theta_diagnostics <- lema_reference_update_theta(
      theta, P, Y, X_label, X_all, graph$L, beta,
      max_iter = admm_max_iter, tolerance = admm_tolerance
    )
    theta <- theta_diagnostics$theta
    graph <- lema_reference_update_graph(
      theta %*% X_all, labels, graph_mass = graph_mass
    )
    P <- lema_reference_update_p(theta, X_label, Y, alpha)
    objective <- lema_reference_objective(
      P, theta, X_label, X_all, Y, graph$W, alpha, beta
    )
    relative <- if (is.finite(previous)) {
      abs(previous - objective[["total"]]) / max(abs(previous), 1e-12)
    } else {
      Inf
    }
    history[[iteration]] <- data.frame(
      iteration = iteration,
      total = objective[["total"]],
      fidelity = objective[["fidelity"]],
      ridge = objective[["ridge"]],
      graph = objective[["graph"]],
      relative_change = relative,
      theta_admm_iterations = theta_diagnostics$iterations,
      theta_admm_residual = theta_diagnostics$residual,
      stringsAsFactors = FALSE
    )
    if (isTRUE(verbose)) {
      message("LeMA reference iteration ", iteration,
              ": objective=", format(objective[["total"]], digits = 7))
    }
    if (is.finite(relative) && relative <= tolerance) {
      status <- "converged"
      history <- history[seq_len(iteration)]
      break
    }
    previous <- objective[["total"]]
  }
  history <- do.call(rbind, history[!vapply(history, is.null, logical(1))])
  H <- theta %*% X_all
  target_index <- n + seq_len(n + nu)
  source_embedding <- t(H[, seq_len(n), drop = FALSE])
  target_embedding <- t(H[, target_index, drop = FALSE])
  score <- t(P %*% H)
  probability <- lema_reference_softmax(score)
  colnames(probability) <- class_levels
  prediction <- class_levels[max.col(probability, ties.method = "first")]

  structure(
    list(
      theta = theta,
      P = P,
      W = graph$W,
      L = graph$L,
      source_embedding = source_embedding,
      target_embedding = target_embedding,
      target_labeled_rows = seq_len(n),
      target_unlabeled_rows = n + seq_len(nu),
      probabilities = probability,
      predictions = prediction,
      class_levels = class_levels,
      history = history,
      objective = utils::tail(history$total, 1L),
      convergence = list(
        status = status,
        iterations = nrow(history),
        theta_admm = theta_diagnostics
      ),
      constraints = list(
        orthogonality_error = max(abs(theta %*% t(theta) - diag(ncomp))),
        graph_symmetry_error = max(abs(graph$W - t(graph$W))),
        graph_minimum = min(graph$W),
        laplacian_row_sum_error = max(abs(rowSums(graph$L))),
        learned_mass = graph$learned_mass
      ),
      preprocessing = list(
        high_center = high_scale$center,
        high_scale = high_scale$scale,
        low_center = low_scale$center,
        low_scale = low_scale$scale,
        transductive_low_scaling = TRUE
      ),
      parameters = list(
        ncomp = ncomp,
        alpha = alpha,
        beta = beta,
        graph_mass = graph_mass,
        seed = seed
      ),
      provenance = list(
        reference = "Hong et al. (2019), Eqs. (2)-(11), Algorithms 1-4",
        doi = "10.1016/j.isprsjprs.2018.10.006",
        implementation = "independent equation-level reimplementation",
        production_backend = FALSE,
        dense_joint_graph = TRUE
      )
    ),
    class = "lema_reference_fit"
  )
}
