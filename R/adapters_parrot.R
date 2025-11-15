#' PARROT adapter for align_many()
#'
#' Provides a pairwise interface for PARROT using internal building blocks.
#' Returns a soft assignment/transport map suitable for permutation sync.

#' Construct a PARROT aligner descriptor
#' @return an object of class c("parrot_aligner", "aligner")
#' @export
parrot_aligner <- function() {
  new_aligner("parrot", group = "perm", supports_multi = FALSE)
}

# Normalize rows helper (avoid division by zero)
.normalize_rows <- function(M) {
  rs <- rowSums(M)
  rs[!is.finite(rs) | rs == 0] <- 1
  M / rs
}

#' Fit PARROT on a pair of domains (adapter)
#' @export
fit_pair.parrot_aligner <- function(algo, X_i, X_j, links = NULL,
                                    ncomp = NULL,
                                    sigma = 0.15,
                                    lambda = 0.1,
                                    lambda_e = NULL,
                                    lambda_n = NULL,
                                    lambda_p = NULL,
                                    tau = 0.05,
                                    alpha = 0.2,
                                    gamma = 0.1,
                                    solver = c("sinkhorn"),
                                    max_iter = 100,
                                    tol = 1e-6,
                                    use_cpp = FALSE,
                                    anchors_policy = c("warn","off","error"),
                                    unsup_strategy = c("mnn_bootstrap","attribute_only","self_training"),
                                    conf_threshold = 0.9,
                                    max_rounds = 1,
                                    seed_strength = 0.05,
                                    laplacian_reg = NULL,
                                    ...) {
  solver <- match.arg(solver)
  anchors_policy <- match.arg(anchors_policy)
  unsup_strategy <- match.arg(unsup_strategy)
  if (is.null(lambda_e)) lambda_e <- lambda * 0.5
  if (is.null(lambda_n)) lambda_n <- if (!is.null(laplacian_reg)) laplacian_reg else lambda * 0.5
  if (is.null(lambda_p)) lambda_p <- lambda

  strata <- list(list(x = as.matrix(X_i)), list(x = as.matrix(X_j)))
  networks <- extract_parrot_networks(strata)

  n1 <- nrow(networks[[1]]$features)
  n2 <- nrow(networks[[2]]$features)

  # Build anchor_info from links, or default to empty/uniform prior
  anchor_info <- list(
    vec1 = rep(NA_real_, n1),
    vec2 = rep(NA_real_, n2),
    idx1 = integer(0),
    idx2 = integer(0),
    n1 = n1,
    n2 = n2
  )
  if (!is.null(links)) {
    if (is.list(links) && !is.null(links$vec1) && !is.null(links$vec2)) {
      v1 <- links$vec1; v2 <- links$vec2
    } else if (is.list(links) && length(links) == 2 && is.vector(links[[1]]) && is.vector(links[[2]])) {
      v1 <- links[[1]]; v2 <- links[[2]]
    } else {
      stop("links must be list(vec1=, vec2=) or list(v1, v2)")
    }
    if (length(v1) != n1 || length(v2) != n2) {
      stop("links vectors must match domain sizes")
    }
    anchor_info$vec1 <- v1
    anchor_info$vec2 <- v2
    anchor_info$idx1 <- which(!is.na(v1))
    anchor_info$idx2 <- which(!is.na(v2))
  }

  # If no anchors were provided, decide policy
  if (length(anchor_info$idx1) == 0 || length(anchor_info$idx2) == 0) {
    if (anchors_policy == "error") {
      stop("PARROT: anchors required but none provided or discovered.")
    }
    if (anchors_policy == "warn") {
      warning("PARROT: no anchors provided; using unsupervised seeding (", unsup_strategy, ")")
    }
    if (unsup_strategy == "mnn_bootstrap") {
      seeds <- .parrot_mnn_seeds(networks[[1]]$features, networks[[2]]$features)
      if (!is.null(seeds)) {
        anchor_info$vec1 <- seeds$vec1
        anchor_info$vec2 <- seeds$vec2
        anchor_info$idx1 <- which(!is.na(seeds$vec1))
        anchor_info$idx2 <- which(!is.na(seeds$vec2))
        # Increase anchor prior modestly if seeds were found
        lambda_p <- max(lambda_p, seed_strength)
      }
    }
  }

  # RWR features: if anchors exist, compute RWR; else constant columns to downweight RWR cost
  if (length(anchor_info$idx1) > 0 && length(anchor_info$idx2) > 0) {
    rwr_features <- compute_parrot_rwr(networks, anchor_info, sigma, max_iter, tol, use_cpp)
  } else {
    rwr_features <- list(matrix(1, n1, 1), matrix(1, n2, 1))
  }

  tr <- solve_parrot_transport(networks, rwr_features, anchor_info,
                               lambda_e, lambda_n, lambda_p,
                               tau, alpha, sigma, gamma, solver,
                               max_iter, tol = tol, use_cpp = use_cpp)

  structure(list(
    transport = tr$transport_plan,
    objective = tr$final_objective %||% NA_real_,
    n1 = n1, n2 = n2
  ), class = "parrot_pair_fit")
}

# Compute mutual-NN seed anchors from attributes (whitened)
.parrot_mnn_seeds <- function(X1, X2, k = NULL, ratio = 0.9) {
  X1 <- scale(as.matrix(X1))
  X2 <- scale(as.matrix(X2))
  n1 <- nrow(X1); n2 <- nrow(X2)
  if (is.null(k)) k <- max(1L, min(10L, floor(sqrt(min(n1, n2)))))
  # Find knn in the other domain
  if (!requireNamespace("RANN", quietly = TRUE)) return(NULL)
  nn12 <- RANN::nn2(X2, X1, k = k)
  nn21 <- RANN::nn2(X1, X2, k = k)

  # Mutual nearest neighbors with simple ratio test on distances
  pairs <- list()
  for (i in seq_len(n1)) {
    j <- nn12$nn.idx[i, 1]
    # Check mutual
    if (i %in% nn21$nn.idx[j, ]) {
      d1 <- nn12$nn.dists[i, 1]
      d2 <- if (k >= 2) nn12$nn.dists[i, 2] else Inf
      if (is.finite(d2) && d1 / d2 > ratio) next
      pairs[[length(pairs) + 1L]] <- c(i, j)
    }
  }
  if (!length(pairs)) return(NULL)
  pairs <- do.call(rbind, pairs)
  # Build anchor label vectors
  vec1 <- rep(NA_real_, n1); vec2 <- rep(NA_real_, n2)
  for (lab in seq_len(nrow(pairs))) {
    i <- pairs[lab, 1]; j <- pairs[lab, 2]
    vec1[i] <- lab; vec2[j] <- lab
  }
  list(vec1 = vec1, vec2 = vec2)
}

#' Extract relative transform from a PARROT pairwise fit
#' @export
relative_transform.parrot_pair_fit <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) stop("from and to must differ")
  S <- fit$transport  # n1 x n2 (i -> j mass)
  op <- if (from == "i" && to == "j") {
    # Row operator acting on rows: want n_j x n_i; use row-normalized transpose for barycentric mapping
    .normalize_rows(t(S))
  } else {
    .normalize_rows(S)
  }
  new_align_transform("perm", op, from = from, to = to, dim = c(fit$n1, fit$n2))
}

#' Pairwise loss for PARROT fit
#' @export
pair_loss.parrot_pair_fit <- function(fit, X_i = NULL, X_j = NULL, ...) {
  fit$objective %||% NA_real_
}
