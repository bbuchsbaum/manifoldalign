#' Token-level Optimal Transport Graph Alignment
#'
#' Align two graphs (node sets) by representing each node as a **bag of context
#' tokens** and aligning nodes via:
#' 1) multi-view node embeddings (late-fusion similarity), and
#' 2) token-level entropic optimal transport (Sinkhorn) between node bags.
#'
#' This implementation is designed to stay sparse and efficient by generating
#' a top-k candidate set per source node (ANN / kNN) and running Sinkhorn only
#' on the candidate bipartite graph.
#'
#' @param data A `hyperdesign` object (2 domains) or a list of two matrices.
#' @param anchors Optional unquoted column name in each domain's design table
#'   providing anchor identifiers (shared values indicate known correspondences).
#'   Use `NULL` for unsupervised alignment.
#' @param preproc Optional preprocessing (default: `multivarious::center()`).
#' @param ncomp Integer number of components stored in `s` (scores) and `v`
#'   (loadings). Used for the baseline embedding view when available.
#' @param control Control settings from [token_ot_graph_align_control()].
#' @param ... Additional arguments (reserved for future extensions).
#'
#' @return A `multiblock_biprojector` with subclass `token_ot_graph_align`.
#'   Extras include `transport_plan` (sparse coupling), `assignment` (greedy or
#'   Hungarian), and `diagnostics`.
#' @export
token_ot_graph_align <- function(data, ...) {
  UseMethod("token_ot_graph_align")
}

#' Control settings for `token_ot_graph_align()`
#'
#' @param n_levels Integer number of hierarchy levels. `1` disables multilevel
#'   refinement.
#' @param coarsen_ratio Numeric factor > 1 controlling the coarse-to-fine size
#'   reduction (roughly `n_next ≈ n_prev / coarsen_ratio`).
#' @param coarsen_method Coarsening strategy for multilevel:
#'   `"kmeans"` (fast, feature-based) or `"louvain"` (graph community-based with
#'   fallback to k-means if the number of communities is unsuitable).
#' @param min_clusters Minimum nodes per level (stops coarsening once the graph
#'   is small).
#' @param prior_mode Prior lifting mode for multilevel: `"none"`, `"hard"` (use
#'   discrete coarse assignment), or `"soft"` (use top-k coarse coupling mass).
#' @param prior_cluster_k Integer number of coarse target clusters kept per
#'   coarse source cluster when `prior_mode="soft"`.
#' @param prior_max_candidates Maximum prior candidates lifted per fine node
#'   (cap to prevent blowups).
#' @param graph_knn Integer k for within-domain kNN graph construction.
#' @param graph_sigma Optional numeric bandwidth for heat-kernel affinities.
#'   Use `NULL` to auto-tune per domain via [choose_sigma()].
#' @param use_laplacian Logical; if TRUE use normalized Laplacian bases.
#' @param views Character vector of view names in `c("raw", "eigenmap", "hks")`.
#' @param view_weights Optional numeric vector of nonnegative weights, same
#'   length as `views` (will be normalized to sum to 1).
#' @param k_embed Optional integer number of eigenpairs used for eigenmap/HKS.
#'   If `NULL`, defaults to `max(3*ncomp, ncomp + 10)` inside the solver.
#' @param hks_q Integer number of HKS time steps.
#' @param lambda_view Mix weight in [0,1] for view-distance vs token-OT cost.
#' @param candidate_k Integer number of ANN candidates per source node.
#' @param ann_backend Nearest-neighbor backend for candidate generation:
#'   `"auto"` (default), `"rann"`, or `"annoy"` (requires RcppAnnoy).
#' @param ann_trees Number of trees used by Annoy when `ann_backend="annoy"`.
#' @param prior_k Integer number of candidates taken from prior (reserved).
#' @param ensure_cols Logical; ensure every target column has at least one
#'   incoming candidate edge (adds reverse-1NN edges if needed).
#' @param token_mode Either `"view_only"` (self token only) or
#'   `"view_plus_neighbors"` (includes sampled neighbor tokens).
#' @param max_hops Integer hop depth for neighborhood tokens.
#' @param samples_per_hop Integer number of neighbors sampled per hop.
#' @param max_tokens Maximum tokens retained per node (cap for runtime).
#' @param token_metric Token distance: `"cosine"` or `"sqeuclidean"`.
#' @param eps_token Entropic regularization for token-level OT (>0).
#' @param iters_token Maximum Sinkhorn iterations for token-level OT.
#' @param tol_token Token-level Sinkhorn convergence tolerance.
#' @param eps_node Entropic regularization for node-level sparse Sinkhorn (>0).
#' @param iters_node Maximum iterations for node-level sparse Sinkhorn.
#' @param tol_node Node-level Sinkhorn convergence tolerance.
#' @param projection Discrete projection mode: `"greedy"` or `"hungarian"`.
#' @param anchor_weight Nonnegative weight applied to an anchor mismatch penalty.
#' @param anchor_penalty Nonnegative penalty added to costs for anchor mismatch
#'   edges when both endpoints are anchored.
#' @param seed Integer seed used for any randomized sampling.
#' @param verbose Logical; print progress messages.
#'
#' @return A list of class `token_ot_graph_align_control`.
#' @export
token_ot_graph_align_control <- function(
  n_levels = 1L,
  coarsen_ratio = 4,
  coarsen_method = c("kmeans", "louvain"),
  min_clusters = 20L,
  prior_mode = c("none", "hard", "soft"),
  prior_cluster_k = 2L,
  prior_max_candidates = 256L,
  graph_knn = 15L,
  graph_sigma = NULL,
  use_laplacian = TRUE,
  views = c("raw", "eigenmap", "hks"),
  view_weights = NULL,
  k_embed = NULL,
  hks_q = 16L,
  lambda_view = 0.5,
  candidate_k = 100L,
  ann_backend = c("auto", "rann", "annoy"),
  ann_trees = 50L,
  prior_k = 0L,
  ensure_cols = TRUE,
  token_mode = c("view_only", "view_plus_neighbors"),
  max_hops = 1L,
  samples_per_hop = 8L,
  max_tokens = 16L,
  token_metric = c("cosine", "sqeuclidean"),
  eps_token = 0.05,
  iters_token = 50L,
  tol_token = 1e-7,
  eps_node = 0.05,
  iters_node = 400L,
  tol_node = 1e-7,
  projection = c("greedy", "hungarian"),
  anchor_weight = 0,
  anchor_penalty = 1,
  seed = 1L,
  verbose = FALSE
) {
  prior_mode <- match.arg(prior_mode)
  coarsen_method <- match.arg(coarsen_method)
  ann_backend <- match.arg(ann_backend)
  token_mode <- match.arg(token_mode)
  token_metric <- match.arg(token_metric)
  projection <- match.arg(projection)

  chk::chk_number(n_levels)
  chk::chk_true(n_levels >= 1)
  chk::chk_number(coarsen_ratio)
  chk::chk_true(coarsen_ratio > 1)
  chk::chk_character(coarsen_method)
  chk::chk_number(min_clusters)
  chk::chk_true(min_clusters >= 2)
  chk::chk_number(prior_cluster_k)
  chk::chk_true(prior_cluster_k >= 1)
  chk::chk_number(prior_max_candidates)
  chk::chk_true(prior_max_candidates >= 1)

  chk::chk_number(graph_knn)
  chk::chk_true(graph_knn >= 1)
  if (!is.null(graph_sigma)) {
    chk::chk_number(graph_sigma)
    chk::chk_true(is.finite(graph_sigma) && graph_sigma > 0)
  }
  chk::chk_logical(use_laplacian)

  chk::chk_character(views)
  if (!length(views)) stop("`views` must be non-empty.", call. = FALSE)
  allowed <- c("raw", "eigenmap", "hks")
  bad <- setdiff(views, allowed)
  if (length(bad)) stop("Unknown view(s): ", paste(bad, collapse = ", "), call. = FALSE)

  if (!is.null(view_weights)) {
    chk::chk_numeric(view_weights)
    chk::chk_true(length(view_weights) == length(views))
    chk::chk_true(all(is.finite(view_weights)))
    chk::chk_true(all(view_weights >= 0))
  }

  if (!is.null(k_embed)) {
    chk::chk_number(k_embed)
    chk::chk_true(k_embed >= 2)
  }
  chk::chk_number(hks_q)
  chk::chk_true(hks_q >= 1)

  chk::chk_number(lambda_view)
  chk::chk_true(lambda_view >= 0 && lambda_view <= 1)

  chk::chk_number(candidate_k)
  chk::chk_true(candidate_k >= 1)
  chk::chk_number(ann_trees)
  chk::chk_true(ann_trees >= 1)
  chk::chk_number(prior_k)
  chk::chk_true(prior_k >= 0)
  chk::chk_logical(ensure_cols)

  chk::chk_number(max_hops)
  chk::chk_true(max_hops >= 0)
  chk::chk_number(samples_per_hop)
  chk::chk_true(samples_per_hop >= 0)
  chk::chk_number(max_tokens)
  chk::chk_true(max_tokens >= 1)

  chk::chk_number(eps_token)
  chk::chk_true(eps_token > 0)
  chk::chk_number(iters_token)
  chk::chk_true(iters_token >= 1)
  chk::chk_number(tol_token)
  chk::chk_true(tol_token > 0)

  chk::chk_number(eps_node)
  chk::chk_true(eps_node > 0)
  chk::chk_number(iters_node)
  chk::chk_true(iters_node >= 1)
  chk::chk_number(tol_node)
  chk::chk_true(tol_node > 0)

  chk::chk_number(anchor_weight)
  chk::chk_true(anchor_weight >= 0)
  chk::chk_number(anchor_penalty)
  chk::chk_true(anchor_penalty >= 0)

  chk::chk_number(seed)
  chk::chk_logical(verbose)

  structure(
    list(
      n_levels = as.integer(n_levels),
      coarsen_ratio = as.numeric(coarsen_ratio),
      coarsen_method = coarsen_method,
      min_clusters = as.integer(min_clusters),
      prior_mode = prior_mode,
      prior_cluster_k = as.integer(prior_cluster_k),
      prior_max_candidates = as.integer(prior_max_candidates),
      graph_knn = as.integer(graph_knn),
      graph_sigma = if (is.null(graph_sigma)) NULL else as.numeric(graph_sigma),
      use_laplacian = isTRUE(use_laplacian),
      views = as.character(views),
      view_weights = if (is.null(view_weights)) NULL else as.numeric(view_weights),
      k_embed = if (is.null(k_embed)) NULL else as.integer(k_embed),
      hks_q = as.integer(hks_q),
      lambda_view = as.numeric(lambda_view),
      candidate_k = as.integer(candidate_k),
      ann_backend = ann_backend,
      ann_trees = as.integer(ann_trees),
      prior_k = as.integer(prior_k),
      ensure_cols = isTRUE(ensure_cols),
      token_mode = token_mode,
      max_hops = as.integer(max_hops),
      samples_per_hop = as.integer(samples_per_hop),
      max_tokens = as.integer(max_tokens),
      token_metric = token_metric,
      eps_token = as.numeric(eps_token),
      iters_token = as.integer(iters_token),
      tol_token = as.numeric(tol_token),
      eps_node = as.numeric(eps_node),
      iters_node = as.integer(iters_node),
      tol_node = as.numeric(tol_node),
      projection = projection,
      anchor_weight = as.numeric(anchor_weight),
      anchor_penalty = as.numeric(anchor_penalty),
      seed = as.integer(seed),
      verbose = isTRUE(verbose)
    ),
    class = "token_ot_graph_align_control"
  )
}

resolve_token_ot_graph_align_control <- function(control) {
  defaults <- token_ot_graph_align_control()
  if (missing(control) || is.null(control)) return(defaults)
  if (inherits(control, "token_ot_graph_align_control")) {
    merged <- modifyList(defaults, control)
    return(do.call(token_ot_graph_align_control, merged))
  }
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or token_ot_graph_align_control().", call. = FALSE)
  }
  if (is.null(names(control)) || any(names(control) == "")) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control argument(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  merged <- modifyList(defaults, control)
  do.call(token_ot_graph_align_control, merged)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.normalize_rows_dense <- function(M) {
  rs <- rowSums(M)
  rs[!is.finite(rs) | rs == 0] <- 1
  M / rs
}

.normalize_tokens <- function(T) {
  T <- as.matrix(T)
  nrm <- sqrt(rowSums(T * T))
  nrm[!is.finite(nrm) | nrm < 1e-12] <- 1
  T / nrm
}

.view_weights_vec <- function(views, view_weights) {
  if (is.null(view_weights)) {
    w <- rep(1, length(views))
  } else {
    w <- as.numeric(view_weights)
  }
  if (sum(w) <= 0) w <- rep(1, length(views))
  w / sum(w)
}

.scale_columns <- function(M) {
  M <- as.matrix(M)
  s <- apply(M, 2, stats::sd)
  s[!is.finite(s) | s < 1e-12] <- 1
  sweep(M, 2, s, "/")
}

.build_knn_graph <- function(X, knn, sigma = NULL, seed = 1L) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 3L) stop("Each domain must have at least 3 rows.", call. = FALSE)
  k_use <- min(as.integer(knn), n - 1L)
  sigma_use <- if (is.null(sigma)) {
    tryCatch({
      old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
      on.exit({
        if (is.null(old)) {
          if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
        } else {
          assign(".Random.seed", old, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(as.integer(seed))
      choose_sigma(X)
    }, error = function(e) 0.5)
  } else {
    as.numeric(sigma)
  }
  gw <- safe_compute(
    neighborweights::graph_weights(X, k = k_use, weight_mode = "heat",
                                   sigma = sigma_use, neighbor_mode = "knn"),
    "Graph construction failed in token_ot_graph_align"
  )
  A <- neighborweights::adjacency(gw)
  A <- (A + Matrix::t(A)) / 2
  Matrix::diag(A) <- 0
  A <- Matrix::drop0(methods::as(A, "dgCMatrix"))
  list(A = A, sigma = sigma_use, knn = k_use)
}

.neighbors_from_sparse <- function(A) {
  n <- nrow(A)
  s <- Matrix::summary(A)
  if (!nrow(s)) return(replicate(n, integer(0), simplify = FALSE))
  split_j <- split(s$j, s$i)
  out <- vector("list", n)
  for (i in seq_len(n)) {
    key <- as.character(i)
    out[[i]] <- if (!is.null(split_j[[key]])) unique(as.integer(split_j[[key]])) else integer(0)
  }
  out
}

.compute_views <- function(X, A, ncomp, ctrl, seed) {
  views <- ctrl$views
  Z <- vector("list", length(views))
  names(Z) <- views

  if ("raw" %in% views) {
    Xc <- scale(as.matrix(X), center = TRUE, scale = FALSE)
    k_raw <- min(as.integer(ncomp), ncol(Xc), nrow(Xc))
    if (k_raw < 1L) stop("ncomp too small for raw view.", call. = FALSE)
    if (ncol(Xc) == k_raw) {
      Z[["raw"]] <- Xc
    } else {
      pc <- tryCatch({
        if (requireNamespace("irlba", quietly = TRUE) && k_raw < min(nrow(Xc), ncol(Xc))) {
          irlba::irlba(Xc, nv = k_raw, nu = 0)
        } else {
          stats::prcomp(Xc, rank. = k_raw)
        }
      }, error = function(e) stats::prcomp(Xc, rank. = k_raw))
      scores <- if (!is.null(pc$x)) pc$x else Xc %*% pc$v
      Z[["raw"]] <- as.matrix(scores[, seq_len(k_raw), drop = FALSE])
    }
  }

  if (any(c("eigenmap", "hks") %in% views)) {
    k_embed <- ctrl$k_embed
    if (is.null(k_embed)) {
      k_embed <- max(3L * as.integer(ncomp), as.integer(ncomp) + 10L)
    }
    k_embed <- min(as.integer(k_embed), nrow(A) - 1L)
    basis <- compute_laplacian_basis(A, k = k_embed, normalized = ctrl$use_laplacian, seed = as.integer(seed))
    k_eig <- min(as.integer(ncomp), ncol(basis$vectors))
    if ("eigenmap" %in% views) {
      Z[["eigenmap"]] <- as.matrix(basis$vectors[, seq_len(k_eig), drop = FALSE])
    }
    if ("hks" %in% views) {
      Z[["hks"]] <- compute_hks_descriptors(
        basis = basis,
        q = ctrl$hks_q,
        time_mode = "auto",
        spacing = "log",
        time_scale = 1,
        normalize = "both"
      )
    }
  }

  Z
}

.fused_embedding <- function(Z_views, ctrl) {
  views <- ctrl$views
  w <- .view_weights_vec(views, ctrl$view_weights)
  mats <- lapply(seq_along(views), function(k) {
    v <- views[[k]]
    M <- .scale_columns(Z_views[[v]])
    sqrt(w[[k]]) * M
  })
  do.call(cbind, mats)
}

.token_bag_for_node <- function(idx, fused, neighbors, degree, ctrl) {
  d <- ncol(fused)
  hop_dim <- as.integer(ctrl$max_hops) + 1L
  deg_feat <- if (!is.null(degree)) 1L else 0L

  make_token <- function(node, hop) {
    hop_vec <- rep(0, hop_dim)
    hop_vec[hop + 1L] <- 1
    base <- fused[node, , drop = TRUE]
    if (deg_feat) {
      c(base, hop_vec, degree[node])
    } else {
      c(base, hop_vec)
    }
  }

  tokens <- list(make_token(idx, hop = 0L))
  if (identical(ctrl$token_mode, "view_only")) {
    out <- do.call(rbind, tokens)
    return(.normalize_tokens(out))
  }

  max_hops <- as.integer(ctrl$max_hops)
  if (max_hops < 1L) {
    out <- do.call(rbind, tokens)
    return(.normalize_tokens(out))
  }

  seen <- rep(FALSE, nrow(fused))
  seen[idx] <- TRUE
  frontier <- idx

  for (hop in seq_len(max_hops)) {
    # Expand frontier
    nbrs <- unique(unlist(neighbors[frontier], use.names = FALSE))
    nbrs <- nbrs[!seen[nbrs]]
    if (!length(nbrs)) break
    seen[nbrs] <- TRUE

    take <- nbrs
    if (length(take) > as.integer(ctrl$samples_per_hop)) {
      # Deterministic cap for reproducibility (avoid RNG breaking equivariance).
      take <- sort(take)[seq_len(as.integer(ctrl$samples_per_hop))]
    }
    for (u in take) {
      tokens[[length(tokens) + 1L]] <- make_token(u, hop = hop)
    }
    frontier <- nbrs
    if (length(tokens) >= as.integer(ctrl$max_tokens)) break
  }

  if (length(tokens) > as.integer(ctrl$max_tokens)) {
    tokens <- tokens[seq_len(as.integer(ctrl$max_tokens))]
  }
  out <- do.call(rbind, tokens)
  .normalize_tokens(out)
}

.token_ot_cost <- function(tokensA, tokensB, ctrl) {
  # Tokens are already normalized when produced by .token_bag_for_node().
  XA <- as.matrix(tokensA)
  YB <- as.matrix(tokensB)
  p <- nrow(XA)
  q <- nrow(YB)
  if (p < 1L || q < 1L) stop("Empty token bag.", call. = FALSE)

  if (identical(ctrl$token_metric, "cosine")) {
    sim <- XA %*% t(YB)
    sim <- pmin(pmax(sim, -1), 1)
    D <- 1 - sim
  } else {
    xa2 <- rowSums(XA^2)
    yb2 <- rowSums(YB^2)
    D <- outer(xa2, yb2, "+") - 2 * (XA %*% t(YB))
    D[D < 0 & D > -1e-12] <- 0
  }

  # If either bag has size 1, OT constraints uniquely determine the plan
  # (independent of epsilon), so skip Sinkhorn for speed and determinism.
  if (p == 1L) return(mean(D))
  if (q == 1L) return(mean(D))

  a <- rep(1 / p, p)
  b <- rep(1 / q, q)

  pi <- sinkhorn_unified(
    D,
    a,
    b,
    epsilon = ctrl$eps_token,
    max_iter = ctrl$iters_token,
    tol = ctrl$tol_token,
    stabilized = TRUE
  )
  sum(pi * D)
}

.node_pair_cost <- function(i, j, Zs_views, Zt_views, Ms, Mt, ctrl, anchors = NULL) {
  views <- ctrl$views
  w <- .view_weights_vec(views, ctrl$view_weights)

  base <- 0
  for (k in seq_along(views)) {
    v <- views[[k]]
    xs <- Zs_views[[v]][i, , drop = TRUE]
    yt <- Zt_views[[v]][j, , drop = TRUE]
    sim <- sum(xs * yt)
    sim <- max(min(sim, 1), -1)
    base <- base + w[[k]] * (1 - sim)
  }

  ot_cost <- .token_ot_cost(Ms[[i]], Mt[[j]], ctrl)

  cost <- ctrl$lambda_view * base + (1 - ctrl$lambda_view) * ot_cost

  if (!is.null(anchors) && ctrl$anchor_weight > 0) {
    a1 <- anchors$vec1[i]
    a2 <- anchors$vec2[j]
    if (!is.na(a1) && !is.na(a2) && a1 != a2) {
      cost <- cost + ctrl$anchor_weight * ctrl$anchor_penalty
    }
  }

  cost
}

.pair_costs_view_only <- function(i_idx, j_idx, Zs_views, Zt_views, Ts, Tt, ctrl, anchors = NULL, chunk_size = 50000L) {
  i_idx <- as.integer(i_idx)
  j_idx <- as.integer(j_idx)
  n_edges <- length(i_idx)
  if (n_edges != length(j_idx)) stop("edge index vectors must have equal length.", call. = FALSE)
  if (!n_edges) return(numeric(0))

  views <- ctrl$views
  w <- .view_weights_vec(views, ctrl$view_weights)
  lam <- ctrl$lambda_view

  chunk_size <- max(1L, as.integer(chunk_size))
  out <- numeric(n_edges)

  for (start in seq.int(1L, n_edges, by = chunk_size)) {
    end <- min(n_edges, start + chunk_size - 1L)
    idx <- start:end
    ii <- i_idx[idx]
    jj <- j_idx[idx]

    base <- numeric(length(idx))
    for (k in seq_along(views)) {
      v <- views[[k]]
      xs <- Zs_views[[v]][ii, , drop = FALSE]
      yt <- Zt_views[[v]][jj, , drop = FALSE]
      sim <- rowSums(xs * yt)
      sim <- pmin(pmax(sim, -1), 1)
      base <- base + w[[k]] * (1 - sim)
    }

    xs_tok <- Ts[ii, , drop = FALSE]
    yt_tok <- Tt[jj, , drop = FALSE]
    if (identical(ctrl$token_metric, "cosine")) {
      sim_tok <- rowSums(xs_tok * yt_tok)
      sim_tok <- pmin(pmax(sim_tok, -1), 1)
      tok <- 1 - sim_tok
    } else {
      dif <- xs_tok - yt_tok
      tok <- rowSums(dif * dif)
    }

    cost <- lam * base + (1 - lam) * tok

    if (!is.null(anchors) && ctrl$anchor_weight > 0) {
      a1 <- anchors$vec1[ii]
      a2 <- anchors$vec2[jj]
      mismatch <- !is.na(a1) & !is.na(a2) & a1 != a2
      if (any(mismatch)) {
        cost[mismatch] <- cost[mismatch] + ctrl$anchor_weight * ctrl$anchor_penalty
      }
    }

    out[idx] <- cost
  }
  out
}

.ann_candidates <- function(src, tgt, k = 50L, ctrl = NULL) {
  src <- as.matrix(src)
  tgt <- as.matrix(tgt)
  k <- min(as.integer(k), nrow(tgt))
  # Use Euclidean kNN on normalized vectors for cosine-like behaviour
  srcn <- .normalize_tokens(src)
  tgtn <- .normalize_tokens(tgt)
  backend <- if (is.null(ctrl)) "rann" else as.character(ctrl$ann_backend %||% "auto")
  if (backend == "auto") {
    backend <- if (requireNamespace("RcppAnnoy", quietly = TRUE) && nrow(tgtn) >= 2000L) "annoy" else "rann"
  }

  if (backend == "annoy") {
    if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
      stop("RcppAnnoy is required for ann_backend='annoy'.", call. = FALSE)
    }
    ann_ctor <- get0("AnnoyAngular", envir = asNamespace("RcppAnnoy"))
    if (is.null(ann_ctor) || !inherits(ann_ctor, "C++Class")) {
      stop("RcppAnnoy::AnnoyAngular constructor not found; please reinstall RcppAnnoy.", call. = FALSE)
    }
    ann <- methods::new(ann_ctor, ncol(tgtn))
    for (j in seq_len(nrow(tgtn))) {
      ann$addItem(j - 1L, tgtn[j, ])
    }
    n_trees <- if (is.null(ctrl)) 50L else as.integer(ctrl$ann_trees %||% 50L)
    ann$build(n_trees)

    idx <- matrix(NA_integer_, nrow(srcn), k)
    for (i in seq_len(nrow(srcn))) {
      idx[i, ] <- ann$getNNsByVector(srcn[i, ], k) + 1L
    }
    return(idx)
  }

  nn <- RANN::nn2(data = tgtn, query = srcn, k = k)
  nn$nn.idx
}

.ensure_target_columns <- function(Cand, src_fused, tgt_fused) {
  n <- nrow(src_fused)
  m <- nrow(tgt_fused)
  incoming <- integer(m)
  for (i in seq_len(n)) {
    jj <- Cand[[i]]
    if (length(jj)) incoming[jj] <- incoming[jj] + 1L
  }
  missing <- which(incoming == 0L)
  if (!length(missing)) return(Cand)

  # Add a reverse 1NN edge for each missing target
  srcn <- .normalize_tokens(src_fused)
  tgtn <- .normalize_tokens(tgt_fused)
  nn <- RANN::nn2(data = srcn, query = tgtn[missing, , drop = FALSE], k = 1)
  for (k in seq_along(missing)) {
    j <- missing[k]
    i <- nn$nn.idx[k, 1]
    Cand[[i]] <- unique(c(Cand[[i]], j))
  }
  Cand
}

.sparse_sinkhorn_items <- function(i_idx, j_idx, cost, n, m, eps, iters, tol) {
  i_idx <- as.integer(i_idx)
  j_idx <- as.integer(j_idx)
  cost <- as.numeric(cost)
  if (length(i_idx) != length(j_idx) || length(i_idx) != length(cost)) {
    stop("Sparse Sinkhorn: edge arrays must have the same length.", call. = FALSE)
  }
  if (!length(i_idx)) stop("Sparse Sinkhorn: no candidate edges.", call. = FALSE)

  row_counts <- tabulate(i_idx, nbins = as.integer(n))
  col_counts <- tabulate(j_idx, nbins = as.integer(m))
  if (any(row_counts == 0L)) {
    stop("At least one source node has zero candidates; increase candidate_k or adjust constraints.", call. = FALSE)
  }
  if (any(col_counts == 0L)) {
    stop("At least one target node has zero incoming candidates; set ensure_cols=TRUE or increase candidate_k.", call. = FALSE)
  }

  a <- rep(1 / n, n)
  b <- rep(1 / m, m)

  # Cost scaling for numeric stability
  c0 <- min(cost)
  if (is.finite(c0) && c0 != 0) cost <- cost - c0
  pos <- cost[cost > 0]
  s <- if (length(pos)) stats::median(pos, na.rm = TRUE) else 1
  if (!is.finite(s) || s <= 0) s <- 1
  cost_s <- cost / s

  K <- exp(-cost_s / eps)
  K[!is.finite(K)] <- 0

  Ksp <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = K, dims = c(n, m))
  Ksp_t <- Matrix::t(Ksp)

  u <- rep(1, n)
  v <- rep(1, m)
  tiny <- 1e-300
  for (it in seq_len(iters)) {
    u_old <- u
    denom_u <- as.numeric(Ksp %*% v)
    denom_u[!is.finite(denom_u) | denom_u <= 0] <- tiny
    u <- a / denom_u

    denom_v <- as.numeric(Ksp_t %*% u)
    denom_v[!is.finite(denom_v) | denom_v <= 0] <- tiny
    v <- b / denom_v

    if (max(abs(u - u_old)) < tol) break
  }

  w <- u[i_idx] * K * v[j_idx]
  w[!is.finite(w)] <- 0
  list(
    i = i_idx,
    j = j_idx,
    w = w,
    dims = c(n, m),
    cost_shift = c0,
    cost_scale = s
  )
}

.greedy_projection <- function(edge_i, edge_j, edge_w, n, m) {
  ord <- order(edge_w, decreasing = TRUE, na.last = NA)
  assign <- rep(NA_integer_, n)
  used <- rep(FALSE, m)
  for (k in ord) {
    i <- edge_i[[k]]
    j <- edge_j[[k]]
    if (!is.na(assign[i]) || used[j]) next
    assign[i] <- j
    used[j] <- TRUE
  }
  assign
}

.hungarian_projection <- function(edge_i, edge_j, edge_w, n, m) {
  # Maximize weights via LSAP on costs = -log(w).
  # Fill missing edges with large cost.
  tiny <- 1e-300
  C <- matrix(1e6, n, m)
  for (k in seq_along(edge_w)) {
    C[edge_i[[k]], edge_j[[k]]] <- -log(edge_w[[k]] + tiny)
  }
  solve_lsap_with_padding(C, method = "dense")
}

.kmeans_membership <- function(X, k, seed = 1L) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (k >= n) {
    return(seq_len(n))
  }
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old, envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  km <- stats::kmeans(scale(X, center = TRUE, scale = FALSE), centers = as.integer(k), nstart = 1, iter.max = 50)
  as.integer(km$cluster)
}

.centroids_from_membership <- function(X, membership) {
  X <- as.matrix(X)
  g <- as.integer(membership)
  k <- max(g)
  counts <- as.numeric(tabulate(g, nbins = k))
  sums <- rowsum(X, group = g, reorder = FALSE)
  sweep(sums, 1, pmax(counts, 1), "/")
}

.build_hierarchy <- function(X, n_levels, ratio, min_clusters, seed = 1L,
                             method = c("kmeans", "louvain"),
                             graph_knn = 15L,
                             graph_sigma = NULL) {
  method <- match.arg(method)
  X_levels <- list(as.matrix(X))
  parent <- list()
  sizes <- integer(0)
  sizes[1] <- nrow(X_levels[[1]])

  for (lvl in seq_len(as.integer(n_levels) - 1L)) {
    X_prev <- X_levels[[lvl]]
    n_prev <- nrow(X_prev)
    if (n_prev <= as.integer(min_clusters)) break
    k_next <- max(as.integer(min_clusters), ceiling(n_prev / as.numeric(ratio)))
    if (k_next >= n_prev) break
    mem <- NULL
    if (identical(method, "louvain")) {
      g <- tryCatch(
        .build_knn_graph(X_prev, knn = as.integer(graph_knn), sigma = graph_sigma, seed = as.integer(seed) + lvl * 1009L),
        error = function(e) NULL
      )
      if (!is.null(g)) {
        ig <- igraph::graph_from_adjacency_matrix(g$A, mode = "undirected", weighted = TRUE, diag = FALSE)
        cl <- tryCatch(igraph::cluster_louvain(ig, weights = igraph::E(ig)$weight), error = function(e) NULL)
        if (!is.null(cl)) {
          mem0 <- igraph::membership(cl)
          mem1 <- as.integer(factor(mem0))
          k_mem <- max(mem1)
          if (is.finite(k_mem) && k_mem >= as.integer(min_clusters) && k_mem < n_prev) {
            # Accept if it's in a plausible range for the requested ratio.
            if (k_mem <= max(ceiling(k_next * 4), ceiling(n_prev / 2))) {
              mem <- mem1
            }
          }
        }
      }
    }
    if (is.null(mem)) {
      mem <- .kmeans_membership(X_prev, k = k_next, seed = as.integer(seed) + lvl * 1009L)
    }
    X_next <- .centroids_from_membership(X_prev, mem)
    parent[[lvl]] <- mem
    X_levels[[lvl + 1L]] <- X_next
    sizes[lvl + 1L] <- nrow(X_next)
  }

  list(
    X_levels = X_levels,
    parent = parent,
    sizes = sizes,
    n_levels = length(X_levels)
  )
}

.prior_candidates_from_coarse <- function(fit_coarse, parent_s, parent_t, ctrl) {
  mode <- ctrl$prior_mode
  if (is.null(mode) || identical(mode, "none")) return(NULL)
  parent_s <- as.integer(parent_s)
  parent_t <- as.integer(parent_t)
  n_fine_s <- length(parent_s)

  # Target children per coarse node
  split_child <- split(seq_along(parent_t), parent_t)
  n_coarse_t <- max(parent_t)
  children_t <- vector("list", n_coarse_t)
  for (ct in seq_len(n_coarse_t)) {
    key <- as.character(ct)
    children_t[[ct]] <- if (!is.null(split_child[[key]])) as.integer(split_child[[key]]) else integer(0)
  }

  # Coarse map: for each coarse source node, list of coarse target nodes
  n_coarse_s <- nrow(fit_coarse$transport)
  coarse_map <- vector("list", n_coarse_s)
  if (identical(mode, "hard")) {
    perm <- as.integer(fit_coarse$assignment)
    for (cs in seq_len(n_coarse_s)) {
      ct <- perm[[cs]]
      coarse_map[[cs]] <- if (!is.na(ct) && ct >= 1L && ct <= n_coarse_t) as.integer(ct) else integer(0)
    }
  } else {
    # soft: top-k target clusters by coupling mass per coarse source
    W <- fit_coarse$transport
    sm <- Matrix::summary(W)
    row_idx <- split(seq_len(nrow(sm)), sm$i)
    k_keep <- as.integer(ctrl$prior_cluster_k)
    for (cs in seq_len(n_coarse_s)) {
      key <- as.character(cs)
      rr <- row_idx[[key]]
      if (is.null(rr) || !length(rr)) {
        coarse_map[[cs]] <- integer(0)
        next
      }
      jj <- sm$j[rr]
      xx <- sm$x[rr]
      ord <- order(xx, decreasing = TRUE)
      jj <- jj[ord]
      jj <- jj[seq_len(min(k_keep, length(jj)))]
      coarse_map[[cs]] <- unique(as.integer(jj))
    }
  }

  # Lift to fine candidates
  prior <- vector("list", n_fine_s)
  max_keep <- if (as.integer(ctrl$prior_k) > 0L) as.integer(ctrl$prior_k) else as.integer(ctrl$prior_max_candidates)
  for (i in seq_len(n_fine_s)) {
    cs <- parent_s[[i]]
    cts <- coarse_map[[cs]]
    if (!length(cts)) {
      prior[[i]] <- integer(0)
      next
    }
    cand <- unique(unlist(children_t[cts], use.names = FALSE))
    if (length(cand) > max_keep) {
      cand <- sort(cand)[seq_len(max_keep)]
    }
    prior[[i]] <- as.integer(cand)
  }
  prior
}

.token_ot_graph_align_level <- function(Xs, Xt, anchors = NULL, ncomp, ctrl, priorCand = NULL, level = 1L) {
  verbose <- isTRUE(ctrl$verbose)
  seed <- as.integer(ctrl$seed)

  if (verbose) message("token_ot_graph_align: building within-domain graphs")
  gs <- .build_knn_graph(Xs, knn = ctrl$graph_knn, sigma = ctrl$graph_sigma, seed = seed + 101)
  gt <- .build_knn_graph(Xt, knn = ctrl$graph_knn, sigma = ctrl$graph_sigma, seed = seed + 202)

  As <- gs$A
  At <- gt$A
  deg_s <- as.numeric(Matrix::rowSums(As))
  deg_t <- as.numeric(Matrix::rowSums(At))
  deg_s <- deg_s / max(deg_s + 1e-12)
  deg_t <- deg_t / max(deg_t + 1e-12)

  if (verbose) message("token_ot_graph_align: computing views")
  Zs <- .compute_views(Xs, As, ncomp = ncomp, ctrl = ctrl, seed = seed + 11)
  Zt <- .compute_views(Xt, At, ncomp = ncomp, ctrl = ctrl, seed = seed + 17)

  # Pre-normalize view embeddings for fast cosine distances in .node_pair_cost
  Zs_norm <- lapply(Zs, .normalize_tokens)
  Zt_norm <- lapply(Zt, .normalize_tokens)

  fused_s <- .fused_embedding(Zs, ctrl)
  fused_t <- .fused_embedding(Zt, ctrl)

  n <- nrow(fused_s)
  m <- nrow(fused_t)

  Ms <- NULL
  Mt <- NULL
  Ts <- NULL
  Tt <- NULL
  if (verbose) message("token_ot_graph_align: building token representations")
  if (identical(ctrl$token_mode, "view_only")) {
    hop_dim <- as.integer(ctrl$max_hops) + 1L
    hop0 <- c(1, rep(0, max(0L, hop_dim - 1L)))
    Hs <- matrix(rep(hop0, n), nrow = n, ncol = hop_dim, byrow = TRUE)
    Ht <- matrix(rep(hop0, m), nrow = m, ncol = hop_dim, byrow = TRUE)
    Ts <- .normalize_tokens(cbind(fused_s, Hs, deg_s))
    Tt <- .normalize_tokens(cbind(fused_t, Ht, deg_t))
  } else {
    nbr_s <- .neighbors_from_sparse(As)
    nbr_t <- .neighbors_from_sparse(At)
    Ms <- lapply(seq_len(n), function(i) .token_bag_for_node(i, fused_s, nbr_s, degree = deg_s, ctrl = ctrl))
    Mt <- lapply(seq_len(m), function(j) .token_bag_for_node(j, fused_t, nbr_t, degree = deg_t, ctrl = ctrl))
  }

  if (verbose) message("token_ot_graph_align: generating candidates (k=", ctrl$candidate_k, ")")
  nn_idx <- .ann_candidates(fused_s, fused_t, k = ctrl$candidate_k, ctrl = ctrl)
  Cand <- lapply(seq_len(n), function(i) as.integer(nn_idx[i, ]))

  # Add anchor-consistent candidates if anchors provided
  if (!is.null(anchors)) {
    vec1 <- anchors$vec1
    vec2 <- anchors$vec2
    for (i in seq_len(n)) {
      lab <- vec1[[i]]
      if (is.na(lab)) next
      idx <- which(vec2 == lab)
      if (length(idx)) Cand[[i]] <- unique(c(Cand[[i]], idx))
    }
  }

  # Add prior candidates from coarser level
  if (!is.null(priorCand)) {
    if (length(priorCand) != n) {
      stop("priorCand must have length equal to the number of source nodes at this level.", call. = FALSE)
    }
    max_keep <- if (as.integer(ctrl$prior_k) > 0L) as.integer(ctrl$prior_k) else as.integer(ctrl$prior_max_candidates)
    for (i in seq_len(n)) {
      pc <- as.integer(priorCand[[i]])
      if (length(pc) > max_keep) pc <- sort(pc)[seq_len(max_keep)]
      if (length(pc)) Cand[[i]] <- unique(c(Cand[[i]], pc))
    }
  }

  if (isTRUE(ctrl$ensure_cols)) {
    Cand <- .ensure_target_columns(Cand, fused_s, fused_t)
  }

  if (verbose) message("token_ot_graph_align: computing sparse costs")
  edge_counts <- vapply(Cand, length, integer(1))
  n_edges <- sum(edge_counts)
  i_idx <- integer(n_edges)
  j_idx <- integer(n_edges)
  ptr <- 0L
  for (i in seq_len(n)) {
    js <- Cand[[i]]
    if (!length(js)) next
    for (j in js) {
      ptr <- ptr + 1L
      i_idx[[ptr]] <- i
      j_idx[[ptr]] <- j
    }
  }
  if (ptr < n_edges) {
    i_idx <- i_idx[seq_len(ptr)]
    j_idx <- j_idx[seq_len(ptr)]
  }
  c_vec <- if (identical(ctrl$token_mode, "view_only")) {
    .pair_costs_view_only(i_idx, j_idx, Zs_norm, Zt_norm, Ts = Ts, Tt = Tt, ctrl = ctrl, anchors = anchors)
  } else {
    vapply(seq_along(i_idx), function(k) {
      .node_pair_cost(i_idx[[k]], j_idx[[k]], Zs_norm, Zt_norm, Ms, Mt, ctrl, anchors = anchors)
    }, numeric(1))
  }

  if (verbose) message("token_ot_graph_align: sparse Sinkhorn on ", length(c_vec), " edges")
  sol <- .sparse_sinkhorn_items(i_idx, j_idx, c_vec, n = n, m = m,
                                eps = ctrl$eps_node, iters = ctrl$iters_node, tol = ctrl$tol_node)

  cost_scaled <- (c_vec - sol$cost_shift) / sol$cost_scale
  objective_scaled <- sum(sol$w * cost_scaled)
  objective <- sol$cost_scale * objective_scaled + sol$cost_shift

  if (verbose) message("token_ot_graph_align: projecting to discrete assignment (", ctrl$projection, ")")
  perm <- if (identical(ctrl$projection, "hungarian") && max(n, m) <= 2000) {
    .hungarian_projection(sol$i, sol$j, sol$w, n = n, m = m)
  } else {
    .greedy_projection(sol$i, sol$j, sol$w, n = n, m = m)
  }

  W <- Matrix::sparseMatrix(i = sol$i, j = sol$j, x = sol$w, dims = c(n, m))

  list(
    transport = W,
    assignment = perm,
    views_s = Zs,
    views_t = Zt,
    fused_s = fused_s,
    fused_t = fused_t,
    graphs = list(source = As, target = At),
    diagnostics = list(
      nnz = length(sol$w),
      n = n, m = m,
      eps_token = ctrl$eps_token,
      eps_node = ctrl$eps_node,
      objective_scaled = objective_scaled,
      objective = objective,
      level = as.integer(level)
    )
  )
}

.token_ot_graph_pair_fit <- function(Xs, Xt, anchors = NULL, ncomp, ctrl) {
  L <- as.integer(ctrl$n_levels %||% 1L)
  if (L <= 1L) {
    return(.token_ot_graph_align_level(Xs, Xt, anchors = anchors, ncomp = ncomp, ctrl = ctrl, priorCand = NULL, level = 1L))
  }

  hs <- .build_hierarchy(
    Xs,
    n_levels = L,
    ratio = ctrl$coarsen_ratio,
    min_clusters = ctrl$min_clusters,
    seed = ctrl$seed + 5001L,
    method = ctrl$coarsen_method %||% "kmeans",
    graph_knn = ctrl$graph_knn,
    graph_sigma = ctrl$graph_sigma
  )
  ht <- .build_hierarchy(
    Xt,
    n_levels = L,
    ratio = ctrl$coarsen_ratio,
    min_clusters = ctrl$min_clusters,
    seed = ctrl$seed + 7001L,
    method = ctrl$coarsen_method %||% "kmeans",
    graph_knn = ctrl$graph_knn,
    graph_sigma = ctrl$graph_sigma
  )
  L_eff <- min(hs$n_levels, ht$n_levels)
  if (L_eff < L && isTRUE(ctrl$verbose)) {
    message("token_ot_graph_align: reduced n_levels from ", L, " to ", L_eff, " due to problem size")
  }

  prior <- NULL
  fits <- vector("list", L_eff)
  for (lvl in seq(from = L_eff, to = 1L, by = -1L)) {
    Xs_l <- hs$X_levels[[lvl]]
    Xt_l <- ht$X_levels[[lvl]]
    # Anchors only at the finest level by default.
    anchors_l <- if (lvl == 1L) anchors else NULL
    ncomp_l <- min(as.integer(ncomp), nrow(Xs_l), nrow(Xt_l))
    fits[[lvl]] <- .token_ot_graph_align_level(
      Xs = Xs_l, Xt = Xt_l,
      anchors = anchors_l,
      ncomp = ncomp_l,
      ctrl = ctrl,
      priorCand = prior,
      level = lvl
    )

    if (lvl > 1L) {
      parent_s <- hs$parent[[lvl - 1L]]
      parent_t <- ht$parent[[lvl - 1L]]
      prior <- .prior_candidates_from_coarse(fits[[lvl]], parent_s = parent_s, parent_t = parent_t, ctrl = ctrl)
    }
  }

  out <- fits[[1L]]
  out$multilevel <- list(
    n_levels = L_eff,
    source_sizes = hs$sizes,
    target_sizes = ht$sizes,
    fits = fits
  )
  out
}

#' @rdname token_ot_graph_align
#' @method token_ot_graph_align hyperdesign
#' @export
token_ot_graph_align.hyperdesign <- function(
  data,
  anchors = NULL,
  preproc = multivarious::center(),
  ncomp = 20L,
  control = token_ot_graph_align_control(),
  ...
) {
  ctrl <- resolve_token_ot_graph_align_control(control)
  verbose <- isTRUE(ctrl$verbose)

  chk::chk_s3_class(data, "hyperdesign")
  resolved <- resolve_hyperdesign(data)
  domains <- resolved$domains
  domain_names <- resolved$domain_names
  if (length(domains) != 2L) {
    stop("token_ot_graph_align currently supports exactly 2 domains.", call. = FALSE)
  }

  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }
  ncomp <- as.integer(ncomp)

  if (verbose) message("token_ot_graph_align: preprocessing domains")
  pre <- .apply_preproc_domains_spectral_mnn(domains, preproc = preproc)
  pdata <- pre$pdata
  proc <- pre$proc
  names(pdata) <- domain_names

  # Extract anchors (optional)
  anchor_info <- NULL
  anchors_quo <- rlang::enquo(anchors)
  if (!rlang::quo_is_null(anchors_quo)) {
    anchor_name <- rlang::as_name(rlang::ensym(anchors))
    v1 <- pdata[[1]]$design[[anchor_name]]
    v2 <- pdata[[2]]$design[[anchor_name]]
    if (is.null(v1) || is.null(v2)) {
      stop("Anchor column `", anchor_name, "` not found in both domain designs.", call. = FALSE)
    }
    if (length(v1) != nrow(pdata[[1]]$x) || length(v2) != nrow(pdata[[2]]$x)) {
      stop("Anchor column length must match nrow(x) for each domain.", call. = FALSE)
    }
    anchor_info <- list(vec1 = v1, vec2 = v2)
  }

  fit <- .token_ot_graph_pair_fit(
    Xs = pdata[[1]]$x,
    Xt = pdata[[2]]$x,
    anchors = anchor_info,
    ncomp = ncomp,
    ctrl = ctrl
  )

  # Store baseline embeddings as scores (for compatibility)
  E1 <- fit$views_s[["eigenmap"]] %||% fit$views_s[[ctrl$views[[1]]]]
  E2 <- fit$views_t[["eigenmap"]] %||% fit$views_t[[ctrl$views[[1]]]]
  k_use <- min(ncomp, ncol(E1), ncol(E2))
  scores <- rbind(as.matrix(E1[, seq_len(k_use), drop = FALSE]),
                  as.matrix(E2[, seq_len(k_use), drop = FALSE]))

  loadings <- do.call(rbind, lapply(seq_along(pdata), function(i) {
    Xi <- as.matrix(pdata[[i]]$x)
    Ei <- if (i == 1L) E1 else E2
    Matrix::crossprod(Xi, Ei[, seq_len(k_use), drop = FALSE])
  }))

  feature_blocks <- feature_block_indices(pdata)

  new_alignment_result(
    scores = scores,
    loadings = as.matrix(loadings),
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "token_ot_graph_align",
    extras = list(
      transport_plan = fit$transport,
      assignment = fit$assignment,
      anchors = anchor_info,
      control = ctrl,
      diagnostics = fit$diagnostics,
      graphs = fit$graphs,
      multilevel = fit$multilevel
    )
  )
}

#' @rdname token_ot_graph_align
#' @method token_ot_graph_align list
#' @export
token_ot_graph_align.list <- function(
  data,
  ...
) {
  if (!is.list(data) || length(data) != 2L) {
    stop("token_ot_graph_align.list expects a list of exactly two matrices.", call. = FALSE)
  }
  if (!all(vapply(data, function(x) is.matrix(x) || is.data.frame(x), logical(1)))) {
    stop("Each element of `data` must be a numeric matrix.", call. = FALSE)
  }
  hd <- as_hyperdesign(data)
  token_ot_graph_align(hd, ...)
}
