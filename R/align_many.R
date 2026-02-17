#' Align multiple datasets by synchronizing pairwise transforms
#'
#' This meta-aligner lifts any pairwise method to N domains by:
#' 1) choosing a sparse dataset graph, 2) fitting pairwise models on edges,
#' 3) extracting relative transforms, and 4) running group-appropriate
#' synchronization to get per-domain transforms into a shared space.
#'
#' If the adapter reports native multi-view support, dispatches to fit_many().
#'
#' @param domains list of domain objects (often matrices), one per dataset
#' @param algo an aligner object created by new_aligner() with adapter methods
#' @param graph "auto" | "mst" | "complete" | "star"
#' @param consensus "auto" | "sync" | "star" | "tree" (currently: "sync"/"star")
#' @param k target latent dimension (optional)
#' @param parallel logical; reserved for future parallel pair fits
#' @param sync list of synchronization options (graph, refine, gl, perm settings)
#' @param ... forwarded to adapter methods
#' @return list with transforms, embeddings (if apply_transform succeeds),
#'   dataset graph, pairwise fits, and diagnostics
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 10, 5)
#' X2 <- matrix(rnorm(50), 10, 5)
#' X3 <- matrix(rnorm(50), 10, 5)
#' domains <- list(X1, X2, X3)
#' # align_many requires an aligner object with adapter methods
#' # See package vignettes for complete examples
#' }
#' @export
align_many <- function(domains, algo,
                       graph = c("auto", "mst", "complete", "star"),
                       consensus = c("auto", "sync", "star", "tree"),
                       k = NULL, parallel = TRUE,
                       sync = NULL, ...) {

  graph <- match.arg(graph)
  consensus <- match.arg(consensus)

  caps <- aligner_capabilities(algo)
  if (isTRUE(caps$supports_multi)) {
    return(fit_many(algo, domains, k = k, ...))
  }

  N <- length(domains)
  if (N < 2) stop("align_many: need at least two domains")

  # Choose dataset graph
  if (graph == "auto") {
    graph <- if (N <= 5) "complete" else if (N <= 8) "knn" else "mst"
  }
  graph_opts <- tryCatch(sync$graph, error = function(e) NULL)
  E <- build_dataset_graph_(N, graph = graph, domains = domains, algo = algo, graph_opts = graph_opts, ...)

  # Pairwise fits
  fits <- fit_pairs_on_edges_(domains, E, algo, parallel = parallel, ...)

  # Relative transforms and weights
  G <- lapply(fits, function(f) relative_transform(f, from = "i", to = "j"))
  raw_w <- vapply(seq_along(fits), function(idx) {
    val <- tryCatch(pair_loss(fits[[idx]], NULL, NULL, ...), error = function(e) NA_real_)
    if (is.na(val) || !is.finite(val)) 1.0 else 1.0 / (1.0 + val)
  }, numeric(1))
  # Robust clipping of weights to reduce outliers' influence
  w_ok <- is.finite(raw_w)
  if (any(w_ok)) {
    q <- stats::quantile(raw_w[w_ok], probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE, type = 7)
    W <- pmax(pmin(raw_w, q[2L]), q[1L])
  } else {
    W <- raw_w
  }

  # Group-specific synchronization
  grp <- caps$group
  if (is.null(k)) {
    k <- tryCatch({
      ks <- vapply(G, function(g) {
        if (inherits(g, "align_transform") && !is.null(g$k)) return(g$k)
        op <- if (inherits(g, "align_transform")) g$op else as.matrix(g)
        ncol(op)
      }, numeric(1))
      if (length(ks)) min(ks) else NULL
    }, error = function(e) NULL)
  }

  T_list <- switch(grp,
    "O" = {
      T_tmp <- rotation_sync(G, W, E, N = N, k = k)
      # Optional one-pass local refinement on O(k)
      ref_opts <- tryCatch(sync$refine, error = function(e) NULL)
      if (is.list(ref_opts) && isTRUE(ref_opts$gpa) && length(G)) {
        T_tmp <- refine_orthogonal_local_(T_tmp, G, W, E, N = N)
      }
      T_tmp
    },
    "GL" = {
      gl_opts <- tryCatch(sync$gl, error = function(e) NULL)
      lambda <- if (is.list(gl_opts) && !is.null(gl_opts$ridge)) gl_opts$ridge else 1e-6
      iters  <- if (is.list(gl_opts) && !is.null(gl_opts$iters)) gl_opts$iters else 20
      gauge  <- if (is.list(gl_opts) && !is.null(gl_opts$gauge)) gl_opts$gauge else 1L
      project<- if (is.list(gl_opts) && !is.null(gl_opts$project)) gl_opts$project else "none"
      gl_sync_least_squares(G, W, E, N = N, k = k, lambda = lambda, iters = iters, gauge_idx = gauge, project = project)
    },
    "perm" = {
      perm_opts <- tryCatch(sync$perm, error = function(e) NULL)
      mode <- if (is.list(perm_opts) && !is.null(perm_opts$mode)) perm_opts$mode else if (N > 4) "spectral" else "star"
      if (mode == "star") {
        permutation_sync_star(G, W, E, N = N)
      } else if (mode == "global") {
        permutation_sync(G, W, E, N = N)
      } else { # spectral
        L <- if (is.list(perm_opts)) perm_opts$k else NULL
        rounding <- if (is.list(perm_opts) && !is.null(perm_opts$rounding)) perm_opts$rounding else "sinkhorn"
        max_sink <- if (is.list(perm_opts) && !is.null(perm_opts$maxit_sinkhorn)) perm_opts$maxit_sinkhorn else 200
        tol_s    <- if (is.list(perm_opts) && !is.null(perm_opts$tol)) perm_opts$tol else 1e-6
        permutation_sync_spectral(G, W, E, N = N, L = L, rounding = rounding, max_sinkhorn = max_sink, tol = tol_s)
      }
    },
    "kernel" = rotation_sync(G, W, E, N = N, k = k),
    stop("Unknown group: ", grp)
  )

  # Try embedding via apply_transform; if adapter uses different return type,
  # user can ignore 'embeddings'
  Z <- lapply(seq_len(N), function(i) {
    tryCatch(apply_transform(T_list[[i]], domains[[i]]), error = function(e) NULL)
  })

  diags <- compute_diagnostics_(grp, T_list, G, E, W)

  list(transforms = T_list,
       embeddings = Z,
       graph = E,
       pair_fits = fits,
       diagnostics = diags)
}

# -- Graph construction and pair fitting -----------------------------------

# Build a dataset graph as a data.frame with columns i, j (1-based indices)
# star: hub = 1; complete: all pairs; mst: TODO (falls back to complete)
build_dataset_graph_ <- function(N, graph = c("complete", "star", "mst", "knn"), domains = NULL, algo = NULL, graph_opts = NULL, ...) {
  graph <- match.arg(graph)
  if (graph == "complete") {
    ij <- which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
    return(data.frame(i = ij[, 1], j = ij[, 2]))
  }
  if (graph == "star") {
    return(data.frame(i = rep(1L, N - 1L), j = 2:N))
  }
  # For mst/knn we need pairwise scores; smaller is better
  scores <- compute_dataset_pair_scores_(domains, algo, ...)
  if (graph == "mst") {
    edges <- mst_from_scores_(scores)
    return(edges)
  }
  # knn
  k <- 2L
  if (is.list(graph_opts) && !is.null(graph_opts$k)) k <- as.integer(graph_opts$k)
  edges <- knn_edges_from_scores_(scores, k)
  edges
}

# Compute a symmetric matrix of dataset-pair scores using light landmarks
compute_dataset_pair_scores_ <- function(domains, algo, n_landmarks = 200, q = 8, seed = 1L, ...) {
  set.seed(seed)
  N <- length(domains)
  S <- matrix(0, N, N)
  # Precompute PCA projections on small random subsets for speed
  pca_list <- vector("list", N)
  for (i in seq_len(N)) {
    Xi <- as.matrix(domains[[i]])
    if (!is.matrix(Xi)) Xi <- as.matrix(domains[[i]]$x)
    ni <- nrow(Xi)
    li <- min(ni, n_landmarks)
    idx <- sample.int(ni, li)
    Xs <- scale(Xi[idx, , drop = FALSE], center = TRUE, scale = FALSE)
    # Low-rank PCA
    qq <- min(q, ncol(Xs), nrow(Xs))
    pc <- tryCatch({
      if (requireNamespace("irlba", quietly = TRUE) && qq < min(nrow(Xs), ncol(Xs))) {
        irlba::irlba(Xs, nv = qq, nu = qq)
      } else {
        svd(Xs, nu = qq, nv = qq)
      }
    }, error = function(e) svd(Xs, nu = qq, nv = qq))
    Zi <- Xs %*% pc$v
    pca_list[[i]] <- Zi
  }
  for (i in seq_len(N - 1L)) {
    for (j in (i + 1L):N) {
      Zi <- pca_list[[i]]; Zj <- pca_list[[j]]
      # Score: mean nearest neighbor distance between Zi and Zj (symmetric)
      d_ij <- tryCatch({
        if (requireNamespace("RANN", quietly = TRUE)) {
          nn <- RANN::nn2(Zj, Zi, k = 1)
          mean(nn$nn.dists)
        } else {
          D <- as.matrix(dist(rbind(Zi, Zj)))
          n_i <- nrow(Zi); n_j <- nrow(Zj)
          minD <- apply(D[seq_len(n_i), (n_i + 1):(n_i + n_j), drop = FALSE], 1, min)
          mean(minD)
        }
      }, error = function(e) {
        # Fallback: norm difference of covariance matrices
        Ci <- cov(Zi); Cj <- cov(Zj)
        sqrt(mean((Ci - Cj)^2))
      })
      S[i, j] <- S[j, i] <- d_ij
    }
  }
  S
}

# Build MST edges from a symmetric score matrix
mst_from_scores_ <- function(S) {
  N <- nrow(S)
  # Prim's algorithm over dense matrix
  inT <- rep(FALSE, N); inT[1] <- TRUE
  edges <- data.frame(i = integer(0), j = integer(0))
  while (sum(inT) < N) {
    best <- Inf; bi <- NA_integer_; bj <- NA_integer_
    for (i in which(inT)) {
      j_candidates <- which(!inT)
      if (!length(j_candidates)) next
      j <- j_candidates[which.min(S[i, j_candidates])]
      if (S[i, j] < best) { best <- S[i, j]; bi <- i; bj <- j }
    }
    edges <- rbind(edges, data.frame(i = bi, j = bj))
    inT[bj] <- TRUE
  }
  edges
}

# Build undirected k-NN edges from a symmetric score matrix
knn_edges_from_scores_ <- function(S, k = 2L) {
  N <- nrow(S)
  adj <- matrix(FALSE, N, N)
  for (i in seq_len(N)) {
    ord <- order(S[i, ], decreasing = FALSE)
    nbrs <- ord[ord != i][seq_len(min(k, N - 1L))]
    adj[i, nbrs] <- TRUE
  }
  # symmetrize
  adj <- adj | t(adj)
  ij <- which(upper.tri(adj) & adj, arr.ind = TRUE)
  data.frame(i = ij[, 1], j = ij[, 2])
}

fit_pairs_on_edges_ <- function(domains, E, algo, parallel = TRUE, ...) {
  # For now, simple sequential fits; adapters decide how to use domains[[i/j]]
  lapply(seq_len(nrow(E)), function(idx) {
    ii <- E$i[idx]; jj <- E$j[idx]
    fit_pair(algo, domains[[ii]], domains[[jj]], links = NULL, i = ii, j = jj, ...)
  })
}

# -- Synchronization backends ----------------------------------------------

#' Rotation synchronization for orthogonal group O(k)
#'
#' @param G list of align_transform objects (group="O"), each a relative g_ij
#' @param W numeric weights per edge in E
#' @param E data.frame with columns i, j (1-based dataset indices)
#' @param N number of datasets
#' @param k latent dimension
#' @return list of per-domain align_transform objects mapping to consensus
#' @examples
#' \donttest{
#' # Rotation synchronization for orthogonal transformations
#' k <- 3
#' N <- 3
#' G <- list(matrix(rnorm(9), k, k), matrix(rnorm(9), k, k))
#' W <- c(1, 1)
#' E <- data.frame(i = c(1, 2), j = c(2, 3))
#' # result <- rotation_sync(G, W, E, N = N, k = k)
#' }
#' @export
rotation_sync <- function(G, W, E, N, k) {
  if (missing(N)) N <- max(c(E$i, E$j))
  if (is.null(k)) {
    k <- tryCatch({
      g1 <- G[[1]]; if (inherits(g1, "align_transform")) ncol(g1$op) else ncol(as.matrix(g1))
    }, error = function(e) NULL)
  }
  if (is.null(k)) stop("rotation_sync: could not infer latent dimension k")

  # Build symmetric block matrix S of size (N*k) x (N*k) with off-diagonal blocks W_ij * g_ij
  nblk <- N * k
  Ii <- integer(0); Jj <- integer(0); Xx <- numeric(0)
  deg <- numeric(N)

  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]
    if (!inherits(gij, "align_transform")) {
      gij <- new_align_transform("O", as.matrix(gij), from = i, to = j, k = k)
    }
    R <- as.matrix(gij$op)
    w <- as.numeric(W[idx] %||% 1)
    deg[i] <- deg[i] + w
    deg[j] <- deg[j] + w

    # block indices
    ri <- (i - 1L) * k + 1L; rj <- (j - 1L) * k + 1L
    # Indices must match R's column-major vectorization (as.vector()).
    ii <- rep(ri:(ri + k - 1L), times = k)
    jj <- rep(rj:(rj + k - 1L), each = k)
    # S_ij += w * R
    Ii <- c(Ii, ii)
    Jj <- c(Jj, jj)
    Xx <- c(Xx, as.vector(w * R))
    # S_ji += w * R^T
    ii2 <- rep(rj:(rj + k - 1L), times = k)
    jj2 <- rep(ri:(ri + k - 1L), each = k)
    Xx2 <- as.vector(w * t(R))
    Ii <- c(Ii, ii2)
    Jj <- c(Jj, jj2)
    Xx <- c(Xx, Xx2)
  }

  S <- Matrix::sparseMatrix(i = Ii, j = Jj, x = Xx, dims = c(nblk, nblk))
  # Normalize by degrees: D^{-1/2} S D^{-1/2}
  blocks <- vector("list", N)
  for (i in seq_len(N)) {
    di <- deg[i]
    if (!is.finite(di) || di <= 0) {
      blocks[[i]] <- Matrix::Diagonal(k, x = 0)
    } else {
      blocks[[i]] <- Matrix::Diagonal(k, x = rep(1 / sqrt(di), k))
    }
  }
  Dm <- Matrix::bdiag(blocks)
  M <- Matrix::forceSymmetric(Dm %*% Matrix::forceSymmetric(S) %*% Dm)

  # Top-k eigenvectors
  V <- tryCatch({
    if (requireNamespace("RSpectra", quietly = TRUE) && k < nrow(M)) {
      RSpectra::eigs_sym(M, k = k, which = "LA")$vectors
    } else if (requireNamespace("PRIMME", quietly = TRUE) && k < nrow(M)) {
      PRIMME::eigs_sym(M, NEig = k, which = "LA")$vectors
    } else {
      eigen(as.matrix(M), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
    }
  }, error = function(e) {
    eigen(as.matrix(M), symmetric = TRUE)$vectors[, seq_len(k), drop = FALSE]
  })

  # Slice blocks and project to O(k)
  T_list <- vector("list", N)
  for (i in seq_len(N)) {
    rows <- ((i - 1L) * k + 1L):(i * k)
    Y <- as.matrix(V[rows, , drop = FALSE])
    sv <- tryCatch(svd(Y), error = function(e) NULL)
    if (is.null(sv)) {
      R <- diag(k)
    } else {
      R <- sv$u %*% t(sv$v)
    }
    T_list[[i]] <- new_align_transform("O", R, from = i, to = "consensus", k = k)
  }

  T_list
}

# Placeholders for other groups ------------------------------------------------

gl_sync_least_squares <- function(G, W, E, N, k, lambda = 1e-6, iters = 20, gauge_idx = 1L, project = c("none","orthogonal")) {
  project <- match.arg(project)
  if (is.null(k)) {
    for (g in G) { op <- if (inherits(g, "align_transform")) g$op else g; k <- nrow(as.matrix(op)); break }
  }
  Tlist <- replicate(N, diag(k), simplify = FALSE)
  Tlist[[gauge_idx]] <- diag(k)

  # adjacency list with oriented ops g_ij mapping i->j
  neigh <- vector("list", N)
  for (i in seq_len(N)) neigh[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]
    if (!inherits(gij, "align_transform")) gij <- new_align_transform("GL", as.matrix(gij), from = i, to = j)
    neigh[[i]][[length(neigh[[i]]) + 1L]] <- list(j = j, G = as.matrix(gij$op), w = as.numeric(W[idx] %||% 1))
    # reverse direction
    Gji <- tryCatch(solve(as.matrix(gij$op)), error = function(e) NULL)
    if (!is.null(Gji)) {
      neigh[[j]][[length(neigh[[j]]) + 1L]] <- list(j = i, G = Gji, w = as.numeric(W[idx] %||% 1))
    }
  }

  for (it in seq_len(iters)) {
    for (i in seq_len(N)) {
      if (i == gauge_idx) next
      Ai <- matrix(0, k, k); Bi <- matrix(0, k, k)
      for (nb in neigh[[i]]) {
        Gij <- nb$G; w <- nb$w; j <- nb$j
        Ai <- Ai + w * crossprod(Gij, Gij)
        Bi <- Bi + w * crossprod(Gij, Tlist[[j]])
      }
      Ai <- Ai + lambda * diag(k)
      Ti <- tryCatch(solve(Ai, Bi), error = function(e) diag(k))
      if (project == "orthogonal") { sv <- svd(Ti); Ti <- sv$u %*% t(sv$v) }
      Tlist[[i]] <- Ti
    }
  }
  lapply(seq_len(N), function(i) new_align_transform("GL", Tlist[[i]], from = i, to = "consensus", k = k))
}

# Star/hub consensus for permutations (privileged domain 1)
permutation_sync_star <- function(G, W, E, N) {
  # Map each domain i to hub 1 using available edges; prefer direct (i,1), else invert (1,i)
  # Build lookup for edges
  key <- paste(E$i, E$j, sep = ":")
  idx_by_pair <- setNames(seq_along(key), key)
  get_map <- function(a, b) {
    if (a == b) return(new_align_transform("perm", diag(1), from = a, to = b))
    k1 <- paste(a, b, sep = ":"); k2 <- paste(b, a, sep = ":")
    if (!is.null(idx_by_pair[[k1]])) {
      gij <- G[[idx_by_pair[[k1]]]]
      return(if (inherits(gij, "align_transform")) gij else new_align_transform("perm", as.matrix(gij), from = a, to = b))
    }
    if (!is.null(idx_by_pair[[k2]])) {
      gji <- G[[idx_by_pair[[k2]]]]
      if (inherits(gji, "align_transform")) return(invert_transform(gji))
      op <- as.matrix(gji)
      rs <- rowSums(t(op)); rs[!is.finite(rs) | rs == 0] <- 1
      inv <- t(op) / rs
      return(new_align_transform("perm", inv, from = a, to = b))
    }
    stop("No path between ", a, " and ", b, " for star consensus")
  }
  # Infer hub size
  n1 <- NULL
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    if (i == 1 || j == 1) {
      gij <- G[[idx]]; op <- if (inherits(gij, "align_transform")) as.matrix(gij$op) else as.matrix(gij)
      if (i == 1) { n1 <- ncol(op); break } else { n1 <- nrow(op); break }
    }
  }
  if (is.null(n1)) n1 <- 1L
  T_list <- vector("list", N)
  T_list[[1]] <- new_align_transform("perm", diag(n1), from = 1, to = "consensus")
  for (i in 2:N) {
    g <- get_map(i, 1)
    T_list[[i]] <- new_align_transform("perm", g$op, from = i, to = "consensus")
  }
  T_list
}

# Spectral consensus for permutations (cycle-aware)
permutation_sync_spectral <- function(G, W, E, N, L = NULL,
                                      rounding = c("sinkhorn", "hungarian"),
                                      max_sinkhorn = 200, tol = 1e-6) {
  rounding <- match.arg(rounding)
  # Infer sizes n_i
  n <- rep(NA_integer_, N)
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]
    op <- if (inherits(gij, "align_transform")) as.matrix(gij$op) else as.matrix(gij)
    n[i] <- n[i] %||% ncol(op)
    n[j] <- n[j] %||% nrow(op)
  }
  if (any(is.na(n))) n[is.na(n)] <- min(n, na.rm = TRUE)
  if (is.null(L)) L <- min(n)

  offs <- c(0L, cumsum(n)); dim_big <- sum(n)
  Ii <- integer(0); Jj <- integer(0); Xx <- numeric(0)
  add_block <- function(ri, rj, B) {
    Ii <<- c(Ii, rep(ri, each = length(rj)))
    Jj <<- c(Jj, rep(rj, times = length(ri)))
    Xx <<- c(Xx, as.vector(B))
  }
  # Assemble symmetric block matrix with P_ji and P_ij
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]
    if (!inherits(gij, "align_transform")) gij <- new_align_transform("perm", as.matrix(gij), from = i, to = j)
    Pij <- as.matrix(gij$op)  # n_j x n_i
    Pji <- t(Pij)             # n_i x n_j
    w <- as.numeric(W[idx] %||% 1)
    ri <- (offs[i] + 1):(offs[i + 1]); rj <- (offs[j] + 1):(offs[j + 1])
    add_block(ri, rj, w * Pji)
    add_block(rj, ri, w * Pij)
  }
  M <- Matrix::sparseMatrix(i = Ii, j = Jj, x = Xx, dims = c(dim_big, dim_big))
  M <- Matrix::forceSymmetric(M)

  # Top-L eigenvectors
  Y <- tryCatch({
    if (requireNamespace("RSpectra", quietly = TRUE) && L < nrow(M)) {
      RSpectra::eigs_sym(M, k = L, which = "LM")$vectors
    } else if (requireNamespace("PRIMME", quietly = TRUE) && L < nrow(M)) {
      PRIMME::eigs_sym(M, NEig = L, which = "LM")$vectors
    } else {
      eigen(as.matrix(M), symmetric = TRUE)$vectors[, seq_len(L), drop = FALSE]
    }
  }, error = function(e) {
    eigen(as.matrix(M), symmetric = TRUE)$vectors[, seq_len(L), drop = FALSE]
  })

  # Rectangular Sinkhorn to match row/col marginals
  sinkhorn_rect <- function(A, r, c, iters = max_sinkhorn, tol = tol) {
    A[A < 0] <- 0
    u <- rep(1, length(r)); v <- rep(1, length(c))
    for (t in seq_len(iters)) {
      Av <- A %*% v; u <- r / pmax(Av, .Machine$double.eps)
      ATu <- t(A) %*% u; v <- c / pmax(ATu, .Machine$double.eps)
      if (t %% 10 == 0) {
        X <- (u * A) %*% diag(v)
        if (max(abs(rowSums(X) - r)) < tol && max(abs(colSums(X) - c)) < tol) break
      }
    }
    (u * A) %*% diag(v)
  }

  T_list <- vector("list", N)
  for (i in seq_len(N)) {
    rows <- (offs[i] + 1):(offs[i + 1])
    A <- abs(Y[rows, , drop = FALSE])
    r <- rep(1 / n[i], n[i]); c <- rep(1 / L, L)
    Q <- if (rounding == "hungarian" && n[i] == L) {
      if (!requireNamespace("clue", quietly = TRUE)) stop("clue package required for Hungarian rounding")
      cost <- 1 - A / (max(A) + 1e-12)
      m <- clue::solve_LSAP(as.matrix(cost))
      Mperm <- matrix(0, nrow = n[i], ncol = L); Mperm[cbind(seq_len(n[i]), m)] <- 1
      Mperm
    } else {
      sinkhorn_rect(A, r, c)
    }
    T_list[[i]] <- new_align_transform("perm", t(Q), from = i, to = "consensus")
  }
  T_list
}

permutation_sync <- function(G, W, E, N, max_iter = 25, tol = 1e-6) {
  # Synchronize soft maps without hard hub dependence using an averaging
  # iteration: T_i ≈ average_j T_j %*% g_{ji}. Fix gauge with T_1 = I.

  # Row-normalization helper
  normalize_rows <- function(M) {
    rs <- rowSums(M)
    rs[!is.finite(rs) | rs == 0] <- 1
    M / rs
  }

  # Infer per-domain sizes n[i] from edge operators (g_ij: n_j x n_i)
  n <- rep(NA_integer_, N)
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]
    op <- if (inherits(gij, "align_transform")) as.matrix(gij$op) else as.matrix(gij)
    n[j] <- n[j] %||% nrow(op)
    n[i] <- n[i] %||% ncol(op)
  }
  # Fallback for any unresolved sizes
  n[is.na(n)] <- n[which(!is.na(n))[1]] %||% 1L

  # Consensus rows set to domain 1 size (gauge)
  m <- n[1]
  if (is.na(m) || m <= 0) m <- 1L

  # Build neighbor list with oriented maps g_{ji}
  neighbors <- vector("list", N)
  for (i in seq_len(N)) neighbors[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    a <- E$i[idx]; b <- E$j[idx]
    gij <- G[[idx]]
    if (!inherits(gij, "align_transform")) {
      gij <- new_align_transform("perm", as.matrix(gij), from = a, to = b)
    }
    # Store g_ab and its inverse g_ba
    neighbors[[a]][[length(neighbors[[a]]) + 1L]] <- list(j = b, op = invert_transform(gij)$op, w = W[idx])
    neighbors[[b]][[length(neighbors[[b]]) + 1L]] <- list(j = a, op = gij$op,                    w = W[idx])
  }

  # Initialize T_i: identity for i=1, uniform for others
  T <- vector("list", N)
  T[[1]] <- diag(m, n[1])
  for (i in 2:N) {
    T[[i]] <- matrix(1 / n[i], nrow = m, ncol = n[i])
  }

  # Iterative averaging
  for (iter in seq_len(max_iter)) {
    T_new <- T
    max_diff <- 0
    for (i in seq_len(N)) {
      if (i == 1) next  # gauge
      acc <- matrix(0, nrow = m, ncol = n[i])
      wsum <- 0
      for (nb in neighbors[[i]]) {
        j <- nb$j
        op_ji <- as.matrix(nb$op) # n_i x n_j
        w <- as.numeric(nb$w %||% 1)
        acc <- acc + w * (T[[j]] %*% op_ji)
        wsum <- wsum + w
      }
      if (wsum > 0) acc <- acc / wsum
      T_new[[i]] <- normalize_rows(acc)
      d <- max(abs(T_new[[i]] - T[[i]]))
      if (is.finite(d) && d > max_diff) max_diff <- d
    }
    T <- T_new
    if (max_diff < tol) break
  }

  # Package as transforms to a shared consensus
  lapply(seq_len(N), function(i) new_align_transform("perm", T[[i]], from = i, to = "consensus"))
}

# -- Diagnostics ------------------------------------------------------------

compute_diagnostics_ <- function(group, T_list, G, E, W) {
  # Edge residuals: measure ||constraint|| under learned transforms
  edge_res <- rep(NA_real_, nrow(E))
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    edge_res[idx] <- tryCatch({
      if (group %in% c("O", "GL")) {
        A <- as.matrix(T_list[[j]]$op)
        B <- if (inherits(G[[idx]], "align_transform")) as.matrix(G[[idx]]$op) else as.matrix(G[[idx]])
        C <- as.matrix(T_list[[i]]$op)
        D <- A - B %*% C
        sqrt(mean(D^2))
      } else if (group == "perm") {
        A <- as.matrix(T_list[[i]]$op) # L x n_i
        B <- as.matrix(T_list[[j]]$op) # L x n_j
        Gop <- if (inherits(G[[idx]], "align_transform")) as.matrix(G[[idx]]$op) else as.matrix(G[[idx]]) # n_j x n_i
        D <- A - B %*% Gop
        sqrt(mean(D^2))
      } else NA_real_
    }, error = function(e) NA_real_)
  }
  edge_df <- data.frame(i = E$i, j = E$j, residual = edge_res)
  edge_wt <- data.frame(i = E$i, j = E$j, weight = as.numeric(W))

  # Cycle residuals on triangles (if present)
  cyc_df <- NULL
  if (nrow(E) >= 3) {
    key <- paste(E$i, E$j, sep = ":")
    idx_map <- setNames(seq_len(nrow(E)), key)
    tri <- list(i = integer(0), j = integer(0), k = integer(0), residual = numeric(0))
    Nmax <- max(c(E$i, E$j))
    for (i in seq_len(Nmax - 2)) for (j in (i + 1):(Nmax - 1)) for (k in (j + 1):Nmax) {
      kij <- idx_map[[paste(i, j, sep = ":")]]
      kjk <- idx_map[[paste(j, k, sep = ":")]]
      kik <- idx_map[[paste(i, k, sep = ":")]]
      if (is.null(kij) || is.null(kjk) || is.null(kik)) next
      Gij <- if (inherits(G[[kij]], "align_transform")) as.matrix(G[[kij]]$op) else as.matrix(G[[kij]])
      Gjk <- if (inherits(G[[kjk]], "align_transform")) as.matrix(G[[kjk]]$op) else as.matrix(G[[kjk]])
      Gik <- if (inherits(G[[kik]], "align_transform")) as.matrix(G[[kik]]$op) else as.matrix(G[[kik]])
      D <- tryCatch(Gik - Gjk %*% Gij, error = function(e) NULL)
      if (!is.null(D)) {
        tri$i <- c(tri$i, i); tri$j <- c(tri$j, j); tri$k <- c(tri$k, k); tri$residual <- c(tri$residual, sqrt(mean(D^2)))
      }
    }
    if (length(tri$residual)) cyc_df <- data.frame(i = tri$i, j = tri$j, k = tri$k, residual = tri$residual)
  }

  # Entropy metrics for permutation group
  entropy_df <- NULL
  if (identical(group, "perm")) {
    ent <- lapply(seq_along(T_list), function(i) {
      Ti <- as.matrix(T_list[[i]]$op) # L x n_i
      L <- nrow(Ti); n <- ncol(Ti); eps <- 1e-12
      cs <- colSums(Ti); Cn <- sweep(Ti, 2, pmax(cs, eps), "/")
      rs <- rowSums(Ti); Rn <- sweep(Ti, 1, pmax(rs, eps), "/")
      Hc <- -colSums(Cn * log(pmax(Cn, eps))) / log(max(L, 2))
      Hr <- -rowSums(Rn * log(pmax(Rn, eps))) / log(max(n, 2))
      data.frame(domain = i, L = L, n = n,
                 mean_col_entropy = mean(Hc),
                 mean_row_entropy = mean(Hr))
    })
    entropy_df <- do.call(rbind, ent)
  }

  list(
    weights = W,
    edge_residuals = edge_df,
    edge_weights = edge_wt,
    mean_edge_residual = if (all(is.na(edge_res))) NA_real_ else mean(edge_res, na.rm = TRUE),
    cycle_residuals = cyc_df,
    mean_cycle_residual = if (is.null(cyc_df)) NA_real_ else mean(cyc_df$residual, na.rm = TRUE),
    entropy = entropy_df
  )
}

# One-pass local refinement for O(k) synchronization
refine_orthogonal_local_ <- function(T_list, G, W, E, N) {
  k <- ncol(as.matrix(T_list[[1]]$op))
  out <- T_list
  # Build neighbor list with oriented g_ij (i->j)
  neigh <- vector("list", N)
  for (i in seq_len(N)) neigh[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    gij <- G[[idx]]; Gij <- if (inherits(gij, "align_transform")) as.matrix(gij$op) else as.matrix(gij)
    w <- as.numeric(W[idx] %||% 1)
    neigh[[i]][[length(neigh[[i]]) + 1L]] <- list(j = j, G = Gij, w = w)
    neigh[[j]][[length(neigh[[j]]) + 1L]] <- list(j = i, G = t(Gij), w = w) # reverse for orthogonal
  }
  for (i in seq_len(N)) {
    Mi <- matrix(0, k, k)
    for (nb in neigh[[i]]) {
      j <- nb$j; Gij <- nb$G; w <- nb$w
      Mi <- Mi + w * t(Gij) %*% as.matrix(out[[j]]$op)
    }
    sv <- tryCatch(svd(Mi), error = function(e) NULL)
    if (!is.null(sv)) {
      Ri <- sv$u %*% t(sv$v)
      out[[i]] <- new_align_transform("O", Ri, from = i, to = "consensus", k = k)
    }
  }
  out
}
