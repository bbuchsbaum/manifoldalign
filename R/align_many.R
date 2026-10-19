#' Align multiple datasets by synchronizing pairwise transforms
#'
#' This meta-aligner lifts any pairwise method to N domains by:
#' 1) choosing a sparse dataset graph, 2) fitting pairwise models on edges,
#' 3) extracting relative transforms, and 4) running group-appropriate
#' synchronization to get per-domain transforms into a shared space.
#'
#' If the adapter reports native multi-view support, dispatches to fit_many().
#'
#' @section Orientation:
#' All transforms follow the contract documented in [new_align_transform()]:
#' the edge transform `g_ij = relative_transform(fit_ij, "i", "j")` maps
#' domain `i` onto domain `j`, and the returned `transforms[[i]]` maps domain
#' `i` onto the consensus space, so that `embeddings[[i]] =
#' apply_transform(transforms[[i]], domains[[i]])` agree across domains.
#' Synchronisation, the optional `sync$refine$gpa` refinement and the
#' diagnostics all solve or measure the same relation
#' `compose_transform(g_ij, T_j) = T_i`.
#'
#' @param domains list of domain objects (often matrices), one per dataset
#' @param algo an aligner object created by new_aligner() with adapter methods
#' @param graph "auto" | "mst" | "complete" | "star" | "knn". `"auto"` picks
#'   `"complete"` for up to 5 domains, `"knn"` for up to 8 and `"mst"` beyond.
#' @param consensus "auto" or "sync" (pairwise synchronisation). `"star"` and
#'   `"tree"` are not implemented and are rejected; for permutation aligners a
#'   hub consensus is selected with `sync = list(perm = list(mode = "star"))`.
#' @param k target latent dimension (optional). Native multi-domain aligners
#'   receive it as their number of components (`ncomp`).
#' @param parallel logical; reserved for future parallel pair fits
#' @param sync list of synchronization options (graph, refine, gl, perm settings).
#'   `sync$gl` accepts `ridge`, `gauge`, `project` (`"none"` solves the
#'   least-squares synchronisation exactly; `"orthogonal"` adds `iters`
#'   projected sweeps). `sync$refine = list(gpa = TRUE)` adds one local O(k)
#'   refinement sweep.
#'   `sync$perm` accepts `mode` (`"star"`, `"global"`, `"spectral"`; default
#'   `"spectral"` for more than four domains, else `"star"`), `hub` (the domain
#'   whose rows form the consensus; default the largest domain, the first on
#'   ties), `k` (number of spectral eigenvectors), `rounding` (`"sinkhorn"` or
#'   `"hungarian"`), `maxit_sinkhorn`, `tol` and `max_iter`. Permutation
#'   consensus transforms are `n_hub x n_i`; for hard (partial) permutations
#'   every row of domain `i` maps to one consensus row and consensus rows
#'   absent from domain `i` stay zero.
#' @param links Optional correspondences for the pairwise fits. Either a list
#'   named by edge keys `"i:j"` (each element is passed as `links` to
#'   [fit_pair()] for the edge from domain `i` to domain `j`), or, for exactly
#'   two domains, a single links object. Keys that do not name an edge of the
#'   dataset graph are rejected. For native multi-domain aligners `links` is
#'   forwarded to [fit_many()].
#' @param ... forwarded to adapter methods
#' @return list with transforms, embeddings (`NULL` for a domain whose
#'   transform cannot be applied to its data, with a warning naming the cause),
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
                       graph = c("auto", "mst", "complete", "star", "knn"),
                       consensus = c("auto", "sync", "star", "tree"),
                       k = NULL, parallel = TRUE,
                       sync = NULL, links = NULL, ...) {

  graph <- match.arg(graph)
  consensus <- match.arg(consensus)
  if (consensus %in% c("star", "tree")) {
    stop("align_many: consensus = \"", consensus, "\" is not implemented; ",
         "use consensus = \"sync\" (for permutation aligners a hub consensus is ",
         "available via sync = list(perm = list(mode = \"star\"))).", call. = FALSE)
  }

  caps <- aligner_capabilities(algo)
  if (isTRUE(caps$supports_multi)) {
    if (is.null(links)) {
      return(fit_many(algo, domains, k = k, ...))
    }
    return(fit_many(algo, domains, k = k, links = links, ...))
  }

  grp <- caps$group
  if (!is.character(grp) || length(grp) != 1L || !grp %in% c("O", "GL", "perm")) {
    stop("align_many: aligners of group '", paste(grp, collapse = ", "),
         "' have no relative transform to synchronise pairwise; only groups ",
         "'O', 'GL' and 'perm' are supported (native multi-domain aligners ",
         "must declare supports_multi = TRUE).", call. = FALSE)
  }

  N <- length(domains)
  if (N < 2) stop("align_many: need at least two domains")

  # Choose dataset graph
  if (graph == "auto") {
    graph <- if (N <= 5) "complete" else if (N <= 8) "knn" else "mst"
  }
  graph_opts <- tryCatch(sync$graph, error = function(e) NULL)
  E <- build_dataset_graph_(N, graph = graph, domains = domains, algo = algo, graph_opts = graph_opts, ...)

  edge_links <- .align_many_edge_links(links, E, N)

  # Pairwise fits
  fits <- fit_pairs_on_edges_(domains, E, algo, parallel = parallel,
                              edge_links = edge_links, ...)

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
  if (is.null(k) && grp %in% c("O", "GL")) {
    k <- tryCatch({
      ks <- vapply(G, function(g) {
        if (inherits(g, "align_transform") && !is.null(g$k)) return(as.numeric(g$k))
        op <- if (inherits(g, "align_transform")) g$op else as.matrix(g)
        as.numeric(ncol(op))
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
      if (is.list(gl_opts) && !is.null(gl_opts$iters) && identical(project, "none")) {
        stop("align_many: sync$gl$iters only applies with sync$gl$project = ",
             "\"orthogonal\"; the unprojected GL synchronisation is solved exactly.",
             call. = FALSE)
      }
      gl_sync_least_squares(G, W, E, N = N, k = k, lambda = lambda, iters = iters, gauge_idx = gauge, project = project)
    },
    "perm" = {
      perm_opts <- tryCatch(sync$perm, error = function(e) NULL)
      if (!is.null(perm_opts) && !is.list(perm_opts)) {
        stop("align_many: sync$perm must be a list.", call. = FALSE)
      }
      mode <- perm_opts$mode %||% (if (N > 4) "spectral" else "star")
      mode <- match.arg(mode, c("star", "global", "spectral"))
      if (mode == "star") {
        permutation_sync_star(G, W, E, N = N, hub = perm_opts$hub)
      } else if (mode == "global") {
        permutation_sync(G, W, E, N = N,
                         max_iter = perm_opts$max_iter %||% 100L,
                         tol = perm_opts$tol %||% 1e-8,
                         gauge = perm_opts$hub)
      } else { # spectral
        permutation_sync_spectral(G, W, E, N = N,
                                  L = perm_opts$k,
                                  rounding = perm_opts$rounding %||% "sinkhorn",
                                  max_sinkhorn = perm_opts$maxit_sinkhorn %||% 200,
                                  tol = perm_opts$tol %||% 1e-6,
                                  ref = perm_opts$hub)
      }
    }
  )

  # Embed each domain with its consensus transform. A transform that does not
  # act on the raw domain data (for example a map between spectral bases)
  # yields NULL with a warning naming the cause.
  Z <- lapply(seq_len(N), function(i) {
    tryCatch(apply_transform(T_list[[i]], domains[[i]]), error = function(e) {
      warning("align_many: could not embed domain ", i,
              " with its consensus transform: ", conditionMessage(e),
              call. = FALSE)
      NULL
    })
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
# star: hub = 1; complete: all pairs i < j; mst: Prim's tree on landmark
# scores; knn: symmetrised k-nearest-neighbour graph on the same scores.
build_dataset_graph_ <- function(N, graph = c("complete", "star", "mst", "knn"), domains = NULL, algo = NULL, graph_opts = NULL, ...) {
  graph <- match.arg(graph)
  if (graph == "complete") {
    ij <- which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
    ij <- ij[order(ij[, 1], ij[, 2]), , drop = FALSE]
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
  .with_local_seed(seed, {
  N <- length(domains)
  S <- matrix(0, N, N)
  # Precompute PCA projections on small random subsets for speed
  pca_list <- vector("list", N)
  for (i in seq_len(N)) {
    Xi <- domains[[i]]
    if (is.list(Xi) && !is.null(Xi$x)) Xi <- Xi$x
    ni <- nrow(Xi)
    li <- min(ni, n_landmarks)
    idx <- sample.int(ni, li)
    # Subset before densifying so sparse domains are never copied whole
    Xs <- scale(as.matrix(Xi[idx, , drop = FALSE]), center = TRUE, scale = FALSE)
    # Low-rank PCA
    qq <- min(q, ncol(Xs), nrow(Xs))
    pc <- tryCatch({
      if (requireNamespace("irlba", quietly = TRUE) && qq < min(nrow(Xs), ncol(Xs))) {
        .irlba_det(Xs, nv = qq, nu = qq)
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
  })
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

# Resolve user-supplied links into one entry per edge of E.
.align_many_edge_links <- function(links, E, N) {
  out <- vector("list", nrow(E))
  if (is.null(links)) return(out)
  edge_keys <- paste(E$i, E$j, sep = ":")
  nms <- names(links)
  keyed <- is.list(links) && length(links) > 0L && !is.null(nms) &&
    all(grepl("^[0-9]+:[0-9]+$", nms))
  if (!keyed) {
    if (N == 2L && nrow(E) == 1L) {
      out[[1L]] <- links
      return(out)
    }
    stop("align_many: `links` must be a list named by edge keys \"i:j\" ",
         "(one entry per dataset-graph edge); a single links object is only ",
         "accepted for two domains.", call. = FALSE)
  }
  unknown <- setdiff(nms, edge_keys)
  if (length(unknown)) {
    reversed <- unknown[vapply(strsplit(unknown, ":", fixed = TRUE), function(p) {
      paste(p[2], p[1], sep = ":") %in% edge_keys
    }, logical(1))]
    msg <- paste0("align_many: `links` names pair(s) that are not edges of the ",
                  "dataset graph: ", paste(unknown, collapse = ", "), ". Edges are: ",
                  paste(edge_keys, collapse = ", "), ".")
    if (length(reversed)) {
      msg <- paste0(msg, " Links for ", paste(reversed, collapse = ", "),
                    " are given in the opposite orientation to the edge.")
    }
    stop(msg, call. = FALSE)
  }
  for (nm in nms) out[[match(nm, edge_keys)]] <- links[[nm]]
  out
}

fit_pairs_on_edges_ <- function(domains, E, algo, parallel = TRUE,
                                edge_links = NULL, ...) {
  if (is.null(edge_links)) edge_links <- vector("list", nrow(E))
  # For now, simple sequential fits; adapters decide how to use domains[[i/j]]
  lapply(seq_len(nrow(E)), function(idx) {
    ii <- E$i[idx]; jj <- E$j[idx]
    fit_pair(algo, domains[[ii]], domains[[jj]], links = edge_links[[idx]],
             i = ii, j = jj, ...)
  })
}

# -- Shared helpers ----------------------------------------------------------

# Operator of a transform given either as an align_transform or a raw matrix.
.transform_op <- function(g) {
  if (inherits(g, "align_transform")) g$op else g
}

# Coerce an edge transform to an align_transform oriented from -> to.
.as_edge_transform <- function(g, group, from, to, k = NULL) {
  if (inherits(g, "align_transform")) return(g)
  new_align_transform(group, as.matrix(g), from = from, to = to, k = k)
}

# Top-k eigenvectors (largest algebraic eigenvalues) of a symmetric matrix.
# Synchronisation matrices have highly degenerate leading eigenvalues, which
# Krylov solvers can resolve incompletely, so moderate problems use a dense
# decomposition. Larger problems use RSpectra, then PRIMME when installed.
.sync_top_eigenvectors <- function(M, k, dense_max = 1000L) {
  n <- nrow(M)
  dense <- function() {
    e <- eigen(as.matrix(M), symmetric = TRUE)
    e$vectors[, seq_len(k), drop = FALSE]
  }
  if (k >= n || n <= dense_max) return(dense())
  tryCatch({
    if (requireNamespace("RSpectra", quietly = TRUE)) {
      RSpectra::eigs_sym(M, k = k, which = "LA")$vectors
    } else if (.primme_available()) {
      PRIMME::eigs_sym(M, NEig = k, which = "LA")$vectors
    } else {
      dense()
    }
  }, error = function(e) dense())
}

# -- Synchronization backends ----------------------------------------------

#' Rotation synchronization for orthogonal group O(k)
#'
#' Solves `T_i ~ g_ij %*% T_j` over all edges (the relation implied by the
#' transform contract in [new_align_transform()]) by the spectral relaxation:
#' the top-`k` eigenvectors of the degree-normalised block matrix with blocks
#' `W_ij g_ij` (and `W_ij t(g_ij)` below the diagonal), projected blockwise to
#' O(k).
#'
#' @param G list of align_transform objects (group="O"), each a relative g_ij
#'   from domain `E$i` to domain `E$j`
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
  if (missing(k) || is.null(k)) {
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
    gij <- .as_edge_transform(G[[idx]], "O", i, j, k = k)
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
  V <- .sync_top_eigenvectors(M, k)

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

# Least-squares synchronisation on GL(k):
#   min sum_{(i,j) in E} w_ij || T_i - g_ij T_j ||_F^2 + lambda sum_i ||T_i||_F^2
# with T_gauge = I fixed. The objective is quadratic, so with
# project = "none" the joint normal equations are solved exactly. For domain i
# the block row reads
#   (sum_out w I + sum_in w g_ai' g_ai + lambda I) T_i
#     - sum_out w g_ij T_j - sum_in w g_ai' T_a = 0,
# where "out" are edges i -> j and "in" are edges a -> i; no edge operator is
# inverted. With project = "orthogonal" the exact solution initialises `iters`
# Gauss-Seidel sweeps of the projected block minimiser.
gl_sync_least_squares <- function(G, W, E, N, k, lambda = 1e-6, iters = 20, gauge_idx = 1L, project = c("none","orthogonal")) {
  project <- match.arg(project)
  if (is.null(k)) {
    for (g in G) { k <- nrow(as.matrix(.transform_op(g))); break }
  }
  gauge_idx <- as.integer(gauge_idx)
  if (length(gauge_idx) != 1L || is.na(gauge_idx) || gauge_idx < 1L || gauge_idx > N) {
    stop("gl_sync_least_squares: gauge must index one of the ", N, " domains.", call. = FALSE)
  }
  Tlist <- replicate(N, diag(k), simplify = FALSE)

  # Incident edges per domain, with the role of the domain on the edge
  inc <- vector("list", N)
  for (i in seq_len(N)) inc[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    Gij <- as.matrix(.transform_op(G[[idx]]))
    if (!all(dim(Gij) == c(k, k))) {
      stop("gl_sync_least_squares: edge ", i, "->", j, " operator is ",
           nrow(Gij), "x", ncol(Gij), " but k = ", k, ".", call. = FALSE)
    }
    w <- as.numeric(W[idx] %||% 1)
    inc[[i]][[length(inc[[i]]) + 1L]] <- list(other = j, G = Gij, w = w, out = TRUE)
    inc[[j]][[length(inc[[j]]) + 1L]] <- list(other = i, G = Gij, w = w, out = FALSE)
  }

  # Exact joint solve of the normal equations for the free domains
  free <- setdiff(seq_len(N), gauge_idx)
  pos <- match(seq_len(N), free)
  nf <- length(free)
  if (nf) {
    A <- matrix(0, nf * k, nf * k)
    B <- matrix(0, nf * k, k)
    blk <- function(p) ((p - 1L) * k + 1L):(p * k)
    for (i in free) {
      ri <- blk(pos[i])
      A[ri, ri] <- A[ri, ri] + lambda * diag(k)
      for (nb in inc[[i]]) {
        o <- nb$other
        if (nb$out) {
          A[ri, ri] <- A[ri, ri] + nb$w * diag(k)
          C <- nb$w * nb$G
        } else {
          A[ri, ri] <- A[ri, ri] + nb$w * crossprod(nb$G)
          C <- nb$w * t(nb$G)
        }
        if (o == gauge_idx) {
          B[ri, ] <- B[ri, ] + C %*% Tlist[[o]]
        } else {
          ro <- blk(pos[o])
          A[ri, ro] <- A[ri, ro] - C
        }
      }
    }
    sol <- tryCatch(solve(A, B), error = function(e) {
      stop("gl_sync_least_squares: the synchronisation system is singular (",
           conditionMessage(e), "); check that the dataset graph is connected ",
           "or increase sync$gl$ridge.", call. = FALSE)
    })
    for (i in free) Tlist[[i]] <- sol[blk(pos[i]), , drop = FALSE]
  }
  if (project == "none") iters <- 0L
  if (project == "orthogonal") {
    for (i in free) { sv <- svd(Tlist[[i]]); Tlist[[i]] <- sv$u %*% t(sv$v) }
  }

  for (it in seq_len(iters)) {
    for (i in seq_len(N)) {
      if (i == gauge_idx) next
      Ai <- lambda * diag(k); Bi <- matrix(0, k, k)
      for (nb in inc[[i]]) {
        To <- Tlist[[nb$other]]
        if (nb$out) {
          # edge i -> other: || T_i - G T_other ||^2
          Ai <- Ai + nb$w * diag(k)
          Bi <- Bi + nb$w * (nb$G %*% To)
        } else {
          # edge other -> i: || T_other - G T_i ||^2
          Ai <- Ai + nb$w * crossprod(nb$G)
          Bi <- Bi + nb$w * crossprod(nb$G, To)
        }
      }
      Ti <- tryCatch(solve(Ai, Bi), error = function(e) Tlist[[i]])
      if (project == "orthogonal") { sv <- svd(Ti); Ti <- sv$u %*% t(sv$v) }
      Tlist[[i]] <- Ti
    }
  }
  lapply(seq_len(N), function(i) new_align_transform("GL", Tlist[[i]], from = i, to = "consensus", k = k))
}

# Sizes n_i of each domain implied by the permutation edge operators
# (g_ij is n_j x n_i). Conflicting or missing sizes are errors.
.perm_domain_sizes <- function(G, E, N) {
  n <- rep(NA_integer_, N)
  set_size <- function(d, val, idx) {
    if (is.na(n[d])) {
      n[d] <<- as.integer(val)
    } else if (n[d] != val) {
      stop("permutation sync: edge ", E$i[idx], "->", E$j[idx], " implies ",
           val, " rows for domain ", d, " but another edge implies ", n[d], ".",
           call. = FALSE)
    }
  }
  for (idx in seq_len(nrow(E))) {
    op <- .transform_op(G[[idx]])
    set_size(E$i[idx], ncol(op), idx)
    set_size(E$j[idx], nrow(op), idx)
  }
  if (anyNA(n)) {
    stop("permutation sync: domain(s) ", paste(which(is.na(n)), collapse = ", "),
         " are not covered by any edge of the dataset graph.", call. = FALSE)
  }
  n
}

# Star/hub consensus for permutations: the hub's rows are the consensus and
# every other domain is chained to the hub along a maximum-weight spanning
# tree of the dataset graph, so no direct edge to the hub is required.
permutation_sync_star <- function(G, W, E, N, hub = NULL) {
  n <- .perm_domain_sizes(G, E, N)
  if (is.null(hub)) hub <- which.max(n)
  hub <- as.integer(hub)
  if (length(hub) != 1L || is.na(hub) || hub < 1L || hub > N) {
    stop("permutation_sync_star: hub must index one of the ", N, " domains.", call. = FALSE)
  }
  tr <- lapply(seq_len(nrow(E)), function(idx) {
    .as_edge_transform(G[[idx]], "perm", E$i[idx], E$j[idx])
  })
  w <- as.numeric(W)
  if (length(w) != nrow(E)) w <- rep(1, nrow(E))
  w[!is.finite(w)] <- 0

  T_list <- vector("list", N)
  T_list[[hub]] <- new_align_transform("perm", diag(n[hub]), from = hub, to = "consensus")
  visited <- rep(FALSE, N); visited[hub] <- TRUE
  while (!all(visited)) {
    cand <- which(xor(visited[E$i], visited[E$j]))
    if (!length(cand)) {
      stop("permutation_sync_star: domain(s) ", paste(which(!visited), collapse = ", "),
           " are not connected to hub ", hub, " in the dataset graph.", call. = FALSE)
    }
    idx <- cand[which.max(w[cand])]
    i <- E$i[idx]; j <- E$j[idx]
    if (visited[i]) {
      # j is new: T_j = T_i o g_ji with g_ji the pseudo-inverse of g_ij
      child <- j; parent <- i
      g_cp <- invert_transform(tr[[idx]])
    } else {
      child <- i; parent <- j
      g_cp <- tr[[idx]]
    }
    Tc <- compose_transform(g_cp, T_list[[parent]])
    T_list[[child]] <- new_align_transform("perm", Tc$op, from = child, to = "consensus")
    visited[child] <- TRUE
  }
  T_list
}

# Column-stochastic rounding of a nonnegative affinity A (L x n): every
# domain row (column of A) is assigned unit mass while no consensus row
# receives more than unit mass (alternating Bregman projections).
.perm_sinkhorn_round <- function(A, max_iter = 200, tol = 1e-6) {
  Q <- pmax(A, 0)
  for (it in seq_len(max_iter)) {
    cs <- colSums(Q); cs[cs <= 0] <- 1
    Q <- sweep(Q, 2, cs, "/")
    rs <- rowSums(Q)
    over <- rs > 1
    if (any(over)) Q[over, ] <- Q[over, , drop = FALSE] / rs[over]
    if (all(rs <= 1 + tol) && max(abs(colSums(Q)[colSums(A > 0) > 0] - 1)) < tol) break
  }
  Q
}

.perm_hungarian_round <- function(A) {
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("clue package required for Hungarian rounding", call. = FALSE)
  }
  L <- nrow(A); n <- ncol(A)
  cost <- max(A) - A
  Q <- matrix(0, L, n)
  if (n <= L) {
    m <- as.integer(clue::solve_LSAP(t(cost)))   # consensus row per domain row
    idx <- cbind(m, seq_len(n))
  } else {
    m <- as.integer(clue::solve_LSAP(cost))      # domain row per consensus row
    idx <- cbind(seq_len(L), m)
  }
  # Assignments forced onto zero affinity are not correspondences.
  idx <- idx[A[idx] > 0, , drop = FALSE]
  Q[idx] <- 1
  Q
}

# Spectral consensus for permutations (cycle-aware). The block matrix with
# blocks M[j, i] = w g_ij (n_j x n_i) equals T' T for consistent maps, so its
# leading eigenvectors Y satisfy Y_ref Y_i' proportional to the map from
# domain i to the reference domain. Rounding that product gives T_i, whose
# rows are the reference domain's rows.
permutation_sync_spectral <- function(G, W, E, N, L = NULL,
                                      rounding = c("sinkhorn", "hungarian"),
                                      max_sinkhorn = 200, tol = 1e-6,
                                      ref = NULL) {
  rounding <- match.arg(rounding)
  n <- .perm_domain_sizes(G, E, N)
  if (is.null(L)) L <- max(n)
  L <- as.integer(L)
  offs <- c(0L, cumsum(n)); dim_big <- sum(n)
  if (L < 1L || L > dim_big) {
    stop("permutation_sync_spectral: L must lie in 1..", dim_big, ".", call. = FALSE)
  }
  if (is.null(ref)) ref <- which.max(n)
  ref <- as.integer(ref)

  Ii <- integer(0); Jj <- integer(0); Xx <- numeric(0)
  add_block <- function(ri, rj, B) {
    # as.vector(B) is column-major: row index varies fastest
    Ii <<- c(Ii, rep(ri, times = length(rj)))
    Jj <<- c(Jj, rep(rj, each = length(ri)))
    Xx <<- c(Xx, as.vector(B))
  }
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    Pij <- as.matrix(.transform_op(G[[idx]]))  # n_j x n_i
    w <- as.numeric(W[idx] %||% 1)
    ri <- (offs[i] + 1):(offs[i + 1]); rj <- (offs[j] + 1):(offs[j + 1])
    add_block(rj, ri, w * Pij)
    add_block(ri, rj, w * t(Pij))
  }
  M <- Matrix::sparseMatrix(i = Ii, j = Jj, x = Xx, dims = c(dim_big, dim_big))
  M <- Matrix::forceSymmetric((M + Matrix::t(M)) / 2)

  Y <- .sync_top_eigenvectors(M, L)
  Y_ref <- Y[(offs[ref] + 1):(offs[ref + 1]), , drop = FALSE]

  T_list <- vector("list", N)
  for (i in seq_len(N)) {
    rows <- (offs[i] + 1):(offs[i + 1])
    A <- Y_ref %*% t(Y[rows, , drop = FALSE])  # n_ref x n_i
    # Entries at round-off level carry no correspondence; dropping them keeps
    # domain rows without any path to the reference unassigned.
    A[A <= 1e-8 * max(abs(A))] <- 0
    Q <- if (rounding == "hungarian") {
      .perm_hungarian_round(A)
    } else {
      .perm_sinkhorn_round(A, max_iter = max_sinkhorn, tol = tol)
    }
    T_list[[i]] <- new_align_transform("perm", Q, from = i, to = "consensus")
  }
  T_list
}

permutation_sync <- function(G, W, E, N, max_iter = 100, tol = 1e-8, gauge = NULL) {
  # Synchronise soft maps without a hard hub path using Gauss-Seidel
  # averaging of the consistency relations T_i ~ T_j g_ij (edge i -> j) and
  # T_i ~ T_a g_ai^+ (edge a -> i, pseudo-inverse), with T_gauge = I fixed.
  # After each update every domain row (column of T_i) carries unit mass and
  # no consensus row more than unit mass, so exact (partial) permutations are
  # fixed points even when domains have different sizes.
  n <- .perm_domain_sizes(G, E, N)
  if (is.null(gauge)) gauge <- which.max(n)
  gauge <- as.integer(gauge)
  if (length(gauge) != 1L || is.na(gauge) || gauge < 1L || gauge > N) {
    stop("permutation_sync: gauge must index one of the ", N, " domains.", call. = FALSE)
  }
  m <- n[gauge]

  # Neighbour list: for domain d, op maps d -> neighbour (n_nb x n_d)
  neighbors <- vector("list", N)
  for (i in seq_len(N)) neighbors[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    a <- E$i[idx]; b <- E$j[idx]
    gab <- .as_edge_transform(G[[idx]], "perm", a, b)
    w <- as.numeric(W[idx] %||% 1)
    neighbors[[a]][[length(neighbors[[a]]) + 1L]] <- list(j = b, op = as.matrix(gab$op), w = w)
    neighbors[[b]][[length(neighbors[[b]]) + 1L]] <- list(j = a, op = as.matrix(invert_transform(gab)$op), w = w)
  }

  # Initialize T_i: identity for the gauge, uniform for others
  T <- vector("list", N)
  for (i in seq_len(N)) T[[i]] <- matrix(1 / m, nrow = m, ncol = n[i])
  T[[gauge]] <- diag(m)

  for (iter in seq_len(max_iter)) {
    max_diff <- 0
    for (i in seq_len(N)) {
      if (i == gauge) next
      acc <- matrix(0, nrow = m, ncol = n[i])
      wsum <- 0
      for (nb in neighbors[[i]]) {
        acc <- acc + nb$w * (T[[nb$j]] %*% nb$op)
        wsum <- wsum + nb$w
      }
      if (wsum > 0) acc <- acc / wsum
      new_i <- .perm_sinkhorn_round(acc, max_iter = 50, tol = tol)
      d <- max(abs(new_i - T[[i]]))
      if (is.finite(d) && d > max_diff) max_diff <- d
      T[[i]] <- new_i
    }
    if (max_diff < tol) break
  }

  # Package as transforms to a shared consensus
  lapply(seq_len(N), function(i) new_align_transform("perm", T[[i]], from = i, to = "consensus"))
}

# -- Diagnostics ------------------------------------------------------------

# Edge residual: RMS entry of compose_transform(g_ij, T_j) - T_i (zero when the
# consensus transforms reproduce the edge exactly). Cycle residual for a
# triangle (i, j, k): RMS entry of compose_transform(g_ij, g_jk) - g_ik, with
# edges stored in the opposite orientation inverted first. For "perm" a
# partial map only constrains the rows of the source domain that it maps
# (columns with nonzero mass), so the comparison is restricted to them.
compute_diagnostics_ <- function(group, T_list, G, E, W) {
  tr <- lapply(seq_len(nrow(E)), function(idx) {
    .as_edge_transform(G[[idx]], group, E$i[idx], E$j[idx])
  })
  mapped_cols <- function(...) {
    if (!identical(group, "perm")) return(NULL)
    keep <- NULL
    for (A in list(...)) {
      k <- colSums(abs(as.matrix(A))) > 1e-12
      keep <- if (is.null(keep)) k else keep & k
    }
    keep
  }
  rms_diff <- function(A, B, cols = NULL) {
    A <- as.matrix(A); B <- as.matrix(B)
    if (!all(dim(A) == dim(B))) return(NA_real_)
    if (!is.null(cols)) {
      if (!any(cols)) return(NA_real_)
      A <- A[, cols, drop = FALSE]; B <- B[, cols, drop = FALSE]
    }
    sqrt(mean((A - B)^2))
  }

  edge_res <- rep(NA_real_, nrow(E))
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    edge_res[idx] <- tryCatch({
      pred_i <- compose_transform(tr[[idx]], T_list[[j]])
      rms_diff(pred_i$op, T_list[[i]]$op, mapped_cols(tr[[idx]]$op))
    }, error = function(e) NA_real_)
  }
  edge_df <- data.frame(i = E$i, j = E$j, residual = edge_res)
  edge_wt <- data.frame(i = E$i, j = E$j, weight = as.numeric(W))

  # Cycle residuals on triangles (if present)
  cyc_df <- NULL
  Nmax <- max(c(E$i, E$j, length(T_list)))
  if (nrow(E) >= 3 && Nmax >= 3) {
    key <- paste(E$i, E$j, sep = ":")
    edge_between <- function(a, b) {
      hit <- match(paste(a, b, sep = ":"), key)
      if (!is.na(hit)) return(tr[[hit]])
      hit <- match(paste(b, a, sep = ":"), key)
      if (!is.na(hit)) {
        return(tryCatch(invert_transform(tr[[hit]]), error = function(e) NULL))
      }
      NULL
    }
    adj <- matrix(FALSE, Nmax, Nmax)
    adj[cbind(E$i, E$j)] <- TRUE
    adj <- adj | t(adj)
    tri <- list(i = integer(0), j = integer(0), k = integer(0), residual = numeric(0))
    for (i in seq_len(Nmax - 2L)) {
      for (j in which(adj[i, ])) {
        if (j <= i) next
        for (k in which(adj[i, ] & adj[j, ])) {
          if (k <= j) next
          gij <- edge_between(i, j); gjk <- edge_between(j, k); gik <- edge_between(i, k)
          if (is.null(gij) || is.null(gjk) || is.null(gik)) next
          res <- tryCatch({
            path <- compose_transform(gij, gjk)$op
            rms_diff(path, gik$op, mapped_cols(path, gik$op))
          }, error = function(e) NA_real_)
          tri$i <- c(tri$i, i); tri$j <- c(tri$j, j); tri$k <- c(tri$k, k)
          tri$residual <- c(tri$residual, res)
        }
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
    mean_cycle_residual = if (is.null(cyc_df) || all(is.na(cyc_df$residual))) NA_real_ else mean(cyc_df$residual, na.rm = TRUE),
    entropy = entropy_df
  )
}

# One-pass local refinement for O(k) synchronization: Gauss-Seidel sweep of
# the exact block minimiser of sum w || T_i - g_ij T_j ||^2 over T_i in O(k),
# i.e. T_i = polar(sum_out w g_ij T_j + sum_in w g_ai' T_a).
refine_orthogonal_local_ <- function(T_list, G, W, E, N) {
  k <- ncol(as.matrix(T_list[[1]]$op))
  out <- T_list
  neigh <- vector("list", N)
  for (i in seq_len(N)) neigh[[i]] <- list()
  for (idx in seq_len(nrow(E))) {
    i <- E$i[idx]; j <- E$j[idx]
    Gij <- as.matrix(.transform_op(G[[idx]]))
    w <- as.numeric(W[idx] %||% 1)
    # T_i ~ g_ij T_j  and, for j, the least-squares term g_ij' T_i
    neigh[[i]][[length(neigh[[i]]) + 1L]] <- list(j = j, G = Gij, w = w)
    neigh[[j]][[length(neigh[[j]]) + 1L]] <- list(j = i, G = t(Gij), w = w)
  }
  for (i in seq_len(N)) {
    Mi <- matrix(0, k, k)
    for (nb in neigh[[i]]) {
      Mi <- Mi + nb$w * nb$G %*% as.matrix(out[[nb$j]]$op)
    }
    sv <- tryCatch(svd(Mi), error = function(e) NULL)
    if (!is.null(sv)) {
      Ri <- sv$u %*% t(sv$v)
      out[[i]] <- new_align_transform("O", Ri, from = i, to = "consensus", k = k)
    }
  }
  out
}
