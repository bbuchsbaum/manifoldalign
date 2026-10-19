#' Graph Alignment by Spectral Corresponding Functions (GRASP) - Corrected Version
#'
#' Performs GRASP alignment on hyperdesign data structures. 
#' Projects graph data from multiple domains into aligned spectral spaces.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param preproc Preprocessing function to apply to the data (default: center)
#' @param ncomp Number of spectral components to extract (default: 30)
#' @param q_descriptors Number of diffusion-time descriptors (default: 100)
#' @param sigma Diffusion parameter for descriptor computation (default: 0.73)
#' @param lambda Weight of the descriptor term in the base-alignment
#'   objective \eqn{\mathrm{off}(M^\top \Lambda_2 M) + \lambda
#'   \|F^\top \Phi_1 - G^\top \Phi_2 M\|_F^2} (default: 0.1).
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param graph_metric Distance used to build each domain's kNN graph:
#'   \code{"euclidean"} (default) uses Euclidean distances between rows, so
#'   the graph is invariant to rotations and translations of the features;
#'   \code{"correlation"} first standardises each row across its features
#'   (needs at least three features) and suits rows whose features are
#'   comparable profiles. The heat-kernel bandwidth of each graph is the
#'   median distance to the k-th neighbour under the same metric.
#' @param solver Assignment method: "linear" (default) or its synonym
#'   "hungarian", both the exact Hungarian solver. "auction" is not available.
#' @param max_iter Maximum Riemannian gradient iterations of the base
#'   alignment (default: 200).
#' @param tol Relative-decrease tolerance of the base alignment
#'   (default: 1e-8).
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the GRASP alignment.
#'   Score rows follow each domain's input row order: block 1 is the spectral
#'   basis of domain 1 and block 2 the basis of domain 2 mapped into it by the
#'   base alignment and the diagonal functional map. \code{assignment[i]} is
#'   the domain-2 row matched to domain-1 row \code{i} (\code{NA} when
#'   unmatched because domain 2 has fewer rows). \code{v} holds least-squares
#'   loadings from each domain's preprocessed features to its scores (exact
#'   only for full-row-rank features; otherwise the best linear
#'   approximation). \code{base_alignment} reports the base-alignment
#'   iterations, convergence and objective.
#'
#' @examples
#' \donttest{
#' # Example with hyperdesign graph data
#' library(multidesign)
#' library(tibble)
#' 
#' # Create synthetic graph domains (node coordinates for graph construction)
#' set.seed(123)
#' X1 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 1
#' X2 <- matrix(rnorm(100), 50, 2)  # Node coordinates for domain 2
#' 
#' # Create design data frames with node indices (GRASP doesn't use features)
#' design1 <- tibble(node_id = 1:50)
#' design2 <- tibble(node_id = 1:50)
#' 
#' # Create multidesign objects
#' md1 <- multidesign(X1, design1)
#' md2 <- multidesign(X2, design2)
#' 
#' # Create hyperdesign from multidesign objects
#' hd <- hyperdesign(list(domain1 = md1, domain2 = md2))
#' 
#' # Run GRASP alignment (purely spectral, uses graph structure from X)
#' result <- grasp(hd, ncomp = 20, q_descriptors = 50)
#' }
#'
#' @rdname grasp
#' @method grasp hyperdesign
#' @export
#' @importFrom multivarious init_transform prep concat_pre_processors
#' @importFrom manifoldalign block_indices
#' @importFrom stats dist
grasp.hyperdesign <- function(data, 
                             preproc = multivarious::center(), 
                             ncomp = 30,
                             q_descriptors = 100,
                             sigma = 0.73,
                             lambda = 0.1,
                             use_laplacian = TRUE,
                             solver = "linear",
                             max_iter = 200,
                             tol = 1e-8,
                             graph_metric = c("euclidean", "correlation"),
                             ...) {
  graph_metric <- match.arg(graph_metric)
  
  # Input validation (following package patterns + your suggestions)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda > 0)
  chk::chk_logical(use_laplacian)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter >= 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)

  # Both names select the exact Hungarian solver (clue::solve_LSAP()).
  solver <- .check_lsap_solver(solver, allowed = c("linear", "hungarian"))

  # Validate input data - should be a hyperdesign object containing 2 domains
  chk::chk_s3_class(data, "hyperdesign")
  
  if (!is.list(data) || length(data) == 0) {
    stop("hyperdesign data must contain domains", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("GRASP currently supports exactly 2 domains. For multi-graph ", 
         "alignment with 3+ domains, use grasp_multiset()", call. = FALSE)
  }
  
  # Apply preprocessing per domain. Stateful pre-processors are fitted on each
  # domain and concatenated by feature block for out-of-sample projection.
  domains <- lapply(unclass(data), function(d) {
    list(x = as.matrix(d$x), design = d$design)
  })
  if (is.list(preproc) && !is.function(preproc) &&
      !inherits(preproc, "prepper") && !inherits(preproc, "pre_processor") &&
      length(preproc) == 1L) {
    preproc <- rep(preproc, length(domains))
  }
  prepped <- .fit_domain_preproc(domains, preproc)
  pdata <- prepped$pdata
  proc <- prepped$proc
  names(pdata) <- names(data)

  grasp_fit(pdata, proc, ncomp, q_descriptors, sigma, lambda,
            use_laplacian, solver, feature_block_indices(pdata),
            max_iter = max_iter, tol = tol, graph_metric = graph_metric)
}

#' @keywords internal
grasp_fit <- function(strata, proc, ncomp, q_descriptors, sigma, lambda,
                     use_laplacian, solver, feature_blocks,
                     max_iter = 200, tol = 1e-8, graph_metric = "euclidean") {
  
  # Block A: Spectral basis construction
  bases <- compute_grasp_basis(strata, ncomp, use_laplacian, graph_metric)
  
  # Block B: Multi-scale descriptors  
  descriptors <- compute_grasp_descriptors(bases, q_descriptors, sigma)
  
  # Block C: Base alignment
  alignment_result <- align_grasp_bases(bases[[1]], bases[[2]],
                                       descriptors[[1]], descriptors[[2]], lambda,
                                       max_iter = max_iter, tol = tol)
  
  # Blocks D & E: Functional map and assignment
  assignment_result <- compute_grasp_assignment(bases[[1]], bases[[2]], 
                                               descriptors[[1]], descriptors[[2]],
                                               alignment_result$rotation,
                                               distance_method = "cosine",
                                               solver_method = solver)
  
  # Scores keep each domain's input row order. Domain 2 is expressed in the
  # basis of domain 1 through the base alignment and the diagonal functional
  # map; the node correspondence is reported separately in `assignment`.
  embed1 <- as.matrix(bases[[1]]$vectors)
  embed2 <- as.matrix(bases[[2]]$vectors %*% alignment_result$rotation %*% assignment_result$mapping_matrix)
  scores <- rbind(embed1, embed2)

  # assignment[i] is the domain-2 row matched to domain-1 row i; rows left
  # unmatched when n1 > n2 are NA.
  assignment <- as.integer(assignment_result$assignment)
  assignment[assignment < 1L | assignment > nrow(embed2)] <- NA_integer_

  # Loadings regress each block's aligned scores on its preprocessed features.
  loadings <- .stack_regression_loadings(strata, list(embed1, embed2))

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "grasp",
    extras = list(
      assignment = assignment,
      rotation = alignment_result$rotation,
      mapping_matrix = assignment_result$mapping_matrix,
      base_alignment = alignment_result[setdiff(names(alignment_result), "rotation")]
    )
  )
}

# Rows used to build a GRASP graph. "correlation" standardises each row across
# its features (as adjoin's "normalized" weights do); with fewer than three
# features every standardised row takes at most two values, so it is refused.
.grasp_graph_rows <- function(X, graph_metric) {
  X <- as.matrix(X)
  if (identical(graph_metric, "euclidean")) return(X)
  if (ncol(X) < 3L) {
    stop("graph_metric = \"correlation\" needs at least three features; with ",
         ncol(X), " every standardised row takes at most two values. Use ",
         "graph_metric = \"euclidean\".", call. = FALSE)
  }
  Xc <- X - rowMeans(X)
  s <- sqrt(rowSums(Xc * Xc) / (ncol(X) - 1L))
  s[!is.finite(s) | s == 0] <- 1
  Xc / s
}

# Heat-kernel bandwidth: median distance from each row to its k-th nearest
# neighbour (exact kd-tree search). Deterministic, invariant to row order and
# to isometries of the rows, and independent of the caller's RNG.
.grasp_knn_bandwidth <- function(X, k) {
  k <- as.integer(min(k, nrow(X) - 1L))
  if (k < 1L) return(1)
  nn <- RANN::nn2(X, k = k + 1L, eps = 0)
  d <- nn$nn.dists[, k + 1L]
  d <- d[is.finite(d) & d > 0]
  sigma <- if (length(d)) stats::median(d) else NA_real_
  if (!is.finite(sigma) || sigma <= 0) 1 else sigma
}

#' Block A: Spectral Basis Construction - Corrected
#' 
#' @param strata List of data domains 
#' @param ncomp Number of spectral components
#' @param use_laplacian Whether to use Laplacian normalization
#' @param graph_metric \code{"euclidean"} or \code{"correlation"}; see
#'   \code{\link{grasp}}.
#' @return List of spectral bases for each domain
#' @keywords internal
compute_grasp_basis <- function(strata, ncomp, use_laplacian = TRUE,
                                graph_metric = c("euclidean", "correlation")) {
  graph_metric <- match.arg(graph_metric)
  # Validate input
  if (!is.list(strata) || length(strata) == 0) {
    stop("strata must be a non-empty list", call. = FALSE)
  }
  
  bases <- lapply(seq_along(strata), function(stratum_idx) {
    stratum <- strata[[stratum_idx]]
    # Ensure stratum$x is a matrix for adjoin
    if (!is.matrix(stratum$x)) {
      stratum$x <- as.matrix(stratum$x)
      if (!is.matrix(stratum$x)) {
        stop("Failed to convert domain data to matrix format", call. = FALSE)
      }
    }
    
    # Validate dimensions
    if (nrow(stratum$x) < 3) {
      stop("GRASP requires at least 3 samples per domain", call. = FALSE)
    }
    
    # Heat-kernel kNN graph on the rows of the chosen metric, with a
    # deterministic bandwidth from the same distances.
    X_graph <- .grasp_graph_rows(stratum$x, graph_metric)
    k_graph <- min(max(5, nrow(X_graph) %/% 3), nrow(X_graph) - 1)
    sigma_g <- .grasp_knn_bandwidth(X_graph, k_graph)
    graph_weights <- tryCatch(
      adjoin::graph_weights(
        X_graph,
        weight_mode = "heat",
        neighbor_mode = "knn",
        k = k_graph,
        type = "normal",
        sigma = sigma_g
      ),
      error = function(e) {
        stop("Graph construction failed: ", e$message, call. = FALSE)
      }
    )

    adj_matrix <- adjoin::adjacency(graph_weights)
    adj_matrix <- (adj_matrix + Matrix::t(adj_matrix)) / 2
    adj_matrix <- Matrix::drop0(adj_matrix)

    total_dim <- nrow(adj_matrix)
    if (total_dim <= 1) {
      stop("Graph must contain at least two nodes", call. = FALSE)
    }

    k_use <- if (use_laplacian) min(as.integer(ncomp), total_dim - 1L) else min(as.integer(ncomp), total_dim)
    if (k_use < 1L) {
      stop("Not enough non-trivial eigenvectors available for requested components", call. = FALSE)
    }

    if (use_laplacian) {
      basis <- compute_laplacian_basis(
        adj_matrix,
        k = k_use,
        normalized = TRUE,
        seed = as.integer(1L + stratum_idx * 1009L)
      )
      vectors <- Matrix::Matrix(basis$vectors, sparse = TRUE)
      values <- basis$values
    } else {
      A <- methods::as(Matrix::Matrix(adj_matrix, sparse = TRUE), "dgCMatrix")
      # Largest-algebraic adjacency eigenpairs: RSpectra for moderate k,
      # a dense eigen() when most of the spectrum is requested (PRIMME was used
      # here before and can miss repeated eigenvalues such as the -1 cluster).
      decomp <- .sym_eigs(A, k_use, which = "LA",
                          backend = if (k_use <= ceiling(0.8 * total_dim)) "rspectra" else "dense",
                          tol = 1e-10, maxitr = 1000L,
                          context = "Adjacency eigen-decomposition in GRASP")

      ord <- order(decomp$values, decreasing = TRUE)
      values <- as.numeric(decomp$values)[ord]
      vecs <- as.matrix(decomp$vectors)[, ord, drop = FALSE]
      keep <- which(is.finite(values) & abs(values) > 1e-10)
      if (!length(keep)) stop("No non-trivial adjacency eigenvalues available", call. = FALSE)
      keep <- keep[seq_len(min(k_use, length(keep)))]
      vectors <- Matrix::Matrix(vecs[, keep, drop = FALSE], sparse = TRUE)
      values <- values[keep]
    }

    col_norms <- sqrt(Matrix::colSums(vectors * vectors))
    vectors <- vectors %*% Matrix::Diagonal(x = 1 / pmax(col_norms, 1e-12))

    list(vectors = vectors, values = values)
  })

  # All domains must share one basis width: the base alignment and the
  # functional map are k x k. Truncate to the narrowest basis (graphs with
  # fewer nodes, or with fewer non-trivial eigenpairs, give narrower bases).
  k_common <- min(vapply(bases, function(b) ncol(b$vectors), integer(1)))
  if (k_common < 1L) {
    stop("No non-trivial spectral components available in at least one domain.", call. = FALSE)
  }
  lapply(bases, function(b) {
    if (ncol(b$vectors) == k_common) return(b)
    list(vectors = b$vectors[, seq_len(k_common), drop = FALSE],
         values = b$values[seq_len(k_common)])
  })
}

#' Block B: Multi-scale Descriptors - Corrected
#' 
#' @param bases List of spectral bases
#' @param q_descriptors Number of descriptor functions
#' @param sigma Diffusion parameter
#' @return List of descriptor matrices
#' @keywords internal
compute_grasp_descriptors <- function(bases, q_descriptors, sigma = 0.73) {
  # Validation
  chk::chk_number(q_descriptors)
  chk::chk_true(q_descriptors > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)

  descriptors <- lapply(bases, function(basis) {
    desc_matrix <- compute_hks_descriptors(
      basis = basis,
      q = as.integer(q_descriptors),
      time_mode = "fixed",
      time_range = c(0.1, 50),
      spacing = "linear",
      time_scale = sigma,
      eigen_tol = 1e-12,
      normalize = "column",
      clamp_exponent = -745
    )
    Matrix::Matrix(desc_matrix, sparse = TRUE)
  })

  descriptors
}

#' Block C: Base Alignment
#'
#' Aligns the spectral basis of domain 2 to that of domain 1 with an orthogonal
#' \eqn{k \times k} matrix \eqn{M} minimising the GRASP base-alignment objective
#' \deqn{f(M) = \mathrm{off}(M^\top \Lambda_2 M) +
#'   \lambda \| F^\top \Phi_1 - G^\top \Phi_2 M \|_F^2,}
#' where \eqn{\mathrm{off}(\cdot)} is the squared Frobenius norm of the
#' off-diagonal part, \eqn{\Lambda_2} holds the eigenvalues of domain 2 and
#' \eqn{F, G} are the descriptor matrices. The optimisation starts from the
#' signed permutation that best matches descriptor projections column by column
#' (at which the off-diagonal term is zero) and then takes Riemannian gradient
#' steps on the orthogonal group with a polar retraction and Armijo
#' backtracking. Every accepted step decreases \eqn{f}.
#'
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis
#' @param desc1 First domain's descriptors
#' @param desc2 Second domain's descriptors
#' @param lambda Non-negative weight of the descriptor term. With
#'   \code{lambda = 0} the signed-permutation start is already optimal.
#' @param max_iter Maximum number of gradient iterations (0 returns the
#'   signed-permutation start, reported as not converged).
#' @param tol Convergence tolerance on the relative decrease of \eqn{f}.
#' @return List with \code{rotation} (\eqn{M}), \code{iterations} (gradient
#'   iterations run), \code{converged} (relative decrease below \code{tol}, a
#'   vanishing gradient, or no further descent possible, before
#'   \code{max_iter}), \code{objective} (final \eqn{f}) and
#'   \code{objective_init} (\eqn{f} at the start).
#' @keywords internal
align_grasp_bases <- function(basis1, basis2, desc1, desc2, lambda = 0.1,
                             max_iter = 200, tol = 1e-8) {
  # Validation
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter >= 0)
  max_iter <- as.integer(max_iter)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)

  # Extract matrices
  phi1 <- basis1$vectors
  phi2 <- basis2$vectors

  k <- min(ncol(phi1), ncol(phi2))
  if (!is.finite(k) || k < 1L) stop("Invalid basis dimension in align_grasp_bases.", call. = FALSE)
  k <- as.integer(k)

  if (ncol(phi1) != k) phi1 <- phi1[, seq_len(k), drop = FALSE]
  if (ncol(phi2) != k) phi2 <- phi2[, seq_len(k), drop = FALSE]
  lam2 <- as.numeric(basis2$values)[seq_len(k)]
  if (length(lam2) != k || anyNA(lam2)) {
    stop("basis2$values must provide one eigenvalue per basis vector.", call. = FALSE)
  }

  # Descriptor projections (q x k) used for basis alignment.
  Fm <- as.matrix(Matrix::crossprod(desc1, phi1))
  Gm <- as.matrix(Matrix::crossprod(desc2, phi2))

  # Start: the signed permutation induced by descriptor correlations. It is
  # invariant to eigenvector sign flips and leaves M' Lambda_2 M diagonal.
  Fn <- sqrt(colSums(Fm * Fm) + 1e-12)
  Gn <- sqrt(colSums(Gm * Gm) + 1e-12)
  Fnorm <- sweep(Fm, 2, Fn, "/")
  Gnorm <- sweep(Gm, 2, Gn, "/")

  sim <- abs(crossprod(Fnorm, Gnorm)) # k x k (F columns x G columns)
  cost <- 1 - sim

  perm <- as.integer(solve_lsap_with_padding(cost))

  dots <- crossprod(Fm, Gm)
  sgn <- sign(dots[cbind(seq_len(k), perm)])
  sgn[sgn == 0] <- 1

  M <- matrix(0, k, k)
  M[cbind(perm, seq_len(k))] <- sgn

  objective <- function(M) {
    S <- crossprod(M, lam2 * M)
    diag(S) <- 0
    sum(S * S) + lambda * sum((Fm - Gm %*% M)^2)
  }
  f_cur <- objective(M)
  f_init <- f_cur
  converged <- FALSE
  iterations <- 0L
  armijo_c <- 1e-4

  for (it in seq_len(max_iter)) {
    iterations <- it
    S <- crossprod(M, lam2 * M)
    diag(S) <- 0
    G <- 4 * ((lam2 * M) %*% S) + 2 * lambda * crossprod(Gm, Gm %*% M - Fm)
    MtG <- crossprod(M, G)
    rgrad <- G - M %*% ((MtG + t(MtG)) / 2)
    gnorm2 <- sum(rgrad * rgrad)
    if (!is.finite(gnorm2) || gnorm2 < 1e-24) {
      converged <- TRUE
      break
    }

    step <- 1
    accepted <- FALSE
    for (ls in seq_len(40L)) {
      sv <- svd(M - step * rgrad)
      M_try <- sv$u %*% t(sv$v)
      f_try <- objective(M_try)
      if (is.finite(f_try) && f_try <= f_cur - armijo_c * step * gnorm2) {
        accepted <- TRUE
        break
      }
      step <- step / 2
    }
    if (!accepted) {
      # No descent step left at machine precision: a stationary point.
      converged <- TRUE
      break
    }
    rel <- (f_cur - f_try) / max(abs(f_cur), 1e-12)
    M <- M_try
    f_cur <- f_try
    if (rel < tol) {
      converged <- TRUE
      break
    }
  }

  list(rotation = Matrix::Matrix(M, sparse = FALSE), iterations = iterations,
       converged = converged, objective = f_cur, objective_init = f_init)
}

#' Blocks D & E: Functional Map and Assignment - Corrected
#' 
#' @param basis1 First domain's spectral basis
#' @param basis2 Second domain's spectral basis
#' @param desc1 First domain's descriptors  
#' @param desc2 Second domain's descriptors
#' @param M Rotation matrix from base alignment
#' @param distance_method Distance method: "cosine" or "euclidean"
#' @param solver_method Assignment solver: "linear" or "hungarian" (synonyms
#'   for the exact Hungarian solver; "auction" is not available)
#' @return List with assignment and related matrices
#' @keywords internal
compute_grasp_assignment <- function(basis1, basis2, desc1, desc2, M, 
                                   distance_method = "cosine", 
                                   solver_method = "linear",
                                   alpha = 0.85) {
  # Validation
  if (!distance_method %in% c("cosine", "euclidean")) {
    stop("distance_method must be 'cosine' or 'euclidean'", call. = FALSE)
  }
  
  # "linear" and "hungarian" both select the exact Hungarian solver.
  solver_method <- .check_lsap_solver(solver_method, allowed = c("linear", "hungarian"))
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1]", call. = FALSE)
  }
  
  # Block D: Functional correspondence with ridge regularization (your suggestion)
  # Disconnected or numerically singular graphs can yield different numbers of
  # non-trivial modes in the two domains.  `align_grasp_bases()` constructs M
  # on their common subspace, so use that same subspace here rather than
  # multiplying a wider basis by M.
  k <- nrow(M)
  if (ncol(M) != k || k < 1L || ncol(basis1$vectors) < k ||
      ncol(basis2$vectors) < k) {
    stop("GRASP rotation is incompatible with the spectral bases", call. = FALSE)
  }
  phi1 <- basis1$vectors[, seq_len(k), drop = FALSE]
  phi2 <- basis2$vectors[, seq_len(k), drop = FALSE]
  phi2_aligned <- phi2 %*% M
  
  # Solve F^T*Phi1 ≈ G^T*Phi2_aligned * C for diagonal C
  A_prime <- Matrix::crossprod(desc2, phi2_aligned)  # (q, k)
  B_prime <- Matrix::crossprod(desc1, phi1)         # (q, k)
  
  # Ridge regularization for numerical stability (your suggestion)
  c_diag <- Matrix::colSums(A_prime * B_prime) / (Matrix::colSums(A_prime^2) + 1e-9)
  C <- Matrix::Diagonal(n = length(c_diag), x = c_diag)
  
  # Block E: Final embeddings and assignment
  embed1 <- as.matrix(phi1)
  embed2 <- as.matrix(phi2_aligned %*% C)

  # Augment with descriptor-driven embeddings for robustness (mirrors multiset features)
  B_proj <- as.matrix(desc1 %*% B_prime)           # n1 x k
  A_proj <- as.matrix(desc2 %*% (A_prime %*% C))   # n2 x k
  # Balance scales to avoid one dominating
  fn1 <- sqrt(sum(embed1^2)); fnB <- sqrt(sum(B_proj^2)); s1 <- if (fnB > 0) fn1 / fnB else 1
  fn2 <- sqrt(sum(embed2^2)); fnA <- sqrt(sum(A_proj^2)); s2 <- if (fnA > 0) fn2 / fnA else 1
  embed1 <- cbind(embed1, B_proj * s1)
  embed2 <- cbind(embed2, A_proj * s2)
  
  # Improved distance computation (cosine + small Euclidean blend)
  if (distance_method == "cosine") {
    # Robust cosine computation without external dependencies
    embed1_norms <- sqrt(rowSums(embed1^2) + 1e-12)
    embed2_norms <- sqrt(rowSums(embed2^2) + 1e-12)
    embed1_norm <- embed1 / embed1_norms
    embed2_norm <- embed2 / embed2_norms

    sim_matrix <- embed1_norm %*% t(embed2_norm)
    cost_cos <- 1 - sim_matrix
    cost_cos[cost_cos < 0] <- 0

    # Small Euclidean component to break cosine ties and stabilize
    e1_sq <- rowSums(embed1^2)
    e2_sq <- rowSums(embed2^2)
    dm <- sqrt(pmax(outer(e1_sq, e2_sq, "+") - 2 * (embed1 %*% t(embed2)), 0))
    if (max(dm) > 0) dm <- dm / max(dm)

    cost_matrix <- alpha * cost_cos + (1 - alpha) * dm
  } else {
    # Efficient Euclidean distance computation
    embed1_sq_norms <- rowSums(embed1^2)
    embed2_sq_norms <- rowSums(embed2^2)
    cross_product <- embed1 %*% t(embed2)

    cost_sq <- outer(embed1_sq_norms, embed2_sq_norms, "+") - 2 * cross_product
    cost_matrix <- sqrt(pmax(cost_sq, 0))
  }
  
  assignment <- solve_lsap_with_padding(cost_matrix)

  list(assignment = assignment,
       cost_matrix = cost_matrix,
       mapping_matrix = C)
}

# Exact linear assignment (Hungarian algorithm, clue::solve_LSAP()) for a
# possibly rectangular cost matrix. Returns, for each row, the matched column;
# when nrow > ncol, rows matched to padding columns get indices > ncol.
# `method` is accepted only for the exact solver ("linear", "hungarian",
# "dense" or NULL): clue has no approximate/auction backend.
solve_lsap_with_padding <- function(cost_matrix, method = NULL) {
  if (!is.null(method)) {
    .check_lsap_solver(method, allowed = c("linear", "hungarian", "dense"))
  }
  cost_matrix <- as.matrix(cost_matrix)
  nr <- nrow(cost_matrix)
  nc <- ncol(cost_matrix)
  if (nr == 0 || nc == 0) {
    stop("Cost matrix must have positive dimension", call. = FALSE)
  }
  solve_core <- function(mat) clue::solve_LSAP(mat)

  if (nr == nc) {
    sol <- solve_core(cost_matrix)
    return(as.integer(sol))
  }

  dim_pad <- max(nr, nc)
  pad_value <- max(cost_matrix)
  if (!is.finite(pad_value)) {
    pad_value <- 0
  }
  padded <- matrix(pad_value, nrow = dim_pad, ncol = dim_pad)
  padded[seq_len(nr), seq_len(nc)] <- cost_matrix

  sol <- solve_core(padded)

  as.integer(sol[seq_len(nr)])
}

# Validate an assignment-solver name. Only exact solvers backed by
# clue::solve_LSAP() exist; "auction" was advertised but never implemented by
# the backend, so it is rejected with an explicit message.
.check_lsap_solver <- function(solver, allowed = "linear") {
  if (!is.character(solver) || length(solver) != 1L || is.na(solver)) {
    stop("`solver` must be a single character string.", call. = FALSE)
  }
  if (identical(solver, "auction")) {
    stop("solver = \"auction\" is not available: clue::solve_LSAP() implements ",
         "only the exact Hungarian algorithm. Use solver = \"linear\".", call. = FALSE)
  }
  if (!solver %in% allowed) {
    stop("`solver` should be one of ", paste0("\"", allowed, "\"", collapse = ", "),
         ".", call. = FALSE)
  }
  solver
}

#' @export
grasp <- function(data, ...) {
  UseMethod("grasp")
}

#' @rdname grasp
#' @method grasp default
#' @export
grasp.default <- function(data, ...) {
  stop("grasp() requires a hyperdesign object. ",
       "Got: ", paste(class(data), collapse = ", "), 
       ". See ?grasp for usage examples.", call. = FALSE)
}
