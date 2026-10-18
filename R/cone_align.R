#' CONE-Align: Consensus Optimization for Node Embedding Alignment
#'
#' Performs CONE-Align on hyperdesign data structures. Aligns graph embeddings
#' by alternating between orthogonal transformations and node assignments.
#'
#' CONE-Align tackles graph alignment by iteratively refining both orthogonal
#' transformations and node-to-node correspondences. The algorithm alternates
#' between two steps: (1) finding optimal orthogonal rotations given current
#' assignments via Procrustes analysis, and (2) updating assignments given
#' current transformations via linear assignment.
#'
#' @param data Input data object containing graph domains
#' @param ... Additional arguments passed to specific methods. See 
#'   \code{\link{cone_align.hyperdesign}} for details on method-specific 
#'   parameters such as \code{preproc}, \code{ncomp}, \code{sigma}, 
#'   \code{lambda}, \code{use_laplacian}, \code{solver}, \code{max_iter}, 
#'   and \code{tol}
#'
#' @details
#' CONE-Align operates through the following algorithmic blocks:
#' \itemize{
#'   \item \strong{Spectral Embedding}: Compute Laplacian eigenmaps for each graph
#'   \item \strong{Orthogonal Alignment}: Find rotation matrices via
#'     Procrustes analysis
#'   \item \strong{Assignment Update}: Solve linear assignment problem
#'     for correspondences
#'   \item \strong{Convergence Check}: Iterate until assignments stabilize
#' }
#'
#' The algorithm minimizes the objective:
#' \deqn{\sum_{i,j} ||Q_i^T Z_i - P_{ij} Q_j^T Z_j||_F^2}
#'
#' where \eqn{Z_i} are the embeddings, \eqn{Q_i} are orthogonal transforms,
#' and \eqn{P_{ij}} are permutation matrices.
#'
#' Currently supports *exactly two* domains. For multi-graph alignment with
#' three or more domains, see \code{\link{cone_align_multiple}}.
#'
#' Key parameters:
#' \itemize{
#'   \item \code{ncomp}: Dimension of spectral embeddings
#'   \item \code{sigma}: Bandwidth for embedding computation
#'   \item \code{lambda}: Regularization for numerical stability
#'   \item \code{solver}: Assignment algorithm; only the exact Hungarian
#'     solver (\code{"linear"}) is available
#'   \item \code{knn}: Number of nearest neighbors for graph construction
#' }
#'
#' @return A \code{multiblock_biprojector} object containing:
#' \itemize{
#'   \item \code{s}: Aligned embeddings. Rows follow the input row order of
#'     each domain (domain 1 rows first, then domain 2 rows); domain 1 is
#'     rotated into the frame of domain 2.
#'   \item \code{v}: Loadings for out-of-sample projection, fitted by least
#'     squares so that the preprocessed \eqn{X_i v_i} best approximates the
#'     aligned scores of block \eqn{i} (minimum-norm solution). This is a
#'     linear bridge to a nonlinear spectral embedding: it reproduces the
#'     training scores exactly only when \eqn{X_i} has full row rank, and is
#'     otherwise the best linear approximation.
#'   \item \code{assignment}: Integer vector of length \eqn{n_1};
#'     \code{assignment[i]} is the row of domain 2 matched to row \code{i} of
#'     domain 1. Score rows are not reordered by it.
#'   \item \code{rotation}: Orthogonal transformation matrices
#'     (\code{list(Q, I)}; domain 1 scores are \code{E1 \%*\% Q})
#'   \item \code{converged}, \code{iterations}: Whether the assignment
#'     stabilised before \code{max_iter}, and the number of iterations run
#'   \item \code{sdev}: Standard deviations of aligned components
#'   \item \code{preproc}: The fitted per-domain preprocessor (concatenated by
#'     feature block), or \code{NULL} when \code{preproc = NULL}
#' }
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
#' # Create design data frames with node indices (CONE-Align doesn't use features)
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
#' # Run CONE-Align (purely spectral, uses graph structure from X)
#' result <- cone_align(hd, ncomp = 10)
#' 
#' # Access alignment results
#' node_assignment <- result$assignment
#' aligned_embeddings <- result$s
#' 
#' # Use exact solver for optimal results
#' result_exact <- cone_align(hd, ncomp = 10, solver = "linear")
#' }
#'
#' @references
#' Heimann, M., Shen, H., Safavi, T., & Koutra, D. (2018). REGAL: Representation 
#' learning-based graph alignment. In Proceedings of the 27th ACM International 
#' Conference on Information and Knowledge Management (pp. 117-126).
#'
#' @seealso \code{\link{cone_align.hyperdesign}}
#' @export
cone_align <- function(data, ...) {
  UseMethod("cone_align")
}

#' CONE-Align for Hyperdesign Objects
#'
#' Performs CONE-Align on hyperdesign data structures containing graph domains.
#'
#' @param data A hyperdesign object containing multiple graph domains
#' @param anchors Optional unquoted column name in each domain's design table
#'   containing anchor/correspondence labels. Shared non-NA values indicate
#'   fixed correspondences used to estimate the rotation and constrain the
#'   assignment update.
#' @param preproc A multivarious pre-processor (default: \code{center()}),
#'   a list of one per domain, a plain function, or \code{NULL}. Pre-processors
#'   are fitted separately on each domain and stored in \code{result$preproc}
#'   for out-of-sample projection. A plain function is applied to the training
#'   data only and is not stored, so \code{result$preproc} is then \code{NULL}.
#' @param ncomp Number of embedding dimensions (default: 10)
#' @param sigma Diffusion parameter for embedding computation (default: 0.73)
#' @param lambda Regularization parameter for numerical stability (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment algorithm. Only \code{"linear"} (the exact
#'   Hungarian solver from \code{clue::solve_LSAP()}) is available;
#'   \code{"auction"} was never implemented by the backend and is rejected.
#' @param max_iter Maximum number of iterations (default: 30)
#' @param tol Convergence tolerance for assignment changes (default: 0.01)
#' @param knn Number of nearest neighbors for graph construction
#'   (default: adaptive)
#' @param verbose Logical; if TRUE, report convergence of the alternating
#'   assignment/rotation updates (default: FALSE).
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the CONE-Align results
#'
#' @rdname cone_align
#' @method cone_align hyperdesign
#' @export
#' @importFrom chk chk_number chk_true chk_logical
#' @importFrom multivarious center concat_pre_processors
cone_align.hyperdesign <- function(data, 
                                  anchors = NULL,
                                  preproc = center(), 
                                  ncomp = 10,
                                  sigma = 0.73,
                                  lambda = 0.1,
                                  use_laplacian = TRUE,
                                  solver = "linear",
                                  max_iter = 30,
                                  tol = 0.01,
                                  knn = NULL,
                                  verbose = FALSE,
                                  ...) {
  
  # Fix S3 dispatch - check class first
  if (!inherits(data, "hyperdesign")) {
    stop("data must be an object of class 'hyperdesign'", call. = FALSE)
  }
  
  # Extract domains from hyperdesign object
  domains <- unclass(data)
  
  # Input validation (following exact package patterns)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_logical(use_laplacian)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  chk::chk_flag(verbose)
  solver <- .check_lsap_solver(solver)

  # Validate input data structure
  if (length(domains) == 0) {
    stop("data must be a non-empty hyperdesign object", call. = FALSE)
  }
  
  if (length(domains) != 2) {
    stop("CONE-Align currently supports exactly 2 domains. For multi-graph ", 
         "alignment with 3+ domains, use cone_align_multiple()", call. = FALSE)
  }
  
  # Fit the per-domain preprocessors (multivarious preppers are fitted and
  # stored for out-of-sample projection).
  if (is.function(preproc) && !inherits(preproc, "scaler")) {
    warning("Consider scaling data before embedding as sigma parameter is ",
            "sensitive to data units", call. = FALSE)
  }
  prepped <- .fit_domain_preproc(
    lapply(domains, function(domain) list(x = domain$x, design = domain$design)),
    preproc
  )
  pdata <- prepped$pdata
  proc <- prepped$proc

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
    anchor_info <- resolve_cone_anchor_map(v1, v2)
  }

  # Call CONE-Align fitting function
  cone_align_fit(pdata, proc, ncomp, sigma, lambda, use_laplacian,
                solver, max_iter, tol, feature_block_indices(pdata), knn,
                anchors = anchor_info, verbose = verbose)
}

#' CONE-Align for List of Feature Matrices
#'
#' Performs CONE-Align on a list containing exactly two feature matrices.
#' This provides a simplified interface for direct matrix input without
#' requiring hyperdesign object construction.
#'
#' @param data A list containing exactly two matrices representing node features
#'   for each graph domain. Each matrix should have format: nodes x features.
#' @param preproc See the hyperdesign method.
#' @param ncomp Number of embedding dimensions (default: 10)
#' @param sigma Diffusion parameter for graph construction (default: 0.73)
#' @param lambda Regularization parameter for numerical stability (default: 0.1)
#' @param use_laplacian Whether to use Laplacian normalization (default: TRUE)
#' @param solver Assignment algorithm; see the hyperdesign method.
#' @param max_iter Maximum number of iterations (default: 30)
#' @param tol Convergence tolerance for assignment changes (default: 0.01)
#' @param knn Number of nearest neighbors for graph construction
#'   (default: adaptive)
#' @param verbose Logical; if TRUE, report convergence of the alternating
#'   assignment/rotation updates (default: FALSE).
#' @param ... Additional arguments (currently unused)
#'
#' @return A multiblock_biprojector object containing the CONE-Align results
#'
#' @examples
#' \donttest{
#' # Simple example with two feature matrices
#' set.seed(123)
#' X1 <- matrix(rnorm(100), 50, 2)  # 50 nodes, 2 features
#' X2 <- matrix(rnorm(100), 50, 2)  # 50 nodes, 2 features
#' 
#' # Run CONE-Align directly on matrices
#' result <- cone_align(list(X1, X2), ncomp = 5)
#' 
#' # Access results
#' assignments <- result$assignment
#' embeddings <- result$s
#' }
#'
#' @rdname cone_align
#' @method cone_align list
#' @export
cone_align.list <- function(data,
                           preproc = center(), 
                           ncomp = 10,
                           sigma = 0.73,
                           lambda = 0.1,
                           use_laplacian = TRUE,
                           solver = "linear",
                           max_iter = 30,
                           tol = 0.01,
                           knn = NULL,
                           verbose = FALSE,
                           ...) {
  
  # Validate input list
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty list", call. = FALSE)
  }
  
  if (length(data) != 2) {
    stop("CONE-Align currently supports exactly 2 domains", call. = FALSE)
  }
  if ("anchors" %in% ...names()) {
    stop("`anchors` is only supported for hyperdesign input; ",
         "wrap the matrices in a hyperdesign with an anchor column.", call. = FALSE)
  }

  # Validate that all elements are matrices
  if (!all(vapply(data, is.matrix, logical(1)))) {
    stop("All elements in data list must be matrices", call. = FALSE)
  }
  
  # Validate matrix dimensions
  nrows <- vapply(data, nrow, integer(1))
  ncols <- vapply(data, ncol, integer(1))
  
  if (any(nrows < 3)) {
    stop("Each matrix must have at least 3 rows (nodes)", call. = FALSE)
  }
  
  if (any(ncols < 1)) {
    stop("Each matrix must have at least 1 column (feature)", call. = FALSE)
  }
  
  # Input validation (following exact package patterns)
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(sigma)
  chk::chk_true(sigma > 0)
  chk::chk_number(lambda)
  chk::chk_true(lambda >= 0)
  chk::chk_logical(use_laplacian)
  chk::chk_number(max_iter)
  chk::chk_true(max_iter > 0)
  chk::chk_number(tol)
  chk::chk_true(tol > 0)
  
  chk::chk_flag(verbose)
  solver <- .check_lsap_solver(solver)

  # Convert list of matrices to stratum format (internal structure)
  strata <- lapply(seq_along(data), function(i) {
    list(x = data[[i]], design = data.frame(node_id = seq_len(nrow(data[[i]]))))
  })
  names(strata) <- paste0("domain_", seq_along(strata))

  if (is.function(preproc) && !inherits(preproc, "scaler")) {
    warning("Consider scaling data before embedding as sigma parameter is ",
            "sensitive to data units", call. = FALSE)
  }
  prepped <- .fit_domain_preproc(strata, preproc)
  pdata <- prepped$pdata
  proc <- prepped$proc

  # Call CONE-Align fitting function
  cone_align_fit(pdata, proc, ncomp, sigma, lambda, use_laplacian,
                solver, max_iter, tol, feature_block_indices(pdata), knn,
                verbose = verbose)
}

#' @rdname cone_align
#' @method cone_align default
#' @export
cone_align.default <- function(data, ...) {
  stop("No applicable method for cone_align. data must be a hyperdesign object or list of matrices.", call. = FALSE)
}



#' Compute block indices for multiblock structure
#' @param pdata Preprocessed data list
#' @return Matrix with start and end indices for each block
#' @keywords internal
compute_block_indices <- function(pdata) {
  lens <- vapply(pdata, function(b) nrow(b$x), integer(1))
  cbind(start = c(1L, head(cumsum(lens) + 1L, -1)),
        end   = cumsum(lens))
}

# Fit per-domain preprocessing for the graph-alignment methods (CONE-Align,
# CONE multiple, GRASP, MMA, coupled diagonalization).
#
# `preproc` may be NULL, a multivarious prepper/pre_processor (broadcast to
# every domain and fitted separately on each), a plain function, or a list with
# one entry per domain. Stateful pre-processors are cloned before fitting so a
# caller's object is never mutated. Fitted pre-processors are concatenated by
# *feature* column blocks, which is the index space
# `multivarious::concat_pre_processors()` and `project()` work in. A plain
# function carries no fitted state, so it is applied to the training data only
# and the returned `proc` is NULL (projection of new data cannot replay it).
# Returns list(pdata = domains with transformed `x`, proc = concat pre-processor
# or NULL).
.fit_domain_preproc <- function(domains, preproc) {
  nb <- length(domains)
  pdata <- domains
  if (is.null(preproc)) {
    return(list(pdata = pdata, proc = NULL))
  }

  is_stateful <- function(p) inherits(p, "prepper") || inherits(p, "pre_processor")
  if (is.function(preproc) || is_stateful(preproc)) {
    preproc <- rep(list(preproc), nb)
  } else if (!is.list(preproc) || length(preproc) != nb) {
    stop("`preproc` must be NULL, a pre-processor, a function, or a list with one ",
         "entry per domain (", nb, ").", call. = FALSE)
  }

  proclist <- vector("list", nb)
  stateless <- FALSE
  for (i in seq_len(nb)) {
    pre_i <- preproc[[i]]
    if (is.null(pre_i)) pre_i <- multivarious::pass()
    Xi <- pdata[[i]]$x

    if (is.function(pre_i)) {
      pdata[[i]]$x <- safe_compute(
        pre_i(Xi),
        "Data preprocessing failed. Check preprocessing function or data format."
      )
      stateless <- TRUE
      next
    }
    if (!is_stateful(pre_i)) {
      stop("`preproc[[", i, "]]` must be NULL, a pre-processor or a function.",
           call. = FALSE)
    }

    pre_i <- unserialize(serialize(pre_i, connection = NULL))
    Xi <- as.matrix(Xi)
    if (exists("fit_transform", envir = asNamespace("multivarious"), mode = "function")) {
      ft <- multivarious::fit_transform(pre_i, Xi)
      pdata[[i]]$x <- ft$transformed
      proclist[[i]] <- ft$preproc
    } else {
      templ <- multivarious::prep(pre_i)
      Xi_tf <- multivarious::init_transform(templ, Xi)
      pdata[[i]]$x <- Xi_tf
      proc_attr <- attr(Xi_tf, "preproc")
      proclist[[i]] <- if (is.null(proc_attr)) templ else proc_attr
    }
  }

  if (stateless) {
    return(list(pdata = pdata, proc = NULL))
  }
  names(proclist) <- names(pdata)
  proc <- multivarious::concat_pre_processors(proclist, feature_block_indices(pdata))
  list(pdata = pdata, proc = proc)
}

# Loadings that map (preprocessed) features onto aligned scores.
#
# The spectral graph methods produce scores that are nonlinear functions of the
# features, so there is no exact linear encoder. This returns the
# minimum-norm least-squares bridge `v = X^+ S`, so `X %*% v` is the orthogonal
# projection of the aligned scores `S` onto the column space of `X`: exact when
# `X` has full row rank, otherwise the best linear approximation. Directions
# with singular value below `rel_tol * max(d)` are treated as null.
.score_regression_loadings <- function(X, S, rel_tol = 1e-8) {
  X <- as.matrix(X)
  S <- as.matrix(S)
  if (nrow(X) != nrow(S)) {
    stop("Loadings regression: X has ", nrow(X), " rows but scores have ",
         nrow(S), ".", call. = FALSE)
  }
  out <- matrix(0, ncol(X), ncol(S))
  if (!ncol(X) || !nrow(X)) return(out)
  sv <- svd(X)
  keep <- sv$d > rel_tol * max(sv$d)
  if (!any(keep)) return(out)
  sv$v[, keep, drop = FALSE] %*%
    (crossprod(sv$u[, keep, drop = FALSE], S) / sv$d[keep])
}

# Stack per-block regression loadings for a list of (preprocessed) strata and
# the matching list of aligned score blocks.
.stack_regression_loadings <- function(strata, score_blocks) {
  do.call(rbind, lapply(seq_along(strata), function(i) {
    .score_regression_loadings(strata[[i]]$x, score_blocks[[i]])
  }))
}

#' @keywords internal
cone_align_fit <- function(strata, proc, ncomp, sigma, lambda, use_laplacian,
                          solver, max_iter, tol, feature_blocks, knn,
                          anchors = NULL, verbose = FALSE) {

  # Compute spectral embeddings (following GRASP pattern)
  embeddings <- compute_cone_embeddings(strata, ncomp, sigma, use_laplacian, knn)

  if (any(vapply(embeddings, is.null, logical(1)))) {
    stop("Embedding computation failed for one or more graphs.", call. = FALSE)
  }

  # Iterative alignment (core CONE-Align algorithm)
  alignment_result <- cone_align_iterate(embeddings, solver, max_iter, tol, lambda,
                                         anchors = anchors, verbose = verbose)

  # Scores keep each domain's input row order: graph 1 is rotated into the
  # frame of graph 2, and the correspondence is reported separately in
  # `assignment` (assignment[i] = graph-2 row matched to graph-1 row i).
  embed1_aligned <- embeddings[[1]] %*% alignment_result$Q
  embed2 <- embeddings[[2]]
  scores <- rbind(embed1_aligned, embed2)

  # Loadings regress the aligned scores of each block on its preprocessed
  # features, so X_i v_i lands in the aligned frame.
  loadings <- .stack_regression_loadings(strata, list(embed1_aligned, embed2))

  # Store rotations used to produce scores (graph 1 rotated, graph 2 identity)
  rotation_list <- list(alignment_result$Q,
                        diag(ncol(embed2)))

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "cone_align",
    extras = list(
      assignment = alignment_result$P,
      rotation = rotation_list,
      converged = alignment_result$converged,
      iterations = alignment_result$iterations
    )
  )
}

#' Compute Spectral Embeddings for CONE-Align
#' 
#' @param strata List of data domains
#' @param ncomp Number of embedding dimensions
#' @param sigma Diffusion parameter
#' @param use_laplacian Whether to use Laplacian normalization
#' @param knn Number of nearest neighbors (NULL for adaptive)
#' @return List of embedding matrices
#' @keywords internal
compute_cone_embeddings <- function(strata, ncomp, sigma, use_laplacian, knn) {
  embeddings <- lapply(strata, function(s) {
    tryCatch({
      compute_embedding(s, ncomp, sigma, use_laplacian, knn)
    }, error = function(e) {
      warning(sprintf("Embedding computation failed for a graph: %s", e$message))
      return(NULL) # Return NULL on failure
    })
  })
  
  embeddings
}

#' Worker function to compute a single cone embedding
#' @param domain Data domain
#' @param ncomp Number of embedding dimensions
#' @param sigma Diffusion parameter
#' @param use_laplacian Whether to use Laplacian normalization
#' @param knn Number of nearest neighbors (NULL for adaptive)
#' @return Embedding matrix
#' @keywords internal
compute_embedding <- function(domain, ncomp, sigma, use_laplacian, knn) {
  x <- domain$x
  n_nodes <- nrow(x)
  
  # Dynamically select knn if not provided, based on graph size
  if (is.null(knn)) {
    knn_val <- min(max(ceiling(log(n_nodes)), 3), n_nodes - 1)
  } else {
    knn_val <- min(knn, n_nodes - 1)
  }

  # Construct k-nearest neighbor graph (robust to data type)
  graph_w <- safe_compute(
    adjoin::graph_weights(x, k = knn_val, weight_mode = "heat",
                                   sigma = sigma, neighbor_mode = "knn"),
    "Graph construction failed in CONE-Align"
  )
  A <- adjoin::adjacency(graph_w)
  
  # Compute graph Laplacian (following GRASP patterns)
  L <- graph_laplacian(A, normalized = use_laplacian)
  
  # Robustness: ensure ncomp is not too large for the graph size.
  # We can compute at most n-1 non-trivial eigenvectors.
  # A safe upper bound is n-2 to avoid issues with very small graphs.
  ncomp_safe <- min(ncomp, n_nodes - 2)
  if (ncomp_safe < 1) {
    stop("Cannot compute embeddings. Graph size (", n_nodes, ") is too small ",
         "for the requested number of components (", ncomp, ").", call. = FALSE)
  }
  
  if (ncomp_safe < ncomp) {
    warning("ncomp reduced to ", ncomp_safe, " (was ", ncomp, ") because graph size is ",
            n_nodes, call. = FALSE)
  }
  
  # Eigen-decomposition with robust solver
  # Request one extra eigenvector to account for the trivial one if using Laplacian
  k_request <- min(ncomp_safe + 1, n_nodes - 1)

  # PRIMME (optional) for very small graphs, RSpectra otherwise; see
  # .sym_eigs() for the dense fallbacks and the convergence check.
  decomp <- .sym_eigs(L, k_request, which = "SA",
                      backend = if (n_nodes <= 50) "primme" else "rspectra",
                      tol = 1e-4, maxitr = 300L,
                      context = "Eigen-decomposition in CONE-Align")

  evals <- decomp$values
  vecs <- decomp$vectors
  
  # Select eigenvectors corresponding to smallest non-zero eigenvalues.
  # The first eigenvector for the Laplacian is often trivial (constant vector).
  # This logic robustly handles cases where fewer eigenvectors are returned
  # than requested, or if some eigenvalues are numerically zero.
  nz <- which(evals > 1e-12)  # Indices of all informative eigenvectors
  
  if (length(nz) < ncomp_safe) {
    warning("Graph has only ", length(nz), " informative eigenvectors, ", 
            "but ", ncomp_safe, " were requested. Using all available.", call. = FALSE)
  }
  
  # Select the best available eigenvectors
  ncomp_actual <- min(ncomp_safe, length(nz))
  selected_idx <- nz[seq_len(ncomp_actual)]

  embedding <- vecs[, selected_idx, drop = FALSE]

  # Normalize embeddings for alignment stability (row-wise unit norm + spectral scaling)
  if (!is.null(dim(embedding)) && ncol(embedding) > 0) {
    row_norms <- sqrt(rowSums(embedding^2))
    row_norms[row_norms < 1e-12] <- 1
    embedding <- embedding / row_norms

    if (nrow(embedding) > 0 && ncol(embedding) > 0) {
      spectral_norm <- tryCatch({
        max(svd(embedding, nu = 0, nv = 0)$d)
      }, error = function(e) NA_real_)

      if (is.finite(spectral_norm) && spectral_norm > 1e-12) {
        embedding <- embedding / spectral_norm
      }
    }
  }

  embedding
}

#' @keywords internal
#' @noRd
graph_laplacian <- function(A, normalized = TRUE) {
  # Symmetrise first: the eigensolvers below assume a symmetric operator and
  # read only one triangle, so a directed kNN adjacency would otherwise yield
  # a silently wrong spectrum. Self-loops are kept, as before.
  compute_graph_laplacian(A, normalized = normalized, symmetrize = TRUE,
                          zero_diagonal = FALSE)
}

#' Iterative CONE-Align Algorithm
#' 
#' @param embeddings List of embedding matrices
#' @param solver Assignment solver
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param lambda Regularization parameter
#' @param verbose Logical; report convergence via message()
#' @return List with final rotation, permutation, iteration count and
#'   convergence flag
#' @keywords internal
cone_align_iterate <- function(embeddings, solver, max_iter, tol, lambda, anchors = NULL,
                               verbose = FALSE) {

  if (nrow(embeddings[[1]]) == 0 || nrow(embeddings[[2]]) == 0) {
    stop("Input embeddings to cone_align_iterate must have non-zero rows.")
  }

  if (nrow(embeddings[[1]]) != nrow(embeddings[[2]])) {
    stop("This implementation requires equal node counts. Pad the smaller graph with singleton nodes before alignment.",
         call. = FALSE)
  }

  # Initialize with identity permutation on the reference set
  n1 <- nrow(embeddings[[1]])
  P <- seq_len(n1)

  if (!is.null(anchors)) {
    anchor_map <- anchors$map
    anchor_idx1 <- which(!is.na(anchor_map))
    anchor_idx2 <- as.integer(anchor_map[anchor_idx1])

    if (length(anchor_idx1) > 0) {
      if (any(anchor_idx2 < 1L | anchor_idx2 > nrow(embeddings[[2]]))) {
        stop("Anchor mapping contains indices outside the target graph.", call. = FALSE)
      }
      if (anyDuplicated(anchor_idx2)) {
        stop("Anchor mapping must be one-to-one (duplicate target indices found).", call. = FALSE)
      }

      P[anchor_idx1] <- anchor_idx2

      remaining_src <- setdiff(seq_len(n1), anchor_idx1)
      remaining_tgt <- setdiff(seq_len(n1), anchor_idx2)
      if (length(remaining_tgt) != length(remaining_src)) {
        stop("Anchor mapping size is incompatible with graph sizes.", call. = FALSE)
      }
      P[remaining_src] <- remaining_tgt
    }
  }
  
  # Pre-allocate for efficiency
  P_old <- integer(n1)
  converged <- FALSE
  
  # Iterative alignment
  for (i in 1:max_iter) {
    P_old <- P
    
    # Given P, find Q (orthogonal transformation)
    Q <- solve_procrustes_cone(embeddings[[1]], embeddings[[2]], P, lambda, anchors = anchors)
    
    # Given Q, find P (assignment)
    P <- solve_assignment_cone(embeddings[[1]], embeddings[[2]], Q, solver, anchors = anchors)
    
    # Check for convergence (permutation matrix is stable)
    if (identical(P, P_old)) {
      converged <- TRUE
      break
    }
    
    # Early termination if assignment changes are minimal
    changes <- sum(P != P_old)
    if (changes / n1 < tol) {
      converged <- TRUE
      break
    }
  }
  
  if (isTRUE(verbose)) {
    if (converged) {
      message(sprintf("Converged after %d iterations.", i))
    } else {
      message(sprintf("Maximum iterations (%d) reached without convergence.", max_iter))
    }
  }

  # Return the permutation, rotation, matched indices, and iteration count
  list(P = P, Q = Q, P_matched = P, iterations = i, converged = converged)
}

#' Solve Procrustes Problem for CONE-Align
#' 
#' @param Z1 Source embedding matrix
#' @param Z2 Target embedding matrix
#' @param P Permutation matrix (or vector of indices)
#' @param lambda Regularization parameter for Procrustes problem
#' @return Orthogonal transformation matrix Q
#' @keywords internal
solve_procrustes_cone <- function(Z1, Z2, P, lambda = 0, anchors = NULL) {
  # If anchors are available, estimate the orthogonal map using only the
  # anchored correspondences. This prevents a poor initial assignment from
  # dominating the rotation estimate in symmetric/noisy settings.
  use_anchors <- FALSE
  if (!is.null(anchors) && !is.null(anchors$map)) {
    anchor_map <- anchors$map
    anchor_idx1 <- which(!is.na(anchor_map))
    if (length(anchor_idx1) >= 2L) {
      anchor_idx2 <- as.integer(anchor_map[anchor_idx1])
      use_anchors <- TRUE
    }
  }

  if (use_anchors) {
    Z2_permuted <- Z2[anchor_idx2, , drop = FALSE]
    Z1_use <- Z1[anchor_idx1, , drop = FALSE]
  } else {
    Z2_permuted <- Z2[P, , drop = FALSE]
    Z1_use <- Z1
  }

  if (ncol(Z1_use) != ncol(Z2_permuted) || nrow(Z1_use) != nrow(Z2_permuted)) {
    stop("Procrustes inputs must have identical dimensions.", call. = FALSE)
  }

  M <- crossprod(Z1_use, Z2_permuted)

  # Regularization can bias the estimate toward the identity. When anchors are
  # present, rely on the anchored correspondences directly.
  if (lambda > 0 && !use_anchors) {
    diag(M) <- diag(M) + lambda
  }

  sv <- svd(M)
  Q <- sv$u %*% t(sv$v)

  Q
}

#' Solve Assignment Problem for CONE-Align
#' 
#' @param Z1 First embedding matrix
#' @param Z2 Second embedding matrix
#' @param Q Orthogonal transformation matrix
#' @param solver Assignment algorithm (only "linear", the exact Hungarian solver)
#' @return Vector of assignment indices
#' @keywords internal
solve_assignment_cone <- function(Z1, Z2, Q, solver, anchors = NULL) {
  
  solver <- .check_lsap_solver(solver)

  # Pre-compute transformed embeddings to avoid repeated multiplication
  Z1_transformed <- Z1 %*% Q
  
  # Compute squared distance matrix
  D <- pairwise_sqdist(Z1_transformed, Z2)
  
  # Solve assignment (following exact GRASP dependency pattern)
  # Check if clue is installed, if not, suggest installation
  if (!requireNamespace("clue", quietly = TRUE)) {
    stop("Package 'clue' needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  n1 <- nrow(Z1)
  
  anchor_map <- NULL
  if (!is.null(anchors) && !is.null(anchors$map)) {
    anchor_map <- anchors$map
  }

  # Solve the assignment problem (optionally with hard anchors).
  if (is.null(anchor_map) || all(is.na(anchor_map))) {
    assignment <- clue::solve_LSAP(D)
    return(as.integer(assignment))
  }

  anchor_idx1 <- which(!is.na(anchor_map))
  anchor_idx2 <- as.integer(anchor_map[anchor_idx1])
  if (any(anchor_idx2 < 1L | anchor_idx2 > nrow(Z2))) {
    stop("Anchor mapping contains indices outside the target graph.", call. = FALSE)
  }
  if (anyDuplicated(anchor_idx2)) {
    stop("Anchor mapping must be one-to-one (duplicate target indices found).", call. = FALSE)
  }

  free_i <- setdiff(seq_len(n1), anchor_idx1)
  free_j <- setdiff(seq_len(nrow(Z2)), anchor_idx2)

  P <- integer(n1)
  P[anchor_idx1] <- anchor_idx2
  if (!length(free_i)) return(P)

  if (length(free_i) != length(free_j)) {
    stop("Anchor mapping size is incompatible with graph sizes.", call. = FALSE)
  }

  D_sub <- D[free_i, free_j, drop = FALSE]
  assignment <- clue::solve_LSAP(D_sub)

  P[free_i] <- free_j[as.integer(assignment)]
  P
}

#' Resolve anchor correspondences for CONE-Align
#'
#' Returns a mapping from source indices to target indices using label vectors
#' where shared non-NA values indicate anchored correspondences.
#'
#' @param vec1 Anchor labels for source domain (length n1)
#' @param vec2 Anchor labels for target domain (length n2)
#' @return A list with `map` (integer vector length n1, NA for unanchored)
#'   and `values` of shared anchors; or NULL if no anchors found.
#' @keywords internal
resolve_cone_anchor_map <- function(vec1, vec2) {
  idx1 <- which(!is.na(vec1))
  idx2 <- which(!is.na(vec2))
  if (length(idx1) == 0 || length(idx2) == 0) return(NULL)

  vals1 <- vec1[idx1]
  vals2 <- vec2[idx2]

  if (anyDuplicated(vals1)) {
    stop("Anchor labels must be unique within the source domain.", call. = FALSE)
  }
  if (anyDuplicated(vals2)) {
    stop("Anchor labels must be unique within the target domain.", call. = FALSE)
  }

  common <- intersect(vals1, vals2)
  if (!length(common)) {
    stop("No shared anchor labels found across domains.", call. = FALSE)
  }

  i_pos <- match(common, vals1)
  j_pos <- match(common, vals2)
  idx1 <- idx1[i_pos]
  idx2 <- idx2[j_pos]

  map <- rep(NA_integer_, length(vec1))
  map[idx1] <- idx2
  list(map = map, values = common)
}

#' Compute Pairwise Squared Distances
#' 
#' Compute squared Euclidean distances between rows of two matrices
#' 
#' @param X First matrix (n x d)
#' @param Y Second matrix (m x d)
#' @return A matrix of squared Euclidean distances (n x m)
#' @keywords internal
#' @noRd
pairwise_sqdist <- function(X, Y = NULL) {
  if (is.null(Y)) {
    Y <- X
  }
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  X_norm_sq <- rowSums(X^2)
  Y_norm_sq <- rowSums(Y^2)

  D <- outer(X_norm_sq, Y_norm_sq, "+") - 2 * tcrossprod(X, Y)

  pmax(D, 0)
}
