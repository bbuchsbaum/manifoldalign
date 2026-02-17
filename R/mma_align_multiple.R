#' Multiset Manifold Alignment (MMA)
#'
#' Aligns 3+ domains by: (i) building spectral embeddings (commute-time or
#' Laplacian), (ii) aligning eigenbases via eigensignatures (histograms +
#' Hungarian), and (iii) refining per-domain alignments with an EM-based
#' probabilistic point-registration on the orthogonal group.
#'
#' Two matching modes are supported:
#' - match_to = "reference": pairwise EM to a chosen reference domain.
#' - match_to = "consensus": joint EM that learns a shared template (consensus centers)
#'   along with per-domain orthogonal transforms; supports unequal node counts naturally.
#'
#' The implementation reuses shared helpers in this package:
#' - neighborweights::graph_weights + adjacency for graph building
#' - graph_laplacian (internal) for Laplacian construction
#' - safe_compute, compute_block_indices, feature_block_indices, new_alignment_result
#'
#' @param data A hyperdesign object containing three or more domains
#' @param ref_idx Integer; which domain acts as reference (default 1)
#' @param preproc Preprocessing function (default multivarious::center())
#' @param ncomp Number of components; set NULL to choose automatically
#' @param sigma Gaussian affinity bandwidth (in data units)
#' @param knn Number of nearest neighbors for graph construction (default adaptive)
#' @param embedding One of "ctd" (commute-time; default) or "lap"
#' @param normalize One of "hypersphere" (default) or "none" for embeddings
#' @param hist_bins Integer number of bins for eigensignature histograms; NULL for auto
#' @param hist_max_bins Clamp for auto-bins (default 128)
#' @param hist_similarity One of "cor" or "cosine" (default "cor")
#' @param signature_method One of "hist", "w2", or "hybrid" for eigenvector matching
#' @param eig_penalty Eigenvalue penalty weight for signature alignment (default 0)
#' @param signature_gate_multiplicity Logical; gate matches by eigenvalue clusters
#' @param signature_gate_tau Threshold for eigenvalue clustering (NULL for auto)
#' @param signature_retry_if_weak Logical; retry with alternative method if weak match
#' @param em_max_iter Integer EM iterations (default 50)
#' @param em_tol Numeric EM tolerance on relative change (default 1e-5)
#' @param em_sigma0 Optional initial stddev in embedding units; NULL for auto
#' @param em_outlier Prior mass for uniform outlier component (default 0.01)
#' @param force_SO Logical; if TRUE, force det(R)=+1 (default FALSE)
#' @param target_var Lower bound criterion for auto K when ncomp=NULL (default 0.95)
#' @param match_to One of "reference" (default) or "consensus"
#' @param consensus_centers Strategy to set number of consensus centers when
#'   match_to="consensus": "ref" (use ref domain size) or "min" (use minimum size)
#' @param n_centers Optional integer overriding consensus_centers; must be <= min sizes
#' @param consensus_init Initialization for consensus centers: "ref" or "kmeans"
#' @param final_assignment Optional post-processing: "none" (default), "hungarian",
#'   or "argmax". Affects only extras$assignment; core EM remains soft.
#' @param ... Unused.
#'
#' @return A multiblock alignment object with stacked scores, block loadings,
#'   and extras (rotations, posteriors, histogram alignment info, and optionally
#'   consensus and final assignment).
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(30), 10, 3)
#' X2 <- matrix(rnorm(30), 10, 3)
#' X3 <- matrix(rnorm(30), 10, 3)
#' data_list <- list(X1, X2, X3)
#' result <- mma_align_multiple(data_list, ncomp = 2)
#' dim(result$scores)
#' }
#' @export
#' @importFrom chk chk_number chk_true chk_logical
#' @importFrom multivarious center concat_pre_processors
#' @importFrom Matrix crossprod
mma_align_multiple.hyperdesign <- function(
  data,
  ref_idx = 1L,
  preproc = center(),
  ncomp = 10,
  sigma = 0.73,
  knn = NULL,
  embedding = c("ctd", "lap"),
  normalize = c("hypersphere", "none"),
  hist_bins = NULL,
  hist_max_bins = 128L,
  hist_similarity = c("cor", "cosine"),
  signature_method = c("hist","w2","hybrid"),
  eig_penalty = 0.0,
  signature_gate_multiplicity = FALSE,
  signature_gate_tau = NULL,
  signature_retry_if_weak = FALSE,
  em_max_iter = 50L,
  em_tol = 1e-5,
  em_sigma0 = NULL,
  em_outlier = 0.01,
  force_SO = FALSE,
  target_var = 0.95,
  match_to = c("reference", "consensus"),
  consensus_centers = c("ref", "min"),
  n_centers = NULL,
  consensus_init = c("ref", "kmeans"),
  final_assignment = c("none", "hungarian", "argmax"),
  ...
){
  # Validate top-level inputs
  if (!inherits(data, "hyperdesign")) {
    stop("data must be a 'hyperdesign' object", call. = FALSE)
  }
  domains <- unclass(data)
  if (length(domains) < 3) {
    stop("mma_align_multiple() requires three or more domains.", call. = FALSE)
  }

  chk::chk_number(ref_idx); chk::chk_true(ref_idx >= 1 && ref_idx <= length(domains))
  if (!is.null(ncomp)) { chk::chk_number(ncomp); chk::chk_true(ncomp > 0) }
  chk::chk_number(sigma); chk::chk_true(sigma > 0)
  chk::chk_number(hist_max_bins); chk::chk_true(hist_max_bins >= 16)
  chk::chk_number(em_max_iter); chk::chk_true(em_max_iter > 0)
  chk::chk_number(em_tol); chk::chk_true(em_tol > 0)
  chk::chk_number(em_outlier); chk::chk_true(em_outlier >= 0 && em_outlier < 1)

  embedding        <- match.arg(embedding)
  normalize        <- match.arg(normalize)
  hist_similarity  <- match.arg(hist_similarity)
  signature_method <- match.arg(signature_method)
  match_to         <- match.arg(match_to)
  consensus_centers<- match.arg(consensus_centers)
  consensus_init   <- match.arg(consensus_init)
  final_assignment <- match.arg(final_assignment)

  # Preprocess per-domain (reuse cone_align pattern)
  pdata <- lapply(domains, function(domain) {
    x_pre <- safe_compute(
      if (is.function(preproc)) preproc(domain$x) else domain$x,
      "Data preprocessing failed. Check preprocessing function or data format."
    )

    if (is.function(preproc) && !inherits(preproc, "scaler")) {
      warning("Consider scaling data before embedding as sigma is in data units", call. = FALSE)
    }
    list(x = x_pre, design = domain$design)
  })
  names(pdata) <- names(domains)

  # Block indices and concat preproc descriptor
  block_idx <- compute_block_indices(pdata)
  proc <- if (is.function(preproc)) {
    proclist <- setNames(rep(list(preproc), length(pdata)), names(pdata))
    block_idx_list <- split(block_idx, row(block_idx))
    safe_compute(multivarious::concat_pre_processors(proclist, block_idx_list),
                 "Preprocessing concatenation failed")
  } else NULL

  mma_align_multiple_fit(
    strata            = pdata,
    proc              = proc,
    ref_idx           = ref_idx,
    ncomp             = ncomp,
    sigma             = sigma,
    knn               = knn,
    embedding         = embedding,
    normalize         = normalize,
    hist_bins         = hist_bins,
    hist_max_bins     = hist_max_bins,
    hist_similarity   = hist_similarity,
    signature_method  = signature_method,
    eig_penalty       = eig_penalty,
    signature_gate_multiplicity = signature_gate_multiplicity,
    signature_gate_tau = signature_gate_tau,
    signature_retry_if_weak = signature_retry_if_weak,
    em_max_iter       = em_max_iter,
    em_tol            = em_tol,
    em_sigma0         = em_sigma0,
    em_outlier        = em_outlier,
    force_SO          = force_SO,
    target_var        = target_var,
    match_to          = match_to,
    consensus_centers = consensus_centers,
    n_centers         = n_centers,
    consensus_init    = consensus_init,
    final_assignment  = final_assignment,
    feature_blocks    = feature_block_indices(pdata)
  )
}

#' @rdname mma_align_multiple
#' @export
mma_align_multiple.list <- function(data, ...) {
  if (!is.list(data) || length(data) < 3) {
    stop("Provide a list with at least three matrices.", call. = FALSE)
  }
  if (!all(vapply(data, is.matrix, logical(1)))) {
    stop("Every element in the list must be a numeric matrix.", call. = FALSE)
  }
  nrows <- vapply(data, nrow, integer(1))
  ncols <- vapply(data, ncol, integer(1))
  if (any(nrows < 3)) stop("Each matrix must have >= 3 rows (nodes)", call. = FALSE)
  if (any(ncols < 1)) stop("Each matrix must have >= 1 column (feature)", call. = FALSE)

  strata <- lapply(seq_along(data), function(i) {
    list(x = data[[i]], design = data.frame(node_id = seq_len(nrow(data[[i]]))))
  })
  names(strata) <- paste0("domain_", seq_along(strata))
  class(strata) <- "hyperdesign"
  mma_align_multiple.hyperdesign(strata, ...)
}

#' @export
mma_align_multiple.default <- function(data, ...) {
  stop("No applicable method: data must be a 'hyperdesign' or list of matrices.", call. = FALSE)
}

# ---- Core fit ----

mma_align_multiple_fit <- function(
  strata, proc, ref_idx,
  ncomp, sigma, knn,
  embedding, normalize,
  hist_bins, hist_max_bins, hist_similarity,
  signature_method, eig_penalty,
  signature_gate_multiplicity, signature_gate_tau, signature_retry_if_weak,
  em_max_iter, em_tol, em_sigma0, em_outlier, force_SO, target_var,
  match_to, consensus_centers, n_centers, consensus_init, final_assignment,
  feature_blocks
){
  m <- length(strata)

  # 1) Spectral embeddings
  embeddings <- compute_mma_embeddings(
    strata = strata,
    ncomp = ncomp,
    sigma = sigma,
    knn = knn,
    embedding = embedding,
    normalize = normalize,
    target_var = target_var
  )

  # 2) Eigenbasis alignment to reference (histograms + Hungarian)
  aligns <- vector("list", m)
  R0     <- vector("list", m)
  refU   <- attr(embeddings[[ref_idx]], "U_raw")
  refDeg <- attr(embeddings[[ref_idx]], "deg")
  refLam <- attr(embeddings[[ref_idx]], "lam_keep")
  K      <- ncol(refU)
  for (g in seq_len(m)) {
    if (g == ref_idx) {
      aligns[[g]] <- list(P = diag(K), S = diag(K))
      R0[[g]]     <- diag(K)
    } else {
      ag <- eigensignature_align(
        Uref = refU,
        Ug   = attr(embeddings[[g]], "U_raw"),
        deg_ref = refDeg,
        deg_g   = attr(embeddings[[g]], "deg"),
        lam_ref = refLam,
        lam_g   = attr(embeddings[[g]], "lam_keep"),
        method = signature_method,
        bins = hist_bins,
        max_bins = hist_max_bins,
        similarity = hist_similarity,
        eig_penalty = eig_penalty,
        gate_multiplicity = signature_gate_multiplicity,
        gate_tau = signature_gate_tau,
        retry_if_weak = signature_retry_if_weak
      )
      R0[[g]]     <- ag$S %*% ag$P
      aligns[[g]] <- ag
    }
  }

  if (match_to == "reference") {
    refE <- embeddings[[ref_idx]]
    rotations  <- vector("list", m)
    posteriors <- vector("list", m)
    em_stats   <- vector("list", m)
    rotations[[ref_idx]]  <- diag(K)
    posteriors[[ref_idx]] <- matrix(NA_real_, nrow = nrow(refE), ncol = nrow(refE))
    em_stats[[ref_idx]]   <- list(iterations = 0L, converged = TRUE, sigma2 = NA_real_, loglik_last = NA_real_)

    for (g in seq_len(m)) {
      if (g == ref_idx) next
      EM <- em_mma_pairwise_R(
        X_ref    = refE,
        X_g      = embeddings[[g]],
        R_init   = R0[[g]],
        sigma2_0 = if (is.null(em_sigma0)) NULL else (em_sigma0^2),
        outlier  = em_outlier,
        max_iter = em_max_iter,
        tol      = em_tol,
        force_SO = force_SO
      )
      rotations[[g]]  <- EM$R
      posteriors[[g]] <- EM$alpha
      em_stats[[g]]   <- list(iterations = EM$iterations, converged = EM$converged,
                               sigma2 = EM$sigma2, loglik_last = EM$loglik_last)
    }

    aligned <- lapply(seq_len(m), function(g) embeddings[[g]] %*% t(rotations[[g]]))
    scores  <- do.call(rbind, aligned)

    assignment <- build_final_assignment(final_assignment, posteriors, mode = "reference")

    loadings <- do.call(rbind, lapply(seq_len(m), function(i) {
      xi <- strata[[i]]$x
      Matrix::crossprod(xi, embeddings[[i]])
    }))

    return(new_alignment_result(
      scores = scores,
      loadings = loadings,
      preproc = proc,
      feature_blocks = feature_blocks,
      subclass = "mma_align_multiple",
      extras = list(
        rotations  = rotations,
        posteriors = posteriors,
        hist_align = aligns,
        signature_align = aligns,
        ref_idx    = ref_idx,
        embedding  = embedding,
        normalize  = normalize,
        K          = K,
        assignment = assignment,
        em_stats   = em_stats
      )
    ))
  }

  # consensus branch
  aligned0 <- lapply(seq_len(m), function(g) embeddings[[g]] %*% t(R0[[g]]))
  n_list   <- vapply(aligned0, nrow, integer(1))
  n0       <- determine_n_centers(n_list, consensus_centers, n_centers, ref_idx)

  C0 <- init_consensus_centers(aligned0, ref_idx = ref_idx, n_centers = n0, method = consensus_init)

  EMc <- em_mma_consensus(
    X_list   = embeddings,
    R_init   = R0,
    C_init   = C0,
    sigma2_0 = if (is.null(em_sigma0)) NULL else (em_sigma0^2),
    outlier  = em_outlier,
    max_iter = em_max_iter,
    tol      = em_tol,
    force_SO = force_SO
  )

  rotations  <- EMc$R_list
  C          <- EMc$C
  posteriors <- EMc$alpha_list
  em_stats   <- list(iterations = EMc$iterations, converged = EMc$converged,
                     sigma2 = EMc$sigma2, loglik_last = EMc$loglik_last)

  aligned <- lapply(seq_len(m), function(g) embeddings[[g]] %*% t(rotations[[g]]))
  scores  <- do.call(rbind, aligned)

  assignment <- build_final_assignment(final_assignment, posteriors, mode = "consensus")

  loadings <- do.call(rbind, lapply(seq_len(m), function(i) {
    xi <- strata[[i]]$x
    Matrix::crossprod(xi, embeddings[[i]])
  }))

  new_alignment_result(
    scores = scores,
    loadings = loadings,
    preproc = proc,
    feature_blocks = feature_blocks,
    subclass = "mma_align_multiple",
    extras = list(
      rotations   = rotations,
      posteriors  = posteriors,
      consensus   = C,
      hist_align  = aligns,
      signature_align = aligns,
      ref_idx     = ref_idx,
      embedding   = embedding,
      normalize   = normalize,
      K           = K,
      n_centers   = n0,
      assignment  = assignment,
      em_stats    = em_stats
    )
  )
}

# ---------- Spectral embedding pipeline ----------

compute_mma_embeddings <- function(strata, ncomp, sigma, knn, embedding, normalize, target_var) {
  {
    # First pass: if ncomp is NULL, compute a common K across domains (min of per-domain Ks)
    if (is.null(ncomp)) {
      Kc <- vapply(strata, function(s) {
        X <- s$x; n <- nrow(X)
        knn_val <- if (is.null(knn)) min(max(ceiling(log(n)), 3), max(n - 1, 1)) else min(knn, max(n - 1, 1))
        gw <- safe_compute(
          neighborweights::graph_weights(X, k = knn_val, weight_mode = "heat", sigma = sigma, neighbor_mode = "knn"),
          "Graph construction failed in MMA (K selection)"
        )
        A <- neighborweights::adjacency(gw)
        L <- graph_laplacian(A, normalized = FALSE)
        max(2L, min(choose_K_theta_min(L, target_var = target_var), max(n - 2L, 1L)))
      }, integer(1))
      K_final <- max(2L, min(Kc))
    } else {
      # fixed K across all domains
      K_final <- max(2L, as.integer(ncomp))
    }

    # Second pass: compute embeddings with K_final for every domain
    embeds <- lapply(strata, function(s) {
      X <- s$x
      n <- nrow(X)
      knn_val <- if (is.null(knn)) min(max(ceiling(log(n)), 3), max(n - 1, 1)) else min(knn, max(n - 1, 1))
      gw <- safe_compute(
        neighborweights::graph_weights(X, k = knn_val, weight_mode = "heat", sigma = sigma, neighbor_mode = "knn"),
        "Graph construction failed in MMA"
      )
      A <- neighborweights::adjacency(gw)
      L <- graph_laplacian(A, normalized = FALSE)

      K <- max(2L, min(K_final, max(n - 2L, 1L)))
      k_req <- min(K + 1L, max(n - 1L, 1L))

      decomp <- NULL
      if (n <= 50 && requireNamespace("PRIMME", quietly = TRUE)) {
        decomp <- safe_compute(PRIMME::eigs_sym(L, NEig = k_req, which = "SA"), "Eigen-decomposition failed in MMA")
      } else if (requireNamespace("RSpectra", quietly = TRUE)) {
        decomp <- safe_compute(RSpectra::eigs_sym(L, k = k_req, which = "SA", opts = list(tol = 1e-4, maxitr = 300)),
                               "Eigen-decomposition failed in MMA")
      } else {
        decomp <- safe_compute(PRIMME::eigs_sym(L, NEig = k_req, which = "SA"), "Eigen-decomposition failed in MMA")
      }

      lam <- as.numeric(decomp$values)
      U   <- decomp$vectors
      # Ensure ascending eigenvalue order and consistent eigenvector ordering
      if (length(lam) > 1) {
        o <- order(lam, decreasing = FALSE)
        lam <- lam[o]
        if (!is.null(dim(U))) U <- U[, o, drop = FALSE]
      }
      nz <- which(lam > 1e-12)
      if (length(nz) < K) {
        warning("Graph has only ", length(nz), " informative eigenvectors; requested ", K, ". Using all available.", call. = FALSE)
      }
      keep <- nz[seq_len(min(K, length(nz)))]
      lam_keep <- lam[keep]
      U_keep   <- U[, keep, drop = FALSE]

      Xemb <- if (embedding == "ctd") {
        sweep(U_keep, 2L, sqrt(lam_keep), "/")
      } else {
        U_keep
      }

      if (normalize == "hypersphere" && ncol(Xemb) > 0) {
        r <- sqrt(rowSums(Xemb^2)); r[r == 0] <- 1
        Xemb <- Xemb / r
      }

      # Attach attributes for signature alignment
      attr(Xemb, "U_raw")   <- U_keep
      attr(Xemb, "lam_keep")<- lam_keep
      # Degrees (row sums of adjacency) for degree-weighted signatures
      deg <- as.numeric(Matrix::rowSums(A))
      deg[!is.finite(deg)] <- 0
      if (sum(deg) <= 0) deg[] <- 1
      attr(Xemb, "deg") <- deg
      Xemb
    })

    # Ensure a common embedding dimension across domains.
    # Disconnected graphs can yield multiple zero eigenvalues, reducing the count of informative eigenvectors.
    K_dom <- vapply(embeds, ncol, integer(1))
    K_common <- min(K_dom)
    if (K_common < max(K_dom)) {
      warning("MMA: using K = ", K_common, " components across all domains (some graphs had fewer informative eigenvectors).",
              call. = FALSE)
    }
    if (K_common <= 0) {
      stop("MMA: no informative eigenvectors were found in at least one domain.", call. = FALSE)
    }

    trim_embed <- function(E, K) {
      if (ncol(E) <= K) return(E)
      U_raw <- attr(E, "U_raw")
      lam_keep <- attr(E, "lam_keep")
      deg <- attr(E, "deg")
      E2 <- E[, seq_len(K), drop = FALSE]
      if (!is.null(U_raw)) attr(E2, "U_raw") <- U_raw[, seq_len(K), drop = FALSE]
      if (!is.null(lam_keep)) attr(E2, "lam_keep") <- lam_keep[seq_len(K)]
      attr(E2, "deg") <- deg
      E2
    }
    lapply(embeds, trim_embed, K = K_common)
  }
}

choose_K_theta_min <- function(L, target_var = 0.95) {
  n <- nrow(L)
  Kcap <- max(2L, min(64L, max(n - 2L, 1L)))
  k_req <- min(Kcap + 1L, max(n - 1L, 1L))
  decomp <- NULL
  if (n <= 50 && requireNamespace("PRIMME", quietly = TRUE)) {
    decomp <- safe_compute(PRIMME::eigs_sym(L, NEig = k_req, which = "SA"), "Eigen-decomposition failed in MMA (K)")
  } else if (requireNamespace("RSpectra", quietly = TRUE)) {
    decomp <- safe_compute(RSpectra::eigs_sym(L, k = k_req, which = "SA", opts = list(tol = 1e-4, maxitr = 300)),
                           "Eigen-decomposition failed in MMA (K)")
  } else {
    decomp <- safe_compute(PRIMME::eigs_sym(L, NEig = k_req, which = "SA"), "Eigen-decomposition failed in MMA (K)")
  }
  lam <- as.numeric(decomp$values)
  # Ensure ascending order and filter out near-zero eigenvalues
  lam <- sort(lam[lam > 1e-12], decreasing = FALSE)
  if (length(lam) <= 2) return(max(2L, min(length(lam), 2L)))

  inv <- 1 / lam
  for (K in 2:min(Kcap, length(lam) - 1L)) {
    num <- sum(inv[1:K])
    denom <- sum(inv[1:(K-1)]) + (n - K + 1L) * inv[K]
    theta_min <- num / denom
    if (theta_min >= target_var) return(K)
  }
  min(Kcap, length(lam) - 1L)
}

# ---------- Eigensignature alignment (histograms + Hungarian) ----------

eigensignature_align <- function(Uref, Ug,
                                 deg_ref = NULL, deg_g = NULL,
                                 lam_ref = NULL, lam_g = NULL,
                                 method = c("hist","w2","hybrid"),
                                 bins = NULL, max_bins = 128L,
                                 similarity = c("cor","cosine"),
                                 eig_penalty = 0.0,
                                 gate_multiplicity = FALSE,
                                 gate_tau = NULL,
                                 retry_if_weak = FALSE) {
  method <- match.arg(method)
  similarity <- match.arg(similarity)
  K <- ncol(Uref)
  if (ncol(Ug) != K) stop("Incompatible eigen-basis sizes in eigensignature_align")

  # Helper: histogram similarity path (original behavior)
  hist_path <- function() {
    n <- nrow(Uref)
    if (is.null(bins)) {
      bins_loc <- max(32L, min(max_bins, as.integer(ceiling(n^(1/3)))))
    } else bins_loc <- bins
    br <- seq(-1, 1, length.out = bins_loc + 1L)

    h1  <- lapply(seq_len(K), function(k) hist(Uref[,k], breaks = br, plot = FALSE)$counts)
    h2  <- lapply(seq_len(K), function(k) hist(Ug[,k],   breaks = br, plot = FALSE)$counts)
    h2m <- lapply(seq_len(K), function(k) hist(-Ug[,k],  breaks = br, plot = FALSE)$counts)

    simfun <- switch(similarity,
                     cor    = function(a,b) if (stats::sd(a) == 0 || stats::sd(b) == 0) 0 else stats::cor(a,b),
                     cosine = function(a,b) {
                       na <- sqrt(sum(a*a)); nb <- sqrt(sum(b*b)); if (na == 0 || nb == 0) 0 else sum(a*b)/(na*nb)
                     })

    A <- matrix(0, K, K); B <- matrix(0, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      s1 <- simfun(h1[[i]], h2[[j]]); s2 <- simfun(h1[[i]], h2m[[j]])
      if (s1 >= s2) { A[i,j] <- s1; B[i,j] <- +1 } else { A[i,j] <- s2; B[i,j] <- -1 }
    }
    # Convert similarity to cost and solve
    cost <- max(A) - A
    perm <- clue::solve_LSAP(cost)
    P <- diag(K)[, as.integer(perm), drop = FALSE]
    svec <- vapply(seq_len(K), function(i) B[i, as.integer(perm)[i]], numeric(1))
    S <- diag(svec)
    # Confidence in [0,1] by rescaling A
    minA <- suppressWarnings(min(A, na.rm = TRUE)); maxA <- suppressWarnings(max(A, na.rm = TRUE))
    if (!is.finite(minA)) minA <- 0; if (!is.finite(maxA)) maxA <- 1
    A_scaled <- if ((maxA - minA) > 1e-12) (A - minA) / (maxA - minA) else matrix(0.5, K, K)
    pair_scores <- vapply(seq_len(K), function(i) A[i, as.integer(perm)[i]], numeric(1))
    pair_conf <- vapply(seq_len(K), function(i) A_scaled[i, as.integer(perm)[i]], numeric(1))
    list(P = P, S = S, scores = A, signs = svec, perm = as.integer(perm),
         pair_scores = pair_scores, pair_confidence = pair_conf,
         method = "hist")
  }

  # Helper: build W2 cost matrix with degree weights and optional eigen penalty
  # Build gating penalty based on eigenvalue clusters (optional)
  gate_info <- NULL
  if (isTRUE(gate_multiplicity) && !is.null(lam_ref) && !is.null(lam_g)) {
    cluster_eigs <- function(lam, tau = NULL) {
      lam <- as.numeric(lam)
      if (length(lam) <= 1) return(rep(1L, length(lam)))
      dif <- diff(lam)
      if (is.null(tau)) {
        md <- stats::median(abs(dif[is.finite(dif)]))
        if (!is.finite(md) || md <= 0) md <- mean(abs(dif[is.finite(dif)]))
        if (!is.finite(md) || md <= 0) md <- 1e-6
        tau <- 1e-3 * md
      }
      grp <- integer(length(lam)); grp[1] <- 1L; g <- 1L
      for (i in 2:length(lam)) {
        if (abs(lam[i] - lam[i-1]) > tau) g <- g + 1L
        grp[i] <- g
      }
      grp
    }
    ref_grp <- cluster_eigs(lam_ref, gate_tau)
    g_grp   <- cluster_eigs(lam_g,   gate_tau)
    gate_info <- list(ref = ref_grp, tgt = g_grp)
  }

  w2_path <- function(hybrid = FALSE, A_hist_raw = NULL) {
    n1 <- nrow(Uref); n2 <- nrow(Ug)
    wr <- if (is.null(deg_ref)) rep(1/n1, n1) else as.numeric(deg_ref)
    wg <- if (is.null(deg_g))   rep(1/n2, n2) else as.numeric(deg_g)
    if (!all(is.finite(wr))) wr[!is.finite(wr)] <- 0
    if (!all(is.finite(wg))) wg[!is.finite(wg)] <- 0
    sr <- sum(wr); sg <- sum(wg)
    if (sr <= 0) wr[] <- 1/n1 else wr <- wr / sr
    if (sg <= 0) wg[] <- 1/n2 else wg <- wg / sg

    pre <- function(U, w) {
      lapply(seq_len(ncol(U)), function(k) {
        o <- order(U[,k])
        list(x = U[o,k], w = w[o])
      })
    }
    # Prefer fast C++ kernel; fallback to R implementation on error
    C_w2 <- NULL; sign_w2 <- NULL
    cw2 <- tryCatch({
      w2_cost_matrix_1d(Uref, Ug, wr, wg)
    }, error = function(e) NULL)
    if (!is.null(cw2) && is.list(cw2) && !is.null(cw2$C) && !is.null(cw2$sign)) {
      C_w2 <- as.matrix(cw2$C)
      sign_w2 <- as.matrix(cw2$sign)
    } else {
      # R fallback
      Rpre <- pre(Uref, wr)
      Gpre <- pre(Ug,   wg)
      Gpre_m <- lapply(Gpre, function(g) {
        nx <- length(g$x)
        list(x = -rev(g$x), w = rev(g$w))
      })
      w2_pair <- function(a, b) {
        i <- 1L; j <- 1L
        wi <- a$w[i]; wj <- b$w[j]
        xi <- a$x[i]; xj <- b$x[j]
        res <- 0.0
        na <- length(a$x); nb <- length(b$x)
        while (i <= na && j <= nb) {
          m <- if (wi < wj) wi else wj
          dx <- xi - xj
          res <- res + (dx * dx) * m
          wi <- wi - m; wj <- wj - m
          if (wi <= 1e-15) { i <- i + 1L; if (i <= na) { xi <- a$x[i]; wi <- a$w[i] } }
          if (wj <= 1e-15) { j <- j + 1L; if (j <= nb) { xj <- b$x[j]; wj <- b$w[j] } }
        }
        res
      }
      C_w2 <- matrix(0, K, K)
      sign_w2 <- matrix(1, K, K)
      for (k in seq_len(K)) for (l in seq_len(K)) {
        c1 <- w2_pair(Rpre[[k]], Gpre[[l]])
        c2 <- w2_pair(Rpre[[k]], Gpre_m[[l]])
        cmin <- if (c1 <= c2) { sign_w2[k,l] <- +1; c1 } else { sign_w2[k,l] <- -1; c2 }
        C_w2[k,l] <- cmin
      }
    }
    # Optional eigenvalue penalty (small nudge)
    if (!is.null(lam_ref) && !is.null(lam_g) && eig_penalty > 0) {
      for (k in seq_len(K)) for (l in seq_len(K)) {
        if (length(lam_ref) >= k && length(lam_g) >= l) {
          C_w2[k,l] <- C_w2[k,l] + eig_penalty * (abs(lam_ref[k] - lam_g[l]) / (lam_ref[k] + lam_g[l] + 1e-12))
        }
      }
    }

    if (!hybrid) {
      # Apply gating as large penalties for cross-cluster matches
      if (!is.null(gate_info)) {
        # penalty scale
        Cscale <- stats::median(C_w2[is.finite(C_w2)])
        if (!is.finite(Cscale) || Cscale <= 0) Cscale <- max(1, mean(C_w2[is.finite(C_w2)]))
        big <- 10 * Cscale + 1
        for (k in seq_len(K)) for (l in seq_len(K)) {
          if (gate_info$ref[k] != gate_info$tgt[l]) C_w2[k,l] <- C_w2[k,l] + big
        }
      }
      perm <- clue::solve_LSAP(C_w2)  # minimize cost
      P <- diag(K)[, as.integer(perm), drop = FALSE]
      svec <- vapply(seq_len(K), function(i) sign_w2[i, as.integer(perm)[i]], numeric(1))
      S <- diag(svec)
      # Similarity from cost for confidence
      medC <- stats::median(C_w2[is.finite(C_w2)])
      if (!is.finite(medC) || medC <= 0) medC <- mean(C_w2[is.finite(C_w2) & C_w2 > 0])
      if (!is.finite(medC) || medC <= 0) medC <- 1.0
      Sim_w2 <- exp(-C_w2 / medC)
      pair_scores <- vapply(seq_len(K), function(i) (-C_w2)[i, as.integer(perm)[i]], numeric(1))
      pair_conf <- vapply(seq_len(K), function(i) Sim_w2[i, as.integer(perm)[i]], numeric(1))
      return(list(P = P, S = S, scores = -C_w2, signs = svec, perm = as.integer(perm),
                  pair_scores = pair_scores, pair_confidence = pair_conf,
                  cost_matrix = C_w2, sign_matrix = sign_w2, method = "w2"))
    } else {
      # Hybrid: combine W2 similarity with histogram similarity
      Cn <- C_w2
      medC <- stats::median(Cn)
      if (!is.finite(medC) || medC <= 0) medC <- mean(Cn[is.finite(Cn) & Cn > 0])
      if (!is.finite(medC) || medC <= 0) medC <- 1.0
      Sim_w2 <- exp(-Cn / medC)

      A_hist <- A_hist_raw
      # Rescale histogram similarity to [0,1] for blending
      minA <- suppressWarnings(min(A_hist, na.rm = TRUE)); maxA <- suppressWarnings(max(A_hist, na.rm = TRUE))
      if (!is.finite(minA)) minA <- 0
      if (!is.finite(maxA)) maxA <- 1
      A_hist_scaled <- if ((maxA - minA) > 1e-12) (A_hist - minA) / (maxA - minA) else matrix(0.5, K, K)

      A_hybrid <- 0.7 * Sim_w2 + 0.3 * A_hist_scaled
      # Convert to cost to apply gating easily: cost = max - sim
      cost_h <- max(A_hybrid) - A_hybrid
      if (!is.null(gate_info)) {
        big <- 10 * stats::median(cost_h[is.finite(cost_h)]) + 1
        if (!is.finite(big) || big <= 0) big <- 10
        for (k in seq_len(K)) for (l in seq_len(K)) {
          if (gate_info$ref[k] != gate_info$tgt[l]) cost_h[k,l] <- cost_h[k,l] + big
        }
      }
      # Minimize cost
      perm <- clue::solve_LSAP(cost_h)
      P <- diag(K)[, as.integer(perm), drop = FALSE]
      svec <- vapply(seq_len(K), function(i) sign_w2[i, as.integer(perm)[i]], numeric(1))
      S <- diag(svec)
      pair_scores <- vapply(seq_len(K), function(i) A_hybrid[i, as.integer(perm)[i]], numeric(1))
      pair_conf <- vapply(seq_len(K), function(i) A_hybrid[i, as.integer(perm)[i]], numeric(1))
      return(list(P = P, S = S, scores = A_hybrid, signs = svec, perm = as.integer(perm),
                  pair_scores = pair_scores, pair_confidence = pair_conf,
                  cost_matrix = C_w2, sign_matrix = sign_w2, method = "hybrid"))
    }
  }

  if (method == "hist") {
    res <- hist_path()
  } else if (method == "hybrid") {
    # Build histogram similarity matrix A as in original path
    n <- nrow(Uref)
    bins_loc <- if (is.null(bins)) max(32L, min(max_bins, as.integer(ceiling(n^(1/3))))) else bins
    br <- seq(-1, 1, length.out = bins_loc + 1L)
    h1  <- lapply(seq_len(K), function(k) hist(Uref[,k], breaks = br, plot = FALSE)$counts)
    h2  <- lapply(seq_len(K), function(k) hist(Ug[,k],   breaks = br, plot = FALSE)$counts)
    h2m <- lapply(seq_len(K), function(k) hist(-Ug[,k],  breaks = br, plot = FALSE)$counts)
    simfun <- switch(similarity,
                     cor    = function(a,b) if (stats::sd(a) == 0 || stats::sd(b) == 0) 0 else stats::cor(a,b),
                     cosine = function(a,b) { na <- sqrt(sum(a*a)); nb <- sqrt(sum(b*b)); if (na == 0 || nb == 0) 0 else sum(a*b)/(na*nb) })
    A <- matrix(0, K, K)
    for (i in seq_len(K)) for (j in seq_len(K)) {
      s1 <- simfun(h1[[i]], h2[[j]]); s2 <- simfun(h1[[i]], h2m[[j]]);
      A[i,j] <- max(s1, s2)
    }
    res <- w2_path(hybrid = TRUE, A_hist_raw = A)
  } else {
    res <- w2_path(hybrid = FALSE)
  }

  # Optional weak-score retry: compare with alternative method and keep the better
  if (isTRUE(retry_if_weak)) {
    # Compute assignment score sum for current res
    score_sum <- sum(res$scores[cbind(seq_len(K), res$perm)])
    # Heuristic threshold based on matrix stats
    sc <- res$scores
    mu <- mean(sc); sdv <- stats::sd(c(sc))
    weak <- is.finite(score_sum) && (score_sum < (mu - sdv))
    if (weak) {
      alt <- if (method == "w2") "hist" else "w2"
      alt_res <- eigensignature_align(Uref, Ug, deg_ref, deg_g, lam_ref, lam_g,
                                      method = alt, bins = bins, max_bins = max_bins,
                                      similarity = similarity, eig_penalty = eig_penalty,
                                      gate_multiplicity = gate_multiplicity,
                                      gate_tau = gate_tau, retry_if_weak = FALSE)
      alt_sum <- sum(alt_res$scores[cbind(seq_len(K), alt_res$perm)])
      if (is.finite(alt_sum) && alt_sum > score_sum) res <- alt_res
    }
  }

  res
}

# ---------- EM: pairwise to reference ----------

em_mma_pairwise_R <- function(X_ref, X_g, R_init, sigma2_0 = NULL, outlier = 0.01,
                              max_iter = 50L, tol = 1e-5, force_SO = FALSE) {
  nr <- nrow(X_ref); ng <- nrow(X_g); K <- ncol(X_ref)
  if (ncol(X_g) != K) stop("Dimension mismatch in EM pairwise")

  R <- R_init
  # init sigma^2 from rough residuals
  if (is.null(sigma2_0)) {
    RX <- X_ref %*% t(R)
    muY <- colMeans(X_g); muRX <- colMeans(RX)
    sigma2 <- max((sum(rowSums((X_g  - rep(1, ng) %o% muRX)^2)) / (ng * K) +
                   sum(rowSums((RX   - rep(1, nr) %o% muY )^2)) / (nr * K)) / 2, 1e-8)
  } else sigma2 <- max(sigma2_0, 1e-8)

  pi_out <- outlier
  pi_in  <- (1 - pi_out) / nr
  c_out  <- (pi_out / pi_in)
  obj_prev <- -Inf
  alpha <- matrix(0, nr, ng)
  ll <- NA_real_

  for (it in seq_len(max_iter)) {
    RX <- X_ref %*% t(R)              # nr x K
    YtRX <- X_g %*% t(RX)             # ng x nr
    YY <- rowSums(X_g^2)              # ng
    RR <- rowSums(RX^2)               # nr
    D2 <- outer(RR, YY, "+") - 2 * t(YtRX)   # nr x ng

    G <- exp(-D2 / (2 * sigma2))      # nr x ng
    denom <- colSums(G) + c_out
    alpha <- sweep(G, 2L, denom, "/")

    M <- crossprod(X_ref, alpha %*% X_g)   # K x K
    sv <- svd(M)
    Rnew <- sv$u %*% t(sv$v)
    if (force_SO && det(Rnew) < 0) { sv$u[,K] <- -sv$u[,K]; Rnew <- sv$u %*% t(sv$v) }
    R <- Rnew

    num <- sum(alpha * D2)
    den <- K * sum(alpha)
    sigma2 <- max(num / max(den, 1e-12), 1e-12)

    ll <- sum(log(colSums(G) + c_out))
    if (is.finite(obj_prev)) {
      rel <- abs(ll - obj_prev) / (abs(obj_prev) + 1e-8)
      if (rel < tol) break
    }
    obj_prev <- ll
  }
  list(R = R, alpha = alpha, sigma2 = sigma2, iterations = it, converged = (it < max_iter), loglik_last = ll)
}

# ---------- EM: consensus (joint template) ----------

em_mma_consensus <- function(X_list, R_init, C_init, sigma2_0 = NULL,
                             outlier = 0.01, max_iter = 50L, tol = 1e-5, force_SO = FALSE) {
  m <- length(X_list)
  K <- ncol(X_list[[1]])
  n0 <- nrow(C_init)
  R_list <- R_init
  C <- C_init

  sigma2 <- if (is.null(sigma2_0)) {
    res <- 0; cnt <- 0
    for (g in seq_len(m)) {
      RX <- C %*% t(R_list[[g]])
      muX <- colMeans(X_list[[g]]); muC <- colMeans(RX)
      res <- res + sum((muX - muC)^2); cnt <- cnt + 1
    }
    max(res / max(cnt, 1) / K, 1e-6)
  } else max(sigma2_0, 1e-8)

  pi_out <- outlier
  obj_prev <- -Inf
  alpha_list <- vector("list", m)
  ll <- NA_real_

  for (it in seq_len(max_iter)) {
    denom_sum <- 0.0
    num_sigma <- 0.0
    sum_alpha <- 0.0

    # E-step
    for (g in seq_len(m)) {
      Xg <- X_list[[g]]
      Rg <- R_list[[g]]
      RC <- C %*% t(Rg)            # n0 x K
      YY <- rowSums(Xg^2)          # n_g
      CC <- rowSums(RC^2)          # n0
      D2 <- outer(CC, YY, "+") - 2 * (RC %*% t(Xg))
      G <- exp(-D2 / (2 * sigma2))
      c_out <- (pi_out / ((1 - pi_out) / n0))
      denom <- colSums(G) + c_out
      alpha <- sweep(G, 2L, denom, "/")
      alpha_list[[g]] <- alpha
      num_sigma <- num_sigma + sum(alpha * D2)
      sum_alpha <- sum_alpha + sum(alpha)
      denom_sum <- denom_sum + sum(log(denom))
    }

    # M-step: rotations
    for (g in seq_len(m)) {
      Xg <- X_list[[g]]
      alpha <- alpha_list[[g]]
      M <- crossprod(C, alpha %*% Xg)
      sv <- svd(M)
      Rg <- sv$u %*% t(sv$v)
      if (force_SO && det(Rg) < 0) { sv$u[,K] <- -sv$u[,K]; Rg <- sv$u %*% t(sv$v) }
      R_list[[g]] <- Rg
    }

    # M-step: consensus centers C
    numC <- matrix(0, nrow = n0, ncol = K)
    denC <- Reduce("+", lapply(alpha_list, rowSums))
    for (g in seq_len(m)) {
      Xg <- X_list[[g]]; Rg <- R_list[[g]]
      numC <- numC + alpha_list[[g]] %*% (Xg %*% Rg)
    }
    denC[denC == 0] <- 1
    C <- numC / denC

    sigma2 <- max(num_sigma / (K * sum_alpha), 1e-12)
    ll <- denom_sum
    if (is.finite(obj_prev)) {
      rel <- abs(ll - obj_prev) / (abs(obj_prev) + 1e-8)
      if (rel < tol) break
    }
    obj_prev <- ll
  }

  list(R_list = R_list, C = C, alpha_list = alpha_list, sigma2 = sigma2,
       iterations = it, converged = (it < max_iter), loglik_last = ll)
}

# ---------- Helpers: centers and final assignments ----------

determine_n_centers <- function(n_list, strategy = c("ref","min"), n_centers = NULL, ref_idx = 1L) {
  strategy <- match.arg(strategy)
  if (!is.null(n_centers)) {
    nc <- as.integer(n_centers)
    if (nc <= 0 || nc > min(n_list)) {
      stop("n_centers must be in 1..min(domain sizes).", call. = FALSE)
    }
    return(nc)
  }
  if (strategy == "ref") return(n_list[[ref_idx]])
  min(n_list)
}

init_consensus_centers <- function(aligned0, ref_idx, n_centers, method = c("ref","kmeans")) {
  method <- match.arg(method)
  if (method == "ref") {
    Eref <- aligned0[[ref_idx]]
    if (nrow(Eref) == n_centers) return(Eref)
    km <- stats::kmeans(Eref, centers = n_centers, nstart = 5)
    return(km$centers)
  } else {
    Xstack <- do.call(rbind, aligned0)
    km <- stats::kmeans(Xstack, centers = n_centers, nstart = 5)
    return(km$centers)
  }
}

build_final_assignment <- function(final_assignment, posteriors, mode = c("reference","consensus")) {
  mode <- match.arg(mode)
  if (final_assignment == "none") return(NULL)

  if (mode == "reference") {
    if (final_assignment == "argmax") {
      return(lapply(posteriors, function(A) if (all(is.na(A))) NA_integer_ else max.col(t(A), ties.method = "first")))
    } else {
      return(lapply(posteriors, function(A) if (all(is.na(A))) NA_integer_ else hungarian_from_posterior(A)))
    }
  } else {
    if (final_assignment == "argmax") {
      return(lapply(posteriors, function(A) max.col(t(A), ties.method = "first")))
    } else {
      return(lapply(posteriors, hungarian_from_posterior))
    }
  }
}

hungarian_from_posterior <- function(A) {
  nr <- nrow(A); nc <- ncol(A); s <- max(nr, nc)
  base <- if (length(A)) max(A, na.rm = TRUE) else 0
  if (!is.finite(base)) base <- 0
  C <- matrix(base, s, s)
  C[1:nr, 1:nc] <- base - A
  perm <- clue::solve_LSAP(C)
  res <- rep(NA_integer_, nc)
  for (j in seq_len(min(s, nc))) {
    if (perm[j] <= nr) res[j] <- perm[j] else res[j] <- NA_integer_
  }
  res
}
