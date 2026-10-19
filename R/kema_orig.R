#' Original Kernel Manifold Alignment (KEMA)
#'
#' `kema_orig()` is a deprecated duplicate of [kema()]. New code should call
#' `kema()`, which runs this same fidelity-gated original formulation.
#' `kema_orig()` is scheduled for removal in manifoldalign 1.0.0.
#'
#' Paper-faithful implementation of Tuia & Camps-Valls (2016) using the
#' generalized eigenproblems:
#'
#' Full KEMA (Eq. 6):
#' \deqn{K (L + \mu L_s) K \Lambda = \lambda K L_d K \Lambda}
#'
#' REKEMA (Eq. 10, reduced-rank):
#' \deqn{K_{rn} (L + \mu L_s) K_{nr} \Lambda = \lambda K_{rn} L_d K_{nr} \Lambda}
#'
#' This function intentionally omits extension layers (regression solver,
#' repulsion term \eqn{L_r}, fail-soft retries, and kernel centering variants).
#'
#' @param data Input object. A `hyperdesign` or `multidesign` object.
#' @param y Label column used for class graphs.
#' @param ... Method-specific arguments.
#' @return A `multiblock_biprojector` object with KEMA embeddings.
#' @references
#' Tuia, D., & Camps-Valls, G. (2016). Kernel manifold alignment for domain
#' adaptation. PLoS ONE, 11(2), e0148655.
#' @export
kema_orig <- function(data, y, ...) {
  .Deprecated(
    "kema",
    package = "manifoldalign",
    msg = paste0(
      "kema_orig() is deprecated because it duplicates the public KEMA API. ",
      "Use kema(); it now runs the same fidelity-gated original formulation. ",
      "kema_orig() is scheduled for removal in manifoldalign 1.0.0."
    )
  )
  UseMethod("kema_orig")
}

#' @rdname kema_orig
#' @method kema_orig multidesign
#' @param subject Subject/domain column for splitting a multidesign object.
#' @param preproc Preprocessing recipe (default `center()`).
#' @param ncomp Number of components.
#' @param knn Number of nearest neighbors for within-domain graph.
#' @param sigma Kernel graph scale for local geometry.
#' @param mu Weight for class-pull Laplacian in \eqn{L + \mu L_s}.
#' @param kernel Kernel function.
#' @param sample_frac Fraction of samples retained as landmarks (<1 => REKEMA).
#'   The selected landmark rows are stored in `fit$encoder$landmarks`.
#' @param seed Optional seed for REKEMA landmark selection; `NULL` draws from
#'   the session RNG. A supplied seed leaves the caller's RNG stream untouched.
#' @param lambda Positive Tikhonov regularization added to the RHS generalized
#'   matrix. A positive value is required because the partial generalized
#'   eigensolver uses that matrix as an inner-product metric.
#' @param backend Backend for the generalized eigensolver. One of `"auto"`,
#'   `"full_exact"`, `"reduced_exact"`, or `"operator_exact"`.
#' @param backend_control Optional list for backend auto-selection and fidelity
#'   diagnostics. Supported keys: `full_exact_max_n`, `reduced_exact_max_r`,
#'   `dense_exact_max_dim`, `dense_fallback_max_dim`,
#'   `fidelity_residual_tol`, `fidelity_orth_tol`, `fidelity_action`,
#'   `eigencore_tol`, `eigencore_maxit`, `eigencore_seed`, `kernel_rank_tol`,
#'   and `max_dense_bytes`. Quotients of dimension at most
#'   `dense_exact_max_dim` (default 400) are solved by the dense generalized
#'   eigensolver; larger ones by shifted-Cholesky preconditioned LOBPCG
#'   (`eigencore_maxit`, `eigencore_tol`). A certified partial solve whose
#'   pairs miss the original-pencil residual target is repeated with
#'   `eigencore_tol` tightened by 100 (at most twice). When no partial solve
#'   succeeds, `full_exact` and `reduced_exact` retry with the dense
#'   generalized eigensolver if the dimension is at most
#'   `dense_fallback_max_dim` (default 2000); `operator_exact` only if it is
#'   at most `dense_exact_max_dim`. Every returned solve must pass the
#'   eigencore certificate and the fidelity gate. `eigencore_seed` (default 1) seeds the iterative
#'   partial eigensolver locally, so fits do not depend on, or change, the
#'   caller's RNG state. Every backend materialises the dense kernel blocks,
#'   their eigen factors and dense label-comparison matrices
#'   (O(sum n_i r_i + n_labelled^2)); `full_exact` and `reduced_exact` also
#'   form dense O(r^2) quotient matrices, which `operator_exact` replaces by
#'   matrix-free products and per-domain O(r_i^2) preconditioner blocks.
#'   `max_dense_bytes` (default 8 GiB) bounds a coarse per-backend estimate of
#'   that footprint: candidates that exceed it are skipped before allocation
#'   and the fit errors when none fits. Use `Inf` to disable the guard.
#' @export
kema_orig.multidesign <- function(data, y,
                                  subject,
                                  preproc = center(),
                                  ncomp = 2,
                                  knn = 5,
                                  sigma = 0.73,
                                  mu = 0.5,
                                  kernel = coskern(),
                                  sample_frac = 1,
                                  lambda = 1e-4,
                                  backend = "auto",
                                  backend_control = NULL,
                                  seed = NULL,
                                  ...) {
  if (!requireNamespace("multidesign", quietly = TRUE)) {
    stop(
      "KEMA on a multidesign object requires the 'multidesign' package ",
      "(a suggested dependency). Install it, or pass a hyperdesign object.",
      call. = FALSE
    )
  }
  subject <- rlang::enquo(subject)
  y <- rlang::enquo(y)
  strata <- multidesign::hyperdesign(split(data, subject))
  kema_orig.hyperdesign(
    data = strata,
    y = !!y,
    preproc = preproc,
    ncomp = ncomp,
    knn = knn,
    sigma = sigma,
    mu = mu,
    kernel = kernel,
    sample_frac = sample_frac,
    lambda = lambda,
    backend = backend,
    backend_control = backend_control,
    seed = seed
  )
}

#' @rdname kema_orig
#' @method kema_orig hyperdesign
#' @export
kema_orig.hyperdesign <- function(data, y,
                                  preproc = center(),
                                  ncomp = 2,
                                  knn = 5,
                                  sigma = 0.73,
                                  mu = 0.5,
                                  kernel = coskern(),
                                  sample_frac = 1,
                                  lambda = 1e-4,
                                  backend = "auto",
                                  backend_control = NULL,
                                  seed = NULL,
                                  ...) {
  chk::chk_number(ncomp)
  chk::chk_true(ncomp > 0)
  chk::chk_number(knn)
  chk::chk_true(knn > 0)
  chk::chk_number(mu)
  chk::chk_true(mu >= 0)
  chk::chk_range(sample_frac, c(0, 1))
  chk::chk_number(lambda)
  if (lambda <= 0) {
    stop(
      "lambda must be strictly positive so the KEMA denominator defines an SPD metric.",
      call. = FALSE
    )
  }
  backend <- match.arg(
    backend,
    c("auto", "full_exact", "reduced_exact", "operator_exact")
  )
  control <- normalize_kema_backend_control(backend_control)
  if (!is.list(data) || length(data) == 0) {
    stop("data must be a non-empty hyperdesign list", call. = FALSE)
  }

  y <- rlang::enquo(y)
  y_char <- rlang::as_name(y)

  # Apply preprocessing per domain and preserve fitted processor state.
  pdata <- data
  n_domains <- length(data)
  is_preproc_object <- inherits(preproc, "prepper") ||
    inherits(preproc, "pre_processor")
  if (!is.list(preproc) || is_preproc_object) {
    preproc_list <- rep(list(preproc), n_domains)
  } else if (length(preproc) == 1L) {
    preproc_list <- rep(preproc, n_domains)
  } else if (length(preproc) == n_domains) {
    preproc_list <- preproc
  } else {
    stop(
      "Length of preproc list (", length(preproc),
      ") must match number of KEMA domains (", n_domains, ").",
      call. = FALSE
    )
  }

  proclist <- vector("list", n_domains)
  for (i in seq_along(data)) {
    Xi <- as.matrix(data[[i]]$x)
    pre_i <- preproc_list[[i]]
    if (is.null(pre_i)) {
      pdata[[i]]$x <- Xi
      proclist[i] <- list(NULL)
    } else if (is.function(pre_i)) {
      pdata[[i]]$x <- pre_i(Xi)
      proclist[[i]] <- pre_i
    } else {
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }
      fitted <- .preproc_fit_transform(pre_i, Xi)
      pdata[[i]]$x <- fitted$transformed
      proclist[[i]] <- fitted$preproc
    }
  }
  names(proclist) <- names(data)

  labels <- unlist(purrr::map(pdata, function(x) {
    if (!"design" %in% names(x)) {
      stop("Each domain in hyperdesign must have a 'design' component.", call. = FALSE)
    }
    if (!y_char %in% names(x$design)) {
      stop("Label column '", y_char, "' not found in domain design frame.", call. = FALSE)
    }
    x$design[[y_char]]
  }))

  if (sum(!is.na(labels)) < 2) {
    stop("kema_orig requires at least 2 labeled samples.", call. = FALSE)
  }
  if (length(unique(labels[!is.na(labels)])) < 2) {
    stop("kema_orig requires at least 2 distinct class labels.", call. = FALSE)
  }

  # Backend candidates and dense-memory guard. Every backend materialises the
  # dense kernel blocks and their eigen factors; the dense quotient backends
  # also form r_total x r_total quotient matrices, which operator_exact
  # applies matrix-free. Candidates whose estimate exceeds the guard are
  # recorded as skipped, before any kernel is allocated.
  n_per_block <- vapply(pdata, function(d) nrow(d$x), integer(1))
  landmarks <- .kema_select_landmarks(n_per_block, sample_frac, seed = seed)
  r_per_block <- if (is.null(landmarks)) {
    n_per_block
  } else {
    vapply(landmarks, length, integer(1))
  }
  auto_backend <- if (identical(backend, "auto")) {
    choose_kema_orig_backend(
      sample_frac, sum(n_per_block), sum(r_per_block), control
    )
  } else {
    backend
  }
  candidates <- kema_orig_backend_candidates(
    requested_backend = backend,
    auto_backend = auto_backend,
    sample_frac = sample_frac,
    fidelity_action = "fallback"
  )
  memory <- .kema_check_dense_memory(
    n_per_block = n_per_block,
    r_per_block = r_per_block,
    n_labeled = sum(!is.na(labels)),
    max_bytes = control$max_dense_bytes,
    candidates = candidates,
    ncomp = ncomp
  )
  attempts <- lapply(setdiff(candidates, memory$admissible), function(b) {
    list(
      backend = b,
      solved = FALSE,
      passed = FALSE,
      skipped = "max_dense_bytes",
      message = paste0(
        "skipped before allocation: estimated dense footprint ",
        format(structure(memory$estimates[[b]], class = "object_size"),
               units = "auto"),
        " exceeds backend_control$max_dense_bytes"
      )
    )
  })
  candidates <- memory$admissible
  if (!identical(control$fidelity_action, "fallback")) {
    candidates <- candidates[1L]
  }

  # concat_pre_processors() expects the feature (column) indices of each block.
  feat_sizes <- vapply(pdata, function(b) ncol(b$x), integer(1))
  feat_ends <- cumsum(feat_sizes)
  feat_starts <- c(1L, utils::head(feat_ends, -1L) + 1L)
  feature_index_list <- lapply(seq_along(feat_sizes), function(i) {
    seq.int(feat_starts[i], feat_ends[i])
  })
  # NULL or plain-function preprocessing has no fitted pre_processor; the
  # combined object then records a pass-through for that block. Prediction
  # never relies on it: the kernel encoder applies `proclist` itself.
  proc_for_concat <- lapply(seq_along(proclist), function(i) {
    p <- proclist[[i]]
    if (inherits(p, "pre_processor")) {
      return(p)
    }
    multivarious::fit_transform(multivarious::pass(), as.matrix(pdata[[i]]$x))$preproc
  })
  proc <- multivarious::concat_pre_processors(proc_for_concat, feature_index_list)

  # Build within-domain graph W (block diagonal).
  W_blocks <- lapply(pdata, function(dom) {
    g <- adjoin::graph_weights(
      dom$x,
      weight_mode = "normalized",
      neighbor_mode = "knn",
      k = knn,
      type = "normal",
      sigma = sigma
    )
    A <- adjoin::adjacency(g)
    if (!methods::is(A, "sparseMatrix")) {
      A <- Matrix::Matrix(A, sparse = TRUE)
    }
    A <- (A + Matrix::t(A)) / 2
    Matrix::diag(A) <- 0
    as(Matrix::drop0(A), "dgCMatrix")
  })
  W <- Matrix::bdiag(W_blocks)
  n_total <- nrow(W)

  # Build class graphs from labels only (unlabeled nodes contribute zero edges).
  class_graphs <- build_original_class_graphs(labels)
  Ws <- class_graphs$Ws
  Wd <- class_graphs$Wd

  # Unnormalized Laplacians (paper form).
  L <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
  Ls <- Matrix::Diagonal(x = Matrix::rowSums(Ws)) - Ws
  Ld <- Matrix::Diagonal(x = Matrix::rowSums(Wd)) - Wd

  component_basis <- kema_graph_component_basis(W + mu * Ws)

  # compute_kernels() returns dgCMatrix blocks without explicit zeros, already
  # symmetrised for full KEMA, so they are used without further copies.
  Ks <- compute_kernels(pdata, kernel, landmarks = landmarks)
  Ks_sparse <- Ks

  # Block-diagonal kernel map Z = blockdiag(K_i), applied block by block so no
  # second (block-diagonal) copy of the kernel entries is formed.
  Z <- .kema_kernel_operator(Ks_sparse)
  if (Z$nrow != n_total) {
    stop("Kernel and graph dimensions do not match.", call. = FALSE)
  }

  quotient <- build_kema_kernel_quotient(
    Ks = Ks_sparse,
    full_form = sample_frac >= 1,
    component_basis = component_basis,
    rank_tol = control$kernel_rank_tol
  )

  eig <- NULL
  chosen_backend <- NULL
  selected_fidelity <- NULL
  for (candidate in candidates) {
    solved <- tryCatch(
      switch(
        candidate,
        full_exact = solve_kema_orig_full(
          Z = Z,
          quotient = quotient,
          L = L,
          Ls = Ls,
          Ld = Ld,
          mu = mu,
          lambda = lambda,
          ncomp = ncomp,
          control = control
        ),
        reduced_exact = solve_kema_orig_rekema(
          Z = Z,
          quotient = quotient,
          L = L,
          Ls = Ls,
          Ld = Ld,
          mu = mu,
          lambda = lambda,
          ncomp = ncomp,
          control = control
        ),
        operator_exact = solve_kema_orig_operator(
          Z = Z,
          quotient = quotient,
          L = L,
          Ls = Ls,
          Ld = Ld,
          mu = mu,
          lambda = lambda,
          ncomp = ncomp,
          control = control
        )
      ),
      error = function(e) {
        structure(list(message = conditionMessage(e)), class = "kema_solver_error")
      }
    )

    if (inherits(solved, "kema_solver_error")) {
      attempts[[length(attempts) + 1L]] <- list(
        backend = candidate,
        solved = FALSE,
        passed = FALSE,
        message = solved$message
      )
      next
    }

    fidelity <- compute_kema_fidelity(
      values = solved$values,
      vectors = solved$vectors,
      A_apply = solved$A_apply,
      B_apply = solved$B_apply,
      residual_tol = control$fidelity_residual_tol,
      orth_tol = control$fidelity_orth_tol
    )
    attempts[[length(attempts) + 1L]] <- c(
      list(backend = candidate, solved = TRUE),
      fidelity
    )
    eig <- solved
    chosen_backend <- candidate
    selected_fidelity <- fidelity

    should_fallback <- identical(control$fidelity_action, "fallback") &&
      identical(backend, "auto") &&
      !isTRUE(fidelity$passed)
    if (!should_fallback) {
      break
    }
  }

  if (is.null(eig)) {
    attempt_msgs <- vapply(attempts, function(x) {
      paste0(x$backend, ": ", x$message)
    }, character(1))
    stop(
      "All KEMA backend attempts failed. ",
      paste(attempt_msgs, collapse = " | "),
      call. = FALSE
    )
  }

  vectors <- eig$vectors
  eigenvalues <- eig$values
  fidelity <- selected_fidelity
  if (!isTRUE(fidelity$passed)) {
    fidelity_msg <- paste0(
      "KEMA fidelity checks failed for backend '", chosen_backend,
      "': max_rel_residual=", signif(fidelity$max_rel_residual, 3),
      ", max_B_orth_offdiag=", signif(fidelity$max_b_orth_offdiag, 3)
    )
    stop(
      fidelity_msg,
      ". No KEMA fit was returned because its generalized-eigen solution ",
      "did not meet the configured fidelity gate.",
      call. = FALSE
    )
  }

  # In both the full (Eq. 6) and reduced (Eq. 10) forms the score of row x of
  # domain i is k(x, R_i) alpha_i, where R_i are the kernel's reference rows
  # (all rows for full KEMA, the landmark rows for REKEMA). The encoder keeps
  # exactly those ingredients so out-of-sample rows use the same map.
  scores <- Z$apply(vectors)
  if (!is.matrix(scores) && !methods::is(scores, "Matrix")) {
    scores <- matrix(scores, ncol = 1L)
  }
  scores <- as.matrix(scores)
  domain_names <- names(pdata)
  if (is.null(domain_names) || any(!nzchar(domain_names))) {
    domain_names <- paste0("domain_", seq_along(pdata))
  }
  encoder <- .kema_build_encoder(
    pdata = pdata,
    proclist = proclist,
    kernel = kernel,
    alpha = vectors,
    landmarks = landmarks,
    domain_names = domain_names
  )
  # Feature loadings exist only for the linear kernel; NA otherwise so that
  # no generic linear projection silently stands in for the kernel encoder.
  v <- .kema_encoder_loadings(encoder)

  feat_per_block <- vapply(pdata, function(b) ncol(b$x), integer(1))
  end_idx <- cumsum(feat_per_block)
  start_idx <- c(1L, head(end_idx, -1L) + 1L)
  feature_block_idx <- lapply(seq_along(feat_per_block), function(i) start_idx[i]:end_idx[i])
  names(feature_block_idx) <- paste0("block_", seq_along(feature_block_idx))

  result <- multivarious::multiblock_biprojector(
    v = v,
    s = scores,
    sdev = apply(scores, 2, stats::sd),
    preproc = proc,
    alpha = vectors,
    block_indices = feature_block_idx,
    Ks = Ks,
    sample_frac = sample_frac,
    labels = labels,
    encoder = encoder,
    landmarks = landmarks,
    loadings_status = if (isTRUE(encoder$linear_kernel)) {
      "exact_linear_kernel"
    } else {
      "unavailable_nonlinear_kernel"
    },
    classes = c("kema_fit", "kema_orig")
  )
  result$eigenvalues <- eigenvalues
  result$formulation <- paste0("kema_orig_", eig$formulation)
  result$backend <- chosen_backend
  result$backend_requested <- backend
  result$backend_candidates <- candidates
  result$fidelity <- list(
    passed = isTRUE(fidelity$passed),
    eigenvalue_zero_tol = eig$zero_tol,
    spectral_scale = eig$spectral_scale,
    n_eigenpairs_examined = eig$n_eigenpairs_examined,
    n_residual_rejected = eig$n_residual_rejected,
    eigensolver = eig$eigensolver,
    quotient_dimension = eig$eigensolver_stats$quotient_dimension,
    kernel_rank_discarded = eig$eigensolver_stats$kernel_rank_discarded,
    graph_nullity = eig$eigensolver_stats$graph_nullity,
    nullity_deflated = eig$eigensolver_stats$nullity_deflated,
    eigencore_certificate_passed = isTRUE(
      eig$eigensolver_stats$certificate$passed
    ),
    dense_quotient_formed = eig$eigensolver_stats$dense_quotient_formed,
    preconditioner = eig$eigensolver_stats$preconditioner,
    lobpcg_block_size = eig$eigensolver_stats$block_size,
    eigensolver_iterations = eig$eigensolver_stats$iterations,
    partial_attempts = eig$eigensolver_stats$partial_attempts,
    dense_memory_estimates = memory$estimates,
    max_rel_residual = fidelity$max_rel_residual,
    max_b_orth_offdiag = fidelity$max_b_orth_offdiag,
    min_b_gram_diag = fidelity$min_b_gram_diag,
    max_b_gram_diag = fidelity$max_b_gram_diag,
    residuals_rel = fidelity$residuals_rel
  )
  result$fidelity_history <- attempts
  result
}

#' @rdname kema_orig
#' @method kema_orig default
#' @export
kema_orig.default <- function(data, ...) {
  stop(
    "kema_orig() requires either a hyperdesign or multidesign object. Got: ",
    paste(class(data), collapse = ", "),
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
build_original_class_graphs <- function(labels) {
  n <- length(labels)
  Ws <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
  Wd <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))

  labeled_idx <- which(!is.na(labels))
  if (length(labeled_idx) <= 1) {
    return(list(Ws = as(Ws, "dgCMatrix"), Wd = as(Wd, "dgCMatrix")))
  }

  labs <- as.character(labels[labeled_idx])
  same <- outer(labs, labs, FUN = "==")
  diff <- outer(labs, labs, FUN = "!=")
  diag(same) <- FALSE
  diag(diff) <- FALSE

  Ws[labeled_idx, labeled_idx] <- Matrix::Matrix(same * 1, sparse = TRUE)
  Wd[labeled_idx, labeled_idx] <- Matrix::Matrix(diff * 1, sparse = TRUE)

  list(Ws = as(Ws, "dgCMatrix"), Wd = as(Wd, "dgCMatrix"))
}

#' @keywords internal
#' @noRd
select_nontrivial_eigenpairs <- function(values, vectors, ncomp,
                                         tol = sqrt(.Machine$double.eps),
                                         spectral_scale = NULL,
                                         residuals_rel = NULL,
                                         residual_tol = Inf) {
  finite_abs <- abs(values[is.finite(values)])
  observed_scale <- if (length(finite_abs)) max(finite_abs) else 0
  supplied_scale <- spectral_scale[
    is.finite(spectral_scale) & spectral_scale >= 0
  ]
  spectral_scale <- max(c(observed_scale, supplied_scale, 0))
  zero_tol <- max(100 * .Machine$double.eps, tol * spectral_scale)
  negative <- which(is.finite(values) & values < -zero_tol)
  if (length(negative)) {
    stop(
      "KEMA generalized eigenproblem returned ", length(negative),
      " materially negative eigenvalue(s); the paper formulation requires a ",
      "positive-semidefinite numerator. Smallest value: ",
      signif(min(values[negative]), 6), ".",
      call. = FALSE
    )
  }
  non_trivial_mask <- is.finite(values) & values > zero_tol
  residual_mask <- rep(TRUE, length(values))
  if (!is.null(residuals_rel)) {
    if (length(residuals_rel) != length(values)) {
      stop("residuals_rel must have one value per eigenpair.", call. = FALSE)
    }
    residual_mask <- is.finite(residuals_rel) & residuals_rel <= residual_tol
  }
  non_trivial <- which(non_trivial_mask & residual_mask)
  non_trivial <- non_trivial[order(values[non_trivial])]
  n_residual_rejected <- sum(non_trivial_mask & !residual_mask)

  if (length(non_trivial) < ncomp) {
    stop(
      "Failed to extract ", ncomp, " non-trivial KEMA eigenpairs: only ",
      length(non_trivial), " of ", length(values),
      " returned eigenpairs exceeded the numerical-zero threshold ",
      signif(zero_tol, 4), " and met the residual tolerance ",
      signif(residual_tol, 4), ". The generalized problem may have a larger ",
      "null space than anticipated; improve graph/class connectivity, add ",
      "regularization, or request fewer components.",
      call. = FALSE
    )
  }

  take <- non_trivial[seq_len(ncomp)]
  list(
    values = values[take],
    vectors = vectors[, take, drop = FALSE],
    zero_tol = zero_tol,
    spectral_scale = spectral_scale,
    n_eigenpairs_examined = length(values),
    n_residual_rejected = n_residual_rejected,
    selected_residuals_rel = if (is.null(residuals_rel)) {
      NULL
    } else {
      residuals_rel[take]
    }
  )
}

#' @keywords internal
#' @noRd
default_kema_backend_control <- function() {
  list(
    full_exact_max_n = 1200L,
    reduced_exact_max_r = 5000L,
    dense_exact_max_dim = 400L,
    dense_fallback_max_dim = 2000L,
    fidelity_residual_tol = 1e-6,
    fidelity_orth_tol = 1e-6,
    fidelity_action = "fallback",
    eigencore_tol = 1e-8,
    eigencore_maxit = 500L,
    kernel_rank_tol = sqrt(.Machine$double.eps),
    max_dense_bytes = 8 * 1024^3,
    eigencore_seed = 1L
  )
}

#' @keywords internal
#' @noRd
normalize_kema_backend_control <- function(control) {
  defaults <- default_kema_backend_control()
  if (is.null(control)) {
    return(defaults)
  }
  if (!is.list(control)) {
    stop("backend_control must be NULL or a named list.", call. = FALSE)
  }

  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) {
    warning(
      "Ignoring unknown backend_control fields: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }

  out <- defaults
  shared <- intersect(names(control), names(defaults))
  for (nm in shared) {
    out[[nm]] <- control[[nm]]
  }

  out$full_exact_max_n <- as.integer(out$full_exact_max_n)
  out$reduced_exact_max_r <- as.integer(out$reduced_exact_max_r)
  out$dense_exact_max_dim <- as.integer(out$dense_exact_max_dim)
  out$dense_fallback_max_dim <- as.integer(out$dense_fallback_max_dim)
  out$eigencore_maxit <- as.integer(out$eigencore_maxit)
  if (length(out$fidelity_action) != 1L ||
      !is.character(out$fidelity_action) ||
      !out$fidelity_action %in% c("fallback", "error")) {
    stop(
      "backend_control$fidelity_action must be one of 'fallback' or 'error'.",
      call. = FALSE
    )
  }
  chk::chk_true(out$full_exact_max_n > 2)
  chk::chk_true(out$reduced_exact_max_r > 2)
  chk::chk_true(out$dense_exact_max_dim > 2)
  chk::chk_true(out$dense_fallback_max_dim > 2)
  chk::chk_number(out$fidelity_residual_tol)
  chk::chk_true(out$fidelity_residual_tol > 0)
  chk::chk_number(out$fidelity_orth_tol)
  chk::chk_true(out$fidelity_orth_tol > 0)
  chk::chk_number(out$eigencore_tol)
  chk::chk_true(out$eigencore_tol > 0)
  chk::chk_true(out$eigencore_maxit > 0)
  chk::chk_number(out$kernel_rank_tol)
  chk::chk_true(out$kernel_rank_tol > 0)
  chk::chk_true(out$kernel_rank_tol < 1)
  chk::chk_number(out$max_dense_bytes)
  chk::chk_true(out$max_dense_bytes > 0)
  chk::chk_number(out$eigencore_seed)
  out$eigencore_seed <- as.integer(out$eigencore_seed)

  out
}

#' @keywords internal
#' @noRd
choose_kema_orig_backend <- function(sample_frac, n_total, r_total, control) {
  if (sample_frac >= 1 && n_total <= control$full_exact_max_n) {
    return("full_exact")
  }
  if (r_total <= control$reduced_exact_max_r) {
    return("reduced_exact")
  }
  "operator_exact"
}

#' @keywords internal
#' @noRd
kema_orig_backend_candidates <- function(requested_backend, auto_backend, sample_frac, fidelity_action) {
  if (sample_frac < 1 && identical(requested_backend, "full_exact")) {
    stop(
      "backend='full_exact' requires sample_frac = 1 (square block kernels).",
      call. = FALSE
    )
  }

  if (!identical(requested_backend, "auto")) {
    return(requested_backend)
  }

  ordered <- unique(c(auto_backend, "full_exact", "reduced_exact", "operator_exact"))
  if (sample_frac < 1) {
    ordered <- ordered[ordered != "full_exact"]
  }
  if (!identical(fidelity_action, "fallback")) {
    ordered <- ordered[1L]
  }
  ordered
}

#' @keywords internal
#' @noRd
compute_kema_relative_residuals <- function(values, vectors, A_apply, B_apply) {
  V <- if (is.matrix(vectors) || methods::is(vectors, "Matrix")) {
    vectors
  } else {
    matrix(vectors, ncol = 1L)
  }
  Av <- A_apply(V)
  Bv <- B_apply(V)
  if (!is.matrix(Av) && !methods::is(Av, "Matrix")) {
    Av <- matrix(Av, ncol = 1L)
  }
  if (!is.matrix(Bv) && !methods::is(Bv, "Matrix")) {
    Bv <- matrix(Bv, ncol = 1L)
  }

  resid <- Av - sweep(Bv, 2, values, "*")
  resid_num <- sqrt(colSums(resid^2))
  resid_den <- sqrt(colSums(Av^2)) + abs(values) * sqrt(colSums(Bv^2))
  resid_den[resid_den < .Machine$double.eps] <- 1
  as.numeric(resid_num / resid_den)
}

#' @keywords internal
#' @noRd
compute_kema_fidelity <- function(values, vectors, A_apply, B_apply, residual_tol, orth_tol) {
  V <- if (is.matrix(vectors) || methods::is(vectors, "Matrix")) {
    vectors
  } else {
    matrix(vectors, ncol = 1L)
  }
  Bv <- B_apply(V)
  if (!is.matrix(Bv) && !methods::is(Bv, "Matrix")) {
    Bv <- matrix(Bv, ncol = 1L)
  }
  residuals_rel <- compute_kema_relative_residuals(
    values,
    V,
    A_apply,
    B_apply
  )

  gram <- as.matrix(Matrix::crossprod(V, Bv))
  gram_diag <- diag(gram)
  valid_diag <- is.finite(gram_diag) & (gram_diag > .Machine$double.eps)

  if (!all(valid_diag)) {
    max_offdiag <- Inf
    min_diag <- suppressWarnings(min(gram_diag, na.rm = TRUE))
    max_diag <- suppressWarnings(max(gram_diag, na.rm = TRUE))
    passed <- FALSE
  } else {
    d <- sqrt(gram_diag)
    gram_norm <- sweep(sweep(gram, 1L, d, "/"), 2L, d, "/")
    offdiag <- gram_norm
    diag(offdiag) <- 0
    max_offdiag <- max(abs(offdiag))
    min_diag <- min(diag(gram_norm))
    max_diag <- max(diag(gram_norm))
    passed <- all(is.finite(residuals_rel)) &&
      max(residuals_rel) <= residual_tol &&
      is.finite(max_offdiag) &&
      max_offdiag <= orth_tol
  }

  list(
    passed = isTRUE(passed),
    max_rel_residual = max(residuals_rel),
    max_b_orth_offdiag = max_offdiag,
    min_b_gram_diag = min_diag,
    max_b_gram_diag = max_diag,
    residuals_rel = residuals_rel
  )
}

#' Exact graph-Laplacian null basis from connected components
#'
#' @keywords internal
#' @noRd
kema_graph_component_basis <- function(adjacency) {
  adjacency <- Matrix::Matrix(adjacency, sparse = TRUE)
  adjacency <- (adjacency + Matrix::t(adjacency)) / 2
  Matrix::diag(adjacency) <- 0
  adjacency <- Matrix::drop0(adjacency)
  n <- nrow(adjacency)
  if (n < 1L) {
    stop("KEMA graph must contain at least one sample.", call. = FALSE)
  }

  graph <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  membership <- igraph::components(graph)$membership
  component_sizes <- tabulate(membership, nbins = max(membership))
  Matrix::sparseMatrix(
    i = seq_len(n),
    j = membership,
    x = 1 / sqrt(component_sizes[membership]),
    dims = c(n, length(component_sizes))
  )
}

#' Compact coefficient basis for the numerical image of the kernel blocks
#'
#' @keywords internal
#' @noRd
kema_kernel_block_basis <- function(K, full_form, rank_tol) {
  K <- as.matrix(K)
  storage.mode(K) <- "double"

  if (isTRUE(full_form)) {
    if (nrow(K) != ncol(K)) {
      stop("Full KEMA requires square kernel blocks.", call. = FALSE)
    }
    K <- (K + t(K)) / 2
    decomp <- eigen(K, symmetric = TRUE)
    scale <- max(abs(decomp$values), 0)
    threshold <- max(100 * .Machine$double.eps, rank_tol * scale)
    materially_negative <- decomp$values < -threshold
    if (any(materially_negative)) {
      stop(
        "KEMA kernel block is not positive semidefinite; smallest eigenvalue is ",
        signif(min(decomp$values), 6), ".",
        call. = FALSE
      )
    }
    keep <- which(decomp$values > threshold)
    if (!length(keep)) {
      stop("KEMA kernel block has zero numerical rank.", call. = FALSE)
    }
    U <- decomp$vectors[, keep, drop = FALSE]
    return(list(
      U = U,
      d = decomp$values[keep],
      V = U,
      rank = length(keep),
      original_rank = ncol(K),
      threshold = threshold
    ))
  }

  decomp <- svd(K, nu = min(dim(K)), nv = min(dim(K)))
  scale <- max(decomp$d, 0)
  threshold <- max(100 * .Machine$double.eps, rank_tol * scale)
  keep <- which(decomp$d > threshold)
  if (!length(keep)) {
    stop("REKEMA kernel block has zero numerical rank.", call. = FALSE)
  }
  list(
    U = decomp$u[, keep, drop = FALSE],
    d = decomp$d[keep],
    V = decomp$v[, keep, drop = FALSE],
    rank = length(keep),
    original_rank = ncol(K),
    threshold = threshold
  )
}

#' Build the KEMA quotient by kernel null directions and graph components
#'
#' @keywords internal
#' @noRd
build_kema_kernel_quotient <- function(Ks, full_form, component_basis,
                                       rank_tol = sqrt(.Machine$double.eps)) {
  pieces <- lapply(
    Ks,
    kema_kernel_block_basis,
    full_form = full_form,
    rank_tol = rank_tol
  )
  # The quotient basis is block diagonal, H = blockdiag(U_i diag(d_i)) and
  # V = blockdiag(V_i). Only the dense per-domain factors are stored (for the
  # full form V_i is the same object as U_i, so it costs no extra memory);
  # products with H and V are applied block by block, and a sparse
  # block-diagonal H is formed only by the dense quotient backends.
  U_blocks <- lapply(pieces, `[[`, "U")
  V_blocks <- lapply(pieces, `[[`, "V")
  d_blocks <- lapply(pieces, `[[`, "d")
  singular_values <- unlist(d_blocks, use.names = FALSE)
  row_sizes <- vapply(U_blocks, nrow, integer(1))

  component_basis <- Matrix::Matrix(component_basis, sparse = TRUE)
  if (nrow(component_basis) != sum(row_sizes)) {
    stop("Kernel and graph-component dimensions do not match.", call. = FALSE)
  }
  row_ends <- cumsum(row_sizes)
  row_starts <- c(1L, utils::head(row_ends, -1L) + 1L)
  UtC <- do.call(rbind, lapply(seq_along(U_blocks), function(i) {
    rows <- seq.int(row_starts[i], length.out = row_sizes[i])
    as.matrix(Matrix::crossprod(
      U_blocks[[i]],
      component_basis[rows, , drop = FALSE]
    ))
  }))
  residual_gram <- as.matrix(
    Matrix::crossprod(component_basis) - Matrix::crossprod(UtC)
  )
  residual_gram <- (residual_gram + t(residual_gram)) / 2
  intersection <- eigen(residual_gram, symmetric = TRUE)
  intersection_scale <- max(abs(intersection$values), 1)
  intersection_tol <- max(
    100 * .Machine$double.eps,
    rank_tol * intersection_scale
  )
  keep_intersection <- which(intersection$values <= intersection_tol)
  constraints <- if (length(keep_intersection)) {
    component_coefficients <- intersection$vectors[
      , keep_intersection, drop = FALSE
    ]
    sweep(
      UtC %*% component_coefficients,
      1L,
      singular_values,
      "/"
    )
  } else {
    matrix(numeric(0), nrow = length(singular_values), ncol = 0L)
  }

  structure(list(
    U_blocks = U_blocks,
    d_blocks = d_blocks,
    V_blocks = V_blocks,
    row_sizes = row_sizes,
    col_sizes = vapply(V_blocks, nrow, integer(1)),
    constraints = constraints,
    rank = length(singular_values),
    original_rank = sum(vapply(pieces, `[[`, integer(1), "original_rank")),
    discarded_rank = sum(vapply(
      pieces,
      function(x) x$original_rank - x$rank,
      integer(1)
    )),
    graph_nullity = ncol(component_basis),
    nullity_deflated = ncol(constraints),
    block_ranks = vapply(pieces, `[[`, integer(1), "rank"),
    block_thresholds = vapply(pieces, `[[`, numeric(1), "threshold")
  ), class = "kema_kernel_quotient")
}

#' Block offsets of a KEMA kernel quotient
#'
#' @keywords internal
#' @noRd
.kema_quotient_offsets <- function(sizes) {
  sizes <- as.integer(sizes)
  ends <- cumsum(sizes)
  starts <- c(1L, utils::head(ends, -1L) + 1L)
  lapply(seq_along(sizes), function(i) seq.int(starts[i], length.out = sizes[i]))
}

#' Sparse block-diagonal quotient basis `H = blockdiag(U_i diag(d_i))`
#'
#' Used by the dense quotient backends and by tests; `operator_exact` never
#' forms it.
#'
#' @keywords internal
#' @noRd
.kema_quotient_H <- function(quotient) {
  Matrix::bdiag(lapply(seq_along(quotient$U_blocks), function(i) {
    Matrix::Matrix(
      sweep(quotient$U_blocks[[i]], 2L, quotient$d_blocks[[i]], "*"),
      sparse = TRUE
    )
  }))
}

#' Sparse block-diagonal coefficient map `V = blockdiag(V_i)`
#'
#' @keywords internal
#' @noRd
.kema_quotient_V <- function(quotient) {
  Matrix::bdiag(lapply(quotient$V_blocks, Matrix::Matrix, sparse = TRUE))
}

#' Blockwise products with the quotient factors
#'
#' `.kema_quotient_H_apply(q, X)` returns `H X` (`n x k`),
#' `.kema_quotient_Ht_apply(q, Y)` returns `t(H) Y` (`rank x k`) and
#' `.kema_quotient_map(q, X)` returns `V X` (original coefficients). All are
#' dense BLAS products over the per-domain factors.
#'
#' @keywords internal
#' @noRd
.kema_quotient_H_apply <- function(quotient, X) {
  X <- as.matrix(X)
  cols <- .kema_quotient_offsets(quotient$block_ranks)
  do.call(rbind, lapply(seq_along(quotient$U_blocks), function(i) {
    quotient$U_blocks[[i]] %*%
      (quotient$d_blocks[[i]] * X[cols[[i]], , drop = FALSE])
  }))
}

#' @keywords internal
#' @noRd
.kema_quotient_Ht_apply <- function(quotient, Y) {
  Y <- as.matrix(Y)
  rows <- .kema_quotient_offsets(quotient$row_sizes)
  do.call(rbind, lapply(seq_along(quotient$U_blocks), function(i) {
    quotient$d_blocks[[i]] *
      crossprod(quotient$U_blocks[[i]], Y[rows[[i]], , drop = FALSE])
  }))
}

#' @keywords internal
#' @noRd
.kema_quotient_map <- function(quotient, X) {
  X <- as.matrix(X)
  cols <- .kema_quotient_offsets(quotient$block_ranks)
  do.call(rbind, lapply(seq_along(quotient$V_blocks), function(i) {
    quotient$V_blocks[[i]] %*% X[cols[[i]], , drop = FALSE]
  }))
}

#' Block-diagonal kernel map applied block by block
#'
#' Represents `Z = blockdiag(K_1, ..., K_m)` through the list of kernel blocks:
#' `apply(X)` returns `Z X` and `crossprod(Y)` returns `t(Z) Y`, both as dense
#' matrices, without forming a block-diagonal copy or a transpose of `Z`.
#'
#' @keywords internal
#' @noRd
.kema_kernel_operator <- function(Ks) {
  row_sizes <- vapply(Ks, nrow, integer(1))
  col_sizes <- vapply(Ks, ncol, integer(1))
  rows <- .kema_quotient_offsets(row_sizes)
  cols <- .kema_quotient_offsets(col_sizes)
  list(
    nrow = sum(row_sizes),
    ncol = sum(col_sizes),
    apply = function(X) {
      X <- as.matrix(X)
      do.call(rbind, lapply(seq_along(Ks), function(i) {
        as.matrix(Ks[[i]] %*% X[cols[[i]], , drop = FALSE])
      }))
    },
    crossprod = function(Y) {
      Y <- as.matrix(Y)
      do.call(rbind, lapply(seq_along(Ks), function(i) {
        as.matrix(Matrix::crossprod(Ks[[i]], Y[rows[[i]], , drop = FALSE]))
      }))
    }
  )
}

#' LOBPCG block width for a KEMA quotient solve
#'
#' The partial solver iterates on a block wider than `ncomp` (guard vectors):
#' convergence of the `ncomp`-th eigenpair is then governed by its gap to the
#' first eigenvalue outside the block rather than to its nearest neighbour,
#' and residual-rejected pairs can be skipped as in the dense solve. The block
#' is capped so the LOBPCG search space (three blocks) fits in the deflated
#' quotient.
#'
#' @keywords internal
#' @noRd
.kema_lobpcg_block_size <- function(ncomp, available) {
  ncomp <- as.integer(ncomp)
  wanted <- max(2L * ncomp, ncomp + 4L)
  cap <- max(ncomp, as.integer(floor(available / 3)))
  as.integer(min(wanted, cap))
}

#' Dense KEMA quotient pencil `(H' M H, H' Ld H + lambda I)`
#'
#' Used by `full_exact` and `reduced_exact` (and by the small-quotient dense
#' fallback of `operator_exact`).
#'
#' @keywords internal
#' @noRd
.kema_dense_quotient_pencil <- function(quotient, M, Ld, lambda) {
  H <- .kema_quotient_H(quotient)
  q <- ncol(H)
  A <- as.matrix(Matrix::crossprod(H, M %*% H))
  B <- as.matrix(
    Matrix::crossprod(H, Ld %*% H) +
      lambda * Matrix::Diagonal(q)
  )
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"
  list(A = A, B = B)
}

#' Shifted-Cholesky preconditioner for the smallest pencil eigenpairs
#'
#' Factors `A_i + shift * B_i` for each diagonal block (a single block for the
#' dense quotient, one block per domain for the operator path) and returns a
#' function applying `blockdiag(A_i + shift B_i)^{-1}`. The shift is a small
#' multiple of `trace(A) / trace(B)`, the B-weighted mean generalized
#' eigenvalue; the preconditioner only changes the LOBPCG search directions,
#' never the pencil or its certification. Returns `NULL` when no block
#' factorisation succeeds, in which case LOBPCG runs unpreconditioned.
#'
#' @param trace_A,trace_B Traces of the full pencil matrices.
#' @param block_matrix `function(i, shift)` returning block `i` of
#'   `A + shift * B` (only its upper triangle is read by `chol()`).
#' @param index List of row indices of each block.
#' @keywords internal
#' @noRd
.kema_shifted_cholesky_preconditioner <- function(trace_A, trace_B,
                                                  block_matrix, index) {
  base <- trace_A / trace_B
  if (!is.finite(base) || base <= 0) {
    base <- 1
  }
  for (shift in base * c(1e-3, 1)) {
    factors <- tryCatch(
      lapply(seq_along(index), function(i) chol(block_matrix(i, shift))),
      error = function(e) NULL
    )
    if (!is.null(factors)) {
      return(list(
        apply = .kema_block_cholesky_apply(factors, index),
        shift = shift
      ))
    }
  }
  NULL
}

#' Closures handed to eigencore
#'
#' eigencore hashes the solver method descriptor, including a preconditioner
#' closure and everything its environment references, when it plans and
#' validates a solve. These factories therefore build closures whose
#' environments hold only the data the closure needs (never the quotient,
#' the pencil or the graph matrices), which keeps that transient
#' serialisation proportional to the preconditioner factors.
#'
#' @keywords internal
#' @noRd
.kema_block_cholesky_apply <- function(factors, index) {
  force(factors)
  force(index)
  function(X) {
    X <- as.matrix(X)
    for (i in seq_along(factors)) {
      idx <- index[[i]]
      X[idx, ] <- backsolve(
        factors[[i]],
        backsolve(factors[[i]], X[idx, , drop = FALSE], transpose = TRUE)
      )
    }
    X
  }
}

#' @keywords internal
#' @noRd
.kema_identity_plus_low_rank_apply <- function(Q, identity_scale, low_rank_scale) {
  force(Q)
  force(identity_scale)
  force(low_rank_scale)
  function(X) {
    X <- as.matrix(X)
    identity_scale * X + Q %*% (low_rank_scale * crossprod(Q, X))
  }
}

#' @keywords internal
#' @noRd
.kema_sandwich_apply <- function(outer, inner) {
  force(outer)
  force(inner)
  function(X) outer(inner(outer(X)))
}

#' Column chunks of a KEMA domain factor
#'
#' @keywords internal
#' @noRd
.kema_column_chunks <- function(r, chunk = 256L) {
  starts <- seq.int(1L, r, by = chunk)
  lapply(starts, function(s) seq.int(s, min(r, s + chunk - 1L)))
}

#' Diagonal of `d * t(U) S U * d` for one domain block
#'
#' @keywords internal
#' @noRd
.kema_quotient_block_diag <- function(U, d, S, chunk = 256L) {
  out <- numeric(ncol(U))
  for (cols in .kema_column_chunks(ncol(U), chunk)) {
    Uc <- U[, cols, drop = FALSE]
    out[cols] <- colSums(Uc * as.matrix(S %*% Uc))
  }
  d^2 * out
}

#' Congruence `d * t(U) S U * d + diag_add * I` of one domain block
#'
#' Computed in column chunks so the transient `S U` product holds at most
#' `n_i x chunk` entries; the result is the only `r_i x r_i` allocation. Only
#' its upper triangle is consumed (by `chol()`).
#'
#' @keywords internal
#' @noRd
.kema_quotient_block_gram <- function(U, d, S, diag_add = 0, chunk = 256L) {
  r <- ncol(U)
  out <- matrix(0, r, r)
  for (cols in .kema_column_chunks(r, chunk)) {
    block <- crossprod(U, as.matrix(S %*% U[, cols, drop = FALSE]))
    out[, cols] <- d * sweep(block, 2L, d[cols], "*")
  }
  diag(out) <- diag(out) + diag_add
  out
}

#' Square root of the KEMA quotient metric without forming it
#'
#' The quotient metric is `B = lambda I + t(H) Ld H`. `Ld` is a graph
#' Laplacian supported on the labelled rows, so `t(H) Ld H = G t(G)` with
#' `G = t(H_L) S` (`H_L` the labelled rows of `H`, `Ld_LL = S t(S)`), which has
#' at most `n_labelled` columns. With the thin SVD `G = Q diag(s) W'`,
#' `B^{1/2} = sqrt(lambda) I + Q diag(sqrt(lambda + s^2) - sqrt(lambda)) Q'`
#' and `B^{-1/2}` is analogous; both cost `O(rank * n_labelled)` per column.
#'
#' @keywords internal
#' @noRd
.kema_quotient_metric_sqrt <- function(quotient, Ld, lambda) {
  q <- quotient$rank
  labelled <- which(Matrix::diag(Ld) > 0)
  Q <- matrix(0, q, 0L)
  s2 <- numeric(0)
  if (length(labelled)) {
    Ld_LL <- as.matrix(Ld[labelled, labelled, drop = FALSE])
    Ld_LL <- (Ld_LL + t(Ld_LL)) / 2
    eig_L <- eigen(Ld_LL, symmetric = TRUE)
    keep <- which(eig_L$values > 0)
    if (length(keep)) {
      S <- sweep(
        eig_L$vectors[, keep, drop = FALSE], 2L,
        sqrt(eig_L$values[keep]), "*"
      )
      rows <- .kema_quotient_offsets(quotient$row_sizes)
      cols <- .kema_quotient_offsets(quotient$block_ranks)
      H_L <- matrix(0, length(labelled), q)
      for (i in seq_along(rows)) {
        local <- match(labelled, rows[[i]])
        hit <- which(!is.na(local))
        if (length(hit)) {
          H_L[hit, cols[[i]]] <- sweep(
            quotient$U_blocks[[i]][local[hit], , drop = FALSE], 2L,
            quotient$d_blocks[[i]], "*"
          )
        }
      }
      G <- crossprod(H_L, S)
      sv <- svd(G, nu = ncol(G), nv = 0L)
      positive <- which(sv$d > 0)
      Q <- sv$u[, positive, drop = FALSE]
      s2 <- sv$d[positive]^2
    }
  }
  root_lambda <- sqrt(lambda)
  list(
    apply = .kema_identity_plus_low_rank_apply(
      Q, root_lambda, sqrt(lambda + s2) - root_lambda
    ),
    apply_inverse = .kema_identity_plus_low_rank_apply(
      Q, 1 / root_lambda, 1 / sqrt(lambda + s2) - 1 / root_lambda
    ),
    rank = ncol(Q)
  )
}

#' Matrix-free KEMA quotient problem for `operator_exact`
#'
#' The quotient matrices `A = H' M H` and `B = H' Ld H + lambda I` are applied
#' as operators through the per-domain kernel factors (`H X`, a sparse product
#' with `M` or `Ld`, then `t(H) Y`), so no `rank x rank` quotient matrix is
#' formed. The pencil is solved in the equivalent standard form
#' `C y = theta y`, `C = B^{-1/2} A B^{-1/2}`, `x = B^{-1/2} y`, using the
#' low-rank metric square root above: eigencore can then verify target
#' completeness (that the returned pairs are the smallest ones) without
#' materialising `B`. Constraints in the `B` inner product map to Euclidean
#' constraints `B^{1/2} c`. The preconditioner is the block-Jacobi
#' (per-domain) shifted Cholesky factor of the pencil, transformed the same
#' way; the omitted off-diagonal blocks come only from the label graphs and
#' have rank at most twice the number of labelled rows.
#'
#' @keywords internal
#' @noRd
.kema_operator_quotient_pencil <- function(quotient, M, Ld, lambda) {
  q <- quotient$rank
  A_apply <- function(X) {
    X <- as.matrix(X)
    HX <- .kema_quotient_H_apply(quotient, X)
    .kema_quotient_Ht_apply(quotient, as.matrix(M %*% HX))
  }
  B_apply <- function(X) {
    X <- as.matrix(X)
    HX <- .kema_quotient_H_apply(quotient, X)
    .kema_quotient_Ht_apply(quotient, as.matrix(Ld %*% HX)) + lambda * X
  }
  metric <- .kema_quotient_metric_sqrt(quotient, Ld, lambda)
  C_apply <- function(Y) {
    metric$apply_inverse(A_apply(metric$apply_inverse(Y)))
  }

  rows <- .kema_quotient_offsets(quotient$row_sizes)
  cols <- .kema_quotient_offsets(quotient$block_ranks)
  M_blocks <- lapply(rows, function(idx) M[idx, idx, drop = FALSE])
  Ld_blocks <- lapply(rows, function(idx) Ld[idx, idx, drop = FALSE])
  trace_A <- sum(vapply(seq_along(rows), function(i) {
    sum(.kema_quotient_block_diag(
      quotient$U_blocks[[i]], quotient$d_blocks[[i]], M_blocks[[i]]
    ))
  }, numeric(1)))
  trace_B <- lambda * q + sum(vapply(seq_along(rows), function(i) {
    sum(.kema_quotient_block_diag(
      quotient$U_blocks[[i]], quotient$d_blocks[[i]], Ld_blocks[[i]]
    ))
  }, numeric(1)))
  pencil_preconditioner <- .kema_shifted_cholesky_preconditioner(
    trace_A,
    trace_B,
    function(i, shift) {
      .kema_quotient_block_gram(
        quotient$U_blocks[[i]],
        quotient$d_blocks[[i]],
        M_blocks[[i]] + shift * Ld_blocks[[i]],
        diag_add = shift * lambda
      )
    },
    cols
  )
  rm(M_blocks, Ld_blocks)
  preconditioner <- if (is.null(pencil_preconditioner)) {
    NULL
  } else {
    list(
      apply = .kema_sandwich_apply(metric$apply, pencil_preconditioner$apply),
      shift = pencil_preconditioner$shift
    )
  }

  list(
    C = eigencore::linear_operator(
      dim = c(q, q),
      apply = function(X, alpha = 1, beta = 0, Y = NULL) {
        Z <- alpha * C_apply(X)
        if (is.null(Y) || beta == 0) Z else Z + beta * Y
      },
      structure = eigencore::hermitian(),
      name = "KEMA quotient pencil in metric-whitened standard form"
    ),
    to_original = metric$apply_inverse,
    constraints_to_standard = metric$apply,
    A_apply = A_apply,
    B_apply = B_apply,
    metric_rank = metric$rank,
    preconditioner = preconditioner
  )
}

#' Form and solve the nullspace-safe KEMA quotient pencil
#'
#' Solver routes. Every returned result must pass the eigencore certificate;
#' its pairs are then screened by the original-pencil relative residual
#' (`fidelity_residual_tol`) and the caller applies the full fidelity gate.
#' * dense quotient (`operator = FALSE`): the full dense generalized solve when
#'   the quotient dimension is at most `dense_exact_max_dim` (and
#'   `force_partial` is `FALSE`), otherwise shifted-Cholesky preconditioned
#'   LOBPCG with guard vectors.
#' * operator (`operator = TRUE`): block-Jacobi preconditioned LOBPCG on the
#'   matrix-free quotient pencil in metric-whitened standard form.
#' A certified partial solve whose pairs miss the original-pencil residual
#' target is repeated with a tolerance tightened by 100 (twice at most, not
#' below 1e-14): eigencore's normwise backward error does not bound the
#' per-pair relative residual of small eigenvalues. When no partial solve
#' succeeds, the dense generalized solve is used if the dimension is at most
#' `dense_fallback_max_dim` (`dense_exact_max_dim` for `operator_exact`);
#' otherwise the last failure is reported.
#'
#' @keywords internal
#' @noRd
solve_kema_quotient <- function(quotient, M, Ld, lambda, ncomp, control,
                                force_partial = FALSE,
                                original_A_apply = NULL,
                                original_B_apply = NULL,
                                operator = FALSE) {
  q <- quotient$rank
  n_constraints <- ncol(quotient$constraints)
  available <- q - n_constraints
  if (available < ncomp) {
    stop(
      "KEMA quotient has only ", available,
      " positive-dimensional direction(s) after nullspace deflation; need ",
      ncomp, ".",
      call. = FALSE
    )
  }
  constraints <- if (n_constraints) quotient$constraints else NULL
  block_size <- .kema_lobpcg_block_size(ncomp, available)
  use_original_residual <- is.function(original_A_apply) &&
    is.function(original_B_apply)

  dense_pencil <- NULL
  get_dense_pencil <- function() {
    if (is.null(dense_pencil)) {
      dense_pencil <<- .kema_dense_quotient_pencil(quotient, M, Ld, lambda)
    }
    dense_pencil
  }
  dense_apply <- function(which) {
    function(x) get_dense_pencil()[[which]] %*% x
  }
  solve_dense_full <- function() {
    pencil <- get_dense_pencil()
    eigencore::eig_full(
      pencil$A,
      B = pencil$B,
      structure = eigencore::hermitian()
    )
  }
  certification_error <- function(fit, cert) {
    paste0(
      "eigencore failed to certify the KEMA quotient solve (method: ",
      fit$method, ", max backward error: ",
      signif(cert$max_backward_error, 4),
      if (!is.null(cert$target_completeness)) {
        paste0(", target completeness: ", cert$target_completeness)
      },
      ")."
    )
  }
  # Certified pairs -> original coefficients, residuals and selection.
  finish <- function(fit, to_original, A_apply, B_apply) {
    values <- as.numeric(Re(eigencore::values(fit)))
    vectors <- as.matrix(Re(eigencore::vectors(fit)))
    if (!is.null(to_original)) {
      vectors <- to_original(vectors)
    }
    mapped_vectors <- .kema_quotient_map(quotient, vectors)
    residuals_rel <- if (use_original_residual) {
      compute_kema_relative_residuals(
        values, mapped_vectors,
        A_apply = original_A_apply, B_apply = original_B_apply
      )
    } else {
      compute_kema_relative_residuals(
        values, vectors, A_apply = A_apply, B_apply = B_apply
      )
    }
    select_nontrivial_eigenpairs(
      values,
      mapped_vectors,
      ncomp,
      spectral_scale = max(abs(values), 0),
      residuals_rel = residuals_rel,
      residual_tol = control$fidelity_residual_tol
    )
  }

  # Partial (LOBPCG) problem for this route.
  preconditioner <- NULL
  partial <- NULL
  if (isTRUE(operator)) {
    operator_pencil <- .kema_operator_quotient_pencil(quotient, M, Ld, lambda)
    preconditioner <- operator_pencil$preconditioner
    partial <- list(
      route = "operator_lobpcg",
      A = operator_pencil$C,
      B = NULL,
      constraints = if (is.null(constraints)) {
        NULL
      } else {
        operator_pencil$constraints_to_standard(constraints)
      },
      to_original = operator_pencil$to_original,
      A_apply = operator_pencil$A_apply,
      B_apply = operator_pencil$B_apply
    )
    dense_fallback_max_dim <- control$dense_exact_max_dim
  } else if (isTRUE(force_partial) || q > control$dense_exact_max_dim) {
    pencil <- get_dense_pencil()
    preconditioner <- .kema_shifted_cholesky_preconditioner(
      sum(diag(pencil$A)),
      sum(diag(pencil$B)),
      function(i, shift) pencil$A + shift * pencil$B,
      list(seq_len(q))
    )
    partial <- list(
      route = "dense_lobpcg",
      A = pencil$A,
      B = pencil$B,
      constraints = constraints,
      to_original = NULL,
      A_apply = dense_apply("A"),
      B_apply = dense_apply("B")
    )
    rm(pencil)
    dense_fallback_max_dim <- control$dense_fallback_max_dim
  }

  partial_attempts <- list()
  selected <- NULL
  fit <- NULL
  cert <- NULL
  route <- NULL
  last_error <- NULL
  if (!is.null(partial)) {
    tolerances <- control$eigencore_tol * c(1, 1e-2, 1e-4)
    tolerances <- c(
      tolerances[1L],
      tolerances[-1L][tolerances[-1L] >= 1e-14]
    )
    for (tol in tolerances) {
      fit_try <- eigencore::eig_partial(
        partial$A,
        B = partial$B,
        k = block_size,
        target = eigencore::smallest(),
        method = eigencore::lobpcg(
          maxit = control$eigencore_maxit,
          constraints = partial$constraints,
          preconditioner = preconditioner$apply
        ),
        tol = tol,
        seed = control$eigencore_seed,
        allow_dense_fallback = "never"
      )
      cert_try <- eigencore::certificate(fit_try)
      attempt <- list(
        tol = tol,
        method = fit_try$method,
        certified = isTRUE(cert_try$passed),
        max_backward_error = cert_try$max_backward_error,
        target_completeness = cert_try$target_completeness
      )
      if (!isTRUE(cert_try$passed)) {
        last_error <- certification_error(fit_try, cert_try)
        attempt$message <- last_error
        partial_attempts[[length(partial_attempts) + 1L]] <- attempt
        break
      }
      outcome <- tryCatch(
        finish(fit_try, partial$to_original, partial$A_apply, partial$B_apply),
        error = function(e) e
      )
      if (inherits(outcome, "error")) {
        last_error <- conditionMessage(outcome)
        attempt$message <- last_error
        partial_attempts[[length(partial_attempts) + 1L]] <- attempt
        next
      }
      partial_attempts[[length(partial_attempts) + 1L]] <- attempt
      selected <- outcome
      fit <- fit_try
      cert <- cert_try
      route <- partial$route
      break
    }
  }

  if (is.null(selected)) {
    dense_allowed <- is.null(partial) || q <= dense_fallback_max_dim
    if (!dense_allowed) {
      stop(
        paste(
          unique(vapply(partial_attempts, `[[`, character(1), "message")),
          collapse = " | "
        ),
        call. = FALSE
      )
    }
    preconditioner <- NULL
    fit <- solve_dense_full()
    cert <- eigencore::certificate(fit)
    if (!isTRUE(cert$passed)) {
      stop(
        certification_error(fit, cert),
        if (length(partial_attempts)) {
          paste0(
            " Partial solve before the dense fallback: ",
            paste(
              unique(vapply(partial_attempts, `[[`, character(1), "message")),
              collapse = " | "
            )
          )
        },
        call. = FALSE
      )
    }
    selected <- finish(fit, NULL, dense_apply("A"), dense_apply("B"))
    route <- "dense_full"
  }

  selected$eigensolver <- switch(
    route,
    dense_full = "eigencore_dense_full",
    dense_lobpcg = "eigencore_lobpcg_deflated",
    operator_lobpcg = "eigencore_lobpcg_operator"
  )
  iterations <- tryCatch(
    eigencore::diagnostics(fit)$iterations,
    error = function(e) NULL
  )
  selected$eigensolver_stats <- list(
    method = fit$method,
    certificate = cert,
    quotient_dimension = q,
    original_dimension = quotient$original_rank,
    kernel_rank_discarded = quotient$discarded_rank,
    graph_nullity = quotient$graph_nullity,
    nullity_deflated = quotient$nullity_deflated,
    full_spectrum = identical(route, "dense_full"),
    dense_quotient_formed = !is.null(dense_pencil),
    block_size = if (identical(route, "dense_full")) NA_integer_ else block_size,
    iterations = if (is.numeric(iterations)) iterations else NA_integer_,
    preconditioner = if (is.null(preconditioner)) {
      "none"
    } else if (isTRUE(operator)) {
      "block_jacobi_shifted_cholesky"
    } else {
      "shifted_cholesky"
    },
    preconditioner_shift = if (is.null(preconditioner)) {
      NA_real_
    } else {
      preconditioner$shift
    },
    partial_attempts = partial_attempts
  )
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_full <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  n <- Z$nrow
  if (n != Z$ncol) {
    stop("full_exact requires a square block-kernel matrix (set sample_frac = 1).", call. = FALSE)
  }
  if (ncomp >= n) {
    stop("ncomp must be less than total sample count for full KEMA.", call. = FALSE)
  }

  M <- L + mu * Ls
  A_apply <- function(x) Z$apply(M %*% Z$crossprod(x))
  B_apply <- function(x) {
    Z$apply(Ld %*% Z$crossprod(x)) + lambda * as.matrix(x)
  }
  selected <- solve_kema_quotient(
    quotient,
    M,
    Ld,
    lambda,
    ncomp,
    control,
    original_A_apply = A_apply,
    original_B_apply = B_apply
  )
  selected$A_apply <- A_apply
  selected$B_apply <- B_apply
  selected$score_apply <- Z$apply
  selected$coefficient_space <- "full"
  selected$formulation <- "eq6_full_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_rekema <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  r <- Z$ncol
  if (r <= ncomp) {
    stop("REKEMA landmark rank (", r, ") must be > ncomp (", ncomp, ").", call. = FALSE)
  }

  M <- L + mu * Ls
  A_apply <- function(x) Z$crossprod(M %*% Z$apply(x))
  B_apply <- function(x) {
    Z$crossprod(Ld %*% Z$apply(x)) + lambda * as.matrix(x)
  }
  selected <- solve_kema_quotient(
    quotient,
    M,
    Ld,
    lambda,
    ncomp,
    control,
    original_A_apply = A_apply,
    original_B_apply = B_apply
  )
  selected$A_apply <- A_apply
  selected$B_apply <- B_apply
  selected$score_apply <- Z$apply
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_reduced_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_operator <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  r <- Z$ncol
  if (r <= ncomp) {
    stop("Operator rank (", r, ") must be > ncomp (", ncomp, ").", call. = FALSE)
  }

  M <- L + mu * Ls
  A_apply <- function(x) Z$crossprod(M %*% Z$apply(x))
  B_apply <- function(x) {
    Z$crossprod(Ld %*% Z$apply(x)) + lambda * as.matrix(x)
  }

  selected <- solve_kema_quotient(
    quotient,
    M,
    Ld,
    lambda,
    ncomp,
    control,
    original_A_apply = A_apply,
    original_B_apply = B_apply,
    operator = TRUE
  )
  selected$A_apply <- A_apply
  selected$B_apply <- B_apply
  selected$score_apply <- Z$apply
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_operator_exact"
  selected
}
