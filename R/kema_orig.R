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
#' @param lambda Positive Tikhonov regularization added to the RHS generalized
#'   matrix. A positive value is required because the partial generalized
#'   eigensolver uses that matrix as an inner-product metric.
#' @param backend Backend for the generalized eigensolver. One of `"auto"`,
#'   `"full_exact"`, `"reduced_exact"`, or `"operator_exact"`.
#' @param backend_control Optional list for backend auto-selection and fidelity
#'   diagnostics. Supported keys: `full_exact_max_n`, `reduced_exact_max_r`,
#'   `dense_exact_max_dim`,
#'   `fidelity_residual_tol`, `fidelity_orth_tol`, `fidelity_action`,
#'   `eigencore_tol`, `eigencore_maxit`, and `kernel_rank_tol`.
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
                                  ...) {
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
    backend_control = backend_control
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
      proclist[[i]] <- NULL
    } else if (is.function(pre_i)) {
      pdata[[i]]$x <- pre_i(Xi)
      proclist[[i]] <- pre_i
    } else {
      if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
        pre_i <- unserialize(serialize(pre_i, connection = NULL))
      }
      if (exists(
        "fit_transform",
        envir = asNamespace("multivarious"),
        mode = "function"
      )) {
        fitted <- multivarious::fit_transform(pre_i, Xi)
        pdata[[i]]$x <- fitted$transformed
        proclist[[i]] <- fitted$preproc
      } else {
        proc_template <- multivarious::prep(pre_i)
        transformed_x <- multivarious::init_transform(proc_template, Xi)
        pdata[[i]]$x <- transformed_x
        proc_attr <- attr(transformed_x, "preproc")
        proclist[[i]] <- if (is.null(proc_attr)) proc_template else proc_attr
      }
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

  sample_block_indices <- block_indices(pdata)
  block_indices_list <- split(sample_block_indices, row(sample_block_indices))
  proc <- multivarious::concat_pre_processors(proclist, block_indices_list)

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

  Ks <- compute_kernels(pdata, kernel, sample_frac = sample_frac, centre_kernel = FALSE)
  Ks_sparse <- lapply(Ks, function(k) {
    if (!methods::is(k, "sparseMatrix")) {
      k <- Matrix::Matrix(k, sparse = TRUE)
    } else {
      k <- k
    }
    if (sample_frac >= 1) {
      k <- (k + Matrix::t(k)) / 2
    }
    Matrix::drop0(k)
  })
  r_total <- sum(vapply(Ks_sparse, ncol, integer(1)))

  Z_cache <- NULL
  get_Z <- function() {
    if (is.null(Z_cache)) {
      Z_cache <<- Matrix::bdiag(Ks_sparse)
      if (nrow(Z_cache) != n_total) {
        stop("Kernel and graph dimensions do not match.", call. = FALSE)
      }
    }
    Z_cache
  }

  quotient <- build_kema_kernel_quotient(
    Ks = Ks_sparse,
    full_form = sample_frac >= 1,
    component_basis = component_basis,
    rank_tol = control$kernel_rank_tol
  )

  auto_backend <- if (identical(backend, "auto")) {
    choose_kema_orig_backend(sample_frac, n_total, r_total, control)
  } else {
    backend
  }
  candidates <- kema_orig_backend_candidates(
    requested_backend = backend,
    auto_backend = auto_backend,
    sample_frac = sample_frac,
    fidelity_action = control$fidelity_action
  )

  attempts <- list()
  eig <- NULL
  chosen_backend <- NULL
  selected_fidelity <- NULL
  for (candidate in candidates) {
    solved <- tryCatch(
      switch(
        candidate,
        full_exact = solve_kema_orig_full(
          Z = get_Z(),
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
          Z = get_Z(),
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
          Z = get_Z(),
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

  if (identical(eig$coefficient_space, "full")) {
    Z <- get_Z()
    scores <- Z %*% vectors
    k_idx <- get_block_indices(Ks, byrow = TRUE)
    v <- do.call(rbind, lapply(seq_along(pdata), function(i) {
      xi <- pdata[[i]]$x
      idx <- seq(k_idx[i, 1], k_idx[i, 2])
      alpha_i <- vectors[idx, , drop = FALSE]
      Matrix::crossprod(xi, alpha_i)
    }))
  } else {
    scores <- if (is.function(eig$score_apply)) {
      eig$score_apply(vectors)
    } else {
      get_Z() %*% vectors
    }
    landmark_counts <- vapply(Ks, ncol, integer(1))
    landmark_end <- cumsum(landmark_counts)
    landmark_start <- c(1L, head(landmark_end, -1L) + 1L)
    v <- do.call(rbind, lapply(seq_along(pdata), function(i) {
      xi <- pdata[[i]]$x
      idx <- seq(landmark_start[i], landmark_end[i])
      alpha_i <- vectors[idx, , drop = FALSE]
      Matrix::crossprod(xi, Ks[[i]] %*% alpha_i)
    }))
  }
  if (!is.matrix(scores) && !methods::is(scores, "Matrix")) {
    scores <- matrix(scores, ncol = 1L)
  }

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
    classes = "kema_orig"
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
    fidelity_residual_tol = 1e-6,
    fidelity_orth_tol = 1e-6,
    fidelity_action = "fallback",
    eigencore_tol = 1e-8,
    eigencore_maxit = 500L,
    kernel_rank_tol = sqrt(.Machine$double.eps)
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
  U <- Matrix::bdiag(lapply(pieces, function(x) Matrix::Matrix(x$U, sparse = TRUE)))
  V <- Matrix::bdiag(lapply(pieces, function(x) Matrix::Matrix(x$V, sparse = TRUE)))
  singular_values <- unlist(lapply(pieces, `[[`, "d"), use.names = FALSE)
  H <- U %*% Matrix::Diagonal(x = singular_values)

  component_basis <- Matrix::Matrix(component_basis, sparse = TRUE)
  if (nrow(component_basis) != nrow(U)) {
    stop("Kernel and graph-component dimensions do not match.", call. = FALSE)
  }
  UtC <- as.matrix(Matrix::crossprod(U, component_basis))
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

  list(
    H = H,
    V = V,
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
  )
}

#' Form and solve the nullspace-safe KEMA quotient pencil
#'
#' @keywords internal
#' @noRd
solve_kema_quotient <- function(quotient, M, Ld, lambda, ncomp, control,
                                force_partial = FALSE,
                                original_A_apply = NULL,
                                original_B_apply = NULL) {
  H <- quotient$H
  q <- ncol(H)
  available <- q - ncol(quotient$constraints)
  if (available < ncomp) {
    stop(
      "KEMA quotient has only ", available,
      " positive-dimensional direction(s) after nullspace deflation; need ",
      ncomp, ".",
      call. = FALSE
    )
  }

  A <- as.matrix(Matrix::crossprod(H, M %*% H))
  B <- as.matrix(
    Matrix::crossprod(H, Ld %*% H) +
      lambda * Matrix::Diagonal(q)
  )
  A <- (A + t(A)) / 2
  B <- (B + t(B)) / 2
  storage.mode(A) <- "double"
  storage.mode(B) <- "double"

  use_full <- !isTRUE(force_partial) && q <= control$dense_exact_max_dim
  fit <- if (use_full) {
    eigencore::eig_full(A, B = B, structure = eigencore::hermitian())
  } else {
    constraints <- quotient$constraints
    if (!ncol(constraints)) {
      constraints <- NULL
    }
    eigencore::eig_partial(
      A,
      B = B,
      k = ncomp,
      target = eigencore::smallest(),
      method = eigencore::lobpcg(
        maxit = control$eigencore_maxit,
        constraints = constraints
      ),
      tol = control$eigencore_tol,
      allow_dense_fallback = "never"
    )
  }

  cert <- eigencore::certificate(fit)
  if (!isTRUE(cert$passed) && !use_full && q <= control$dense_exact_max_dim) {
    fit <- eigencore::eig_full(
      A,
      B = B,
      structure = eigencore::hermitian()
    )
    cert <- eigencore::certificate(fit)
    use_full <- TRUE
  }
  if (!isTRUE(cert$passed)) {
    stop(
      "eigencore failed to certify the KEMA quotient solve (method: ",
      fit$method, ", max backward error: ",
      signif(cert$max_backward_error, 4), ").",
      call. = FALSE
    )
  }
  values <- as.numeric(Re(eigencore::values(fit)))
  vectors <- as.matrix(Re(eigencore::vectors(fit)))
  mapped_vectors <- quotient$V %*% vectors
  use_original_residual <- is.function(original_A_apply) &&
    is.function(original_B_apply)
  residuals_rel <- if (use_original_residual) {
    compute_kema_relative_residuals(
      values,
      mapped_vectors,
      A_apply = original_A_apply,
      B_apply = original_B_apply
    )
  } else {
    compute_kema_relative_residuals(
      values,
      vectors,
      A_apply = function(x) A %*% x,
      B_apply = function(x) B %*% x
    )
  }
  selected <- select_nontrivial_eigenpairs(
    values,
    if (use_original_residual) mapped_vectors else vectors,
    ncomp,
    spectral_scale = max(abs(values), 0),
    residuals_rel = residuals_rel,
    residual_tol = control$fidelity_residual_tol
  )
  if (!use_original_residual) {
    selected$vectors <- quotient$V %*% selected$vectors
  }
  selected$eigensolver <- if (use_full) {
    "eigencore_dense_full"
  } else {
    "eigencore_lobpcg_deflated"
  }
  selected$eigensolver_stats <- list(
    method = fit$method,
    certificate = cert,
    quotient_dimension = q,
    original_dimension = quotient$original_rank,
    kernel_rank_discarded = quotient$discarded_rank,
    graph_nullity = quotient$graph_nullity,
    nullity_deflated = quotient$nullity_deflated,
    full_spectrum = use_full
  )
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_full <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  n <- nrow(Z)
  if (n != ncol(Z)) {
    stop("full_exact requires a square block-kernel matrix (set sample_frac = 1).", call. = FALSE)
  }
  if (ncomp >= n) {
    stop("ncomp must be less than total sample count for full KEMA.", call. = FALSE)
  }

  M <- L + mu * Ls
  A_apply <- function(x) Z %*% (M %*% (Matrix::t(Z) %*% x))
  B_apply <- function(x) {
    Z %*% (Ld %*% (Matrix::t(Z) %*% x)) + lambda * x
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
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "full"
  selected$formulation <- "eq6_full_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_rekema <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  r <- ncol(Z)
  if (r <= ncomp) {
    stop("REKEMA landmark rank (", r, ") must be > ncomp (", ncomp, ").", call. = FALSE)
  }

  M <- L + mu * Ls
  A_apply <- function(x) Matrix::crossprod(Z, M %*% (Z %*% x))
  B_apply <- function(x) {
    Matrix::crossprod(Z, Ld %*% (Z %*% x)) + lambda * x
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
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_reduced_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_operator <- function(Z, quotient, L, Ls, Ld, mu, lambda, ncomp, control) {
  r <- ncol(Z)
  if (r <= ncomp) {
    stop("Operator rank (", r, ") must be > ncomp (", ncomp, ").", call. = FALSE)
  }

  M <- L + mu * Ls
  to_dense_input <- function(x, nrow_target) {
    if (is.matrix(x)) {
      return(matrix(as.numeric(x), nrow = nrow_target))
    }
    if (methods::is(x, "Matrix")) {
      return(matrix(as.numeric(as.matrix(x)), nrow = nrow_target))
    }
    matrix(as.numeric(x), nrow = nrow_target)
  }
  restore_shape <- function(x, template) {
    if (is.matrix(template) || methods::is(template, "Matrix")) {
      x
    } else {
      drop(x)
    }
  }

  A_apply <- function(x) {
    xin <- x
    X <- to_dense_input(xin, nrow_target = r)
    out <- Matrix::crossprod(Z, M %*% (Z %*% X))
    restore_shape(matrix(as.numeric(out), nrow = r), xin)
  }
  B_apply <- function(x) {
    xin <- x
    X <- to_dense_input(xin, nrow_target = r)
    out <- Matrix::crossprod(Z, Ld %*% (Z %*% X)) + lambda * X
    restore_shape(matrix(as.numeric(out), nrow = r), xin)
  }

  selected <- solve_kema_quotient(
    quotient,
    M,
    Ld,
    lambda,
    ncomp,
    control,
    force_partial = TRUE,
    original_A_apply = A_apply,
    original_B_apply = B_apply
  )
  selected$A_apply <- A_apply
  selected$B_apply <- B_apply
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_operator_exact"
  selected
}
