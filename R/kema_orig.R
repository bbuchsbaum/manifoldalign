#' Original Kernel Manifold Alignment (KEMA)
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
#' @param lambda Non-negative regularization added to RHS generalized matrix.
#'   Set to `0` for strict paper formulation.
#' @param backend Backend for the generalized eigensolver. One of `"auto"`,
#'   `"full_exact"`, `"reduced_exact"`, or `"operator_exact"`.
#' @param backend_control Optional list for backend auto-selection and fidelity
#'   diagnostics. Supported keys: `full_exact_max_n`, `reduced_exact_max_r`,
#'   `fidelity_residual_tol`, `fidelity_orth_tol`, `fidelity_action`,
#'   `primme_tol`, and `primme_method`.
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
                                  lambda = 0,
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
                                  lambda = 0,
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
  chk::chk_true(lambda >= 0)
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
  proclist <- vector("list", length(data))
  is_stateful_preproc <- inherits(preproc, "prepper") || inherits(preproc, "pre_processor")
  if (is_stateful_preproc) {
    preproc_list <- replicate(
      length(data),
      unserialize(serialize(preproc, connection = NULL)),
      simplify = FALSE
    )
  } else {
    preproc_list <- rep(list(preproc), length(data))
  }

  for (i in seq_along(data)) {
    pre_i <- preproc_list[[i]]
    if (inherits(pre_i, "prepper") || inherits(pre_i, "pre_processor")) {
      pre_i <- unserialize(serialize(pre_i, connection = NULL))
    }
    proc_template <- multivarious::prep(pre_i)
    transformed_x <- multivarious::init_transform(proc_template, data[[i]]$x)
    pdata[[i]]$x <- transformed_x
    proc_attr <- attr(transformed_x, "preproc")
    proclist[[i]] <- if (is.null(proc_attr)) proc_template else proc_attr
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
    g <- neighborweights::graph_weights(
      dom$x,
      weight_mode = "normalized",
      neighbor_mode = "knn",
      k = knn,
      type = "normal",
      sigma = sigma
    )
    A <- neighborweights::adjacency(g)
    if (!methods::is(A, "sparseMatrix")) {
      A <- Matrix::Matrix(A, sparse = TRUE)
    }
    as(A, "dgCMatrix")
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

  Ks <- compute_kernels(pdata, kernel, sample_frac = sample_frac, centre_kernel = FALSE)
  Ks_sparse <- lapply(Ks, function(k) {
    if (!methods::is(k, "sparseMatrix")) {
      Matrix::Matrix(k, sparse = TRUE)
    } else {
      k
    }
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
    if (identical(control$fidelity_action, "error")) {
      stop(fidelity_msg, call. = FALSE)
    }
    if (identical(control$fidelity_action, "warn") ||
        identical(control$fidelity_action, "fallback")) {
      warning(fidelity_msg, call. = FALSE)
    }
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
select_nontrivial_eigenpairs <- function(values, vectors, ncomp, tol = 1e-10) {
  non_trivial <- which(abs(values) > tol)
  if (length(non_trivial) >= ncomp) {
    take <- non_trivial[seq_len(ncomp)]
  } else {
    take <- seq_len(min(ncomp, length(values)))
  }
  if (length(take) < ncomp) {
    stop("Failed to extract enough eigenvectors (got ", length(take), ", need ", ncomp, ").", call. = FALSE)
  }
  list(values = values[take], vectors = vectors[, take, drop = FALSE])
}

#' @keywords internal
#' @noRd
default_kema_backend_control <- function() {
  list(
    full_exact_max_n = 1200L,
    reduced_exact_max_r = 5000L,
    fidelity_residual_tol = 5e-4,
    fidelity_orth_tol = 1e-3,
    fidelity_action = "warn",
    primme_tol = 1e-6,
    primme_method = "PRIMME_DEFAULT_MIN_MATVECS"
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
  out$fidelity_action <- match.arg(
    out$fidelity_action,
    c("warn", "error", "fallback", "ignore")
  )
  chk::chk_true(out$full_exact_max_n > 2)
  chk::chk_true(out$reduced_exact_max_r > 2)
  chk::chk_number(out$fidelity_residual_tol)
  chk::chk_true(out$fidelity_residual_tol > 0)
  chk::chk_number(out$fidelity_orth_tol)
  chk::chk_true(out$fidelity_orth_tol > 0)
  chk::chk_number(out$primme_tol)
  chk::chk_true(out$primme_tol > 0)
  if (!is.null(out$primme_method) && !is.character(out$primme_method)) {
    stop("backend_control$primme_method must be NULL or character.", call. = FALSE)
  }

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
compute_kema_fidelity <- function(values, vectors, A_apply, B_apply, residual_tol, orth_tol) {
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
  residuals_rel <- as.numeric(resid_num / resid_den)

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

#' @keywords internal
#' @noRd
kema_primme_options <- function(control) {
  opts <- list(tol = control$primme_tol)
  if (!is.null(control$primme_method) && nzchar(control$primme_method)) {
    opts$method <- control$primme_method
  }
  opts
}

#' @keywords internal
#' @noRd
solve_kema_orig_full <- function(Z, L, Ls, Ld, mu, lambda, ncomp, control) {
  n <- nrow(Z)
  if (n != ncol(Z)) {
    stop("full_exact requires a square block-kernel matrix (set sample_frac = 1).", call. = FALSE)
  }
  if (ncomp >= n) {
    stop("ncomp must be less than total sample count for full KEMA.", call. = FALSE)
  }

  A <- Z %*% (L + mu * Ls) %*% Matrix::t(Z)
  B <- Z %*% Ld %*% Matrix::t(Z) + lambda * Matrix::Diagonal(n)
  decomp <- do.call(
    PRIMME::eigs_sym,
    c(
      list(
        A = A,
        NEig = min(ncomp + 1, n - 1),
        which = "SA",
        B = B
      ),
      kema_primme_options(control)
    )
  )
  selected <- select_nontrivial_eigenpairs(decomp$values, decomp$vectors, ncomp)
  selected$A_apply <- function(x) A %*% x
  selected$B_apply <- function(x) B %*% x
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "full"
  selected$formulation <- "eq6_full_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_rekema <- function(Z, L, Ls, Ld, mu, lambda, ncomp, control) {
  r <- ncol(Z)
  if (r <= ncomp) {
    stop("REKEMA landmark rank (", r, ") must be > ncomp (", ncomp, ").", call. = FALSE)
  }

  A <- Matrix::crossprod(Z, (L + mu * Ls) %*% Z)
  B <- Matrix::crossprod(Z, Ld %*% Z) + lambda * Matrix::Diagonal(r)
  decomp <- do.call(
    PRIMME::eigs_sym,
    c(
      list(
        A = A,
        NEig = min(ncomp + 1, r - 1),
        which = "SA",
        B = B
      ),
      kema_primme_options(control)
    )
  )
  selected <- select_nontrivial_eigenpairs(decomp$values, decomp$vectors, ncomp)
  selected$A_apply <- function(x) A %*% x
  selected$B_apply <- function(x) B %*% x
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_reduced_exact"
  selected
}

#' @keywords internal
#' @noRd
solve_kema_orig_operator <- function(Z, L, Ls, Ld, mu, lambda, ncomp, control) {
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

  decomp <- do.call(
    PRIMME::eigs_sym,
    c(
      list(
        A = A_apply,
        NEig = min(ncomp + 1, r - 1),
        which = "SA",
        B = B_apply,
        isreal = TRUE,
        n = r
      ),
      kema_primme_options(control)
    )
  )

  selected <- select_nontrivial_eigenpairs(decomp$values, decomp$vectors, ncomp)
  selected$A_apply <- A_apply
  selected$B_apply <- B_apply
  selected$score_apply <- function(x) Z %*% x
  selected$coefficient_space <- "reduced"
  selected$formulation <- "eq10_operator_exact"
  selected
}
