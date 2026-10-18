# Stimulus-weighted SUMCOR alignment ---------------------------------------
#
# Mathematical contract (all quantities are per subject s unless stated):
#   u_ser   = phi_s(x_ser)                      fitted linear brain basis
#   ubar_se = mean_r u_ser                      event means (numerator)
#   C_s     = (1 - rho) Chat_s + rho tr(Chat_s)/p I
#   Chat_s  = (1/m_s) sum_e (w_e / R_se) sum_r u_ser u_ser'   (presentations)
#   H       = I + beta sum_m theta_m K_m,  K_m = F_m (F_m'F_m + eta_m I)^-1 F_m'
#   B_st    = ubar_s' W^1/2 H W^1/2 ubar_t / m_st              (overlap-normalised)
#   maximise sum_{s<t} tr(A_s' B_st A_t)  subject to  A_s' C_s A_s = I
#   encoder f_s(x) = ((x - mu_s) / sd_s) P_s A_s
# Event weights w sum to one over the event universe; m_s is the weight of
# the events subject s observed and m_st the weight of the shared events.

# -- Validation helpers ------------------------------------------------------

.ssc_number <- function(x, name, lower = -Inf, upper = Inf,
                        strict_lower = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop("`", name, "` must be a single finite number.", call. = FALSE)
  }
  if (x < lower || (strict_lower && x <= lower) || x > upper) {
    bound <- if (strict_lower) paste0("> ", lower) else paste0(">= ", lower)
    if (is.finite(upper)) bound <- paste0(bound, " and <= ", upper)
    stop("`", name, "` must be ", bound, ".", call. = FALSE)
  }
  as.numeric(x)
}

.ssc_int <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > 1e-8 || x < lower) {
    stop("`", name, "` must be a single integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(round(x))
}

.ssc_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.ssc_with_preserved_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  expr
}

# -- Control ------------------------------------------------------------------

#' Control settings for stimulus-weighted SUMCOR alignment
#'
#' @param beta Nonnegative weight of the stimulus-feature metric relative to
#'   exact event correspondence. `beta = 0` gives the correspondence-only
#'   model exactly; supplied features are then not used.
#' @param family_weights Optional named nonnegative weights over feature
#'   families (`theta_m`). They are normalised to sum to one over the families
#'   supplied. `NULL` gives equal weights.
#' @param feature_ridge Positive ridge for each feature kernel, relative to the
#'   mean squared singular value of the preprocessed feature block:
#'   `eta = feature_ridge * mean(sigma^2)`.
#' @param feature_rank Optional maximum number of singular components kept per
#'   feature block. `NULL` keeps every component with nonnegligible weight.
#' @param feature_scale Logical; scale feature columns to unit variance after
#'   centring (constant columns are dropped).
#' @param kernel_normalization `"max_eigenvalue"` rescales each family kernel so
#'   its largest eigenvalue is one; `"none"` keeps `K_m` as defined.
#' @param cross_modal Logical; if `FALSE` (default) each family kernel is block
#'   diagonal by event modality, so features never relate events from
#'   different modalities.
#' @param covariance_shrinkage Shrinkage `rho` in `[0, 1]` of each subject's
#'   presentation covariance toward a scaled identity.
#' @param brain_basis Linear brain representation: `"pca"` (weighted PCA over
#'   all features), `"parcel_pca"` (PCA within supplied parcels), or
#'   `"identity"` (centring/scaling only).
#' @param brain_rank Maximum number of PCA components per subject for
#'   `"pca"`.
#' @param parcel_rank Maximum number of PCA components per parcel for
#'   `"parcel_pca"`.
#' @param brain_scale Logical; scale brain features to unit (event-weighted)
#'   variance before the basis is fitted.
#' @param event_weights Event weighting: `"uniform"`, `"modality"` (equal total
#'   weight per modality), `"clip"` (equal total weight per clip), or a named
#'   nonnegative numeric vector over event identifiers.
#' @param min_overlap Minimum shared event weight (a fraction of the total) for
#'   a pair of subjects to contribute to the objective.
#' @param max_iter Maximum number of sequential polar sweeps.
#' @param tol Relative objective-change tolerance for convergence.
#' @param grad_tol Relative projected-gradient (stationarity) tolerance for
#'   convergence.
#' @param rank_tol Relative singular-value threshold below which a polar update
#'   direction is treated as unsupported and completed deterministically.
#' @param identity_max_dim Largest feature count allowed with
#'   `brain_basis = "identity"`.
#' @param store_basis Logical; keep each subject's fitted basis matrix in the
#'   result. Prediction does not need it.
#' @param verbose Logical; print sweep progress. Logging never changes the
#'   computation.
#' @return An object of class `stimulus_sumcor_control`.
#' @export
stimulus_sumcor_control <- function(
  beta = 0.5,
  family_weights = NULL,
  feature_ridge = 1,
  feature_rank = NULL,
  feature_scale = TRUE,
  kernel_normalization = c("max_eigenvalue", "none"),
  cross_modal = FALSE,
  covariance_shrinkage = 0.1,
  brain_basis = c("pca", "parcel_pca", "identity"),
  brain_rank = 50L,
  parcel_rank = 5L,
  brain_scale = FALSE,
  event_weights = "uniform",
  min_overlap = 0.25,
  max_iter = 500L,
  tol = 1e-9,
  grad_tol = 1e-5,
  rank_tol = 1e-8,
  identity_max_dim = 2000L,
  store_basis = TRUE,
  verbose = FALSE
) {
  kernel_normalization <- match.arg(kernel_normalization)
  brain_basis <- match.arg(brain_basis)
  beta <- .ssc_number(beta, "beta", lower = 0)
  feature_ridge <- .ssc_number(feature_ridge, "feature_ridge", lower = 0,
                               strict_lower = TRUE)
  if (!is.null(feature_rank)) feature_rank <- .ssc_int(feature_rank, "feature_rank")
  feature_scale <- .ssc_flag(feature_scale, "feature_scale")
  cross_modal <- .ssc_flag(cross_modal, "cross_modal")
  covariance_shrinkage <- .ssc_number(covariance_shrinkage,
                                      "covariance_shrinkage", 0, 1)
  brain_rank <- .ssc_int(brain_rank, "brain_rank")
  parcel_rank <- .ssc_int(parcel_rank, "parcel_rank")
  brain_scale <- .ssc_flag(brain_scale, "brain_scale")
  min_overlap <- .ssc_number(min_overlap, "min_overlap", 0, 1)
  max_iter <- .ssc_int(max_iter, "max_iter")
  tol <- .ssc_number(tol, "tol", lower = 0)
  grad_tol <- .ssc_number(grad_tol, "grad_tol", lower = 0)
  rank_tol <- .ssc_number(rank_tol, "rank_tol", 0, 1, strict_lower = TRUE)
  identity_max_dim <- .ssc_int(identity_max_dim, "identity_max_dim")
  store_basis <- .ssc_flag(store_basis, "store_basis")
  verbose <- .ssc_flag(verbose, "verbose")

  if (!is.null(family_weights)) {
    if (!is.numeric(family_weights) || !length(family_weights) ||
        any(!is.finite(family_weights)) || any(family_weights < 0) ||
        is.null(names(family_weights)) || any(!nzchar(names(family_weights))) ||
        anyDuplicated(names(family_weights))) {
      stop("`family_weights` must be a named, finite, nonnegative numeric ",
           "vector with unique names.", call. = FALSE)
    }
    if (sum(family_weights) <= 0) {
      stop("`family_weights` must not all be zero.", call. = FALSE)
    }
  }

  if (is.character(event_weights)) {
    if (length(event_weights) != 1L ||
        !event_weights %in% c("uniform", "modality", "clip")) {
      stop("`event_weights` must be \"uniform\", \"modality\", \"clip\", or a ",
           "named numeric vector.", call. = FALSE)
    }
  } else if (is.numeric(event_weights)) {
    if (!length(event_weights) || any(!is.finite(event_weights)) ||
        any(event_weights < 0) || is.null(names(event_weights)) ||
        anyNA(names(event_weights)) || any(!nzchar(names(event_weights))) ||
        anyDuplicated(names(event_weights))) {
      stop("Numeric `event_weights` must be finite, nonnegative, and uniquely ",
           "named by event identifier.", call. = FALSE)
    }
  } else {
    stop("`event_weights` must be a string or a named numeric vector.",
         call. = FALSE)
  }

  structure(
    list(
      beta = beta,
      family_weights = family_weights,
      feature_ridge = feature_ridge,
      feature_rank = feature_rank,
      feature_scale = feature_scale,
      kernel_normalization = kernel_normalization,
      cross_modal = cross_modal,
      covariance_shrinkage = covariance_shrinkage,
      brain_basis = brain_basis,
      brain_rank = brain_rank,
      parcel_rank = parcel_rank,
      brain_scale = brain_scale,
      event_weights = event_weights,
      min_overlap = min_overlap,
      max_iter = max_iter,
      tol = tol,
      grad_tol = grad_tol,
      rank_tol = rank_tol,
      identity_max_dim = identity_max_dim,
      store_basis = store_basis,
      verbose = verbose
    ),
    class = "stimulus_sumcor_control"
  )
}

# Every path re-runs the constructor, so nested or list input is always
# validated and unknown fields are always rejected.
.ssc_resolve_control <- function(control) {
  if (is.null(control)) return(stimulus_sumcor_control())
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or stimulus_sumcor_control().",
         call. = FALSE)
  }
  values <- unclass(control)
  if (length(values) &&
      (is.null(names(values)) || any(!nzchar(names(values))))) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(values), names(formals(stimulus_sumcor_control)))
  if (length(unknown)) {
    stop("Unknown stimulus_sumcor control field(s): ",
         paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  do.call(stimulus_sumcor_control, values)
}

# -- Input resolution ----------------------------------------------------------

.ssc_as_numeric_matrix <- function(X, name) {
  if (is.data.frame(X)) X <- as.matrix(X)
  if (inherits(X, "Matrix")) X <- as.matrix(X)
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`", name, "` must be a numeric matrix, data.frame, or Matrix.",
         call. = FALSE)
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("`", name, "` must have at least one row and one column.",
         call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop("`", name, "` must contain only finite values.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  X
}

.ssc_design_column <- function(design, column, name, domain) {
  if (is.null(column)) return(NULL)
  if (!is.character(column) || length(column) != 1L || !nzchar(column)) {
    stop("`", name, "` must name one design column.", call. = FALSE)
  }
  if (!column %in% names(design)) {
    stop("Domain `", domain, "`: design column `", column, "` (", name,
         ") is missing.", call. = FALSE)
  }
  design[[column]]
}

# Resolve one subject's observation table. Correspondence is carried only by
# keyed identifiers; row counts and row order are never used to infer it.
.ssc_resolve_domain <- function(X, design, ids, domain) {
  X <- .ssc_as_numeric_matrix(X, paste0("data[[\"", domain, "\"]]"))
  n <- nrow(X)
  if (is.null(design)) {
    stop("Domain `", domain, "`: a design table with event identifiers is ",
         "required.", call. = FALSE)
  }
  design <- as.data.frame(design, stringsAsFactors = FALSE)
  if (nrow(design) != n) {
    stop("Domain `", domain, "`: design has ", nrow(design), " rows but the ",
         "data have ", n, ".", call. = FALSE)
  }

  event <- .ssc_design_column(design, ids$event_id, "event_id", domain)
  if (is.null(event)) {
    stop("`event_id` must name the design column holding event identifiers.",
         call. = FALSE)
  }
  event <- as.character(event)
  if (any(!is.na(event) & !nzchar(event))) {
    stop("Domain `", domain, "`: event identifiers must not be empty strings.",
         call. = FALSE)
  }

  obs <- .ssc_design_column(design, ids$observation_id, "observation_id", domain)
  if (is.null(obs)) {
    obs <- rownames(X)
    if (is.null(obs) || anyDuplicated(obs) || anyNA(obs)) {
      obs <- as.character(seq_len(n))
    }
  }
  obs <- as.character(obs)
  if (anyNA(obs) || anyDuplicated(obs)) {
    stop("Domain `", domain, "`: observation identifiers must be unique and ",
         "non-missing.", call. = FALSE)
  }

  rep_raw <- .ssc_design_column(design, ids$repeat_id, "repeat_id", domain)
  if (is.null(rep_raw)) {
    repeat_id <- rep(NA_character_, n)
    labelled <- which(!is.na(event))
    occurrence <- stats::ave(seq_along(labelled), event[labelled],
                             FUN = seq_along)
    repeat_id[labelled] <- as.character(occurrence)
  } else {
    repeat_id <- as.character(rep_raw)
    if (any(!is.na(event) & is.na(repeat_id))) {
      stop("Domain `", domain, "`: every event-labelled observation needs a ",
           "repeat identifier.", call. = FALSE)
    }
    keys <- paste(event, repeat_id, sep = "\r")[!is.na(event)]
    if (anyDuplicated(keys)) {
      stop("Domain `", domain, "`: (event_id, repeat_id) pairs must be unique.",
           call. = FALSE)
    }
  }

  run <- .ssc_design_column(design, ids$run_id, "run_id", domain)
  modality <- .ssc_design_column(design, ids$modality, "modality", domain)
  clip <- .ssc_design_column(design, ids$clip, "clip", domain)

  table <- data.frame(
    row = seq_len(n),
    observation_id = obs,
    event_id = event,
    repeat_id = repeat_id,
    run_id = if (is.null(run)) NA_character_ else as.character(run),
    used = !is.na(event),
    stringsAsFactors = FALSE
  )
  list(
    X = X,
    observations = table,
    modality = if (is.null(modality)) NULL else as.character(modality),
    clip = if (is.null(clip)) NULL else as.character(clip)
  )
}

.ssc_event_attribute <- function(domains, events_df, field, universe) {
  value <- rep(NA_character_, length(universe))
  names(value) <- universe
  assign_value <- function(ev, val, source) {
    keep <- !is.na(ev) & !is.na(val)
    if (!any(keep)) return(invisible(NULL))
    ev <- ev[keep]
    val <- as.character(val[keep])
    pairs <- unique(data.frame(ev = ev, val = val, stringsAsFactors = FALSE))
    if (anyDuplicated(pairs$ev)) {
      bad <- pairs$ev[duplicated(pairs$ev)][1L]
      stop("Event `", bad, "` has more than one ", field, " in ", source, ".",
           call. = FALSE)
    }
    idx <- match(pairs$ev, universe)
    pairs <- pairs[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
    clash <- !is.na(value[idx]) & value[idx] != pairs$val
    if (any(clash)) {
      stop("Event `", pairs$ev[clash][1L], "` has inconsistent ", field,
           " values across inputs.", call. = FALSE)
    }
    value[idx] <<- pairs$val
    invisible(NULL)
  }
  for (dom in names(domains)) {
    observed <- domains[[dom]][[field]]
    if (!is.null(observed)) {
      assign_value(domains[[dom]]$observations$event_id, observed,
                   paste0("domain `", dom, "`"))
    }
  }
  if (!is.null(events_df) && field %in% names(events_df)) {
    assign_value(as.character(events_df$event_id), events_df[[field]],
                 "`events`")
  }
  unname(value)
}

.ssc_event_weights <- function(spec, events) {
  n <- nrow(events)
  if (is.numeric(spec)) {
    missing <- setdiff(events$event_id, names(spec))
    if (length(missing)) {
      stop("Numeric `event_weights` lacks ", length(missing), " event(s), e.g. `",
           missing[1L], "`.", call. = FALSE)
    }
    w <- as.numeric(spec[events$event_id])
  } else if (identical(spec, "uniform")) {
    w <- rep(1, n)
  } else {
    group <- events[[spec]]
    if (anyNA(group)) {
      stop("`event_weights = \"", spec, "\"` requires a ", spec,
           " for every event.", call. = FALSE)
    }
    counts <- table(group)
    w <- 1 / (length(counts) * as.numeric(counts[group]))
  }
  if (sum(w) <= 0) stop("Event weights sum to zero.", call. = FALSE)
  w / sum(w)
}

# -- Brain basis ---------------------------------------------------------------

.ssc_svd <- function(M, k) {
  k <- min(as.integer(k), nrow(M), ncol(M))
  if (k < 1L) return(list(d = numeric(0), u = matrix(0, nrow(M), 0L),
                          v = matrix(0, ncol(M), 0L)))
  small <- min(dim(M))
  if (small <= 300L) {
    s <- svd(M, nu = k, nv = k)
    out <- list(d = s$d[seq_len(k)], u = s$u, v = s$v)
  } else if (small <= 6000L) {
    # Eigendecomposition of the smaller Gram matrix: one level-3 BLAS product
    # plus a small symmetric eigenproblem. The long-side vectors are recovered
    # by projection and re-orthonormalised (span preserved).
    rows_short <- nrow(M) <= ncol(M)
    G <- if (rows_short) tcrossprod(M) else crossprod(M)
    e <- eigen(G, symmetric = TRUE)
    lam <- pmax(e$values[seq_len(k)], 0)
    d <- sqrt(lam)
    short <- e$vectors[, seq_len(k), drop = FALSE]
    positive <- d > 0
    long <- matrix(0, if (rows_short) ncol(M) else nrow(M), k)
    if (any(positive)) {
      proj <- if (rows_short) crossprod(M, short[, positive, drop = FALSE]) else
        M %*% short[, positive, drop = FALSE]
      long[, positive] <- qr.Q(qr(sweep(proj, 2L, d[positive], "/")))
    }
    out <- if (rows_short) list(d = d, u = short, v = long) else
      list(d = d, u = long, v = short)
  } else {
    start <- rep(1 / sqrt(ncol(M)), ncol(M))
    s <- .ssc_with_preserved_rng(
      irlba::irlba(M, nv = k, nu = k, v = start, tol = 1e-8, maxit = 2000L)
    )
    out <- list(d = s$d, u = s$u, v = s$v)
  }
  # Deterministic signs: largest-magnitude loading of each component positive.
  if (k > 0L) {
    pivot <- apply(abs(out$v), 2L, which.max)
    flip <- sign(out$v[cbind(pivot, seq_len(k))])
    flip[flip == 0] <- 1
    out$v <- sweep(out$v, 2L, flip, "*")
    out$u <- sweep(out$u, 2L, flip, "*")
  }
  out
}

.ssc_fit_basis <- function(X, weights, control, parcels, domain) {
  V <- ncol(X)
  center <- colSums(X * weights)
  Xc <- sweep(X, 2L, center, "-")
  scale <- NULL
  if (isTRUE(control$brain_scale)) {
    scale <- sqrt(colSums(Xc^2 * weights))
    scale[!is.finite(scale) | scale <= 0] <- 1
    Xc <- sweep(Xc, 2L, scale, "/")
  }
  Xw <- Xc * sqrt(weights)
  n_eff <- sum(weights > 0)

  if (identical(control$brain_basis, "identity")) {
    if (V > control$identity_max_dim) {
      stop("Domain `", domain, "`: brain_basis = \"identity\" with ", V,
           " features exceeds `identity_max_dim = ", control$identity_max_dim,
           "`.", call. = FALSE)
    }
    P <- diag(1, V)
    sv <- svd(Xw, nu = 0L, nv = 0L)$d
    info <- list(type = "identity", singular_values = sv)
  } else if (identical(control$brain_basis, "pca")) {
    k <- min(control$brain_rank, V, max(n_eff - 1L, 1L))
    s <- .ssc_svd(Xw, k)
    keep <- s$d > control$rank_tol * max(s$d, 0)
    P <- s$v[, keep, drop = FALSE]
    info <- list(type = "pca", singular_values = s$d[keep],
                 dropped_degenerate = sum(!keep))
  } else {
    if (is.null(parcels)) {
      stop("Domain `", domain, "`: brain_basis = \"parcel_pca\" requires ",
           "`parcels`.", call. = FALSE)
    }
    if (length(parcels) != V) {
      stop("Domain `", domain, "`: `parcels` must have one entry per feature.",
           call. = FALSE)
    }
    labels <- as.character(parcels)
    levels <- sort(unique(labels[!is.na(labels)]), method = "radix")
    if (!length(levels)) {
      stop("Domain `", domain, "`: `parcels` assigns no feature to a parcel.",
           call. = FALSE)
    }
    rows <- integer(0)
    cols <- integer(0)
    vals <- numeric(0)
    svals <- numeric(0)
    parcel_of <- character(0)
    offset <- 0L
    for (lev in levels) {
      idx <- which(labels == lev)
      k <- min(control$parcel_rank, length(idx), max(n_eff - 1L, 1L))
      s <- .ssc_svd(Xw[, idx, drop = FALSE], k)
      keep <- s$d > control$rank_tol * max(s$d, 0)
      if (!any(keep)) next
      v <- s$v[, keep, drop = FALSE]
      kk <- ncol(v)
      rows <- c(rows, rep(idx, kk))
      cols <- c(cols, rep(offset + seq_len(kk), each = length(idx)))
      vals <- c(vals, as.numeric(v))
      svals <- c(svals, s$d[keep])
      parcel_of <- c(parcel_of, rep(lev, kk))
      offset <- offset + kk
    }
    if (!offset) {
      stop("Domain `", domain, "`: no parcel has a nondegenerate component.",
           call. = FALSE)
    }
    P <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(V, offset))
    info <- list(type = "parcel_pca", singular_values = svals,
                 parcel = parcel_of, n_parcels = length(levels))
  }
  U <- as.matrix(Xc %*% P)
  list(center = center, scale = scale, P = P, U = U, info = info)
}

# -- Stimulus kernels -----------------------------------------------------------

.ssc_kernel_factor <- function(features, events, control) {
  universe <- events$event_id
  n <- length(universe)
  empty <- list(L = NULL, families = NULL, dropped_rows = integer(0))
  if (is.null(features) || control$beta == 0) return(empty)
  if (!is.list(features) || is.data.frame(features) || !length(features)) {
    stop("`stimulus_features` must be a named list of feature matrices.",
         call. = FALSE)
  }
  fam_names <- names(features)
  if (is.null(fam_names) || anyNA(fam_names) || any(!nzchar(fam_names)) ||
      anyDuplicated(fam_names)) {
    stop("`stimulus_features` must have unique, non-empty family names.",
         call. = FALSE)
  }
  theta <- control$family_weights
  if (is.null(theta)) {
    theta <- stats::setNames(rep(1, length(fam_names)), fam_names)
  } else {
    if (!setequal(names(theta), fam_names)) {
      stop("`family_weights` names must match the stimulus feature families.",
           call. = FALSE)
    }
    theta <- theta[fam_names]
  }
  theta <- theta / sum(theta)

  factors <- list()
  diagnostics <- list()
  dropped <- integer(0)
  for (fam in fam_names) {
    F <- features[[fam]]
    if (is.data.frame(F)) F <- as.matrix(F)
    if (inherits(F, "Matrix")) F <- as.matrix(F)
    if (!is.matrix(F) || !is.numeric(F)) {
      stop("Feature family `", fam, "` must be a numeric matrix.", call. = FALSE)
    }
    ids <- rownames(F)
    if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
      stop("Feature family `", fam, "` needs unique row names giving event ",
           "identifiers.", call. = FALSE)
    }
    if (any(!is.finite(F))) {
      stop("Feature family `", fam, "` contains non-finite values.",
           call. = FALSE)
    }
    keep_rows <- ids %in% universe
    dropped[[fam]] <- sum(!keep_rows)
    F <- F[keep_rows, , drop = FALSE]
    if (nrow(F) < 2L) {
      stop("Feature family `", fam, "` describes fewer than two fitted events.",
           call. = FALSE)
    }
    idx <- match(rownames(F), universe)
    F <- sweep(F, 2L, colMeans(F), "-")
    sds <- sqrt(colSums(F^2) / max(nrow(F) - 1L, 1L))
    informative <- sds > 1e-12 * max(sds, 0)
    if (!any(informative)) {
      stop("Feature family `", fam, "` has no non-constant column.",
           call. = FALSE)
    }
    F <- F[, informative, drop = FALSE]
    if (isTRUE(control$feature_scale)) {
      F <- sweep(F, 2L, sds[informative], "/")
    }

    block_of <- if (!isTRUE(control$cross_modal) && !all(is.na(events$modality))) {
      mod <- events$modality[idx]
      if (anyNA(mod)) {
        stop("Feature family `", fam, "`: every described event needs a ",
             "modality when `cross_modal = FALSE`.", call. = FALSE)
      }
      mod
    } else {
      rep("all", length(idx))
    }
    blocks <- list()
    for (b in sort(unique(block_of), method = "radix")) {
      rows_b <- which(block_of == b)
      if (length(rows_b) < 2L) next
      Fb <- F[rows_b, , drop = FALSE]
      s <- svd(Fb)
      positive <- s$d > 1e-10 * max(s$d, 0)
      if (!any(positive)) next
      sig2 <- s$d[positive]^2
      eta <- control$feature_ridge * mean(sig2)
      kappa <- sig2 / (sig2 + eta)
      U <- s$u[, positive, drop = FALSE]
      if (!is.null(control$feature_rank) && length(kappa) > control$feature_rank) {
        kappa <- kappa[seq_len(control$feature_rank)]
        U <- U[, seq_len(control$feature_rank), drop = FALSE]
      }
      nb <- length(rows_b)
      cbar <- sum(kappa) / nb
      isotropy <- sqrt((sum((kappa - cbar)^2) + (nb - length(kappa)) * cbar^2) /
                         sum(kappa^2))
      blocks[[b]] <- list(rows = idx[rows_b], U = U, kappa = kappa, eta = eta,
                          n_events = nb, isotropy_deviation = isotropy)
    }
    if (!length(blocks)) {
      stop("Feature family `", fam, "` yields no usable kernel block.",
           call. = FALSE)
    }
    max_kappa <- max(vapply(blocks, function(b) max(b$kappa), numeric(1)))
    norm <- if (identical(control$kernel_normalization, "max_eigenvalue")) {
      max_kappa
    } else {
      1
    }
    for (b in names(blocks)) {
      blk <- blocks[[b]]
      Lb <- matrix(0, n, length(blk$kappa))
      Lb[blk$rows, ] <- sweep(blk$U, 2L, sqrt(blk$kappa / norm), "*")
      factors[[length(factors) + 1L]] <- sqrt(control$beta * theta[[fam]]) * Lb
      diagnostics[[length(diagnostics) + 1L]] <- data.frame(
        family = fam,
        block = b,
        n_events = blk$n_events,
        n_features = ncol(F),
        rank = length(blk$kappa),
        eta = blk$eta,
        effective_rank = sum(blk$kappa),
        isotropy_deviation = blk$isotropy_deviation,
        theta = theta[[fam]],
        stringsAsFactors = FALSE
      )
      # Centring alone leaves a deviation of about 1/sqrt(n) for an isotropic
      # kernel, so the threshold allows for it.
      if (blk$isotropy_deviation < max(0.1, 1.5 / sqrt(blk$n_events))) {
        warning(
          "Feature family `", fam, "` (block `", b, "`) gives a kernel close ",
          "to a multiple of the identity; it will change the fit very little. ",
          "Reduce its effective rank or increase `feature_ridge`.",
          call. = FALSE
        )
      }
    }
  }
  diag_df <- do.call(rbind, diagnostics)
  rownames(diag_df) <- NULL
  list(L = do.call(cbind, factors), families = diag_df, dropped_rows = dropped,
       theta = theta)
}

.ssc_apply_H <- function(L, Y) {
  if (is.null(L)) return(Y)
  Y + L %*% crossprod(L, Y)
}

# -- Polar update with deterministic completion ---------------------------------

.ssc_polar <- function(G, Q_prev = NULL, rank_tol = 1e-8) {
  p <- nrow(G)
  d <- ncol(G)
  s <- svd(G)
  smax <- if (length(s$d)) max(s$d) else 0
  k <- if (smax > 0) sum(s$d > rank_tol * smax) else 0L
  if (k == d) {
    return(list(Q = s$u %*% t(s$v), rank = d, singular_values = s$d))
  }
  Uk <- s$u[, seq_len(k), drop = FALSE]
  Vk <- s$v[, seq_len(k), drop = FALSE]
  Vperp <- s$v[, seq.int(k + 1L, d), drop = FALSE]
  complement <- function(M) M - Uk %*% crossprod(Uk, M)
  candidate <- NULL
  if (!is.null(Q_prev)) {
    candidate <- complement(Q_prev %*% Vperp)
    cs <- svd(candidate, nu = 0L, nv = 0L)$d
    if (!length(cs) || min(cs) < 1e-8) candidate <- NULL
  }
  if (is.null(candidate)) {
    B <- complement(diag(1, p))
    qrB <- qr(B, tol = 1e-10)
    candidate <- B[, qrB$pivot[seq_len(d - k)], drop = FALSE]
  }
  cs <- svd(candidate)
  N <- cs$u %*% t(cs$v)
  list(Q = Uk %*% t(Vk) + N %*% t(Vperp), rank = k, singular_values = s$d)
}

# -- Solver ----------------------------------------------------------------------

.ssc_pair_weights <- function(masks, w, min_overlap) {
  S <- length(masks)
  overlap <- matrix(0, S, S)
  for (s in seq_len(S)) for (t in seq_len(S)) {
    overlap[s, t] <- sum(w[masks[[s]] & masks[[t]]])
  }
  lambda <- ifelse(overlap >= min_overlap & overlap > 0, 1 / overlap, 0)
  diag(lambda) <- 0
  list(overlap = overlap, lambda = lambda)
}

.ssc_objective <- function(Y, HY, lambda) {
  S <- length(Y)
  total <- 0
  for (s in seq_len(S - 1L)) for (t in seq.int(s + 1L, S)) {
    if (lambda[s, t] > 0) total <- total + lambda[s, t] * sum(Y[[s]] * HY[[t]])
  }
  total
}

.ssc_gradient_input <- function(s, HY, lambda) {
  acc <- NULL
  for (t in seq_along(HY)) {
    if (t == s || lambda[s, t] == 0) next
    term <- lambda[s, t] * HY[[t]]
    acc <- if (is.null(acc)) term else acc + term
  }
  acc
}

.ssc_solve <- function(Ut, L, lambda, d, control, domain_names) {
  S <- length(Ut)
  n_ev <- nrow(Ut[[1L]])
  off <- lambda[upper.tri(lambda)]
  full <- all(off > 0) && (max(off) - min(off)) <= 1e-12 * max(off)

  # Order-invariant MAXVAR-style start: leading left singular vectors of the
  # concatenated whitened event scores (column order does not change them).
  stacked <- do.call(cbind, Ut)
  T0 <- .ssc_svd(stacked, d)$u
  if (ncol(T0) < d) T0 <- cbind(T0, matrix(0, n_ev, d - ncol(T0)))
  HT0 <- .ssc_apply_H(L, T0)
  Q <- vector("list", S)
  supported <- integer(S)
  for (s in seq_len(S)) {
    pol <- .ssc_polar(crossprod(Ut[[s]], HT0), NULL, control$rank_tol)
    Q[[s]] <- pol$Q
    supported[s] <- pol$rank
  }
  Y <- Map(function(U, Qs) U %*% Qs, Ut, Q)
  HY <- lapply(Y, function(y) .ssc_apply_H(L, y))
  objective <- .ssc_objective(Y, HY, lambda)

  history <- list()
  converged <- FALSE
  monotone <- TRUE
  status <- "max_iter"
  it <- 0L
  for (it in seq_len(control$max_iter)) {
    previous <- objective
    if (full) {
      HYsum <- Reduce(`+`, HY)
    }
    for (s in seq_len(S)) {
      R <- if (full) off[1L] * (HYsum - HY[[s]]) else .ssc_gradient_input(s, HY, lambda)
      if (is.null(R)) next
      pol <- .ssc_polar(crossprod(Ut[[s]], R), Q[[s]], control$rank_tol)
      Q[[s]] <- pol$Q
      supported[s] <- pol$rank
      Y[[s]] <- Ut[[s]] %*% Q[[s]]
      HYs <- .ssc_apply_H(L, Y[[s]])
      if (full) HYsum <- HYsum + HYs - HY[[s]]
      HY[[s]] <- HYs
    }
    objective <- .ssc_objective(Y, HY, lambda)
    if (objective < previous - 1e-10 * max(1, abs(previous))) monotone <- FALSE

    stationarity <- numeric(S)
    constraint <- numeric(S)
    for (s in seq_len(S)) {
      R <- if (full) off[1L] * (HYsum - HY[[s]]) else .ssc_gradient_input(s, HY, lambda)
      constraint[s] <- max(abs(crossprod(Q[[s]]) - diag(1, d)))
      if (is.null(R)) next
      G <- crossprod(Ut[[s]], R)
      QtG <- crossprod(Q[[s]], G)
      grad <- G - Q[[s]] %*% (0.5 * (QtG + t(QtG)))
      stationarity[s] <- sqrt(sum(grad^2)) / max(sqrt(sum(G^2)), 1e-300)
    }
    rel_change <- abs(objective - previous) / max(1, abs(objective))
    history[[it]] <- data.frame(
      iteration = it,
      objective = objective,
      relative_change = rel_change,
      max_stationarity = max(stationarity),
      max_constraint_error = max(constraint),
      min_supported_rank = min(supported)
    )
    if (isTRUE(control$verbose)) {
      message(sprintf(
        "stimulus_sumcor: sweep %d objective=%.10g rel_change=%.3g stationarity=%.3g",
        it, objective, rel_change, max(stationarity)
      ))
    }
    if (rel_change <= control$tol && max(stationarity) <= control$grad_tol) {
      converged <- TRUE
      status <- "converged"
      break
    }
  }
  history <- do.call(rbind, history)

  # Common rotation ordering components by cross-subject agreement. It leaves
  # the objective and every constraint unchanged.
  M <- matrix(0, d, d)
  for (s in seq_len(S - 1L)) for (t in seq.int(s + 1L, S)) {
    if (lambda[s, t] > 0) M <- M + lambda[s, t] * crossprod(Y[[s]], HY[[t]])
  }
  eg <- eigen(0.5 * (M + t(M)), symmetric = TRUE)
  rot <- eg$vectors
  Ysum <- Reduce(`+`, Y) %*% rot
  pivot <- apply(abs(Ysum), 2L, which.max)
  flip <- sign(Ysum[cbind(pivot, seq_len(d))])
  flip[flip == 0] <- 1
  rot <- sweep(rot, 2L, flip, "*")
  Q <- lapply(Q, function(q) q %*% rot)
  Y <- lapply(Y, function(y) y %*% rot)
  n_pairs <- sum(lambda[upper.tri(lambda)] > 0)

  list(
    Q = Q,
    Y = Y,
    objective = objective,
    history = history,
    converged = converged,
    status = status,
    iterations = it,
    monotone = monotone,
    supported_rank = stats::setNames(supported, domain_names),
    component_agreement = eg$values / max(n_pairs, 1L),
    n_pairs = n_pairs
  )
}

# -- Per-subject preparation ------------------------------------------------------

.ssc_prepare_domain <- function(dom, events, w, control, parcels, domain) {
  obs <- dom$observations
  used <- which(obs$used)
  ev_idx <- match(obs$event_id[used], events$event_id)
  counts <- tabulate(ev_idx, nbins = nrow(events))
  mask <- counts > 0
  m_s <- sum(w[mask])
  if (m_s <= 0) {
    stop("Domain `", domain, "`: its observed events carry zero total weight.",
         call. = FALSE)
  }
  omega <- w[ev_idx] / (counts[ev_idx] * m_s)
  basis <- .ssc_fit_basis(dom$X[used, , drop = FALSE], omega, control,
                          parcels, domain)
  U <- basis$U
  p <- ncol(U)
  # Event means through a sparse averaging operator (events x used rows).
  avg <- Matrix::sparseMatrix(i = ev_idx, j = seq_along(ev_idx),
                              x = 1 / counts[ev_idx],
                              dims = c(nrow(events), length(ev_idx)))
  Ubar <- as.matrix(avg %*% U)
  Chat <- crossprod(U * sqrt(omega))
  rho <- control$covariance_shrinkage
  C <- (1 - rho) * Chat + rho * sum(diag(Chat)) / p * diag(1, p)
  C <- 0.5 * (C + t(C))
  eg <- eigen(C, symmetric = TRUE)
  lam <- eg$values
  if (!length(lam) || min(lam) <= 1e-12 * max(lam, 0) || max(lam) <= 0) {
    stop("Domain `", domain, "`: the presentation covariance is singular; ",
         "increase `covariance_shrinkage` or reduce the brain rank.",
         call. = FALSE)
  }
  C_inv_half <- eg$vectors %*% (t(eg$vectors) / sqrt(lam))
  list(
    basis = basis,
    Ubar = Ubar,
    Ut = (Ubar * sqrt(w)) %*% C_inv_half,
    C = C,
    C_eigenvalues = lam,
    C_inv_half = C_inv_half,
    mask = mask,
    counts = counts,
    weight_mass = m_s,
    p = p
  )
}

.ssc_encoder <- function(prep, A, feature_names, control) {
  P <- prep$basis$P
  W <- as.matrix(P %*% A)
  if (!is.null(prep$basis$scale)) W <- W / prep$basis$scale
  list(
    center = prep$basis$center,
    weights = W,
    feature_names = feature_names,
    n_features = length(prep$basis$center),
    A = A,
    C = prep$C,
    C_eigenvalues = prep$C_eigenvalues,
    condition_number = max(prep$C_eigenvalues) / min(prep$C_eigenvalues),
    basis = list(
      type = prep$basis$info$type,
      rank = prep$p,
      singular_values = prep$basis$info$singular_values,
      scale = prep$basis$scale,
      P = if (isTRUE(control$store_basis)) P else NULL,
      details = prep$basis$info[setdiff(names(prep$basis$info),
                                        c("type", "singular_values"))]
    )
  )
}

.ssc_encode <- function(encoder, X) {
  Xc <- sweep(X, 2L, encoder$center, "-")
  Z <- Xc %*% encoder$weights
  colnames(Z) <- paste0("ssc_", seq_len(ncol(Z)))
  Z
}

.ssc_event_scores <- function(prep, A, event_ids) {
  Z <- prep$Ubar %*% A
  Z[!prep$mask, ] <- NA_real_
  dimnames(Z) <- list(event_ids, paste0("ssc_", seq_len(ncol(Z))))
  Z
}

# -- Fit ----------------------------------------------------------------------------

.ssc_resolve_events_table <- function(events) {
  if (is.null(events)) return(NULL)
  events <- as.data.frame(events, stringsAsFactors = FALSE)
  if (!"event_id" %in% names(events)) {
    stop("`events` must have an `event_id` column.", call. = FALSE)
  }
  if (anyNA(events$event_id) || anyDuplicated(as.character(events$event_id))) {
    stop("`events$event_id` must be unique and non-missing.", call. = FALSE)
  }
  events
}

.ssc_fit <- function(X_list, design_list, domain_names, ids, events,
                     stimulus_features, parcels, ncomp, control, call) {
  control <- .ssc_resolve_control(control)
  ncomp <- .ssc_int(ncomp, "ncomp")
  S <- length(X_list)
  if (S < 2L) stop("stimulus_sumcor_align() needs at least two domains.",
                   call. = FALSE)
  if (anyDuplicated(domain_names)) {
    stop("Domain names must be unique.", call. = FALSE)
  }
  events_in <- .ssc_resolve_events_table(events)
  if (!is.null(parcels)) {
    if (!is.list(parcels) || length(parcels) != S) {
      stop("`parcels` must be a list with one entry per domain.", call. = FALSE)
    }
    if (!is.null(names(parcels)) && all(domain_names %in% names(parcels))) {
      parcels <- parcels[domain_names]
    }
  }

  domains <- lapply(seq_len(S), function(i) {
    .ssc_resolve_domain(X_list[[i]], design_list[[i]], ids, domain_names[[i]])
  })
  names(domains) <- domain_names

  observed <- unlist(lapply(domains, function(d) d$observations$event_id[d$observations$used]),
                     use.names = FALSE)
  universe <- sort(unique(observed), method = "radix")
  if (length(universe) < 2L) {
    stop("At least two distinct events are needed.", call. = FALSE)
  }
  ev <- data.frame(event_id = universe, stringsAsFactors = FALSE)
  ev$modality <- .ssc_event_attribute(domains, events_in, "modality", universe)
  ev$clip <- .ssc_event_attribute(domains, events_in, "clip", universe)
  w <- .ssc_event_weights(control$event_weights, ev)
  ev$weight <- w

  preps <- lapply(seq_len(S), function(i) {
    .ssc_prepare_domain(domains[[i]], ev, w, control,
                        if (is.null(parcels)) NULL else parcels[[i]],
                        domain_names[[i]])
  })
  names(preps) <- domain_names
  ev$n_domains <- Reduce(`+`, lapply(preps, function(p) as.integer(p$mask)))

  ranks <- vapply(preps, `[[`, integer(1), "p")
  if (ncomp > min(ranks)) {
    stop("`ncomp = ", ncomp, "` exceeds the brain rank of domain `",
         names(ranks)[which.min(ranks)], "` (", min(ranks), ").", call. = FALSE)
  }

  pairs <- .ssc_pair_weights(lapply(preps, `[[`, "mask"), w, control$min_overlap)
  dimnames(pairs$overlap) <- dimnames(pairs$lambda) <- list(domain_names, domain_names)
  isolated <- domain_names[rowSums(pairs$lambda > 0) == 0]
  if (length(isolated)) {
    stop("Domain(s) ", paste0("`", isolated, "`", collapse = ", "),
         " share less than `min_overlap` event weight with every other domain.",
         call. = FALSE)
  }

  kernel <- .ssc_kernel_factor(stimulus_features, ev, control)
  solved <- .ssc_solve(lapply(preps, `[[`, "Ut"), kernel$L, pairs$lambda,
                       ncomp, control, domain_names)

  encoders <- list()
  scores <- list()
  event_scores <- list()
  for (i in seq_len(S)) {
    nm <- domain_names[[i]]
    A <- preps[[i]]$C_inv_half %*% solved$Q[[i]]
    enc <- .ssc_encoder(preps[[i]], A, colnames(domains[[i]]$X), control)
    Z <- .ssc_encode(enc, domains[[i]]$X)
    rownames(Z) <- domains[[i]]$observations$observation_id
    encoders[[nm]] <- enc
    scores[[nm]] <- Z
    event_scores[[nm]] <- .ssc_event_scores(preps[[i]], A, universe)
  }
  stack <- simplify2array(event_scores)
  reference <- apply(stack, c(1L, 2L), function(v) {
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  })
  dimnames(reference) <- dimnames(event_scores[[1L]])

  if (!solved$converged) {
    warning(
      "stimulus_sumcor_align() did not converge in ", control$max_iter,
      " sweeps (relative change ",
      signif(utils::tail(solved$history$relative_change, 1L), 3),
      ", stationarity ",
      signif(utils::tail(solved$history$max_stationarity, 1L), 3), ").",
      call. = FALSE
    )
  }
  if (min(solved$supported_rank) < ncomp) {
    warning(
      "Only ", min(solved$supported_rank), " of ", ncomp, " latent directions ",
      "are supported by cross-subject agreement in at least one domain; the ",
      "remaining directions are a deterministic completion and are not learned ",
      "correspondence.", call. = FALSE
    )
  }

  structure(
    list(
      call = call,
      domain_names = domain_names,
      domain_roles = stats::setNames(rep("reference", S), domain_names),
      ncomp = ncomp,
      scores = scores,
      event_scores = event_scores,
      reference_scores = reference,
      events = ev,
      observations = lapply(domains, `[[`, "observations"),
      encoders = encoders,
      overlap = pairs$overlap,
      pair_weights = pairs$lambda,
      kernel = list(
        beta = control$beta,
        factor = kernel$L,
        families = kernel$families,
        family_weights = kernel$theta,
        dropped_feature_rows = kernel$dropped_rows
      ),
      convergence = list(
        converged = solved$converged,
        status = solved$status,
        iterations = solved$iterations,
        objective = solved$objective,
        monotone = solved$monotone,
        history = solved$history
      ),
      supported_rank = solved$supported_rank,
      component_agreement = solved$component_agreement,
      calibration = list(),
      ids = ids,
      control = control,
      provenance = list(
        package_version = as.character(utils::packageVersion("manifoldalign")),
        r_version = R.version.string,
        n_domains = S,
        n_events = length(universe),
        observations_used = vapply(domains, function(d) sum(d$observations$used),
                                   integer(1)),
        observations_unlabelled = vapply(domains,
                                         function(d) sum(!d$observations$used),
                                         integer(1))
      )
    ),
    class = c("stimulus_sumcor", "manifold_alignment")
  )
}

#' Stimulus-weighted SUMCOR alignment of repeated naturalistic data
#'
#' Learns one linear encoder per subject (domain) that maps brain observations
#' into a shared latent space. The encoders maximise reproducible
#' cross-subject agreement of event-locked activity, where an event is a
#' stimulus moment (for example a time point within a music excerpt or video)
#' identified by a keyed `event_id`. Repeated presentations of an event are
#' averaged in the agreement term but enter the per-subject normalisation
#' individually, so directions dominated by presentation-to-presentation
#' variability are down-weighted. Optional stimulus features add a metric
#' over events, so that events with similar descriptors also contribute to
#' agreement.
#'
#' @section Model:
#' For subject `s`, a fitted linear basis gives `u = phi_s(x)`. With event
#' weights `w` (summing to one), event means `ubar_s`, presentation covariance
#' `C_s` (shrunk toward a scaled identity), and stimulus metric
#' `H = I + beta * sum_m theta_m K_m` with
#' `K_m = F_m (F_m' F_m + eta_m I)^{-1} F_m'`, the encoders solve
#' \deqn{\max \sum_{s<t} \mathrm{tr}(A_s^\top B_{st} A_t)
#'   \quad \mathrm{s.t.} \quad A_s^\top C_s A_s = I,}
#' with `B_st = ubar_s' W^{1/2} H W^{1/2} ubar_t / m_st`, where `m_st` is the
#' event weight both subjects observed. Unobserved events are masked by key.
#' The problem is solved by sequential polar updates, which never decrease the
#' objective; joint global optimality is not guaranteed for more than two
#' subjects. The encoder `f_s(x) = phi_s(x) A_s` is stored explicitly, and the
#' training scores are produced by that same function.
#'
#' @section Identities:
#' Correspondence between subjects is carried only by event identifiers. Row
#' order and row counts are never used to infer it. Observation, repeat, run,
#' modality and clip identifiers are kept in the result. Rows whose event is
#' `NA` (for example censored or non-stimulus volumes) are excluded from
#' fitting but still receive scores.
#'
#' @section Features:
#' Each element of `stimulus_features` is a matrix with one row per event and
#' row names equal to event identifiers. Rows for events no subject observed
#' are dropped. Apply any hemodynamic convolution or lag before matching
#' features to fMRI observations. With `cross_modal = FALSE` (default) each
#' feature kernel is block diagonal by event modality. A feature kernel close
#' to a multiple of the identity changes nothing and triggers a warning.
#'
#' @param data A list of numeric observation-by-feature matrices (one per
#'   subject) or a `hyperdesign`.
#' @param ... Method arguments. Unknown arguments are rejected.
#' @param design For list input, a list of data frames (one per subject, rows
#'   matching `data`) holding the identifier columns.
#' @param event_id Name of the design column with event identifiers.
#' @param repeat_id Optional design column with repeat identifiers. Without it,
#'   repeats are numbered by occurrence within each event.
#' @param run_id Optional design column with run identifiers (stored only).
#' @param observation_id Optional design column with unique observation
#'   identifiers. Defaults to row names or row positions.
#' @param modality Optional design column giving each observation's event
#'   modality.
#' @param clip Optional design column giving each observation's clip.
#' @param events Optional data frame with an `event_id` column and optional
#'   `modality` and `clip` columns describing events.
#' @param stimulus_features Optional named list of event-by-feature matrices
#'   (row names are event identifiers).
#' @param parcels For `brain_basis = "parcel_pca"`, a list with one parcel
#'   label vector per subject (one entry per brain feature; `NA` excludes a
#'   feature).
#' @param ncomp Latent dimension.
#' @param control Settings from [stimulus_sumcor_control()] or a named list of
#'   its arguments.
#' @return A `stimulus_sumcor` object with per-subject `scores` (all
#'   observations, original order, rows named by observation identifier),
#'   `event_scores`, frozen `reference_scores`, `encoders`, identity tables,
#'   kernel diagnostics, `convergence` (objective history, stationarity,
#'   constraint error, supported rank), and `component_agreement`.
#' @seealso [calibrate_domain()] to add a new subject against the frozen
#'   reference, [oos_predict()] to encode new observations.
#' @examples
#' set.seed(1)
#' n_events <- 40
#' latent <- matrix(rnorm(n_events * 3), n_events, 3)
#' make_subject <- function(p) {
#'   mix <- matrix(rnorm(3 * p), 3, p)
#'   rows <- rep(seq_len(n_events), times = 2)
#'   X <- latent[rows, ] %*% mix + matrix(rnorm(length(rows) * p, sd = 0.5),
#'                                        length(rows), p)
#'   list(X = X, design = data.frame(event_id = paste0("e", rows),
#'                                   repeat_id = rep(1:2, each = n_events)))
#' }
#' subjects <- list(a = make_subject(20), b = make_subject(25),
#'                  c = make_subject(18))
#' fit <- stimulus_sumcor_align(
#'   lapply(subjects, `[[`, "X"),
#'   design = lapply(subjects, `[[`, "design"),
#'   event_id = "event_id", repeat_id = "repeat_id",
#'   ncomp = 3,
#'   control = stimulus_sumcor_control(beta = 0, brain_rank = 10)
#' )
#' fit$component_agreement
#' @export
stimulus_sumcor_align <- function(data, ...) {
  UseMethod("stimulus_sumcor_align")
}

.ssc_reject_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    nm[!nzchar(nm)] <- "<unnamed>"
    stop("Unknown argument(s): ", paste(nm, collapse = ", "), ".",
         call. = FALSE)
  }
}

#' @rdname stimulus_sumcor_align
#' @export
stimulus_sumcor_align.list <- function(
  data,
  design,
  event_id = "event_id",
  repeat_id = NULL,
  run_id = NULL,
  observation_id = NULL,
  modality = NULL,
  clip = NULL,
  events = NULL,
  stimulus_features = NULL,
  parcels = NULL,
  ncomp = 10L,
  control = stimulus_sumcor_control(),
  ...
) {
  .ssc_reject_dots(...)
  domain_names <- names(data)
  if (is.null(domain_names) || anyNA(domain_names) || any(!nzchar(domain_names))) {
    domain_names <- paste0("domain", seq_along(data))
  }
  if (missing(design) || !is.list(design) || is.data.frame(design) ||
      length(design) != length(data)) {
    stop("`design` must be a list with one data frame per domain.",
         call. = FALSE)
  }
  if (!is.null(names(design)) && all(domain_names %in% names(design))) {
    design <- design[domain_names]
  }
  ids <- list(event_id = event_id, repeat_id = repeat_id, run_id = run_id,
              observation_id = observation_id, modality = modality, clip = clip)
  .ssc_fit(unname(data), unname(design), domain_names, ids, events,
           stimulus_features, parcels, ncomp, control, match.call())
}

#' @rdname stimulus_sumcor_align
#' @export
stimulus_sumcor_align.hyperdesign <- function(
  data,
  event_id = "event_id",
  repeat_id = NULL,
  run_id = NULL,
  observation_id = NULL,
  modality = NULL,
  clip = NULL,
  events = NULL,
  stimulus_features = NULL,
  parcels = NULL,
  ncomp = 10L,
  control = stimulus_sumcor_control(),
  ...
) {
  .ssc_reject_dots(...)
  resolved <- resolve_hyperdesign(data)
  X_list <- lapply(resolved$domains, `[[`, "x")
  design_list <- lapply(resolved$domains, `[[`, "design")
  ids <- list(event_id = event_id, repeat_id = repeat_id, run_id = run_id,
              observation_id = observation_id, modality = modality, clip = clip)
  .ssc_fit(unname(X_list), unname(design_list), resolved$domain_names, ids,
           events, stimulus_features, parcels, ncomp, control, match.call())
}

#' @export
stimulus_sumcor_align.default <- function(data, ...) {
  stop("stimulus_sumcor_align() needs a list of matrices or a hyperdesign.",
       call. = FALSE)
}

#' @method print stimulus_sumcor
#' @export
print.stimulus_sumcor <- function(x, ...) {
  cat("Stimulus-weighted SUMCOR alignment\n")
  roles <- x$domain_roles
  cat("  domains: ", length(x$domain_names), " (",
      sum(roles == "reference"), " reference, ",
      sum(roles == "calibrated"), " calibrated)\n", sep = "")
  cat("  events: ", nrow(x$events), "\n", sep = "")
  cat("  components: ", x$ncomp, "\n", sep = "")
  fam <- if (is.null(x$kernel$families)) "none" else
    paste(unique(x$kernel$families$family), collapse = ", ")
  cat("  beta: ", format(x$kernel$beta), " (features: ", fam, ")\n", sep = "")
  cat("  status: ", x$convergence$status, " after ", x$convergence$iterations,
      " sweep(s); objective ", format(x$convergence$objective, digits = 7),
      "\n", sep = "")
  cat("  min supported rank: ", min(x$supported_rank), "\n", sep = "")
  agree <- x$component_agreement
  cat("  component agreement: ",
      paste(format(utils::head(agree, 6L), digits = 3), collapse = " "),
      if (length(agree) > 6L) " ..." else "", "\n", sep = "")
  invisible(x)
}
