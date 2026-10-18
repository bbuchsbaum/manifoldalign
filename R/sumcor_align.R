# Event-keyed SUMCOR alignment ---------------------------------------------
#
# Mathematical contract (all quantities are per domain s unless stated):
#   u_ser   = phi_s(x_ser)                      fitted linear basis
#   ubar_se = mean_r u_ser                      event means (numerator)
#   C_s     = (1 - rho) Chat_s + rho tr(Chat_s)/p I
#   Chat_s  = (1/m_s) sum_e (w_e / R_se) sum_r u_ser u_ser'   (repeats)
#   H       = I + beta sum_m theta_m K_m,  K_m = F_m (F_m'F_m + eta_m I)^-1 F_m'
#   B_st    = ubar_s' W^1/2 H W^1/2 ubar_t / m_st              (overlap-normalised)
#   maximise sum_{s<t} tr(A_s' B_st A_t)  subject to  A_s' C_s A_s = I
#   encoder f_s(x) = ((x - mu_s) / sd_s) P_s A_s
# Event weights w sum to one over the event universe; m_s is the weight of
# the events domain s observed and m_st the weight of the shared events.

# -- Validation helpers ------------------------------------------------------

.sumcor_number <- function(x, name, lower = -Inf, upper = Inf,
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

.sumcor_int <- function(x, name, lower = 1L) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      abs(x - round(x)) > 1e-8 || x < lower) {
    stop("`", name, "` must be a single integer >= ", lower, ".", call. = FALSE)
  }
  as.integer(round(x))
}

.sumcor_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.sumcor_name <- function(x, name, allow_null = TRUE) {
  if (is.null(x) && allow_null) return(NULL)
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be one non-empty string.", call. = FALSE)
  }
  x
}

.sumcor_names <- function(x, name) {
  if (is.null(x)) return(character(0))
  if (!is.character(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) {
    stop("`", name, "` must be a character vector of unique column names.",
         call. = FALSE)
  }
  x
}

.sumcor_with_preserved_rng <- function(expr) {
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

#' Control settings for event-keyed SUMCOR alignment
#'
#' @param beta Nonnegative weight of the event-feature metric relative to exact
#'   event correspondence. `beta = 0` gives the correspondence-only model
#'   exactly; supplied event features are then not used.
#' @param kernel_weights Optional named nonnegative weights over event-feature
#'   families (`theta_m`), normalised to sum to one over the families
#'   supplied. `NULL` gives equal weights.
#' @param kernel_ridge Positive ridge for each event-feature kernel, relative
#'   to the mean squared singular value of the preprocessed feature block:
#'   `eta = kernel_ridge * mean(sigma^2)`.
#' @param kernel_rank Optional maximum number of singular components kept per
#'   kernel block. `NULL` keeps every component with nonnegligible weight.
#' @param kernel_scale Logical; scale event-feature columns to unit variance
#'   after centring (constant columns are dropped).
#' @param kernel_normalization `"max_eigenvalue"` rescales each family kernel so
#'   its largest eigenvalue is one; `"none"` keeps `K_m` as defined.
#' @param kernel_blocks Optional name of an event attribute. When given, each
#'   kernel is block diagonal by that attribute, so event features never
#'   relate events with different attribute values. `NULL` uses one block.
#' @param covariance_shrinkage Shrinkage `rho` in `[0, 1]` of each domain's
#'   repeat-level covariance toward a scaled identity.
#' @param basis Linear per-domain representation: `"pca"` (weighted PCA over
#'   all columns), `"grouped_pca"` (PCA within supplied column groups), or
#'   `"identity"` (centring and optional scaling only).
#' @param basis_rank Maximum number of components per domain for `"pca"`.
#' @param group_rank Maximum number of components per column group for
#'   `"grouped_pca"`.
#' @param basis_scale Logical; scale columns to unit (event-weighted) variance
#'   before the basis is fitted.
#' @param event_weights Event weighting: `"uniform"`, the name of an event
#'   attribute (each attribute value then receives equal total weight, split
#'   uniformly among its events), or a named nonnegative numeric vector over
#'   event identifiers.
#' @param min_overlap Minimum shared event weight (a fraction of the total) for
#'   a pair of domains to contribute to the objective.
#' @param max_iter Maximum number of sequential polar sweeps.
#' @param tol Relative objective-change tolerance for convergence.
#' @param grad_tol Relative projected-gradient (stationarity) tolerance for
#'   convergence.
#' @param rank_tol Relative singular-value threshold below which a polar update
#'   direction is treated as unsupported and completed deterministically.
#' @param identity_max_dim Largest column count allowed with
#'   `basis = "identity"`.
#' @param store_basis Logical; keep each domain's fitted basis matrix in the
#'   result. Prediction does not need it.
#' @param verbose Logical; print sweep progress. Logging never changes the
#'   computation.
#' @return An object of class `sumcor_control`.
#' @export
sumcor_control <- function(
  beta = 0.5,
  kernel_weights = NULL,
  kernel_ridge = 1,
  kernel_rank = NULL,
  kernel_scale = TRUE,
  kernel_normalization = c("max_eigenvalue", "none"),
  kernel_blocks = NULL,
  covariance_shrinkage = 0.1,
  basis = c("pca", "grouped_pca", "identity"),
  basis_rank = 50L,
  group_rank = 5L,
  basis_scale = FALSE,
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
  basis <- match.arg(basis)
  beta <- .sumcor_number(beta, "beta", lower = 0)
  kernel_ridge <- .sumcor_number(kernel_ridge, "kernel_ridge", lower = 0,
                                 strict_lower = TRUE)
  if (!is.null(kernel_rank)) kernel_rank <- .sumcor_int(kernel_rank, "kernel_rank")
  kernel_scale <- .sumcor_flag(kernel_scale, "kernel_scale")
  kernel_blocks <- .sumcor_name(kernel_blocks, "kernel_blocks")
  covariance_shrinkage <- .sumcor_number(covariance_shrinkage,
                                         "covariance_shrinkage", 0, 1)
  basis_rank <- .sumcor_int(basis_rank, "basis_rank")
  group_rank <- .sumcor_int(group_rank, "group_rank")
  basis_scale <- .sumcor_flag(basis_scale, "basis_scale")
  min_overlap <- .sumcor_number(min_overlap, "min_overlap", 0, 1)
  max_iter <- .sumcor_int(max_iter, "max_iter")
  tol <- .sumcor_number(tol, "tol", lower = 0)
  grad_tol <- .sumcor_number(grad_tol, "grad_tol", lower = 0)
  rank_tol <- .sumcor_number(rank_tol, "rank_tol", 0, 1, strict_lower = TRUE)
  identity_max_dim <- .sumcor_int(identity_max_dim, "identity_max_dim")
  store_basis <- .sumcor_flag(store_basis, "store_basis")
  verbose <- .sumcor_flag(verbose, "verbose")

  if (!is.null(kernel_weights)) {
    if (!is.numeric(kernel_weights) || !length(kernel_weights) ||
        any(!is.finite(kernel_weights)) || any(kernel_weights < 0) ||
        is.null(names(kernel_weights)) || any(!nzchar(names(kernel_weights))) ||
        anyDuplicated(names(kernel_weights))) {
      stop("`kernel_weights` must be a named, finite, nonnegative numeric ",
           "vector with unique names.", call. = FALSE)
    }
    if (sum(kernel_weights) <= 0) {
      stop("`kernel_weights` must not all be zero.", call. = FALSE)
    }
  }

  if (is.character(event_weights)) {
    .sumcor_name(event_weights, "event_weights", allow_null = FALSE)
  } else if (is.numeric(event_weights)) {
    if (!length(event_weights) || any(!is.finite(event_weights)) ||
        any(event_weights < 0) || is.null(names(event_weights)) ||
        anyNA(names(event_weights)) || any(!nzchar(names(event_weights))) ||
        anyDuplicated(names(event_weights))) {
      stop("Numeric `event_weights` must be finite, nonnegative, and uniquely ",
           "named by event identifier.", call. = FALSE)
    }
  } else {
    stop("`event_weights` must be \"uniform\", an event attribute name, or a ",
         "named numeric vector.", call. = FALSE)
  }

  structure(
    list(
      beta = beta,
      kernel_weights = kernel_weights,
      kernel_ridge = kernel_ridge,
      kernel_rank = kernel_rank,
      kernel_scale = kernel_scale,
      kernel_normalization = kernel_normalization,
      kernel_blocks = kernel_blocks,
      covariance_shrinkage = covariance_shrinkage,
      basis = basis,
      basis_rank = basis_rank,
      group_rank = group_rank,
      basis_scale = basis_scale,
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
    class = "sumcor_control"
  )
}

# Every path re-runs the constructor, so list input is always validated and
# unknown fields are always rejected.
.sumcor_resolve_control <- function(control) {
  if (is.null(control)) return(sumcor_control())
  if (!is.list(control)) {
    stop("`control` must be NULL, a named list, or sumcor_control().",
         call. = FALSE)
  }
  values <- unclass(control)
  if (length(values) &&
      (is.null(names(values)) || any(!nzchar(names(values))))) {
    stop("`control` must be a named list.", call. = FALSE)
  }
  unknown <- setdiff(names(values), names(formals(sumcor_control)))
  if (length(unknown)) {
    stop("Unknown sumcor control field(s): ", paste(unknown, collapse = ", "),
         ".", call. = FALSE)
  }
  do.call(sumcor_control, values)
}

# -- Input resolution ----------------------------------------------------------

.sumcor_as_numeric_matrix <- function(X, name) {
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

.sumcor_design_column <- function(design, column, role, domain) {
  if (is.null(column)) return(NULL)
  if (!column %in% names(design)) {
    stop("Domain `", domain, "`: design column `", column, "` (", role,
         ") is missing.", call. = FALSE)
  }
  design[[column]]
}

# Resolve one domain's observation table. Correspondence is carried only by
# keyed identifiers; row counts and row order are never used to infer it.
.sumcor_resolve_domain <- function(X, design, ids, domain) {
  X <- .sumcor_as_numeric_matrix(X, paste0("data[[\"", domain, "\"]]"))
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

  event <- as.character(.sumcor_design_column(design, ids$event_id, "event_id",
                                              domain))
  if (any(!is.na(event) & !nzchar(event))) {
    stop("Domain `", domain, "`: event identifiers must not be empty strings.",
         call. = FALSE)
  }

  obs <- .sumcor_design_column(design, ids$observation_id, "observation_id",
                               domain)
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

  rep_raw <- .sumcor_design_column(design, ids$repeat_id, "repeat_id", domain)
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

  table <- data.frame(
    row = seq_len(n),
    observation_id = obs,
    event_id = event,
    repeat_id = repeat_id,
    used = !is.na(event),
    stringsAsFactors = FALSE
  )
  for (col in ids$observation_attributes) {
    table[[col]] <- .sumcor_design_column(design, col, "observation attribute",
                                          domain)
  }
  attributes <- list()
  for (col in ids$event_attributes) {
    attributes[[col]] <- as.character(
      .sumcor_design_column(design, col, "event attribute", domain)
    )
  }
  list(X = X, observations = table, event_attributes = attributes)
}

# Collapse an event attribute from observation-level design columns and the
# optional `events` table, requiring one consistent value per event.
.sumcor_event_attribute <- function(domains, events_df, field, universe) {
  value <- rep(NA_character_, length(universe))
  assign_value <- function(ev, val, source) {
    keep <- !is.na(ev) & !is.na(val)
    if (!any(keep)) return(invisible(NULL))
    pairs <- unique(data.frame(ev = ev[keep], val = as.character(val[keep]),
                               stringsAsFactors = FALSE))
    if (anyDuplicated(pairs$ev)) {
      bad <- pairs$ev[duplicated(pairs$ev)][1L]
      stop("Event `", bad, "` has more than one `", field, "` value in ",
           source, ".", call. = FALSE)
    }
    idx <- match(pairs$ev, universe)
    pairs <- pairs[!is.na(idx), , drop = FALSE]
    idx <- idx[!is.na(idx)]
    clash <- !is.na(value[idx]) & value[idx] != pairs$val
    if (any(clash)) {
      stop("Event `", pairs$ev[clash][1L], "` has inconsistent `", field,
           "` values across inputs.", call. = FALSE)
    }
    value[idx] <<- pairs$val
    invisible(NULL)
  }
  for (dom in names(domains)) {
    observed <- domains[[dom]]$event_attributes[[field]]
    if (!is.null(observed)) {
      assign_value(domains[[dom]]$observations$event_id, observed,
                   paste0("domain `", dom, "`"))
    }
  }
  if (!is.null(events_df) && field %in% names(events_df)) {
    assign_value(as.character(events_df$event_id), events_df[[field]],
                 "`events`")
  }
  value
}

.sumcor_event_weights <- function(spec, events) {
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
    if (!spec %in% names(events) || identical(spec, "event_id")) {
      stop("`event_weights = \"", spec, "\"` names no event attribute.",
           call. = FALSE)
    }
    group <- events[[spec]]
    if (anyNA(group)) {
      stop("`event_weights = \"", spec, "\"` requires that attribute for ",
           "every event.", call. = FALSE)
    }
    counts <- table(group)
    w <- 1 / (length(counts) * as.numeric(counts[group]))
  }
  if (sum(w) <= 0) stop("Event weights sum to zero.", call. = FALSE)
  w / sum(w)
}

# -- Per-domain basis ------------------------------------------------------------

.sumcor_svd <- function(M, k) {
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
    d <- sqrt(pmax(e$values[seq_len(k)], 0))
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
    s <- .sumcor_with_preserved_rng(
      irlba::irlba(M, nv = k, nu = k, v = start, tol = 1e-8, maxit = 2000L)
    )
    out <- list(d = s$d, u = s$u, v = s$v)
  }
  # Deterministic signs: largest-magnitude loading of each component positive.
  pivot <- apply(abs(out$v), 2L, which.max)
  flip <- sign(out$v[cbind(pivot, seq_len(k))])
  flip[flip == 0] <- 1
  out$v <- sweep(out$v, 2L, flip, "*")
  out$u <- sweep(out$u, 2L, flip, "*")
  out
}

.sumcor_fit_basis <- function(X, weights, control, groups, domain) {
  V <- ncol(X)
  center <- colSums(X * weights)
  Xc <- sweep(X, 2L, center, "-")
  scale <- NULL
  if (isTRUE(control$basis_scale)) {
    scale <- sqrt(colSums(Xc^2 * weights))
    scale[!is.finite(scale) | scale <= 0] <- 1
    Xc <- sweep(Xc, 2L, scale, "/")
  }
  Xw <- Xc * sqrt(weights)
  n_eff <- sum(weights > 0)

  if (identical(control$basis, "identity")) {
    if (V > control$identity_max_dim) {
      stop("Domain `", domain, "`: basis = \"identity\" with ", V,
           " columns exceeds `identity_max_dim = ", control$identity_max_dim,
           "`.", call. = FALSE)
    }
    P <- diag(1, V)
    info <- list(type = "identity",
                 singular_values = svd(Xw, nu = 0L, nv = 0L)$d)
  } else if (identical(control$basis, "pca")) {
    k <- min(control$basis_rank, V, max(n_eff - 1L, 1L))
    s <- .sumcor_svd(Xw, k)
    keep <- s$d > control$rank_tol * max(s$d, 0)
    P <- s$v[, keep, drop = FALSE]
    info <- list(type = "pca", singular_values = s$d[keep],
                 dropped_degenerate = sum(!keep))
  } else {
    if (is.null(groups)) {
      stop("Domain `", domain, "`: basis = \"grouped_pca\" requires ",
           "`feature_groups`.", call. = FALSE)
    }
    if (length(groups) != V) {
      stop("Domain `", domain, "`: `feature_groups` must have one entry per ",
           "column.", call. = FALSE)
    }
    labels <- as.character(groups)
    levels <- sort(unique(labels[!is.na(labels)]), method = "radix")
    if (!length(levels)) {
      stop("Domain `", domain, "`: `feature_groups` assigns no column to a ",
           "group.", call. = FALSE)
    }
    rows <- integer(0)
    cols <- integer(0)
    vals <- numeric(0)
    svals <- numeric(0)
    group_of <- character(0)
    offset <- 0L
    for (lev in levels) {
      idx <- which(labels == lev)
      k <- min(control$group_rank, length(idx), max(n_eff - 1L, 1L))
      s <- .sumcor_svd(Xw[, idx, drop = FALSE], k)
      keep <- s$d > control$rank_tol * max(s$d, 0)
      if (!any(keep)) next
      v <- s$v[, keep, drop = FALSE]
      kk <- ncol(v)
      rows <- c(rows, rep(idx, kk))
      cols <- c(cols, rep(offset + seq_len(kk), each = length(idx)))
      vals <- c(vals, as.numeric(v))
      svals <- c(svals, s$d[keep])
      group_of <- c(group_of, rep(lev, kk))
      offset <- offset + kk
    }
    if (!offset) {
      stop("Domain `", domain, "`: no column group has a nondegenerate ",
           "component.", call. = FALSE)
    }
    P <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(V, offset))
    info <- list(type = "grouped_pca", singular_values = svals,
                 group = group_of, n_groups = length(levels))
  }
  U <- as.matrix(Xc %*% P)
  list(center = center, scale = scale, P = P, U = U, info = info)
}

# -- Event-feature kernels ------------------------------------------------------

.sumcor_kernel_factor <- function(features, events, control) {
  universe <- events$event_id
  n <- length(universe)
  empty <- list(L = NULL, families = NULL, dropped_rows = integer(0),
                theta = NULL)
  if (is.null(features) || control$beta == 0) return(empty)
  if (!is.list(features) || is.data.frame(features) || !length(features)) {
    stop("`event_features` must be a named list of feature matrices.",
         call. = FALSE)
  }
  fam_names <- names(features)
  if (is.null(fam_names) || anyNA(fam_names) || any(!nzchar(fam_names)) ||
      anyDuplicated(fam_names)) {
    stop("`event_features` must have unique, non-empty family names.",
         call. = FALSE)
  }
  theta <- control$kernel_weights
  if (is.null(theta)) {
    theta <- stats::setNames(rep(1, length(fam_names)), fam_names)
  } else {
    if (!setequal(names(theta), fam_names)) {
      stop("`kernel_weights` names must match the event-feature families.",
           call. = FALSE)
    }
    theta <- theta[fam_names]
  }
  theta <- theta / sum(theta)

  blocks_by <- control$kernel_blocks
  if (!is.null(blocks_by) &&
      (!blocks_by %in% names(events) || identical(blocks_by, "event_id"))) {
    stop("`kernel_blocks = \"", blocks_by, "\"` names no event attribute.",
         call. = FALSE)
  }

  factors <- list()
  diagnostics <- list()
  dropped <- integer(0)
  for (fam in fam_names) {
    F <- features[[fam]]
    if (is.data.frame(F)) F <- as.matrix(F)
    if (inherits(F, "Matrix")) F <- as.matrix(F)
    if (!is.matrix(F) || !is.numeric(F)) {
      stop("Event-feature family `", fam, "` must be a numeric matrix.",
           call. = FALSE)
    }
    ids <- rownames(F)
    if (is.null(ids) || anyNA(ids) || anyDuplicated(ids)) {
      stop("Event-feature family `", fam, "` needs unique row names giving ",
           "event identifiers.", call. = FALSE)
    }
    if (any(!is.finite(F))) {
      stop("Event-feature family `", fam, "` contains non-finite values.",
           call. = FALSE)
    }
    keep_rows <- ids %in% universe
    dropped[[fam]] <- sum(!keep_rows)
    F <- F[keep_rows, , drop = FALSE]
    if (nrow(F) < 2L) {
      stop("Event-feature family `", fam, "` describes fewer than two fitted ",
           "events.", call. = FALSE)
    }
    idx <- match(rownames(F), universe)
    F <- sweep(F, 2L, colMeans(F), "-")
    sds <- sqrt(colSums(F^2) / max(nrow(F) - 1L, 1L))
    informative <- sds > 1e-12 * max(sds, 0)
    if (!any(informative)) {
      stop("Event-feature family `", fam, "` has no non-constant column.",
           call. = FALSE)
    }
    F <- F[, informative, drop = FALSE]
    if (isTRUE(control$kernel_scale)) {
      F <- sweep(F, 2L, sds[informative], "/")
    }

    block_of <- if (is.null(blocks_by)) {
      rep("all", length(idx))
    } else {
      b <- events[[blocks_by]][idx]
      if (anyNA(b)) {
        stop("Event-feature family `", fam, "`: every described event needs ",
             "a `", blocks_by, "` value.", call. = FALSE)
      }
      b
    }
    blocks <- list()
    for (b in sort(unique(block_of), method = "radix")) {
      rows_b <- which(block_of == b)
      if (length(rows_b) < 2L) next
      s <- svd(F[rows_b, , drop = FALSE])
      positive <- s$d > 1e-10 * max(s$d, 0)
      if (!any(positive)) next
      sig2 <- s$d[positive]^2
      eta <- control$kernel_ridge * mean(sig2)
      kappa <- sig2 / (sig2 + eta)
      U <- s$u[, positive, drop = FALSE]
      if (!is.null(control$kernel_rank) && length(kappa) > control$kernel_rank) {
        kappa <- kappa[seq_len(control$kernel_rank)]
        U <- U[, seq_len(control$kernel_rank), drop = FALSE]
      }
      nb <- length(rows_b)
      cbar <- sum(kappa) / nb
      isotropy <- sqrt((sum((kappa - cbar)^2) + (nb - length(kappa)) * cbar^2) /
                         sum(kappa^2))
      blocks[[b]] <- list(rows = idx[rows_b], U = U, kappa = kappa, eta = eta,
                          n_events = nb, isotropy_deviation = isotropy)
    }
    if (!length(blocks)) {
      stop("Event-feature family `", fam, "` yields no usable kernel block.",
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
          "Event-feature family `", fam, "` (block `", b, "`) gives a kernel ",
          "close to a multiple of the identity; it will change the fit very ",
          "little. Reduce its effective rank or increase `kernel_ridge`.",
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

.sumcor_apply_H <- function(L, Y) {
  if (is.null(L)) return(Y)
  Y + L %*% crossprod(L, Y)
}

# -- Polar update with deterministic completion ---------------------------------

.sumcor_polar <- function(G, Q_prev = NULL, rank_tol = 1e-8) {
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

.sumcor_pair_weights <- function(masks, w, min_overlap) {
  S <- length(masks)
  overlap <- matrix(0, S, S)
  for (s in seq_len(S)) for (t in seq_len(S)) {
    overlap[s, t] <- sum(w[masks[[s]] & masks[[t]]])
  }
  lambda <- ifelse(overlap >= min_overlap & overlap > 0, 1 / overlap, 0)
  diag(lambda) <- 0
  list(overlap = overlap, lambda = lambda)
}

.sumcor_objective <- function(Y, HY, lambda) {
  S <- length(Y)
  total <- 0
  for (s in seq_len(S - 1L)) for (t in seq.int(s + 1L, S)) {
    if (lambda[s, t] > 0) total <- total + lambda[s, t] * sum(Y[[s]] * HY[[t]])
  }
  total
}

.sumcor_gradient_input <- function(s, HY, lambda) {
  acc <- NULL
  for (t in seq_along(HY)) {
    if (t == s || lambda[s, t] == 0) next
    term <- lambda[s, t] * HY[[t]]
    acc <- if (is.null(acc)) term else acc + term
  }
  acc
}

.sumcor_solve <- function(Ut, L, lambda, d, control, domain_names) {
  S <- length(Ut)
  n_ev <- nrow(Ut[[1L]])
  off <- lambda[upper.tri(lambda)]
  full <- all(off > 0) && (max(off) - min(off)) <= 1e-12 * max(off)

  # Order-invariant MAXVAR-style start: leading left singular vectors of the
  # concatenated whitened event scores (column order does not change them).
  T0 <- .sumcor_svd(do.call(cbind, Ut), d)$u
  if (ncol(T0) < d) T0 <- cbind(T0, matrix(0, n_ev, d - ncol(T0)))
  HT0 <- .sumcor_apply_H(L, T0)
  Q <- vector("list", S)
  supported <- integer(S)
  for (s in seq_len(S)) {
    pol <- .sumcor_polar(crossprod(Ut[[s]], HT0), NULL, control$rank_tol)
    Q[[s]] <- pol$Q
    supported[s] <- pol$rank
  }
  Y <- Map(function(U, Qs) U %*% Qs, Ut, Q)
  HY <- lapply(Y, function(y) .sumcor_apply_H(L, y))
  objective <- .sumcor_objective(Y, HY, lambda)

  history <- list()
  converged <- FALSE
  monotone <- TRUE
  status <- "max_iter"
  it <- 0L
  for (it in seq_len(control$max_iter)) {
    previous <- objective
    if (full) HYsum <- Reduce(`+`, HY)
    for (s in seq_len(S)) {
      R <- if (full) off[1L] * (HYsum - HY[[s]]) else
        .sumcor_gradient_input(s, HY, lambda)
      if (is.null(R)) next
      pol <- .sumcor_polar(crossprod(Ut[[s]], R), Q[[s]], control$rank_tol)
      Q[[s]] <- pol$Q
      supported[s] <- pol$rank
      Y[[s]] <- Ut[[s]] %*% Q[[s]]
      HYs <- .sumcor_apply_H(L, Y[[s]])
      if (full) HYsum <- HYsum + HYs - HY[[s]]
      HY[[s]] <- HYs
    }
    objective <- .sumcor_objective(Y, HY, lambda)
    if (objective < previous - 1e-10 * max(1, abs(previous))) monotone <- FALSE

    stationarity <- numeric(S)
    constraint <- numeric(S)
    for (s in seq_len(S)) {
      R <- if (full) off[1L] * (HYsum - HY[[s]]) else
        .sumcor_gradient_input(s, HY, lambda)
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
        "sumcor_align: sweep %d objective=%.10g rel_change=%.3g stationarity=%.3g",
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

  # Common rotation ordering components by cross-domain agreement. It leaves
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
  n_pairs <- sum(lambda[upper.tri(lambda)] > 0)

  list(
    Q = Q,
    objective = objective,
    history = history,
    converged = converged,
    status = status,
    iterations = it,
    monotone = monotone,
    supported_rank = stats::setNames(supported, domain_names),
    component_agreement = eg$values / max(n_pairs, 1L)
  )
}

# -- Per-domain preparation ------------------------------------------------------

.sumcor_prepare_domain <- function(dom, events, w, control, groups, domain) {
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
  basis <- .sumcor_fit_basis(dom$X[used, , drop = FALSE], omega, control,
                             groups, domain)
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
  if (!length(lam) || max(lam) <= 0 || min(lam) <= 1e-12 * max(lam)) {
    stop("Domain `", domain, "`: the repeat-level covariance is singular; ",
         "increase `covariance_shrinkage` or reduce the basis rank.",
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
    p = p
  )
}

.sumcor_encoder <- function(prep, A, column_names, control) {
  P <- prep$basis$P
  W <- as.matrix(P %*% A)
  if (!is.null(prep$basis$scale)) W <- W / prep$basis$scale
  list(
    center = prep$basis$center,
    weights = W,
    column_names = column_names,
    n_columns = length(prep$basis$center),
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

.sumcor_encode <- function(encoder, X) {
  Z <- sweep(X, 2L, encoder$center, "-") %*% encoder$weights
  colnames(Z) <- paste0("sumcor_", seq_len(ncol(Z)))
  Z
}

.sumcor_event_scores <- function(prep, A, event_ids) {
  Z <- prep$Ubar %*% A
  Z[!prep$mask, ] <- NA_real_
  dimnames(Z) <- list(event_ids, paste0("sumcor_", seq_len(ncol(Z))))
  Z
}

# -- Fit ----------------------------------------------------------------------------

.sumcor_resolve_events_table <- function(events) {
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

.sumcor_fit <- function(X_list, design_list, domain_names, ids, events,
                        event_features, feature_groups, ncomp, control, call) {
  control <- .sumcor_resolve_control(control)
  ncomp <- .sumcor_int(ncomp, "ncomp")
  S <- length(X_list)
  if (S < 2L) stop("sumcor_align() needs at least two domains.", call. = FALSE)
  if (anyDuplicated(domain_names)) {
    stop("Domain names must be unique.", call. = FALSE)
  }
  events_in <- .sumcor_resolve_events_table(events)
  if (!is.null(feature_groups)) {
    if (!is.list(feature_groups) || length(feature_groups) != S) {
      stop("`feature_groups` must be a list with one entry per domain.",
           call. = FALSE)
    }
    if (!is.null(names(feature_groups)) &&
        all(domain_names %in% names(feature_groups))) {
      feature_groups <- feature_groups[domain_names]
    }
  }

  domains <- lapply(seq_len(S), function(i) {
    .sumcor_resolve_domain(X_list[[i]], design_list[[i]], ids,
                           domain_names[[i]])
  })
  names(domains) <- domain_names

  observed <- unlist(lapply(domains, function(d) {
    d$observations$event_id[d$observations$used]
  }), use.names = FALSE)
  universe <- sort(unique(observed), method = "radix")
  if (length(universe) < 2L) {
    stop("At least two distinct events are needed.", call. = FALSE)
  }
  ev <- data.frame(event_id = universe, stringsAsFactors = FALSE)
  attribute_names <- unique(c(ids$event_attributes,
                              setdiff(names(events_in), "event_id")))
  for (field in attribute_names) {
    ev[[field]] <- .sumcor_event_attribute(domains, events_in, field, universe)
  }
  w <- .sumcor_event_weights(control$event_weights, ev)
  ev$weight <- w

  preps <- lapply(seq_len(S), function(i) {
    .sumcor_prepare_domain(domains[[i]], ev, w, control,
                           if (is.null(feature_groups)) NULL else
                             feature_groups[[i]],
                           domain_names[[i]])
  })
  names(preps) <- domain_names
  ev$n_domains <- Reduce(`+`, lapply(preps, function(p) as.integer(p$mask)))

  ranks <- vapply(preps, `[[`, integer(1), "p")
  if (ncomp > min(ranks)) {
    stop("`ncomp = ", ncomp, "` exceeds the basis rank of domain `",
         names(ranks)[which.min(ranks)], "` (", min(ranks), ").", call. = FALSE)
  }

  pairs <- .sumcor_pair_weights(lapply(preps, `[[`, "mask"), w,
                                control$min_overlap)
  dimnames(pairs$overlap) <- dimnames(pairs$lambda) <-
    list(domain_names, domain_names)
  isolated <- domain_names[rowSums(pairs$lambda > 0) == 0]
  if (length(isolated)) {
    stop("Domain(s) ", paste0("`", isolated, "`", collapse = ", "),
         " share less than `min_overlap` event weight with every other domain.",
         call. = FALSE)
  }

  kernel <- .sumcor_kernel_factor(event_features, ev, control)
  solved <- .sumcor_solve(lapply(preps, `[[`, "Ut"), kernel$L, pairs$lambda,
                          ncomp, control, domain_names)

  encoders <- list()
  scores <- list()
  event_scores <- list()
  for (i in seq_len(S)) {
    nm <- domain_names[[i]]
    A <- preps[[i]]$C_inv_half %*% solved$Q[[i]]
    enc <- .sumcor_encoder(preps[[i]], A, colnames(domains[[i]]$X), control)
    Z <- .sumcor_encode(enc, domains[[i]]$X)
    rownames(Z) <- domains[[i]]$observations$observation_id
    encoders[[nm]] <- enc
    scores[[nm]] <- Z
    event_scores[[nm]] <- .sumcor_event_scores(preps[[i]], A, universe)
  }
  stack <- simplify2array(event_scores)
  reference <- apply(stack, c(1L, 2L), function(v) {
    if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
  })
  dimnames(reference) <- dimnames(event_scores[[1L]])

  if (!solved$converged) {
    warning(
      "sumcor_align() did not converge in ", control$max_iter,
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
      "are supported by cross-domain agreement in at least one domain; the ",
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
        observations_used = vapply(domains,
                                   function(d) sum(d$observations$used),
                                   integer(1)),
        observations_unlabelled = vapply(domains,
                                         function(d) sum(!d$observations$used),
                                         integer(1))
      )
    ),
    class = c("sumcor_align", "manifold_alignment")
  )
}

#' Event-keyed SUMCOR alignment with optional event-feature metric
#'
#' Learns one linear encoder per domain that maps observations into a shared
#' latent space. Domains are linked only through shared *events*: keyed items,
#' conditions or time points that several domains observed, possibly more than
#' once. The encoders maximise reproducible cross-domain agreement of
#' event-level means. Repeated observations of an event are averaged in the
#' agreement term but enter each domain's normalisation individually, so
#' directions dominated by repeat-to-repeat variability are down-weighted.
#' Optional event features add a metric over events, so that events with
#' similar descriptors also contribute to agreement.
#'
#' @section Model:
#' For domain `s`, a fitted linear basis gives `u = phi_s(x)`. With event
#' weights `w` (summing to one), event means `ubar_s`, repeat-level covariance
#' `C_s` (shrunk toward a scaled identity), and event metric
#' `H = I + beta * sum_m theta_m K_m` with
#' `K_m = F_m (F_m' F_m + eta_m I)^{-1} F_m'`, the encoders solve
#' \deqn{\max \sum_{s<t} \mathrm{tr}(A_s^\top B_{st} A_t)
#'   \quad \mathrm{s.t.} \quad A_s^\top C_s A_s = I,}
#' with `B_st = ubar_s' W^{1/2} H W^{1/2} ubar_t / m_st`, where `m_st` is the
#' event weight both domains observed. Unobserved events are masked by key.
#' The problem is solved by sequential polar updates, which never decrease the
#' objective; joint global optimality is not guaranteed for more than two
#' domains. The encoder `f_s(x) = phi_s(x) A_s` is stored explicitly, and the
#' training scores are produced by that same function.
#'
#' @section Identities:
#' Correspondence between domains is carried only by event identifiers. Row
#' order and row counts are never used to infer it. Observation and repeat
#' identifiers, and any requested observation attributes, are kept in the
#' result. Rows whose event is `NA` are excluded from fitting but still
#' receive scores.
#'
#' @section Event features:
#' Each element of `event_features` is a matrix with one row per event and
#' row names equal to event identifiers. Rows for events no domain observed
#' are dropped; events without a row in a family take no part in that
#' family's kernel. Feature preprocessing is fitted on the events supplied.
#' With `kernel_blocks` set, each kernel is block diagonal by that event
#' attribute. A kernel close to a multiple of the identity changes nothing and
#' triggers a warning.
#'
#' @param data A list of numeric observation-by-column matrices (one per
#'   domain) or a `hyperdesign`.
#' @param ... Method arguments. Unknown arguments are rejected.
#' @param design For list input, a list of data frames (one per domain, rows
#'   matching `data`) holding the identifier columns.
#' @param event_id Name of the design column with event identifiers.
#' @param repeat_id Optional design column with repeat identifiers. Without it,
#'   repeats are numbered by occurrence within each event.
#' @param observation_id Optional design column with unique observation
#'   identifiers. Defaults to row names or row positions.
#' @param event_attributes Optional design columns holding event-level
#'   attributes (one consistent value per event across all domains). They can
#'   be named by `event_weights` and `kernel_blocks` in [sumcor_control()].
#' @param observation_attributes Optional design columns copied into each
#'   domain's observation table (stored only).
#' @param events Optional data frame with an `event_id` column and further
#'   event-attribute columns.
#' @param event_features Optional named list of event-by-feature matrices
#'   (row names are event identifiers).
#' @param feature_groups For `basis = "grouped_pca"`, a list with one group
#'   label vector per domain (one entry per column; `NA` excludes a column).
#' @param ncomp Latent dimension.
#' @param control Settings from [sumcor_control()] or a named list of its
#'   arguments.
#' @return A `sumcor_align` object with per-domain `scores` (all observations,
#'   original order, rows named by observation identifier), `event_scores`,
#'   frozen `reference_scores`, `encoders`, identity tables, kernel
#'   diagnostics, `convergence` (objective history, stationarity, constraint
#'   error, supported rank), and `component_agreement`.
#' @seealso [calibrate_domain()] to add a domain against the frozen reference,
#'   [oos_predict()] to encode new observations.
#' @examples
#' set.seed(1)
#' n_events <- 40
#' latent <- matrix(rnorm(n_events * 3), n_events, 3)
#' make_domain <- function(p) {
#'   mix <- matrix(rnorm(3 * p), 3, p)
#'   rows <- rep(seq_len(n_events), times = 2)
#'   X <- latent[rows, ] %*% mix + matrix(rnorm(length(rows) * p, sd = 0.5),
#'                                        length(rows), p)
#'   list(X = X, design = data.frame(event_id = paste0("e", rows),
#'                                   repeat_id = rep(1:2, each = n_events)))
#' }
#' domains <- list(a = make_domain(20), b = make_domain(25),
#'                 c = make_domain(18))
#' fit <- sumcor_align(
#'   lapply(domains, `[[`, "X"),
#'   design = lapply(domains, `[[`, "design"),
#'   event_id = "event_id", repeat_id = "repeat_id",
#'   ncomp = 3,
#'   control = sumcor_control(beta = 0, basis_rank = 10)
#' )
#' fit$component_agreement
#' @export
sumcor_align <- function(data, ...) {
  UseMethod("sumcor_align")
}

.sumcor_reject_dots <- function(...) {
  dots <- list(...)
  if (length(dots)) {
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    nm[!nzchar(nm)] <- "<unnamed>"
    stop("Unknown argument(s): ", paste(nm, collapse = ", "), ".",
         call. = FALSE)
  }
}

.sumcor_ids <- function(event_id, repeat_id, observation_id, event_attributes,
                        observation_attributes) {
  list(
    event_id = .sumcor_name(event_id, "event_id", allow_null = FALSE),
    repeat_id = .sumcor_name(repeat_id, "repeat_id"),
    observation_id = .sumcor_name(observation_id, "observation_id"),
    event_attributes = .sumcor_names(event_attributes, "event_attributes"),
    observation_attributes = .sumcor_names(observation_attributes,
                                           "observation_attributes")
  )
}

#' @rdname sumcor_align
#' @export
sumcor_align.list <- function(
  data,
  design,
  event_id = "event_id",
  repeat_id = NULL,
  observation_id = NULL,
  event_attributes = NULL,
  observation_attributes = NULL,
  events = NULL,
  event_features = NULL,
  feature_groups = NULL,
  ncomp = 10L,
  control = sumcor_control(),
  ...
) {
  .sumcor_reject_dots(...)
  domain_names <- names(data)
  if (is.null(domain_names) || anyNA(domain_names) ||
      any(!nzchar(domain_names))) {
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
  ids <- .sumcor_ids(event_id, repeat_id, observation_id, event_attributes,
                     observation_attributes)
  .sumcor_fit(unname(data), unname(design), domain_names, ids, events,
              event_features, feature_groups, ncomp, control, match.call())
}

#' @rdname sumcor_align
#' @export
sumcor_align.hyperdesign <- function(
  data,
  event_id = "event_id",
  repeat_id = NULL,
  observation_id = NULL,
  event_attributes = NULL,
  observation_attributes = NULL,
  events = NULL,
  event_features = NULL,
  feature_groups = NULL,
  ncomp = 10L,
  control = sumcor_control(),
  ...
) {
  .sumcor_reject_dots(...)
  resolved <- resolve_hyperdesign(data)
  ids <- .sumcor_ids(event_id, repeat_id, observation_id, event_attributes,
                     observation_attributes)
  .sumcor_fit(unname(lapply(resolved$domains, `[[`, "x")),
              unname(lapply(resolved$domains, `[[`, "design")),
              resolved$domain_names, ids, events, event_features,
              feature_groups, ncomp, control, match.call())
}

#' @export
sumcor_align.default <- function(data, ...) {
  stop("sumcor_align() needs a list of matrices or a hyperdesign.",
       call. = FALSE)
}

#' @method print sumcor_align
#' @export
print.sumcor_align <- function(x, ...) {
  cat("Event-keyed SUMCOR alignment\n")
  roles <- x$domain_roles
  cat("  domains: ", length(x$domain_names), " (",
      sum(roles == "reference"), " reference, ",
      sum(roles == "calibrated"), " calibrated)\n", sep = "")
  cat("  events: ", nrow(x$events), "\n", sep = "")
  cat("  components: ", x$ncomp, "\n", sep = "")
  fam <- if (is.null(x$kernel$families)) "none" else
    paste(unique(x$kernel$families$family), collapse = ", ")
  cat("  beta: ", format(x$kernel$beta), " (event features: ", fam, ")\n",
      sep = "")
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
