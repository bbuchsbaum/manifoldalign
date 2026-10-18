# RAPID-MA dual (feature-side) coordinates ---------------------------------

.rapid_feature_weights <- function(X, weight) {
  if (inherits(X, "Matrix")) {
    if (identical(weight, "shifted")) {
      X <- as.matrix(X)
    } else {
      X <- .rapid_force_dgC(X)
      X@x <- if (identical(weight, "absolute")) abs(X@x) else pmax(X@x, 0)
      return(Matrix::drop0(X))
    }
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  switch(
    weight,
    positive = {
      X[X < 0] <- 0
      X
    },
    absolute = abs(X),
    shifted = {
      mins <- apply(X, 2L, min)
      sweep(X, 2L, mins, "-")
    }
  )
}

.rapid_feature_barycentres <- function(X, Z, weight, name) {
  .rapid_validate_numeric_matrix(X, name, n_rows = nrow(Z))
  vals <- if (inherits(X, "Matrix")) X@x else X
  if (any(!is.finite(vals))) {
    stop("`", name, "` must contain only finite values.", call. = FALSE)
  }
  W <- .rapid_feature_weights(X, weight)
  mass <- as.numeric(Matrix::colSums(W))
  numerator <- as.matrix(Matrix::crossprod(W, Z))
  coords <- numerator / mass
  zero_mass <- which(mass <= 0)
  if (length(zero_mass)) coords[zero_mass, ] <- NA_real_
  feature_names <- colnames(X)
  if (is.null(feature_names) || any(!nzchar(feature_names)) ||
      anyDuplicated(feature_names)) {
    feature_names <- paste0("feature_", seq_len(ncol(X)))
  }
  dimnames(coords) <- list(feature_names, colnames(Z))
  names(mass) <- feature_names
  list(coords = coords, mass = mass, zero_mass = as.integer(zero_mass))
}

.rapid_feature_scores_data <- function(fit, data) {
  domain_names <- fit$domain_names
  M <- length(domain_names)
  if (inherits(data, "hyperdesign")) {
    resolved <- resolve_hyperdesign(data)
    matrices <- lapply(resolved$domains, `[[`, "x")
    names(matrices) <- resolved$domain_names
    data <- matrices
  }
  if (!is.list(data) || inherits(data, "Matrix")) {
    stop("`data` must be a list with one matrix per fitted domain or a ",
         "`hyperdesign`.", call. = FALSE)
  }
  data_names <- names(data)
  if (!is.null(data_names) && all(nzchar(data_names))) {
    unknown <- setdiff(data_names, domain_names)
    if (length(unknown)) {
      stop("`data` names do not match fitted domains: ",
           paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    out <- rep(list(NULL), M)
    names(out) <- domain_names
    out[data_names] <- data
    return(out)
  }
  if (length(data) != M) {
    stop("`data` must have one entry per fitted domain (", M, ").",
         call. = FALSE)
  }
  names(data) <- domain_names
  data
}

#' Feature-side (dual) coordinates for a RAPID-MA embedding
#'
#' RAPID-MA embeds samples (for example voxels) into a shared latent space but
#' is not a linear projector, so it has no feature loadings. This function
#' supplies the dual view instead: every feature column (for example an
#' activation map over voxels) is placed at the weighted barycentre of the
#' latent coordinates of the samples it loads on. Because the latent space is
#' shared across domains, feature coordinates from different domains are
#' directly comparable, and features that distribute similarly over similarly
#' embedded samples land close together.
#'
#' For feature \eqn{j} in domain \eqn{d} with sample coordinates
#' \eqn{Z_d} and nonnegative weights \eqn{w_{ij}} derived from the data,
#' the coordinate is \eqn{m_{dj} = \sum_i w_{ij} z_{di} / \sum_i w_{ij}}.
#' Weights are derived from the data column by column according to `weight`:
#' `"positive"` keeps nonnegative values and zeroes negative ones,
#' `"absolute"` uses magnitudes, and `"shifted"` subtracts each column's
#' minimum so every entry is nonnegative. Sparse `Matrix` inputs stay sparse
#' for `"positive"` and `"absolute"`; `"shifted"` densifies.
#'
#' Columns with zero total weight have no barycentre and are returned as `NA`
#' rows; their indices are reported in `zero_mass` and a warning is raised.
#'
#' The matrices in `data` need not be the training matrices: any matrix with
#' one row per fitted sample of that domain can be embedded, so new maps over
#' the same samples can be placed after fitting. Pooling requires every
#' supplied domain to have the same number of columns (and identical column
#' names when both sides are named).
#'
#' @param fit A fitted `rapid_ma` object or a `rapid_ma_pair_fit`.
#' @param data A list of sample-by-feature matrices (base or `Matrix`), one per
#'   fitted domain in fitting order or named by domain, or a `hyperdesign`.
#'   Entries may be `NULL` to skip a domain. For a `rapid_ma_pair_fit`, supply
#'   `list(X_i, X_j)` or `list(i = X_i, j = X_j)`.
#' @param weight Nonnegative weighting rule applied column by column:
#'   `"positive"`, `"absolute"`, or `"shifted"`.
#' @param pool Logical; also return the equal-weight mean of each feature's
#'   per-domain coordinates across the supplied domains, plus a dispersion
#'   diagnostic.
#' @return An object of class `rapid_ma_feature_scores` with elements
#'   `scores` (per-domain feature-by-component matrices), `mass` (per-domain
#'   total weight per feature), `zero_mass` (per-domain indices of features
#'   without a barycentre), `pooled` (feature-by-component matrix, or `NULL`
#'   when `pool = FALSE`), `dispersion` (per-feature root-mean-square distance
#'   of the per-domain coordinates from the pooled coordinate; `NA` when fewer
#'   than two domains contribute), `n_domains` (per-feature count of
#'   contributing domains), `weight`, `domain_names`, and `ncomp`.
#' @examples
#' set.seed(1)
#' X <- list(
#'   a = matrix(abs(rnorm(40 * 6)), 40, 6),
#'   b = matrix(abs(rnorm(35 * 6)), 35, 6)
#' )
#' fit <- rapid_ma(X, ncomp = 3L,
#'                 control = rapid_ma_control(max_iter = 1L,
#'                                            prototype = list(unknown_prototypes = 4L)))
#' dual <- rapid_ma_feature_scores(fit, X)
#' dim(dual$pooled)
#' @export
rapid_ma_feature_scores <- function(
  fit,
  data,
  weight = c("positive", "absolute", "shifted"),
  pool = TRUE
) {
  weight <- match.arg(weight)
  if (length(pool) != 1L || !is.logical(pool) || is.na(pool)) {
    stop("`pool` must be TRUE or FALSE.", call. = FALSE)
  }
  if (inherits(fit, "rapid_ma_pair_fit")) fit <- fit$core
  if (!inherits(fit, "rapid_ma")) {
    stop("`fit` must be a `rapid_ma` or `rapid_ma_pair_fit` object.",
         call. = FALSE)
  }
  data <- .rapid_feature_scores_data(fit, data)
  supplied <- which(!vapply(data, is.null, logical(1)))
  if (!length(supplied)) {
    stop("`data` must supply at least one domain matrix.", call. = FALSE)
  }

  scores <- list()
  mass <- list()
  zero_mass <- list()
  for (m in supplied) {
    name <- fit$domain_names[[m]]
    Z <- fit$scores[[m]]
    label <- paste0("data[[\"", name, "\"]]")
    .rapid_validate_numeric_matrix(data[[m]], label)
    if (nrow(data[[m]]) != fit$domain_sizes[[m]]) {
      stop("Domain `", name, "`: `data` must have ", fit$domain_sizes[[m]],
           " rows (one per fitted sample).", call. = FALSE)
    }
    out <- .rapid_feature_barycentres(data[[m]], Z, weight, label)
    scores[[name]] <- out$coords
    mass[[name]] <- out$mass
    zero_mass[[name]] <- out$zero_mass
  }
  dropped <- vapply(zero_mass, length, integer(1))
  if (any(dropped > 0L)) {
    warning(
      "Features with zero total weight have no barycentre and are NA: ",
      paste0(names(dropped)[dropped > 0L], " (", dropped[dropped > 0L], ")",
             collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  pooled <- NULL
  dispersion <- NULL
  n_domains <- NULL
  if (isTRUE(pool)) {
    widths <- vapply(scores, nrow, integer(1))
    if (length(unique(widths)) != 1L) {
      stop("Pooling requires the same number of features in every supplied ",
           "domain; widths are ", paste(widths, collapse = ", "),
           ". Use `pool = FALSE` for per-domain coordinates only.",
           call. = FALSE)
    }
    named <- Filter(function(x) !is.null(colnames(x)), data[supplied])
    if (length(named) > 1L) {
      reference <- colnames(named[[1L]])
      for (x in named[-1L]) {
        if (!identical(colnames(x), reference)) {
          stop("Pooling requires identical feature names across domains.",
               call. = FALSE)
        }
      }
    }
    stacked <- simplify2array(scores)
    if (length(dim(stacked)) == 2L) {
      dim(stacked) <- c(dim(stacked), 1L)
    }
    n_domains <- as.integer(apply(!is.na(stacked[, 1L, , drop = FALSE]), 1L, sum))
    pooled <- apply(stacked, c(1L, 2L), mean, na.rm = TRUE)
    pooled[n_domains == 0L, ] <- NA_real_
    squared <- sweep(stacked, c(1L, 2L), pooled, "-")^2
    dispersion <- sqrt(apply(squared, 1L, sum, na.rm = TRUE) /
                         pmax(n_domains, 1L))
    dispersion[n_domains < 2L] <- NA_real_
    dimnames(pooled) <- dimnames(scores[[1L]])
    names(dispersion) <- names(n_domains) <- rownames(pooled)
  }

  structure(
    list(
      scores = scores,
      mass = mass,
      zero_mass = zero_mass,
      pooled = pooled,
      dispersion = dispersion,
      n_domains = n_domains,
      weight = weight,
      domain_names = names(scores),
      ncomp = fit$ncomp
    ),
    class = "rapid_ma_feature_scores"
  )
}

#' @method print rapid_ma_feature_scores
#' @export
print.rapid_ma_feature_scores <- function(x, ...) {
  cat("RAPID-MA feature-side coordinates\n")
  cat("  domains: ", length(x$domain_names), " (",
      paste(x$domain_names, collapse = ", "), ")\n", sep = "")
  cat("  features per domain: ",
      paste(vapply(x$scores, nrow, integer(1)), collapse = ", "), "\n",
      sep = "")
  cat("  components: ", x$ncomp, "\n", sep = "")
  cat("  weight: ", x$weight, "\n", sep = "")
  dropped <- sum(vapply(x$zero_mass, length, integer(1)))
  if (dropped > 0L) cat("  features without barycentre: ", dropped, "\n", sep = "")
  if (!is.null(x$pooled)) {
    cat("  pooled: ", nrow(x$pooled), " features, median dispersion ",
        format(stats::median(x$dispersion, na.rm = TRUE), digits = 4), "\n",
        sep = "")
  } else {
    cat("  pooled: no\n")
  }
  invisible(x)
}
