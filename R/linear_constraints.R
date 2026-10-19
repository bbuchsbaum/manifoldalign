#' Linear correspondence constraints for alignment coefficients
#'
#' Describe soft agreement between linear summaries of two domains' projection
#' coefficients. For maps H1 and H2, coefficients W1 and W2, and effective row
#' weights d, the added objective is
#' `strength * sum(d * rowSums((H1 %*% W1 - H2 %*% W2)^2))`.
#'
#' @param domains Character vector of two distinct domain names, matching the
#'   names in the alignment input.
#' @param maps List of two finite numeric matrices (base or `Matrix`), with equal
#'   positive row counts. Columns index each domain's processed input features.
#'   Different domains may have different numbers of features. If row names
#'   are supplied, both maps must have identical, unique row names in order.
#'   Named columns must exactly match the processed data's feature names.
#' @param weights Nonnegative finite confidence weights, one per matched row,
#'   or a scalar recycled across rows. Zero weights omit individual links.
#' @param strength Finite nonnegative multiplier for this constraint.
#' @param reduction `"weighted_mean"` divides weights by their sum;
#'   `"sum"` retains their absolute magnitude. Neither option rescales the maps
#'   or harmonizes their units. An all-zero weighted mean is undefined.
#' @param name Nonempty label for diagnostics; labels must be unique within a fit.
#' @param space Coordinate contract. Currently only `"processed"` is supported:
#'   maps act on coefficients after input preprocessing but before SSMA's
#'   internal compression. Maps are never centered or standardized as samples.
#'
#' @details
#' These are homogeneous quadratic soft constraints, not additional samples or
#' hard equalities. They do not modify sample neighborhoods or their degree
#' normalization. Selecting feature coefficients, averaging groups, and applying
#' linear measurement maps all use the same interface. For processed data X,
#' maps `A %*% X` express correspondence between weighted sample summaries A.
#'
#' With externally reduced input `Z = X %*% P`, use maps `H %*% P` and compose
#' returned coefficients as `P %*% W`. With covariance maps, calculate
#' `G %*% crossprod(X, Z) / (nrow(X) - 1)` for centered data instead of forming
#' the full feature covariance. Covariance estimation and preprocessing belong
#' to the caller and must respect validation folds.
#'
#' Unnamed columns are interpreted positionally. Raw-space constraints are not
#' automatically transformed through arbitrary preprocessing. For a known
#' linear preprocessing map S, a raw-coefficient functional H becomes H S;
#' intercepts are outside this homogeneous coefficient contract.
#'
#' @return A `linear_alignment_constraint` object. Pass a list of these objects
#'   to [ssma_align()] via `constraints`. Currently SSMA is the consumer.
#' @md
#' @export
#' @examples
#' h <- diag(3)
#' link <- linear_correspondence(c("a", "b"), list(h, h), name = "features")
#' smooth <- linear_penalty("a", diff(diag(3)), name = "smoothness")
linear_correspondence <- function(domains, maps, weights = 1, strength = 1,
                                  reduction = c("weighted_mean", "sum"),
                                  name = "correspondence", space = "processed") {
  .new_linear_constraint(domains, maps, weights, strength, match.arg(reduction),
                         name, space, "correspondence")
}

#' Linear regularization of alignment coefficients
#'
#' A single-domain counterpart to [linear_correspondence()]. Adds
#' `strength * sum(d * rowSums((map %*% W)^2))` to the objective. For example,
#' a feature-graph incidence matrix penalizes differences between neighboring
#' feature coefficients. This does not establish preservation of the complete
#' feature covariance geometry.
#'
#' @param domain One domain name.
#' @param map A finite numeric matrix whose columns index processed features.
#' @inheritParams linear_correspondence
#' @return A `linear_alignment_constraint` object.
#' @md
#' @export
linear_penalty <- function(domain, map, weights = 1, strength = 1,
                           reduction = c("weighted_mean", "sum"),
                           name = "penalty", space = "processed") {
  .new_linear_constraint(domain, list(map), weights, strength, match.arg(reduction),
                         name, space, "penalty")
}

.linear_finite_matrix <- function(x) {
  numeric_matrix <- (is.matrix(x) && is.numeric(x) && !is.complex(x)) ||
    inherits(x, "dMatrix") || inherits(x, "iMatrix")
  if (!numeric_matrix || length(dim(x)) != 2L || any(dim(x) < 1L)) return(FALSE)
  values <- if (inherits(x, "sparseMatrix") && "x" %in% methods::slotNames(x)) {
    methods::slot(x, "x")
  } else x
  all(is.finite(values))
}

.new_linear_constraint <- function(domains, maps, weights, strength, reduction,
                                    name, space, type) {
  nd <- if (identical(type, "correspondence")) 2L else 1L
  if (!is.character(domains) || length(domains) != nd || anyNA(domains) ||
      any(!nzchar(domains)) || anyDuplicated(domains)) {
    stop("`domains` must contain ", nd, " distinct, nonempty domain name(s).", call. = FALSE)
  }
  if (!is.list(maps) || length(maps) != nd ||
      !all(vapply(maps, .linear_finite_matrix, logical(1)))) {
    stop("`maps` must contain finite numeric matrices with positive dimensions.", call. = FALSE)
  }
  q <- nrow(maps[[1L]])
  if (any(vapply(maps, nrow, integer(1)) != q)) {
    stop("Constraint maps must have the same number of rows.", call. = FALSE)
  }
  rn <- lapply(maps, rownames)
  if (any(vapply(rn, function(x) !is.null(x), logical(1)))) {
    if (!all(vapply(rn, function(x) identical(x, rn[[1L]]), logical(1))) ||
        is.null(rn[[1L]]) || anyNA(rn[[1L]]) || any(!nzchar(rn[[1L]])) ||
        anyDuplicated(rn[[1L]])) {
      stop("Map row names must be identical, unique, nonempty and in matching order.", call. = FALSE)
    }
  }
  for (h in maps) {
    cn <- colnames(h)
    if (!is.null(cn) && (anyNA(cn) || any(!nzchar(cn)) || anyDuplicated(cn))) {
      stop("Map column names must be unique and nonempty.", call. = FALSE)
    }
  }
  if (!is.numeric(weights) || is.complex(weights) || !is.null(dim(weights)) || !length(weights) ||
      !(length(weights) %in% c(1L, q)) || any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite, nonnegative and scalar or one per map row.", call. = FALSE)
  }
  weights <- rep_len(as.numeric(weights), q)
  if (!is.numeric(strength) || is.complex(strength) || !is.null(dim(strength)) || length(strength) != 1L || !is.finite(strength) || strength < 0) {
    stop("`strength` must be one finite nonnegative number.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be one nonempty diagnostic label.", call. = FALSE)
  }
  if (!identical(space, "processed")) {
    stop("Only space='processed' is supported; transform raw-space maps explicitly.", call. = FALSE)
  }
  reduction <- match.arg(reduction, c("weighted_mean", "sum"))
  if (reduction == "weighted_mean") {
    if (!any(weights > 0)) stop("A weighted mean requires positive total weight.", call. = FALSE)
    # Normalize before summing to avoid overflow for large confidence weights.
    effective <- weights / max(weights)
    effective <- effective / sum(effective)
  } else {
    effective <- weights
  }
  if (any(!is.finite(strength * effective))) {
    stop("Effective constraint weights overflow; rescale weights or strength.", call. = FALSE)
  }
  structure(list(type = type, domains = domains, maps = maps, weights = weights,
                 effective_weights = effective, strength = as.numeric(strength),
                 reduction = reduction, name = name, space = space),
            class = "linear_alignment_constraint")
}

.resolve_linear_constraints <- function(constraints, domains, domain_names) {
  if (is.null(constraints)) return(list())
  if (!is.list(constraints) || inherits(constraints, "linear_alignment_constraint")) {
    stop("`constraints` must be a list of linear constraint objects.", call. = FALSE)
  }
  out <- lapply(constraints, function(x) {
    if (!inherits(x, "linear_alignment_constraint") ||
        length(x$type) != 1L || is.na(x$type) ||
        !x$type %in% c("correspondence", "penalty")) {
      stop("Each constraint must be made by linear_correspondence() or linear_penalty().", call. = FALSE)
    }
    # Revalidate mutable R objects at the fit boundary.
    x <- .new_linear_constraint(x$domains, x$maps, x$weights, x$strength,
                                x$reduction, x$name, x$space, x$type)
    ids <- match(x$domains, domain_names)
    if (anyNA(ids)) stop("Unknown constraint domain: ", paste(x$domains[is.na(ids)], collapse = ", "), call. = FALSE)
    for (j in seq_along(ids)) {
      h <- x$maps[[j]]
      d <- domains[[ids[j]]]$x
      if (ncol(h) != ncol(d)) stop("Constraint map columns must match processed features for ", x$domains[j], ".", call. = FALSE)
      if (!is.null(colnames(h)) && !identical(colnames(h), colnames(d))) {
        stop("Constraint map column names must match processed feature names in order for ", x$domains[j], ".", call. = FALSE)
      }
    }
    x$ids <- ids
    x
  })
  labels <- vapply(out, `[[`, character(1), "name")
  if (anyDuplicated(labels)) stop("Constraint names must be unique within a fit.", call. = FALSE)
  out
}

.compile_linear_constraints <- function(constraints, basis) {
  lapply(constraints, function(x) {
    maps <- lapply(seq_along(x$ids), function(j) {
      q <- basis[[x$ids[j]]]$P
      # NULL is the identity when compression is disabled; no p-by-p allocation.
      if (is.null(q)) x$maps[[j]] else x$maps[[j]] %*% q
    })
    if (!all(vapply(maps, .linear_finite_matrix, logical(1)))) {
      stop("Nonfinite reduced constraint map; rescale the supplied maps.", call. = FALSE)
    }
    loss <- vapply(seq_along(maps), function(j) {
      original <- sum(x$maps[[j]]^2)
      kept <- sum(maps[[j]]^2)
      if (!is.finite(original) || !is.finite(kept)) {
        stop("Constraint map energy overflow; rescale the supplied maps.", call. = FALSE)
      }
      if (original == 0) 0 else max(0, min(1, 1 - kept / original))
    }, numeric(1))
    list(ids = x$ids, maps = maps, signs = if (length(maps) == 2L) c(1, -1) else 1,
         weights = x$effective_weights, strength = x$strength, name = x$name,
         type = x$type, domains = x$domains, reduction = x$reduction,
         n_rows = nrow(maps[[1L]]), active_rows = sum(x$effective_weights > 0),
         discarded_energy_fraction = stats::setNames(loss, x$domains))
  })
}

.linear_constraint_apply <- function(terms, blocks) {
  out <- lapply(blocks, function(x) matrix(0, nrow(x), ncol(x)))
  for (term in terms) {
    if (term$strength == 0) next
    residual <- matrix(0, term$n_rows, ncol(blocks[[1L]]))
    for (j in seq_along(term$ids)) {
      residual <- residual + term$signs[j] * as.matrix(term$maps[[j]] %*% blocks[[term$ids[j]]])
    }
    weighted <- residual * (term$strength * term$weights)
    for (j in seq_along(term$ids)) {
      i <- term$ids[j]
      out[[i]] <- out[[i]] + term$signs[j] * as.matrix(crossprod(term$maps[[j]], weighted))
    }
  }
  if (!all(vapply(out, function(x) all(is.finite(x)), logical(1)))) {
    stop("Linear constraint product overflow; rescale maps or strengths.", call. = FALSE)
  }
  out
}

.add_linear_constraint_matrix <- function(A, terms, ranks) {
  offsets <- c(0L, cumsum(ranks))
  for (term in terms) {
    if (term$strength == 0) next
    for (j in seq_along(term$ids)) for (k in seq_along(term$ids)) {
      i <- term$ids[j]
      l <- term$ids[k]
      ri <- seq.int(offsets[i] + 1L, offsets[i + 1L])
      rl <- seq.int(offsets[l] + 1L, offsets[l + 1L])
      weighted <- Matrix::Diagonal(x = term$strength * term$weights) %*% term$maps[[k]]
      block <- term$signs[j] * term$signs[k] * crossprod(term$maps[[j]], weighted)
      A[ri, rl] <- A[ri, rl] + as.matrix(block)
    }
  }
  if (any(!is.finite(A))) stop("Linear constraint matrix overflow; rescale maps or strengths.", call. = FALSE)
  A
}

.linear_constraint_diagnostics <- function(terms, coefficient_blocks) {
  lapply(terms, function(term) {
    residual <- matrix(0, term$n_rows, ncol(coefficient_blocks[[1L]]))
    for (j in seq_along(term$ids)) {
      residual <- residual + term$signs[j] * as.matrix(term$maps[[j]] %*% coefficient_blocks[[term$ids[j]]])
    }
    per_component <- colSums(residual^2 * term$weights)
    list(name = term$name, type = term$type, domains = term$domains,
         reduction = term$reduction, strength = term$strength,
         n_rows = term$n_rows, active_rows = term$active_rows,
         effective_weight_sum = sum(term$weights),
         discarded_energy_fraction = term$discarded_energy_fraction,
         residual_per_component = per_component,
         residual = sum(per_component), objective = term$strength * sum(per_component))
  })
}
