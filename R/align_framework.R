#' Generic alignment adapter interface and transform helpers
#'
#' This file defines a minimal, method-agnostic interface for lifting
#' pairwise manifold alignment algorithms to the multi-domain setting.
#' Adapters for concrete methods (e.g., GRASP, PARROT, GPA, KEMA) should
#' implement these generics for their aligner classes.

#' Create a simple aligner descriptor
#'
#' @param name Character name of the algorithm
#' @param group One of "O" (orthogonal), "GL" (linear), "perm" (assignment),
#'   or "kernel" (embedding-only)
#' @param supports_multi Logical; if TRUE, adapter should also implement
#'   a native multi-domain entry point via fit_many()
#' @param ... Additional fields stored on the object
#' @return An object of class c("<name>_aligner", "aligner")
#' @export
new_aligner <- function(name, group, supports_multi = FALSE, ...) {
  structure(
    list(name = name, group = group, supports_multi = supports_multi, ...),
    class = c(paste0(name, "_aligner"), "aligner")
  )
}

#' Capabilities for an aligner
#' @export
aligner_capabilities <- function(algo) UseMethod("aligner_capabilities")

#' @export
aligner_capabilities.aligner <- function(algo) {
  list(group = algo$group, supports_multi = isTRUE(algo$supports_multi))
}

# -- Adapter generics -------------------------------------------------------

#' Fit a pairwise alignment model
#' @export
fit_pair <- function(algo, X_i, X_j, links = NULL, ...) UseMethod("fit_pair")

#' @export
fit_pair.default <- function(algo, X_i, X_j, links = NULL, ...) {
  stop("fit_pair() is not implemented for this aligner class.")
}

#' Extract a relative transform from a pairwise fit
#' @export
relative_transform <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  UseMethod("relative_transform")
}

#' @export
relative_transform.default <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  stop("relative_transform() is not implemented for object of class ",
       paste(class(fit), collapse = ", "))
}

#' Apply a transform to data
#' @export
apply_transform <- function(transform, X, ...) UseMethod("apply_transform")

#' @export
apply_transform.default <- function(transform, X, ...) {
  if (inherits(transform, "align_transform")) {
    return(apply_transform.align_transform(transform, X, ...))
  }
  stop("apply_transform() is not implemented for object of class ",
       paste(class(transform), collapse = ", "))
}

#' Pairwise loss (for weighting/diagnostics)
#' @export
pair_loss <- function(fit, X_i = NULL, X_j = NULL, ...) UseMethod("pair_loss")

#' @export
pair_loss.default <- function(fit, X_i = NULL, X_j = NULL, ...) NA_real_

#' Invert a transform
#' @export
invert_transform <- function(transform, ...) UseMethod("invert_transform")

#' @export
invert_transform.default <- function(transform, ...) {
  if (inherits(transform, "align_transform")) {
    return(invert_transform.align_transform(transform, ...))
  }
  stop("invert_transform() is not implemented for object of class ",
       paste(class(transform), collapse = ", "))
}

#' Compose two transforms (a->b and b->c)
#' @export
compose_transform <- function(T_ab, T_bc, ...) UseMethod("compose_transform")

#' @export
compose_transform.default <- function(T_ab, T_bc, ...) {
  if (inherits(T_ab, "align_transform") && inherits(T_bc, "align_transform")) {
    return(compose_transform.align_transform(T_ab, T_bc, ...))
  }
  stop("compose_transform() is not implemented for provided objects")
}

#' Out-of-sample projection/prediction
#' @export
oos_predict <- function(fit_or_transform, newX, side, ...) UseMethod("oos_predict")

#' @export
oos_predict.default <- function(fit_or_transform, newX, side, ...) {
  stop("oos_predict() is not implemented for object of class ",
       paste(class(fit_or_transform), collapse = ", "))
}

#' Latent dimension selected by the method (if applicable)
#' @export
latent_dim <- function(fit, ...) UseMethod("latent_dim")

#' @export
latent_dim.default <- function(fit, ...) NULL

#' Native multi-domain entry point (optional)
#' @export
fit_many <- function(algo, domains, ...) UseMethod("fit_many")

#' @export
fit_many.default <- function(algo, domains, ...) {
  stop("fit_many() is not implemented for this aligner class.")
}

# -- Transform S3 -----------------------------------------------------------

#' Construct an alignment transform object
#'
#' @param group Group identifier: "O", "GL", or "perm" (or "kernel" surrogate)
#' @param op The underlying operator (e.g., kxk matrix, or assignment matrix)
#' @param from Source index/label
#' @param to Target index/label
#' @param k Latent dimension (if applicable)
#' @param dim Optional problem dimensions metadata
#' @return An object of class c("align_transform_<group>", "align_transform")
#' @export
new_align_transform <- function(group, op, from, to, k = NULL, dim = NULL) {
  structure(
    list(group = group, op = op, from = from, to = to, k = k, dim = dim),
    class = c(paste0("align_transform_", group), "align_transform")
  )
}

#' @export
apply_transform.align_transform <- function(transform, X, ...) {
  grp <- transform$group
  op <- transform$op
  if (grp %in% c("O", "GL")) {
    X <- as.matrix(X)
    op <- as.matrix(op)
    if (ncol(X) != nrow(op)) {
      stop("apply_transform: dimension mismatch: ncol(X) != nrow(op)")
    }
    return(X %*% op)
  } else if (grp == "perm") {
    # Permutation/assignment acts on rows: op (n_to x n_from) times X (n_from x d)
    X <- as.matrix(X)
    op <- as.matrix(op)
    if (ncol(op) != nrow(X)) {
      stop("apply_transform[perm]: dimension mismatch: ncol(op) != nrow(X)")
    }
    return(op %*% X)
  }
  stop("apply_transform: not implemented for group '", grp, "'")
}

#' @export
invert_transform.align_transform <- function(transform, ...) {
  grp <- transform$group
  op <- transform$op
  if (grp == "O") {
    inv <- t(op)
  } else if (grp == "GL") {
    inv <- tryCatch(solve(op), error = function(e) {
      stop("invert_transform: GL operator is not invertible: ", e$message)
    })
  } else if (grp == "perm") {
    # For soft maps, use transpose with row-normalization as pseudo-inverse
    A <- as.matrix(op)
    At <- t(A)
    rs <- rowSums(At)
    rs[!is.finite(rs) | rs == 0] <- 1
    inv <- At / rs
  } else {
    stop("invert_transform: not implemented for group '", grp, "'")
  }
  new_align_transform(grp, inv, from = transform$to, to = transform$from, k = transform$k, dim = transform$dim)
}

#' @export
compose_transform.align_transform <- function(T_ab, T_bc, ...) {
  if (T_ab$group != T_bc$group) stop("compose_transform: group mismatch")
  grp <- T_ab$group
  if (grp %in% c("O", "GL")) {
    op <- as.matrix(T_ab$op) %*% as.matrix(T_bc$op)
    return(new_align_transform(grp, op, from = T_ab$from, to = T_bc$to, k = T_ab$k %||% T_bc$k))
  } else if (grp == "perm") {
    # Row-operator composition: if T_ab is (b x a) and T_bc is (c x b), result is (c x a)
    A <- as.matrix(T_ab$op)
    B <- as.matrix(T_bc$op)
    if (ncol(B) != nrow(A)) {
      stop("compose_transform[perm]: inner dimensions do not match")
    }
    op <- B %*% A
    return(new_align_transform("perm", op, from = T_ab$from, to = T_bc$to))
  }
  stop("compose_transform: not implemented for group '", grp, "'")
}

#' Lightweight infix for NULL coalescing
#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
