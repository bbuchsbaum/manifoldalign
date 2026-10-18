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
#'   "kernel" (embedding-only), or "encoder" (a native multi-domain method that
#'   learns per-domain encoders and is never synchronised pairwise). Only
#'   "O", "GL" and "perm" have relative transforms that [align_many()] can
#'   synchronise; "kernel" and "encoder" aligners must be native
#'   (`supports_multi = TRUE`).
#' @param supports_multi Logical; if TRUE, adapter should also implement
#'   a native multi-domain entry point via fit_many()
#' @param ... Additional fields stored on the object
#' @return An object of class c("<name>_aligner", "aligner")
#' @examples
#' algo <- new_aligner("my_method", group = "O")
#' @export
new_aligner <- function(name, group, supports_multi = FALSE, ...) {
  structure(
    list(name = name, group = group, supports_multi = supports_multi, ...),
    class = c(paste0(name, "_aligner"), "aligner")
  )
}

#' Capabilities for an aligner
#'
#' @param algo An aligner object created by new_aligner or an aligner constructor
#' @return A list with elements group (character) and supports_multi (logical).
#'   Aligners that declare a fuller contract also report, in this order,
#'   `encoder` (`"linear"`, `"nonlinear"`, or `"none"` for row transports),
#'   `supports_predict`, `supports_calibration`,
#'   `requires_auxiliary_at_predict`, and `invertible`. Fields an aligner does
#'   not declare are omitted rather than guessed. Their meaning:
#'   \describe{
#'     \item{`encoder`}{how a row of a fitted domain reaches the shared space:
#'       a fixed `"linear"` map of its features, a `"nonlinear"` (for example
#'       kernel) map, or `"none"` when the method only produces row
#'       correspondences, transports or maps between derived coordinates.}
#'     \item{`supports_predict`}{new rows of a fitted domain can be projected
#'       with a validated out-of-sample method.}
#'     \item{`supports_calibration`}{a new domain can be added to a frozen fit.}
#'     \item{`requires_auxiliary_at_predict`}{prediction needs inputs beyond the
#'       new rows' features (for example graphs, labels or training data that
#'       the fit does not store).}
#'     \item{`invertible`}{the map from a domain to the shared space can be
#'       inverted to recover the domain's features.}
#'   }
#' @examples
#' algo <- new_aligner("test", group = "O")
#' caps <- aligner_capabilities(algo)
#' @export
aligner_capabilities <- function(algo) UseMethod("aligner_capabilities")

#' @export
aligner_capabilities.aligner <- function(algo) {
  caps <- list(group = algo$group, supports_multi = isTRUE(algo$supports_multi))
  contract <- c("encoder", "supports_predict", "supports_calibration",
                "requires_auxiliary_at_predict", "invertible")
  for (field in contract) {
    if (!is.null(algo[[field]])) caps[[field]] <- algo[[field]]
  }
  caps
}

# -- Adapter generics -------------------------------------------------------

#' Fit a pairwise alignment model
#' @param algo An aligner object created by an aligner constructor
#' @param X_i First domain data matrix (samples x features)
#' @param X_j Second domain data matrix (samples x features)
#' @param links Optional correspondence links between domains
#' @param ... Additional arguments passed to methods
#' @return A pairwise fit object whose class depends on the method (e.g.,
#'   grasp_pair_fit, parrot_pair_fit). Contains alignment results including
#'   transform parameters and loss.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 25, 2)
#' X2 <- matrix(rnorm(50), 25, 2)
#' algo <- grasp_aligner()
#' fit <- fit_pair(algo, X1, X2)
#' }
#' @export
fit_pair <- function(algo, X_i, X_j, links = NULL, ...) UseMethod("fit_pair")

#' @export
fit_pair.default <- function(algo, X_i, X_j, links = NULL, ...) {
  stop("fit_pair() is not implemented for this aligner class.")
}

#' Extract a relative transform from a pairwise fit
#' @param fit A pairwise fit object returned by fit_pair
#' @param from Source domain identifier ("i" or "j")
#' @param to Target domain identifier ("j" or "i")
#' @param ... Additional arguments passed to methods
#' @return An align_transform object (see new_align_transform) containing the
#'   transformation operator, group type, and source/target indices.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 25, 2)
#' X2 <- matrix(rnorm(50), 25, 2)
#' algo <- grasp_aligner()
#' fit <- fit_pair(algo, X1, X2)
#' transform <- relative_transform(fit, from = "i", to = "j")
#' }
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
#'
#' @param transform An alignment transform object (typically from relative_transform
#'   or new_align_transform)
#' @param X A data matrix to transform (samples x features)
#' @param ... Additional arguments passed to methods
#' @return A transformed data matrix. For orthogonal/linear groups, returns
#'   `X %*% op` (action on the feature axis); for the permutation group,
#'   returns `op %*% X` (action on the sample axis). See the transform
#'   contract in [new_align_transform()]. A `Matrix` input (for example a
#'   sparse `dgCMatrix`) yields a `Matrix` result; base matrices yield base
#'   matrices.
#' @examples
#' # Create a simple orthogonal transform
#' R <- diag(3)
#' tr <- new_align_transform("O", R, from = 1, to = 2)
#' X <- matrix(rnorm(12), 4, 3)
#' X_transformed <- apply_transform(tr, X)
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
#' @param fit A pairwise fit object returned by fit_pair
#' @param X_i Optional first domain data for loss computation
#' @param X_j Optional second domain data for loss computation
#' @param ... Additional arguments passed to methods
#' @return A single numeric value representing the pairwise alignment loss, or
#'   NA_real_ if not available.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 25, 2)
#' X2 <- matrix(rnorm(50), 25, 2)
#' algo <- grasp_aligner()
#' fit <- fit_pair(algo, X1, X2)
#' loss <- pair_loss(fit)
#' }
#' @export
pair_loss <- function(fit, X_i = NULL, X_j = NULL, ...) UseMethod("pair_loss")

#' @export
pair_loss.default <- function(fit, X_i = NULL, X_j = NULL, ...) NA_real_

#' Invert a transform
#'
#' @param transform An alignment transform object to invert
#' @param ... Additional arguments passed to methods
#' @return An align_transform object representing the inverse transform, with
#'   from and to swapped.
#' @examples
#' # Create and invert an orthogonal transform
#' R <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), 2, 2)
#' tr <- new_align_transform("O", R, from = 1, to = 2)
#' tr_inv <- invert_transform(tr)
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
#'
#' @param T_ab An alignment transform from domain a to domain b
#' @param T_bc An alignment transform from domain b to domain c
#' @param ... Additional arguments passed to methods
#' @return An align_transform object representing the composed transform from
#'   domain a to domain c.
#' @examples
#' # Compose two orthogonal transforms
#' R1 <- diag(3)
#' R2 <- diag(3)
#' T_ab <- new_align_transform("O", R1, from = 1, to = 2)
#' T_bc <- new_align_transform("O", R2, from = 2, to = 3)
#' T_ac <- compose_transform(T_ab, T_bc)
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
#' @param fit_or_transform A fit object or transform object
#' @param newX New data matrix to project (samples x features)
#' @param side Which domain side the new data belongs to ("i" or "j")
#' @param ... Additional arguments passed to methods
#' @return A matrix of predicted/projected coordinates for the new data in the
#'   aligned space.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 25, 2)
#' X2 <- matrix(rnorm(50), 25, 2)
#' newX <- matrix(rnorm(20), 10, 2)
#' # Create a simple 2-block projector (loadings stacked by feature block)
#' v1 <- prcomp(X1)$rotation[, 1, drop = FALSE]
#' v2 <- prcomp(X2)$rotation[, 1, drop = FALSE]
#' v <- rbind(v1, v2)
#' s <- rbind(X1 %*% v1, X2 %*% v2)
#' fit <- structure(
#'   list(
#'     v = v,
#'     s = s,
#'     sdev = apply(s, 2, sd),
#'     preproc = NULL,
#'     block_indices = list(
#'       i = seq_len(ncol(X1)),
#'       j = (ncol(X1) + 1L):(ncol(X1) + ncol(X2))
#'     )
#'   ),
#'   class = "multiblock_biprojector"
#' )
#' pred <- oos_predict(fit, newX, side = "i")
#' }
#' @export
oos_predict <- function(fit_or_transform, newX, side, ...) UseMethod("oos_predict")

#' @export
oos_predict.default <- function(fit_or_transform, newX, side, ...) {
  stop("oos_predict() is not implemented for object of class ",
       paste(class(fit_or_transform), collapse = ", "))
}

#' @export
oos_predict.multiblock_biprojector <- function(fit_or_transform, newX, side = c("i", "j"), ...) {
  blocks <- fit_or_transform$block_indices
  if (is.null(blocks)) {
    stop("oos_predict.multiblock_biprojector: missing block_indices in fit object.", call. = FALSE)
  }

  side_idx <- side
  if (is.character(side_idx)) {
    side_idx <- side_idx[1]
    if (side_idx %in% c("i", "j")) {
      side_idx <- if (side_idx == "i") 1L else 2L
    } else if (is.list(blocks) && !is.null(names(blocks)) && side_idx %in% names(blocks)) {
      side_idx <- match(side_idx, names(blocks))
    } else {
      stop("side must be one of 'i' or 'j' (or a valid block name).", call. = FALSE)
    }
  } else {
    side_idx <- as.integer(side_idx[1])
  }

  if (side_idx < 1L) stop("Invalid side index.", call. = FALSE)

  # Resolve feature indices for this block
  feat_idx <- NULL
  if (is.list(blocks)) {
    if (side_idx > length(blocks)) stop("Invalid side index for block_indices.", call. = FALSE)
    feat_idx <- blocks[[side_idx]]
  } else if (is.matrix(blocks) && ncol(blocks) >= 2L) {
    if (side_idx > nrow(blocks)) stop("Invalid side index for block_indices.", call. = FALSE)
    feat_idx <- seq.int(blocks[side_idx, 1], blocks[side_idx, 2])
  } else {
    stop("Unsupported block_indices format; expected list or 2-column matrix.", call. = FALSE)
  }

  V <- fit_or_transform$v
  if (is.null(V)) stop("Fit object is missing loadings matrix 'v'.", call. = FALSE)

  V_block <- V[feat_idx, , drop = FALSE]
  X <- as.matrix(newX)

  # Apply preprocessing for this block if available
  pre <- fit_or_transform$preproc
  if (!is.null(pre)) {
    if (inherits(pre, "concat_pre_processor")) {
      X <- multivarious::transform(pre, X, colind = feat_idx)
    } else {
      if (is.list(pre) && length(pre) >= side_idx) pre <- pre[[side_idx]]
      if (inherits(pre, "prepper") || inherits(pre, "pre_processor")) {
        X <- multivarious::init_transform(pre, X)
      } else if (is.function(pre)) {
        X <- pre(X)
      }
    }
  }

  if (ncol(X) != nrow(V_block)) {
    stop("newX has ", ncol(X), " columns but this block expects ", nrow(V_block), ".", call. = FALSE)
  }
  X %*% V_block
}

#' Latent dimension selected by the method (if applicable)
#' @param fit A fit object returned by fit_pair or fit_many
#' @param ... Additional arguments passed to methods
#' @return An integer giving the latent dimension selected by the method, or
#'   NULL if not applicable.
#' @examples
#' \donttest{
#' set.seed(1)
#' X1 <- matrix(rnorm(50), 25, 2)
#' X2 <- matrix(rnorm(50), 25, 2)
#' algo <- grasp_aligner()
#' fit <- fit_pair(algo, X1, X2)
#' dim <- latent_dim(fit)
#' }
#' @export
latent_dim <- function(fit, ...) UseMethod("latent_dim")

#' @export
latent_dim.default <- function(fit, ...) NULL

#' Native multi-domain entry point (optional)
#' @param algo An aligner object that supports multi-domain fitting
#' @param domains A hyperdesign object or list of domain data
#' @param ... Additional arguments passed to methods
#' @return The return value depends on the aligner method. Typically a
#'   multiblock_biprojector or similar alignment result object.
#' @examples
#' \donttest{
#' set.seed(1)
#' library(multidesign)
#' X1 <- matrix(rnorm(60), 30, 2)
#' X2 <- matrix(rnorm(60), 30, 2)
#' labels <- sample(c("A", "B"), 30, TRUE)
#' design1 <- data.frame(labels = labels)
#' design2 <- data.frame(labels = labels)
#' md1 <- multidesign(X1, design1)
#' md2 <- multidesign(X2, design2)
#' hd <- hyperdesign(list(d1 = md1, d2 = md2))
#' algo <- kema_aligner()
#' result <- fit_many(algo, hd, y = labels, ncomp = 2)
#' }
#' @export
fit_many <- function(algo, domains, ...) UseMethod("fit_many")

#' @export
fit_many.default <- function(algo, domains, ...) {
  stop("fit_many() is not implemented for this aligner class.")
}

# -- Transform S3 -----------------------------------------------------------

#' Construct an alignment transform object
#'
#' @section Transform contract:
#' Every relative and per-domain transform in the package follows one
#' contract, from which [apply_transform()], [compose_transform()],
#' [invert_transform()], the synchronisation backends of [align_many()] and
#' its diagnostics are all derived. A transform `T` from `from = a` to
#' `to = b` maps observations of domain `a` onto observations of domain `b`:
#' \describe{
#'   \item{`"O"`, `"GL"` (action on columns)}{`op` is a `p_a x p_b` matrix
#'     acting on the feature axis from the right: `apply_transform(T, X_a) =
#'     X_a %*% op ~ X_b`, row by row.}
#'   \item{`"perm"` (action on rows)}{`op` is an `n_b x n_a` (soft)
#'     assignment acting on the sample axis from the left:
#'     `apply_transform(T, X_a) = op %*% X_a ~ X_b`.}
#' }
#' Composition applies the first transform first:
#' `apply_transform(compose_transform(T_ab, T_bc), X_a) =
#' apply_transform(T_bc, apply_transform(T_ab, X_a))`, so the composed operator
#' is `op_ab %*% op_bc` for `"O"`/`"GL"` and `op_bc %*% op_ab` for `"perm"`.
#'
#' In [align_many()] a pairwise edge transform `g_ij` goes from domain `i` to
#' domain `j` (`X_i g_ij ~ X_j`), and each per-domain transform `T_i` goes from
#' domain `i` to the shared consensus space (`Z ~ apply_transform(T_i, X_i)`
#' for every `i`). Consistency of an edge with the consensus therefore reads
#' `compose_transform(g_ij, T_j) = T_i`, that is `T_i = g_ij %*% T_j` for
#' `"O"`/`"GL"` and `T_i = T_j %*% g_ij` for `"perm"`. Synchronisation solves
#' these relations in the least-squares sense, and the reported edge residual
#' is the root-mean-square entry of their difference. A triangle `(i, j, k)` is
#' cycle-consistent when `compose_transform(g_ij, g_jk) = g_ik`.
#'
#' @param group Group identifier: "O", "GL", or "perm" (or "kernel" surrogate)
#' @param op The underlying operator (see section 'Transform contract').
#' @param from Source index/label
#' @param to Target index/label
#' @param k Latent dimension (if applicable)
#' @param dim Optional problem dimensions metadata
#' @return An object of class c("align_transform_<group>", "align_transform")
#' @examples
#' # Create a simple orthogonal transform
#' R <- diag(3)
#' tr <- new_align_transform("O", R, from = 1, to = 2)
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
  # Sparse `Matrix` inputs are kept as `Matrix` objects; base inputs are
  # returned as base matrices.
  x_is_Matrix <- inherits(X, "Matrix")
  if (!x_is_Matrix) X <- as.matrix(X)
  if (grp %in% c("O", "GL")) {
    # Column action: X (n x p_from) %*% op (p_from x p_to)
    if (!inherits(op, "Matrix")) op <- as.matrix(op)
    if (ncol(X) != nrow(op)) {
      stop("apply_transform: dimension mismatch: ncol(X) != nrow(op)")
    }
    out <- X %*% op
  } else if (grp == "perm") {
    # Row action: op (n_to x n_from) %*% X (n_from x d)
    if (!inherits(op, "Matrix")) op <- as.matrix(op)
    if (ncol(op) != nrow(X)) {
      stop("apply_transform[perm]: dimension mismatch: ncol(op) != nrow(X)")
    }
    out <- op %*% X
  } else {
    stop("apply_transform: not implemented for group '", grp, "'")
  }
  if (x_is_Matrix) return(out)
  if (inherits(out, "Matrix")) out <- as.matrix(out)
  out
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

#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b
