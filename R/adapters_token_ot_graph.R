#' Token-OT graph aligner adapter for align_many()
#'
#' Provides a pairwise interface for token_ot_graph_align() using the unified
#' aligner adapter framework (fit_pair/relative_transform).

#' Construct a Token-OT Graph aligner descriptor
#'
#' @return an object of class c("token_ot_graph_aligner", "aligner")
#' @examples
#' algo <- token_ot_graph_aligner()
#' aligner_capabilities(algo)
#' @export
token_ot_graph_aligner <- function() {
  new_aligner("token_ot_graph", group = "perm", supports_multi = FALSE)
}

.normalize_rows_sparse <- function(M) {
  if (!inherits(M, "Matrix")) {
    return(.normalize_rows_dense(as.matrix(M)))
  }
  rs <- Matrix::rowSums(M)
  rs[!is.finite(rs) | rs == 0] <- 1
  Dinv <- Matrix::Diagonal(x = 1 / rs)
  Dinv %*% M
}

#' Fit Token-OT Graph alignment on a pair of domains
#'
#' @param algo A token_ot_graph_aligner object.
#' @param X_i First domain data matrix (nodes x features).
#' @param X_j Second domain data matrix (nodes x features).
#' @param links Optional anchor links. Accepts either list(vec1=, vec2=) or
#'   list(v1, v2), where shared values indicate anchors and NAs denote unknowns.
#' @param ncomp Latent dimension used for baseline embeddings (default: 20).
#' @param control Control settings from token_ot_graph_align_control().
#' @param ... Additional arguments (unused).
#' @return An object of class token_ot_graph_pair_fit with fields transport
#'   (sparse coupling), assignment, objective, and dimensions.
#' @export
fit_pair.token_ot_graph_aligner <- function(algo, X_i, X_j, links = NULL,
                                           ncomp = 20L,
                                           control = token_ot_graph_align_control(),
                                           ...) {
  ctrl <- resolve_token_ot_graph_align_control(control)
  if (!is.numeric(ncomp) || length(ncomp) != 1L || is.na(ncomp) || ncomp < 1) {
    stop("`ncomp` must be a positive scalar.", call. = FALSE)
  }
  ncomp <- as.integer(ncomp)

  anchor_info <- NULL
  if (!is.null(links)) {
    if (is.list(links) && !is.null(links$vec1) && !is.null(links$vec2)) {
      v1 <- links$vec1
      v2 <- links$vec2
    } else if (is.list(links) && length(links) == 2 && is.vector(links[[1]]) && is.vector(links[[2]])) {
      v1 <- links[[1]]
      v2 <- links[[2]]
    } else {
      stop("links must be list(vec1=, vec2=) or list(v1, v2)", call. = FALSE)
    }
    if (length(v1) != nrow(X_i) || length(v2) != nrow(X_j)) {
      stop("links vectors must match nrow(X_i) and nrow(X_j).", call. = FALSE)
    }
    anchor_info <- list(vec1 = v1, vec2 = v2)
  }

  fit <- .token_ot_graph_pair_fit(
    Xs = as.matrix(X_i),
    Xt = as.matrix(X_j),
    anchors = anchor_info,
    ncomp = ncomp,
    ctrl = ctrl
  )

  objective <- NA_real_
  if (!is.null(fit$diagnostics$objective_scaled)) {
    objective <- as.numeric(fit$diagnostics$objective_scaled)
  }

  structure(
    list(
      transport = fit$transport,
      assignment = fit$assignment,
      objective = objective,
      n1 = nrow(X_i),
      n2 = nrow(X_j),
      diagnostics = fit$diagnostics
    ),
    class = "token_ot_graph_pair_fit"
  )
}

#' @rdname relative_transform
#' @export
relative_transform.token_ot_graph_pair_fit <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) stop("from and to must differ", call. = FALSE)

  W <- fit$transport
  op <- if (from == "i" && to == "j") {
    # Want (n2 x n1): normalize columns of W == normalize rows of t(W)
    .normalize_rows_sparse(Matrix::t(W))
  } else {
    .normalize_rows_sparse(W)
  }
  new_align_transform("perm", op, from = from, to = to, dim = c(fit$n1, fit$n2))
}

#' @rdname pair_loss
#' @export
pair_loss.token_ot_graph_pair_fit <- function(fit, X_i = NULL, X_j = NULL, ...) {
  fit$objective %||% NA_real_
}

