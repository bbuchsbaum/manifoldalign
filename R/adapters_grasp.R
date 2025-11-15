#' GRASP adapter for align_many()
#'
#' Provides the minimal pairwise interface required by align_many() for the
#' GRASP algorithm. Uses the existing GRASP building blocks to estimate a
#' relative orthogonal transform between two domains.

#' Construct a GRASP aligner descriptor
#' @return an object of class c("grasp_aligner", "aligner")
#' @export
grasp_aligner <- function() {
  new_aligner("grasp", group = "O", supports_multi = FALSE)
}

#' Fit GRASP on a pair of domains (adapter)
#' @export
fit_pair.grasp_aligner <- function(algo, X_i, X_j, links = NULL,
                                   ncomp = 30,
                                   q_descriptors = 100,
                                   sigma = 0.73,
                                   lambda = 0.1,
                                   use_laplacian = TRUE,
                                   solver = c("linear", "hungarian", "auction"),
                                   ...) {
  solver <- match.arg(solver)

  strata <- list(list(x = as.matrix(X_i)), list(x = as.matrix(X_j)))
  bases <- compute_grasp_basis(strata, ncomp = ncomp, use_laplacian = use_laplacian)
  descs <- compute_grasp_descriptors(bases, q_descriptors = q_descriptors, sigma = sigma)

  align_res <- align_grasp_bases(bases[[1]], bases[[2]], descs[[1]], descs[[2]], lambda = lambda)
  assign_res <- compute_grasp_assignment(bases[[1]], bases[[2]], descs[[1]], descs[[2]],
                                         M = align_res$rotation,
                                         distance_method = "cosine",
                                         solver_method = solver)

  # Pair loss: descriptor alignment residual
  phi1 <- bases[[1]]$vectors
  phi2 <- bases[[2]]$vectors
  F_T_Phi1 <- Matrix::crossprod(descs[[1]], phi1)
  G_T_Phi2 <- Matrix::crossprod(descs[[2]], phi2)
  resid <- F_T_Phi1 - G_T_Phi2 %*% align_res$rotation
  loss <- as.numeric(Matrix::norm(resid, "F")^2)
  k <- ncol(phi1)

  structure(list(
    rotation = as.matrix(align_res$rotation),
    mapping_matrix = as.matrix(assign_res$mapping_matrix),
    assignment = assign_res$assignment,
    loss = loss,
    k = k
  ), class = "grasp_pair_fit")
}

#' Extract relative transform from a GRASP pairwise fit
#' @export
relative_transform.grasp_pair_fit <- function(fit, from = c("i", "j"), to = c("j", "i"), ...) {
  from <- match.arg(from)
  to <- match.arg(to)
  if (from == to) stop("from and to must differ")
  op <- if (from == "i" && to == "j") {
    fit$rotation
  } else {
    t(fit$rotation) # inverse for orthogonal
  }
  new_align_transform("O", op, from = from, to = to, k = fit$k)
}

#' Pairwise loss for GRASP fit
#' @export
pair_loss.grasp_pair_fit <- function(fit, X_i = NULL, X_j = NULL, ...) {
  fit$loss
}

#' Latent dimension for GRASP fit
#' @export
latent_dim.grasp_pair_fit <- function(fit, ...) fit$k

