# Shared helpers for alignment algorithm tests --------------------------------

create_alignment_hyperdesign <- function(X1, X2) {
  stopifnot(nrow(X1) == nrow(X2))
  domains <- list(
    list(x = X1, design = data.frame(node_id = seq_len(nrow(X1)))),
    list(x = X2, design = data.frame(node_id = seq_len(nrow(X2))))
  )
  names(domains) <- c("domain1", "domain2")
  class(domains) <- "hyperdesign"
  domains
}

create_multi_alignment_hyperdesign <- function(mats) {
  domains <- lapply(seq_along(mats), function(i) {
    Xi <- mats[[i]]
    list(
      x = Xi,
      design = data.frame(node_id = seq_len(nrow(Xi)))
    )
  })
  names(domains) <- paste0("domain", seq_along(domains))
  class(domains) <- "hyperdesign"
  domains
}

skip_if_missing_cd_deps <- function() {
  for (pkg in c("multidesign", "multivarious", "neighborweights", "tibble")) {
    testthat::skip_if_not_installed(pkg)
  }
}

alignment_error <- function(E1, E2, P, Q) {
  diff <- E1 %*% Q - E2[P, , drop = FALSE]
  sum(diff^2)
}

permutation_accuracy <- function(P, target_perm) {
  stopifnot(length(P) == length(target_perm))
  mean(P == target_perm)
}

orthogonality_residual <- function(Q) {
  max(abs(crossprod(Q) - diag(ncol(Q))))
}

inverse_permutation <- function(perm) {
  inv <- integer(length(perm))
  inv[perm] <- seq_along(perm)
  inv
}

generate_identity_alignment <- function(n = 20, d = 3, noise_sd = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X1 <- matrix(rnorm(n * d), n, d)
  X2 <- X1 + matrix(rnorm(n * d, sd = noise_sd), n, d)
  perm <- seq_len(n)
  list(X1 = X1, X2 = X2, perm = perm, map = perm)
}

generate_permuted_alignment <- function(n = 20, d = 3, shift = 1L, noise_sd = 0.01, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base_angles <- 2 * pi * seq_len(n) / n
  X1 <- cbind(cos(base_angles), sin(base_angles), sin(2 * base_angles))
  if (d > 3) {
    extra <- matrix(rnorm(n * (d - 3)), n, d - 3)
    X1 <- cbind(X1, extra)
  }
  perm <- ((seq_len(n) + shift - 1L) %% n) + 1L
  X2 <- X1[perm, , drop = FALSE] + matrix(rnorm(n * ncol(X1), sd = noise_sd), n)
  list(X1 = X1, X2 = X2, perm = perm, map = inverse_permutation(perm))
}

compute_test_embeddings <- function(hd, ncomp, sigma, use_laplacian = TRUE, knn = NULL) {
  suppressWarnings(
    manifoldalign:::compute_cone_embeddings(unclass(hd), ncomp = ncomp, sigma = sigma,
                                            use_laplacian = use_laplacian, knn = knn)
  )
}
