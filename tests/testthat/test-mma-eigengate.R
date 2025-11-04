test_that("eigenvalue gating constrains matches within clusters", {
  skip_if_missing_cd_deps()
  set.seed(101)
  # Circle points yield near-multiplicity in Laplacian eigenvalues
  n <- 80; d <- 2
  theta <- seq(0, 2*pi, length.out = n + 1)[- (n+1)]
  X1 <- cbind(cos(theta), sin(theta))
  X2 <- X1 # identical structure

  # Run alignment with large K to expose clustered eigenvalues
  K <- 8
  res <- mma_align_multiple(list(X1, X2, X2),
                            ref_idx = 1,
                            ncomp = K,
                            sigma = 0.5,
                            embedding = "ctd",
                            normalize = "hypersphere",
                            signature_method = "w2",
                            signature_gate_multiplicity = TRUE,
                            em_max_iter = 5)

  expect_true(inherits(res, "mma_align_multiple"))
  # Extract eigenvalues from attributes
  # Using internal function attributes available on embeddings via extras$K and ref idx
  # We can recompute lam from eigensignature alignment stored in extras$hist_align (pairing info present)
  aligns <- res$hist_align
  expect_equal(length(aligns), 3L)
  ag <- aligns[[2]]
  perm <- ag$perm

  # Recompute embeddings to extract lam_keep attributes
  emb <- suppressWarnings(manifoldalign:::compute_mma_embeddings(
    strata = list(list(x = X1), list(x = X2)), ncomp = K, sigma = 0.5, knn = NULL,
    embedding = "ctd", normalize = "hypersphere", target_var = 0.95))
  lam_ref <- attr(emb[[1]], "lam_keep")
  lam_g   <- attr(emb[[2]], "lam_keep")

  cluster_eigs <- function(lam, tau = NULL) {
    lam <- as.numeric(lam)
    if (length(lam) <= 1) return(rep(1L, length(lam)))
    dif <- diff(lam)
    if (is.null(tau)) {
      md <- stats::median(abs(dif[is.finite(dif)]))
      if (!is.finite(md) || md <= 0) md <- mean(abs(dif[is.finite(dif)]))
      if (!is.finite(md) || md <= 0) md <- 1e-6
      tau <- 1e-3 * md
    }
    grp <- integer(length(lam)); grp[1] <- 1L; g <- 1L
    for (i in 2:length(lam)) { if (abs(lam[i] - lam[i-1]) > tau) g <- g + 1L; grp[i] <- g }
    grp
  }
  ref_grp <- cluster_eigs(lam_ref)
  tgt_grp <- cluster_eigs(lam_g)
  # Assert that gating enforced within-cluster matches
  expect_true(all(ref_grp == tgt_grp[perm]))
})

test_that("pair_confidence is present and bounded", {
  skip_if_missing_cd_deps()
  set.seed(102)
  n <- 40; d <- 3; K <- 6
  X <- matrix(rnorm(n*d), n, d)
  X1 <- X
  X2 <- X %*% {M <- matrix(rnorm(d*d), d, d); svd(M)$u %*% t(svd(M)$v)}
  X3 <- X %*% {M <- matrix(rnorm(d*d), d, d); svd(M)$u %*% t(svd(M)$v)}
  res <- mma_align_multiple(list(X1, X2, X3), ncomp = K, signature_method = "hybrid",
                            sigma = 0.6, embedding = "ctd", normalize = "hypersphere",
                            em_max_iter = 10)
  ag <- res$hist_align[[2]]
  expect_true(is.numeric(ag$pair_confidence))
  expect_equal(length(ag$pair_confidence), K)
  expect_true(all(ag$pair_confidence >= 0 - 1e-12))
  expect_true(all(ag$pair_confidence <= 1 + 1e-12))
})
