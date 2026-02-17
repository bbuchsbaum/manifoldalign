test_that("createSimFun preserves sparsity and symmetry", {
  S <- matrix(c(1, 0.2, 0.3,
                0.2, 1, 0.4,
                0.3, 0.4, 1), nrow = 3, byrow = TRUE,
              dimnames = list(letters[1:3], letters[1:3]))

  expect_warning(simfun_warn <- createSimFun(S, na_value = 5),
                 "ignores non-zero na_value", fixed = TRUE)
  expect_type(simfun_warn, "closure")

  simfun <- createSimFun(S)
  labels <- c("a", "b", NA, "c", "ghost")
  M <- simfun(labels)

  expect_s4_class(M, "dgCMatrix")
  expect_true(Matrix::isSymmetric(M))
  expect_equal(as.numeric(Matrix::diag(M)), rep(0, length(labels)))
  expect_gt(Matrix::nnzero(M[1:3, 1:3]), 0)
  expect_equal(Matrix::nnzero(M[5, ]), 0)
  expect_equal(Matrix::nnzero(M[, 5]), 0)
})

test_that("lowrank_align agrees across solver backends", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("RSpectra")
  skip_if_not_installed("tibble")

  set.seed(42)
  X1 <- matrix(rnorm(5 * 3), 5, 3)
  X2 <- matrix(rnorm(4 * 3), 4, 3)

  labels1 <- c("L1", "L2", NA, "L1", "L3")
  labels2 <- c("L1", NA, "L3", "L2")

  if (!("package:tibble" %in% search())) {
    suppressWarnings(library(tibble))
  }

  d1 <- multidesign::multidesign(X1, tibble::tibble(y = labels1))
  d2 <- multidesign::multidesign(X2, tibble::tibble(y = labels2))
  hd <- multidesign::hyperdesign(list(block1 = d1, block2 = d2))

  label_pool <- sort(unique(stats::na.omit(c(labels1, labels2))))
  S <- diag(length(label_pool))
  dimnames(S) <- list(label_pool, label_pool)
  simfun <- createSimFun(S)

  total_samples <- length(labels1) + length(labels2)
  total_features <- ncol(X1) + ncol(X2)

  res_exp <- lowrank_align(hd, y,
                           simfun = simfun,
                           ncomp = 2,
                           solver = "explicit",
                           sv_thresh = 0.5,
                           lambda = 0.01)

  res_op <- lowrank_align(hd, y,
                          simfun = simfun,
                          ncomp = 2,
                          solver = "operator",
                          sv_thresh = 0.5,
                          lambda = 0.01)

  expect_s3_class(res_exp, c("lowrank_align", "multiblock_biprojector"))
  expect_s3_class(res_op, c("lowrank_align", "multiblock_biprojector"))

  expect_equal(dim(res_exp$s), c(total_samples, 2))
  expect_equal(dim(res_op$s), c(total_samples, 2))
  expect_equal(nrow(res_exp$v), total_features)
  expect_equal(nrow(res_op$v), total_features)

  expect_equal(sum(is.na(res_exp$labels)), sum(is.na(c(labels1, labels2))))
  expect_equal(sum(is.na(res_op$labels)), sum(is.na(c(labels1, labels2))))

  overlap <- svd(crossprod(res_exp$s, res_op$s))
  expect_gt(min(overlap$d), 0.8)
})

test_that("lowrank_align recovers a linear latent space on synthetic data", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("tibble")

  set.seed(20240501)
  n <- 80
  d_lat <- 2
  d1 <- 5
  d2 <- 6

  # Latent 2D Gaussian blobs with two classes
  n1 <- n %/% 2
  n2 <- n - n1
  latent_c1 <- matrix(rnorm(n1 * d_lat, mean = -0.7, sd = 0.4), n1, d_lat)
  latent_c2 <- matrix(rnorm(n2 * d_lat, mean = 0.7, sd = 0.4), n2, d_lat)
  Z <- rbind(latent_c1, latent_c2)
  labels <- factor(c(rep("A", n1), rep("B", n2)))

  # Two linear views with different projections and noise
  W1 <- matrix(rnorm(d_lat * d1), d_lat, d1)
  W2 <- matrix(rnorm(d_lat * d2), d_lat, d2)
  X1 <- Z %*% W1 + matrix(rnorm(n * d1, sd = 0.05), n, d1)
  X2 <- Z %*% W2 + matrix(rnorm(n * d2, sd = 0.05), n, d2)

  if (!("package:tibble" %in% search())) {
    suppressWarnings(library(tibble))
  }

  d1_md <- multidesign::multidesign(X1, tibble::tibble(y = labels))
  d2_md <- multidesign::multidesign(X2, tibble::tibble(y = labels))
  hd <- multidesign::hyperdesign(list(view1 = d1_md, view2 = d2_md))

  # Simple similarity: same label = 1, different = 0
  lab_all <- c(labels, labels)
  lab_pool <- sort(unique(stats::na.omit(lab_all)))
  S <- diag(length(lab_pool))
  dimnames(S) <- list(lab_pool, lab_pool)
  simfun <- createSimFun(S)

  # Fit lowrank_align and extract 2D embedding
  fit <- lowrank_align(hd, y,
                       simfun = simfun,
                       ncomp = 2,
                       solver = "explicit",
                       sv_thresh = 0.2,
                       lambda = 0.01,
                       mu = 0)

  expect_s3_class(fit, c("lowrank_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$s), 2 * n)
  expect_equal(ncol(fit$s), 2)
  expect_true(all(is.finite(fit$s)))

  # Align embedding from first view back to latent space via Procrustes
  emb1 <- fit$s[seq_len(n), , drop = FALSE]

  if (!requireNamespace("vegan", quietly = TRUE)) {
    skip("vegan package required for Procrustes analysis - install with: install.packages('vegan')")
  }

  proc <- vegan::procrustes(as.matrix(Z), as.matrix(emb1), symmetric = FALSE)
  err <- proc$ss
  cor_lat <- sqrt(1 - proc$ss / sum(Z^2))

  # On clean linear data, lowrank_align should recover the latent space up to
  # an orthogonal transform, though PCA can still outperform depending on the
  # similarity structure and regularization choices.
  pca_scores <- stats::prcomp(X1, rank. = 2)$x
  proc_pca <- vegan::procrustes(as.matrix(Z), as.matrix(pca_scores), symmetric = FALSE)
  err_pca <- proc_pca$ss

  expect_lt(err, 10 * err_pca)         # avoid pathological regressions vs PCA
  expect_gt(cor_lat, 0.8)              # correlation with latent coordinates should be strong
})
