test_that("spectral_mnn_align aligns near-identical domains (unsupervised anchors)", {
  set.seed(123)
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.02, seed = 123)

  fit <- spectral_mnn_align(
    list(domain1 = dat$X1, domain2 = dat$X2),
    preproc = multivarious::center(),
    ncomp = 3,
    control = spectral_mnn_align_control(
      knn = 10,
      q_descriptors = 12,
      mnn_k = 10,
      max_anchors = 200,
      dist_quantile = 0.3,
      refine_rounds = 1,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("spectral_mnn_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(ncol(fit$s), 3)
  expect_equal(nrow(fit$s), 2 * nrow(dat$X1))
  expect_equal(nrow(fit$v), ncol(dat$X1) + ncol(dat$X2))

  n <- nrow(dat$X1)
  E1 <- fit$s[seq_len(n), , drop = FALSE]
  E2 <- fit$s[(n + 1):(2 * n), , drop = FALSE]

  # Aligned embeddings should be close for corresponding samples
  mse <- mean(rowSums((E1 - E2)^2))
  expect_lt(mse, 0.15)

  # Rotation should be orthogonal
  R2 <- fit$rotations[[2]]
  expect_true(is.matrix(R2))
  expect_lt(max(abs(crossprod(R2) - diag(ncol(R2)))), 1e-6)
})

test_that("spectral_mnn_align respects explicit correspondences", {
  set.seed(321)
  dat <- generate_identity_alignment(n = 30, d = 5, noise_sd = 0.03, seed = 321)
  corr <- data.frame(
    domain_i = 1L, index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L, index_j = seq_len(nrow(dat$X2))
  )

  fit <- spectral_mnn_align(
    list(d1 = dat$X1, d2 = dat$X2),
    correspondences = corr,
    preproc = multivarious::center(),
    ncomp = 2,
    control = spectral_mnn_align_control(
      knn = 8,
      q_descriptors = 8,
      max_anchors = 128,
      refine_rounds = 0,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("spectral_mnn_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(ncol(fit$s), 2)

  expect_true(is.data.frame(fit$anchors[[2]]))
  expect_true(all(grepl("correspondences", fit$anchors[[2]]$mode)))
  expect_gte(nrow(fit$anchors[[2]]), 10)
})

test_that("spectral_mnn_align aligns 3+ domains to a reference", {
  set.seed(77)
  n <- 28
  z <- matrix(rnorm(n * 3), n, 3)
  X1 <- z %*% matrix(rnorm(3 * 7), 3, 7) + matrix(rnorm(n * 7, sd = 0.02), n, 7)
  X2 <- z %*% matrix(rnorm(3 * 6), 3, 6) + matrix(rnorm(n * 6, sd = 0.02), n, 6)
  X3 <- z %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(n * 5, sd = 0.02), n, 5)

  fit <- spectral_mnn_align(
    list(a = X1, b = X2, c = X3),
    preproc = multivarious::center(),
    ncomp = 3,
    ref_idx = 1L,
    control = spectral_mnn_align_control(
      knn = 9,
      q_descriptors = 10,
      max_anchors = 200,
      dist_quantile = 0.3,
      refine_rounds = 1,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("spectral_mnn_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$s), 3 * n)
  expect_equal(ncol(fit$s), 3)
  expect_equal(length(fit$rotations), 3)
  expect_true(all(vapply(fit$rotations, is.matrix, logical(1))))
})

test_that("spectral_mnn_align accepts unquoted label columns (tidy-style)", {
  set.seed(20260209)
  n <- 34
  X1 <- matrix(rnorm(n * 6), n, 6)
  X2 <- X1 + matrix(rnorm(n * 6, sd = 0.03), n, 6)
  condition_vec <- factor(rep(c("A", "B"), length.out = n))

  hd <- as_hyperdesign(
    list(domain1 = X1, domain2 = X2),
    labels = list(condition_vec, condition_vec),
    label_name = "condition"
  )

  fit <- spectral_mnn_align(
    hd,
    y = condition, # should NOT require a `condition` variable in scope
    preproc = multivarious::center(),
    ncomp = 2,
    control = spectral_mnn_align_control(
      knn = 10,
      q_descriptors = 8,
      max_pairs_per_label = 20,
      refine_rounds = 0,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("spectral_mnn_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_true(any(grepl("labels", fit$anchors[[2]]$mode)))
})
