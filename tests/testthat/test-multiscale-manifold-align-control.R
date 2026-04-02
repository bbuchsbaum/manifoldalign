library(testthat)

test_that("multiscale_manifold_align_control validates core fields", {
  ctrl <- multiscale_manifold_align_control(
    enabled = TRUE,
    backend = "spectral",
    rank_per_domain = 64,
    seed = 7
  )

  expect_s3_class(ctrl, "multiscale_manifold_align_control")
  expect_true(isTRUE(ctrl$enabled))
  expect_equal(ctrl$backend, "spectral")
  expect_equal(ctrl$rank_per_domain, 64L)

  expect_error(multiscale_manifold_align_control(rank_per_domain = 0), ">= 1")
  expect_error(multiscale_manifold_align_control(eps_rank = 2), "eps_rank")
  expect_error(multiscale_manifold_align_control(cg_tol = 0), "positive")
})

test_that("resolve_multiscale_manifold_align_control merges and rejects unknown fields", {
  resolved <- manifoldalign:::resolve_multiscale_manifold_align_control(list(
    enabled = TRUE,
    backend = "randomized_dwt",
    randomized_sketch_size = 128
  ))

  expect_s3_class(resolved, "multiscale_manifold_align_control")
  expect_true(isTRUE(resolved$enabled))
  expect_equal(resolved$backend, "randomized_dwt")
  expect_equal(resolved$randomized_sketch_size, 128L)

  expect_error(
    manifoldalign:::resolve_multiscale_manifold_align_control(list(not_a_field = 1)),
    "Unknown multiscale control argument"
  )
})

test_that("multiscale_manifold_align spectral backend returns valid alignment", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(101)
  dat <- generate_identity_alignment(n = 28, d = 5, noise_sd = 0.03, seed = 101)
  corr <- data.frame(
    domain_i = 1L, index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L, index_j = seq_len(nrow(dat$X2)),
    weight = 1
  )

  fit <- suppressWarnings(
    multiscale_manifold_align(
      list(d1 = dat$X1, d2 = dat$X2),
      correspondences = corr,
      preproc = multivarious::center(),
      ncomp = 3,
      control = multiscale_manifold_align_control(
        enabled = TRUE,
        backend = "spectral",
        rank_per_domain = 16L,
        knn = 8L,
        eigen_k_max = 32L,
        verbose = FALSE
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("multiscale_manifold_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(ncol(fit$s), 3L)
  expect_equal(nrow(fit$s), 2L * nrow(dat$X1))
  expect_equal(fit$backend, "spectral")
  expect_true(length(fit$level_dims) >= 1L)
  expect_true(all(diff(fit$level_dims) <= 0))
  expect_true(ncol(fit$v) == 3L)
})

test_that("multiscale randomized backend approximates the reduced diffusion spectrum", {
  set.seed(17)
  A <- crossprod(matrix(rnorm(40 * 25), 40, 25)) + diag(0.1, 25)
  ctrl <- multiscale_manifold_align_control(
    enabled = TRUE,
    backend = "randomized_dwt",
    randomized_sketch_size = 10L,
    randomized_oversample = 4L,
    randomized_power_iters = 2L,
    seed = 17L,
    verbose = FALSE
  )

  exact <- manifoldalign:::.msa_diffusion_eig_exact(A, 6L)
  randomized <- manifoldalign:::.msa_diffusion_eig_randomized(A, 6L, ctrl = ctrl)

  rel_err <- abs(exact$values - randomized$values) / pmax(abs(exact$values), 1e-8)
  expect_true(all(is.finite(rel_err)))
  expect_lt(max(rel_err), 0.05)
})

test_that("multiscale randomized backend yields matched cross-domain neighbors", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(707)
  dat <- generate_identity_alignment(n = 32, d = 6, noise_sd = 0.02, seed = 707)
  corr <- data.frame(
    domain_i = 1L, index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L, index_j = seq_len(nrow(dat$X2)),
    weight = 1
  )

  fit <- suppressWarnings(
    multiscale_manifold_align(
      list(d1 = dat$X1, d2 = dat$X2),
      correspondences = corr,
      preproc = multivarious::center(),
      ncomp = 3,
      control = multiscale_manifold_align_control(
        enabled = TRUE,
        backend = "randomized_dwt",
        rank_per_domain = 16L,
        knn = 8L,
        eigen_k_max = 24L,
        max_levels = 4L,
        randomized_sketch_size = 12L,
        randomized_oversample = 4L,
        randomized_power_iters = 1L,
        seed = 11L,
        verbose = FALSE
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("multiscale_manifold_align", "multiblock_biprojector"))
  expect_equal(fit$backend, "randomized_dwt")
  expect_true(all(is.finite(fit$s)))

  n <- nrow(dat$X1)
  scores_1 <- fit$s[seq_len(n), , drop = FALSE]
  scores_2 <- fit$s[n + seq_len(n), , drop = FALSE]
  cross_dist <- as.matrix(stats::dist(rbind(scores_1, scores_2)))[seq_len(n), n + seq_len(n), drop = FALSE]
  matched_mean <- mean(diag(cross_dist))
  mismatched_mean <- mean(cross_dist[row(cross_dist) != col(cross_dist)])

  expect_lt(matched_mean, mismatched_mean)
})

test_that("multiscale_manifold_align retains preprocessors for high-dimensional OOS projection", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(919)
  n <- 22L
  latent <- matrix(rnorm(n * 3L), n, 3L)
  W1 <- matrix(rnorm(3L * 160L), 3L, 160L)
  W2 <- matrix(rnorm(3L * 190L), 3L, 190L)
  X1 <- latent %*% W1 + matrix(rnorm(n * 160L, sd = 0.03), n, 160L)
  X2 <- latent %*% W2 + matrix(rnorm(n * 190L, sd = 0.03), n, 190L)
  corr <- data.frame(
    domain_i = 1L, index_i = seq_len(n),
    domain_j = 2L, index_j = seq_len(n),
    weight = 1
  )

  fit <- suppressWarnings(
    multiscale_manifold_align(
      list(d1 = X1, d2 = X2),
      correspondences = corr,
      preproc = multivarious::center(),
      ncomp = 2,
      control = multiscale_manifold_align_control(
        enabled = TRUE,
        backend = "spectral",
        rank_per_domain = 12L,
        knn = 8L,
        eigen_k_max = 20L,
        verbose = FALSE
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("multiscale_manifold_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(nrow(fit$v), ncol(X1) + ncol(X2))
  expect_s3_class(fit$preproc, "concat_pre_processor")

  new1 <- latent[1:4, , drop = FALSE] %*% W1
  new2 <- latent[1:4, , drop = FALSE] %*% W2
  pred1 <- oos_predict(fit, new1, side = 1L)
  pred2 <- oos_predict(fit, new2, side = 2L)

  expect_equal(dim(pred1), c(4L, 2L))
  expect_equal(dim(pred2), c(4L, 2L))
  expect_true(all(is.finite(pred1)))
  expect_true(all(is.finite(pred2)))
})

test_that("multiscale_manifold_align accepts unquoted label columns", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")

  set.seed(20260215)
  n <- 24
  X1 <- matrix(rnorm(n * 4), n, 4)
  X2 <- X1 + matrix(rnorm(n * 4, sd = 0.05), n, 4)
  condition_vec <- factor(rep(c("A", "B"), length.out = n))
  hd <- as_hyperdesign(
    list(domain1 = X1, domain2 = X2),
    labels = list(condition_vec, condition_vec),
    label_name = "condition"
  )

  fit <- multiscale_manifold_align(
    hd,
    y = condition,
    preproc = NULL,
    ncomp = 2,
    control = multiscale_manifold_align_control(
      enabled = TRUE,
      backend = "spectral",
      rank_per_domain = 12L,
      knn = 8L,
      max_pairs_per_label = 20L,
      verbose = FALSE
    ),
    verbose = FALSE
  )

  expect_s3_class(fit, c("multiscale_manifold_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(ncol(fit$s), 2L)
  expect_true(is.data.frame(fit$correspondences))
  expect_true(nrow(fit$correspondences) > 0L)
})

test_that("multiscale_manifold_align can tune core diffusion controls", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")

  set.seed(303)
  dat <- generate_identity_alignment(n = 24, d = 5, noise_sd = 0.03, seed = 303)
  corr <- data.frame(
    domain_i = 1L, index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L, index_j = seq_len(nrow(dat$X2)),
    weight = 1
  )

  fit <- suppressWarnings(
    multiscale_manifold_align(
      list(d1 = dat$X1, d2 = dat$X2),
      correspondences = corr,
      preproc = multivarious::center(),
      ncomp = 2,
      control = multiscale_manifold_align_control(
        enabled = TRUE,
        backend = "spectral",
        rank_per_domain = 10L,
        max_levels = 3L,
        knn = 8L,
        tune = TRUE,
        candidate_cross_edge_weight = c(0.5, 1, 2),
        candidate_rank_per_domain = c(8L, 10L),
        candidate_max_levels = c(2L, 3L),
        verbose = FALSE
      ),
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("multiscale_manifold_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_true(is.list(fit$tuning))
  expect_true(isTRUE(fit$tuning$enabled))
  expect_true(is.finite(fit$tuning$objective))
  expect_true(fit$control$cross_edge_weight %in% c(0.5, 1, 2))
  expect_true(fit$control$rank_per_domain %in% c(8L, 10L))
  expect_true(fit$control$max_levels %in% c(2L, 3L))
})
