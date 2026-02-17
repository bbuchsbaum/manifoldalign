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
