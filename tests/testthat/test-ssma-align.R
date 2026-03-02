library(testthat)

.ar1_noise <- function(n, p, phi = 0.9, sd = 1) {
  E <- matrix(rnorm(n * p, sd = sd), n, p)
  X <- matrix(0, n, p)
  X[1, ] <- E[1, ]
  for (t in 2:n) {
    X[t, ] <- phi * X[t - 1, ] + E[t, ]
  }
  X
}

test_that("ssma_align aligns near-identical domains with explicit correspondences", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(20260227)
  dat <- generate_identity_alignment(n = 36, d = 6, noise_sd = 0.03, seed = 20260227)
  corr <- data.frame(
    domain_i = 1L,
    index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L,
    index_j = seq_len(nrow(dat$X2))
  )

  fit <- ssma_align(
    list(d1 = dat$X1, d2 = dat$X2),
    correspondences = corr,
    preproc = multivarious::center(),
    ncomp = 3,
    control = ssma_align_control(
      knn = 10,
      rank_per_domain = 12,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("ssma_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_equal(ncol(fit$s), 3)

  n <- nrow(dat$X1)
  E1 <- fit$s[seq_len(n), , drop = FALSE]
  E2 <- fit$s[(n + 1):(2 * n), , drop = FALSE]
  mse <- mean(rowSums((E1 - E2)^2))

  expect_lt(mse, 0.6)
  expect_true(all(grepl("provided", fit$correspondences$source)))
})

test_that("ssma_align operator solver matches reduced-solver subspace on same fixture", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")
  skip_if_not_installed("RSpectra")

  set.seed(333)
  dat <- generate_identity_alignment(n = 44, d = 8, noise_sd = 0.04, seed = 333)
  corr <- data.frame(
    domain_i = 1L,
    index_i = seq_len(nrow(dat$X1)),
    domain_j = 2L,
    index_j = seq_len(nrow(dat$X2))
  )

  fit_reduced <- ssma_align(
    list(d1 = dat$X1, d2 = dat$X2),
    correspondences = corr,
    preproc = multivarious::center(),
    ncomp = 3,
    control = ssma_align_control(
      knn = 10,
      rank_per_domain = 16,
      solver = "reduced",
      verbose = FALSE
    )
  )

  fit_operator <- ssma_align(
    list(d1 = dat$X1, d2 = dat$X2),
    correspondences = corr,
    preproc = multivarious::center(),
    ncomp = 3,
    control = ssma_align_control(
      knn = 10,
      rank_per_domain = 16,
      solver = "operator",
      operator_tol = 1e-6,
      operator_maxiter = 3000,
      operator_power_iter = 20,
      verbose = FALSE
    )
  )

  expect_equal(fit_reduced$solver, "reduced")
  expect_equal(fit_operator$solver, "operator")
  expect_equal(ncol(fit_operator$s), 3)
  expect_true(all(is.finite(fit_operator$s)))

  Qr <- qr.Q(qr(fit_reduced$s))
  Qo <- qr.Q(qr(fit_operator$s))
  k <- min(ncol(Qr), ncol(Qo))
  sv <- svd(crossprod(Qr[, seq_len(k), drop = FALSE], Qo[, seq_len(k), drop = FALSE]), nu = 0, nv = 0)$d
  subspace_dist <- sqrt(sum(1 - pmin(sv, 1)^2))

  expect_lt(subspace_dist, 0.5)
})

test_that("ssma_align reports serial decontamination diagnostics", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("Matrix")

  set.seed(101)
  n <- 70
  z <- matrix(rnorm(n * 3), n, 3)
  X1 <- z %*% matrix(rnorm(3 * 7), 3, 7) + 0.25 * .ar1_noise(n, 7, phi = 0.95)
  X2 <- z %*% matrix(rnorm(3 * 7), 3, 7) + 0.25 * .ar1_noise(n, 7, phi = 0.92)

  corr <- data.frame(
    domain_i = 1L,
    index_i = seq_len(n),
    domain_j = 2L,
    index_j = seq_len(n)
  )

  fit <- ssma_align(
    list(a = X1, b = X2),
    correspondences = corr,
    serial_index = list(seq_len(n), seq_len(n)),
    preproc = multivarious::center(),
    ncomp = 2,
    control = ssma_align_control(
      knn = 12,
      rank_per_domain = 20,
      verbose = FALSE,
      serial_control = ssma_serial_control(
        enabled = TRUE,
        row_whiten = "ar1",
        lag_mode = "hard",
        lag_exclusion = 2L,
        preserve_degree = TRUE
      )
    )
  )

  expect_s3_class(fit, c("ssma_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_true(is.list(fit$serial))
  expect_equal(length(fit$serial), 2)
  expect_true(isTRUE(fit$serial[[1]]$enabled))
  expect_true(is.finite(fit$serial[[1]]$ar1_phi))
  expect_gt(fit$serial[[1]]$masked_edges, 0)
})

test_that("ssma_align accepts unquoted labels and design-column serial_index", {
  skip_if_not_installed("neighborweights")
  skip_if_not_installed("multivarious")

  set.seed(20260301)
  n <- 40
  X1 <- matrix(rnorm(n * 5), n, 5)
  X2 <- X1 + matrix(rnorm(n * 5, sd = 0.05), n, 5)
  cond <- factor(rep(c("A", "B"), length.out = n))

  hd <- list(
    d1 = list(x = X1, design = data.frame(condition = cond, time = seq_len(n))),
    d2 = list(x = X2, design = data.frame(condition = cond, time = seq_len(n)))
  )
  class(hd) <- "hyperdesign"

  fit <- ssma_align(
    hd,
    y = condition,
    serial_index = "time",
    preproc = multivarious::center(),
    ncomp = 2,
    control = ssma_align_control(
      knn = 8,
      rank_per_domain = 10,
      max_pairs_per_label = 25,
      verbose = FALSE,
      serial_control = ssma_serial_control(enabled = TRUE, lag_mode = "hard", lag_exclusion = 1)
    )
  )

  expect_s3_class(fit, c("ssma_align", "multiblock_biprojector"))
  expect_true(all(is.finite(fit$s)))
  expect_true(is.data.frame(fit$correspondences))
  expect_gt(nrow(fit$correspondences), 0)
  expect_true(any(grepl("label", fit$correspondences$source)))
})
