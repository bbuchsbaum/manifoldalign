# Spectral MNN: out-of-sample ridge bridge, common HKS time grid, rigid
# Procrustes with translation, RNG hygiene.

smnn_fit <- function(dat, ...) {
  spectral_mnn_align(
    list(d1 = dat$X1, d2 = dat$X2),
    preproc = multivarious::center(),
    ncomp = 3,
    control = spectral_mnn_align_control(
      knn = 8, q_descriptors = 8, max_anchors = 128, refine_rounds = 1,
      store_basis = TRUE, store_descriptors = TRUE, verbose = FALSE, ...
    )
  )
}

test_that("the out-of-sample bridge is the ridge least-squares fit to the scores", {
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.03, seed = 21)
  fit <- smnn_fit(dat, bridge_ridge = 1e-3)
  n <- nrow(dat$X1)
  p <- ncol(dat$X1)
  for (i in 1:2) {
    X <- scale(if (i == 1) dat$X1 else dat$X2, scale = FALSE)
    S <- as.matrix(fit$s)[(i - 1) * n + seq_len(n), , drop = FALSE]
    B <- fit$v[(i - 1) * p + seq_len(p), , drop = FALSE]
    lambda <- fit$bridge$lambda[[i]]
    expect_equal(lambda, 1e-3 * sum(X^2) / p)
    # Ridge normal equations: X'(X B - S) + lambda B = 0.
    grad <- crossprod(X, X %*% B - S) + lambda * B
    expect_lt(max(abs(grad)), 1e-8 * max(1, max(abs(crossprod(X, S)))))
    # oos_predict applies preprocessing and the bridge, not the training map.
    pred <- oos_predict(fit, if (i == 1) dat$X1 else dat$X2, side = i)
    expect_equal(unname(as.matrix(pred)), unname(X %*% B), tolerance = 1e-10)
  }
  expect_identical(fit$bridge$type, "linear_ridge")
})

test_that("a zero-ridge bridge is the minimum-norm least-squares solution", {
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.03, seed = 22)
  fit <- smnn_fit(dat, bridge_ridge = 0)
  X <- scale(dat$X1, scale = FALSE)
  S <- as.matrix(fit$s)[seq_len(40), , drop = FALSE]
  B <- fit$v[seq_len(6), , drop = FALSE]
  expect_equal(unname(B), unname(qr.solve(X, S)), tolerance = 1e-8)
})

test_that("training scores are the aligned eigenmaps, kept separate from the bridge", {
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.03, seed = 23)
  fit <- smnn_fit(dat)
  n <- nrow(dat$X1)
  E_ref <- scale(fit$basis[[1]]$vectors[, 1:3], center = TRUE, scale = FALSE)
  expect_equal(unname(as.matrix(fit$s)[seq_len(n), ]), E_ref, tolerance = 1e-12, ignore_attr = TRUE)
  E2 <- scale(fit$basis[[2]]$vectors[, 1:3], center = TRUE, scale = FALSE)
  aligned2 <- sweep(E2 %*% fit$rotations[[2]], 2L, fit$translations[[2]], "+")
  expect_equal(unname(as.matrix(fit$s)[n + seq_len(n), ]), aligned2, tolerance = 1e-12, ignore_attr = TRUE)
  # The bridge is an approximation of the eigenmap on 6 features: not exact.
  pred <- oos_predict(fit, dat$X1, side = 1)
  expect_gt(max(abs(pred - as.matrix(fit$s)[seq_len(n), ])), 1e-6)
})

test_that("HKS descriptors use one time grid shared by all domains", {
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.03, seed = 24)
  # Domain 2 lives at a different scale, so its spectrum differs.
  dat$X2 <- dat$X2 * 3
  fit <- smnn_fit(dat)
  pooled <- unlist(lapply(fit$basis, `[[`, "values"))
  expected <- manifoldalign:::hks_time_grid(pooled, q = 8, mode = "auto",
                                            time_range = c(0.1, 50), spacing = "log")
  expect_equal(fit$hks_times, expected)
  for (i in 1:2) {
    d_i <- compute_hks_descriptors(basis = fit$basis[[i]], times = fit$hks_times,
                                   normalize = "both")
    expect_equal(fit$descriptors[[i]], d_i)
  }
})

test_that("Procrustes estimates and applies the translation", {
  set.seed(31)
  B <- matrix(stats::rnorm(30), 10, 3)
  R0 <- qr.Q(qr(matrix(stats::rnorm(9), 3, 3)))
  if (det(R0) < 0) R0[, 3] <- -R0[, 3]
  t0 <- c(2, -1, 0.5)
  A <- sweep(B %*% R0, 2L, t0, "+")
  map <- manifoldalign:::.procrustes_rotation_spectral_mnn(A, B)
  expect_equal(map$rotation, R0, tolerance = 1e-10)
  expect_equal(map$translation, t0, tolerance = 1e-10)
  expect_equal(manifoldalign:::.apply_rigid_spectral_mnn(B, map), A, tolerance = 1e-10)
  # Weighted anchors: an exact rigid relation is still recovered.
  w <- stats::runif(10)
  map_w <- manifoldalign:::.procrustes_rotation_spectral_mnn(A, B, weights = w)
  expect_equal(manifoldalign:::.apply_rigid_spectral_mnn(B, map_w), A, tolerance = 1e-10)
})

test_that("spectral_mnn_align does not disturb the caller's RNG stream", {
  dat <- generate_identity_alignment(n = 40, d = 6, noise_sd = 0.03, seed = 25)
  set.seed(17)
  before <- .Random.seed
  fit <- spectral_mnn_align(
    list(d1 = dat$X1, d2 = dat$X2), ncomp = 2,
    control = spectral_mnn_align_control(knn = 8, verbose = FALSE)
  )
  expect_identical(.Random.seed, before)
})
