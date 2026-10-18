# KEMA prediction contract: a fit is a per-domain kernel encoder
# f_i(x) = k(phi_i(x), R_i) alpha_i, with R_i every training row (full KEMA)
# or the stored landmark rows (REKEMA).

kema_encoder_hd <- function(n = 30, p = 4, seed = 11) {
  set.seed(seed)
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.5, n, p)
  y <- rep(1:3, length.out = n)
  md <- list(
    d1 = multidesign::multidesign(X1, data.frame(lbl = factor(y))),
    d2 = multidesign::multidesign(X2, data.frame(lbl = factor(y)))
  )
  list(hd = multidesign::hyperdesign(md), X1 = X1, X2 = X2, n = n, p = p)
}

test_that("re-encoding training rows reproduces full RBF KEMA scores", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd()
  fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
              kernel = kernlab::rbfdot(sigma = 0.5))

  expect_s3_class(fit, "kema_fit")
  expect_identical(fit$encoder$form, "full")
  s <- as.matrix(fit$s)
  pred <- rbind(
    oos_predict(fit, dat$X1, side = "d1"),
    oos_predict(fit, dat$X2, side = 2)
  )
  expect_equal(dim(pred), dim(s))
  expect_lt(max(abs(pred - s)), 1e-10 * max(1, max(abs(s))))

  # Independent oracle: k(X_new, X_train) alpha computed by hand.
  Xc <- scale(dat$X1, center = TRUE, scale = FALSE)
  K <- exp(-0.5 * as.matrix(dist(Xc))^2)
  alpha1 <- as.matrix(fit$alpha)[seq_len(dat$n), , drop = FALSE]
  expect_lt(max(abs(K %*% alpha1 - s[seq_len(dat$n), ])), 1e-10)

  # The former linear surrogate X v cannot stand in for the kernel map.
  expect_true(all(is.na(fit$v)))
  expect_identical(fit$loadings_status, "unavailable_nonlinear_kernel")
})

test_that("reduced linear KEMA reproduces its own training scores", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd(n = 40, seed = 5)
  fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
              kernel = kernlab::vanilladot(), sample_frac = 0.5, seed = 3)

  expect_identical(fit$encoder$form, "reduced")
  s <- as.matrix(fit$s)
  pred <- rbind(
    oos_predict(fit, dat$X1, side = 1),
    oos_predict(fit, dat$X2, side = 2)
  )
  expect_lt(max(abs(pred - s)), 1e-10 * max(1, max(abs(s))))

  # For the linear kernel the loadings are exact: X_c v = scores.
  X1c <- scale(dat$X1, center = TRUE, scale = FALSE)
  expect_identical(fit$loadings_status, "exact_linear_kernel")
  expect_lt(max(abs(X1c %*% fit$v[seq_len(dat$p), ] - s[seq_len(dat$n), ])), 1e-10)
})

test_that("KEMA prediction of a row is independent of the rest of the batch", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd()
  fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
              kernel = kernlab::rbfdot(sigma = 0.5))
  set.seed(99)
  newX <- matrix(rnorm(8 * dat$p), 8, dat$p)
  batch <- oos_predict(fit, newX, side = 2)
  single <- t(vapply(seq_len(nrow(newX)), function(i) {
    oos_predict(fit, newX[i, , drop = FALSE], side = 2)[1, ]
  }, numeric(2)))
  expect_equal(batch, single, tolerance = 1e-12)
  shuffled <- oos_predict(fit, newX[8:1, ], side = 2)
  expect_equal(shuffled[8:1, ], batch, tolerance = 1e-12)
  # A numeric vector is a single row.
  expect_equal(oos_predict(fit, newX[3, ], side = 2)[1, ], batch[3, ], tolerance = 1e-12)
})

test_that("REKEMA landmarks are seeded, stored and leave the caller RNG alone", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd(n = 40, seed = 6)
  set.seed(2024)
  before <- .Random.seed
  f1 <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
             kernel = kernlab::rbfdot(sigma = 0.5), sample_frac = 0.5, seed = 17)
  expect_identical(.Random.seed, before)
  f2 <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
             kernel = kernlab::rbfdot(sigma = 0.5), sample_frac = 0.5, seed = 17)
  expect_identical(f1$encoder$landmarks, f2$encoder$landmarks)
  expect_identical(lengths(f1$encoder$landmarks), c(d1 = 20L, d2 = 20L))
  expect_equal(f1$encoder$reference$d1, scale(dat$X1, scale = FALSE)[f1$encoder$landmarks$d1, ],
               ignore_attr = TRUE)
  f3 <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
             kernel = kernlab::rbfdot(sigma = 0.5), sample_frac = 0.5, seed = 18)
  expect_false(identical(f1$encoder$landmarks, f3$encoder$landmarks))
  # The stored landmarks reproduce the training scores.
  pred <- oos_predict(f1, dat$X1, side = 1)
  expect_lt(max(abs(pred - as.matrix(f1$s)[seq_len(dat$n), ])), 1e-10)
})

test_that("generic linear projection of a KEMA fit is an explicit error", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd()
  fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
              kernel = kernlab::rbfdot(sigma = 0.5))
  expect_error(multivarious::project(fit, cbind(dat$X1, dat$X2)),
               "not defined for KEMA")
  expect_error(multivarious::partial_project(fit, dat$X1, colind = 1:4),
               "not defined for KEMA")
  expect_equal(
    multivarious::project_block(fit, dat$X2, block = 2),
    oos_predict(fit, dat$X2, side = 2)
  )
  expect_error(oos_predict(fit, dat$X1[, 1:3], side = 1), "columns")
  expect_error(oos_predict(fit, dat$X1, side = "nope"), "Unknown KEMA domain")
  expect_error(oos_predict(fit, dat$X1, side = 1, extra = 1), "no additional")
})

test_that("KEMA encoder honours NULL and function preprocessing", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd()
  for (pp in list(NULL, function(x) 2 * x)) {
    fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
                kernel = kernlab::rbfdot(sigma = 0.5), preproc = pp)
    pred <- oos_predict(fit, dat$X1, side = 1)
    expect_lt(max(abs(pred - as.matrix(fit$s)[seq_len(dat$n), ])), 1e-10)
  }
})

test_that("every KEMA backend returns a consistent encoder", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd(n = 24, seed = 8)
  for (backend in c("full_exact", "reduced_exact", "operator_exact")) {
    fit <- kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
                kernel = kernlab::rbfdot(sigma = 0.5), backend = backend)
    pred <- rbind(oos_predict(fit, dat$X1, side = 1),
                  oos_predict(fit, dat$X2, side = 2))
    expect_lt(max(abs(pred - as.matrix(fit$s))), 1e-10, label = backend)
  }
})

test_that("KEMA rejects problems whose dense footprint exceeds the guard", {
  skip_if_not_installed("multidesign")
  dat <- kema_encoder_hd()
  expect_error(
    kema(dat$hd, y = lbl, ncomp = 2, knn = 4,
         backend_control = list(max_dense_bytes = 1000)),
    "dense O\\(n\\^2\\).*max_dense_bytes"
  )
  est <- manifoldalign:::.kema_dense_memory_estimate(c(20000, 20000), c(20000, 20000), 0)
  expect_gt(est, 8 * 1024^3)
  expect_lt(
    manifoldalign:::.kema_dense_memory_estimate(c(100, 100), c(100, 100), 200),
    8 * 1024^3
  )
})

test_that("block_indices counts rows of sparse matrices", {
  blocks <- list(
    list(x = Matrix::Matrix(matrix(1, 3, 2), sparse = TRUE)),
    list(x = matrix(1, 4, 2))
  )
  bi <- block_indices(blocks)
  expect_equal(unname(bi[, "end"]), c(3, 7))
})
