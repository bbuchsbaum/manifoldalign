# Contracts of global_geo_align(): per-domain preprocessing lists, the
# concatenated preprocessor, and RNG hygiene.

gga_fixture <- function(seed = 42, n = 30) {
  set.seed(seed)
  z <- matrix(rnorm(n * 3), n, 3)
  X <- list(
    d1 = z %*% matrix(rnorm(3 * 8), 3, 8) + matrix(rnorm(n * 8, sd = 0.05), n, 8) + 2,
    d2 = z %*% matrix(rnorm(3 * 7), 3, 7) + matrix(rnorm(n * 7, sd = 0.05), n, 7) - 1,
    d3 = z %*% matrix(rnorm(3 * 6), 3, 6) + matrix(rnorm(n * 6, sd = 0.05), n, 6)
  )
  corr <- rbind(
    data.frame(domain_i = 1L, index_i = seq_len(n), domain_j = 2L, index_j = seq_len(n)),
    data.frame(domain_i = 1L, index_i = seq_len(n), domain_j = 3L, index_j = seq_len(n))
  )
  list(X = X, corr = corr)
}

gga_control <- function() {
  global_geo_align_control(knn = 8, n_landmarks_total = 24, k_embed = 10,
                           ridge_eps = 1e-4, verbose = FALSE, seed = 11)
}

test_that("global_geo_align keeps one preprocessor per domain for mixed lists", {
  fx <- gga_fixture()
  # Regression: `proclist[[i]] <- NULL` deleted list entries, so a mixed list
  # shifted the fitted preprocessors onto the wrong domains.
  fit <- global_geo_align(fx$X, correspondences = fx$corr,
                          preproc = list(NULL, multivarious::center(), NULL),
                          ncomp = 3, control = gga_control())
  expect_true(all(is.finite(fit$s)))

  # The concatenated preprocessor acts on feature blocks (8 + 7 + 6 columns).
  Xnew <- do.call(cbind, lapply(fx$X, function(x) x[1:4, , drop = FALSE]))
  tf <- multivarious::transform(fit$preproc, Xnew)
  expect_equal(dim(tf), dim(Xnew))
  expect_equal(unname(tf[, 1:8]), unname(Xnew[, 1:8]))
  expect_equal(unname(tf[, 9:15]),
               unname(sweep(Xnew[, 9:15], 2, colMeans(fx$X$d2))), tolerance = 1e-12)
  expect_equal(unname(tf[, 16:21]), unname(Xnew[, 16:21]))
})

test_that("broadcast preprocessing yields a feature-block preprocessor", {
  fx <- gga_fixture()
  fit <- global_geo_align(fx$X, correspondences = fx$corr,
                          preproc = multivarious::center(),
                          ncomp = 3, control = gga_control())
  Xnew <- do.call(cbind, lapply(fx$X, function(x) x[1:4, , drop = FALSE]))
  # Regression: blocks were indexed by sample rows (90) instead of features (21).
  tf <- multivarious::transform(fit$preproc, Xnew)
  centred <- do.call(cbind, lapply(fx$X, function(x) {
    sweep(x[1:4, , drop = FALSE], 2, colMeans(x))
  }))
  expect_equal(unname(tf), unname(centred), tolerance = 1e-12)

  # All-NULL preprocessing still records no preprocessor.
  fit_null <- global_geo_align(fx$X, correspondences = fx$corr, preproc = NULL,
                               ncomp = 3, control = gga_control())
  expect_null(fit_null$preproc)
})

test_that("global_geo_align leaves the caller's RNG stream untouched", {
  fx <- gga_fixture(n = 24)
  labs <- rep(c("A", "B", "C"), each = 8)
  hd <- create_multi_alignment_hyperdesign(unname(fx$X))
  for (i in seq_along(hd)) hd[[i]]$design$lbl <- labs
  coords <- list(cbind(seq_len(8), 0, 0), cbind(seq_len(7), 0, 0), cbind(seq_len(6), 0, 0))
  names(coords) <- names(hd)

  run <- function() {
    global_geo_align(
      hd, y = lbl, preproc = NULL, ncomp = 2,
      feature_correspondence = list(coords = coords, num_anchors = 4, nn_anchors = 2,
                                    anchor_method = "kmeans", seed = 5),
      control = global_geo_align_control(knn = 6, n_landmarks_total = 20, k_embed = 8,
                                         max_pairs_per_label = 4, verbose = FALSE,
                                         seed = 4)
    )
  }

  set.seed(123)
  expected <- runif(3)
  set.seed(123)
  fit1 <- run()
  expect_equal(runif(3), expected)

  fit2 <- run()
  expect_equal(fit1$s, fit2$s)
  expect_equal(fit1$correspondences, fit2$correspondences)
})
