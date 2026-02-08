test_that("global_geo_align fits from explicit correspondences", {
  set.seed(42)
  n <- 36
  z <- matrix(rnorm(n * 3), n, 3)

  W1 <- matrix(rnorm(3 * 8), 3, 8)
  W2 <- matrix(rnorm(3 * 7), 3, 7)
  W3 <- matrix(rnorm(3 * 6), 3, 6)

  X1 <- z %*% W1 + matrix(rnorm(n * 8, sd = 0.05), n, 8)
  X2 <- z %*% W2 + matrix(rnorm(n * 7, sd = 0.05), n, 7)
  X3 <- z %*% W3 + matrix(rnorm(n * 6, sd = 0.05), n, 6)

  corr <- rbind(
    data.frame(domain_i = 1L, index_i = seq_len(n), domain_j = 2L, index_j = seq_len(n)),
    data.frame(domain_i = 1L, index_i = seq_len(n), domain_j = 3L, index_j = seq_len(n))
  )

  fit <- global_geo_align(
    list(d1 = X1, d2 = X2, d3 = X3),
    correspondences = corr,
    preproc = NULL,
    ncomp = 3,
    control = global_geo_align_control(
      knn = 8,
      n_landmarks_total = 30,
      k_embed = 10,
      ridge_eps = 1e-4,
      verbose = FALSE,
      seed = 11
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$s), 3 * n)
  expect_equal(ncol(fit$s), 3)
  expect_equal(nrow(fit$v), ncol(X1) + ncol(X2) + ncol(X3))
  expect_true(all(is.finite(fit$s)))
  expect_true(length(fit$landmarks) > 0)
  expect_equal(nrow(fit$correspondences), 2 * n)
})

test_that("global_geo_align builds correspondences from labels", {
  set.seed(101)
  n <- 30
  labs <- rep(c("A", "B", "C"), each = n / 3)
  z <- matrix(rnorm(n * 2), n, 2)

  X1 <- z %*% matrix(rnorm(2 * 6), 2, 6) + matrix(rnorm(n * 6, sd = 0.03), n, 6)
  X2 <- z %*% matrix(rnorm(2 * 5), 2, 5) + matrix(rnorm(n * 5, sd = 0.03), n, 5)
  X3 <- z %*% matrix(rnorm(2 * 7), 2, 7) + matrix(rnorm(n * 7, sd = 0.03), n, 7)

  hd <- create_multi_alignment_hyperdesign(list(X1, X2, X3))
  for (i in seq_along(hd)) {
    hd[[i]]$design$lbl <- labs
  }

  fit <- global_geo_align(
    hd,
    y = lbl,
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      knn = 6,
      n_landmarks_total = 24,
      k_embed = 8,
      max_pairs_per_label = 8,
      verbose = FALSE,
      seed = 4
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$s), 3 * n)
  expect_equal(ncol(fit$s), 2)
  expect_true(nrow(fit$correspondences) > 0)
  expect_true(all(fit$domain_scales > 0))
})

test_that("global_geo_aligner dispatches native fit_many", {
  set.seed(7)
  X1 <- matrix(rnorm(24 * 4), 24, 4)
  X2 <- X1 + matrix(rnorm(24 * 4, sd = 0.02), 24, 4)
  corr <- data.frame(domain_i = 1L, index_i = 1:24, domain_j = 2L, index_j = 1:24)

  algo <- global_geo_aligner()
  fit <- fit_many(
    algo,
    list(X1, X2),
    correspondences = corr,
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      knn = 5,
      n_landmarks_total = 16,
      k_embed = 6,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_equal(ncol(fit$s), 2)
})

test_that("global_geo_align supports spatial anchor feature correspondence", {
  set.seed(21)
  n <- 28
  z <- matrix(rnorm(n * 3), n, 3)
  W1 <- matrix(rnorm(3 * 9), 3, 9)
  W2 <- matrix(rnorm(3 * 7), 3, 7)
  X1 <- z %*% W1 + matrix(rnorm(n * 9, sd = 0.03), n, 9)
  X2 <- z %*% W2 + matrix(rnorm(n * 7, sd = 0.03), n, 7)

  coords1 <- cbind(seq_len(9), rep(0, 9), rep(0, 9))
  coords2 <- cbind(seq_len(7), rep(0, 7), rep(0, 7))
  corr <- data.frame(domain_i = 1L, index_i = seq_len(n), domain_j = 2L, index_j = seq_len(n))

  fit <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = list(
      voxel_coords = list(d1 = coords1, d2 = coords2),
      num_anchors = 6,
      nn_anchors = 3,
      sigma = 1.5,
      anchor_method = "kmeans",
      transform = "identity",
      store_B = TRUE,
      seed = 9
    ),
    preproc = NULL,
    ncomp = 3,
    control = global_geo_align_control(
      knn = 6,
      n_landmarks_total = 20,
      k_embed = 10,
      verbose = FALSE,
      seed = 9
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_equal(nrow(fit$v), ncol(X1) + ncol(X2))
  expect_equal(ncol(fit$s), 3)
  expect_true(!is.null(fit$feature_correspondence))
  expect_equal(nrow(fit$feature_correspondence$anchors), 6)
  expect_equal(length(fit$feature_correspondence$B), 2)
})

test_that("global_geo_align supports anchor smoothness regularization", {
  set.seed(77)
  n <- 26
  z <- matrix(rnorm(n * 2), n, 2)
  X1 <- z %*% matrix(rnorm(2 * 8), 2, 8) + matrix(rnorm(n * 8, sd = 0.03), n, 8)
  X2 <- z %*% matrix(rnorm(2 * 8), 2, 8) + matrix(rnorm(n * 8, sd = 0.03), n, 8)
  corr <- data.frame(domain_i = 1L, index_i = 1:n, domain_j = 2L, index_j = 1:n)
  coords <- cbind(seq_len(8), rep(0, 8), rep(0, 8))

  fit <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = list(
      voxel_coords = list(d1 = coords, d2 = coords),
      num_anchors = 6,
      nn_anchors = 3,
      sigma = 1.2,
      smoothness = list(
        enabled = TRUE,
        lambda = 0.2,
        knn = 3,
        weight_mode = "heat",
        store_laplacian = TRUE
      )
    ),
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      knn = 5,
      n_landmarks_total = 18,
      k_embed = 8,
      verbose = FALSE,
      seed = 5
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_true(!is.null(fit$feature_correspondence$smoothness))
  expect_equal(fit$feature_correspondence$smoothness$lambda, 0.2)
  expect_true(methods::is(fit$feature_correspondence$smoothness$laplacian, "Matrix"))
})

test_that("global_geo_align validates SAB voxel coordinate dimensions", {
  set.seed(9)
  X1 <- matrix(rnorm(20 * 5), 20, 5)
  X2 <- matrix(rnorm(20 * 4), 20, 4)
  corr <- data.frame(domain_i = 1L, index_i = 1:20, domain_j = 2L, index_j = 1:20)

  expect_error(
    global_geo_align(
      list(d1 = X1, d2 = X2),
      correspondences = corr,
      feature_correspondence = list(
        voxel_coords = list(
          d1 = matrix(rnorm(4 * 3), 4, 3), # wrong: should be 5 rows
          d2 = matrix(rnorm(4 * 3), 4, 3)
        ),
        num_anchors = 4,
        nn_anchors = 2
      ),
      preproc = NULL,
      ncomp = 2,
      control = global_geo_align_control(
        n_landmarks_total = 12,
        k_embed = 6,
        verbose = FALSE
      )
    ),
    "Voxel coordinate rows must match number of features"
  )
})

test_that("global_geo_align ignores SAB when explicitly disabled", {
  set.seed(12)
  X1 <- matrix(rnorm(18 * 5), 18, 5)
  X2 <- X1 + matrix(rnorm(18 * 5, sd = 0.02), 18, 5)
  corr <- data.frame(domain_i = 1L, index_i = 1:18, domain_j = 2L, index_j = 1:18)

  fit <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = list(enabled = FALSE),
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      n_landmarks_total = 12,
      k_embed = 6,
      verbose = FALSE
    )
  )

  expect_s3_class(fit, c("global_geo_align", "multiblock_biprojector"))
  expect_true(is.null(fit$feature_correspondence))
})

test_that("SAB transfer matrices are row-stochastic and sparse", {
  set.seed(123)
  n <- 24
  p1 <- 15
  p2 <- 12
  z <- matrix(rnorm(n * 2), n, 2)
  X1 <- z %*% matrix(rnorm(2 * p1), 2, p1) + matrix(rnorm(n * p1, sd = 0.05), n, p1)
  X2 <- z %*% matrix(rnorm(2 * p2), 2, p2) + matrix(rnorm(n * p2, sd = 0.05), n, p2)
  corr <- data.frame(domain_i = 1L, index_i = 1:n, domain_j = 2L, index_j = 1:n)

  coords1 <- cbind(runif(p1), runif(p1), runif(p1))
  coords2 <- cbind(runif(p2), runif(p2), runif(p2))
  nn_anchors <- 4L

  fit <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = list(
      voxel_coords = list(d1 = coords1, d2 = coords2),
      num_anchors = 10,
      nn_anchors = nn_anchors,
      sigma = 0.25,
      anchor_method = "kmeans",
      store_B = TRUE,
      seed = 13
    ),
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      knn = 6,
      n_landmarks_total = 16,
      k_embed = 8,
      verbose = FALSE,
      seed = 13
    )
  )

  B1 <- fit$feature_correspondence$B[[1]]
  B2 <- fit$feature_correspondence$B[[2]]
  expect_true(methods::is(B1, "sparseMatrix"))
  expect_true(methods::is(B2, "sparseMatrix"))

  rs1 <- as.numeric(Matrix::rowSums(B1))
  rs2 <- as.numeric(Matrix::rowSums(B2))
  expect_lt(max(abs(rs1 - 1)), 1e-8)
  expect_lt(max(abs(rs2 - 1)), 1e-8)

  nz1 <- Matrix::rowSums(B1 != 0)
  nz2 <- Matrix::rowSums(B2 != 0)
  expect_true(all(nz1 <= nn_anchors))
  expect_true(all(nz2 <= nn_anchors))
})

test_that("SAB gives identical transfer matrices for identical coordinates", {
  set.seed(55)
  n <- 20
  p <- 11
  z <- matrix(rnorm(n * 2), n, 2)
  X1 <- z %*% matrix(rnorm(2 * p), 2, p) + matrix(rnorm(n * p, sd = 0.04), n, p)
  X2 <- z %*% matrix(rnorm(2 * p), 2, p) + matrix(rnorm(n * p, sd = 0.04), n, p)
  corr <- data.frame(domain_i = 1L, index_i = 1:n, domain_j = 2L, index_j = 1:n)

  coords <- cbind(seq_len(p), rep(0, p), rep(0, p))
  fit <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = list(
      voxel_coords = list(d1 = coords, d2 = coords),
      num_anchors = 7,
      nn_anchors = 3,
      sigma = 1.0,
      anchor_method = "grid",
      reference_dataset = "d1",
      transform = "identity",
      store_B = TRUE
    ),
    preproc = NULL,
    ncomp = 2,
    control = global_geo_align_control(
      knn = 5,
      n_landmarks_total = 14,
      k_embed = 6,
      verbose = FALSE,
      seed = 2
    )
  )

  B1 <- as.matrix(fit$feature_correspondence$B[[1]])
  B2 <- as.matrix(fit$feature_correspondence$B[[2]])
  expect_equal(B1, B2, tolerance = 1e-10)
})

test_that("SAB smoothness lowers anchor roughness on synthetic smooth signal", {
  set.seed(808)
  n <- 40
  p <- 24
  ncomp <- 1
  corr <- data.frame(domain_i = 1L, index_i = 1:n, domain_j = 2L, index_j = 1:n)

  # Smooth latent signal across voxel coordinates in 1D (embedded in 3D).
  coords <- cbind(seq(0, 1, length.out = p), rep(0, p), rep(0, p))
  w_true <- exp(-((coords[, 1] - 0.5)^2) / (2 * 0.12^2))
  w_true <- w_true / sqrt(sum(w_true^2))
  z <- matrix(rnorm(n), n, 1)
  X1 <- z %*% t(w_true) + matrix(rnorm(n * p, sd = 0.20), n, p)
  X2 <- z %*% t(w_true) + matrix(rnorm(n * p, sd = 0.20), n, p)

  fc_base <- list(
    voxel_coords = list(d1 = coords, d2 = coords),
    num_anchors = 10,
    nn_anchors = 3,
    sigma = 0.10,
    anchor_method = "grid",
    transform = "identity",
    store_B = TRUE
  )

  fit_no <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = fc_base,
    preproc = NULL,
    ncomp = ncomp,
    control = global_geo_align_control(
      knn = 6,
      n_landmarks_total = 20,
      k_embed = 8,
      verbose = FALSE,
      seed = 17
    )
  )

  fit_sm <- global_geo_align(
    list(d1 = X1, d2 = X2),
    correspondences = corr,
    feature_correspondence = modifyList(
      fc_base,
      list(smoothness = list(
        enabled = TRUE,
        lambda = 2.0,
        knn = 3,
        weight_mode = "heat",
        store_laplacian = TRUE
      ))
    ),
    preproc = NULL,
    ncomp = ncomp,
    control = global_geo_align_control(
      knn = 6,
      n_landmarks_total = 20,
      k_embed = 8,
      verbose = FALSE,
      seed = 17
    )
  )

  K <- nrow(fit_sm$feature_correspondence$anchors)
  L <- as.matrix(fit_sm$feature_correspondence$smoothness$laplacian)

  split_anchor <- function(Astack) {
    list(
      Astack[1:K, , drop = FALSE],
      Astack[(K + 1):(2 * K), , drop = FALSE]
    )
  }
  norm_cols <- function(A) {
    den <- sqrt(colSums(A^2))
    den[den <= .Machine$double.eps] <- 1
    sweep(A, 2, den, "/")
  }
  roughness <- function(Ablocks) {
    sum(vapply(Ablocks, function(A) {
      An <- norm_cols(A)
      sum(diag(crossprod(An, L %*% An)))
    }, numeric(1)))
  }

  r_no <- roughness(split_anchor(fit_no$feature_correspondence$anchor_loadings))
  r_sm <- roughness(split_anchor(fit_sm$feature_correspondence$anchor_loadings))
  expect_lt(r_sm, r_no)
})
