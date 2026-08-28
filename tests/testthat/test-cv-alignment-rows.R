library(testthat)

.make_cv_alignment_fixture <- function(seed = 20260410, n = 18, latent_dim = 3, p = 6) {
  set.seed(seed)

  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  A1 <- matrix(rnorm(latent_dim * p), latent_dim, p)
  A2 <- matrix(rnorm(latent_dim * p), latent_dim, p)

  X1 <- Z %*% A1 + matrix(rnorm(n * p, sd = 0.03), n, p)
  X2 <- Z %*% A2 + matrix(rnorm(n * p, sd = 0.03), n, p)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(d1 = md1, d2 = md2), block_names = c("d1", "d2"))

  features <- list(
    d1 = Z + matrix(rnorm(n * latent_dim, sd = 0.01), n, latent_dim),
    d2 = Z + matrix(rnorm(n * latent_dim, sd = 0.01), n, latent_dim)
  )

  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    ids1 <- resolved$domains[[1]]$design$row_id
    ids2 <- resolved$domains[[2]]$design$row_id
    shared <- intersect(ids1, ids2)

    corr <- data.frame(
      domain_i = rep(1L, length(shared)),
      index_i = match(shared, ids1),
      domain_j = rep(2L, length(shared)),
      index_j = match(shared, ids2)
    )

    ssma_align(
      analysis,
      correspondences = corr,
      preproc = multivarious::center(),
      ncomp = 2,
      control = ssma_align_control(
        knn = 8,
        rank_per_domain = 10,
        verbose = FALSE
      )
    )
  }

  list(hd = hd, features = features, fit_fn = fit_fn, Z = Z)
}

.make_dual_mode_fixture <- function(seed = 20260411, n = 16, latent_dim = 3) {
  set.seed(seed)

  Za <- matrix(rnorm(n * latent_dim), n, latent_dim)
  Zb <- matrix(rnorm(n * latent_dim), n, latent_dim)

  md1 <- multidesign::multidesign(Za, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(Zb, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    Xa <- as.matrix(resolved$domains[[1]]$x)
    Xb <- as.matrix(resolved$domains[[2]]$x)
    d <- ncol(Xa)

    structure(
      list(
        s = rbind(Xa, Xb),
        v = rbind(diag(d), diag(d)),
        sdev = apply(rbind(Xa, Xb), 2, stats::sd),
        preproc = NULL,
        block_indices = list(a = seq_len(d), b = d + seq_len(d)),
        domain_names = resolved$domain_names,
        targets = list(a = Xa, b = Xb)
      ),
      class = c("mock_dual_fit", "multiblock_biprojector")
    )
  }

  list(
    hd = hd,
    features = list(a = Za, b = Zb),
    fit_fn = fit_fn,
    Za = Za,
    Zb = Zb
  )
}

test_that("cv_alignment_rows scores synchronized SSMA row holdouts", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("adjoin")

  fixture <- .make_cv_alignment_fixture(seed = 20260402)

  rows <- list(
    list(d1 = 1:3, d2 = 1:3),
    list(d1 = 4:6, d2 = 4:6),
    list(d1 = 7:9, d2 = 7:9)
  )

  res <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 2,
    target_pool = "assessment"
  )

  expect_s3_class(res, "cv_result")
  expect_equal(nrow(res$scores), length(rows))
  expect_true(all(c(
    ".fold",
    "mean_top1_similarity",
    "mean_topk_similarity",
    "oracle_top1_similarity",
    "oracle_topk_similarity",
    "top1_gap",
    "topk_gap",
    "n_queries",
    "n_pairs"
  ) %in% names(res$scores)))
  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
  expect_true(all(is.finite(res$scores$mean_topk_similarity)))
  expect_true(all(res$scores$oracle_top1_similarity >= res$scores$mean_top1_similarity))
  expect_true(all(res$scores$oracle_topk_similarity >= res$scores$mean_topk_similarity))
  expect_gt(mean(res$scores$mean_top1_similarity), 0.7)
  expect_gt(mean(res$scores$mean_topk_similarity), 0.45)
})

test_that("cv_alignment_rows supports one-sided holdouts against analysis targets", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("adjoin")

  fixture <- .make_cv_alignment_fixture(seed = 20260403)

  rows <- list(
    list(d1 = 1:3),
    list(d1 = 4:6),
    list(d1 = 7:9)
  )

  res <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 2,
    target_pool = "analysis"
  )

  expect_s3_class(res, "cv_result")
  expect_equal(nrow(res$scores), length(rows))
  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
  expect_true(all(res$scores$n_pairs == 1))
  expect_gt(mean(res$scores$mean_top1_similarity), 0.6)
})

test_that("cv_alignment_rows treats target_pool = both like analysis for one-sided folds", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("adjoin")

  fixture <- .make_cv_alignment_fixture(seed = 20260405)
  rows <- list(
    list(d1 = 1:3),
    list(d1 = 4:6)
  )

  res_analysis <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 2,
    target_pool = "analysis"
  )

  res_both <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 2,
    target_pool = "both"
  )

  expect_equal(
    res_both$scores$mean_top1_similarity,
    res_analysis$scores$mean_top1_similarity,
    tolerance = 1e-10
  )
  expect_equal(
    res_both$scores$mean_topk_similarity,
    res_analysis$scores$mean_topk_similarity,
    tolerance = 1e-10
  )
  expect_equal(res_both$scores$n_pairs, res_analysis$scores$n_pairs)
})

test_that("cv_alignment_rows supports weight-based prediction dispatch", {
  skip_if_not_installed("multidesign")

  predict.mock_weight_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_weight_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_weight_fit", predict.mock_weight_fit)

  set.seed(20260404)
  n <- 14
  latent_dim <- 2

  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  X1 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)
  X2 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    structure(
      list(
        domain_names = resolved$domain_names,
        targets = lapply(resolved$domains, function(dom) as.matrix(dom$x))
      ),
      class = "mock_weight_fit"
    )
  }

  res <- cv_alignment_rows(
    hd,
    rows = list(list(a = 1:4), list(a = 5:8)),
    fit_fn = fit_fn,
    features = list(a = Z, b = Z),
    k = 2,
    target_pool = "analysis",
    prediction_mode = "weights"
  )

  expect_s3_class(res, "cv_result")
  expect_equal(nrow(res$scores), 2)
  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
  expect_true(all(res$scores$n_pairs == 1))
  expect_gt(mean(res$scores$mean_top1_similarity), 0.6)
})

test_that("cv_alignment_rows embedding and weight modes agree when rankings match", {
  skip_if_not_installed("multidesign")

  predict.mock_dual_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_dual_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_dual_fit", predict.mock_dual_fit)

  fixture <- .make_dual_mode_fixture(seed = 20260408)
  rows <- list(
    list(a = 1:3),
    list(a = 4:6)
  )

  res_embedding <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  res_weights <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "weights"
  )

  expect_equal(
    res_embedding$scores$mean_top1_similarity,
    res_weights$scores$mean_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_embedding$scores$mean_topk_similarity,
    res_weights$scores$mean_topk_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_embedding$scores$oracle_top1_similarity,
    res_weights$scores$oracle_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_embedding$scores$oracle_topk_similarity,
    res_weights$scores$oracle_topk_similarity,
    tolerance = 1e-12
  )
})

test_that("cv_alignment_rows auto mode prefers weights for transport-style fits", {
  skip_if_not_installed("multidesign")

  predict.mock_transport_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_transport_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_transport_fit", predict.mock_transport_fit)

  set.seed(20260406)
  n <- 12
  latent_dim <- 2
  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  X1 <- Z + matrix(rnorm(n * latent_dim, sd = 0.02), n, latent_dim)
  X2 <- Z + matrix(rnorm(n * latent_dim, sd = 0.02), n, latent_dim)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    structure(
      list(
        domain_names = resolved$domain_names,
        targets = lapply(resolved$domains, function(dom) as.matrix(dom$x))
      ),
      class = c("mock_transport_fit", "fpgw")
    )
  }

  res <- cv_alignment_rows(
    hd,
    rows = list(list(a = 1:3), list(a = 4:6)),
    fit_fn = fit_fn,
    features = list(a = Z, b = Z),
    k = 2,
    target_pool = "analysis",
    prediction_mode = "auto"
  )

  expect_s3_class(res, "cv_result")
  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
})

test_that("cv_alignment_rows auto and weights agree for real gromov_wasserstein fits", {
  skip_if_not_installed("multidesign")

  set.seed(20260415)
  n <- 10
  latent_dim <- 2
  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  X1 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)
  X2 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))
  rows <- list(list(a = 1:2), list(a = 3:4))
  fit_fn <- function(analysis) {
    gromov_wasserstein(analysis, epsilon = 0.05, max_iter = 50, verbose = FALSE)
  }

  res_auto <- cv_alignment_rows(
    hd,
    rows = rows,
    fit_fn = fit_fn,
    features = list(a = Z, b = Z),
    k = 2,
    target_pool = "analysis",
    prediction_mode = "auto"
  )

  res_weights <- cv_alignment_rows(
    hd,
    rows = rows,
    fit_fn = fit_fn,
    features = list(a = Z, b = Z),
    k = 2,
    target_pool = "analysis",
    prediction_mode = "weights"
  )

  expect_equal(
    res_auto$scores$mean_top1_similarity,
    res_weights$scores$mean_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_auto$scores$mean_topk_similarity,
    res_weights$scores$mean_topk_similarity,
    tolerance = 1e-12
  )
  expect_true(all(is.finite(res_auto$scores$mean_top1_similarity)))
})

test_that("cv_alignment_rows auto falls back to embedding when weights fail", {
  skip_if_not_installed("multidesign")

  predict.mock_fallback_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (identical(type, "weights")) {
      stop("synthetic weight failure", call. = FALSE)
    }
    stop("mock_fallback_fit only implements failing weights", call. = FALSE)
  }

  registerS3method("predict", "mock_fallback_fit", predict.mock_fallback_fit)

  fixture <- .make_dual_mode_fixture(seed = 20260417)
  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    Xa <- as.matrix(resolved$domains[[1]]$x)
    Xb <- as.matrix(resolved$domains[[2]]$x)
    d <- ncol(Xa)

    structure(
      list(
        s = rbind(Xa, Xb),
        v = rbind(diag(d), diag(d)),
        sdev = apply(rbind(Xa, Xb), 2, stats::sd),
        preproc = NULL,
        block_indices = list(a = seq_len(d), b = d + seq_len(d)),
        domain_names = resolved$domain_names,
        targets = list(a = Xa, b = Xb)
      ),
      class = c("mock_fallback_fit", "fpgw", "multiblock_biprojector")
    )
  }

  rows <- list(
    list(a = 1:3),
    list(a = 4:6)
  )

  res_embedding <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  res_auto <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "auto"
  )

  expect_equal(
    res_auto$scores$mean_top1_similarity,
    res_embedding$scores$mean_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_auto$scores$mean_topk_similarity,
    res_embedding$scores$mean_topk_similarity,
    tolerance = 1e-12
  )
})

test_that("cv_alignment_rows runs with real fpgw fits in auto mode", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("lpSolve")

  set.seed(20260418)
  n <- 10
  latent_dim <- 2
  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  X1 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)
  X2 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  res <- cv_alignment_rows(
    hd,
    rows = list(list(a = 1:2), list(a = 3:4)),
    fit_fn = function(analysis) {
      fpgw(analysis, omega1 = 0.5, verbose = FALSE)
    },
    features = list(a = Z, b = Z),
    k = 2,
    target_pool = "analysis",
    prediction_mode = "auto"
  )

  expect_s3_class(res, "cv_result")
  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
  expect_true(all(is.finite(res$scores$mean_topk_similarity)))
})

test_that("cv_alignment_rows is permutation-invariant under consistent row remapping", {
  skip_if_not_installed("multidesign")

  predict.mock_dual_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_dual_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_dual_fit", predict.mock_dual_fit)

  fixture <- .make_dual_mode_fixture(seed = 20260409)
  rows <- list(
    list(a = 1:3),
    list(a = 4:6)
  )

  res_original <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  perm <- sample(seq_len(nrow(fixture$Za)))
  rows_perm <- lapply(rows, function(fold) {
    lapply(fold, function(idx) match(idx, perm))
  })

  md1_perm <- multidesign::multidesign(fixture$Za[perm, , drop = FALSE], data.frame(row_id = seq_along(perm)))
  md2_perm <- multidesign::multidesign(fixture$Zb[perm, , drop = FALSE], data.frame(row_id = seq_along(perm)))
  hd_perm <- multidesign::hyperdesign(list(a = md1_perm, b = md2_perm), block_names = c("a", "b"))
  features_perm <- list(
    a = fixture$features$a[perm, , drop = FALSE],
    b = fixture$features$b[perm, , drop = FALSE]
  )

  fit_fn_perm <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    Xa <- as.matrix(resolved$domains[[1]]$x)
    Xb <- as.matrix(resolved$domains[[2]]$x)
    d <- ncol(Xa)

    structure(
      list(
        s = rbind(Xa, Xb),
        v = rbind(diag(d), diag(d)),
        sdev = apply(rbind(Xa, Xb), 2, stats::sd),
        preproc = NULL,
        block_indices = list(a = seq_len(d), b = d + seq_len(d)),
        domain_names = resolved$domain_names,
        targets = list(a = Xa, b = Xb)
      ),
      class = c("mock_dual_fit", "multiblock_biprojector")
    )
  }

  res_perm <- cv_alignment_rows(
    hd_perm,
    rows = rows_perm,
    fit_fn = fit_fn_perm,
    features = features_perm,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  expect_equal(
    res_original$scores$mean_top1_similarity,
    res_perm$scores$mean_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_original$scores$mean_topk_similarity,
    res_perm$scores$mean_topk_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_original$scores$oracle_top1_similarity,
    res_perm$scores$oracle_top1_similarity,
    tolerance = 1e-12
  )
  expect_equal(
    res_original$scores$oracle_topk_similarity,
    res_perm$scores$oracle_topk_similarity,
    tolerance = 1e-12
  )
})

test_that("cv_alignment_rows keeps zero-norm feature rows finite", {
  skip_if_not_installed("multidesign")

  predict.mock_weight_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_weight_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_weight_fit", predict.mock_weight_fit)

  set.seed(20260407)
  n <- 10
  latent_dim <- 2
  Z <- matrix(rnorm(n * latent_dim), n, latent_dim)
  Z[c(1, 7), ] <- 0
  X1 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)
  X2 <- Z + matrix(rnorm(n * latent_dim, sd = 0.03), n, latent_dim)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = seq_len(n)))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = seq_len(n)))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  fit_fn <- function(analysis) {
    resolved <- manifoldalign:::resolve_hyperdesign(analysis)
    structure(
      list(
        domain_names = resolved$domain_names,
        targets = lapply(resolved$domains, function(dom) as.matrix(dom$x))
      ),
      class = "mock_weight_fit"
    )
  }

  res <- cv_alignment_rows(
    hd,
    rows = list(list(a = 1:3), list(a = 4:5)),
    fit_fn = fit_fn,
    features = list(a = Z, b = Z),
    k = 20,
    target_pool = "analysis",
    prediction_mode = "weights"
  )

  expect_true(all(is.finite(res$scores$mean_top1_similarity)))
  expect_true(all(is.finite(res$scores$mean_topk_similarity)))
  expect_true(all(is.finite(res$scores$oracle_top1_similarity)))
  expect_true(all(is.finite(res$scores$oracle_topk_similarity)))
})

test_that("cv_alignment_rows scores degrade when target features are shuffled", {
  skip_if_not_installed("multidesign")

  predict.mock_dual_fit <- function(object, newdata, from, to, type = c("weights", "transport"), ...) {
    type <- match.arg(type)
    if (!identical(type, "weights")) {
      stop("mock_dual_fit only supports type = 'weights'", call. = FALSE)
    }

    target_name <- if (is.character(to)) to else object$domain_names[[to]]
    target <- object$targets[[target_name]]
    dists <- as.matrix(stats::dist(rbind(newdata, target)))
    nq <- nrow(newdata)
    nt <- nrow(target)
    sim <- exp(-dists[seq_len(nq), nq + seq_len(nt), drop = FALSE])
    sim / rowSums(sim)
  }

  registerS3method("predict", "mock_dual_fit", predict.mock_dual_fit)

  fixture <- .make_dual_mode_fixture(seed = 20260414)
  rows <- list(
    list(a = 1:3),
    list(a = 4:6)
  )

  res_signal <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = fixture$features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  set.seed(20260416)
  shuffled_features <- list(
    a = fixture$features$a,
    b = fixture$features$b[sample(seq_len(nrow(fixture$features$b))), , drop = FALSE]
  )

  res_shuffled <- cv_alignment_rows(
    fixture$hd,
    rows = rows,
    fit_fn = fixture$fit_fn,
    features = shuffled_features,
    k = 3,
    target_pool = "analysis",
    prediction_mode = "embedding"
  )

  expect_gt(
    mean(res_signal$scores$mean_top1_similarity),
    mean(res_shuffled$scores$mean_top1_similarity) + 0.25
  )
  expect_gt(
    mean(res_signal$scores$mean_topk_similarity),
    mean(res_shuffled$scores$mean_topk_similarity) + 0.15
  )
})

test_that("cv_alignment_rows errors for one-sided assessment-only target pools", {
  skip_if_not_installed("multidesign")
  skip_if_not_installed("multivarious")
  skip_if_not_installed("adjoin")

  fixture <- .make_cv_alignment_fixture(seed = 20260412)

  expect_error(
    cv_alignment_rows(
      fixture$hd,
      rows = list(list(d1 = 1:3), list(d1 = 4:6)),
      fit_fn = fixture$fit_fn,
      features = fixture$features,
      k = 2,
      target_pool = "assessment"
    ),
    "No target rows were available for block `d2` under `target_pool = \"assessment\"`"
  )
})

test_that("cv_alignment_rows validates incompatible weights configuration", {
  skip_if_not_installed("multidesign")

  X1 <- matrix(rnorm(20), 5, 4)
  X2 <- matrix(rnorm(20), 5, 4)
  md1 <- multidesign::multidesign(X1, data.frame(row_id = 1:5))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = 1:5))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  expect_error(
    cv_alignment_rows(
      hd,
      rows = list(list(a = 1:2)),
      fit_fn = function(analysis) analysis,
      features = list(a = X1, b = X2),
      target_pool = "assessment",
      prediction_mode = "weights"
    ),
    "requires `target_pool = \"analysis\"`"
  )
})

test_that("cv_alignment_rows validates external feature row counts", {
  skip_if_not_installed("multidesign")

  X1 <- matrix(rnorm(20), 5, 4)
  X2 <- matrix(rnorm(20), 5, 4)

  md1 <- multidesign::multidesign(X1, data.frame(row_id = 1:5))
  md2 <- multidesign::multidesign(X2, data.frame(row_id = 1:5))
  hd <- multidesign::hyperdesign(list(a = md1, b = md2), block_names = c("a", "b"))

  bad_features <- list(
    a = matrix(rnorm(12), 4, 3),
    b = matrix(rnorm(15), 5, 3)
  )

  expect_error(
    cv_alignment_rows(
      hd,
      rows = list(list(a = 1:2, b = 1:2)),
      fit_fn = function(analysis) analysis,
      features = bad_features
    ),
    "Feature block `a` has 4 rows"
  )
})
