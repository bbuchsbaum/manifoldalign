test_that("compute_grasp_assignment respects alpha blend", {
  set.seed(101)

  basis1 <- list(
    vectors = matrix(rnorm(24), nrow = 8, ncol = 3),
    values = c(0.9, 1.3, 1.8)
  )
  basis2 <- list(
    vectors = matrix(rnorm(24), nrow = 8, ncol = 3),
    values = c(1.0, 1.2, 1.9)
  )
  desc1 <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
  desc2 <- matrix(rnorm(8 * 5), nrow = 8, ncol = 5)
  M <- Matrix::Diagonal(3)

  out_cos <- manifoldalign:::compute_grasp_assignment(
    basis1, basis2, desc1, desc2, M,
    distance_method = "cosine",
    solver_method = "linear",
    alpha = 1
  )
  out_mix <- manifoldalign:::compute_grasp_assignment(
    basis1, basis2, desc1, desc2, M,
    distance_method = "cosine",
    solver_method = "linear",
    alpha = 0
  )

  expect_false(isTRUE(all.equal(out_cos$cost_matrix, out_mix$cost_matrix)))
})

test_that("compute_grasp_assignment euclidean cost is finite and non-negative", {
  set.seed(202)

  basis1 <- list(
    vectors = matrix(rnorm(30), nrow = 10, ncol = 3),
    values = c(0.8, 1.1, 1.6)
  )
  basis2 <- list(
    vectors = matrix(rnorm(30), nrow = 10, ncol = 3),
    values = c(0.7, 1.0, 1.5)
  )
  desc1 <- matrix(rnorm(10 * 4), nrow = 10, ncol = 4)
  desc2 <- matrix(rnorm(10 * 4), nrow = 10, ncol = 4)
  M <- Matrix::Diagonal(3)

  out <- manifoldalign:::compute_grasp_assignment(
    basis1, basis2, desc1, desc2, M,
    distance_method = "euclidean",
    solver_method = "linear"
  )

  expect_true(all(is.finite(out$cost_matrix)))
  expect_true(all(out$cost_matrix >= 0))
})

test_that("GRASP adapter uses the common non-trivial spectral subspace", {
  set.seed(1)
  # Domains of different sizes give spectral bases of different widths.
  X1 <- matrix(rnorm(50), 25, 2)
  X2 <- matrix(rnorm(24), 12, 2)

  bases <- manifoldalign:::compute_grasp_basis(
    list(list(x = X1), list(x = X2)),
    ncomp = 30
  )
  # The 12-node graph has fewer non-trivial Laplacian eigenvectors than the
  # 25-node one, and both bases are truncated to that common width.
  widths <- vapply(bases, function(b) ncol(b$vectors), integer(1))
  alone <- manifoldalign:::compute_grasp_basis(
    list(list(x = X1), list(x = X1)),
    ncomp = 30
  )
  expect_identical(widths[[1]], widths[[2]])
  expect_lt(widths[[1]], ncol(alone[[1]]$vectors))

  fit <- fit_pair(grasp_aligner(), X1, X2)
  expect_s3_class(fit, "grasp_pair_fit")
  expect_equal(fit$k, min(vapply(bases, function(x) ncol(x$vectors), integer(1))))
  expect_true(is.finite(fit$loss))
  expect_equal(dim(fit$rotation), c(fit$k, fit$k))
})

test_that("grasp_multiset stores anchor-aligned embeddings consistent with permutations", {
  set.seed(303)

  n <- 40
  base <- matrix(rnorm(n * 5), n, 5)
  p2 <- sample(n)
  p3 <- sample(n)
  graphs <- list(
    base,
    base[p2, ] + matrix(rnorm(n * 5, sd = 0.03), n, 5),
    base[p3, ] + matrix(rnorm(n * 5, sd = 0.03), n, 5)
  )

  fit <- manifoldalign::grasp_multiset(
    graphs,
    ncomp = 8,
    q_descriptors = 20,
    max_iter = 20
  )

  expect_true(!is.null(fit$embeddings_raw))
  expect_equal(length(fit$embeddings), length(fit$embeddings_raw))

  anchor_idx <- if (is.numeric(fit$anchor)) as.integer(fit$anchor) else 1L
  for (s in seq_along(fit$embeddings)) {
    if (s == anchor_idx) next
    perm <- as.integer(fit$permutations[[s]])
    if (length(perm) != nrow(fit$embeddings[[anchor_idx]])) next
    expect_equal(
      fit$embeddings[[s]],
      fit$embeddings_raw[[s]][perm, , drop = FALSE]
    )
  }
})

grasp_pair_hd <- function(X, Y) {
  multidesign::hyperdesign(list(
    a = multidesign::multidesign(X, data.frame(i = seq_len(nrow(X)))),
    b = multidesign::multidesign(Y, data.frame(i = seq_len(nrow(Y))))
  ))
}

test_that("GRASP graphs are invariant to rotations and row order of the features", {
  skip_if_not_installed("multidesign")
  # Regression: graphs were built on row-standardised features, which a
  # rotation changes and which collapse to two values with two features, so
  # GRASP matched a rotated copy at chance level.
  for (p in c(2L, 10L)) {
    set.seed(40 + p)
    n <- 45
    X <- matrix(rnorm(n * p), n) %*% matrix(rnorm(p * p), p)
    Q <- qr.Q(qr(matrix(rnorm(p * p), p)))
    perm <- sample(n)
    Y <- (X %*% Q)[perm, , drop = FALSE] + matrix(rnorm(n * p, sd = 1e-3), n)

    fit <- grasp(grasp_pair_hd(X, Y), ncomp = 6, q_descriptors = 8)
    truth <- match(seq_len(n), perm)
    expect_gte(mean(fit$assignment == truth), 0.95, label = paste("accuracy, p =", p))
  }
})

test_that("GRASP graph bandwidth is deterministic and leaves the caller's RNG alone", {
  skip_if_not_installed("multidesign")
  set.seed(7)
  X <- matrix(rnorm(30 * 4), 30, 4)
  Y <- X[sample(30), ] + matrix(rnorm(120, sd = 0.01), 30)
  hd <- grasp_pair_hd(X, Y)

  set.seed(99)
  before <- .Random.seed
  fit1 <- grasp(hd, ncomp = 5, q_descriptors = 8)
  expect_identical(.Random.seed, before)
  set.seed(12345)
  fit2 <- grasp(hd, ncomp = 5, q_descriptors = 8)
  expect_identical(fit1$assignment, fit2$assignment)
  expect_equal(fit1$s, fit2$s)

  # The bandwidth is the median k-th-neighbour distance, so it scales with
  # the data and does not depend on row order.
  bw <- manifoldalign:::.grasp_knn_bandwidth
  expect_equal(bw(3 * X, 5), 3 * bw(X, 5))
  expect_equal(bw(X[30:1, ], 5), bw(X, 5))
  kth <- apply(as.matrix(dist(X)), 1, function(d) sort(d)[6])
  expect_equal(bw(X, 5), stats::median(kth))
})

test_that("graph_metric = 'correlation' standardises rows and needs three features", {
  skip_if_not_installed("multidesign")
  set.seed(3)
  X2 <- matrix(rnorm(40), 20, 2)
  expect_error(
    grasp(grasp_pair_hd(X2, X2[20:1, ]), ncomp = 4, q_descriptors = 6,
          graph_metric = "correlation"),
    "at least three features"
  )

  rows <- manifoldalign:::.grasp_graph_rows(matrix(c(1, 2, 3, 2, 4, 9), 2, 3, byrow = TRUE),
                                            "correlation")
  expect_equal(rowMeans(rows), c(0, 0))
  expect_equal(apply(rows, 1, stats::sd), c(1, 1))

  # Row-wise affine changes leave the correlation graph unchanged.
  n <- 30
  X <- matrix(rnorm(n * 8), n, 8)
  perm <- sample(n)
  Y <- (2 * X + 5)[perm, ]
  fit <- grasp(grasp_pair_hd(X, Y), ncomp = 5, q_descriptors = 8,
               graph_metric = "correlation")
  expect_gte(mean(fit$assignment == match(seq_len(n), perm)), 0.95)
})
