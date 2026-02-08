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
