# SSMA numerical contract: generalized eigenpairs of a diagonal pencil, the
# documented degree-matrix estimand (joint by default), explicit component deficits, and sampling.

ssma_pencil_oracle <- function(a, b, ncomp, regularization) {
  r <- length(a)
  ridge <- max(regularization, 1e-10 * mean(b + 1e-12))
  sort(a / (b + ridge))[seq_len(ncomp)]
}

test_that("reduced SSMA solver matches a diagonal-pencil oracle across the dense/iterative boundary", {
  ncomp <- 2L
  # k_req = ncomp + 5 = 7: r <= 8 uses the dense path, r >= 9 the iterative one.
  for (r in c(6L, 7L, 8L, 9L, 10L, 25L)) {
    set.seed(r)
    a <- sample(seq_len(r))            # spectrum 1..r in shuffled order
    b <- stats::runif(r, 0.5, 2)
    sol <- manifoldalign:::.ssma_solve_reduced_geigen(
      diag(as.numeric(a)), diag(b), ncomp = ncomp, eig_tol = 1e-8,
      regularization = 1e-12
    )
    expected <- ssma_pencil_oracle(a, b, ncomp, 1e-12)
    expect_equal(sol$eigenvalues, expected, tolerance = 1e-8,
                 label = paste0("eigenvalues r=", r))
    # Generalized eigen-equation A alpha = lambda B alpha.
    A <- diag(as.numeric(a))
    B <- diag(b)
    resid <- A %*% sol$alpha - B %*% sol$alpha %*% diag(sol$eigenvalues)
    expect_lt(max(abs(resid)), 1e-6, label = paste0("residual r=", r))
  }
})

test_that("dense SSMA fallback keeps the smallest eigenvalue of diag(1..8)", {
  sol <- manifoldalign:::.ssma_solve_reduced_geigen(
    diag(as.numeric(1:8)), diag(8), ncomp = 2L, eig_tol = 1e-8,
    regularization = 1e-12
  )
  expect_equal(sol$eigenvalues, c(1, 2), tolerance = 1e-8)
  expect_equal(min(sol$eigenvalues_all), 1, tolerance = 1e-8)
})

test_that("SSMA never substitutes trivial eigenvectors and reports deficits", {
  expect_error(
    manifoldalign:::.ssma_solve_reduced_geigen(
      matrix(0, 6, 6), diag(6), ncomp = 2L, eig_tol = 1e-8,
      regularization = 1e-12
    ),
    "no non-trivial"
  )
  expect_warning(
    sol <- manifoldalign:::.ssma_solve_reduced_geigen(
      diag(c(0, 0, 0, 0, 1, 2)), diag(6), ncomp = 4L, eig_tol = 1e-8,
      regularization = 1e-12
    ),
    "returned 2 of the 4 requested"
  )
  expect_equal(sol$eigenvalues, c(1, 2), tolerance = 1e-8)
  expect_equal(ncol(sol$alpha), 2L)

  # Iterative path with many trivial pairs falls back to the full spectrum.
  a <- c(rep(0, 12), 3, 1, 2, 5, 4, 6, 7, 8)
  sol2 <- manifoldalign:::.ssma_solve_reduced_geigen(
    diag(a), diag(length(a)), ncomp = 3L, eig_tol = 1e-8,
    regularization = 1e-12
  )
  expect_equal(sol2$eigenvalues, c(1, 2, 3), tolerance = 1e-8)
})

test_that("ssma_align records requested and returned component counts", {
  set.seed(7)
  dat <- generate_identity_alignment(n = 30, d = 5, noise_sd = 0.03, seed = 7)
  corr <- data.frame(domain_i = 1L, index_i = 1:30, domain_j = 2L, index_j = 1:30)
  fit <- ssma_align(
    list(d1 = dat$X1, d2 = dat$X2), correspondences = corr, ncomp = 3,
    control = ssma_align_control(knn = 8, sigma = 1, rank_per_domain = 5)
  )
  expect_identical(fit$ncomp_requested, 3L)
  expect_identical(fit$ncomp_returned, 3L)
  expect_identical(fit$degree_matrix, "joint")

  # Rank 2 per domain gives at most 4 reduced directions; asking for more
  # than exist is a recorded, warned deficit rather than a silent truncation.
  expect_warning(
    fit_short <- ssma_align(
      list(d1 = dat$X1, d2 = dat$X2), correspondences = corr, ncomp = 6,
      control = ssma_align_control(knn = 8, sigma = 1, rank_per_domain = 2)
    ),
    "requested components"
  )
  expect_identical(fit_short$ncomp_requested, 6L)
  expect_lt(fit_short$ncomp_returned, 6L)
  expect_equal(ncol(fit_short$s), fit_short$ncomp_returned)
})

ssma_estimand_fixture <- function(degree_matrix) {
  set.seed(11)
  dat <- generate_identity_alignment(n = 25, d = 4, noise_sd = 0.05, seed = 11)
  corr <- data.frame(domain_i = 1L, index_i = 1:25, domain_j = 2L, index_j = 1:25,
                     weight = 5)
  ctrl <- ssma_align_control(knn = 6, sigma = 1, rank_per_domain = 4,
                             regularization = 1e-12, corr_edge_weight = 5,
                             degree_matrix = degree_matrix)
  fit <- ssma_align(list(d1 = dat$X1, d2 = dat$X2), correspondences = corr,
                    preproc = multivarious::center(), ncomp = 3, control = ctrl)

  X <- list(scale(dat$X1, scale = FALSE), scale(dat$X2, scale = FALSE))
  within <- lapply(X, function(x) {
    manifoldalign:::.ssma_build_within_graph(x, knn = 6, sigma = 1)
  })
  L_joint <- manifoldalign:::.ssma_build_joint_laplacian(
    within, fit$correspondences, domain_sizes = c(25L, 25L)
  )
  # Independent joint-degree oracle: row sums of the dense joint adjacency
  # (within-domain kNN blocks plus symmetric correspondence edges).
  W <- as.matrix(Matrix::bdiag(lapply(within, `[[`, "A")))
  cr <- fit$correspondences
  gi <- (cr$domain_i - 1L) * 25L + cr$index_i
  gj <- (cr$domain_j - 1L) * 25L + cr$index_j
  for (r in seq_along(gi)) {
    W[gi[r], gj[r]] <- W[gi[r], gj[r]] + cr$weight[r]
    W[gj[r], gi[r]] <- W[gj[r], gi[r]] + cr$weight[r]
  }
  list(
    fit = fit, L_joint = L_joint,
    D_within = unlist(lapply(within, `[[`, "degree")),
    D_joint = rowSums(W)
  )
}

test_that("SSMA solves the documented estimand: joint Laplacian and joint degrees", {
  fx <- ssma_estimand_fixture("joint")
  expect_identical(fx$fit$degree_matrix, "joint")
  S <- as.matrix(fx$fit$s)

  # The joint degree used by the fit is the row sum of the joint adjacency.
  deg <- manifoldalign:::.ssma_degree_list(
    list(list(degree = fx$D_within[1:25]), list(degree = fx$D_within[26:50])),
    fx$fit$correspondences, c(25L, 25L), "joint"
  )
  expect_equal(unlist(deg), fx$D_joint, tolerance = 1e-12)

  # Generalized eigenvectors are D-orthonormal for the joint degree matrix...
  G_joint <- crossprod(S, S * fx$D_joint)
  expect_lt(max(abs(G_joint - diag(ncol(S)))), 1e-6)
  # ...and not for the within-domain one.
  G_within <- crossprod(S, S * fx$D_within)
  expect_gt(max(abs(G_within - diag(ncol(S)))), 1e-2)
  # Rayleigh quotients with the joint Laplacian equal the eigenvalues.
  rq <- diag(crossprod(S, as.matrix(fx$L_joint %*% S)))
  expect_equal(rq, fx$fit$generalized_eigenvalues, tolerance = 1e-6)
})

test_that("degree_matrix = 'within_domain' reproduces the within-domain constraint", {
  fx <- ssma_estimand_fixture("within_domain")
  expect_identical(fx$fit$degree_matrix, "within_domain")
  S <- as.matrix(fx$fit$s)

  G_within <- crossprod(S, S * fx$D_within)
  expect_lt(max(abs(G_within - diag(ncol(S)))), 1e-6)
  G_joint <- crossprod(S, S * fx$D_joint)
  expect_gt(max(abs(G_joint - diag(ncol(S)))), 1e-2)
  rq <- diag(crossprod(S, as.matrix(fx$L_joint %*% S)))
  expect_equal(rq, fx$fit$generalized_eigenvalues, tolerance = 1e-6)

  expect_error(ssma_align_control(degree_matrix = "joint_graph"), "should be one of")
})

test_that("reduced and operator SSMA solvers agree", {
  set.seed(5)
  dat <- generate_identity_alignment(n = 30, d = 5, noise_sd = 0.05, seed = 5)
  corr <- data.frame(domain_i = 1L, index_i = 1:30, domain_j = 2L, index_j = 1:30)
  fits <- lapply(c("reduced", "operator"), function(solver) {
    ssma_align(list(d1 = dat$X1, d2 = dat$X2), correspondences = corr, ncomp = 3,
               control = ssma_align_control(knn = 8, sigma = 1, rank_per_domain = 5,
                                            regularization = 1e-10, solver = solver))
  })
  expect_equal(fits[[1]]$generalized_eigenvalues, fits[[2]]$generalized_eigenvalues,
               tolerance = 1e-5)
})

test_that("label correspondences index the labelled rows, including singletons", {
  labels <- list(c(NA, NA, NA, NA, "a", "b", "b"), c("a", NA, "b", NA, "b", "b"))
  set.seed(1)
  before <- .Random.seed
  corr <- manifoldalign:::.ssma_build_label_correspondences(labels, 10L, 1, seed = 3)
  expect_identical(.Random.seed, before)
  a_row <- corr[corr$index_i == 5L | corr$index_j == 1L, ]
  expect_equal(nrow(a_row), 1L)
  expect_identical(a_row$index_i, 5L)
  expect_identical(a_row$index_j, 1L)
  expect_true(all(labels[[1]][corr$index_i] == labels[[2]][corr$index_j]))
  # Repeat with many seeds: a singleton label never pairs an unrelated row.
  for (sd in 1:25) {
    cc <- manifoldalign:::.ssma_build_label_correspondences(
      list(c(rep(NA, 9), "z"), c("z", rep(NA, 4))), 5L, 1, seed = sd
    )
    expect_identical(cc$index_i, 10L)
    expect_identical(cc$index_j, 1L)
  }
})

test_that("ssma_align does not disturb the caller's RNG stream", {
  dat <- generate_identity_alignment(n = 40, d = 5, noise_sd = 0.05, seed = 3)
  set.seed(99)
  before <- .Random.seed
  fit <- ssma_align(list(d1 = dat$X1, d2 = dat$X2),
                    y = NULL, ncomp = 2,
                    control = ssma_align_control(knn = 6, rank_per_domain = 5,
                                                 solver = "operator"))
  expect_identical(.Random.seed, before)
})
