# Regression tests for multi-set UOT template features, TI-Sinkhorn iteration
# counts, and the sparse cost builder.

# Dense entropic KL-UOT coupling from translation-invariant potentials
uot_dense_coupling <- function(cost, alpha, beta, fbar, gbar, epsilon) {
  outer(as.numeric(alpha), as.numeric(beta)) *
    exp((outer(as.numeric(fbar), as.numeric(gbar), "+") - cost) / epsilon)
}

test_that("learned template features are the weighted mean of transported features", {
  set.seed(42)
  K <- 3
  m <- 7
  D <- 2
  datasets <- lapply(seq_len(K), function(k) {
    n <- 5 + k
    list(X = matrix(rnorm(n * 2), n, 2),
         alpha = runif(n, 0.5, 2) * k,       # total masses differ from one
         F = matrix(rnorm(n * D, mean = k), n, D))
  })
  template <- list(Y = matrix(rnorm(m * 2), m, 2), beta = runif(m, 0.5, 1.5))
  omega <- c(0.5, 0.3, 0.2)                  # non-uniform dataset weights

  for (rescale in c("unit", "alpha")) {
    fit <- multiset_uot_align(datasets, template, omega = omega,
                              epsilon = 0.5, rho1 = 1, rho2 = 1,
                              lambda_coord = 1, neighbor_mode = "dense",
                              max_outer = 1, max_inner = 5000, tol_inner = 1e-12,
                              rescale_pi2 = rescale, learn_template_features = TRUE,
                              delta = 0)
    num <- matrix(0, m, D)
    den <- numeric(m)
    for (k in seq_len(K)) {
      fk <- fit$fits[[k]]
      Pk <- uot_dense_coupling(fk$cost, datasets[[k]]$alpha, template$beta,
                               fk$sol$fbar, fk$sol$gbar, 0.5)
      pi2 <- colSums(Pk)
      mass_k <- if (rescale == "alpha") sum(datasets[[k]]$alpha) else 1
      pi2 <- pi2 / sum(pi2) * mass_k
      mean_F <- crossprod(Pk, datasets[[k]]$F) / colSums(Pk)
      expect_equal(fk$pi2, pi2, tolerance = 1e-10)
      num <- num + omega[k] * pi2 * mean_F
      den <- den + omega[k] * pi2
    }
    expect_equal(fit$template$G, num / den, tolerance = 1e-10,
                 label = paste("template G, rescale_pi2 =", rescale))
    # beta is the pi2 mixture rescaled to the mean dataset mass
    m0 <- mean(vapply(datasets, function(d) sum(d$alpha), numeric(1)))
    expect_equal(fit$template$beta, den / sum(den) * m0, tolerance = 1e-10)
  }

  # With identical features in every dataset, G must reproduce them exactly
  # whatever the masses (the old code returned G * sum(pi2_mix) / m0).
  same <- lapply(datasets, function(d) {
    d$F <- matrix(c(2, -1), nrow(d$X), D, byrow = TRUE)
    d
  })
  fit <- multiset_uot_align(same, template, omega = omega, epsilon = 0.5,
                            rho1 = 1, rho2 = 1, neighbor_mode = "dense",
                            max_outer = 1, rescale_pi2 = "unit",
                            learn_template_features = TRUE, delta = 0)
  expect_equal(fit$template$G, matrix(c(2, -1), m, D, byrow = TRUE), tolerance = 1e-10)
})

test_that("TI-Sinkhorn reports the number of iterations actually run", {
  set.seed(1)
  C <- matrix(runif(30), 5, 6)
  a <- rep(1, 5)
  b <- rep(1, 6)
  sol <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.05, rho1 = 1, rho2 = 1,
                            max_iter = 3, tol = 1e-300)
  expect_false(sol$converged)
  expect_equal(sol$iterations, 3)

  sol <- uot_ti_sinkhorn_kl(C, a, b, epsilon = 0.5, rho1 = 1, rho2 = 1,
                            max_iter = 5000, tol = 1e-10)
  expect_true(sol$converged)
  expect_lte(sol$iterations, 5000)

  # Sparse kernels (CSR only and CSR + CSC)
  nn_idx <- matrix(c(1L, 2L, 3L, 4L, 5L, 2L, 3L, 4L, 5L, 6L), 5, 2)
  nn_d <- matrix(runif(10), 5, 2)
  sp <- uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = 6L)
  sol <- uot_ti_sinkhorn_kl(sp, a, b, epsilon = 0.05, rho1 = 1, rho2 = 1,
                            max_iter = 4, tol = 1e-300)
  expect_equal(sol$iterations, 4)
  sp_csr <- sp[c("row_ptr", "col_idx", "cost", "n_rows", "n_cols")]
  sol <- uot_ti_sinkhorn_kl(sp_csr, a, b, epsilon = 0.05, rho1 = 1, rho2 = 1,
                            max_iter = 4, tol = 1e-300)
  expect_equal(sol$iterations, 4)
})

test_that("sparse cost builder merges duplicate edges and defaults max_dist", {
  set.seed(7)
  n <- 4
  m <- 5
  # Row 1 lists template node 2 twice; row 3 gets template node 4 from both
  # the kNN matrix and the extra edges.
  nn_idx <- matrix(c(2L, 1L, 4L, 5L,
                     2L, 3L, 1L, 4L), n, 2)
  nn_d <- matrix(c(0.5, 0.2, 0.3, 0.4,
                   0.5, 0.6, 0.7, 0.8), n, 2)
  dup <- uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = m,
                                       extra_i = c(3L, 2L), extra_j = c(4L, 5L),
                                       extra_dists = c(0.3, 0.9))
  # unique edges: (1,2) (2,1) (2,3) (2,5) (3,4) (3,1) (4,5) (4,4)
  expect_equal(dup$nnz, 8L)
  expect_equal(diff(dup$row_ptr), c(1L, 3L, 2L, 2L))
  expect_equal(sum(diff(dup$col_ptr)), 8L)
  # within-row columns are sorted and unique
  for (i in seq_len(n)) {
    cols <- dup$col_idx[(dup$row_ptr[i] + 1):dup$row_ptr[i + 1]]
    expect_false(anyDuplicated(cols) > 0)
    expect_equal(cols, sort(cols))
  }

  # Same graph built without duplicates gives identical solver output; with
  # double-counted edges the softmins would differ.
  clean <- uot_build_sparse_cost_knn_cpp(
    matrix(c(2L, 1L, 4L, 5L, NA, 3L, 1L, 4L), n, 2),
    matrix(c(0.5, 0.2, 0.3, 0.4, NA, 0.6, 0.7, 0.8), n, 2),
    n_cols = m, extra_i = 2L, extra_j = 5L, extra_dists = 0.9)
  expect_equal(clean$nnz, 8L)
  alpha <- rep(1, n)
  beta <- rep(1, m)
  s_dup <- uot_ti_sinkhorn_kl(dup, alpha, beta, epsilon = 0.3, rho1 = 1, rho2 = 1,
                              max_iter = 2000, tol = 1e-12)
  s_clean <- uot_ti_sinkhorn_kl(clean, alpha, beta, epsilon = 0.3, rho1 = 1, rho2 = 1,
                                max_iter = 2000, tol = 1e-12)
  expect_equal(s_dup$fbar, s_clean$fbar, tolerance = 1e-10)
  expect_equal(s_dup$gbar, s_clean$gbar, tolerance = 1e-10)

  # max_dist now has a default (NA = no threshold), equivalent to Inf
  a <- uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = m)
  b <- uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = m, max_dist = Inf)
  expect_identical(a, b)
  c <- uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = m, max_dist = 0.45)
  expect_equal(c$nnz, 3L)  # distances 0.2, 0.3, 0.4 survive
  expect_error(uot_build_sparse_cost_knn_cpp(nn_idx, nn_d, n_cols = m, max_dist = -1),
               "max_dist")
})
