library(testthat)
library(multidesign)

quick_hd_orig <- function(Xlist, ylist) {
  md_list <- Map(function(x, y) {
    multidesign::multidesign(x, data.frame(lbl = factor(y)))
  }, Xlist, ylist)
  names(md_list) <- paste0("domain", seq_along(md_list))
  multidesign::hyperdesign(md_list)
}

test_that("kema_orig public entry point is deprecated in favor of kema", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(7)
  X1 <- matrix(rnorm(36), 12, 3)
  X2 <- matrix(rnorm(36) + 0.2, 12, 3)
  y <- rep(1:3, length.out = 12)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  expect_warning(
    fit <- kema_orig(
      hd,
      y = lbl,
      ncomp = 1,
      knn = 3,
      lambda = 1e-3
    ),
    "kema_orig\\(\\) is deprecated.*Use kema.*1\\.0\\.0"
  )
  expect_true(fit$fidelity$passed)
})

test_that("KEMA excludes numerical zero modes using a spectral-scale threshold", {
  values <- c(-1.3e-10, 0.2, 0.3)
  selected <- manifoldalign:::select_nontrivial_eigenpairs(
    values,
    diag(3),
    ncomp = 2
  )

  expect_equal(selected$values, c(0.2, 0.3))
  expect_gt(selected$zero_tol, abs(values[[1L]]))
  expect_equal(selected$vectors, diag(3)[, 2:3, drop = FALSE])
})

test_that("KEMA rejects materially negative modes", {
  expect_error(
    manifoldalign:::select_nontrivial_eigenpairs(
      c(-0.1, 0.2, 0.3),
      diag(3),
      ncomp = 2
    ),
    "materially negative.*positive-semidefinite"
  )
})

test_that("KEMA quotient deflates graph nullity beyond the old eigenpair budget", {
  set.seed(20260827)
  pair_adjacency <- matrix(c(0, 1, 1, 0), 2, 2)
  W <- Matrix::bdiag(replicate(10, pair_adjacency, simplify = FALSE))
  M <- Matrix::Diagonal(x = Matrix::rowSums(W)) - W
  component_basis <- manifoldalign:::kema_graph_component_basis(W)
  quotient <- manifoldalign:::build_kema_kernel_quotient(
    Ks = list(diag(20)),
    full_form = TRUE,
    component_basis = component_basis
  )
  control <- manifoldalign:::default_kema_backend_control()
  control$dense_exact_max_dim <- 3L

  fit <- manifoldalign:::solve_kema_quotient(
    quotient = quotient,
    M = M,
    Ld = Matrix::Matrix(0, 20, 20, sparse = TRUE),
    lambda = 1,
    ncomp = 3,
    control = control
  )

  expect_equal(ncol(component_basis), 10L)
  expect_equal(
    as.matrix(Matrix::crossprod(component_basis)),
    diag(10),
    tolerance = 1e-12
  )
  expect_equal(as.matrix(M %*% component_basis), matrix(0, 20, 10))
  expect_equal(fit$values, rep(2, 3), tolerance = 1e-8)
  expect_identical(fit$eigensolver, "eigencore_lobpcg_deflated")
  expect_equal(fit$eigensolver_stats$graph_nullity, 10L)
  expect_equal(fit$eigensolver_stats$nullity_deflated, 10L)
  expect_true(fit$eigensolver_stats$certificate$passed)

  permutation <- sample.int(20)
  W_permuted <- W[permutation, permutation]
  M_permuted <- M[permutation, permutation]
  quotient_permuted <- manifoldalign:::build_kema_kernel_quotient(
    Ks = list(diag(20)),
    full_form = TRUE,
    component_basis = manifoldalign:::kema_graph_component_basis(W_permuted)
  )
  fit_permuted <- manifoldalign:::solve_kema_quotient(
    quotient = quotient_permuted,
    M = M_permuted,
    Ld = Matrix::Matrix(0, 20, 20, sparse = TRUE),
    lambda = 1,
    ncomp = 3,
    control = control
  )
  expect_equal(fit_permuted$values, fit$values, tolerance = 1e-8)
})

test_that("KEMA kernel quotient removes numerical kernel null directions", {
  x <- c(1, 2, 3)
  K <- tcrossprod(x)
  W <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3)
  quotient <- manifoldalign:::build_kema_kernel_quotient(
    Ks = list(K),
    full_form = TRUE,
    component_basis = manifoldalign:::kema_graph_component_basis(W)
  )

  expect_equal(quotient$rank, 1L)
  expect_equal(quotient$discarded_rank, 2L)
  expect_equal(
    as.matrix(quotient$H %*% Matrix::t(quotient$V)),
    K,
    tolerance = 1e-12
  )
})

test_that("KEMA requires positive metric regularization", {
  X1 <- matrix(seq_len(18), 6, 3)
  X2 <- X1 + 0.1
  y <- rep(1:2, each = 3)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  expect_error(
    manifoldalign:::kema_orig.hyperdesign(hd, y = lbl, lambda = 0),
    "lambda must be strictly positive.*SPD metric"
  )
})

test_that("KEMA fidelity matches a tiny generalized-eigen oracle", {
  A <- diag(c(1, 4, 9))
  B <- diag(c(1, 2, 3))
  vectors <- diag(3)
  values <- c(1, 2, 3)

  fidelity <- manifoldalign:::compute_kema_fidelity(
    values,
    vectors,
    A_apply = function(x) A %*% x,
    B_apply = function(x) B %*% x,
    residual_tol = 1e-12,
    orth_tol = 1e-12
  )
  expect_true(fidelity$passed)
  expect_lte(fidelity$max_rel_residual, 1e-12)
  expect_lte(fidelity$max_b_orth_offdiag, 1e-12)

  wrong <- manifoldalign:::compute_kema_fidelity(
    values + c(0, 0.1, 0),
    vectors,
    A_apply = function(x) A %*% x,
    B_apply = function(x) B %*% x,
    residual_tol = 1e-6,
    orth_tol = 1e-12
  )
  expect_false(wrong$passed)
  expect_gt(wrong$max_rel_residual, 1e-3)
})

test_that("KEMA defaults to a hard fidelity gate", {
  control <- manifoldalign:::default_kema_backend_control()
  expect_equal(control$fidelity_residual_tol, 1e-6)
  expect_equal(control$fidelity_orth_tol, 1e-6)
  expect_equal(control$fidelity_action, "fallback")
  expect_error(
    manifoldalign:::normalize_kema_backend_control(list(fidelity_action = "warn")),
    "must be one of.*fallback.*error"
  )
})

test_that("kema_orig full solver returns expected structure", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  n <- 30
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.5, n, p)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::vanilladot(),
    sample_frac = 1,
    mu = 0.5,
    lambda = 1e-4
  )

  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(ncol(fit$s), 2)
  expect_equal(ncol(fit$alpha), 2)
  expect_true(is.character(fit$formulation))
  expect_match(fit$formulation, "kema_orig")
  expect_true(fit$fidelity$passed)
  expect_lte(fit$fidelity$max_rel_residual, 1e-6)
  expect_lte(fit$fidelity$max_b_orth_offdiag, 1e-6)
  expect_true(all(fit$eigenvalues > fit$fidelity$eigenvalue_zero_tol))
  expect_true(is.null(fit$retry_info))
  expect_true(is.null(fit$regression_quality))
})

test_that("kema_orig never returns a fit after fidelity failure", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(19)
  n <- 18
  X1 <- matrix(rnorm(n * 3), n, 3)
  X2 <- matrix(rnorm(n * 3) + 0.2, n, 3)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  testthat::local_mocked_bindings(
    compute_kema_fidelity = function(...) {
      list(
        passed = FALSE,
        max_rel_residual = 0.25,
        max_b_orth_offdiag = 0.1,
        min_b_gram_diag = 1,
        max_b_gram_diag = 1,
        residuals_rel = c(0.25, 0.2)
      )
    },
    .package = "manifoldalign"
  )

  expect_error(
    manifoldalign:::kema_orig.hyperdesign(
      hd,
      y = lbl,
      ncomp = 2,
      knn = 4,
      kernel = kernlab::vanilladot(),
      sample_frac = 1,
      lambda = 1e-4,
      backend = "full_exact",
      backend_control = list(fidelity_action = "error")
    ),
    "No KEMA fit was returned"
  )
})

test_that("kema_orig REKEMA path uses reduced rank coefficients", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  n <- 50
  p <- 5
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) - 0.25, n, p)
  y <- rep(1:4, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 5,
    kernel = kernlab::rbfdot(sigma = 0.5),
    sample_frac = 0.4,
    mu = 0.5,
    lambda = 1e-4
  )

  n_total <- 2 * n
  expect_true(inherits(fit, "multiblock_biprojector"))
  expect_equal(ncol(fit$s), 2)
  expect_equal(ncol(fit$alpha), 2)
  expect_lt(nrow(fit$alpha), n_total)
})

test_that("kema_orig supports explicit backend selection with fidelity diagnostics", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(101)
  n <- 24
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.2, n, p)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit_full <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::vanilladot(),
    sample_frac = 1,
    lambda = 1e-4,
    backend = "full_exact"
  )
  expect_identical(fit_full$backend, "full_exact")
  expect_true(is.list(fit_full$fidelity))
  expect_true(is.logical(fit_full$fidelity$passed))
  expect_true(all(is.finite(fit_full$fidelity$residuals_rel)))

  fit_reduced <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::rbfdot(sigma = 0.5),
    sample_frac = 0.5,
    lambda = 1e-4,
    backend = "reduced_exact"
  )
  expect_identical(fit_reduced$backend, "reduced_exact")

  fit_operator <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::rbfdot(sigma = 0.5),
    sample_frac = 0.5,
    lambda = 1e-4,
    backend = "operator_exact"
  )
  expect_identical(fit_operator$backend, "operator_exact")
  expect_true(is.numeric(fit_operator$fidelity$max_rel_residual))
})

test_that("KEMA matrix and operator backends agree up to subspace indeterminacy", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(91)
  n <- 18
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.25, n, p)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fits <- lapply(
    c("full_exact", "reduced_exact", "operator_exact"),
    function(backend) {
      manifoldalign:::kema_orig.hyperdesign(
        hd,
        y = lbl,
        ncomp = 2,
        knn = 4,
        kernel = kernlab::rbfdot(sigma = 0.5),
        sample_frac = 1,
        lambda = 1e-3,
        backend = backend
      )
    }
  )
  names(fits) <- c("full", "reduced", "operator")

  for (fit in fits) {
    expect_true(fit$fidelity$passed)
    expect_lte(fit$fidelity$max_rel_residual, 1e-6)
    expect_lte(fit$fidelity$max_b_orth_offdiag, 1e-6)
  }
  expect_equal(
    fits$reduced$eigenvalues,
    fits$full$eigenvalues,
    tolerance = 1e-8
  )
  expect_equal(
    fits$operator$eigenvalues,
    fits$full$eigenvalues,
    tolerance = 1e-7
  )

  Q_full <- qr.Q(qr(as.matrix(fits$full$s)))
  for (name in c("reduced", "operator")) {
    Q_other <- qr.Q(qr(as.matrix(fits[[name]]$s)))
    canonical_correlations <- svd(
      crossprod(Q_full, Q_other),
      nu = 0,
      nv = 0
    )$d
    expect_gt(min(canonical_correlations), 1 - 1e-7)
  }
})

test_that("kema_orig auto backend policy and backend constraints are enforced", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  set.seed(202)
  n <- 28
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) - 0.3, n, p)
  y <- rep(1:4, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit_auto_reduced <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::vanilladot(),
    sample_frac = 1,
    lambda = 1e-4,
    backend = "auto",
    backend_control = list(
      full_exact_max_n = 10L,
      reduced_exact_max_r = 200L,
      fidelity_action = "fallback"
    )
  )
  expect_identical(fit_auto_reduced$backend, "reduced_exact")

  fit_auto_operator <- manifoldalign:::kema_orig.hyperdesign(
    hd,
    y = lbl,
    ncomp = 2,
    knn = 4,
    kernel = kernlab::rbfdot(sigma = 0.6),
    sample_frac = 0.4,
    lambda = 1e-4,
    backend = "auto",
    backend_control = list(
      reduced_exact_max_r = 8L,
      fidelity_action = "fallback"
    )
  )
  expect_identical(fit_auto_operator$backend, "operator_exact")

  expect_error(
    manifoldalign:::kema_orig.hyperdesign(
      hd,
      y = lbl,
      ncomp = 2,
      knn = 4,
      sample_frac = 0.5,
      lambda = 1e-4,
      backend = "full_exact"
    ),
    "requires sample_frac = 1"
  )
})
