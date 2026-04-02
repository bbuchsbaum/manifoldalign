library(testthat)
library(multidesign)

quick_hd_orig <- function(Xlist, ylist) {
  md_list <- Map(function(x, y) {
    multidesign::multidesign(x, data.frame(lbl = factor(y)))
  }, Xlist, ylist)
  names(md_list) <- paste0("domain", seq_along(md_list))
  multidesign::hyperdesign(md_list)
}

test_that("kema_orig full solver returns expected structure", {
  skip_if_not_installed("PRIMME")

  n <- 30
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.5, n, p)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit <- kema_orig(
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
  expect_true(is.null(fit$retry_info))
  expect_true(is.null(fit$regression_quality))
})

test_that("kema_orig REKEMA path uses reduced rank coefficients", {
  skip_if_not_installed("PRIMME")

  n <- 50
  p <- 5
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) - 0.25, n, p)
  y <- rep(1:4, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit <- kema_orig(
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
  skip_if_not_installed("PRIMME")

  set.seed(101)
  n <- 24
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) + 0.2, n, p)
  y <- rep(1:3, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit_full <- kema_orig(
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

  fit_reduced <- kema_orig(
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

  fit_operator <- kema_orig(
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

test_that("kema_orig auto backend policy and backend constraints are enforced", {
  skip_if_not_installed("PRIMME")

  set.seed(202)
  n <- 28
  p <- 4
  X1 <- matrix(rnorm(n * p), n, p)
  X2 <- matrix(rnorm(n * p) - 0.3, n, p)
  y <- rep(1:4, length.out = n)
  hd <- quick_hd_orig(list(X1, X2), list(y, y))

  fit_auto_reduced <- kema_orig(
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
      fidelity_action = "warn"
    )
  )
  expect_identical(fit_auto_reduced$backend, "reduced_exact")

  fit_auto_operator <- kema_orig(
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
      fidelity_action = "warn"
    )
  )
  expect_identical(fit_auto_operator$backend, "operator_exact")

  expect_error(
    kema_orig(
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
