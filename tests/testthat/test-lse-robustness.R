test_that("linear_sim_embed has robust guards for edge cases", {
  set.seed(3)
  X <- matrix(rnorm(80 * 6), 80, 6)

  # (a) Zero-mask should abort optimisation
  M0  <- matrix(0, 80, 80)
  expect_error(
    linear_sim_embed(X, M = M0, ncomp = 2),
    "Mask matrix M sums to zero"
  )

  # (b) Pathologically tiny target similarity – triggers scaling warning
  Ttiny <- matrix(1e-12, 80, 80); diag(Ttiny) <- 1
  expect_warning(
    expect_warning(
      linear_sim_embed(X, T = Ttiny, ncomp = 2, sigma_P = 1),
      "poorly scaled"
    ),
    "too close to zero"
  )

  # (c) Predict with wrong dimensionality must stop
  fit <- linear_sim_embed(X, ncomp = 2, maxit = 300)
  bad  <- matrix(rnorm(10 * 5), 10, 5)              # 5 ≠ 6 columns
  expect_error(
    predict(fit, bad),
    "expects 6 columns"
  )
})

test_that("ADAM reevaluates the current objective independently of verbosity", {
  set.seed(812)
  X <- scale(matrix(rnorm(60), 20, 3))
  W0 <- diag(3)[, 1:2, drop = FALSE]
  T <- exp(-as.matrix(dist(X))^2 / 3)
  M <- matrix(1, 20, 20) - diag(20)
  run <- function(verbose) {
    expect_warning(
      fit <- manifoldalign:::.optimize_W_enhanced(
        X, T, M, W0, sigma = 3, alpha_p = .1,
        maxit = 2, tol = 1e-14, verbose = verbose
      ),
      "did not converge"
    )
    fit
  }
  quiet <- run(FALSE)
  invisible(capture.output(loud <- run(TRUE)))
  expect_identical(quiet$convergence, 1)
  expect_equal(quiet, loud)
  expect_gt(abs(diff(quiet$trace)), 1e-8)
  # Independent scalar objective at the initial iterate.
  P <- exp(-as.matrix(dist(X %*% W0))^2 / 3)
  expected <- .9 * sum(M * (P - T)^2) / (2 * sum(M)) +
    .1 * sum((crossprod(W0) - diag(2))^2) / 8
  expect_equal(quiet$trace[1], expected, tolerance = 1e-12)
})

test_that("the documented default ADAM fit converges with identical logging modes", {
  set.seed(42)
  X <- matrix(rnorm(100 * 10), 100, 10)
  quiet <- linear_sim_embed(X, ncomp = 3)
  invisible(capture.output(loud <- linear_sim_embed(X, ncomp = 3, verbose = TRUE)))
  expect_identical(quiet$convergence$convergence, 0L)
  expect_equal(quiet$weights, loud$weights)
  expect_equal(quiet$objective_trace, loud$objective_trace)
  expect_gt(length(quiet$objective_trace), 2L)
  expect_lt(tail(quiet$objective_trace, 1), quiet$objective_trace[1])
})
