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
    linear_sim_embed(X, T = Ttiny, ncomp = 2, sigma_P = 1),
    "poorly scaled"
  )

  # (c) Predict with wrong dimensionality must stop
  fit <- linear_sim_embed(X, ncomp = 2, maxit = 300)
  bad  <- matrix(rnorm(10 * 5), 10, 5)              # 5 ≠ 6 columns
  expect_error(
    predict(fit, bad),
    "expects 6 columns"
  )
}) 