test_that("linear_sim_embed produces near-orthogonal weights with high alpha_p", {
  set.seed(2)
  X <- matrix(rnorm(300 * 15), 300, 15)     # noisy, no obvious structure

  fit <- linear_sim_embed(X, ncomp = 5, alpha_p = 0.9,
                          sigma_P = "auto", alpha_schedule = TRUE,
                          maxit = 800, tol = 1e-7, verbose = FALSE)

  WtW <- crossprod(fit$weights)
  diag(WtW) <- diag(WtW) - 1
  orth_dev <- sqrt(sum(WtW^2))              # Frobenius ‖WᵀW - I‖

  expect_lt(orth_dev, 0.05)
}) 