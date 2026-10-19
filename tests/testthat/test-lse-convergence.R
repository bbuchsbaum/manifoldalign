# Stopping semantics of linear_sim_embed(): the stopping decision must not
# depend on logging, and an exhausted budget must be reported as such.

lse_fixture <- function(seed = 1, n = 60, d = 5) {
  set.seed(seed)
  matrix(rnorm(n * d), n, d)
}

test_that("verbose logging does not change the R optimizer's answer", {
  X <- lse_fixture()

  quiet <- linear_sim_embed(X, ncomp = 2, sigma_P = 1, verbose = FALSE)
  out <- capture.output(
    loud <- linear_sim_embed(X, ncomp = 2, sigma_P = 1, verbose = TRUE)
  )
  expect_true(length(out) > 0)

  expect_identical(loud$weights, quiet$weights)
  expect_identical(loud$objective_trace, quiet$objective_trace)
  expect_identical(loud$convergence, quiet$convergence)
})

test_that("default R run does not stop after a single update", {
  X <- lse_fixture()
  fit <- linear_sim_embed(X, ncomp = 2, sigma_P = 1)

  expect_equal(fit$convergence$convergence, 0)
  expect_identical(fit$convergence$message, "SUCCESS")
  # Regression: the objective used to be evaluated only when logging, so the
  # default run reported convergence after one update.
  expect_gt(fit$convergence$iterations, 5L)
  # The trace holds the objective at the start and after every update, and the
  # stopping rule (5 consecutive changes below tol) holds at the end.
  tr <- fit$objective_trace
  expect_length(tr, fit$convergence$iterations + 1L)
  expect_true(all(abs(diff(tail(tr, 6))) < 1e-6))
  expect_lt(tail(tr, 1), tr[1])
  expect_equal(fit$convergence$final_value, tail(tr, 1))

  # The reported trace is the documented objective at the returned weights.
  ns <- asNamespace("manifoldalign")
  Xs <- scale(X)
  obj <- ns$.lse_objective(Xs, fit$weights, fit$target_sim, fit$mask,
                           fit$sigma_P, fit$alpha_p)$objective
  expect_equal(tail(tr, 1), obj, tolerance = 1e-12)
})

test_that("an exhausted iteration budget is reported as non-convergence", {
  X <- lse_fixture()
  # The public function refuses to return unconverged weights, with or
  # without logging.
  for (verbose in c(FALSE, TRUE)) {
    capture.output(
      expect_error(
        expect_warning(
          linear_sim_embed(X, ncomp = 2, sigma_P = 1, maxit = 3, tol = 1e-12,
                           verbose = verbose),
          "did not converge within 3 iterations"
        ),
        "failed to converge: MAXITER_REACHED"
      )
    )
  }

  # The optimizer itself reports the exhausted budget and a full trace.
  ns <- asNamespace("manifoldalign")
  Xs <- scale(X)
  T <- ns$gaussian_sim(Xs, 1)
  M <- matrix(1, nrow(Xs), nrow(Xs))
  diag(M) <- 0
  W0 <- svd(Xs, nu = 0, nv = 2)$v
  expect_warning(
    res <- ns$.optimize_W_enhanced(Xs, T, M, W0, sigma = 1, maxit = 3,
                                   tol = 1e-12),
    "did not converge within 3 iterations"
  )
  expect_equal(res$convergence, 1)
  expect_identical(res$message, "MAXITER_REACHED")
  expect_identical(res$iterations, 3L)
  expect_length(res$trace, 4L)
})

test_that("weights are returned as a base matrix", {
  X <- lse_fixture()
  fit <- linear_sim_embed(X, ncomp = 2, sigma_P = 1)
  expect_true(is.matrix(fit$weights))
  expect_false(isS4(fit$weights))
})

test_that("C++ backend maps tol to pgtol and reports iterations without RNG warnings", {
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))
  set.seed(456)
  n <- 15
  X <- scale(matrix(rnorm(n * 3), n, 3))
  labels <- rep(c("A", "B", "C"), each = n / 3)
  T <- outer(labels, labels, "==") * 1
  M <- matrix(1, n, n)
  diag(M) <- 0

  # Regression: arma_rng::set_seed_random() warned under RcppArmadillo.
  expect_no_warning(
    tight <- linear_sim_embed(X, T = T, M = M, sigma_P = 1, ncomp = 2,
                              alpha_p = 0.5, tol = 1e-10, use_cpp = TRUE)
  )
  loose <- linear_sim_embed(X, T = T, M = M, sigma_P = 1, ncomp = 2,
                            alpha_p = 0.5, tol = 1e-1, use_cpp = TRUE)

  expect_equal(tight$convergence$convergence, 0)
  expect_equal(loose$convergence$convergence, 0)
  expect_true(is.numeric(tight$convergence$iterations))
  expect_false(is.na(tight$convergence$iterations))
  # tol used to be passed as reltol, which L-BFGS-B ignores.
  expect_lt(loose$convergence$iterations, tight$convergence$iterations)
  expect_match(loose$convergence$message, "PGTOL")
})
