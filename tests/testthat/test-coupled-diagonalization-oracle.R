library(testthat)

cd_oracle_orthonormal <- function(n, k) {
  qr.Q(qr(matrix(seq_len(n * n) / (n * n + 1), n, n)))[, seq_len(k), drop = FALSE]
}

test_that("coupled diagonalization match cost is zero at its exact target", {
  A <- list(cd_oracle_orthonormal(5, 3), cd_oracle_orthonormal(5, 3)[5:1, ])
  Lambda <- list(
    c(0.15, 0.35, 0.55, 0.75, 0.95),
    c(0.20, 0.40, 0.60, 0.80, 1.00)
  )
  FiUbar <- list(diag(5), diag(5))
  target <- Map(function(Ai, lambda_i) {
    diag(crossprod(Ai, Ai * lambda_i))
  }, A, Lambda)

  baseline <- manifoldalign:::compute_cost_cd(
    A, Lambda, FiUbar, mu_coupling = 0, alpha_match = 0
  )
  matched <- manifoldalign:::compute_cost_cd(
    A, Lambda, FiUbar, mu_coupling = 0,
    lambda_target = target, alpha_match = 17
  )

  expect_equal(as.numeric(matched), as.numeric(baseline), tolerance = 1e-12)
  expect_equal(attr(matched, "components")$match, 0, tolerance = 1e-24)
})

test_that("coupled diagonalization gradient matches finite differences without optional packages", {
  set.seed(240722)
  A <- list(cd_oracle_orthonormal(5, 3), cd_oracle_orthonormal(5, 3)[5:1, ])
  Lambda <- list(
    c(0.15, 0.35, 0.55, 0.75, 0.95),
    c(0.20, 0.40, 0.60, 0.80, 1.00)
  )
  FiUbar <- list(diag(5), diag(5))
  direction <- matrix(rnorm(15), 5, 3)
  direction <- direction / sqrt(sum(direction^2))

  analytical <- manifoldalign:::compute_gradient_cd(
    idx = 1,
    A = A,
    Lambda = Lambda,
    FiUbar = FiUbar,
    mu_coupling = 0
  )
  objective <- function(A1) {
    A_step <- A
    A_step[[1]] <- A1
    as.numeric(manifoldalign:::compute_cost_cd(
      A_step, Lambda, FiUbar, mu_coupling = 0
    ))
  }
  epsilon <- 1e-6
  numerical <- (
    objective(A[[1]] + epsilon * direction) -
      objective(A[[1]] - epsilon * direction)
  ) / (2 * epsilon)

  expect_equal(
    sum(analytical * direction),
    numerical,
    tolerance = 2e-6
  )
})
