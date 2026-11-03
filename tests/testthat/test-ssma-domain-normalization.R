.domain_reduced <- function(A, X) list(A_r = A, ranks = vapply(X, ncol, integer(1)), basis = lapply(X, function(x) list(T = x)))

test_that("domain score solver matches a planted orthogonal correspondence oracle", {
  set.seed(37)
  O <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  A <- rbind(cbind(diag(5), -O), cbind(-t(O), diag(5)))
  z <- .ssma_solve_domain_scores(.domain_reduced(A, list(diag(5), diag(5))), 3, 1e-8, 1000)
  a <- z$alpha[1:5, , drop = FALSE]; b <- z$alpha[6:10, , drop = FALSE]
  expect_equal(a, O %*% b, tolerance = 1e-7)
  expect_equal(crossprod(a), diag(3), tolerance = 1e-10)
  expect_equal(crossprod(b), diag(3), tolerance = 1e-10)
  expect_lt(abs(sum(z$alpha * (A %*% z$alpha))), 1e-10)
  expect_true(z$domain_diagnostics$converged)
  expect_true(all(diff(z$domain_diagnostics$objective_history) <= 1e-10))
  expect_null(z$eigenvalues)
})

test_that("domain normalization prevents collapse for disconnected unequal domains", {
  X <- list(diag(c(2, 1, .5)), diag(c(3, 2, 1, .5)))
  A <- diag(c(100, 80, 60, 1, 3, 2, 4))
  z <- .ssma_solve_domain_scores(.domain_reduced(A, X), 2, 1e-8, 1000)
  b <- .ssma_split_blocks(z$alpha, c(3, 4))
  for (i in 1:2) expect_equal(crossprod(X[[i]] %*% b[[i]]), diag(2), tolerance = 1e-9)
  # Exact independent oracle: choose each domain's two smallest whitened costs.
  oracle <- sum(sort(diag(A)[1:3] / diag(X[[1]])^2)[1:2]) + sum(sort(diag(A)[4:7] / diag(X[[2]])^2)[1:2])
  expect_equal(z$domain_diagnostics$objective, oracle, tolerance = 1e-8)
})

test_that("domain solver handles score-null coordinates and rejects infeasible rank", {
  X <- list(cbind(diag(3), 0), diag(4))
  z <- .ssma_solve_domain_scores(.domain_reduced(diag(8), X), 2, 1e-7, 1000)
  expect_equal(z$domain_diagnostics$retained_score_ranks, c(3L, 4L))
  expect_error(.ssma_solve_domain_scores(.domain_reduced(diag(8), X), 4, 1e-7, 1000), "independent score directions")
})

test_that("domain optimizer descends the independent dense quadratic and rejects nonconvergence", {
  set.seed(71)
  U <- matrix(rnorm(81), 9, 9); A <- crossprod(U)
  X <- list(matrix(rnorm(36), 12, 3), matrix(rnorm(54), 9, 6))
  z <- .ssma_solve_domain_scores(.domain_reduced(A, X), 2, 1e-7, 5000)
  expect_equal(z$domain_diagnostics$objective, sum(z$alpha * (A %*% z$alpha)), tolerance = 1e-10)
  expect_lte(z$domain_diagnostics$objective, z$domain_diagnostics$initial_objective)
  expect_true(all(diff(z$domain_diagnostics$objective_history) <= 1e-10))
  expect_error(.ssma_solve_domain_scores(.domain_reduced(A, X), 2, 1e-14, 1), "did not converge")
})

test_that("public domain score SSMA composes sparse linear constraints and preserves predictions", {
  set.seed(17)
  x <- list(a = matrix(rnorm(180), 30, 6), b = matrix(rnorm(228), 38, 6))
  h <- lapply(x, function(z) Matrix::Matrix(diag(6), sparse = TRUE))
  link <- linear_correspondence(names(x), h, strength = 10)
  control <- ssma_align_control(normalization = 'domain_scores', compression = 'none', auto_correspondences = FALSE, sigma = 1, domain_tol = 1e-6)
  fit <- ssma_align(x, constraints = list(link), preproc = NULL, ncomp = 2, control = control)
  expect_s3_class(fit, 'ssma_align')
  expect_identical(fit$normalization, 'domain_scores')
  expect_null(fit$generalized_eigenvalues)
  scores <- list(fit$s[1:30, , drop = FALSE], fit$s[31:68, , drop = FALSE])
  for (i in 1:2) expect_equal(crossprod(scores[[i]]), diag(2), tolerance = 1e-8)
  expect_equal(x$a %*% fit$v[1:6, , drop = FALSE], scores[[1]], tolerance = 1e-10)
  expect_true(fit$domain_diagnostics$converged)
  expect_equal(length(fit$constraint_diagnostics), 1)
  expect_error(ssma_align_control(normalization = 'domain_scores', solver = 'operator'), 'does not support')
  expect_error(ssma_align_control(domain_tol = NA_real_), 'positive finite')
  expect_error(ssma_align_control(domain_maxiter = 1.5), 'positive integer')
})

test_that("domain score result respects feature reparameterization and domain identity", {
  set.seed(101)
  x <- list(a = matrix(rnorm(150), 30, 5), b = matrix(rnorm(175), 35, 5))
  h <- list(matrix(rnorm(40), 8, 5), matrix(rnorm(40), 8, 5))
  ctl <- ssma_align_control(normalization = 'domain_scores', compression = 'none', auto_correspondences = FALSE, sigma = 1, domain_tol = 1e-7)
  link <- linear_correspondence(names(x), h, strength = 3)
  fit <- ssma_align(x, constraints = list(link), preproc = NULL, ncomp = 2, control = ctl)
  # Domain permutation must preserve the feasible objective; common rotation is arbitrary.
  revfit <- ssma_align(rev(x), constraints = list(link), preproc = NULL, ncomp = 2, control = ctl)
  expect_equal(fit$domain_diagnostics$objective, revfit$domain_diagnostics$objective, tolerance = 1e-5)
  expect_identical(revfit$domain_names, rev(names(x)))
  O <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  rotated <- lapply(x, function(z) z %*% O)
  rotated_link <- linear_correspondence(names(x), lapply(h, function(z) z %*% O), strength = 3)
  transformed <- ssma_align(rotated, constraints = list(rotated_link), preproc = NULL, ncomp = 2, control = ctl)
  expect_equal(fit$domain_diagnostics$objective, transformed$domain_diagnostics$objective, tolerance = 1e-5)
  expect_lt(transformed$domain_diagnostics$score_gram_max_error, 1e-8)
})


test_that("domain tangent gradient matches independent retraction finite differences", {
  set.seed(89)
  ii <- list(1:4, 5:10)
  Y <- rbind(.ssma_retract(matrix(rnorm(8), 4, 2)), .ssma_retract(matrix(rnorm(12), 6, 2)))
  A <- matrix(rnorm(100), 10, 10); A <- crossprod(A)
  G <- .ssma_domain_gradient(Y, 2 * A %*% Y, ii)
  D <- .ssma_domain_gradient(Y, matrix(rnorm(20), 10, 2), ii)
  for (j in ii) expect_equal(crossprod(Y[j, ], G[j, ]) + crossprod(G[j, ], Y[j, ]), matrix(0, 2, 2), tolerance = 1e-10)
  value <- function(t) {
    Z <- do.call(rbind, lapply(ii, function(j) .ssma_retract((Y + t * D)[j, , drop = FALSE])))
    sum(Z * (A %*% Z))
  }
  expect_equal((value(1e-6) - value(-1e-6)) / 2e-6, sum(G * D), tolerance = 1e-6)
})


test_that("Stiefel QR retraction is scale safe and has the tangent derivative", {
  set.seed(126)
  Y <- .ssma_retract(matrix(rnorm(64 * 30), 64, 30))
  D <- .ssma_domain_gradient(Y, matrix(rnorm(64 * 30), 64, 30), list(1:64))
  expect_equal(.ssma_retract(Y), Y, tolerance = 1e-10)
  expect_equal(crossprod(.ssma_retract(1e150 * (Y + D))), diag(30), tolerance = 1e-10)
  expect_equal((.ssma_retract(Y + 1e-6 * D) - .ssma_retract(Y - 1e-6 * D)) / 2e-6, D, tolerance = 1e-7)
  expect_error(.ssma_retract(matrix(Inf, 4, 2)), 'Nonfinite')
})
