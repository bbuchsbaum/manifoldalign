# test-kema.R  —  unit tests for the KEMA implementation
# -------------------------------------------------------------------------
# These tests assume that the package exports
#   kema.hyperdesign(), coskern(), center(), and eigencore are available.
#
# All three tests run with n ≤ 60 and avoid any large eigendecompositions.
#
# NOTE: These tests use the kema() generic function which dispatches to kema.hyperdesign()
# All tests should now run successfully with proper function exports.
#
# -------------------------------------------------------------------------
library(testthat)
library(Matrix)
library(multidesign)
library(manifoldalign)     # Load the manifoldalign package
set.seed(42)

# -------------------------------------------------------------------------
# helper that fabricates a minimal hyper-design object
# -------------------------------------------------------------------------
quick_hd <- function(Xlist, ylist) {
  md_list <- Map(function(x, y) {
    multidesign::multidesign(x, data.frame(lbl = factor(y)))
  }, Xlist, ylist)
  
  # Create names for the domains
  names(md_list) <- paste0("domain", seq_along(md_list))
  
  multidesign::hyperdesign(md_list)
}

# -------------------------------------------------------------------------
# helper for rotation-invariant subspace comparison
# -------------------------------------------------------------------------
principal_angle_distance <- function(A, B) {
  # Orthonormalize both matrices
  qa <- qr.Q(qr(A))
  qb <- qr.Q(qr(B))
  # Singular values of QA^T QB are cosines of principal angles
  sv <- svd(t(qa) %*% qb, nu = 0, nv = 0)$d
  max_angle <- acos(min(sv))  # Largest principal angle in radians
  sin(max_angle)              # Distance in [0,1], where 0 = perfect match
}

# -------------------------------------------------------------------------
# 1.  Exact solver really solves the generalized eigen-problem
#     (Eq. 6 on page 7 of Tuia & Camps-Valls 2016)
# -------------------------------------------------------------------------
test_that("exact solver returns true generalized eigenvectors", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  n  <- 30;  p <- 4
  X1 <- matrix(rnorm(n * p), n)
  X2 <- matrix(rnorm(n * p) + 0.6, n)
  y  <- rep(1:3, length.out = n)

  hd <- quick_hd(list(X1, X2), list(y, y))
  fit <- kema(hd, y = lbl,
                          ncomp        = 2,
                          kernel       = kernlab::vanilladot(), # linear ⇒ match SSMA
                          solver       = "exact",
                          sample_frac  = 1,
                          u            = .5,
                          lambda       = 1e-3)

  # The fit object is a multiblock_biprojector, so we need to access scores differently
  # fit$s contains the scores matrix, fit$alpha contains the eigenvectors
  expect_true(is.matrix(fit$s) || methods::is(fit$s, "Matrix"))
  expect_true(is.matrix(fit$alpha) || methods::is(fit$alpha, "Matrix"))
  expect_equal(ncol(fit$s), 2)  # ncomp = 2
  expect_equal(ncol(fit$alpha), 2)  # ncomp = 2
  
  expect_true(fit$fidelity$passed)
  expect_lte(fit$fidelity$max_rel_residual, 1e-6)
  expect_lte(fit$fidelity$max_b_orth_offdiag, 1e-6)
  expect_true(all(
    fit$eigenvalues$values > fit$fidelity$eigenvalue_zero_tol
  ))
  expect_identical(fit$eigenvalues$solver, "exact")
})

# -------------------------------------------------------------------------
# 2.  Withdrawn compatibility controls cannot masquerade as fitted behavior
# -------------------------------------------------------------------------
test_that("withdrawn KEMA controls are rejected", {
  n  <- 20; p <- 4
  X1 <- matrix(rnorm(n * p), n)
  X2 <- matrix(rnorm(n * p) - 0.3, n)
  y  <- rep(1:2, length.out = n)

  hd <- quick_hd(list(X1, X2), list(y, y))

  expect_error(
    kema(hd, y = lbl, solver = "regression"),
    "solver='regression' has been withdrawn"
  )
  expect_error(
    kema(hd, y = lbl, dweight = 0),
    "Unsupported KEMA argument.*dweight"
  )
  expect_error(
    kema(hd, y = lbl, centre_kernel = TRUE),
    "Unsupported KEMA argument.*centre_kernel"
  )
  expect_error(
    kema(hd, y = lbl, imaginary_control = 1),
    "Unsupported KEMA argument.*imaginary_control"
  )
})

# -------------------------------------------------------------------------
# 3.  Landmark solver (REKEMA) is consistent with the full solver
# -------------------------------------------------------------------------
test_that("REKEMA scores closely match full KEMA scores", {
  skip_if_not_installed("eigencore", minimum_version = "1.0.3")

  n  <- 60;  p <- 6
  X1 <- matrix(rnorm(n * p), n)
  X2 <- matrix(rnorm(n * p) + 1, n)
  y  <- rep(1:3, length.out = n)

  hd <- quick_hd(list(X1, X2), list(y, y))

  full <- kema(hd, y = lbl,
                             ncomp       = 3,
                             sample_frac = 1,
                             solver      = "exact",
                             kernel      = kernlab::rbfdot(sigma = 0.5))

  set.seed(42)
  rekem <- kema(hd, y = lbl,
                             ncomp       = 3,
                             sample_frac = .5,
                             solver      = "exact",
                             kernel      = kernlab::rbfdot(sigma = 0.5))

  # Subspace similarity: principal angles between score spaces
  # Check actual correlation values - REKEMA is an approximation method
  U <- qr.Q(qr(as.matrix(full$s)))
  V <- qr.Q(qr(as.matrix(rekem$s)))
  d <- svd(crossprod(U, V))$d

  expect_true(full$fidelity$passed)
  expect_true(rekem$fidelity$passed)
  expect_lte(full$fidelity$max_rel_residual, 1e-6)
  expect_lte(rekem$fidelity$max_rel_residual, 1e-6)
  expect_gt(min(d[1:2]), 0.7)
  expect_gt(d[[3L]], 0.6)
})
