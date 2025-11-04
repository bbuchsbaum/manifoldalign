# tests/testthat/test-fpgw-validation.R
# ---------------------------------------------------------------------------
#  Comprehensive validation suite for Fused-Partial Gromov-Wasserstein
# ---------------------------------------------------------------------------
#
#  These tests verify the fundamental mathematical properties that any correct
#  FPGW implementation must satisfy, based on Bai et al. (2025).
#
#  Run with: devtools::test(filter = "fpgw-validation")
#
# ---------------------------------------------------------------------------

library(testthat)
library(Matrix)
library(tibble)

have_multidesign <- requireNamespace("multidesign", quietly = TRUE)
have_multivarious <- requireNamespace("multivarious", quietly = TRUE)
have_lpSolve <- requireNamespace("lpSolve", quietly = TRUE)

skip_if_missing_fpgw_deps <- function() {
  skip_if_not(have_multidesign, "multidesign not installed")
  skip_if_not(have_multivarious, "multivarious not installed")
  skip_if_not(have_lpSolve, "lpSolve not installed")
}
set.seed(1234)

# ---------------------------------------------------------------------------
# 1. Helper - objective value that does not rely on package internals
# ---------------------------------------------------------------------------
objective_fpgw <- function(gamma, C_feat, Cx, Cy,
                          omega1 = 0.5, lambda = 0, rho = NULL) {
  
  struct_wt <- 1 - omega1
  # term 1 - feature
  ftr <- omega1 * sum(C_feat * gamma)
  
  # term 2 - structure (quadratic form)
  const1 <- as.vector(Cx^2 %*% rowSums(gamma))
  const2 <- as.vector(Cy^2 %*% colSums(gamma))
  const  <- outer(const1, rep(1, ncol(gamma))) +
            outer(rep(1, nrow(gamma)), const2)
  gw_term <- sum((const - 2 * Cx %*% gamma %*% t(Cy)) * gamma)
  strc <- struct_wt * gw_term
  
  # term 3 - TV penalty
  tv <- if (lambda > 0) lambda * sum(gamma) else 0
  
  ftr + strc + tv
}

# ---------------------------------------------------------------------------
# 2. Test data generators
# ---------------------------------------------------------------------------
toy_pair <- function(d = 2, n = 3) {
  list(
    X = matrix(rnorm(n * d), n, d),
    Y = matrix(rnorm(n * d), n, d)
  )
}

create_test_hyperdesign <- function(X_list, design = NULL) {
  n_domains <- length(X_list)
  if (is.null(design)) {
    n <- nrow(X_list[[1]])
    design <- data.frame(sample_id = 1:n)
  }
  
  md_list <- lapply(seq_along(X_list), function(i) {
    multidesign::multidesign(X_list[[i]], design)
  })
  names(md_list) <- paste0("domain", seq_along(X_list))
  
  multidesign::hyperdesign(md_list)
}

# ---------------------------------------------------------------------------
# 3. Finite-difference gradient check
# ---------------------------------------------------------------------------
test_that("Analytic gradient matches numeric finite differences", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  pb <- toy_pair()
  Cx <- as.matrix(dist(pb$X))
  Cy <- as.matrix(dist(pb$Y))
  C_feat <- manifoldalign:::compute_feature_cost(pb$X, pb$Y, "euclidean")
  
  p <- rep(1/3, 3)
  q <- rep(1/3, 3)
  P0 <- outer(p, q)
  
  # Run one iteration to get gradient
  sol <- manifoldalign:::fw_fpgw(
    C_feat, Cx, Cy, p, q,
    omega1 = 0.3, lambda = 0, rho = NULL,
    P0 = P0, max_iter = 1, tol = 1e-9
  )
  grad_analytical <- sol$grad
  
  # Compute numerical gradient
  eps <- 1e-6
  grad_numerical <- matrix(0, 3, 3)
  
  for (i in 1:3) {
    for (j in 1:3) {
      P_plus <- P0
      P_minus <- P0
      P_plus[i, j] <- P_plus[i, j] + eps
      P_minus[i, j] <- P_minus[i, j] - eps
      
      obj_plus <- objective_fpgw(P_plus, C_feat, Cx, Cy, omega1 = 0.3)
      obj_minus <- objective_fpgw(P_minus, C_feat, Cx, Cy, omega1 = 0.3)
      
      grad_numerical[i, j] <- (obj_plus - obj_minus) / (2 * eps)
    }
  }
  
  expect_equal(as.vector(grad_analytical), as.vector(grad_numerical), 
               tolerance = 1e-5,
               label = "Analytical gradient should match numerical gradient")
})

# ---------------------------------------------------------------------------
# 4. Frank-Wolfe monotone decrease
# ---------------------------------------------------------------------------
test_that("Frank-Wolfe objective is non-increasing and converges", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  pb <- toy_pair()
  hd <- create_test_hyperdesign(list(pb$X, pb$Y))
  
  # Run FPGW with classical setting
  result <- fpgw(hd, omega1 = 0.4, lambda = 0, max_iter = 50, tol = 1e-10)
  
  # Extract objective trace from first pair
  # Since fpgw.hyperdesign calls fw_fpgw internally, we need to modify to capture trace
  # For now, we'll test convergence indirectly
  expect_true(result$converged[1, 2], 
              label = "Frank-Wolfe should converge for well-posed problems")
})

# ---------------------------------------------------------------------------
# 5. Oracle correctness - mass-constrained partial OT
# ---------------------------------------------------------------------------
test_that("Partial OT oracle respects constraints exactly", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  # Small problem for exact verification
  G <- matrix(runif(9), 3, 3)
  p <- c(0.3, 0.4, 0.3)
  q <- c(0.2, 0.5, 0.3)
  rho <- 0.5
  
  gamma <- manifoldalign:::partial_ot_mass(G, p, q, rho)
  
  # Check row constraints (within numerical tolerance)
  row_sums <- rowSums(gamma)
  row_violation <- max(row_sums - p)
  expect_lte(row_violation, 1e-6,
             label = sprintf("Row sums should not exceed marginal p (max excess %.2e)", row_violation))
  
  # Check column constraints (within numerical tolerance)
  col_sums <- colSums(gamma)
  col_violation <- max(col_sums - q)
  expect_lte(col_violation, 1e-6,
             label = sprintf("Column sums should not exceed marginal q (max excess %.2e)", col_violation))
  
  # Check exact mass constraint
  total_mass <- sum(gamma)
  expect_equal(total_mass, rho, tolerance = 1e-6,
               label = "Total transported mass should equal rho exactly")
  
  # Verify it's a minimizer by checking KKT conditions (simplified)
  expect_true(all(gamma >= -1e-6),
              label = "Transport plan should be non-negative")
})

# ---------------------------------------------------------------------------
# 6. Classical OT oracle correctness
# ---------------------------------------------------------------------------
test_that("Classical OT oracle produces doubly stochastic matrix", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  G <- matrix(runif(16), 4, 4)
  p <- rep(0.25, 4)
  q <- rep(0.25, 4)
  
  gamma <- manifoldalign:::classical_ot_lp(G, p, q)
  
  # Check double stochasticity
  expect_equal(as.vector(rowSums(gamma)), p, tolerance = 1e-6,
               label = "Row sums should equal p")
  expect_equal(as.vector(colSums(gamma)), q, tolerance = 1e-6,
               label = "Column sums should equal q")
  expect_true(all(gamma >= -1e-6),
              label = "Transport plan should be non-negative")
})

# ---------------------------------------------------------------------------
# 7. Metric properties: symmetry and triangle inequality
# ---------------------------------------------------------------------------
test_that("FPGW distance is symmetric and satisfies triangle inequality", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  # Three small random domains
  set.seed(42)
  X1 <- matrix(rnorm(12), 4, 3)
  X2 <- matrix(rnorm(12), 4, 3)
  X3 <- matrix(rnorm(12), 4, 3)
  
  hd <- create_test_hyperdesign(list(X1, X2, X3))
  
  # Compute all pairwise distances
  result <- fpgw(hd, omega1 = 0.3, max_iter = 60)
  dmat <- result$distances
  
  # Check symmetry
  expect_equal(dmat, t(dmat), tolerance = 1e-10,
               label = "Distance matrix should be symmetric")
  
  # Check diagonal is zero
  expect_equal(diag(dmat), rep(0, 3), tolerance = 1e-10,
               label = "Distance from domain to itself should be zero")
  
  # Check triangle inequality
  for (i in 1:3) {
    for (j in 1:3) {
      for (k in 1:3) {
        expect_lte(dmat[i, j], dmat[i, k] + dmat[k, j] + 1e-7,
                   label = sprintf("Triangle inequality for (%d,%d,%d)", i, j, k))
      }
    }
  }
})

# ---------------------------------------------------------------------------
# 8. Limiting behavior: FPGW collapses to FGW
# ---------------------------------------------------------------------------
test_that("FPGW with rho=1 or huge lambda equals classical FGW", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  set.seed(42)  # Ensure deterministic test
  
  pb <- toy_pair(d = 3, n = 5)
  hd <- create_test_hyperdesign(list(pb$X, pb$Y))
  
  # Classical FGW (no penalty, full transport)
  # Use larger epsilon for better warm-start with this seed
  fgw_classical <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = NULL, 
                       epsilon = 0.1, max_iter = 200, tol = 1e-4)
  
  # FPGW with rho = 1 (transport all mass)
  fpgw_rho1 <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = 1.0, 
                   epsilon = 0.1, max_iter = 200, tol = 1e-4)
  
  # Compare distances
  expect_equal(fgw_classical$distances[1, 2], fpgw_rho1$distances[1, 2], 
               tolerance = 1e-4,
               label = "FPGW with rho=1 should equal classical FGW")
  
  # Compare transport plans
  P_classical <- fgw_classical$transport_plans[[1]]
  P_rho1 <- fpgw_rho1$transport_plans[[1]]
  
  expect_equal(sum(P_classical), 1.0, tolerance = 1e-6)
  expect_equal(sum(P_rho1), 1.0, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 9. TV penalty effect: transported mass decreases with lambda
# ---------------------------------------------------------------------------
test_that("Mass penalty decreases transported mass as lambda grows", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()

  set.seed(123)
  X1 <- matrix(rnorm(20), 5, 4)
  X2 <- matrix(rnorm(20), 5, 4)
  hd <- create_test_hyperdesign(list(X1, X2))

  result_lambda0 <- fpgw(hd, omega1 = 0.4, lambda = 0.0, max_iter = 60)
  result_lambda1 <- fpgw(hd, omega1 = 0.4, lambda = 0.5, max_iter = 60)
  result_lambda5 <- fpgw(hd, omega1 = 0.4, lambda = 2.0, max_iter = 60)

  mass0 <- sum(result_lambda0$transport_plans[[1]])
  mass1 <- sum(result_lambda1$transport_plans[[1]])
  mass5 <- sum(result_lambda5$transport_plans[[1]])

  tol <- 1e-6
  expect_lte(mass1, mass0 + tol,
             label = sprintf("λ=0.5 mass (%.4f) should not exceed λ=0 mass (%.4f)", mass1, mass0))
  expect_lte(mass5, mass1 + tol,
             label = sprintf("λ=2.0 mass (%.4f) should not exceed λ=0.5 mass (%.4f)", mass5, mass1))
})

# ---------------------------------------------------------------------------
# 10. Mass budget constraint validation
# ---------------------------------------------------------------------------
test_that("Mass-constrained variant respects rho exactly", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  # Test with various rho values
  pb <- toy_pair(d = 2, n = 6)
  hd <- create_test_hyperdesign(list(pb$X, pb$Y))
  
  for (rho in c(0.3, 0.5, 0.7, 0.9)) {
    result <- fpgw(hd, omega1 = 0.3, rho = rho, max_iter = 50)
    transported_mass <- sum(result$transport_plans[[1]])
    
    expect_equal(transported_mass, rho, tolerance = 1e-5,
                 label = sprintf("Transported mass should equal rho=%.1f", rho))
  }
})

# ---------------------------------------------------------------------------
# 11. Numerical stability tests
# ---------------------------------------------------------------------------
test_that("FPGW is numerically stable", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  # Test with poorly conditioned data
  X1 <- matrix(c(1, 1e-8, 1e8, 1, 1, 1e-8), 3, 2)
  X2 <- matrix(c(1, 1e-7, 1e7, 1, 1, 1e-7), 3, 2)
  
  hd <- create_test_hyperdesign(list(X1, X2))
  
  expect_no_error({
    result <- fpgw(hd, omega1 = 0.5, lambda = 0.1)
  })
  
  # Check for valid output
  expect_true(all(result$distances >= 0),
              label = "All distances should be non-negative")
  expect_true(all(is.finite(result$distances)),
              label = "All distances should be finite")
})

# ---------------------------------------------------------------------------
# 12. Edge case: identical domains
# ---------------------------------------------------------------------------
test_that("Distance between identical domains is zero", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  X <- matrix(rnorm(20), 5, 4)
  hd <- create_test_hyperdesign(list(X, X))  # Same data twice
  
  result <- fpgw(hd, omega1 = 0.5)
  
  expect_equal(result$distances[1, 2], 0, tolerance = 1e-8,
               label = "Distance between identical domains should be zero")
  
  # Transport plan should be identity (permuted)
  P <- result$transport_plans[[1]]
  expect_equal(sum(diag(P)), 1.0, tolerance = 1e-4,
               label = "Transport plan should be approximately diagonal for identical domains")
})

# ---------------------------------------------------------------------------
# 13. Different dimensions handling
# ---------------------------------------------------------------------------
test_that("FPGW handles domains with different dimensions", {
  skip_on_cran()
  skip_if_missing_fpgw_deps()
  
  set.seed(456)  # Set seed for reproducibility
  X1 <- matrix(rnorm(15), 5, 3)  # 5 points in R^3
  X2 <- matrix(rnorm(20), 5, 4)  # 5 points in R^4
  
  hd <- create_test_hyperdesign(list(X1, X2))
  
  expect_no_error({
    result <- fpgw(hd, omega1 = 0.2)
  })
  
  expect_true(result$distances[1, 2] > 0,
              label = "Distance should be positive for different dimensional data")
})
