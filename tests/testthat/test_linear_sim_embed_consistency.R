lse_reference_objective <- function(W, X, T, M, sigma_P, alpha_p) {
  Y <- X %*% W
  distances <- as.matrix(stats::dist(Y))^2
  P <- exp(-distances / sigma_P)
  Js <- sum(M * (P - T)^2) / (2 * sum(M))

  orthogonality <- crossprod(W) - diag(ncol(W))
  Jp <- sum(orthogonality^2) / (2 * ncol(W)^2)
  (1 - alpha_p) * Js + alpha_p * Jp
}

lse_central_difference <- function(W, objective, step = 1e-6) {
  gradient <- matrix(0, nrow(W), ncol(W))
  for (j in seq_len(ncol(W))) {
    for (i in seq_len(nrow(W))) {
      plus <- minus <- W
      plus[i, j] <- plus[i, j] + step
      minus[i, j] <- minus[i, j] - step
      gradient[i, j] <- (objective(plus) - objective(minus)) / (2 * step)
    }
  }
  gradient
}

test_that("C++ LSE agrees with an independent R L-BFGS-B reference", {
  skip_if_not_installed("RcppArmadillo")
  
  # Skip if C++ backend not available
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))
  
  set.seed(42)
  
  # Create simple toy dataset
  n <- 20
  d <- 4
  ncomp <- 2
  
  X <- matrix(rnorm(n * d), n, d)
  X <- scale(X)  # Standardize
  
  # Create simple supervised target
  labels <- rep(c("A", "B"), each = n/2)
  T <- outer(labels, labels, "==") * 1
  
  # Create mask (exclude diagonal)
  M <- matrix(1, n, n)
  diag(M) <- 0
  
  # Test parameters
  sigma_P <- 1.0
  alpha_p <- 0.1
  maxit <- 200
  tol <- 1e-8

  cpp <- manifoldalign:::linear_sim_embed_cpp(
    X, T, M, sigma_P, ncomp, alpha_p, maxit, tol
  )
  reference <- manifoldalign:::.optimize_W_reference(
    X, T, M, cpp$initial_W, sigma_P, alpha_p, maxit, tol
  )

  expect_equal(cpp$convergence, 0L, info = cpp$message)
  expect_equal(reference$convergence, 0L, info = reference$message)
  expect_equal(cpp$value, reference$final_value, tolerance = 1e-8)

  D_R <- as.matrix(dist(X %*% reference$W))
  D_cpp <- as.matrix(dist(X %*% cpp$W))
  dist_correlation <- cor(as.vector(D_R), as.vector(D_cpp))
  expect_gt(dist_correlation, 0.9999)

  class_A_idx <- which(labels == "A")
  class_B_idx <- which(labels == "B")
  within_A_R <- mean(D_R[class_A_idx, class_A_idx][upper.tri(D_R[class_A_idx, class_A_idx])])
  between_R <- mean(D_R[class_A_idx, class_B_idx])
  class_sep_R <- between_R / within_A_R
  within_A_cpp <- mean(D_cpp[class_A_idx, class_A_idx][upper.tri(D_cpp[class_A_idx, class_A_idx])])
  between_cpp <- mean(D_cpp[class_A_idx, class_B_idx])
  class_sep_cpp <- between_cpp / within_A_cpp

  expect_equal(class_sep_cpp, class_sep_R, tolerance = 1e-5)
})

test_that("LSE parity holds for rank-deficient and ill-conditioned designs", {
  skip_if_not_installed("RcppArmadillo")
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))

  set.seed(2026)
  n <- 24
  z <- matrix(rnorm(n * 4), n, 4)
  designs <- list(
    rank_deficient = scale(cbind(z[, 1], z[, 2], z[, 1] + z[, 2], 2 * z[, 1])),
    ill_conditioned = scale(cbind(
      z[, 1],
      z[, 1] + 1e-4 * z[, 2],
      z[, 1] + 1e-6 * z[, 3],
      z[, 4]
    ))
  )
  expect_lt(qr(designs$rank_deficient)$rank, ncol(designs$rank_deficient))
  expect_equal(qr(designs$ill_conditioned)$rank, ncol(designs$ill_conditioned))
  expect_gt(kappa(designs$ill_conditioned), 1e6)

  labels <- rep(1:3, each = 8)
  T <- outer(labels, labels, "==") * 1
  M <- matrix(1, n, n)
  diag(M) <- 0

  for (fixture in names(designs)) {
    X <- designs[[fixture]]
    cpp <- manifoldalign:::linear_sim_embed_cpp(
      X, T, M, sigma_P = 1, ncomp = 2, alpha_p = 0.2, maxit = 300, tol = 1e-8
    )
    reference <- manifoldalign:::.optimize_W_reference(
      X, T, M, cpp$initial_W,
      sigma_P = 1, alpha_p = 0.2, maxit = 300, tol = 1e-8
    )

    expect_equal(cpp$convergence, 0L, info = fixture)
    expect_equal(reference$convergence, 0L, info = fixture)
    expect_equal(cpp$value, reference$final_value, tolerance = 1e-8, info = fixture)
    expect_true(all(is.finite(cpp$W)), info = fixture)
    expect_true(is.finite(cpp$gradient_norm), info = fixture)
  }
})

test_that("C++ objective diagnostics remain finite at an iteration limit", {
  skip_if_not_installed("RcppArmadillo")
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))
  
  set.seed(123)
  
  # Small test case for gradient verification
  n <- 10
  d <- 3  
  ncomp <- 2
  
  X <- matrix(rnorm(n * d), n, d)
  X <- scale(X)
  
  T <- matrix(runif(n * n), n, n)
  T <- 0.5 * (T + t(T))  # Make symmetric
  diag(T) <- 1
  
  M <- matrix(1, n, n)
  diag(M) <- 0
  
  sigma_P <- 0.5
  alpha_p <- 0.1
  
  # Test with same initial weights
  W_init <- matrix(rnorm(d * ncomp), d, ncomp)
  W_init <- qr.Q(qr(W_init))  # Orthogonalize
  
  # Helper function to compute objective using R formulation
  compute_objective_R <- function(W) {
    Y <- X %*% W
    D2 <- outer(rowSums(Y^2), rowSums(Y^2), `+`) - 2 * tcrossprod(Y)
    D2[D2 < 0] <- 0
    P <- exp(-D2 / sigma_P)
    P <- 0.5 * (P + t(P))
    
    Js <- sum(M * (P - T)^2) / (2 * sum(M))
    
    WtW <- crossprod(W)
    WtWmI <- WtW
    diag(WtWmI) <- diag(WtWmI) - 1
    Jp <- sum(WtWmI^2) / (2 * ncomp^2)
    
    (1 - alpha_p) * Js + alpha_p * Jp
  }
  
  # Test that initial objectives are similar
  obj_R <- compute_objective_R(W_init)
  
  result_cpp <- manifoldalign:::linear_sim_embed_cpp(
    X, T, M, sigma_P, ncomp, alpha_p, maxit = 1, tol = 1e-8
  )

  expect_true(is.finite(obj_R), info = "R objective should be finite")
  expect_true(is.finite(result_cpp$value),
              info = "C++ objective should be finite")
  expect_equal(result_cpp$convergence, 1L)
  expect_gte(result_cpp$function_evaluations, 1L)
})

test_that("R and C++ LSE gradients match an independent finite-difference oracle", {
  skip_if_not_installed("RcppArmadillo")
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))

  set.seed(321)
  X <- scale(matrix(rnorm(8 * 3), 8, 3))
  labels <- rep(c("A", "B"), each = 4)
  T <- outer(labels, labels, "==") * 1
  M <- matrix(1, 8, 8)
  diag(M) <- 0
  sigma_P <- 1.3
  alpha_p <- 0.35

  cpp <- manifoldalign:::linear_sim_embed_cpp(
    X, T, M, sigma_P, ncomp = 2, alpha_p = alpha_p, maxit = 1, tol = 1e-8
  )
  W <- cpp$initial_W
  objective <- function(candidate) {
    lse_reference_objective(candidate, X, T, M, sigma_P, alpha_p)
  }
  numerical <- lse_central_difference(W, objective)
  analytic_R <- (1 - alpha_p) * manifoldalign:::grad_Js(
    X, W, T, M, sigma_P, sum_M_batch = sum(M)
  ) + alpha_p * manifoldalign:::grad_Jp(W)

  expect_equal(cpp$initial_value, objective(W), tolerance = 1e-12)
  expect_equal(cpp$initial_gradient, numerical, tolerance = 2e-6)
  expect_equal(analytic_R, numerical, tolerance = 2e-6)
})

test_that("LSE maps optimizer status codes and rejects unresolved non-convergence", {
  expect_equal(manifoldalign:::.lse_cpp_status(0L), "converged")
  expect_equal(manifoldalign:::.lse_cpp_status(1L), "iteration_limit")
  expect_equal(manifoldalign:::.lse_cpp_status(51L), "optimizer_warning")
  expect_equal(manifoldalign:::.lse_cpp_status(52L), "optimizer_error")
  expect_equal(manifoldalign:::.lse_cpp_status(999L), "unknown_failure")

  set.seed(99)
  X <- matrix(rnorm(12 * 3), 12, 3)
  labels <- rep(c("A", "B", "C"), each = 4)
  T <- outer(labels, labels, "==") * 1

  expect_error(
    suppressWarnings(
      linear_sim_embed(
        X,
        T = T,
        sigma_P = 1,
        ncomp = 2,
        maxit = 1,
        tol = 1e-12,
        use_cpp = TRUE
      )
    ),
    "failed to converge"
  )
})

test_that("LSE rejects invalid structural inputs before optimization", {
  set.seed(7)
  X <- matrix(rnorm(10 * 3), 10, 3)
  T <- diag(10)
  M <- matrix(1, 10, 10)
  diag(M) <- 0

  asymmetric_T <- T
  asymmetric_T[1, 2] <- 0.5
  expect_error(
    linear_sim_embed(X, T = asymmetric_T, M = M, sigma_P = 1),
    "T must be symmetric"
  )

  asymmetric_M <- M
  asymmetric_M[1, 2] <- 0
  expect_error(
    linear_sim_embed(X, T = T, M = asymmetric_M, sigma_P = 1),
    "M must be symmetric"
  )

  expect_error(
    linear_sim_embed(cbind(X, constant = 1), T = T, M = M, sigma_P = 1),
    "constant column"
  )
})

test_that("LSE exposes a successful exact-objective fallback", {
  set.seed(88)
  X <- matrix(rnorm(20 * 4), 20, 4)
  labels <- rep(c("A", "B"), each = 10)
  T <- outer(labels, labels, "==") * 1

  with_mocked_bindings(
    linear_sim_embed_cpp = function(...) {
      list(convergence = 52L, message = "synthetic line-search failure")
    },
    {
      expect_warning(
        fit <- linear_sim_embed(
          X,
          T = T,
          sigma_P = 1,
          ncomp = 2,
          alpha_p = 0.2,
          maxit = 200,
          tol = 1e-8,
          use_cpp = TRUE
        ),
        "Retrying with the R optimizer"
      )
      expect_equal(fit$convergence$convergence, 0L)
      expect_equal(fit$convergence$status, "fallback_converged")
      expect_equal(fit$convergence$backend, "r")
      expect_equal(fit$convergence$fallback_from$convergence, 52L)
      expect_equal(fit$convergence$fallback_from$status, "optimizer_error")
      expect_equal(fit$optimizer, "L-BFGS-B (R fallback)")
    },
    .package = "manifoldalign"
  )
})

test_that("Alpha parameter scaling works consistently", {
  skip_if_not_installed("RcppArmadillo") 
  skip_if(!exists("linear_sim_embed_cpp", envir = asNamespace("manifoldalign"), mode = "function"))
  
  set.seed(456)
  
  # Test different alpha values
  n <- 15
  d <- 3
  ncomp <- 2
  
  X <- matrix(rnorm(n * d), n, d)
  X <- scale(X)
  
  labels <- rep(c("A", "B", "C"), each = n/3)
  T <- outer(labels, labels, "==") * 1
  M <- matrix(1, n, n)
  diag(M) <- 0
  
  sigma_P <- 1.0
  # ADAM must optimize the measured objective, not stop on a stale value.
  maxit <- 500
  
  # Test with different alpha values
  alpha_values <- c(0.01, 0.1, 0.5, 0.9)
  
  for (alpha in alpha_values) {
    result_R <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                ncomp = ncomp, alpha_p = alpha,
                                maxit = maxit, use_cpp = FALSE, verbose = FALSE)
    
    result_cpp <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                  ncomp = ncomp, alpha_p = alpha,
                                  maxit = maxit, use_cpp = TRUE, verbose = FALSE)
    
    # Both should converge
    expect_equal(result_R$convergence$convergence, 0,
                 info = sprintf("R should converge with alpha=%.2f", alpha))
    expect_equal(result_cpp$convergence$convergence, 0,
                 info = sprintf("C++ should converge with alpha=%.2f", alpha))
    expect_equal(result_cpp$convergence$status, "converged")
    expect_equal(result_cpp$convergence$backend, "cpp")
    
    # Orthogonality should increase with higher alpha
    WtW_R <- crossprod(result_R$weights)
    WtW_cpp <- crossprod(result_cpp$weights)
    
    orthog_error_R <- norm(WtW_R - diag(ncomp), "F")
    orthog_error_cpp <- norm(WtW_cpp - diag(ncomp), "F")
    
    expect_true(orthog_error_R < 1.0, 
                info = sprintf("R orthogonality error too large with alpha=%.2f: %.3f", 
                              alpha, orthog_error_R))
    expect_true(orthog_error_cpp < 1.0,
                info = sprintf("C++ orthogonality error too large with alpha=%.2f: %.3f", 
                              alpha, orthog_error_cpp))
  }
})
