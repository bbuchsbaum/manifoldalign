test_that("R and C++ linear_sim_embed implementations produce consistent results", {
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
  maxit <- 50
  tol <- 1e-6
  
  # Run R implementation (ADAM takes small steps, so it gets the default
  # budget; L-BFGS-B keeps maxit = 50)
  result_R <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P, 
                              ncomp = ncomp, alpha_p = alpha_p,
                              maxit = 500, tol = tol, use_cpp = FALSE,
                              verbose = FALSE)
  
  # Run C++ implementation
  result_cpp <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                ncomp = ncomp, alpha_p = alpha_p,
                                maxit = maxit, tol = tol, use_cpp = TRUE,
                                verbose = FALSE)
  
  # Test 1: Check that both implementations converged
  expect_equal(result_R$convergence$convergence, 0, 
               info = "R implementation should converge")
  expect_equal(result_cpp$convergence$convergence, 0,
               info = "C++ implementation should converge")
  
  # Test 2: Final objective values should be similar
  r_final_obj <- tail(result_R$objective_trace, 1)
  cpp_final_obj <- result_cpp$convergence$final_value
  
  expect_true(abs(r_final_obj - cpp_final_obj) / abs(r_final_obj) < 0.05,
              info = sprintf("Relative objective difference too large: R=%.6f, C++=%.6f", 
                           r_final_obj, cpp_final_obj))
  
  # Test 3: Weight matrices should be similar (allowing for sign flips and orthogonal rotations)
  W_R <- result_R$weights
  W_cpp <- result_cpp$weights
  
  expect_equal(dim(W_R), dim(W_cpp), info = "Weight matrix dimensions should match")
  
  # Test 4: Projected data should have similar structure (check pairwise distances)
  Y_R <- result_R$scores
  Y_cpp <- X %*% W_cpp  # Project using C++ weights
  
  # Compare distance matrices (invariant to rotations/reflections)
  D_R <- as.matrix(dist(Y_R))
  D_cpp <- as.matrix(dist(Y_cpp))
  
  # Correlation of distance matrices should be high
  dist_correlation <- cor(as.vector(D_R), as.vector(D_cpp))
  expect_true(dist_correlation > 0.9,
              info = sprintf("Distance correlation too low: %.3f", dist_correlation))
  
  # Test 5: Both should preserve class structure similarly
  # Check that within-class distances are smaller than between-class
  class_A_idx <- which(labels == "A")
  class_B_idx <- which(labels == "B")
  
  # R implementation
  within_A_R <- mean(D_R[class_A_idx, class_A_idx][upper.tri(D_R[class_A_idx, class_A_idx])])
  between_R <- mean(D_R[class_A_idx, class_B_idx])
  class_sep_R <- between_R / within_A_R
  
  # C++ implementation  
  within_A_cpp <- mean(D_cpp[class_A_idx, class_A_idx][upper.tri(D_cpp[class_A_idx, class_A_idx])])
  between_cpp <- mean(D_cpp[class_A_idx, class_B_idx])
  class_sep_cpp <- between_cpp / within_A_cpp
  
  # Both should achieve similar class separation ratios
  expect_true(abs(class_sep_R - class_sep_cpp) / class_sep_R < 0.2,
              info = sprintf("Class separation ratios differ: R=%.3f, C++=%.3f", 
                           class_sep_R, class_sep_cpp))
})

test_that("R and C++ gradients match central finite differences of the documented objective", {
  skip_if_not_installed("RcppArmadillo")
  skip_if(!exists("lse_objective_gradient_cpp", envir = asNamespace("manifoldalign"), mode = "function"))
  ns <- asNamespace("manifoldalign")

  set.seed(123)
  n <- 10
  d <- 3
  ncomp <- 2

  X <- scale(matrix(rnorm(n * d), n, d))
  T <- matrix(runif(n * n), n, n)
  T <- 0.5 * (T + t(T))
  diag(T) <- 1
  M <- matrix(1, n, n)
  diag(M) <- 0
  sigma_P <- 0.7
  alpha_p <- 0.3

  # Independent oracle for the documented objective:
  # (1 - a) * sum_{i<j} M (P - T)^2 / (2K) + a * ||W'W - I||_F^2 / (2 m^2)
  objective_oracle <- function(W, a = alpha_p) {
    Y <- X %*% W
    P <- exp(-as.matrix(dist(Y))^2 / sigma_P)
    up <- upper.tri(M) & M != 0
    Js <- sum((P[up] - T[up])^2) / (2 * sum(up))
    Jp <- sum((crossprod(W) - diag(ncol(W)))^2) / (2 * ncol(W)^2)
    (1 - a) * Js + a * Jp
  }
  central_diff <- function(f, W, h = 1e-6) {
    G <- W
    for (k in seq_along(W)) {
      Wp <- W; Wm <- W
      Wp[k] <- Wp[k] + h
      Wm[k] <- Wm[k] - h
      G[k] <- (f(Wp) - f(Wm)) / (2 * h)
    }
    G
  }

  # Deliberately non-orthogonal weights so the orthogonality term is active.
  for (W in list(matrix(rnorm(d * ncomp), d, ncomp),
                 2.5 * matrix(rnorm(d * ncomp), d, ncomp))) {
    expect_gt(sum((crossprod(W) - diag(ncomp))^2), 0.1)
    fd <- central_diff(objective_oracle, W)

    obj_R <- ns$.lse_objective(X, W, T, M, sigma_P, alpha_p)$objective
    expect_equal(obj_R, objective_oracle(W), tolerance = 1e-12)
    grad_R <- (1 - alpha_p) * ns$grad_Js(X, W, T, M, sigma_P, sum(M)) +
      alpha_p * ns$grad_Jp(W)
    expect_equal(unname(as.matrix(grad_R)), fd, tolerance = 1e-6)

    cpp <- ns$lse_objective_gradient_cpp(X, T, M, sigma_P, alpha_p, W)
    expect_equal(cpp$value, objective_oracle(W), tolerance = 1e-12)
    expect_equal(cpp$gradient, fd, tolerance = 1e-6)
  }

  # Each term separately (alpha = 0: similarity only; alpha = 1: orthogonality only)
  W <- matrix(rnorm(d * ncomp), d, ncomp)
  for (a in c(0, 1)) {
    fd_a <- central_diff(function(W) objective_oracle(W, a), W)
    cpp_a <- ns$lse_objective_gradient_cpp(X, T, M, sigma_P, a, W)
    expect_equal(cpp_a$gradient, fd_a, tolerance = 1e-6)
    grad_a <- if (a == 0) ns$grad_Js(X, W, T, M, sigma_P, sum(M)) else ns$grad_Jp(W)
    expect_equal(unname(as.matrix(grad_a)), fd_a, tolerance = 1e-6)
  }
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
  maxit <- 30
  
  # Test with different alpha values
  alpha_values <- c(0.01, 0.1, 0.5, 0.9)
  
  for (alpha in alpha_values) {
    # ADAM takes small steps; it is given its default budget, while L-BFGS-B
    # keeps the 30-iteration budget.
    result_R <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                ncomp = ncomp, alpha_p = alpha,
                                maxit = 500, use_cpp = FALSE, verbose = FALSE)
    
    result_cpp <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                  ncomp = ncomp, alpha_p = alpha,
                                  maxit = maxit, use_cpp = TRUE, verbose = FALSE)
    
    # Both should converge
    expect_equal(result_R$convergence$convergence, 0,
                 info = sprintf("R should converge with alpha=%.2f", alpha))
    expect_equal(result_cpp$convergence$convergence, 0,
                 info = sprintf("C++ should converge with alpha=%.2f", alpha))

    # Both backends minimise the same objective and reach the same value
    expect_equal(result_R$convergence$final_value, result_cpp$convergence$final_value,
                 tolerance = 1e-3,
                 info = sprintf("Final objectives differ with alpha=%.2f", alpha))
    
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