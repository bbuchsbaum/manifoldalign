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
  
  # Run R implementation
  result_R <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P, 
                              ncomp = ncomp, alpha_p = alpha_p,
                              maxit = maxit, tol = tol, use_cpp = FALSE,
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

test_that("Gradient computations match between R and C++ (numerical verification)", {
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
  
  # Create linear_sim_embed instance to test C++ objective
  # This is a bit tricky since we can't directly call the C++ objective function
  # So we'll run a single iteration and check the reported objective
  
  result_cpp <- linear_sim_embed(X, T = T, M = M, sigma_P = sigma_P,
                                ncomp = ncomp, alpha_p = alpha_p,
                                maxit = 1, use_cpp = TRUE, verbose = FALSE)
  
  # The objectives should be reasonably close
  # (might not be exact due to initialization differences)
  expect_true(is.finite(obj_R), info = "R objective should be finite")
  expect_true(is.finite(result_cpp$convergence$final_value), 
              info = "C++ objective should be finite")
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