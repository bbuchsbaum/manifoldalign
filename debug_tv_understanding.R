library(manifoldalign)
library(multidesign)
library(tibble)

# Create test data with outliers
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

# Create hyperdesign
design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

cat("Understanding TV penalty behavior:\n\n")

# Test with longer iterations
results <- list()
for (lambda in c(0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0)) {
  cat("\nLambda =", lambda, "\n")
  
  # Run with more iterations
  result <- fpgw(hd, omega1 = 0.4, lambda = lambda, max_iter = 100, 
                 tol = 1e-8, verbose = FALSE)
  
  P <- result$transport_plans[[1]]
  mass <- sum(P)
  converged <- result$converged[1,2]
  
  cat("  Mass:", mass, "\n")
  cat("  Converged:", converged, "\n")
  cat("  Sparsity (zeros):", sum(P < 1e-10), "/", length(P), "\n")
  
  # Check what the gradient looks like at convergence
  C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
  Cx <- as.matrix(dist(X1))
  Cy <- as.matrix(dist(X2))
  p <- rep(1/7, 7)
  q <- rep(1/8, 8)
  
  # Compute final gradient
  Cx2 <- Cx^2
  Cy2 <- Cy^2
  const1_P <- as.vector(Cx2 %*% rowSums(P))
  const2_P <- as.vector(Cy2 %*% colSums(P))
  const_P <- outer(const1_P, rep(1, 8)) + outer(rep(1, 7), const2_P)
  M_P <- const_P - 2 * Cx %*% P %*% t(Cy)
  
  grad <- 0.4 * C_feat + 0.6 * 2 * M_P
  if (lambda > 0) {
    grad <- grad - 4 * lambda * P
  }
  
  cat("  Gradient norm:", norm(grad, "F"), "\n")
  
  results[[as.character(lambda)]] <- list(
    lambda = lambda,
    mass = mass,
    converged = converged,
    P = P
  )
}

# The issue might be that the test's expectation is wrong
# TV penalty doesn't necessarily reduce mass - it encourages sparsity
# Let's check what the actual behavior should be