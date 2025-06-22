library(manifoldalign)
library(multidesign)
library(tibble)

# Test with TV penalty
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

# Compute costs
C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
Cx <- as.matrix(dist(X1))
Cy <- as.matrix(dist(X2))
p <- rep(1/7, 7)
q <- rep(1/8, 8)

cat("Testing Frank-Wolfe with TV penalty lambda=1:\n\n")

# Initialize
P <- outer(p, q)
omega1 <- 0.4
lambda <- 1.0
max_iter <- 10

for (iter in 1:max_iter) {
  # Compute gradient
  Cx2 <- Cx^2
  Cy2 <- Cy^2
  const1_P <- as.vector(Cx2 %*% rowSums(P))
  const2_P <- as.vector(Cy2 %*% colSums(P))
  const_P <- outer(const1_P, rep(1, 8)) + outer(rep(1, 7), const2_P)
  M_P <- const_P - 2 * Cx %*% P %*% t(Cy)
  
  grad <- omega1 * C_feat + (1 - omega1) * 2 * M_P
  
  # Add TV penalty gradient
  grad <- grad - 4 * lambda * P
  
  cat("Iter", iter, "\n")
  cat("  Current mass:", sum(P), "\n")
  cat("  Gradient range:", range(grad), "\n")
  
  # Shift gradient
  G <- grad - min(grad)
  
  # Linear oracle
  gamma_star <- manifoldalign:::partial_ot_tv(G, p, q, lambda)
  cat("  Oracle mass:", sum(gamma_star), "\n")
  
  # Direction
  dP <- gamma_star - P
  
  # Compute step size coefficients
  dr <- rowSums(dP)
  dc <- colSums(dP)
  
  term1 <- sum((Cx2 %*% dr) * dr)
  term2 <- sum((Cy2 %*% dc) * dc)
  term3 <- -2 * sum((Cx %*% dP %*% t(Cy)) * dP)
  
  a <- (1 - omega1) * (term1 + term2 + term3)
  b <- sum(grad * dP)
  
  cat("  a =", a, ", b =", b, "\n")
  
  # Compute optimal step size
  if (abs(a) < 1e-10) {
    alpha <- if (b < 0) 1 else 0
  } else if (a > 0) {
    alpha_star <- -b / (2 * a)
    alpha <- min(1, max(0, alpha_star))
  } else {
    # a < 0
    alpha <- if (a + b < 0) 1 else 0
  }
  
  cat("  Step size:", alpha, "\n")
  
  # Update
  P <- P + alpha * dP
  
  # Check gap
  gap <- abs(b)
  cat("  FW gap:", gap, "\n")
  
  if (gap < 1e-6) {
    cat("  Converged!\n")
    break
  }
  
  cat("\n")
}