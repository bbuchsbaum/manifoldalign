library(manifoldalign)
library(multidesign)
library(tibble)

# Test setup
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Get cost matrices
C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
Cx <- as.matrix(dist(X1))
Cy <- as.matrix(dist(X2))
p <- rep(1/7, 7)
q <- rep(1/8, 8)

# Initialize
P0 <- matrix(1/(7*8), 7, 8)

cat("Testing Frank-Wolfe with lambda = 5.0:\n")

# Manual Frank-Wolfe iteration
P <- P0
lambda <- 5.0
omega1 <- 0.4

for (iter in 1:5) {
  # Compute gradient (without TV term)
  Cx2 <- Cx^2
  Cy2 <- Cy^2
  const1_P <- as.vector(Cx2 %*% rowSums(P))
  const2_P <- as.vector(Cy2 %*% colSums(P))
  const_P <- outer(const1_P, rep(1, 8)) + outer(rep(1, 7), const2_P)
  M_P <- const_P - 2 * Cx %*% P %*% t(Cy)
  
  grad <- omega1 * C_feat + (1 - omega1) * 2 * M_P
  
  cat("\nIteration", iter, "\n")
  cat("Current mass:", sum(P), "\n")
  cat("Gradient range:", range(grad), "\n")
  
  # Add lambda and shift
  G <- grad + lambda
  G <- G - min(G)
  cat("After adding lambda and shifting, G range:", range(G), "\n")
  
  # Call oracle
  gamma_star <- manifoldalign:::partial_ot_tv(G, p, q, lambda)
  cat("Oracle returned mass:", sum(gamma_star), "\n")
  
  # If oracle returns 0 mass, check why
  if (sum(gamma_star) < 0.01) {
    cat("Oracle returned near-zero mass. Checking augmented solution...\n")
    # Check what the augmented problem looks like
    n <- length(p)
    m <- length(q)
    C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
    C_aug[1:n, 1:m] <- G
    p_aug <- c(p, sum(q))
    q_aug <- c(q, sum(p))
    
    cat("Min cost in G:", min(G), "\n")
    cat("Max cost in G:", max(G), "\n")
    cat("Cost to dummy (0) vs real transport (", min(G), ")\n")
  }
  
  # Direction
  dP <- gamma_star - P
  
  if (max(abs(dP)) < 1e-10) {
    cat("No movement possible, stopping.\n")
    break
  }
  
  # Simple step
  alpha <- 0.1
  P <- P + alpha * dP
}