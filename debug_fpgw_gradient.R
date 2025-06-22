library(manifoldalign)
library(multidesign)

# Reproduce FPGW test case to examine gradients
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

cat("Domain shapes: X1 =", dim(X1), ", X2 =", dim(X2), "\n\n")

# Compute feature cost
C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
cat("Feature cost matrix shape:", dim(C_feat), "\n")
cat("Feature cost range:", range(C_feat), "\n\n")

# Create a simple Frank-Wolfe iteration to see the gradient
p <- rep(1/nrow(X1), nrow(X1))
q <- rep(1/nrow(X2), nrow(X2))

# Initial transport plan
P0 <- outer(p, q)

# Compute within-domain distances
Cx <- as.matrix(dist(X1))
Cy <- as.matrix(dist(X2))
Cx <- Cx / (max(Cx) + 1e-10)
Cy <- Cy / (max(Cy) + 1e-10)

# Set parameters
omega1 <- 0.4
feature_weight <- omega1
struct_weight <- 1 - omega1

# Compute gradient at P0
Cx2 <- Cx^2
Cy2 <- Cy^2
const1_P <- as.vector(Cx2 %*% rowSums(P0))
const2_P <- as.vector(Cy2 %*% colSums(P0))
const_P <- outer(const1_P, rep(1, ncol(P0))) + outer(rep(1, nrow(P0)), const2_P)
M_P <- const_P - 2 * Cx %*% P0 %*% t(Cy)

grad <- feature_weight * C_feat + struct_weight * 2 * M_P

cat("Gradient at initial P0:\n")
cat("  Shape:", dim(grad), "\n")
cat("  Range:", range(grad), "\n")
cat("  Number of negative entries:", sum(grad < 0), "\n")
cat("  Min value:", min(grad), "\n\n")

# Test what happens with TV penalty
for (lambda in c(0, 1, 5)) {
  cat("--- Lambda =", lambda, "---\n")
  
  # Shift gradient for classical OT if lambda = 0
  if (lambda == 0) {
    G <- grad - min(grad)
    gamma <- manifoldalign:::classical_ot_lp(G, p, q)
  } else {
    # For TV penalty, pass gradient directly
    gamma <- manifoldalign:::partial_ot_tv(grad, p, q, lambda)
  }
  
  cat("  Transported mass:", sum(gamma), "\n")
  
  # Check the augmented cost matrix values
  G_plus_lambda <- grad + lambda
  cat("  Range of (grad + lambda):", range(G_plus_lambda), "\n")
  cat("  Negative entries in (grad + lambda):", sum(G_plus_lambda < 0), "\n")
  
  # For lambda > 0, check if any transport costs are negative
  if (lambda > 0 && any(G_plus_lambda < 0)) {
    cat("  WARNING: Some transport costs are still negative!\n")
    cat("  This encourages transport instead of escape.\n")
  }
  cat("\n")
}

cat("INSIGHT: The issue occurs when grad + lambda still has negative values.\n")
cat("In this case, those negative-cost routes are cheaper than the zero-cost escape!\n")