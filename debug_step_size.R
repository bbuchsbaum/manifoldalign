library(manifoldalign)
library(multidesign)
library(tibble)

# Debug the step size computation
set.seed(123)
n <- 10  # Smaller for debugging

X1 <- matrix(rnorm(n * 2), n, 2)
X2 <- matrix(rnorm(n * 2), n, 2)

design <- data.frame(id = 1:n)
hd <- hyperdesign(list(
  d1 = multidesign(X1, design),
  d2 = multidesign(X2, design)
))

# Extract data
prep <- manifoldalign:::prepare_ot_data(hd)
X_list <- prep$X_list

# Compute cost matrices
C_feat <- manifoldalign:::compute_feature_cost(X_list[[1]], X_list[[2]], "euclidean")
Cx <- as.matrix(dist(X_list[[1]]))
Cy <- as.matrix(dist(X_list[[2]]))

# Normalize
Cx <- Cx / (max(Cx) + 1e-10)
Cy <- Cy / (max(Cy) + 1e-10)

# Marginals
p <- rep(1/n, n)
q <- rep(1/n, n)

# Initialize
P <- matrix(1/(n*n), n, n)
omega1 <- 0.5
struct_weight <- 1 - omega1

# Pre-compute
Cx2 <- Cx^2
Cy2 <- Cy^2
const1 <- as.vector(Cx2 %*% p)
const2 <- as.vector(Cy2 %*% q)
const <- outer(const1, rep(1, n)) + outer(rep(1, n), const2)

# Compute gradient
M_P <- const - 2 * Cx %*% P %*% t(Cy)
grad <- omega1 * C_feat + struct_weight * M_P

# Get linear oracle solution (simplified - uniform for debugging)
gamma_star <- matrix(1/n, n, n)  # Vertex of transport polytope

# Direction
dP <- gamma_star - P

# Step size computation (OLD WAY)
a_old <- struct_weight * sum((const - 2 * Cx %*% dP %*% t(Cy)) * dP)
b <- sum(grad * dP)

cat("OLD step size computation:\n")
cat("  a =", a_old, "\n")
cat("  b =", b, "\n")
cat("  alpha = ", if(a_old <= 0) "edge case" else -b/(2*a_old), "\n\n")

# Step size computation (NEW WAY)
A <- Cx %*% dP %*% t(Cy)
a_new <- 2 * struct_weight * sum(A * dP)

cat("NEW step size computation:\n")
cat("  a =", a_new, "\n")
cat("  b =", b, "\n")
cat("  alpha = ", if(abs(a_new) < 1e-10) "linear case" else -b/(2*a_new), "\n\n")

# Check the term that should give positive a
cat("Debugging a computation:\n")
cat("  sum(A * dP) =", sum(A * dP), "\n")
cat("  This should be >= 0 (it's <L, dP ⊗ dP>)\n")

# Verify gradient is correct
cat("\nGradient check:\n")
cat("  ||grad||_F =", norm(grad, "F"), "\n")
cat("  <grad, dP> =", sum(grad * dP), "\n")