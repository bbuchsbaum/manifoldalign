library(manifoldalign)

# Test with different initializations
set.seed(1)
n <- 4
X <- matrix(rnorm(n*2), n)
Y <- matrix(rnorm(n*2), n)
Cx <- as.matrix(dist(X))
Cy <- as.matrix(dist(Y))
C <- as.matrix(dist(rbind(X,Y)))[1:n, n+seq_len(n)]
p <- q <- rep(1/n, n)

# Try different initial transport plans
cat("=== Testing different initializations ===\n\n")

# 1. Uniform initialization (default)
P0_uniform <- outer(p, q)
cat("Uniform initialization - mass:", sum(P0_uniform), "\n")

# 2. Sparse initialization 
P0_sparse <- diag(p)
cat("Sparse initialization - mass:", sum(P0_sparse), "\n")

# 3. Partial initialization
P0_partial <- 0.5 * outer(p, q)
cat("Partial initialization - mass:", sum(P0_partial), "\n")

# Test each initialization
for (init_name in c("uniform", "sparse", "partial")) {
  P0 <- switch(init_name,
               uniform = P0_uniform,
               sparse = P0_sparse,
               partial = P0_partial)
  
  cat("\n--- Testing", init_name, "initialization ---\n")
  
  for (lambda in c(0, 5)) {
    result <- manifoldalign:::fw_fpgw(
      C, Cx, Cy, p, q,
      omega1 = 0.4,
      lambda = lambda,
      rho = NULL,
      P0 = P0,
      max_iter = 40,
      tol = 1e-6,
      verbose = FALSE
    )
    
    mass <- sum(result$P)
    cat(sprintf("  Lambda = %.1f: mass = %.4f, converged = %s\n", 
                lambda, mass, result$converged))
  }
}

# Also check what the gradient looks like
cat("\n=== Gradient structure ===\n")
P <- outer(p, q)
const1_P <- as.vector(Cx^2 %*% rowSums(P))
const2_P <- as.vector(Cy^2 %*% colSums(P))
const_P <- outer(const1_P, rep(1, length(q))) + outer(rep(1, length(p)), const2_P)
M_P <- const_P - 2 * Cx %*% P %*% t(Cy)
grad <- 0.4 * C + 0.6 * 2 * M_P

cat("Gradient range:", range(grad), "\n")
cat("Gradient matrix:\n")
print(round(grad, 3))

# Check if gradient structure prevents partial transport
G <- grad - min(grad)
gamma_star <- manifoldalign:::classical_ot_lp(G, p, q)
cat("\nOptimal transport for shifted gradient:\n")
print(round(gamma_star, 3))
cat("Mass:", sum(gamma_star), "\n")