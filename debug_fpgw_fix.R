library(manifoldalign)
library(multidesign)
library(tibble)

# Test the step size fix
set.seed(123)
n <- 20

# Create simple clustered data
X1 <- matrix(rnorm(n * 3), n, 3)
X1[1:10, ] <- X1[1:10, ] + 3  # Stronger separation

X2 <- matrix(rnorm(n * 5), n, 5) 
X2[1:10, ] <- X2[1:10, ] + 3

design <- data.frame(id = 1:n)
hd <- hyperdesign(list(
  d1 = multidesign(X1, design),
  d2 = multidesign(X2, design)
))

# Test with different omega values
cat("Testing FPGW with step size fix:\n")
cat("================================\n\n")

# More structural emphasis (should converge better)
cat("1. Structural emphasis (omega1 = 0.01):\n")
result1 <- fpgw(hd, omega1 = 0.01, verbose = TRUE, max_iter = 30)
cat("\nConverged:", result1$converged[1,2], "\n")
cat("Final distance:", result1$distances[1,2], "\n\n")

# Balanced
cat("2. Balanced (omega1 = 0.5):\n")
result2 <- fpgw(hd, omega1 = 0.5, verbose = TRUE, max_iter = 30)
cat("\nConverged:", result2$converged[1,2], "\n") 
cat("Final distance:", result2$distances[1,2], "\n\n")

# Feature emphasis
cat("3. Feature emphasis (omega1 = 0.99):\n")
result3 <- fpgw(hd, omega1 = 0.99, verbose = TRUE, max_iter = 30)
cat("\nConverged:", result3$converged[1,2], "\n")
cat("Final distance:", result3$distances[1,2], "\n\n")

# Test mass-constrained with exact equality
cat("4. Mass-constrained (rho = 0.8):\n")
result4 <- fpgw(hd, omega1 = 0.3, rho = 0.8, verbose = TRUE, max_iter = 30)
P <- result4$transport_plans[[1]]
cat("\nConverged:", result4$converged[1,2], "\n")
cat("Transported mass:", sum(P), "(should be exactly 0.8)\n")
cat("Final distance:", result4$distances[1,2], "\n")