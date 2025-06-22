library(manifoldalign)

# Simple test case
set.seed(123)
n <- 3
m <- 3
p <- rep(1/n, n)
q <- rep(1/m, m)

# Create a simple cost matrix where diagonal is cheap
C <- matrix(10, n, m)
diag(C) <- 1

cat("Cost matrix:\n")
print(C)

# Test with different lambdas
for (lambda in c(0, 1, 5)) {
  cat("\n=== Testing lambda =", lambda, "===\n")
  
  # Call partial_ot_tv directly
  if (lambda > 0) {
    gamma <- manifoldalign:::partial_ot_tv(C, p, q, lambda)
  } else {
    gamma <- manifoldalign:::classical_ot_lp(C, p, q)
  }
  
  cat("Transport plan:\n")
  print(round(gamma, 4))
  cat("Total mass:", sum(gamma), "\n")
  cat("Row sums:", rowSums(gamma), "\n")
  cat("Col sums:", colSums(gamma), "\n")
}

# Now test what happens in the augmented problem
cat("\n=== Augmented problem details ===\n")
C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
C_aug[1:n, 1:m] <- C
p_aug <- c(p, sum(q))
q_aug <- c(q, sum(p))

cat("Augmented cost matrix:\n")
print(C_aug)
cat("\nAugmented marginals:\n")
cat("p_aug:", p_aug, "\n")
cat("q_aug:", q_aug, "\n")

# Solve augmented problem
gamma_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)
cat("\nAugmented solution:\n")
print(round(gamma_aug, 4))
cat("Top-left block mass:", sum(gamma_aug[1:n, 1:m]), "\n")
cat("Dummy source row mass:", sum(gamma_aug[n+1, ]), "\n")
cat("Dummy sink col mass:", sum(gamma_aug[, m+1]), "\n")