library(manifoldalign)

# Debug augmented solution
n <- 3
m <- 3
p <- rep(1/3, 3)
q <- rep(1/3, 3)

# Simple cost matrix
C <- matrix(5, n, m)

cat("=== Debugging Augmented Solution ===\n")
cat("Original cost matrix:\n")
print(C)

for (lambda in c(0, 5, 20)) {
  cat("\n--- Lambda =", lambda, "---\n")
  
  # Create augmented cost matrix manually
  C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
  C_aug[1:n, 1:m] <- C + lambda
  
  cat("Augmented cost matrix:\n")
  print(C_aug)
  
  # Augmented marginals
  p_aug <- c(p, sum(q))
  q_aug <- c(q, sum(p))
  
  cat("\nAugmented marginals:\n")
  cat("p_aug:", p_aug, "\n")
  cat("q_aug:", q_aug, "\n")
  
  # Solve
  gamma_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)
  
  cat("\nAugmented solution:\n")
  print(round(gamma_aug, 4))
  
  cat("\nMass distribution:\n")
  cat("  Top-left (actual transport):", sum(gamma_aug[1:n, 1:m]), "\n")
  cat("  To dummy sink (row n+1):", sum(gamma_aug[n+1, 1:m]), "\n")
  cat("  From dummy source (col m+1):", sum(gamma_aug[1:n, m+1]), "\n")
  cat("  Dummy-to-dummy:", gamma_aug[n+1, m+1], "\n")
  
  # Check if the issue is with uniform costs
  cat("\nAnalysis: When C+lambda =", C[1,1] + lambda, "vs escape cost = 0\n")
  cat("The solver should prefer escape when C+lambda > 0\n")
}

# Test with gradient-like matrix (some negative values)
cat("\n\n=== Testing with gradient-like cost matrix ===\n")
G <- matrix(c(-2, 3, 1, 
              4, -1, 2,
              0, 2, -3), nrow=3, byrow=TRUE)
cat("Gradient matrix G:\n")
print(G)

for (lambda in c(0, 1, 5)) {
  cat("\n--- Lambda =", lambda, "---\n")
  
  if (lambda == 0) {
    # Shift to make non-negative for classical OT
    G_shifted <- G - min(G)
    gamma <- manifoldalign:::classical_ot_lp(G_shifted, p, q)
  } else {
    gamma <- manifoldalign:::partial_ot_tv(G, p, q, lambda)
  }
  
  cat("Transport plan:\n")
  print(round(gamma, 4))
  cat("Total transported mass:", sum(gamma), "\n")
}