library(manifoldalign)

# Simple test of augmentation
n <- 3
m <- 4
p <- rep(1/n, n)
q <- rep(1/m, m)

# Create a cost matrix where transport is expensive
cost <- matrix(10, n, m)  # High cost

cat("Testing augmentation with high cost:\n")
cat("Original cost matrix (all entries = 10):\n")

# Manual augmentation
C_aug <- matrix(0, n+1, m+1)
C_aug[1:n, 1:m] <- cost
p_aug <- c(p, sum(q))
q_aug <- c(q, sum(p))

cat("\nAugmented cost matrix:\n")
print(C_aug)
cat("\np_aug:", p_aug, "\n")
cat("q_aug:", q_aug, "\n")

# Solve
gamma_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)

cat("\nAugmented solution:\n")
print(round(gamma_aug, 4))

cat("\nMass in original variables:", sum(gamma_aug[1:n, 1:m]), "\n")
cat("Mass to dummy sink (from real sources):", sum(gamma_aug[1:n, m+1]), "\n")
cat("Mass from dummy source (to real sinks):", sum(gamma_aug[n+1, 1:m]), "\n")

# Now test with mixed costs
cat("\n\nTesting with mixed costs:\n")
cost2 <- matrix(c(1, 10, 5, 8,
                  10, 2, 10, 6,
                  5, 10, 3, 10), nrow=3, byrow=TRUE)
                  
C_aug2 <- matrix(0, n+1, m+1)
C_aug2[1:n, 1:m] <- cost2

gamma_aug2 <- manifoldalign:::classical_ot_lp(C_aug2, p_aug, q_aug)

cat("Mixed cost solution:\n")
print(round(gamma_aug2, 4))
cat("\nMass in original variables:", sum(gamma_aug2[1:n, 1:m]), "\n")