library(manifoldalign)

# Test the partial_ot_tv oracle directly
set.seed(123)
n <- 4
m <- 5
cost <- matrix(runif(n * m), n, m)
p <- rep(1/n, n)
q <- rep(1/m, m)

cat("Testing partial_ot_tv oracle:\n\n")

# Without augmentation (classical OT)
P_classical <- manifoldalign:::classical_ot_lp(cost, p, q)
cat("Classical OT mass:", sum(P_classical), "\n")

# With augmentation (partial OT)
P_partial <- manifoldalign:::partial_ot_tv(cost, p, q, lambda = 1)
cat("Partial OT mass:", sum(P_partial), "\n")

# Try with different costs
cat("\nWith higher cost:\n")
high_cost <- cost + 10
P_high <- manifoldalign:::partial_ot_tv(high_cost, p, q, lambda = 1)
cat("Partial OT mass (high cost):", sum(P_high), "\n")

# Check augmented problem directly
cat("\nAugmented problem details:\n")
C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
C_aug[1:n, 1:m] <- cost
p_aug <- c(p, sum(q))
q_aug <- c(q, sum(p))

cat("Augmented cost matrix shape:", dim(C_aug), "\n")
cat("p_aug:", p_aug, "\n")
cat("q_aug:", q_aug, "\n")
cat("sum(p_aug):", sum(p_aug), "\n")
cat("sum(q_aug):", sum(q_aug), "\n")

# Solve augmented
P_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)
cat("\nAugmented solution:\n")
print(round(P_aug, 4))
cat("\nTotal mass in top-left block:", sum(P_aug[1:n, 1:m]), "\n")