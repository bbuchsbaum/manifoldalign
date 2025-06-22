library(manifoldalign)

# Test augmentation manually
n <- 3
m <- 3
p <- rep(1/n, n)
q <- rep(1/m, m)
cost <- matrix(100, n, m)

cat("Manual augmentation test:\n")
cat("Original cost (all 100):\n")

# Create augmented problem
C_aug <- matrix(0, n+1, m+1)
C_aug[1:n, 1:m] <- cost
p_aug <- c(p, sum(q))
q_aug <- c(q, sum(p))

cat("\nAugmented cost matrix:\n")
print(C_aug)

# Check if classical_ot_lp works
cat("\nTesting classical_ot_lp on original problem:\n")
gamma_orig <- manifoldalign:::classical_ot_lp(cost, p, q)
cat("Mass =", sum(gamma_orig), "\n")

cat("\nTesting classical_ot_lp on augmented problem:\n")
gamma_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)
print(round(gamma_aug, 4))
cat("Mass in top-left:", sum(gamma_aug[1:n, 1:m]), "\n")

# Now test our partial_ot_tv function
cat("\nTesting partial_ot_tv:\n")
gamma_partial <- manifoldalign:::partial_ot_tv(cost, p, q, lambda = 1)
cat("Mass from partial_ot_tv:", sum(gamma_partial), "\n")

# Check if network_simplex is being used
cat("\nChecking for network_simplex_ot_cpp:\n")
exists_ns <- exists("network_simplex_ot_cpp", where = asNamespace("manifoldalign"))
cat("network_simplex_ot_cpp exists:", exists_ns, "\n")