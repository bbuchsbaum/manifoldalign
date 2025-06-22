library(manifoldalign)

# Temporarily disable LEMON to test fallback
# Let's test the classical OT LP first
set.seed(123)

# Simple test case
n <- 5
m <- 5
cost <- matrix(runif(n * m), n, m)
p <- rep(1/n, n)
q <- rep(1/m, m)

cat("Testing classical OT LP...\n")
result <- manifoldalign:::classical_ot_lp(cost, p, q)
cat("Solution dimensions:", dim(result), "\n")
cat("Total transported mass:", sum(result), "\n")
cat("Row sums match p:", all(abs(rowSums(result) - p) < 1e-6), "\n")
cat("Col sums match q:", all(abs(colSums(result) - q) < 1e-6), "\n\n")

# Test partial OT
cat("Testing partial OT mass...\n")
p2 <- rep(0.3, n)
q2 <- rep(0.2, m)
mass <- 0.5
result2 <- manifoldalign:::partial_ot_mass(cost, p2, q2, mass)
cat("Solution dimensions:", dim(result2), "\n")
cat("Total transported mass:", sum(result2), "\n")
cat("Mass constraint satisfied:", abs(sum(result2) - mass) < 1e-6, "\n")