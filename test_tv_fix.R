library(manifoldalign)

cat("=== Testing TV Penalty Fix ===\n\n")

# Test 1: Simple uniform cost matrix
cat("Test 1: Uniform cost matrix with cheap diagonal\n")
n <- 3
m <- 3
p <- rep(1/n, n)
q <- rep(1/m, m)

# Create cost matrix with cheap diagonal
C <- matrix(10, n, m)
diag(C) <- 1

cat("Cost matrix:\n")
print(C)
cat("\n")

# Test with different lambda values
cat("Lambda effect on transported mass:\n")
for (lambda in c(0, 0.5, 1, 2, 5, 10, 20)) {
  if (lambda == 0) {
    gamma <- manifoldalign:::classical_ot_lp(C, p, q)
  } else {
    gamma <- manifoldalign:::partial_ot_tv(C, p, q, lambda)
  }
  mass <- sum(gamma)
  cat(sprintf("  Lambda = %4.1f: transported mass = %.4f\n", lambda, mass))
}

# Test 2: FPGW test case with outliers
cat("\n\nTest 2: FPGW with outliers (reproducing test case)\n")
library(multidesign)
library(tibble)

set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

# Create hyperdesign
design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Test with different lambda values
lambdas <- c(0.0, 1.0, 5.0)
masses <- numeric(length(lambdas))

for (i in seq_along(lambdas)) {
  lambda <- lambdas[i]
  cat(sprintf("\nTesting lambda = %.1f...\n", lambda))
  
  result <- fpgw(hd, omega1 = 0.4, lambda = lambda, max_iter = 40, verbose = FALSE)
  masses[i] <- sum(result$transport_plans[[1]])
  
  cat(sprintf("  Transported mass: %.4f\n", masses[i]))
  cat(sprintf("  Converged: %s\n", result$converged[1,2]))
}

# Check monotonic decrease
cat("\n=== Verification ===\n")
cat(sprintf("mass(λ=0) = %.4f\n", masses[1]))
cat(sprintf("mass(λ=1) = %.4f\n", masses[2]))
cat(sprintf("mass(λ=5) = %.4f\n", masses[3]))
cat(sprintf("\nmass(λ=0) > mass(λ=1)? %s\n", masses[1] > masses[2]))
cat(sprintf("mass(λ=1) > mass(λ=5)? %s\n", masses[2] > masses[3]))

if (masses[1] > masses[2] && masses[2] > masses[3]) {
  cat("\n✓ SUCCESS: TV penalty correctly reduces transported mass!\n")
} else {
  cat("\n✗ FAILED: TV penalty still not working correctly\n")
}

# Test 3: Extreme lambda values
cat("\n\nTest 3: Extreme lambda values\n")
C_simple <- matrix(5, 3, 3)
p_simple <- rep(1/3, 3)
q_simple <- rep(1/3, 3)

for (lambda in c(0, 100, 1000)) {
  if (lambda == 0) {
    gamma <- manifoldalign:::classical_ot_lp(C_simple, p_simple, q_simple)
  } else {
    gamma <- manifoldalign:::partial_ot_tv(C_simple, p_simple, q_simple, lambda)
  }
  cat(sprintf("  Lambda = %4.0f: transported mass = %.6f\n", lambda, sum(gamma)))
}