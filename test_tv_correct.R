library(manifoldalign)
library(multidesign)
library(tibble)

cat("=== Testing Corrected TV Penalty Implementation ===\n\n")

# Test 1: Simple synthetic example from the memo
cat("Test 1: Minimal reproducible example\n")
set.seed(1)
n <- 4
X <- matrix(rnorm(n*2), n)
Y <- matrix(rnorm(n*2), n)
Cx <- as.matrix(dist(X))
Cy <- as.matrix(dist(Y))
C <- as.matrix(dist(rbind(X,Y)))[1:n, n+seq_len(n)]
p <- q <- rep(1/n, n)

# Test with different lambda values
lambdas <- c(0, 1, 5, 10)
masses <- numeric(length(lambdas))

for (i in seq_along(lambdas)) {
  lambda <- lambdas[i]
  cat(sprintf("\nLambda = %.1f:\n", lambda))
  
  result <- manifoldalign:::fw_fpgw(
    C, Cx, Cy, p, q,
    omega1 = 0.4,
    lambda = lambda,
    rho = NULL,
    P0 = outer(p, q),
    max_iter = 40,
    tol = 1e-6,
    verbose = FALSE
  )
  
  masses[i] <- sum(result$P)
  cat(sprintf("  Transported mass: %.4f\n", masses[i]))
  cat(sprintf("  Converged: %s\n", result$converged))
  
  # Check objective trace
  obj_trace <- attr(result, "obj_trace")
  if (length(obj_trace) > 1) {
    cat(sprintf("  Objective decreased from %.4f to %.4f\n", 
                obj_trace[1], obj_trace[length(obj_trace)]))
  }
}

# Verify monotonic decrease
cat("\n=== Mass Transport Summary ===\n")
for (i in seq_along(lambdas)) {
  cat(sprintf("Lambda = %.1f: mass = %.4f\n", lambdas[i], masses[i]))
}

cat("\nMonotonic decrease check:\n")
for (i in 2:length(masses)) {
  decrease <- masses[i-1] > masses[i]
  cat(sprintf("  mass(λ=%.1f) > mass(λ=%.1f)? %s\n", 
              lambdas[i-1], lambdas[i], decrease))
}

# Test 2: FPGW with outliers
cat("\n\nTest 2: FPGW with outliers (original test case)\n")
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Test with different lambda values
test_lambdas <- c(0.0, 1.0, 5.0)
test_masses <- numeric(length(test_lambdas))

for (i in seq_along(test_lambdas)) {
  lambda <- test_lambdas[i]
  cat(sprintf("\nTesting lambda = %.1f...\n", lambda))
  
  result <- fpgw(hd, omega1 = 0.4, lambda = lambda, max_iter = 40, verbose = FALSE)
  test_masses[i] <- sum(result$transport_plans[[1]])
  
  cat(sprintf("  Transported mass: %.4f\n", test_masses[i]))
  cat(sprintf("  Converged: %s\n", result$converged[1,2]))
}

# Check monotonic decrease
cat("\n=== FPGW Test Results ===\n")
cat(sprintf("mass(λ=0) = %.4f\n", test_masses[1]))
cat(sprintf("mass(λ=1) = %.4f\n", test_masses[2]))
cat(sprintf("mass(λ=5) = %.4f\n", test_masses[3]))
cat(sprintf("\nmass(λ=0) > mass(λ=1)? %s\n", test_masses[1] > test_masses[2]))
cat(sprintf("mass(λ=1) > mass(λ=5)? %s\n", test_masses[2] > test_masses[3]))

if (all(diff(test_masses) < 0)) {
  cat("\n✓ SUCCESS: TV penalty correctly reduces transported mass!\n")
} else {
  cat("\n✗ Still not working correctly\n")
}