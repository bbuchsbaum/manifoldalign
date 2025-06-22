library(manifoldalign)

# Force using the fallback solver by testing functions that don't trigger LEMON
set.seed(123)

# Test 1: Use lpSolve directly  
cat("Test 1: Classical OT with lpSolve\n")
n <- 5
m <- 5
cost <- matrix(runif(n * m), n, m)
p <- rep(1/n, n)
q <- rep(1/m, m)

# This should fall back to lpSolve when network_simplex_ot_cpp fails
tryCatch({
  result1 <- manifoldalign:::classical_ot_lp(cost, p, q)
  cat("Solution found with sum:", sum(result1), "\n")
}, error = function(e) {
  cat("Error:", e$message, "\n")
})

# Test 2: Test the debug function
cat("\nTest 2: Debug network simplex (northwest corner)\n")
result2 <- manifoldalign:::debug_network_simplex(cost, p, q)
cat("Northwest corner sum:", sum(result2), "\n")

# Test 3: Test Sinkhorn
cat("\nTest 3: Sinkhorn solver\n")
epsilon <- 0.1
result3 <- manifoldalign:::sinkhorn_unified(cost, p, q, epsilon, max_iter = 100)
cat("Sinkhorn solution sum:", sum(result3), "\n")
cat("Row sums match p:", all(abs(rowSums(result3) - p) < 1e-6), "\n")
cat("Col sums match q:", all(abs(colSums(result3) - q) < 1e-6), "\n")

# Test 4: TV proximal
cat("\nTest 4: TV proximal operator\n")
Y <- matrix(rnorm(25), 5, 5)
Y_tv <- manifoldalign:::apply_tv_proximal_cpp(Y, lambda = 0.1)
cat("TV applied successfully, dimensions:", dim(Y_tv), "\n")