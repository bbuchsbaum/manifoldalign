library(manifoldalign)

# Test the network simplex solver
set.seed(123)

# Test 1: Small even-sized problem
cat("Test 1: Even-sized problem (10x10)\n")
n <- 10
m <- 10
cost <- matrix(runif(n * m), n, m)
a <- rep(1/n, n)
b <- rep(1/m, m)

result1 <- manifoldalign:::network_simplex_ot_cpp(cost, a, b)
cat("Solution dimensions:", dim(result1), "\n")
cat("Total transported mass:", sum(result1), "\n")
cat("Row sums match a:", all(abs(rowSums(result1) - a) < 1e-9), "\n")
cat("Col sums match b:", all(abs(colSums(result1) - b) < 1e-9), "\n\n")

# Test 2: Odd-sized problem (should trigger padding)
cat("Test 2: Odd-sized problem (11x9)\n")
n <- 11
m <- 9
cost <- matrix(runif(n * m), n, m)
a <- rep(1/n, n)
b <- rep(1/m, m)

result2 <- manifoldalign:::network_simplex_ot_cpp(cost, a, b)
cat("Solution dimensions:", dim(result2), "\n")
cat("Total transported mass:", sum(result2), "\n")
cat("Row sums match a:", all(abs(rowSums(result2) - a) < 1e-9), "\n")
cat("Col sums match b:", all(abs(colSums(result2) - b) < 1e-9), "\n\n")

# Test 3: Partial OT with mass constraint
cat("Test 3: Partial OT with mass constraint\n")
n <- 10
m <- 10
cost <- matrix(runif(n * m), n, m)
a <- rep(0.2, n)  # Total mass = 2
b <- rep(0.15, m) # Total mass = 1.5
mass <- 1.0  # Transport only 1 unit

result3 <- manifoldalign:::partial_ot_mass_rcpp(cost, a, b, mass)
cat("Solution dimensions:", dim(result3), "\n")
cat("Total transported mass:", sum(result3), "\n")
cat("Mass equals constraint:", abs(sum(result3) - mass) < 1e-9, "\n")
cat("Row sums <= a:", all(rowSums(result3) <= a + 1e-9), "\n")
cat("Col sums <= b:", all(colSums(result3) <= b + 1e-9), "\n\n")

# Test 4: TV proximal operator
cat("Test 4: TV proximal operator\n")
Y <- matrix(rnorm(10 * 10), 10, 10)
lambda <- 0.5
Y_tv <- manifoldalign:::apply_tv_proximal_cpp(Y, lambda)
cat("Input dimensions:", dim(Y), "\n")
cat("Output dimensions:", dim(Y_tv), "\n")
cat("TV reduction achieved:", sum(abs(Y)) - sum(abs(Y_tv)), "\n")