library(manifoldalign)

# Test case from failing test
set.seed(123)
G <- matrix(runif(9), 3, 3)
p <- c(0.3, 0.4, 0.3)
q <- c(0.2, 0.5, 0.3)
rho <- 0.5

cat("Cost matrix G:\n")
print(G)
cat("\nMarginals: p =", p, ", q =", q, "\n")
cat("Mass constraint rho =", rho, "\n")
cat("Sum of p =", sum(p), ", Sum of q =", sum(q), "\n\n")

# Test C++ implementation
cat("Testing C++ implementation:\n")
gamma_cpp <- manifoldalign:::partial_ot_mass_rcpp(G, p, q, rho)
cat("Transport plan:\n")
print(gamma_cpp)
cat("Row sums:", rowSums(gamma_cpp), "\n")
cat("Col sums:", colSums(gamma_cpp), "\n")
cat("Total mass:", sum(gamma_cpp), "\n")
cat("Row sums <= p?", all(rowSums(gamma_cpp) <= p + 1e-6), "\n")
cat("Col sums <= q?", all(colSums(gamma_cpp) <= q + 1e-6), "\n\n")

# Test R implementation  
cat("Testing R implementation:\n")
gamma_r <- manifoldalign:::partial_ot_mass(G, p, q, rho, method = "lpSolve")
cat("Transport plan:\n")
print(gamma_r)
cat("Row sums:", rowSums(gamma_r), "\n")
cat("Col sums:", colSums(gamma_r), "\n")
cat("Total mass:", sum(gamma_r), "\n")
cat("Row sums <= p?", all(rowSums(gamma_r) <= p + 1e-6), "\n")
cat("Col sums <= q?", all(colSums(gamma_r) <= q + 1e-6), "\n")