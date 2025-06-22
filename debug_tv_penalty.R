library(manifoldalign)
library(multidesign)
library(tibble)

# Test with TV penalty
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

# Create appropriate design frames for each domain
design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))

# Create multidesign objects manually since domains have different sizes
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

cat("Testing TV penalty effect on transported mass:\n")
cat("Domain 1 has", nrow(X1), "samples\n")
cat("Domain 2 has", nrow(X2), "samples\n\n")

# Run with different lambda values
for (lambda in c(0.0, 1.0, 5.0)) {
  cat("\nLambda =", lambda, "\n")
  result <- fpgw(hd, omega1 = 0.4, lambda = lambda, max_iter = 40, verbose = FALSE)
  P <- result$transport_plans[[1]]
  mass <- sum(P)
  
  cat("Transport plan shape:", nrow(P), "x", ncol(P), "\n")
  cat("Total transported mass:", mass, "\n")
  cat("Converged:", result$converged[1,2], "\n")
  
  # Check marginals
  row_sums <- rowSums(P)
  col_sums <- colSums(P)
  cat("Row sums range:", range(row_sums), "\n")
  cat("Col sums range:", range(col_sums), "\n")
  
  # Number of sparse entries
  cat("Entries < 1e-10:", sum(P < 1e-10), "out of", length(P), "\n")
}

# Test the partial_ot_tv function directly
cat("\n\nTesting partial_ot_tv directly:\n")
cost <- matrix(runif(7*8), 7, 8)
p <- rep(1/7, 7)
q <- rep(1/8, 8)

for (lambda in c(0.0, 1.0, 5.0)) {
  P_tv <- manifoldalign:::partial_ot_tv(cost, p, q, lambda)
  cat("\nLambda =", lambda, "-> Mass =", sum(P_tv), "\n")
}