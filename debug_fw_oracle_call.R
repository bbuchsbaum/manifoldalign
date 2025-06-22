library(manifoldalign)

# Patch the partial_ot_tv function to add debug output
partial_ot_tv_debug <- function(cost, p, q, lambda) {
  cat("\n=== partial_ot_tv called ===\n")
  cat("Lambda:", lambda, "\n")
  cat("Cost matrix shape:", dim(cost), "\n")
  cat("Cost range:", range(cost), "\n")
  
  n <- length(p)
  m <- length(q)
  
  # Create augmented cost matrix
  C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
  C_aug[1:n, 1:m] <- cost + lambda
  
  cat("Augmented cost range:", range(C_aug), "\n")
  cat("Top-left block range:", range(C_aug[1:n, 1:m]), "\n")
  
  # Augmented marginals
  p_aug <- c(p, sum(q))
  q_aug <- c(q, sum(p))
  
  # Solve
  gamma_aug <- manifoldalign:::classical_ot_lp_r(C_aug, p_aug, q_aug)
  
  # Extract result
  gamma <- gamma_aug[1:n, 1:m]
  
  cat("Result: transported mass =", sum(gamma), "\n")
  cat("Escape to dummy:", sum(gamma_aug[1:n, m+1]), "\n")
  cat("From dummy:", sum(gamma_aug[n+1, 1:m]), "\n")
  
  return(gamma)
}

# Replace the function temporarily
assignInNamespace("partial_ot_tv", partial_ot_tv_debug, "manifoldalign")

# Run a simple test
library(multidesign)
library(tibble)
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)
)

design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Test with lambda = 5
cat("Running FPGW with lambda = 5\n")
result <- fpgw(hd, omega1 = 0.4, lambda = 5.0, max_iter = 5, verbose = TRUE)