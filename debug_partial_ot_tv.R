library(manifoldalign)

# Add debug output to understand what's happening
debug_partial_ot_tv <- function(cost, p, q, lambda) {
  n <- length(p)
  m <- length(q)
  
  cat("In partial_ot_tv:\n")
  cat("  Input cost range:", range(cost), "\n")
  
  # Create augmented cost matrix (n+1) x (m+1)
  C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
  C_aug[1:n, 1:m] <- cost
  
  cat("  Augmented cost shape:", dim(C_aug), "\n")
  cat("  C_aug[1,1]:", C_aug[1,1], ", C_aug[n+1,m+1]:", C_aug[n+1,m+1], "\n")
  
  # Augmented marginals
  p_aug <- c(p, sum(q))
  q_aug <- c(q, sum(p))
  
  # Solve using classical_ot_lp
  cat("  Calling classical_ot_lp...\n")
  gamma_aug <- manifoldalign:::classical_ot_lp(C_aug, p_aug, q_aug)
  
  cat("  Solution shape:", dim(gamma_aug), "\n")
  cat("  Total mass:", sum(gamma_aug), "\n")
  cat("  Mass in top-left:", sum(gamma_aug[1:n, 1:m]), "\n")
  
  # Extract the partial transport plan
  gamma <- gamma_aug[1:n, 1:m]
  
  return(gamma)
}

# Test
n <- 3
m <- 3
p <- rep(1/n, n)
q <- rep(1/m, m)
cost <- matrix(100, n, m)

result <- debug_partial_ot_tv(cost, p, q, lambda = 1)
cat("\nFinal result mass:", sum(result), "\n")