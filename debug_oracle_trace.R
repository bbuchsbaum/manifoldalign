library(manifoldalign)

# Trace through partial_ot_tv
n <- 3
m <- 3
p <- rep(1/n, n)
q <- rep(1/m, m)
cost <- matrix(100, n, m)

cat("Inputs to partial_ot_tv:\n")
cat("cost:\n")
print(cost)
cat("p:", p, "\n")
cat("q:", q, "\n")

# Manually run through the function
C_aug <- matrix(0, nrow = n + 1, ncol = m + 1)
C_aug[1:n, 1:m] <- cost

cat("\nAugmented cost matrix C_aug:\n")
print(C_aug)

p_aug <- c(p, sum(q))
q_aug <- c(q, sum(p))

cat("\np_aug:", p_aug, "\n")
cat("q_aug:", q_aug, "\n")

# Now solve manually with lpSolve
library(lpSolve)

n_aug <- length(p_aug)
m_aug <- length(q_aug)
f.obj <- as.vector(C_aug)

cat("\nObjective coefficients (first 16):\n")
print(f.obj)

# Build constraints
row_constraints <- matrix(0, n_aug, n_aug * m_aug)
for (i in 1:n_aug) {
  idx <- ((i-1)*m_aug + 1):(i*m_aug)
  row_constraints[i, idx] <- 1
}

col_constraints <- matrix(0, m_aug, n_aug * m_aug)
for (j in 1:m_aug) {
  idx <- seq(j, n_aug*m_aug, by = m_aug)
  col_constraints[j, idx] <- 1
}

Amat <- rbind(row_constraints, col_constraints)
rhs <- c(p_aug, q_aug)

cat("\nConstraint matrix shape:", dim(Amat), "\n")
cat("RHS length:", length(rhs), "\n")

# Solve
lp_result <- lp("min", f.obj, Amat, "=", rhs)

cat("\nLP status:", lp_result$status, "\n")
cat("Objective value:", lp_result$objval, "\n")

gamma_aug <- matrix(lp_result$solution, n_aug, m_aug)
cat("\nSolution:\n")
print(round(gamma_aug, 4))

cat("\nMass in top-left block:", sum(gamma_aug[1:n, 1:m]), "\n")