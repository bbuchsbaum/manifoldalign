library(manifoldalign)

# Debug step size computation with TV penalty
set.seed(1)
n <- 4
X <- matrix(rnorm(n*2), n)
Y <- matrix(rnorm(n*2), n)
Cx <- as.matrix(dist(X))
Cy <- as.matrix(dist(Y))
C <- as.matrix(dist(rbind(X,Y)))[1:n, n+seq_len(n)]
p <- q <- rep(1/n, n)

# Patch fw_fpgw to add debug output
fw_fpgw_debug <- function(C, Cx, Cy, p, q, omega1, lambda, rho, P0, max_iter, ...) {
  cat("\n=== Debug fw_fpgw with lambda =", lambda, "===\n")
  
  n1 <- length(p)
  n2 <- length(q)
  P <- P0
  
  # Pre-compute constants
  Cx2 <- Cx^2
  Cy2 <- Cy^2
  feature_weight <- omega1
  struct_weight <- 1 - omega1
  
  for (iter in 1:min(3, max_iter)) {  # Just first 3 iterations
    cat("\nIteration", iter, ":\n")
    
    # Gradient computation
    const1_P <- as.vector(Cx2 %*% rowSums(P))
    const2_P <- as.vector(Cy2 %*% colSums(P))
    const_P <- outer(const1_P, rep(1, n2)) + outer(rep(1, n1), const2_P)
    M_P <- const_P - 2 * Cx %*% P %*% t(Cy)
    grad <- feature_weight * C + struct_weight * 2 * M_P
    
    # Current objective
    current_obj <- sum(P * (feature_weight * C + struct_weight * M_P))
    if (lambda > 0) {
      tv_penalty <- lambda * (sum(p)^2 + sum(q)^2 - 2 * sum(P)^2)
      current_obj <- current_obj + tv_penalty
      cat("  TV penalty term:", tv_penalty, "\n")
    }
    cat("  Current objective:", current_obj, "\n")
    cat("  Current mass:", sum(P), "\n")
    
    # Oracle
    G <- grad - min(grad)
    gamma_star <- manifoldalign:::classical_ot_lp(G, p, q)
    
    # Direction
    dP <- gamma_star - P
    cat("  mass(gamma_star):", sum(gamma_star), "\n")
    cat("  mass(dP):", sum(dP), "\n")
    
    # Step size computation
    dr <- rowSums(dP)
    dc <- colSums(dP)
    term1 <- sum((Cx2 %*% dr) * dr)
    term2 <- sum((Cy2 %*% dc) * dc)
    term3 <- -2 * sum((Cx %*% dP %*% t(Cy)) * dP)
    
    a <- struct_weight * (term1 + term2 + term3)
    b <- sum(grad * dP)
    
    cat("  Base a:", a, ", b:", b, "\n")
    
    if (lambda > 0) {
      mass_P <- sum(P)
      mass_dP <- sum(dP)
      a_tv <- 2 * lambda * mass_dP^2
      b_tv <- -4 * lambda * mass_P * mass_dP
      a <- a + a_tv
      b <- b + b_tv
      cat("  TV contribution: a_tv =", a_tv, ", b_tv =", b_tv, "\n")
      cat("  Total a:", a, ", b:", b, "\n")
    }
    
    # Compute step size
    if (abs(a) < 1e-10) {
      alpha <- if (b < 0) 1 else 0
    } else if (a > 0) {
      alpha_star <- -b / (2 * a)
      alpha <- min(1, max(0, alpha_star))
    } else {
      alpha <- if (a + b < 0) 1 else 0
    }
    
    cat("  Step size alpha:", alpha, "\n")
    
    # Update
    P <- P + alpha * dP
    cat("  New mass:", sum(P), "\n")
  }
  
  return(sum(P))
}

# Test with different lambdas
for (lambda in c(0, 1, 5)) {
  mass <- fw_fpgw_debug(C, Cx, Cy, p, q, omega1 = 0.4, lambda = lambda, 
                        rho = NULL, P0 = outer(p, q), max_iter = 40)
  cat("\nFinal mass with lambda =", lambda, ":", mass, "\n")
  cat(paste(rep("=", 50), collapse=""), "\n")
}