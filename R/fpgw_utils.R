#' FPGW utility functions
#' @keywords internal

#' Project onto doubly stochastic matrices
#' @return A doubly stochastic matrix with row marginals `p` and column marginals `q`.
#' @keywords internal
project_onto_simplex_doubly_stochastic <- function(P, p, q, max_iter = 100, tol = 1e-9) {
  # Use Sinkhorn-Knopp algorithm to project onto doubly stochastic matrices
  # with specified marginals p and q
  
  n <- length(p)
  m <- length(q)
  
  # Ensure P is non-negative
  P[P < 0] <- 0
  
  # Add small epsilon to avoid division by zero
  P <- P + 1e-20
  
  # Sinkhorn-Knopp iterations
  for (iter in seq_len(max_iter)) {
    # Row normalization
    row_sums <- rowSums(P)
    P <- sweep(P, 1, p / row_sums, "*")
    
    # Column normalization
    col_sums <- colSums(P)
    P <- sweep(P, 2, q / col_sums, "*")
    
    # Check convergence
    if (max(abs(rowSums(P) - p)) < tol && max(abs(colSums(P) - q)) < tol) {
      break
    }
  }
  
  return(P)
}
