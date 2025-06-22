library(manifoldalign)
library(multidesign)
library(tibble)

toy_pair <- function(d = 2, n = 3) {
  list(
    X = matrix(rnorm(n * d), n, d),
    Y = matrix(rnorm(n * d), n, d)
  )
}

create_test_hyperdesign <- function(X_list, design = NULL) {
  n_domains <- length(X_list)
  if (is.null(design)) {
    n <- nrow(X_list[[1]])
    design <- data.frame(sample_id = 1:n)
  }
  
  md_list <- lapply(seq_along(X_list), function(i) {
    multidesign::multidesign(X_list[[i]], design)
  })
  names(md_list) <- paste0("domain", seq_along(X_list))
  
  hyperdesign(md_list)
}

# Test with problematic seed
set.seed(42)
pb <- toy_pair(d = 3, n = 5)
hd <- create_test_hyperdesign(list(pb$X, pb$Y))

# Run with different iteration counts
for (max_iter in c(10, 30, 60, 100, 200)) {
  result <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = NULL, max_iter = max_iter, tol = 1e-6, verbose = FALSE)
  P <- result$transport_plans[[1]]
  
  cat("\nmax_iter =", max_iter, "\n")
  cat("Transport plan sum:", sum(P), "\n")
  cat("Converged:", result$converged[1,2], "\n")
  
  # Check if it's a valid transport plan
  row_sums <- rowSums(P)
  col_sums <- colSums(P)
  cat("Row sums valid:", all(abs(row_sums - 0.2) < 1e-6), "\n")
  cat("Col sums valid:", all(abs(col_sums - 0.2) < 1e-6), "\n")
  
  # Look at the transport plan structure
  cat("Number of non-zero entries:", sum(P > 1e-10), "\n")
}

# Try with better initialization
cat("\n\nTrying with entropic warm-start:\n")
result_warm <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = NULL, 
                    epsilon = 0.1, max_iter = 100, tol = 1e-6, verbose = FALSE)
P_warm <- result_warm$transport_plans[[1]]
cat("Transport plan sum:", sum(P_warm), "\n")
cat("Converged:", result_warm$converged[1,2], "\n")