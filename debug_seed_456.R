library(manifoldalign)
library(multidesign)
library(tibble)

# Test with seed 456
set.seed(456)
X1 <- matrix(rnorm(15), 5, 3)
X2 <- matrix(rnorm(20), 5, 4)

cat("X1:\n")
print(round(X1, 3))
cat("\nX2:\n")
print(round(X2, 3))

# Check if data is degenerate
cat("\nX1 column variances:", apply(X1, 2, var), "\n")
cat("X2 column variances:", apply(X2, 2, var), "\n")

# Create hyperdesign
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

hd <- create_test_hyperdesign(list(X1, X2))

# Test FPGW
result <- fpgw(hd, omega1 = 0.2, lambda = 0.1)
cat("\nFPGW distance:", result$distances[1,2], "\n")

# Check transport plan
P <- result$transport_plans[[1]]
cat("Transport plan sum:", sum(P), "\n")
cat("Transport plan shape:", dim(P), "\n")

# Check if all entries are zero
if (sum(P) < 1e-10) {
  cat("WARNING: Transport plan is essentially zero!\n")
}