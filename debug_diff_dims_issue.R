library(manifoldalign)
library(multidesign)
library(tibble)

# Reproduce the test case
set.seed(42)
X1 <- matrix(rnorm(15), 5, 3)  # 5 points in R^3
X2 <- matrix(rnorm(20), 5, 4)  # 5 points in R^4

cat("X1 dimensions:", dim(X1), "\n")
cat("X2 dimensions:", dim(X2), "\n")

# Create test hyperdesign
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

# Test compute_feature_cost directly
cat("\nTesting compute_feature_cost with different dimensions:\n")
tryCatch({
  C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
  cat("Feature cost computed successfully\n")
  cat("C_feat shape:", dim(C_feat), "\n")
  cat("C_feat range:", range(C_feat), "\n")
}, error = function(e) {
  cat("Error in compute_feature_cost:", e$message, "\n")
})

# Try FPGW
cat("\nTesting FPGW:\n")
tryCatch({
  result <- fpgw(hd, omega1 = 0.2, lambda = 0.1)
  cat("FPGW completed\n")
  cat("Distance:", result$distances[1,2], "\n")
}, error = function(e) {
  cat("Error in FPGW:", e$message, "\n")
})