library(manifoldalign)
library(multidesign)
library(tibble)

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

# Test different dimensions
set.seed(42)  # Use same seed as validation test
X1 <- matrix(rnorm(15), 5, 3)  # 5 points in R^3
X2 <- matrix(rnorm(20), 5, 4)  # 5 points in R^4

cat("X1 shape:", dim(X1), "\n")
cat("X2 shape:", dim(X2), "\n")

hd <- create_test_hyperdesign(list(X1, X2))

result <- fpgw(hd, omega1 = 0.2, lambda = 0.1)

cat("\nDistance matrix:\n")
print(result$distances)

cat("\nDistance [1,2]:", result$distances[1, 2], "\n")
cat("Is positive?", result$distances[1, 2] > 0, "\n")