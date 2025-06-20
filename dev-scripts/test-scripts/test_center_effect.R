#!/usr/bin/env Rscript
# Test if centering actually changes the data

library(multivarious)

# Create data with non-zero mean
X <- matrix(rnorm(20, mean = 10), 10, 2)
cat("Original mean:", mean(X), "\n")

# Apply centering
center_fn <- center()
X_centered <- center_fn(X)
cat("Centered mean:", mean(X_centered), "\n")
cat("Are they different?", !all(X == X_centered), "\n")