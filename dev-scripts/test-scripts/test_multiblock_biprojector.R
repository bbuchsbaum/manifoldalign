#!/usr/bin/env Rscript
# Test multivarious::multiblock_biprojector

library(multivarious)

# Create dummy data
v <- matrix(rnorm(20), 10, 2)
s <- matrix(rnorm(10), 5, 2)
sdev <- c(1.5, 0.5)
preproc <- center()
block_indices <- list(block_1 = 1:5, block_2 = 6:10)

# Create biprojector
result <- multiblock_biprojector(
  v = v,
  s = s,
  sdev = sdev,
  preproc = preproc,
  block_indices = block_indices
)

cat("Result class:", class(result), "\n")
cat("Result names:", names(result), "\n")
cat("Has preproc field:", "preproc" %in% names(result), "\n")
cat("preproc value:", class(result$preproc), "\n")