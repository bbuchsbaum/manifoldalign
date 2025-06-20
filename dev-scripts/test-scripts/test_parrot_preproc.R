#!/usr/bin/env Rscript
# Test PARROT preprocessing

library(manifoldalign)
library(multidesign)
library(tibble)

# Create test data with non-zero mean
set.seed(123)
n <- 30
p <- 5

# Create data with different means
X1 <- matrix(rnorm(n * p, mean = 10), n, p)
X2 <- matrix(rnorm(n * p, mean = 10), n, p)

cat("X1 mean before preprocessing:", mean(X1), "\n")
cat("X2 mean before preprocessing:", mean(X2), "\n")

# Create designs
design1 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))
design2 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))

# Create hyperdesign
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)
hd <- hyperdesign(list(d1 = md1, d2 = md2))

# Test with centering
cat("\nRunning PARROT with centering...\n")
result_center <- parrot(hd, anchors = anchors, preproc = multivarious::center(), 
                       ncomp = 4, max_iter = 20)

# Test without preprocessing
cat("\nRunning PARROT without preprocessing...\n")
result_null <- parrot(hd, anchors = anchors, preproc = NULL, 
                     ncomp = 4, max_iter = 20)

# Check results
cat("\nChecking results:\n")
cat("result_center has preproc:", !is.null(result_center$preproc), "\n")
cat("result_null has preproc:", !is.null(result_null$preproc), "\n")

# Check if scores differ
cat("\nScores identical?", all(result_center$s == result_null$s), "\n")
cat("Max difference in scores:", max(abs(result_center$s - result_null$s)), "\n")

# Also check the transport plans
cat("Transport plans identical?", all(result_center$transport_plan == result_null$transport_plan), "\n")
cat("Max difference in transport plans:", max(abs(result_center$transport_plan - result_null$transport_plan)), "\n")