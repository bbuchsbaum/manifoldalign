#!/usr/bin/env Rscript
# Test PARROT preprocessing with data that will show differences

library(manifoldalign)
library(multidesign)
library(tibble)

# Create test data with very different means
set.seed(123)
n <- 30
p <- 5

# Create data with large mean offset
X1 <- matrix(rnorm(n * p), n, p) + 100  # Large offset
X2 <- matrix(rnorm(n * p), n, p) + 100  # Large offset

cat("X1 mean:", mean(X1), "\n")
cat("X2 mean:", mean(X2), "\n")

# Create designs
design1 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))
design2 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))

# Create hyperdesign
md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2)
hd <- hyperdesign(list(d1 = md1, d2 = md2))

# Test with centering - use a pre_processor that works
cat("\nTesting with center preprocessing...\n")
center_proc <- multivarious::center()
cat("center_proc class:", class(center_proc), "\n")

# Apply preprocessing manually to see effect
X1_centered <- center_proc(X1)
cat("X1 centered mean:", mean(X1_centered), "\n")

# Now test PARROT
result_center <- parrot(hd, anchors = anchors, preproc = center_proc, 
                       ncomp = 4, max_iter = 20, tau = 0.1)

result_null <- parrot(hd, anchors = anchors, preproc = NULL, 
                     ncomp = 4, max_iter = 20, tau = 0.1)

# Check results
cat("\nResult with centering:\n")
cat("  Class:", class(result_center), "\n")
cat("  Has preproc:", !is.null(result_center$preproc), "\n")
cat("  Score range:", range(result_center$s), "\n")

cat("\nResult without preprocessing:\n")
cat("  Class:", class(result_null), "\n")
cat("  Has preproc:", !is.null(result_null$preproc), "\n")
cat("  Score range:", range(result_null$s), "\n")

cat("\nAre scores different?", !all(result_center$s == result_null$s), "\n")
if (all(result_center$s == result_null$s)) {
  cat("Max absolute difference:", max(abs(result_center$s - result_null$s)), "\n")
}