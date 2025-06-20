#!/usr/bin/env Rscript
# Fix PARROT preprocessing test

library(manifoldalign)

# The key issue is that the test expects:
# 1. result_center$preproc to be non-NULL when preproc is provided
# 2. The results to be different with/without preprocessing

# Looking at the code, when proc is NULL (which happens for functions),
# we create the result manually without preproc field.
# When proc is not NULL, we use multiblock_biprojector which should preserve it.

# The fix is to ensure center() creates a proper pre_processor that gets preserved.
# Let's check what happens with the current test data:

library(multidesign)
library(tibble)

set.seed(42)
n <- 25

# Create test data
X1 <- matrix(rnorm(n * 3), n, 3)
X2 <- matrix(rnorm(n * 3), n, 3) 

# Add significant mean shift to make preprocessing matter
X1 <- X1 + 50
X2 <- X2 + 50

design1 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))
design2 <- tibble(node_id = 1:n, anchors = c(1:5, rep(NA, n-5)))

md1 <- multidesign(X1, design1)
md2 <- multidesign(X2, design2) 
hd <- hyperdesign(list(d1 = md1, d2 = md2))

# Test preprocessing
preproc_obj <- multivarious::center()
cat("Preprocessor class:", class(preproc_obj), "\n")
cat("Is it a pre_processor?", inherits(preproc_obj, "pre_processor"), "\n")

# The issue is that center() returns a "prepper" not "pre_processor"
# Let's see what multivarious expects
result <- parrot(hd, anchors = anchors, preproc = preproc_obj, ncomp = 4, max_iter = 20)
cat("\nResult preproc field:", class(result$preproc), "\n")