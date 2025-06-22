# Debug the failing tests
library(testthat)
devtools::load_all(".")
library(multidesign)
library(multivarious)
library(tibble)
library(Matrix)

# Test 1: Basic functionality
cat("=== Test 1: Basic hyperdesign input ===\n")
set.seed(123)
n <- 50
X1 <- matrix(rnorm(n * 3), n, 3)
X2 <- matrix(rnorm(n * 4), n, 4)

design1 <- data.frame(sample_id = 1:n)
design2 <- data.frame(sample_id = 1:n)

md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

tryCatch({
  result <- coupled_diagonalization(hd, ncomp = 3, verbose = TRUE)
  cat("✓ Test 1 PASSED\n")
}, error = function(e) {
  cat("✗ Test 1 FAILED:", e$message, "\n")
})

# Test 2: Different domain sizes
cat("\n=== Test 2: Different domain sizes ===\n")
n1 <- 40
n2 <- 60
X1 <- matrix(rnorm(n1 * 3), n1, 3)
X2 <- matrix(rnorm(n2 * 4), n2, 4)

design1 <- data.frame(sample_id = 1:n1)
design2 <- data.frame(sample_id = 1:n2)

md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

tryCatch({
  suppressWarnings({
    result <- coupled_diagonalization(hd, ncomp = 2, verbose = FALSE)
  })
  cat("✓ Test 2 PASSED\n")
}, error = function(e) {
  cat("✗ Test 2 FAILED:", e$message, "\n")
})

# Test 3: Input validation - single domain
cat("\n=== Test 3: Single domain validation ===\n")
X1 <- matrix(rnorm(30), 10, 3)
md1 <- multidesign::multidesign(X1, data.frame(id = 1:10))
hd_single <- multidesign::hyperdesign(list(domain1 = md1))

tryCatch({
  result <- coupled_diagonalization(hd_single)
  cat("✗ Test 3 FAILED: Should have errored\n")
}, error = function(e) {
  if (grepl("at least 2 domains", e$message)) {
    cat("✓ Test 3 PASSED\n")
  } else {
    cat("✗ Test 3 FAILED: Wrong error -", e$message, "\n")
  }
})

# Test 4: Invalid correspondence
cat("\n=== Test 4: Invalid correspondence ===\n")
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(40), 10, 4)
md1 <- multidesign::multidesign(X1, data.frame(id = 1:10))
md2 <- multidesign::multidesign(X2, data.frame(id = 1:10))
hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

tryCatch({
  result <- coupled_diagonalization(hd, correspondence = list(matrix(1, 10, 10)))
  cat("✗ Test 4 FAILED: Should have errored\n")
}, error = function(e) {
  if (grepl("must be a list of length", e$message)) {
    cat("✓ Test 4 PASSED\n")
  } else {
    cat("✗ Test 4 FAILED: Wrong error -", e$message, "\n")
  }
})