library(manifoldalign)
library(multidesign)
library(tibble)

# Simple test case
set.seed(123)
X1 <- matrix(rnorm(12), 4, 3)
X2 <- matrix(rnorm(12), 4, 3)
design <- data.frame(id = 1:4)

md1 <- multidesign(X1, design)
md2 <- multidesign(X2, design)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Test warning
cat("Testing lambda = 1.0:\n")
tryCatch({
  result <- fpgw(hd, omega1 = 0.4, lambda = 1.0)
  cat("No warning produced!\n")
}, warning = function(w) {
  cat("Warning produced:", conditionMessage(w), "\n")
})