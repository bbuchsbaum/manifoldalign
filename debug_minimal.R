# Minimal debug for coupled_diagonalization
devtools::load_all(".")
library(multidesign)
library(multivarious)
library(tibble)

# Create tiny test case
set.seed(123)
n <- 20  # Small dataset
X1 <- matrix(rnorm(n * 3), n, 3)
X2 <- matrix(rnorm(n * 3), n, 3)

design1 <- data.frame(sample_id = 1:n)
design2 <- data.frame(sample_id = 1:n)

md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

# Run with extreme debugging
result <- tryCatch({
  coupled_diagonalization(hd, 
                         ncomp = 2,           # Very few components
                         ncomp_per_domain = 5, # Not too many eigenvectors
                         knn = 5,             # Reasonable neighbors for n=20
                         max_iter = 3,        # Just a few iterations
                         verbose = TRUE)
}, error = function(e) {
  cat("\nError occurred:\n")
  print(e)
  cat("\nTraceback:\n")
  traceback()
  NULL
}, warning = function(w) {
  cat("\nWarning:", w$message, "\n")
  invokeRestart("muffleWarning")
})

if (!is.null(result)) {
  cat("\nSUCCESS! Result obtained.\n")
  cat("Class:", class(result), "\n")
  cat("Components:", ncol(result$s), "\n")
}