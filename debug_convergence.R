# Debug convergence test
devtools::load_all(".")
library(multidesign)
library(multivarious)
library(tibble)

set.seed(123)

# Create well-separated data for good convergence (from test)
n <- 30
X1 <- rbind(
  matrix(rnorm(15 * 3, mean = -2), 15, 3),
  matrix(rnorm(15 * 3, mean = 2), 15, 3)
)
X2 <- rbind(
  matrix(rnorm(15 * 4, mean = -2), 15, 4),
  matrix(rnorm(15 * 4, mean = 2), 15, 4)
)

design1 <- data.frame(sample_id = 1:n)
design2 <- data.frame(sample_id = 1:n)

md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)

hd <- multidesign::hyperdesign(list(domain1 = md1, domain2 = md2))

# Run with tight tolerance and verbose output
cat("Running coupled_diagonalization with verbose output...\n")
result <- coupled_diagonalization(hd, 
                                 ncomp = 2,
                                 max_iter = 100,
                                 tol = 1e-8,
                                 verbose = TRUE)

cat("\n\nResults:\n")
cat("Converged:", result$converged, "\n")
cat("Iterations:", result$iterations, "\n")
cat("Final cost:", result$final_cost, "\n")

# Let's try with different parameters
cat("\n\nTrying with looser tolerance...\n")
result2 <- coupled_diagonalization(hd, 
                                  ncomp = 2,
                                  max_iter = 100,
                                  tol = 1e-4,  # Looser tolerance
                                  verbose = TRUE)

cat("\n\nResults with looser tolerance:\n")
cat("Converged:", result2$converged, "\n")
cat("Iterations:", result2$iterations, "\n")
cat("Final cost:", result2$final_cost, "\n")