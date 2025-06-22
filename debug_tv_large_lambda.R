library(manifoldalign)
library(multidesign)
library(tibble)

# Test setup
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)
)

design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

cat("Testing with very large lambda values:\n\n")

for (lambda in c(0, 10, 50, 100, 500, 1000)) {
  result <- fpgw(hd, omega1 = 0.4, lambda = lambda, max_iter = 100, verbose = FALSE)
  P <- result$transport_plans[[1]]
  mass <- sum(P)
  
  cat("Lambda =", lambda, "-> Mass =", mass, ", Converged =", result$converged[1,2], "\n")
}