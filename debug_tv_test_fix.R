library(manifoldalign)
library(multidesign)
library(tibble)

# Reproduce the exact test scenario
set.seed(123)
X1 <- rbind(
  matrix(rnorm(15, mean = 0), 5, 3),
  matrix(rnorm(6, mean = 5), 2, 3)  # outliers
)
X2 <- rbind(
  matrix(rnorm(15, mean = 0.2), 5, 3),
  matrix(rnorm(9, mean = -5), 3, 3)  # different outliers
)

# Create appropriate design frames for each domain
design1 <- data.frame(sample_id = 1:nrow(X1))
design2 <- data.frame(sample_id = 1:nrow(X2))

# Create multidesign objects manually since domains have different sizes
md1 <- multidesign::multidesign(X1, design1)
md2 <- multidesign::multidesign(X2, design2)
hd <- hyperdesign(list(domain1 = md1, domain2 = md2))

# Run with different lambda values
result_lambda0 <- fpgw(hd, omega1 = 0.4, lambda = 0.0, max_iter = 40)
result_lambda1 <- fpgw(hd, omega1 = 0.4, lambda = 1.0, max_iter = 40)
result_lambda5 <- fpgw(hd, omega1 = 0.4, lambda = 5.0, max_iter = 40)

mass0 <- sum(result_lambda0$transport_plans[[1]])
mass1 <- sum(result_lambda1$transport_plans[[1]])
mass5 <- sum(result_lambda5$transport_plans[[1]])

cat("mass0 =", mass0, "\n")
cat("mass1 =", mass1, "\n")
cat("mass5 =", mass5, "\n")
cat("\n")

# Check the conditions
cat("abs(mass0 - mass1) < 0.1:", abs(mass0 - mass1) < 0.1, "\n")
cat("mass5 > mass1:", mass5 > mass1, "\n")