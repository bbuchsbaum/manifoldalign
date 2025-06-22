library(manifoldalign)
library(multidesign)
library(tibble)

# Reproduce the test case from test-fpgw-validation.R
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
cat("Testing with lambda = 0.0...\n")
result_lambda0 <- fpgw(hd, omega1 = 0.4, lambda = 0.0, max_iter = 40)
mass0 <- sum(result_lambda0$transport_plans[[1]])
cat("Mass with lambda=0:", mass0, "\n\n")

cat("Testing with lambda = 1.0...\n")
result_lambda1 <- fpgw(hd, omega1 = 0.4, lambda = 1.0, max_iter = 40)
mass1 <- sum(result_lambda1$transport_plans[[1]])
cat("Mass with lambda=1:", mass1, "\n\n")

cat("Testing with lambda = 5.0...\n")
result_lambda5 <- fpgw(hd, omega1 = 0.4, lambda = 5.0, max_iter = 40)
mass5 <- sum(result_lambda5$transport_plans[[1]])
cat("Mass with lambda=5:", mass5, "\n\n")

cat("Summary:\n")
cat("mass0 =", mass0, "\n")
cat("mass1 =", mass1, "\n")
cat("mass5 =", mass5, "\n")
cat("mass0 > mass1?", mass0 > mass1, "\n")
cat("mass1 > mass5?", mass1 > mass5, "\n")

# Let's also check what the gradient looks like in fw_fpgw
cat("\nChecking gradient for lambda=1 case:\n")
p <- rep(1/nrow(X1), nrow(X1))
q <- rep(1/nrow(X2), nrow(X2))
C_feat <- manifoldalign:::compute_feature_cost(X1, X2, "euclidean")
cat("Feature cost range:", range(C_feat), "\n")

# Check if the gradient + lambda is making transport too expensive
grad_range <- range(result_lambda1$grad)
cat("Gradient range:", grad_range, "\n")
cat("Gradient + lambda range:", range(result_lambda1$grad + 1.0), "\n")