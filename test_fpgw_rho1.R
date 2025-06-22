library(manifoldalign)
library(multidesign)
library(tibble)

set.seed(1234)

toy_pair <- function(d = 2, n = 3) {
  list(
    X = matrix(rnorm(n * d), n, d),
    Y = matrix(rnorm(n * d), n, d)
  )
}

create_test_hyperdesign <- function(X_list, design = NULL) {
  n_domains <- length(X_list)
  if (is.null(design)) {
    n <- nrow(X_list[[1]])
    design <- data.frame(sample_id = 1:n)
  }
  
  md_list <- lapply(seq_along(X_list), function(i) {
    multidesign::multidesign(X_list[[i]], design)
  })
  names(md_list) <- paste0("domain", seq_along(X_list))
  
  hyperdesign(md_list)
}

pb <- toy_pair(d = 3, n = 5)
hd <- create_test_hyperdesign(list(pb$X, pb$Y))

cat("Test data dimensions:\n")
cat("X:", dim(pb$X), "\n")
cat("Y:", dim(pb$Y), "\n\n")

# Classical FGW
fgw_classical <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = NULL, max_iter = 60, verbose = TRUE)
P_classical <- fgw_classical$transport_plans[[1]]
cat("\nClassical FGW:\n")
cat("Transport plan sum:", sum(P_classical), "\n")
cat("Row sums:", rowSums(P_classical), "\n")
cat("Col sums:", colSums(P_classical), "\n")

# FPGW with rho = 1
fpgw_rho1 <- fpgw(hd, omega1 = 0.5, lambda = 0, rho = 1.0, max_iter = 60, verbose = TRUE)
P_rho1 <- fpgw_rho1$transport_plans[[1]]
cat("\nFPGW with rho=1:\n")
cat("Transport plan sum:", sum(P_rho1), "\n")
cat("Row sums:", rowSums(P_rho1), "\n")
cat("Col sums:", colSums(P_rho1), "\n")

# Check marginals used
p <- rep(1/5, 5)
q <- rep(1/5, 5)
cat("\nExpected marginals: p =", p, ", q =", q, "\n")
cat("Sum p:", sum(p), ", Sum q:", sum(q), "\n")