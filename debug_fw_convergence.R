library(manifoldalign)
library(multidesign)
library(tibble)

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

# Test with problematic seed
set.seed(42)
pb <- toy_pair(d = 3, n = 5)

# Look at the data
cat("X data:\n")
print(pb$X)
cat("\nY data:\n")
print(pb$Y)

# Compute distance matrices
Dx <- as.matrix(dist(pb$X))
Dy <- as.matrix(dist(pb$Y))
cat("\nDistance matrix X (normalized):\n")
print(round(Dx / max(Dx), 3))
cat("\nDistance matrix Y (normalized):\n")
print(round(Dy / max(Dy), 3))

# Compute feature cost
C_feat <- manifoldalign:::compute_feature_cost(pb$X, pb$Y, "euclidean")
cat("\nFeature cost matrix:\n")
print(round(C_feat, 3))

# Try with very simple initialization
p <- rep(1/5, 5)
q <- rep(1/5, 5)
P0 <- outer(p, q)

# Run FW with debugging
cat("\nRunning Frank-Wolfe with debugging...\n")
sol <- manifoldalign:::fw_fpgw(
  C_feat, Dx, Dy, p, q,
  omega1 = 0.5, lambda = 0, rho = NULL,
  P0 = P0, max_iter = 10, tol = 1e-6, verbose = TRUE
)

cat("\nFinal transport plan:\n")
print(sol$P)
cat("Sum:", sum(sol$P), "\n")

# Check objective trace
obj_trace <- attr(sol, "obj_trace")
cat("\nObjective values:\n")
print(obj_trace)