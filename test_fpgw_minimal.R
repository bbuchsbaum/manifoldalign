library(manifoldalign)
library(multidesign)
library(tibble)

# Minimal test with clear structure
set.seed(42)
n <- 6

# Domain 1: 2D with clear clusters
X1 <- rbind(
  matrix(c(0, 0, 0.1, 0, 0, 0.1), ncol = 2, byrow = TRUE),  # Cluster 1
  matrix(c(3, 3, 3.1, 3, 3, 3.1), ncol = 2, byrow = TRUE)   # Cluster 2
)

# Domain 2: 2D with same structure
X2 <- rbind(
  matrix(c(0, 0, 0, 0.1, 0.1, 0), ncol = 2, byrow = TRUE),  # Cluster 1
  matrix(c(2, 2, 2, 2.1, 2.1, 2), ncol = 2, byrow = TRUE)   # Cluster 2  
)

design <- data.frame(id = 1:n, cluster = rep(1:2, each = 3))

hd <- hyperdesign(list(
  visual = multidesign(X1, design),
  semantic = multidesign(X2, design)
))

# Run with pure structural alignment
cat("Testing FPGW on minimal structured data:\n")
cat("======================================\n\n")

# Should converge quickly with good structure
result <- fpgw(hd, omega1 = 0.1, verbose = TRUE, max_iter = 20, tol = 1e-4)

P <- result$transport_plans[[1]]
cat("\nTransport plan:\n")
print(round(P, 3))

cat("\nRow sums:", round(rowSums(P), 3), "\n")
cat("Col sums:", round(colSums(P), 3), "\n")
cat("Total mass:", sum(P), "\n")

# Check if it correctly matched clusters
cat("\nCluster matching (sum of transport within same cluster):\n")
within_cluster_mass <- sum(P[1:3, 1:3]) + sum(P[4:6, 4:6])
cat("Within-cluster mass:", within_cluster_mass, "\n")
cat("Between-cluster mass:", 1 - within_cluster_mass, "\n")