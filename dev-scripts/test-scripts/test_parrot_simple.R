# Simple test of the improved PARROT implementation
library(manifoldalign)

# Create simple test networks
set.seed(123)
n_nodes <- 20

# Create structured node features
angles <- 2 * pi * (1:n_nodes) / n_nodes
X1 <- cbind(
  2 * cos(angles),
  2 * sin(angles),
  0.5 * (1:n_nodes) / n_nodes
)

# Create similar second network with small perturbation
X2 <- X1 + matrix(rnorm(n_nodes * 3, 0, 0.1), n_nodes, 3)

# Create anchor correspondences (first 5 nodes)
anchors1 <- c(1:5, rep(NA, n_nodes - 5))
anchors2 <- c(1:5, rep(NA, n_nodes - 5))

# Create hyperdesign structure
domain1 <- list(
  x = X1,
  design = data.frame(
    node_id = 1:n_nodes,
    anchors = anchors1
  )
)

domain2 <- list(
  x = X2,
  design = data.frame(
    node_id = 1:n_nodes,
    anchors = anchors2
  )
)

hd <- list(domain1 = domain1, domain2 = domain2)
class(hd) <- c("hyperdesign", "list")

# Test PARROT with improved algorithm and new parameters
cat("Testing PARROT with second-round improvements...\n")

result <- tryCatch({
  parrot(hd, anchors = anchors, 
         ncomp = 5, 
         sigma = 0.2,    # RWR restart probability
         lambda = 0.1,   # Consistency regularization
         tau = 0.05,     # Entropy regularization
         alpha = 0.6,    # Weight for feature vs RWR cost
         gamma = 0.15,   # Cross-graph mixing parameter
         solver = "sinkhorn",
         max_iter = 15,  # Conservative iteration limit
         tol = 1e-4)
}, error = function(e) {
  cat("Error occurred:", e$message, "\n")
  cat("Traceback:\n")
  traceback()
  return(NULL)
})

if (!is.null(result)) {
  cat("✓ PARROT completed successfully!\n")
  cat("Transport plan dimensions:", dim(result$transport_plan), "\n")
  cat("Transport plan row sums (should be ~1):", range(rowSums(result$transport_plan)), "\n")
  cat("Transport plan col sums (should be ~1):", range(colSums(result$transport_plan)), "\n")
  
  # Check if anchors have high transport mass
  anchor_mass <- numeric(5)
  for (i in 1:5) {
    anchor_mass[i] <- result$transport_plan[i, i]
  }
  cat("Anchor diagonal masses:", round(anchor_mass, 3), "\n")
  cat("Average anchor mass:", round(mean(anchor_mass), 3), "\n")
  
  # Check return structure
  cat("Result components:", names(result), "\n")
  cat("Embedding dimensions:", dim(result$s), "\n")
  
} else {
  cat("✗ PARROT test failed\n")
}