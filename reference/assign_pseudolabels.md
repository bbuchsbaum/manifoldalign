# Assign pseudolabels based on sparse similarity matrix clustering

Takes a sparse similarity matrix and identifies clusters of highly
similar samples. Focuses on finding diverse, high-confidence cluster
representatives.

## Usage

``` r
assign_pseudolabels(
  sim_matrix,
  min_clusters = 10,
  max_clusters = 100,
  sim_threshold = NULL,
  min_cluster_size = 2,
  max_cluster_size = Inf,
  diversity_weight = 0.3,
  adaptive_threshold_quantile = 0.8,
  seed = NULL,
  use_advanced_diversity = TRUE,
  verbose = FALSE
)
```

## Arguments

- sim_matrix:

  Sparse similarity matrix (dgCMatrix or similar). Should be symmetric
  with values between 0 and 1, where higher values indicate greater
  similarity.

- min_clusters:

  Minimum number of clusters to find (default: 10)

- max_clusters:

  Maximum number of clusters to find (default: 100)

- sim_threshold:

  Minimum similarity threshold for cluster membership. If NULL, will be
  adaptively determined from the distribution of non-zero similarities.

- min_cluster_size:

  Minimum number of samples required to form a cluster (default: 2)

- max_cluster_size:

  Maximum number of samples allowed in a cluster (default: Inf)

- diversity_weight:

  Weight for promoting diversity among cluster representatives \[0,1\].
  Controls the trade-off between diversity and cluster confidence: 0 =
  select largest/most confident clusters, 1 = maximize diversity
  (minimize pairwise similarities), 0.5 = balanced approach (default:
  0.3)

- adaptive_threshold_quantile:

  Quantile of non-zero similarities to use as adaptive threshold
  (default: 0.8)

- seed:

  Random seed for reproducibility

- use_advanced_diversity:

  Whether to use advanced farthest-first traversal for diversity
  selection on large candidate sets (default: TRUE)

- verbose:

  Whether to print progress information (default: FALSE)

## Value

A list containing:

- labels:

  Factor vector of pseudolabels, with NA for unassigned samples

- representatives:

  Integer vector of row indices serving as cluster representatives

- cluster_info:

  Data frame with cluster statistics including cluster ID,
  representative index, size, and average similarity

- threshold_used:

  The similarity threshold that was applied

- n_clusters:

  Number of clusters found

## Details

The algorithm works in several steps:

1.  Determine similarity threshold (adaptive or user-specified)

2.  Build a graph of high-similarity connections

3.  Find connected components as initial clusters

4.  Select diverse representatives from clusters

5.  Optionally merge small clusters or split large ones

The diversity constraint helps ensure that cluster representatives are
spread across the similarity space, making them better anchors for
domain alignment.

## Examples

``` r
# \donttest{
# Create example sparse similarity matrix
library(Matrix)
n <- 1000
# Generate random sparse matrix with values in [0,1]
sim_matrix <- Matrix::rsparsematrix(n, n, density = 0.05, rand.x = runif)
sim_matrix <- (sim_matrix + Matrix::t(sim_matrix)) / 2  # Make symmetric
Matrix::diag(sim_matrix) <- 1  # Self-similarity = 1

# Find pseudolabels
result <- assign_pseudolabels(sim_matrix, min_clusters = 20, verbose = TRUE)
#> Processing 1000 samples with 97248 non-zero similarities
#> Adaptive threshold: 0.4069 ( 0.8 quantile of 97248 similarities)
#> After thresholding: 19450 connections remain
#> Found 1 initial clusters
#> Found 1 valid clusters after size filtering
#> Warning: Final number of representatives (1) is below min_clusters (20). Consider lowering sim_threshold or min_cluster_size.
#> Final result: 1 clusters, 1000 samples assigned
table(result$labels, useNA = "always")
#> 
#> anchor_001       <NA> 
#>       1000          0 
# }
```
