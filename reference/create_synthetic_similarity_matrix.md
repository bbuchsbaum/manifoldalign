# Create synthetic similarity matrix for testing pseudolabeling

Generates a synthetic sparse similarity matrix with known cluster
structure. Used for testing and demonstrating the pseudolabeling
functions.

## Usage

``` r
create_synthetic_similarity_matrix(
  n_samples = 1000,
  n_clusters = 20,
  within_cluster_sim = 0.8,
  between_cluster_sim = 0.1,
  sparsity = 0.1,
  noise_level = 0.1
)
```

## Arguments

- n_samples:

  Total number of samples

- n_clusters:

  Number of true clusters

- within_cluster_sim:

  Average similarity within clusters (default: 0.8)

- between_cluster_sim:

  Average similarity between clusters (default: 0.1)

- sparsity:

  Overall sparsity level (proportion of non-zero entries) (default: 0.1)

- noise_level:

  Amount of noise to add to similarities (default: 0.1)

## Value

A list containing:

- sim_matrix:

  Sparse similarity matrix (dgCMatrix) with synthetic cluster structure

- true_labels:

  Factor vector of true cluster assignments for validation

- cluster_centers:

  Integer vector of cluster center indices

## Examples

``` r
# \donttest{
# Create synthetic data
synthetic <- create_synthetic_similarity_matrix(n_samples = 500, 
                                                n_clusters = 10)

# Apply pseudolabeling
result <- assign_pseudolabels(synthetic$sim_matrix, verbose = TRUE)
#> Processing 500 samples with 9290 non-zero similarities
#> Adaptive threshold: 0.8682 ( 0.8 quantile of 9290 similarities)
#> 'as(<dsCMatrix>, "dgCMatrix")' is deprecated.
#> Use 'as(., "generalMatrix")' instead.
#> See help("Deprecated") and help("Matrix-deprecated").
#> After thresholding: 1858 connections remain
#> Found 17 initial clusters
#> Found 10 valid clusters after size filtering
#> Final result: 10 clusters, 493 samples assigned

# Compare with true labels
table(result$labels, synthetic$true_labels, useNA = "always")
#>             
#>              cluster_1 cluster_10 cluster_2 cluster_3 cluster_4 cluster_5
#>   anchor_001        49          0         0         0         0         0
#>   anchor_002         0          0        49         0         0         0
#>   anchor_003         0          0         0        50         0         0
#>   anchor_004         0          0         0         0        50         0
#>   anchor_005         0          0         0         0         0        48
#>   anchor_006         0          0         0         0         0         0
#>   anchor_007         0          0         0         0         0         0
#>   anchor_008         0          0         0         0         0         0
#>   anchor_009         0          0         0         0         0         0
#>   anchor_010         0         49         0         0         0         0
#>   <NA>               1          1         1         0         0         2
#>             
#>              cluster_6 cluster_7 cluster_8 cluster_9 <NA>
#>   anchor_001         0         0         0         0    0
#>   anchor_002         0         0         0         0    0
#>   anchor_003         0         0         0         0    0
#>   anchor_004         0         0         0         0    0
#>   anchor_005         0         0         0         0    0
#>   anchor_006        50         0         0         0    0
#>   anchor_007         0        50         0         0    0
#>   anchor_008         0         0        50         0    0
#>   anchor_009         0         0         0        48    0
#>   anchor_010         0         0         0         0    0
#>   <NA>               0         0         0         2    0
# }
```
