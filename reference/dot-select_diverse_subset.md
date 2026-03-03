# Select diverse subset of representatives

ENHANCED: Use a greedy algorithm with continuous weighting between
diversity and cluster confidence scores.

## Usage

``` r
.select_diverse_subset(
  sim_matrix,
  candidates,
  target_size,
  diversity_weight,
  cluster_assignments = NULL,
  valid_clusters = NULL
)
```

## Arguments

- sim_matrix:

  Similarity matrix

- candidates:

  Vector of candidate representative indices

- target_size:

  Target number of representatives

- diversity_weight:

  Weight for diversity vs confidence \[0,1\]. 0 = pure confidence
  (largest clusters), 1 = pure diversity, 0.5 = balanced

- cluster_assignments:

  Vector of cluster assignments (for confidence scores)

- valid_clusters:

  Vector of valid cluster IDs (for confidence scores)

## Value

Vector of selected representative indices
