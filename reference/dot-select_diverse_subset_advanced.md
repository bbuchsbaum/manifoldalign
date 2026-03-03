# Advanced diversity selection with farthest-first traversal

OPTIMIZED: Implements efficient farthest-first traversal (FFS) for
diversity selection with improved O(k^2) complexity through early
termination and sparse matrix optimizations.

## Usage

``` r
.select_diverse_subset_advanced(
  sim_matrix,
  candidates,
  target_size,
  diversity_weight,
  cluster_assignments = NULL,
  valid_clusters = NULL,
  use_advanced = TRUE
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

  Weight for diversity vs confidence \[0,1\]

- cluster_assignments:

  Vector of cluster assignments (for confidence scores)

- valid_clusters:

  Vector of valid cluster IDs (for confidence scores)

- use_advanced:

  Whether to use advanced FFS algorithm (default: TRUE)

## Value

Vector of selected representative indices
