# Select diverse cluster representatives

Internal function to select cluster representatives that are diverse
across the similarity space.

## Usage

``` r
.select_diverse_representatives(
  sim_matrix,
  cluster_assignments,
  valid_clusters,
  diversity_weight,
  min_clusters,
  max_clusters,
  use_advanced_diversity = TRUE,
  verbose
)
```

## Arguments

- sim_matrix:

  Original similarity matrix

- cluster_assignments:

  Vector of cluster assignments

- valid_clusters:

  Vector of valid cluster IDs

- diversity_weight:

  Weight for diversity constraint

- min_clusters:

  Minimum number of clusters

- max_clusters:

  Maximum number of clusters

- use_advanced_diversity:

  Whether to use advanced FFS algorithm

- verbose:

  Whether to print progress

## Value

Vector of representative indices
