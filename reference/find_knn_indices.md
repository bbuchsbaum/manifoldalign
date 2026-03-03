# Compute K-nearest neighbours under supported metrics

Compute K-nearest neighbours under supported metrics

## Usage

``` r
find_knn_indices(query, data, k = 5, metric = "euclidean")
```

## Value

A list with elements \`idx\` (integer matrix of neighbor indices) and
\`dists\` (numeric matrix of distances).
