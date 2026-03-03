# Compute Distance Matrix

Memory-efficient distance computation that avoids creating full distance
matrix

## Usage

``` r
compute_distance_matrix(X, metric = "euclidean", chunk_size = 1000)
```

## Arguments

- X:

  Data matrix (samples x features)

- metric:

  Distance metric to use (default: "euclidean")

- chunk_size:

  Chunk size for memory-efficient computation

## Value

A symmetric numeric distance matrix of dimension n x n.
