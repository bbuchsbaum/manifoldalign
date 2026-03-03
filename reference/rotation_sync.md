# Rotation synchronization for orthogonal group O(k)

Rotation synchronization for orthogonal group O(k)

## Usage

``` r
rotation_sync(G, W, E, N, k)
```

## Arguments

- G:

  list of align_transform objects (group="O"), each a relative g_ij

- W:

  numeric weights per edge in E

- E:

  data.frame with columns i, j (1-based dataset indices)

- N:

  number of datasets

- k:

  latent dimension

## Value

list of per-domain align_transform objects mapping to consensus

## Examples

``` r
# \donttest{
# Rotation synchronization for orthogonal transformations
k <- 3
N <- 3
G <- list(matrix(rnorm(9), k, k), matrix(rnorm(9), k, k))
W <- c(1, 1)
E <- data.frame(i = c(1, 2), j = c(2, 3))
# result <- rotation_sync(G, W, E, N = N, k = k)
# }
```
