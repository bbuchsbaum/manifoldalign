# Worker function to compute a single cone embedding

Worker function to compute a single cone embedding

## Usage

``` r
compute_embedding(domain, ncomp, sigma, use_laplacian, knn)
```

## Arguments

- domain:

  Data domain

- ncomp:

  Number of embedding dimensions

- sigma:

  Diffusion parameter

- use_laplacian:

  Whether to use Laplacian normalization

- knn:

  Number of nearest neighbors (NULL for adaptive)

## Value

Embedding matrix
