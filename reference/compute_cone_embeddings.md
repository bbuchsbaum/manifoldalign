# Compute Spectral Embeddings for CONE-Align

Compute Spectral Embeddings for CONE-Align

## Usage

``` r
compute_cone_embeddings(strata, ncomp, sigma, use_laplacian, knn)
```

## Arguments

- strata:

  List of data domains

- ncomp:

  Number of embedding dimensions

- sigma:

  Diffusion parameter

- use_laplacian:

  Whether to use Laplacian normalization

- knn:

  Number of nearest neighbors (NULL for adaptive)

## Value

List of embedding matrices
