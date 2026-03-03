# Convenience: compute diffusion coordinates from an adjacency matrix

Convenience: compute diffusion coordinates from an adjacency matrix

## Usage

``` r
compute_diffusion_coordinates_from_adjacency(
  A,
  k_embed = 30L,
  time = 1,
  use_normalized_laplacian = TRUE,
  store_basis = FALSE,
  ...
)
```

## Arguments

- A:

  Square adjacency/affinity matrix (dense or sparse).

- k_embed:

  Number of eigenpairs used for the basis.

- time:

  Positive scalar or vector of times.

- use_normalized_laplacian:

  Logical; use the normalized Laplacian (default TRUE).

- store_basis:

  Logical; if TRUE, return \`list(coords=, basis=)\`.

- ...:

  Passed to \[compute_diffusion_coordinates()\].

## Value

A numeric matrix (n x k) or list of matrices, or a list when
\`store_basis=TRUE\`.
