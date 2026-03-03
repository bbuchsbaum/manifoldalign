# Convenience: compute HKS descriptors from an adjacency matrix

Convenience: compute HKS descriptors from an adjacency matrix

## Usage

``` r
compute_hks_from_adjacency(
  A,
  k_embed = 30L,
  q = 16L,
  use_normalized_laplacian = TRUE,
  store_basis = FALSE,
  ...
)
```

## Arguments

- A:

  Square adjacency/affinity matrix (dense or sparse).

- k_embed:

  Number of eigenpairs used for the descriptor basis.

- q:

  Number of HKS time steps.

- use_normalized_laplacian:

  Logical; use the normalized Laplacian (default TRUE).

- store_basis:

  Logical; if TRUE, return \`list(descriptors=, basis=)\`.

- ...:

  Passed to \[compute_hks_descriptors()\] (e.g., time grid settings).

## Value

A numeric matrix (n x q), or a list when \`store_basis=TRUE\`.
