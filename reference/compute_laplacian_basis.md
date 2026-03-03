# Compute a Laplacian spectral basis (eigenpairs) from an adjacency matrix

Compute a Laplacian spectral basis (eigenpairs) from an adjacency matrix

## Usage

``` r
compute_laplacian_basis(
  A,
  k,
  normalized = TRUE,
  which = "SA",
  eig_tol = 1e-12,
  extra = 4L,
  seed = 1L,
  tol = 1e-06,
  maxitr = 1000
)
```

## Arguments

- A:

  Square adjacency/affinity matrix (dense or sparse).

- k:

  Integer number of non-trivial eigenpairs to return.

- normalized:

  Logical; use the normalized Laplacian (default TRUE).

- which:

  Which eigenpairs to request (passed to the eigensolver). Default is
  \`"SA"\` (smallest algebraic) for Laplacian bases.

- eig_tol:

  Positive threshold; eigenvalues \`\<= eig_tol\` are treated as trivial
  and dropped.

- extra:

  Integer number of extra eigenpairs requested to survive multiple
  near-zero modes on disconnected graphs.

- seed:

  Integer seed used to build a deterministic initial vector for RSpectra
  (helps reproducibility).

- tol:

  Numeric tolerance passed to RSpectra/PRIMME when available.

- maxitr:

  Maximum iterations passed to RSpectra when available.

## Value

A list with fields \`values\` (length k) and \`vectors\` (n x k).
