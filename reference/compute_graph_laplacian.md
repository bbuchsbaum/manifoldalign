# Compute a (normalized) graph Laplacian from an adjacency matrix

Compute a (normalized) graph Laplacian from an adjacency matrix

## Usage

``` r
compute_graph_laplacian(
  A,
  normalized = TRUE,
  symmetrize = TRUE,
  zero_diagonal = TRUE
)
```

## Arguments

- A:

  Square adjacency/affinity matrix (dense or sparse). Values should be
  nonnegative; the diagonal is set to 0 by default.

- normalized:

  Logical; if TRUE (default), compute the normalized Laplacian \\L = I -
  D^{-1/2} A D^{-1/2}\\. If FALSE, compute the unnormalized Laplacian
  \\L = D - A\\.

- symmetrize:

  Logical; if TRUE (default), symmetrize via \`(A + t(A))/2\`.

- zero_diagonal:

  Logical; if TRUE (default), set \`diag(A)=0\` before forming the
  Laplacian.

## Value

A sparse symmetric matrix (CsparseMatrix).
