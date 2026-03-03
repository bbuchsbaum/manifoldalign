# Joint Basis Alignment for Multiple Graphs

Performs joint diagonalization to find rotations that simultaneously
align all graph spectral bases.

## Usage

``` r
align_bases_multiset(bases, descriptors, lambda, max_iter = 100, tol = 1e-06)
```

## Arguments

- bases:

  List of spectral bases from compute_grasp_basis

- descriptors:

  List of descriptor matrices

- lambda:

  Regularization parameter

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

## Value

List with rotations and convergence info
