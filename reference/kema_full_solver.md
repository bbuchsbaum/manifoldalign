# Full KEMA Solver (sample_frac == 1)

Implements the full KEMA algorithm for complete kernel matrices. This
solver handles the n x n generalized eigenvalue problem.

## Usage

``` r
kema_full_solver(
  strata,
  Z,
  Ks,
  Lap,
  kernel_indices,
  solver,
  ncomp,
  u,
  dweight,
  rweight,
  lambda
)
```

## Arguments

- strata:

  List of data strata

- Z:

  Block diagonal kernel matrix (n x n)

- Ks:

  List of individual kernel matrices

- Lap:

  List of Laplacian matrices

- kernel_indices:

  Block indices for kernels

- solver:

  Solver method: "regression" for fast approximation or "exact" for
  precise solution

- ncomp:

  Number of components to extract

- u:

  Trade-off parameter between manifold and class terms

- dweight:

  Weight for dissimilarity terms

- rweight:

  Weight for repulsion terms

- lambda:

  Regularization parameter

## Value

List with coef, scores, and vectors
