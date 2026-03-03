# Extract Eigenvalues from KEMA Solver

This function extracts eigenvalues from the KEMA generalized eigenvalue
problem to enable validation against paper specifications.

## Usage

``` r
extract_kema_eigenvalues(
  strata,
  labels,
  kernel = kernlab::rbfdot(sigma = 0.1),
  knn = 5,
  u = 0.5,
  lambda = 0,
  ncomp = 3,
  solver = "exact"
)
```

## Arguments

- strata:

  List of data strata

- labels:

  Vector of class labels

- kernel:

  Kernel function

- knn:

  Number of nearest neighbors

- u:

  Trade-off parameter

- lambda:

  Regularization parameter

- ncomp:

  Number of components

- solver:

  Solver method

## Value

List containing eigenvalues and related validation metrics
