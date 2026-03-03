# Block C: Base Alignment - Corrected

Block C: Base Alignment - Corrected

## Usage

``` r
align_grasp_bases(
  basis1,
  basis2,
  desc1,
  desc2,
  lambda = 0.1,
  max_iter = 200,
  tol = 1e-08
)
```

## Arguments

- basis1:

  First domain's spectral basis

- basis2:

  Second domain's spectral basis

- desc1:

  First domain's descriptors

- desc2:

  Second domain's descriptors

- lambda:

  Regularization parameter

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

## Value

List with rotation matrix and convergence info
