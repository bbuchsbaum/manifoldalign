# Core Multi-Graph CONE-Align Fitting

Internal function that implements the iterative multi-graph alignment
algorithm. Uses the helper functions from the pairwise implementation.

## Usage

``` r
cone_align_multiple_fit(
  strata,
  proc,
  ref_idx,
  ncomp,
  sigma,
  lambda,
  use_laplacian,
  solver,
  max_iter,
  tol,
  knn,
  feature_blocks
)
```

## Arguments

- strata:

  List of preprocessed data domains

- proc:

  Preprocessing object

- ref_idx:

  Reference domain index

- ncomp:

  Number of components

- sigma:

  Diffusion parameter

- lambda:

  Regularization parameter

- use_laplacian:

  Whether to use Laplacian normalization

- solver:

  Assignment solver

- max_iter:

  Maximum iterations

- tol:

  Convergence tolerance

- knn:

  Number of nearest neighbors

- feature_blocks:

  Feature block index structure

## Value

multiblock_biprojector object with multi-graph alignment results
