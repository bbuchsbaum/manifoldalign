# Multiple Graph CONE-Align: Consensus Optimization for Multi-Domain Node Embedding Alignment

Performs CONE-Align on three or more graph domains simultaneously.
Extends the pairwise CONE-Align algorithm to handle multiple graphs
through iterative alignment to a common reference frame.

## Usage

``` r
cone_align_multiple(data, ...)

# S3 method for class 'hyperdesign'
cone_align_multiple(
  data,
  ref_idx = 1L,
  preproc = center(),
  ncomp = 10,
  sigma = 0.73,
  lambda = 0.1,
  use_laplacian = TRUE,
  solver = c("linear", "auction"),
  max_iter = 30,
  tol = 0.01,
  knn = NULL,
  ...
)

# S3 method for class 'list'
cone_align_multiple(data, ...)

# Default S3 method
cone_align_multiple(data, ...)
```

## Arguments

- data:

  A list containing three or more matrices (nodes x features)

- ...:

  Additional arguments passed to hyperdesign method

- ref_idx:

  Which domain acts as initial reference (default: 1)

- preproc:

  Preprocessing function to apply to the data (default: center())

- ncomp:

  Number of embedding dimensions (default: 10)

- sigma:

  Diffusion parameter for graph construction (default: 0.73)

- lambda:

  Regularization parameter for numerical stability (default: 0.1)

- use_laplacian:

  Whether to use Laplacian normalization (default: TRUE)

- solver:

  Assignment algorithm: "linear" for exact assignment (default),
  "auction" for large-scale approximation

- max_iter:

  Maximum number of iterations (default: 30)

- tol:

  Convergence tolerance for assignment changes (default: 0.01)

- knn:

  Number of nearest neighbors for graph construction (default: adaptive)

## Value

A `multiblock_biprojector` object containing alignment results for all
graphs

## See also

[`cone_align`](https://bbuchsbaum.github.io/manifoldalign/reference/cone_align.md)
for pairwise alignment, `cone_align_multiple.hyperdesign` for detailed
documentation

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(60), 30, 2)
X2 <- matrix(rnorm(60), 30, 2)
X3 <- matrix(rnorm(60), 30, 2)
result <- cone_align_multiple(list(X1, X2, X3), ncomp = 2)
#> Converged after 3 iterations.
#> Converged after 2 iterations.
# }
```
