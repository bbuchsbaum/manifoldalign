# Fit GRASP on a pair of domains

GRASP-specific method for pairwise alignment using spectral bases and
descriptor alignment with an orthogonal rotation matrix.

## Usage

``` r
# S3 method for class 'grasp_aligner'
fit_pair(
  algo,
  X_i,
  X_j,
  links = NULL,
  ncomp = 30,
  q_descriptors = 100,
  sigma = 0.73,
  lambda = 0.1,
  use_laplacian = TRUE,
  solver = c("linear", "hungarian", "auction"),
  ...
)
```

## Arguments

- algo:

  A grasp_aligner object

- X_i:

  First domain data matrix (samples x features)

- X_j:

  Second domain data matrix (samples x features)

- links:

  Optional correspondence links (unused by GRASP)

- ncomp:

  Number of spectral components to compute

- q_descriptors:

  Number of descriptor vectors for alignment

- sigma:

  Bandwidth parameter for descriptor weighting

- lambda:

  Regularization parameter for rotation alignment

- use_laplacian:

  Logical; if TRUE use graph Laplacian, else adjacency

- solver:

  Assignment solver: "linear", "hungarian", or "auction"

- ...:

  Additional arguments (unused)

## Value

An object of class grasp_pair_fit containing rotation (kxk matrix),
mapping_matrix (assignment matrix), assignment vector, loss (numeric),
and k (latent dimension).

## Examples

``` r
# \donttest{
set.seed(123)
X1 <- matrix(rnorm(100), 50, 2)
X2 <- matrix(rnorm(100), 50, 2)
algo <- grasp_aligner()
fit <- fit_pair(algo, X1, X2, ncomp = 5, q_descriptors = 10)
class(fit)
#> [1] "grasp_pair_fit"
# }
```
