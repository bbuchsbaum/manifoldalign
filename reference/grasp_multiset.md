# Multi-Graph GRASP: Graph Alignment by Spectral Corresponding Functions

Performs GRASP alignment on three or more graph domains simultaneously.
Extends the pairwise GRASP algorithm to handle multiple graphs through
joint diagonalization and shared latent basis alignment.

## Usage

``` r
grasp_multiset(data, ...)

grasp_multiset(data, ...)

# S3 method for class 'hyperdesign'
grasp_multiset(
  data,
  preproc = multivarious::center(),
  ncomp = 30,
  q_descriptors = 100,
  sigma = 0.73,
  lambda = 0.1,
  use_laplacian = TRUE,
  anchor = 1L,
  solver = c("linear", "auction"),
  assignment_alpha = 0.9,
  max_iter = 100,
  tol = 1e-06,
  ...
)

# S3 method for class 'list'
grasp_multiset(data, ...)

# Default S3 method
grasp_multiset(data, ...)
```

## Arguments

- data:

  A list containing three or more matrices (nodes x features)

- ...:

  Additional arguments passed to hyperdesign method

- preproc:

  Preprocessing step applied via multivarious::init_transform

- ncomp:

  Number of spectral components (default: 30)

- q_descriptors:

  Number of descriptors (default: 100)

- sigma:

  Diffusion parameter (default: 0.73)

- lambda:

  Regularization parameter (default: 0.1)

- use_laplacian:

  Whether to use Laplacian normalization (default: TRUE)

- anchor:

  Index of reference graph (1-based) or "mean" (default: 1)

- solver:

  Assignment solver: "auction" or "linear" (default: "auction")

- assignment_alpha:

  Weight for cosine vs. Euclidean assignment cost (0-1)

- max_iter:

  Maximum iterations for joint optimization (default: 100)

- tol:

  Convergence tolerance (default: 1e-6)

## Value

A `grasp_multiset` object containing alignment results for all graphs

A grasp_multiset object

## See also

[`grasp`](https://bbuchsbaum.github.io/manifoldalign/reference/grasp.md)
for pairwise alignment, `grasp_multiset.hyperdesign` for detailed
documentation

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
X3 <- matrix(rnorm(30), 10, 3)
data_list <- list(X1, X2, X3)
result <- grasp_multiset(data_list, ncomp = 5)
#> Warning: Joint alignment reached maximum iterations without convergence
str(result$embeddings)
#> List of 3
#>  $ : num [1:10, 1:5] -0.167 -0.271 -0.182 0.356 -0.057 ...
#>  $ : num [1:10, 1:5] 1.8 -1.84 -1.55 -0.44 1.73 ...
#>  $ : num [1:10, 1:5] -1.3831 -1.3723 -0.1632 1.4176 0.0708 ...
# }
```
