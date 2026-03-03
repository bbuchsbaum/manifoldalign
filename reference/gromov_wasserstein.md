# Gromov-Wasserstein Distance

Generic function for computing Gromov-Wasserstein distance between
domains

## Usage

``` r
gromov_wasserstein(data, ...)

# S3 method for class 'hyperdesign'
gromov_wasserstein(
  data,
  epsilon = 0.1,
  metric = c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
  init = "uniform",
  max_iter = 100,
  tol = 1e-09,
  verbose = FALSE,
  marginals = NULL,
  scale = c("per_domain", "global", "none"),
  inner_max_iter = 30,
  inner_tol = 1e-09,
  dist_chunk_size = 1000,
  ...
)

# Default S3 method
gromov_wasserstein(data, ...)
```

## Arguments

- data:

  The input data (typically a hyperdesign object)

- ...:

  Additional arguments passed to methods:

  - `epsilon`: Entropic regularization parameter (default: 0.1)

  - `metric`: Distance metric for within-domain distances (default:
    "euclidean")

  - `init`: Initialization method: "uniform" (default), "random", or
    initial transport plans

  - `max_iter`: Maximum iterations for optimization (default: 100)

  - `tol`: Convergence tolerance (default: 1e-9)

  - `verbose`: Print convergence information (default: FALSE)

  - `marginals`: Optional list of marginal distributions

  - `scale`: How to scale distance matrices (default: "per_domain")

  - `inner_max_iter`: Maximum iterations for inner Sinkhorn loop
    (default: 30)

  - `inner_tol`: Convergence tolerance for inner Sinkhorn iterations
    (default: 1e-9)

  - `dist_chunk_size`: Chunk size for memory-efficient distance
    computation (default: 1000)

- epsilon:

  Entropic regularization parameter. Default: 0.1

- metric:

  Distance metric for within-domain distances. Default: "euclidean"

- init:

  Initialization method: "uniform" or "random". Default: "uniform"

- max_iter:

  Maximum iterations for outer loop. Default: 100

- tol:

  Convergence tolerance for outer loop. Default: 1e-9

- verbose:

  Print convergence information. Default: FALSE

- marginals:

  Optional list of marginal distributions per domain

- scale:

  Scaling strategy: "per_domain", "global", or "none". Default:
  "per_domain"

- inner_max_iter:

  Maximum Sinkhorn iterations. Default: 30

- inner_tol:

  Sinkhorn convergence tolerance. Default: 1e-9

- dist_chunk_size:

  Chunk size for distance computation. Default: 1000

## Value

A \`gromov_wasserstein\` object containing transport plans, pairwise
distances, convergence status, and domain metadata. See
\`gromov_wasserstein.hyperdesign\` for full details.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
library(multidesign)
md1 <- multidesign(X1, data.frame(id = 1:10))
md2 <- multidesign(X2, data.frame(id = 1:10))
hd <- hyperdesign(list(d1 = md1, d2 = md2))
result <- gromov_wasserstein(hd, epsilon = 0.5, max_iter = 10)
# }
```
