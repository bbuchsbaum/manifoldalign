# Fit PARROT on a pair of domains

PARROT-specific method for pairwise alignment using optimal transport
with random walk regularization. Returns a soft assignment/transport
plan.

## Usage

``` r
# S3 method for class 'parrot_aligner'
fit_pair(
  algo,
  X_i,
  X_j,
  links = NULL,
  ncomp = NULL,
  sigma = 0.15,
  lambda = 0.1,
  lambda_e = NULL,
  lambda_n = NULL,
  lambda_p = NULL,
  tau = 0.05,
  alpha = 0.2,
  gamma = 0.1,
  solver = c("sinkhorn"),
  max_iter = 100,
  tol = 1e-06,
  use_cpp = FALSE,
  anchors_policy = c("warn", "off", "error"),
  unsup_strategy = c("mnn_bootstrap", "attribute_only", "self_training"),
  conf_threshold = 0.9,
  max_rounds = 1,
  seed_strength = 0.05,
  laplacian_reg = NULL,
  ...
)
```

## Arguments

- algo:

  A parrot_aligner object

- X_i:

  First domain data matrix (samples x features)

- X_j:

  Second domain data matrix (samples x features)

- links:

  Optional correspondence links; list(vec1, vec2) with label vectors or
  NA for unmatched samples

- ncomp:

  Number of spectral components (currently unused)

- sigma:

  RWR restart probability

- lambda:

  Overall regularization weight

- lambda_e:

  Edge/attribute cost weight (defaults to lambda \* 0.5)

- lambda_n:

  Network/Laplacian regularization (defaults to lambda \* 0.5)

- lambda_p:

  Anchor prior weight (defaults to lambda)

- tau:

  Entropic regularization parameter for Sinkhorn

- alpha:

  Mixing parameter for cost matrix

- gamma:

  Network structure penalty weight

- solver:

  Transport solver method ("sinkhorn")

- max_iter:

  Maximum iterations for RWR and Sinkhorn

- tol:

  Convergence tolerance

- use_cpp:

  Logical; use C++ implementation if available

- anchors_policy:

  Policy when no anchors provided: "warn", "off", "error"

- unsup_strategy:

  Unsupervised seeding strategy: "mnn_bootstrap", "attribute_only", or
  "self_training"

- conf_threshold:

  Confidence threshold for self-training

- max_rounds:

  Maximum rounds for iterative strategies

- seed_strength:

  Anchor prior strength when using MNN seeds

- laplacian_reg:

  Alternative to lambda_n for Laplacian regularization

- ...:

  Additional arguments (unused)

## Value

An object of class parrot_pair_fit containing transport (n1 x n2
matrix), objective (final cost), n1 (samples in domain i), and n2
(samples in domain j).

## Examples

``` r
# \donttest{
set.seed(456)
X1 <- matrix(rnorm(80), 40, 2)
X2 <- matrix(rnorm(80), 40, 2)
algo <- parrot_aligner()
fit <- fit_pair(algo, X1, X2, sigma = 0.2, tau = 0.1)
#> Warning: PARROT: no anchors provided; using unsupervised seeding (mnn_bootstrap)
class(fit)
#> [1] "parrot_pair_fit"
# }
```
