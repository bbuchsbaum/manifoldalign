# Fit a pairwise alignment model

Fit a pairwise alignment model

## Usage

``` r
fit_pair(algo, X_i, X_j, links = NULL, ...)
```

## Arguments

- algo:

  An aligner object created by an aligner constructor

- X_i:

  First domain data matrix (samples x features)

- X_j:

  Second domain data matrix (samples x features)

- links:

  Optional correspondence links between domains

- ...:

  Additional arguments passed to methods

## Value

A pairwise fit object whose class depends on the method (e.g.,
grasp_pair_fit, parrot_pair_fit). Contains alignment results including
transform parameters and loss.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 25, 2)
X2 <- matrix(rnorm(50), 25, 2)
algo <- grasp_aligner()
fit <- fit_pair(algo, X1, X2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': non-conformable matrix dimensions in .Arith.Csparse(e1, e2, .Generic, class. = "dgCMatrix")
# }
```
