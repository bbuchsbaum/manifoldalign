# Pairwise loss (for weighting/diagnostics)

Pairwise loss (for weighting/diagnostics)

## Usage

``` r
# S3 method for class 'grasp_pair_fit'
pair_loss(fit, X_i = NULL, X_j = NULL, ...)

# S3 method for class 'ot_procrustes_pair_fit'
pair_loss(fit, X_i = NULL, X_j = NULL, ...)

# S3 method for class 'parrot_pair_fit'
pair_loss(fit, X_i = NULL, X_j = NULL, ...)

# S3 method for class 'token_ot_graph_pair_fit'
pair_loss(fit, X_i = NULL, X_j = NULL, ...)

pair_loss(fit, X_i = NULL, X_j = NULL, ...)
```

## Arguments

- fit:

  A pairwise fit object returned by fit_pair

- X_i:

  Optional first domain data for loss computation

- X_j:

  Optional second domain data for loss computation

- ...:

  Additional arguments passed to methods

## Value

A single numeric value representing the pairwise alignment loss, or
NA_real\_ if not available.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 25, 2)
X2 <- matrix(rnorm(50), 25, 2)
algo <- grasp_aligner()
fit <- fit_pair(algo, X1, X2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': non-conformable matrix dimensions in .Arith.Csparse(e1, e2, .Generic, class. = "dgCMatrix")
loss <- pair_loss(fit)
# }
```
