# Latent dimension selected by the method (if applicable)

Latent dimension selected by the method (if applicable)

## Usage

``` r
# S3 method for class 'grasp_pair_fit'
latent_dim(fit, ...)

# S3 method for class 'ot_procrustes_pair_fit'
latent_dim(fit, ...)

latent_dim(fit, ...)
```

## Arguments

- fit:

  A fit object returned by fit_pair or fit_many

- ...:

  Additional arguments passed to methods

## Value

An integer giving the latent dimension selected by the method, or NULL
if not applicable.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 25, 2)
X2 <- matrix(rnorm(50), 25, 2)
algo <- grasp_aligner()
fit <- fit_pair(algo, X1, X2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': non-conformable matrix dimensions in .Arith.Csparse(e1, e2, .Generic, class. = "dgCMatrix")
dim <- latent_dim(fit)
# }
```
