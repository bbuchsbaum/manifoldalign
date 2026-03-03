# Extract a relative transform from a pairwise fit

Extract a relative transform from a pairwise fit

## Usage

``` r
# S3 method for class 'grasp_pair_fit'
relative_transform(fit, from = c("i", "j"), to = c("j", "i"), ...)

# S3 method for class 'ot_procrustes_pair_fit'
relative_transform(fit, from = c("i", "j"), to = c("j", "i"), ...)

# S3 method for class 'parrot_pair_fit'
relative_transform(fit, from = c("i", "j"), to = c("j", "i"), ...)

# S3 method for class 'token_ot_graph_pair_fit'
relative_transform(fit, from = c("i", "j"), to = c("j", "i"), ...)

relative_transform(fit, from = c("i", "j"), to = c("j", "i"), ...)
```

## Arguments

- fit:

  A pairwise fit object returned by fit_pair

- from:

  Source domain identifier ("i" or "j")

- to:

  Target domain identifier ("j" or "i")

- ...:

  Additional arguments passed to methods

## Value

An align_transform object (see new_align_transform) containing the
transformation operator, group type, and source/target indices.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 25, 2)
X2 <- matrix(rnorm(50), 25, 2)
algo <- grasp_aligner()
fit <- fit_pair(algo, X1, X2)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colSums': non-conformable matrix dimensions in .Arith.Csparse(e1, e2, .Generic, class. = "dgCMatrix")
transform <- relative_transform(fit, from = "i", to = "j")
#> Error in relative_transform.default(fit, from = "i", to = "j"): relative_transform() is not implemented for object of class function
# }
```
