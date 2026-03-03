# Apply a transform to data

Apply a transform to data

## Usage

``` r
apply_transform(transform, X, ...)
```

## Arguments

- transform:

  An alignment transform object (typically from relative_transform or
  new_align_transform)

- X:

  A data matrix to transform (samples x features)

- ...:

  Additional arguments passed to methods

## Value

A transformed data matrix. For orthogonal/linear groups, returns X For
permutation group, returns op

## Examples

``` r
# Create a simple orthogonal transform
R <- diag(3)
tr <- new_align_transform("O", R, from = 1, to = 2)
X <- matrix(rnorm(12), 4, 3)
X_transformed <- apply_transform(tr, X)
```
