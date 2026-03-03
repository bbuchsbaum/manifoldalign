# Multiset Manifold Alignment (MMA)

Align 3+ domains via spectral embeddings and EM-based probabilistic
registration with optional consensus template. See
mma_align_multiple.hyperdesign for details.

## Usage

``` r
mma_align_multiple(data, ...)

# S3 method for class 'list'
mma_align_multiple(data, ...)
```

## Arguments

- data:

  Input object; hyperdesign or list of matrices (see methods)

- ...:

  Passed to method

## Value

A list containing aligned embeddings, rotation matrices, and convergence
information for all domains.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(60), 30, 2)
X2 <- matrix(rnorm(60), 30, 2)
X3 <- matrix(rnorm(60), 30, 2)
result <- mma_align_multiple(list(X1, X2, X3), ncomp = 2, max_iter = 10)
# }
```
