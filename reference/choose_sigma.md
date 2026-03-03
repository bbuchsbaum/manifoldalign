# Choose optimal sigma for RBF kernel

Estimates a reasonable sigma parameter for RBF kernels using the median
pairwise distance heuristic, which often works well in practice.

## Usage

``` r
choose_sigma(X, sample_size = 1000)
```

## Arguments

- X:

  Data matrix (samples x features)

- sample_size:

  Maximum number of pairs to sample for distance computation (default:
  1000 for efficiency)

## Value

Suggested sigma value

## Examples

``` r
X <- matrix(rnorm(100), 20, 5)
sigma <- choose_sigma(X)
sigma
#> [1] 1.845579
```
