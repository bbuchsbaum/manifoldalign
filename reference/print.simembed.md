# Print Method for Similarity Embedding

Print Method for Similarity Embedding

## Usage

``` r
# S3 method for class 'simembed'
print(x, ...)
```

## Arguments

- x:

  A simembed object

- ...:

  Additional arguments

## Value

The simembed object x, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(100), 20, 5)
model <- linear_sim_embed(X, ncomp = 2)
print(model)
#> Linear Similarity Embedding
#> ===========================
#> Embedding dimension: 2
#> Original dimension: 5
#> Number of samples: 20
#> Optimizer: ADAM (R)
#> Final sigma_P: 3.162278
#> Alpha_p: 0.100
#> Convergence: SUCCESS
#> Iterations: 2
# }
```
