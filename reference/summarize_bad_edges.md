# Summarize the worst edges by residual

Returns a data.frame of the top-K edges with largest residuals; includes
edge weights when available.

## Usage

``` r
summarize_bad_edges(diagnostics, top = 10)
```

## Arguments

- diagnostics:

  List returned by align_many(...)\[\["diagnostics"\]\]

- top:

  Number of worst edges to return

## Value

data.frame with columns i, j, residual, weight

## Examples

``` r
# \donttest{
# Requires align_many result with diagnostics
# X1 <- matrix(rnorm(100), 50, 2)
# X2 <- matrix(rnorm(100), 50, 2)
# X3 <- matrix(rnorm(100), 50, 2)
# algo <- grasp_aligner()
# result <- align_many(list(X1, X2, X3), algo)
# bad_edges <- summarize_bad_edges(result$diagnostics, top = 5)
# }
```
