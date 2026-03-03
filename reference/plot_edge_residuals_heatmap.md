# Plot edge residuals as a heatmap

Uses diagnostics\$edge_residuals (i,j,residual) to build an N x N matrix
of symmetric residuals (min over directions when both i-\>j and j-\>i
are given), then draws a heatmap.

## Usage

``` r
plot_edge_residuals_heatmap(diagnostics, N = NULL, na_color = "#f0f0f0")
```

## Arguments

- diagnostics:

  List returned by align_many(...)\[\["diagnostics"\]\]

- N:

  Optional number of domains; inferred when NULL

- na_color:

  Fill color for missing cells

## Value

ggplot object

## Examples

``` r
# \donttest{
# Requires align_many result with diagnostics
# X1 <- matrix(rnorm(100), 50, 2)
# X2 <- matrix(rnorm(100), 50, 2)
# X3 <- matrix(rnorm(100), 50, 2)
# algo <- grasp_aligner()
# result <- align_many(list(X1, X2, X3), algo)
# plot_edge_residuals_heatmap(result$diagnostics)
# }
```
