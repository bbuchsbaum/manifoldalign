# Plot cycle-consistency residuals as a heatmap

Aggregates triangle (i,j,k) residuals into a per-pair matrix and draws a
heatmap. Useful for spotting inconsistent domains or weak edges after
synchronization.

## Usage

``` r
plot_cycle_consistency(
  diagnostics,
  N = NULL,
  aggregate = c("min", "mean", "max"),
  na_color = "#f0f0f0"
)
```

## Arguments

- diagnostics:

  List returned in \`align_many(...)\[\["diagnostics"\]\]\`

- N:

  Optional number of domains; inferred from diagnostics when NULL

- aggregate:

  How to aggregate multiple j-paths per (i,k): "min", "mean", or "max"

- na_color:

  Fill color for missing cells (no triangles)

## Value

A ggplot object

## Examples

``` r
# \donttest{
# Requires align_many result with diagnostics
# X1 <- matrix(rnorm(100), 50, 2)
# X2 <- matrix(rnorm(100), 50, 2)
# X3 <- matrix(rnorm(100), 50, 2)
# algo <- grasp_aligner()
# result <- align_many(list(X1, X2, X3), algo)
# plot_cycle_consistency(result$diagnostics)
# }
```
