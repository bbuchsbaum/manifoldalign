# Convert pseudolabels to data.frame

Convert pseudolabels to data.frame

## Usage

``` r
# S3 method for class 'pseudolabels'
as.data.frame(x, ...)
```

## Arguments

- x:

  A pseudolabels object

- ...:

  Additional arguments (ignored)

## Value

A data.frame with columns for sample index, label, and representative
status.

## Examples

``` r
# \donttest{
sim_mat <- Matrix::sparseMatrix(
  i = c(1,2,1,3,2,3), j = c(2,1,3,1,3,2),
  x = c(0.9,0.9,0.8,0.8,0.7,0.7), dims = c(5,5)
)
pl <- assign_pseudolabels(sim_mat, min_clusters = 1, min_cluster_size = 1)
df <- as.data.frame(pl)
head(df)
#>   sample_id      label is_representative
#> 1         1 anchor_001              TRUE
#> 2         2 anchor_001             FALSE
#> 3         3 anchor_002              TRUE
#> 4         4 anchor_003              TRUE
#> 5         5 anchor_004              TRUE
# }
```
