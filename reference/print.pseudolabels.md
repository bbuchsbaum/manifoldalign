# Print method for pseudolabels objects

Print method for pseudolabels objects

## Usage

``` r
# S3 method for class 'pseudolabels'
print(x, ...)
```

## Arguments

- x:

  A pseudolabels object

- ...:

  Additional arguments (ignored)

## Value

The pseudolabels object x, invisibly.

## Examples

``` r
# \donttest{
sim_mat <- Matrix::sparseMatrix(
  i = c(1,2,1,3,2,3), j = c(2,1,3,1,3,2),
  x = c(0.9,0.9,0.8,0.8,0.7,0.7), dims = c(5,5)
)
pl <- assign_pseudolabels(sim_mat, min_clusters = 1, min_cluster_size = 1)
print(pl)
#> Pseudolabeling Results
#> ======================
#> Algorithm:         assign_pseudolabels 
#> Diversity method: greedy (weight: 0.3)
#> Samples:           5 
#> Assigned:         5 (100%)
#> Clusters found:    4 
#> Threshold used:    0.9 
#> 
#> Cluster Summary:
#>   Size range:      1 - 2 
#>   Mean size:       1.2 
#>   Observed similarity: 0.9-0.9 (mean: 0.9)
# }
```
