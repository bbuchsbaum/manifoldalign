# Summary method for pseudolabels objects

Summary method for pseudolabels objects

## Usage

``` r
# S3 method for class 'pseudolabels'
summary(object, ...)
```

## Arguments

- object:

  A pseudolabels object

- ...:

  Additional arguments (ignored)

## Value

The pseudolabels object, invisibly.

## Examples

``` r
# \donttest{
sim_mat <- Matrix::sparseMatrix(
  i = c(1,2,1,3,2,3), j = c(2,1,3,1,3,2),
  x = c(0.9,0.9,0.8,0.8,0.7,0.7), dims = c(5,5)
)
pl <- assign_pseudolabels(sim_mat, min_clusters = 1, min_cluster_size = 1)
summary(pl)
#> Pseudolabeling Summary
#> ======================
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
#> 
#> Detailed Cluster Information:
#>   cluster_id representative size avg_similarity observed_similarity
#> 1          1              1    2            0.9                 0.9
#> 2          2              3    1             NA                  NA
#> 3          3              4    1             NA                  NA
#> 4          4              5    1             NA                  NA
#>   edge_density
#> 1            1
#> 2           NA
#> 3           NA
#> 4           NA
# }
```
