# High-similarity pseudolabels from multi-domain data

Generates pseudolabels based on high cosine similarity across different
domains. Specialized version for multi-domain case.

## Usage

``` r
high_sim_pseudolabels(
  strata,
  k = 5,
  cos_thresh = 0.97,
  min_size = 2,
  ann_trees = 50,
  verbose = FALSE
)
```

## Arguments

- strata:

  List produced by multidesign::hyperdesign. Each element should have
  \$x (samples x features matrix) and \$design (metadata data.frame)

- k:

  Number of cross-domain neighbors to examine (default: 5)

- cos_thresh:

  Minimum cosine similarity threshold (default: 0.97)

- min_size:

  Minimum cluster size to retain (default: 2)

- ann_trees:

  Number of trees for RcppAnnoy index (default: 50)

- verbose:

  Whether to print progress information (default: FALSE)

## Value

Factor vector of pseudolabels with length equal to total number of
samples across all domains. NA indicates samples with no confident
match.

## Details

This function implements the approach described in the pseudolabel
documentation:

1.  Pool and L2-normalize all samples across domains

2.  Build approximate nearest neighbor index using angular distance

3.  Find high-similarity cross-domain pairs

4.  Use connected components to form clusters

5.  Filter clusters by minimum size

## Examples

``` r
# \donttest{
if (requireNamespace("RcppAnnoy", quietly = TRUE) &&
    requireNamespace("igraph", quietly = TRUE)) {
  set.seed(1)
  strata <- list(
    list(x = matrix(rnorm(30), 10, 3), design = data.frame(row = 1:10)),
    list(x = matrix(rnorm(30), 10, 3), design = data.frame(row = 1:10))
  )
  plabs <- high_sim_pseudolabels(strata, k = 3, cos_thresh = 0.95, ann_trees = 10)
  table(plabs, useNA = "always")
}
#> plabs
#> anchor_0005        <NA> 
#>           2          18 
# }
```
