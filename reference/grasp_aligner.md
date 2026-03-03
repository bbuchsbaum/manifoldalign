# GRASP adapter for align_many()

Provides the minimal pairwise interface required by align_many() for the
GRASP algorithm. Uses the existing GRASP building blocks to estimate a
relative orthogonal transform between two domains. Construct a GRASP
aligner descriptor

## Usage

``` r
grasp_aligner()
```

## Value

an object of class c("grasp_aligner", "aligner")

## Examples

``` r
algo <- grasp_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "O"
#> 
#> $supports_multi
#> [1] FALSE
#> 
```
