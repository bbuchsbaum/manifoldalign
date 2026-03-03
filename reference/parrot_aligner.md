# PARROT adapter for align_many()

Provides a pairwise interface for PARROT using internal building blocks.
Returns a soft assignment/transport map suitable for permutation sync.
Construct a PARROT aligner descriptor

## Usage

``` r
parrot_aligner()
```

## Value

an object of class c("parrot_aligner", "aligner")

## Examples

``` r
algo <- parrot_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "perm"
#> 
#> $supports_multi
#> [1] FALSE
#> 
```
