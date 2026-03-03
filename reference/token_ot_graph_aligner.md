# Token-OT graph aligner adapter for align_many()

Provides a pairwise interface for token_ot_graph_align() using the
unified aligner adapter framework (fit_pair/relative_transform).
Construct a Token-OT Graph aligner descriptor

## Usage

``` r
token_ot_graph_aligner()
```

## Value

an object of class c("token_ot_graph_aligner", "aligner")

## Examples

``` r
algo <- token_ot_graph_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "perm"
#> 
#> $supports_multi
#> [1] FALSE
#> 
```
