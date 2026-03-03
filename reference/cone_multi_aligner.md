# CONE-Align multiple aligner (native multi-view)

CONE-Align multiple aligner (native multi-view)

## Usage

``` r
cone_multi_aligner()
```

## Value

An object of class \`c("cone_multi_aligner", "aligner")\` representing a
CONE-Align multiple alignment algorithm descriptor.

## Examples

``` r
algo <- cone_multi_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "O"
#> 
#> $supports_multi
#> [1] TRUE
#> 
```
