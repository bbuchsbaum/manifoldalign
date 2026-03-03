# SSMA aligner (native multi-view)

SSMA aligner (native multi-view)

## Usage

``` r
ssma_aligner()
```

## Value

An object of class \`c("ssma_aligner", "aligner")\` representing a
semi-supervised manifold alignment (SSMA) algorithm descriptor.

## Examples

``` r
algo <- ssma_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "GL"
#> 
#> $supports_multi
#> [1] TRUE
#> 
```
