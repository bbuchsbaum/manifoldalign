# Spectral MNN aligner (native multi-view)

Spectral MNN aligner (native multi-view)

## Usage

``` r
spectral_mnn_aligner()
```

## Value

An object of class \`c("spectral_mnn_aligner", "aligner")\` representing
a spectral MNN alignment algorithm descriptor.

## Examples

``` r
algo <- spectral_mnn_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "O"
#> 
#> $supports_multi
#> [1] TRUE
#> 
```
