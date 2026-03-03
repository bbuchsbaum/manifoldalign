# Multi-view aligner adapters (native delegates)

These adapters declare native multi-domain support and delegate to
existing multi-view functions in this package. KEMA aligner (native
multi-view)

## Usage

``` r
kema_aligner()
```

## Value

An object of class \`c("kema_aligner", "aligner")\` representing a KEMA
alignment algorithm descriptor.

## Examples

``` r
algo <- kema_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "kernel"
#> 
#> $supports_multi
#> [1] TRUE
#> 
```
