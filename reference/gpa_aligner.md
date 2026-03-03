# Generalized Procrustes aligner (native multi-view)

Generalized Procrustes aligner (native multi-view)

## Usage

``` r
gpa_aligner()
```

## Value

An object of class \`c("gpa_aligner", "aligner")\` representing a
Generalized Procrustes alignment algorithm descriptor.

## Examples

``` r
algo <- gpa_aligner()
aligner_capabilities(algo)
#> $group
#> [1] "O"
#> 
#> $supports_multi
#> [1] TRUE
#> 
```
