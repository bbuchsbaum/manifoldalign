# Profile PARROT Memory Usage

Estimates memory usage for different network sizes

## Usage

``` r
profile_parrot_memory(n_nodes = 100)
```

## Arguments

- n_nodes:

  Network size to profile

## Value

List with memory usage statistics

## Examples

``` r
# \donttest{
mem_profile <- profile_parrot_memory(n_nodes = 50)
#> Warning: pryr package not available for memory profiling
# }
```
