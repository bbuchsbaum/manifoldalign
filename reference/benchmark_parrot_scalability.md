# Benchmark PARROT Scalability

Tests PARROT performance across different network sizes

## Usage

``` r
benchmark_parrot_scalability(
  sizes = c(50, 100, 200, 500, 1000),
  n_reps = 3,
  sparse_graph = TRUE
)
```

## Arguments

- sizes:

  Vector of network sizes to test

- n_reps:

  Number of repetitions per size

- sparse_graph:

  Whether to use sparse graph structure

## Value

Data frame with timing results

## Examples

``` r
# \donttest{
results <- benchmark_parrot_scalability(sizes = c(20, 50), n_reps = 2)
#> Benchmarking PARROT scalability...
#> 
#> Testing n=20 nodes:
#>   Rep 1/2... 0.35s
#>   Rep 2/2... 0.31s
#> 
#> Testing n=50 nodes:
#>   Rep 1/2... 0.31s
#>   Rep 2/2... 0.31s
#> 
#> 
#> SCALABILITY SUMMARY
#> ===================
#>   n_nodes total_time.mean total_time.sd
#> 1      20     0.332709551   0.029944493
#> 2      50     0.309624910   0.002591525
#> 
#> Estimated complexity: O(n^-0.08)
# }
```
