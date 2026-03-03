# Print Method for FPGW Objects

Print Method for FPGW Objects

## Usage

``` r
# S3 method for class 'fpgw'
print(x, ...)
```

## Arguments

- x:

  An fpgw object

- ...:

  Additional arguments (currently unused)

## Value

The fpgw object x, invisibly.

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
library(multidesign)
md1 <- multidesign(X1, data.frame(id = 1:10))
md2 <- multidesign(X2, data.frame(id = 1:10))
hd <- hyperdesign(list(d1 = md1, d2 = md2))
res <- fpgw(hd)
print(res)
#> Fused-Partial Gromov-Wasserstein
#> ================================
#> Number of domains: 2 
#> Domain names: d1, d2 
#> Feature weight (omega1): 0.001 
#> Mode: Classical Fused GW
#> 
#> Pairwise distances:
#>        [,1]   [,2]
#> [1,] 0.0000 0.2806
#> [2,] 0.2806 0.0000
# }
```
