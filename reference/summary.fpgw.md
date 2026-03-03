# Summary method for FPGW

Summary method for FPGW

## Usage

``` r
# S3 method for class 'fpgw'
summary(object, ...)
```

## Arguments

- object:

  An fpgw object

- ...:

  Additional arguments

## Value

The fpgw object, invisibly.

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
summary(res)
#> Fused-Partial Gromov-Wasserstein Summary
#> ========================================
#> 
#> Number of domains: 2 
#> Domain names: d1, d2 
#> Samples per domain: 10, 10 
#> 
#> Parameters:
#>   Feature weight (omega1): 0.001 
#>   Structure weight (omega2): 0.999 
#>   Variant: Classical FGW
#> 
#> Distance Matrix:
#>        [,1]   [,2]
#> [1,] 0.0000 0.2806
#> [2,] 0.2806 0.0000
#> 
#> Transport Plan Statistics:
#>   Plan 1: mass = 1.0000, sparsity = 90.00%
# }
```
