# Extract coefficients from FPGW

Extract coefficients from FPGW

## Usage

``` r
# S3 method for class 'fpgw'
coef(object, type = c("transport", "parameters"), ...)
```

## Arguments

- object:

  An fpgw object

- type:

  Type of coefficients to extract

- ...:

  Additional arguments

## Value

Coefficients (transport plans or parameters)

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(30), 10, 3)
X2 <- matrix(rnorm(30), 10, 3)
library(multidesign)
#> 
#> Attaching package: ‘multidesign’
#> The following object is masked from ‘package:manifoldalign’:
#> 
#>     block_indices
md1 <- multidesign(X1, data.frame(id = 1:10))
md2 <- multidesign(X2, data.frame(id = 1:10))
hd <- hyperdesign(list(d1 = md1, d2 = md2))
res <- fpgw(hd)
transport_plans <- coef(res, type = "transport")
params <- coef(res, type = "parameters")
# }
```
