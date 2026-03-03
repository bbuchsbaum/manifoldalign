# Extract fitted values

Extract fitted values

## Usage

``` r
# S3 method for class 'fpgw'
fitted(object, ...)
```

## Arguments

- object:

  An fpgw object

- ...:

  Additional arguments

## Value

List of transport plans

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
plans <- fitted(res)
# }
```
