# Transform new data using FPGW alignment

Transform new data using FPGW alignment

## Usage

``` r
# S3 method for class 'fpgw'
transform(`_data`, newdata, source_index = 1, target_index = 2, ...)
```

## Arguments

- \_data:

  An fpgw object

- newdata:

  New data to transform

- source_index:

  Index of source domain

- target_index:

  Index of target domain

- ...:

  Additional arguments

## Value

Transformed samples in the target domain feature space.

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
X_new <- matrix(rnorm(15), 5, 3)
transformed <- transform(res, X_new, source_index = 1, target_index = 2)
#> Error in as.data.frame.default(`_data`): cannot coerce class ‘c("fpgw", "multiblock_biprojector")’ to a data.frame
# }
```
