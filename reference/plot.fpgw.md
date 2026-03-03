# Plot method for FPGW

Plot method for FPGW

## Usage

``` r
# S3 method for class 'fpgw'
plot(x, which = 1, pair = 1, ...)
```

## Arguments

- x:

  An fpgw object

- which:

  Which plot to create (1: distance matrix, 2: transport plan)

- pair:

  For transport plan, which pair to plot (default: 1)

- ...:

  Additional arguments passed to plotting functions

## Value

NULL, invisibly. Called for its side effect of producing a plot.

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
plot(res, which = 1)

plot(res, which = 2, pair = 1)

# }
```
