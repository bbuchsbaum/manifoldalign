# Compute Z^T Z without forming Z explicitly.

Compute Z^T Z without forming Z explicitly.

## Usage

``` r
gram_ZZ(Zop, block_size = 128L)
```

## Arguments

- Zop:

  Operator as returned by \`make_Zop_from_Ks\` or similar.

- block_size:

  Number of identity columns processed per batch.

## Value

Symmetric dense matrix of size r x r.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
K2 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1, K2))
G <- gram_ZZ(Zop)
# }
```
