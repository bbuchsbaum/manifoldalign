# Generic helper for Z^T A Z using an operator.

Generic helper for Z^T A Z using an operator.

## Usage

``` r
gram_Z_A_Z(Zop, A, block_size = 64L)
```

## Arguments

- Zop:

  Matrix-free operator.

- A:

  Sparse or diagonal matrix compatible with Z columns.

- block_size:

  Number of columns processed per batch.

## Value

Symmetric dense matrix.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1))
A <- Diagonal(10, 1)
G <- gram_Z_A_Z(Zop, A)
# }
```
