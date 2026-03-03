# Fast path for diagonal A in Z^T A Z.

Fast path for diagonal A in Z^T A Z.

## Usage

``` r
gram_Z_diag_Z(Zop, diag_entries, block_size = 128L)
```

## Arguments

- Zop:

  Operator for Z.

- diag_entries:

  Numeric vector of length n (diagonal of A).

- block_size:

  Batch size for columns.

## Value

Symmetric dense matrix.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1))
diag_vals <- rep(1, 10)
G <- gram_Z_diag_Z(Zop, diag_vals)
# }
```
