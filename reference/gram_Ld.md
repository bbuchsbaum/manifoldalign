# Build Gram for between-class Laplacian using label factors.

Z^T L_d Z = Z^T D_d Z + (Z^T C)(Z^T C)^T - (Z^T ell)(Z^T ell)^T

## Usage

``` r
gram_Ld(Zop, F)
```

## Arguments

- Zop:

  Operator for kernel blocks.

- F:

  Label factor list.

## Value

Symmetric dense matrix.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1))
# F <- label_factors(...)  # Requires label factors
# G <- gram_Ld(Zop, F)
# }
```
