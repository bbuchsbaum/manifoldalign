# Build Gram for same-class Laplacian using label factors.

Z^T L_s Z = Z^T D_s Z - (Z^T C)(Z^T C)^T

## Usage

``` r
gram_Ls(Zop, F)
```

## Arguments

- Zop:

  Operator for kernel blocks.

- F:

  Label factor list from \`label_factors()\`.

## Value

Symmetric dense matrix of size r x r.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1))
# F <- label_factors(...)  # Requires label factors
# G <- gram_Ls(Zop, F)
# }
```
