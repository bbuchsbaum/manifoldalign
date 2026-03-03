# Apply Z^T to arbitrary multi-vectors without forming Z.

Apply Z^T to arbitrary multi-vectors without forming Z.

## Usage

``` r
Zt_times(Zop, S)
```

## Arguments

- Zop:

  Operator for Z.

- S:

  Dense or sparse matrix with n rows.

## Value

Dense matrix with r rows.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1))
S <- matrix(rnorm(30), 10, 3)
result <- Zt_times(Zop, S)
# }
```
