# Build matrix-free operator for block-diagonal kernel matrix constructed from pre-computed kernel blocks.

Build matrix-free operator for block-diagonal kernel matrix constructed
from pre-computed kernel blocks.

## Usage

``` r
make_Zop_from_Ks(Ks)
```

## Arguments

- Ks:

  List of \`Matrix\` objects (each n_i x r_i) representing kernel
  evaluations for a block.

## Value

List with dimensions and forward/transpose operators.

## Examples

``` r
# \donttest{
library(Matrix)
K1 <- Matrix(rnorm(20), 10, 2)
K2 <- Matrix(rnorm(20), 10, 2)
Zop <- make_Zop_from_Ks(list(K1, K2))
str(Zop)
#> List of 5
#>  $ n         : int 20
#>  $ r         : int 4
#>  $ mv        :function (x)  
#>  $ tmv       :function (y)  
#>  $ block_cols: int [1:2] 2 2
# }
```
