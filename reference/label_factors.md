# Low-rank class graph helpers and matrix-free kernel operators for scalable KEMA

These helpers avoid materializing dense class similarity matrices and
block diagonal kernel matrices when solving the KEMA generalized
eigenproblem. They expose compact representations that can be combined
to build reduced Grams in O(n + rC) memory, where r is the aggregate
kernel rank and C the number of classes.

## Usage

``` r
label_factors(labels)
```

## Arguments

- labels:

  Character/factor vector with NA for unlabeled samples

## Value

List with sparse label factor matrix \`C\`, indicator vector \`ell\`,
and diagonal entries for the same/different-class Laplacians.
