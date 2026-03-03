# Compose two transforms (a-\>b and b-\>c)

Compose two transforms (a-\>b and b-\>c)

## Usage

``` r
compose_transform(T_ab, T_bc, ...)
```

## Arguments

- T_ab:

  An alignment transform from domain a to domain b

- T_bc:

  An alignment transform from domain b to domain c

- ...:

  Additional arguments passed to methods

## Value

An align_transform object representing the composed transform from
domain a to domain c.

## Examples

``` r
# Compose two orthogonal transforms
R1 <- diag(3)
R2 <- diag(3)
T_ab <- new_align_transform("O", R1, from = 1, to = 2)
T_bc <- new_align_transform("O", R2, from = 2, to = 3)
T_ac <- compose_transform(T_ab, T_bc)
```
