# Create a similarity function from a label similarity matrix

Builds a function that maps a vector of labels to a sparse similarity
matrix. Unknown/NA labels receive all-zero rows/columns, and the result
is guaranteed symmetric with zero diagonal.

## Usage

``` r
createSimFun(S, na_value = 0)
```

## Arguments

- S:

  A square similarity matrix whose row/column names define the label
  vocabulary. Only labels present in \`rownames(S)\` are assigned
  similarities.

- na_value:

  Deprecated; ignored (kept for backward compatibility). NA and unknown
  labels always map to zero similarity to preserve sparsity.

## Value

A function \`f(labels)\` that returns an \`n x n\` sparse similarity
matrix.

## Details

This is useful as \`simfun\` for \[\`lowrank_align()\`\], where labels
may include missing values in semi-supervised settings.

## Examples

``` r
S <- diag(3)
dimnames(S) <- list(letters[1:3], letters[1:3])
simfun <- createSimFun(S)
M <- simfun(c("a", NA, "c"))
Matrix::nnzero(M)
#> [1] 0
```
