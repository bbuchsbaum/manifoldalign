# Construct an alignment transform object

Construct an alignment transform object

## Usage

``` r
new_align_transform(group, op, from, to, k = NULL, dim = NULL)
```

## Arguments

- group:

  Group identifier: "O", "GL", or "perm" (or "kernel" surrogate)

- op:

  The underlying operator (e.g., kxk matrix, or assignment matrix)

- from:

  Source index/label

- to:

  Target index/label

- k:

  Latent dimension (if applicable)

- dim:

  Optional problem dimensions metadata

## Value

An object of class c("align_transform\_\<group\>", "align_transform")

## Examples

``` r
# Create a simple orthogonal transform
R <- diag(3)
tr <- new_align_transform("O", R, from = 1, to = 2)
```
