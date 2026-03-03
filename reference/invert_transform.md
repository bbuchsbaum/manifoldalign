# Invert a transform

Invert a transform

## Usage

``` r
invert_transform(transform, ...)
```

## Arguments

- transform:

  An alignment transform object to invert

- ...:

  Additional arguments passed to methods

## Value

An align_transform object representing the inverse transform, with from
and to swapped.

## Examples

``` r
# Create and invert an orthogonal transform
R <- matrix(c(cos(pi/4), -sin(pi/4), sin(pi/4), cos(pi/4)), 2, 2)
tr <- new_align_transform("O", R, from = 1, to = 2)
tr_inv <- invert_transform(tr)
```
