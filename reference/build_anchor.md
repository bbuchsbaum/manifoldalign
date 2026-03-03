# Build Anchor Representation

Constructs the anchor coefficient matrix either from a specific graph or
as the mean of all aligned representations.

## Usage

``` r
build_anchor(desc, bases, Mlist, anchor)
```

## Arguments

- desc:

  List of descriptor matrices

- bases:

  List of spectral bases

- Mlist:

  List of rotation matrices

- anchor:

  Index of anchor graph or "mean"

## Value

Anchor coefficient matrix
