# Linear regularization of alignment coefficients

A single-domain counterpart to
[`linear_correspondence()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_correspondence.md).
Adds `strength * sum(d * rowSums((map %*% W)^2))` to the objective. For
example, a feature-graph incidence matrix penalizes differences between
neighboring feature coefficients. This does not establish preservation
of the complete feature covariance geometry.

## Usage

``` r
linear_penalty(
  domain,
  map,
  weights = 1,
  strength = 1,
  reduction = c("weighted_mean", "sum"),
  name = "penalty",
  space = "processed"
)
```

## Arguments

- domain:

  One domain name.

- map:

  A finite numeric matrix whose columns index processed features.

- weights:

  Nonnegative finite confidence weights, one per matched row, or a
  scalar recycled across rows. Zero weights omit individual links.

- strength:

  Finite nonnegative multiplier for this constraint.

- reduction:

  `"weighted_mean"` divides weights by their sum; `"sum"` retains their
  absolute magnitude. Neither option rescales the maps or harmonizes
  their units. An all-zero weighted mean is undefined.

- name:

  Nonempty label for diagnostics; labels must be unique within a fit.

- space:

  Coordinate contract. Currently only `"processed"` is supported: maps
  act on coefficients after input preprocessing but before SSMA's
  internal compression. Maps are never centered or standardized as
  samples.

## Value

A `linear_alignment_constraint` object.
