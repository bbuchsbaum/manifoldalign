# Linear correspondence constraints for alignment coefficients

Describe soft agreement between linear summaries of two domains'
projection coefficients. For maps H1 and H2, coefficients W1 and W2, and
effective row weights d, the added objective is
`strength * sum(d * rowSums((H1 %*% W1 - H2 %*% W2)^2))`.

## Usage

``` r
linear_correspondence(
  domains,
  maps,
  weights = 1,
  strength = 1,
  reduction = c("weighted_mean", "sum"),
  name = "correspondence",
  space = "processed"
)
```

## Arguments

- domains:

  Character vector of two distinct domain names, matching the names in
  the alignment input.

- maps:

  List of two finite numeric matrices (base or `Matrix`), with equal
  positive row counts. Columns index each domain's processed input
  features. Different domains may have different numbers of features. If
  row names are supplied, both maps must have identical, unique row
  names in order. Named columns must exactly match the processed data's
  feature names.

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

A `linear_alignment_constraint` object. Pass a list of these objects to
[`ssma_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/ssma_align.md)
via `constraints`. Currently SSMA is the consumer.

## Details

These are homogeneous quadratic soft constraints, not additional samples
or hard equalities. They do not modify sample neighborhoods or their
degree normalization. Selecting feature coefficients, averaging groups,
and applying linear measurement maps all use the same interface. For
processed data X, maps `A %*% X` express correspondence between weighted
sample summaries A.

With externally reduced input `Z = X %*% P`, use maps `H %*% P` and
compose returned coefficients as `P %*% W`. With covariance maps,
calculate `G %*% crossprod(X, Z) / (nrow(X) - 1)` for centered data
instead of forming the full feature covariance. Covariance estimation
and preprocessing belong to the caller and must respect validation
folds.

Unnamed columns are interpreted positionally. Raw-space constraints are
not automatically transformed through arbitrary preprocessing. For a
known linear preprocessing map S, a raw-coefficient functional H becomes
H S; intercepts are outside this homogeneous coefficient contract.

## Examples

``` r
h <- diag(3)
link <- linear_correspondence(c("a", "b"), list(h, h), name = "features")
smooth <- linear_penalty("a", diff(diag(3)), name = "smoothness")
```
