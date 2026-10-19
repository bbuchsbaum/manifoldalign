# Out-of-sample RAPID-MA embedding and prediction

Out-of-sample RAPID-MA embedding and prediction

## Usage

``` r
# S3 method for class 'rapid_ma'
oos_predict(
  fit_or_transform,
  newX,
  side,
  type = c("embedding", "probabilities", "class", "all"),
  k = 8L,
  zero_tolerance = sqrt(.Machine$double.eps),
  new_positions = NULL,
  new_attributes = NULL,
  view_weights = c(feature = 1),
  ...
)
```

## Arguments

- fit_or_transform:

  A fitted \`rapid_ma\` or \`rapid_ma_pair_fit\`.

- newX:

  New rows in the original feature space of \`side\`.

- side:

  Domain name, index, or pairwise alias \`"i"\`/\`"j"\`.

- type:

  Return embedding coordinates, class probabilities, class labels, or
  all diagnostics.

- k:

  Number of bounded interpolation neighbours.

- zero_tolerance:

  Distance treated as an exact in-sample reapplication.

- new_positions:

  Optional position matrix for structure-aware OOS interpolation.

- new_attributes:

  Optional attribute data frame for structure-aware OOS interpolation.

- view_weights:

  Named nonnegative weights over \`feature\`, \`position\`, and
  \`attribute\`. The default is feature-only for backward compatibility.

- ...:

  Unused.

## Value

A matrix/vector, or a diagnostic list when \`type = "all"\`.
