# Fit RAPID-MA to a pair of domains

Fit RAPID-MA to a pair of domains

## Usage

``` r
# S3 method for class 'rapid_ma_aligner'
fit_pair(
  algo,
  X_i,
  X_j,
  links = NULL,
  labels = NULL,
  relations = NULL,
  positions = NULL,
  attributes = NULL,
  ncomp = 16L,
  control = rapid_ma_control(),
  fine_match = TRUE,
  match_control = list(),
  ...
)
```

## Arguments

- algo:

  A \`rapid_ma_aligner\` descriptor.

- X_i, X_j:

  Numeric domain matrices.

- links:

  Optional one-to-one anchor pairs or paired anchor-ID vectors.

- labels, relations, positions, attributes:

  Optional two-domain inputs.

- ncomp:

  Shared embedding dimension.

- control:

  Settings from \[rapid_ma_control()\].

- fine_match:

  Whether to compute bounded node-level matches.

- match_control:

  Named arguments forwarded to \[rapid_ma_match()\].

- ...:

  Unused.

## Value

A \`rapid_ma_pair_fit\` containing the shared core fit and fine match.
