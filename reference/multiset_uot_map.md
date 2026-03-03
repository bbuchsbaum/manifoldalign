# Map multiple subjects into template space from a multi-set fit

Convenience helper for \[multiset_uot_align()\]. It applies each
subject's fitted UOT map to the corresponding signal, using the
\*\*final\*\* template weights and the pairwise potentials stored in
\`fit\$fits\`.

## Usage

``` r
multiset_uot_map(fit, datasets, signals, epsilon = NULL, delta = 1e-08)
```

## Arguments

- fit:

  A result from \[multiset_uot_align()\].

- datasets:

  The datasets list originally passed to \[multiset_uot_align()\]. Each
  dataset must contain \`alpha\`.

- signals:

  A list of length K of signals, each a numeric vector (length N_k) or
  matrix (T x N_k).

- epsilon:

  Optional override of the entropic parameter. If NULL, uses
  \`fit\$params\$epsilon\`.

- delta:

  Stabilizer for the barycentric denominator.

## Value

A list of mapped signals (length K), each in template space.
