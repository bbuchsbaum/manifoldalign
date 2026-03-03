# Map signals into template space using a \`uot_pair_fit\`

Map signals into template space using a \`uot_pair_fit\`

## Usage

``` r
uot_map(fit, signal, delta = 1e-08)
```

## Arguments

- fit:

  An object returned by \[uot_fit_pair()\].

- signal:

  Numeric vector (length n) or matrix (T x n).

- delta:

  Stabilizer for the barycentric denominator.

## Value

A numeric vector (length m) or matrix (T x m).
