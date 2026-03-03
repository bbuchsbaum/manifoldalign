# Compute Heat Kernel Signatures (HKS) from a spectral basis

Compute Heat Kernel Signatures (HKS) from a spectral basis

## Usage

``` r
compute_hks_descriptors(
  basis = NULL,
  values = NULL,
  vectors = NULL,
  q = 16L,
  times = NULL,
  time_mode = c("auto", "fixed"),
  time_range = c(0.1, 50),
  spacing = c("log", "linear"),
  time_scale = 1,
  eigen_tol = 1e-12,
  normalize = c("both", "column", "row", "none"),
  clamp_exponent = -745
)
```

## Arguments

- basis:

  Optional basis list with fields \`values\` and \`vectors\`.

- values:

  Optional numeric vector of eigenvalues (overrides \`basis\$values\`).

- vectors:

  Optional matrix (or sparse Matrix) of eigenvectors (overrides
  \`basis\$vectors\`).

- q:

  Integer number of time points.

- times:

  Optional numeric time grid. If provided, \`time\_\*\` arguments are
  ignored.

- time_mode:

  Either \`"auto"\` (derive a time range from the spectrum) or
  \`"fixed"\` (use \`time_range\`).

- time_range:

  Length-2 positive numeric vector giving min/max times when
  \`time_mode="fixed"\`.

- spacing:

  Either \`"log"\` (default; typical for HKS) or \`"linear"\`.

- time_scale:

  Positive scalar multiplier applied to the time grid.

- eigen_tol:

  Positive threshold; eigenpairs with \`value \<= eigen_tol\` are
  dropped.

- normalize:

  Descriptor normalization: \`"both"\` (default), \`"column"\`,
  \`"row"\`, or \`"none"\`.

- clamp_exponent:

  Lower clamp for the exponent (default \`-745\`), used to avoid
  denormal underflow in \`exp(-λ t)\`.

## Value

A numeric matrix of dimension N x q.
