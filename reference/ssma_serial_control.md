# Control settings for serial-correlation handling in SSMA graph construction

Control settings for serial-correlation handling in SSMA graph
construction

## Usage

``` r
ssma_serial_control(
  enabled = FALSE,
  row_whiten = c("none", "ar1"),
  ar1_shrink = 0,
  ar1_clip = 0.98,
  lag_mode = c("none", "hard", "exponential", "gaussian"),
  lag_exclusion = 0L,
  lag_scale = 2,
  preserve_degree = TRUE
)
```

## Arguments

- enabled:

  Logical; enable serial-correlation adjustments when \`serial_index\`
  is provided.

- row_whiten:

  One of \`"none"\` or \`"ar1"\`. Controls optional row-wise whitening
  applied before within-domain kNN graph construction.

- ar1_shrink:

  Scalar in \`\[0, 1)\` used to shrink estimated AR(1) coefficient
  toward zero.

- ar1_clip:

  Maximum absolute AR(1) coefficient allowed for stability.

- lag_mode:

  One of \`"none"\`, \`"hard"\`, \`"exponential"\`, or \`"gaussian"\`.
  Controls how short-lag neighbor edges are downweighted.

- lag_exclusion:

  Non-negative integer lag exclusion radius used when \`lag_mode =
  "hard"\`.

- lag_scale:

  Positive scalar lag scale used for soft lag kernels.

- preserve_degree:

  Logical; if \`TRUE\`, rescale filtered adjacency to keep row-degree
  magnitudes near pre-filter values.

## Value

A list of class \`ssma_serial_control\`.
