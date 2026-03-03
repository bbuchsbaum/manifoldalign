# Construct a diffusion-time grid for HKS/diffusion descriptors

Construct a diffusion-time grid for HKS/diffusion descriptors

## Usage

``` r
hks_time_grid(
  values,
  q,
  mode = c("auto", "fixed"),
  time_range = c(0.1, 50),
  spacing = c("log", "linear"),
  time_scale = 1,
  eigen_tol = 1e-12
)
```

## Arguments

- values:

  Numeric vector of eigenvalues (typically Laplacian eigenvalues).

- q:

  Integer number of time points.

- mode:

  Either \`"auto"\` (derive a time range from the spectrum) or
  \`"fixed"\` (use \`time_range\`).

- time_range:

  Length-2 positive numeric vector giving min/max times when
  \`mode="fixed"\`.

- spacing:

  Either \`"log"\` (default; typical for HKS) or \`"linear"\`.

- time_scale:

  Positive scalar multiplier applied to the resulting times.

- eigen_tol:

  Positive threshold; eigenvalues \`\<= eigen_tol\` are ignored when
  \`mode="auto"\`.

## Value

Numeric vector of length \`q\`.
