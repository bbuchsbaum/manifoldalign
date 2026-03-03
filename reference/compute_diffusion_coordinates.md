# Compute diffusion coordinates from a spectral basis

Diffusion coordinates at time \\t\\ are given by \\\Phi
\mathrm{diag}(\exp(-\lambda t))\\, where \\(\Phi,\lambda)\\ are
Laplacian eigenvectors/eigenvalues.

## Usage

``` r
compute_diffusion_coordinates(
  basis = NULL,
  values = NULL,
  vectors = NULL,
  time = 1,
  k = NULL,
  eigen_tol = 1e-12,
  normalize = c("none", "row")
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

- time:

  Positive numeric scalar or vector of times.

- k:

  Optional integer number of coordinates to keep (defaults to all).

- eigen_tol:

  Positive threshold; eigenpairs with \`value \<= eigen_tol\` are
  dropped.

- normalize:

  Optional normalization: \`"none"\` (default) or \`"row"\`.

## Value

If \`time\` is length-1, a numeric matrix N x k. If \`time\` has length
\> 1, a list of matrices (one per time).
