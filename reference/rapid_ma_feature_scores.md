# Feature-side (dual) coordinates for a RAPID-MA embedding

RAPID-MA embeds samples (for example voxels) into a shared latent space
but is not a linear projector, so it has no feature loadings. This
function supplies the dual view instead: every feature column (for
example an activation map over voxels) is placed at the weighted
barycentre of the latent coordinates of the samples it loads on. Because
the latent space is shared across domains, feature coordinates from
different domains are directly comparable, and features that distribute
similarly over similarly embedded samples land close together.

## Usage

``` r
rapid_ma_feature_scores(
  fit,
  data,
  weight = c("positive", "absolute", "shifted"),
  pool = TRUE
)
```

## Arguments

- fit:

  A fitted \`rapid_ma\` object or a \`rapid_ma_pair_fit\`.

- data:

  A list of sample-by-feature matrices (base or \`Matrix\`), one per
  fitted domain in fitting order or named by domain, or a
  \`hyperdesign\`. Entries may be \`NULL\` to skip a domain. For a
  \`rapid_ma_pair_fit\`, supply \`list(X_i, X_j)\` or \`list(i = X_i, j
  = X_j)\`.

- weight:

  Nonnegative weighting rule applied column by column: \`"positive"\`,
  \`"absolute"\`, or \`"shifted"\`.

- pool:

  Logical; also return the equal-weight mean of each feature's
  per-domain coordinates across the supplied domains, plus a dispersion
  diagnostic.

## Value

An object of class \`rapid_ma_feature_scores\` with elements \`scores\`
(per-domain feature-by-component matrices), \`mass\` (per-domain total
weight per feature), \`zero_mass\` (per-domain indices of features
without a barycentre), \`pooled\` (feature-by-component matrix, or
\`NULL\` when \`pool = FALSE\`), \`dispersion\` (per-feature
root-mean-square distance of the per-domain coordinates from the pooled
coordinate; \`NA\` when fewer than two domains contribute),
\`n_domains\` (per-feature count of contributing domains), \`weight\`,
\`domain_names\`, and \`ncomp\`.

## Details

For feature \\j\\ in domain \\d\\ with sample coordinates \\Z_d\\ and
nonnegative weights \\w\_{ij}\\ derived from the data, the coordinate is
\\m\_{dj} = \sum_i w\_{ij} z\_{di} / \sum_i w\_{ij}\\. Weights are
derived from the data column by column according to \`weight\`:
\`"positive"\` keeps nonnegative values and zeroes negative ones,
\`"absolute"\` uses magnitudes, and \`"shifted"\` subtracts each
column's minimum so every entry is nonnegative. Sparse \`Matrix\` inputs
stay sparse for \`"positive"\` and \`"absolute"\`; \`"shifted"\`
densifies.

Columns with zero total weight have no barycentre and are returned as
\`NA\` rows; their indices are reported in \`zero_mass\` and a warning
is raised.

The matrices in \`data\` need not be the training matrices: any matrix
with one row per fitted sample of that domain can be embedded, so new
maps over the same samples can be placed after fitting. Pooling requires
every supplied domain to have the same number of columns (and identical
column names when both sides are named).

## Examples

``` r
set.seed(1)
X <- list(
  a = matrix(abs(rnorm(40 * 6)), 40, 6),
  b = matrix(abs(rnorm(35 * 6)), 35, 6)
)
fit <- rapid_ma(X, ncomp = 3L,
                control = rapid_ma_control(max_iter = 1L,
                                           prototype = list(unknown_prototypes = 4L)))
dual <- rapid_ma_feature_scores(fit, X)
dim(dual$pooled)
#> [1] 6 3
```
