# SSMA with linear correspondence constraints

SSMA aligns observations through domain-specific linear projections.
Additional knowledge may concern features, groups of features, or
averages of observations.
[`linear_correspondence()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_correspondence.md)
represents these as weighted linear summaries of the projection
coefficients. It works with unequal feature counts and partial links; no
one-to-one correspondence between entire feature spaces is assumed. SSMA
remains experimental: numerical agreement is not evidence of predictive
benefit.

## Objective and coordinate contract

For processed observations $`X_s`$, let $`Y_s=X_sW_s`$. A correspondence
adds

``` math
\lambda\sum_q d_q\|(H_sW_s)_q-(H_tW_t)_q\|^2
```

to SSMA’s existing objective. It changes neither the observation graph
nor its degree-based normalization. A single-domain
[`linear_penalty()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_penalty.md)
instead adds $`\lambda\sum_q d_q\|(H_sW_s)_q\|^2`$, useful for
smoothness or coefficient contrasts. These are soft, homogeneous
quadratic penalties, not hard equalities.

Maps act on **processed input features, before internal compression**.
They are never centered or standardized as observations. For arbitrary
preprocessing, construct maps in its resulting feature coordinates
explicitly; the package does not guess a raw-to-processed conversion. If
$`X_{processed}=(X-\mu)S`$, a functional $`H`$ of raw coefficients
becomes $`HS`$ in processed coordinates. Intercepts are not part of this
coefficient contract. Input-domain names and optional map row/column
names are checked; unnamed map columns have positional semantics.

## Matching feature summaries alongside observations

This generic example links summaries of two measurement systems with
five and seven features. The maps could select homologous genes, average
sensor groups, or describe linear calibration measurements.

``` r

x <- list(a = matrix(rnorm(24 * 5), 24, 5),
          b = matrix(rnorm(27 * 7), 27, 7))
h <- list(matrix(rnorm(4 * 5), 4, 5), matrix(rnorm(4 * 7), 4, 7))
link <- linear_correspondence(
  domains = c("a", "b"), maps = h,
  weights = c(1, 0, 2, 3), strength = 2, name = "measurements"
)
rows <- data.frame(domain_i = 1L, index_i = 1:12,
                   domain_j = 2L, index_j = 2:13)
fit <- ssma_align(
  x, correspondences = rows, constraints = list(link),
  preproc = NULL, ncomp = 3,
  control = ssma_align_control(knn = 8, sigma = 2, compression = "none")
)
fit$constraint_diagnostics
#> [[1]]
#> [[1]]$name
#> [1] "measurements"
#> 
#> [[1]]$type
#> [1] "correspondence"
#> 
#> [[1]]$domains
#> [1] "a" "b"
#> 
#> [[1]]$reduction
#> [1] "weighted_mean"
#> 
#> [[1]]$strength
#> [1] 2
#> 
#> [[1]]$n_rows
#> [1] 4
#> 
#> [[1]]$active_rows
#> [1] 3
#> 
#> [[1]]$effective_weight_sum
#> [1] 1
#> 
#> [[1]]$discarded_energy_fraction
#> a b 
#> 0 0 
#> 
#> [[1]]$residual_per_component
#> [1] 0.002227325 0.001564937 0.003796108
#> 
#> [[1]]$residual
#> [1] 0.00758837
#> 
#> [[1]]$objective
#> [1] 0.01517674
```

`weights` describe confidence in each matched summary. Zero confidence
excludes a row. The default `reduction="weighted_mean"` normalizes
confidence weights; `"sum"` preserves their magnitude. `strength` scales
the resulting objective. Neither operation equalizes map units.
Normalize measurement units deliberately, and choose strengths within
the permitted training folds. Multiple objects have additive objectives.
If a family spans many domain pairs, allocate its overall strength
across those pairs explicitly; more pairs otherwise add more weight.

Rows with missing measurements should be omitted from the maps and
weights, not filled with zeros. Named map rows must match in order.
Feature names, when supplied on a map, must match the corresponding
processed input in order.

## Observation averages and smooth correspondences

Given averaging matrices $`A_s`$ and $`A_t`$, use $`H_s=A_sX_s`$ and
$`H_t=A_tX_t`$. Rows can encode weighted windows, interpolation, or a
declared temporal smoothing kernel. One-hot selectors recover ordinary
sample-correspondence penalties when weights and reductions match.
Helpers outside the alignment method should choose the windows and
handle boundaries, missing samples, timing, and domain semantics.

Existing `correspondences` continue to support direct observation links.
Set `auto_correspondences=FALSE` when the intended sample graph has no
cross-domain links. This disables label/MNN discovery when explicit
links are absent, which is particularly important when differently sized
feature spaces have no common Euclidean coordinate system. Explicit
links still apply.

## Composition with PCA and covariance maps

With $`Z_s=X_sP_s`$, fit SSMA on $`Z_s`$, supply maps $`H_sP_s`$, and
recover original feature coefficients as $`W_s=P_sB_s`$. This is exact
in the retained subspace; truncation discards information and can change
the observation-neighbor graph. Preserve any scaling or whitening in the
composed mapping.

For spatial or grouped-feature summaries $`G_s`$, $`G_sP_s`$ constrains
averages of projection coefficients. To constrain averages of covariance
fingerprints $`F_s=C_sP_sB_s`$, instead supply $`G_sC_sP_s`$. For
centered real observations this can be computed without a full
feature-by-feature covariance matrix:

``` r

z <- x %*% p
h <- (g %*% t(x)) %*% z / (nrow(x) - 1)
fingerprint <- crossprod(x, z %*% b) / (nrow(x) - 1)
```

These are covariance fingerprints, not automatically standardized
correlations or a guarantee that all feature-connectivity distances are
preserved. Estimate covariances from actual observations, respecting
nuisance processing, weights, missingness and validation folds.
Projection weights and covariance patterns have different meanings.

`compression="none"` preserves already reduced features without
allocating an identity basis. With the default internal SVD, the package
projects maps through its basis and returns `compression_bases`; each
NULL entry denotes identity. `discarded_energy_fraction` reports the
unweighted Frobenius-energy fraction of each map outside the retained
basis. It is a map-support diagnostic, not neural variance explained.
Rank compression can remove useful constraints: inspect this quantity
rather than assuming anatomical or low-variance directions survived.

## Sparse maps, solvers, diagnostics and validation

Base matrices and numeric `Matrix` maps are accepted. Sparse maps remain
sparse until reduced products require dense values. The operator solver
applies maps and their transposes without assembling the global
constraint matrix. The dense solver forms only reduced blocks. With
compression disabled, the caller is responsible for choosing a
manageable input dimension. This is not a streaming data loader.

Diagnostics report active row counts, reduction, strength, effective
weight sum, per-component residuals, total residual and weighted
objective contribution for each named constraint. Small residuals alone
do not imply useful alignment: large strengths can suppress directions
and contradictory constraints cannot generally all be satisfied. Inspect
out-of-sample performance and appropriate shuffled-link controls
separately.

The tests include independent outer-product and generalized-eigenproblem
oracles, dense/operator agreement, sample-link equivalence, sparse and
unequal feature spaces, null penalties, truncation, confidence
splitting, contradictions, and permutation/rotation invariants. Fitted
domains can be transformed using their own learned projections. This
extension does not create a new-domain adaptation protocol; any use of
an unlabeled target during fitting must be declared.
