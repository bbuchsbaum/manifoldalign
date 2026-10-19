# manifoldalign

**manifoldalign** is an R package for aligning data from multiple
domains into a shared latent space. It implements a broad family of
manifold alignment and domain adaptation methods — from kernel-based and
spectral approaches to optimal transport — preserving local manifold
structure while matching samples across domains.

## Alignment methods

No method is currently classified as stable or recommended. The table
below states the input model instead of treating every method as an
observed-graph aligner. See the [complete capability
inventory](https://bbuchsbaum.github.io/manifoldalign/METHOD_COMPARISON_TABLE.html)
and [lifecycle decision
record](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/ALIGNMENT_METHOD_SUNSET_REVIEW.md).

| Method | Function | Input contract | Lifecycle |
|----|----|----|----|
| **RAPID-MA** | [`rapid_ma()`](https://bbuchsbaum.github.io/manifoldalign/reference/rapid_ma.md), [`rapid_ma_feature_scores()`](https://bbuchsbaum.github.io/manifoldalign/reference/rapid_ma_feature_scores.md) | Supplied sparse weighted relations, positions, and mixed attributes | Experimental |
| **KEMA** | [`kema()`](https://bbuchsbaum.github.io/manifoldalign/reference/kema.md) | Supervised feature-derived kernel manifolds | Experimental |
| **GPCA** | [`gpca_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/gpca_align.md) | Supervised feature manifolds | Experimental |
| **GPA** | [`generalized_procrustes()`](https://bbuchsbaum.github.io/manifoldalign/reference/generalized_procrustes.md) | Orthogonal feature-task alignment with partial observations | Experimental |
| **CONE** | [`cone_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/cone_align.md), [`cone_align_multiple()`](https://bbuchsbaum.github.io/manifoldalign/reference/cone_align_multiple.md) | Spectral graphs derived from features | Legacy/reference |
| **GlobalGeo** | [`global_geo_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/global_geo_align.md) | Correspondence-linked geometry derived from features | Legacy/reference |
| **GRASP** | [`grasp()`](https://bbuchsbaum.github.io/manifoldalign/reference/grasp.md), [`grasp_multiset()`](https://bbuchsbaum.github.io/manifoldalign/reference/grasp_multiset.md) | Functional maps on feature-derived graphs | Experimental |
| **PARROT** | [`parrot()`](https://bbuchsbaum.github.io/manifoldalign/reference/parrot.md) | Anchor-guided dense OT on feature-derived networks | Experimental |
| **FPGW** | [`fpgw()`](https://bbuchsbaum.github.io/manifoldalign/reference/fpgw.md) | Partial metric-measure OT | Experimental |
| **Gromov-Wasserstein** | [`gromov_wasserstein()`](https://bbuchsbaum.github.io/manifoldalign/reference/gromov_wasserstein.md) | Balanced metric-measure OT | Experimental |
| **Low-Rank** | [`lowrank_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/lowrank_align.md) | Supervised low-rank feature manifolds | Experimental |
| **MMA** | [`mma_align_multiple()`](https://bbuchsbaum.github.io/manifoldalign/reference/mma_align_multiple.md) | Multiset feature-derived manifolds | Experimental |
| **Coupled Diag.** | [`coupled_diagonalization()`](https://bbuchsbaum.github.io/manifoldalign/reference/coupled_diagonalization.md) | Joint feature-derived spectral bases | Experimental |
| **Spectral MNN** | [`spectral_mnn_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/spectral_mnn_align.md) | Feature-derived spectral graphs | Experimental |
| **Token-OT Graph** | [`token_ot_graph_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/token_ot_graph_align.md) | Sparse OT on feature-derived graphs | Experimental |
| **SSMA** | [`ssma_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/ssma_align.md) | Feature-derived manifolds with sample and linear feature-summary constraints | Experimental |
| **Multiset UOT** | [`multiset_uot_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/multiset_uot_align.md) | Unbalanced metric-measure/template OT | Experimental |
| **LSE** | [`linear_sim_embed()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_sim_embed.md) | Single-domain target-similarity embedding | Experimental utility |
| **Multiscale** | [`multiscale_manifold_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/multiscale_manifold_align.md) | Incoherent experimental hierarchy retained for reproduction | Deprecated; removal 1.0.0 |

## Linear constraints in SSMA

`ssma_align(..., constraints = list(...))` accepts
[`linear_correspondence()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_correspondence.md)
objects for weighted agreement between linear feature summaries and
[`linear_penalty()`](https://bbuchsbaum.github.io/manifoldalign/reference/linear_penalty.md)
objects for within-domain regularization. These add quadratic penalties
without inserting synthetic observations or changing sample-graph
normalization. Unequal feature counts, sparse maps, and external PCA
composition are supported. Maps refer to processed feature coefficients;
rank-truncation loss and fitted residuals are reported per constraint.
See the [worked
vignette](https://bbuchsbaum.github.io/manifoldalign/vignettes/ssma-linear-constraints.Rmd)
for coordinate, scaling, and validation contracts. SSMA remains
experimental.

## Current evidence

The frozen alignment tournament contains 588 typed rows: 42
attributed-graph cases by 14 methods, with zero missing result rows and
zero validation problems. It records failures, unsupported tasks,
runtime, peak RSS, controls, paired effects, and checksums for every
published receipt file. Resumable binary checkpoints and the debug log
remain local operational state. Read the [checked
interpretation](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/INTERPRETATION.md)
or inspect the [raw
results](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/raw_results.csv).

RAPID-MA completed 42/42 and is the only tournament method that consumed
the supplied weighted relations together with positions and mixed
attributes. It is not a universal accuracy winner: a position-only
control beat its node matching in every non-oracle regime. Its frozen
downstream classification probe did beat source-only transfer in every
non-oracle regime. The current full-scoring scale receipt reaches 22,600
total nodes; an older runner reached 220,000, while the one-million-node
and blocked-inductive promotion gates failed. Those limitations are part
of the public evidence contract.

## Installation

Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("bbuchsbaum/manifoldalign")
```

## Quick Start

``` r

library(manifoldalign)

# Create two domains with shared structure
set.seed(1)
X1 <- matrix(rnorm(100 * 20), 100, 20)
X2 <- matrix(rnorm(80 * 15), 80, 15)
labels1 <- rep(1:4, each = 25)
labels2 <- rep(1:4, each = 20)

# Build a hyperdesign (multi-domain container)
hd <- as_hyperdesign(
  list(X1, X2),
  labels = list(labels1, labels2)
)

# Align with KEMA (an experimental supervised feature-manifold method)
result <- kema(hd, y = label, ncomp = 5)

# Or align with Generalized Procrustes
result <- generalized_procrustes(hd, ncomp = 5)
```

## Dual (feature-side) coordinates

Embedding methods place samples in the shared latent space. For
RAPID-MA, which is not a linear projector and therefore has no loadings,
[`rapid_ma_feature_scores()`](https://bbuchsbaum.github.io/manifoldalign/reference/rapid_ma_feature_scores.md)
places each feature column (for example an activation map over voxels)
at the weighted barycentre of the samples it loads on, so features from
different domains can be compared in the same space:

``` r

fit  <- rapid_ma(X, relations = relations, ncomp = 16)
dual <- rapid_ma_feature_scores(fit, X)
dual$scores$domain1   # feature x component coordinates for one domain
dual$pooled           # equal-weight mean across domains
dual$dispersion       # how consistently each feature is placed across domains
```

## Unified aligner interface

Adapter-backed methods share a consistent interface:

``` r

# Create an aligner
aligner <- grasp_aligner()

# Fit a pair of domains
fit <- fit_pair(aligner, X1, X2, ncomp = 5)

# Use a native multi-domain adapter when aligning more than two domains
X3 <- matrix(rnorm(90 * 10), 90, 10)
multi_aligner <- spectral_mnn_aligner()
fit_many(multi_aligner, list(X1, X2, X3), ncomp = 5)
```

## Articles

Detailed tutorials are available as package vignettes:

- [Quickstart: Building a hyperdesign and running
  KEMA/GPA](https://bbuchsbaum.github.io/manifoldalign/articles/quickstart-hyperdesign.html)
- [Features, Correspondences, and Predictive
  Performance](https://bbuchsbaum.github.io/manifoldalign/articles/features-and-predictive-performance.html)
- [Benchmark Guidance by Problem
  Regime](https://bbuchsbaum.github.io/manifoldalign/articles/benchmark-guidance.html)
- [Kernel Manifold Alignment with
  KEMA](https://bbuchsbaum.github.io/manifoldalign/articles/kema-overview.html)
- [Generalized PCA
  Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/gpca-align.html)
- [Pairwise Graph Alignment with
  CONE](https://bbuchsbaum.github.io/manifoldalign/articles/cone-align.html)
- [Multi-Domain Graph Alignment with
  CONE](https://bbuchsbaum.github.io/manifoldalign/articles/cone-align-multiple.html)
- [Spectral Graph Alignment with
  GRASP](https://bbuchsbaum.github.io/manifoldalign/articles/grasp-align.html)
- [Position-Aware OT Alignment with
  PARROT](https://bbuchsbaum.github.io/manifoldalign/articles/parrot-align.html)
- [Fused-Partial
  Gromov-Wasserstein](https://bbuchsbaum.github.io/manifoldalign/articles/fpgw_tutorial.html)
- [Low-Rank
  Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/lowrank-align.html)
- [Multiset Manifold
  Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/mma-align.html)
- [Coupled
  Diagonalization](https://bbuchsbaum.github.io/manifoldalign/articles/coupled-diagonalization.html)
- [Performance
  Overview](https://bbuchsbaum.github.io/manifoldalign/articles/performance-overview.html)

## Documentation

Full API reference and articles at
<https://bbuchsbaum.github.io/manifoldalign/>.

## License

MIT

## Albers theme

This package uses the albersdown theme. Existing vignette theme hooks
are replaced so `albers.css` and local `albers.js` render consistently
on CRAN and GitHub Pages. The palette family is provided via
`params$family` (default ‘red’). The pkgdown site uses
`template: { package: albersdown }`.
