# manifoldalign

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/manifoldalign/)
<!-- badges: end -->

**manifoldalign** is an R package for aligning data from multiple domains into a shared latent space. It implements a broad family of manifold alignment and domain adaptation methods — from kernel-based and spectral approaches to optimal transport — preserving local manifold structure while matching samples across domains.

## Alignment Methods

| Method | Function | Description |
|--------|----------|-------------|
| **KEMA** | `kema()` | Kernel Manifold Alignment with exact and regression-based solvers |
| **GPCA** | `gpca_align()` | Generalized PCA alignment via shared subspace estimation |
| **GPA** | `generalized_procrustes()` | Generalized Procrustes Analysis for partial observations |
| **CONE** | `cone_align()`, `cone_align_multiple()` | Spectral graph alignment with anchor-based matching |
| **GRASP** | `grasp()`, `grasp_multiset()` | Graph Spectral alignment using structural descriptors |
| **PARROT** | `parrot()` | Position-Aware Regularized OT alignment |
| **FPGW** | `fpgw()` | Fused Partial Gromov-Wasserstein transport |
| **Gromov-Wasserstein** | `gromov_wasserstein()` | Structure-preserving optimal transport |
| **Low-Rank** | `lowrank_align()` | Memory-efficient alignment via low-rank factorization |
| **MMA** | `mma_align_multiple()` | Multiset Manifold Alignment for >2 domains |
| **Coupled Diag.** | `coupled_diagonalization()` | Joint diagonalization across modalities |
| **Spectral MNN** | `spectral_mnn_align()` | Mutual nearest neighbor spectral alignment |
| **Token-OT Graph** | `token_ot_graph_align()` | Token-level OT with graph structure |
| **SSMA** | `ssma_align()` | Semi-Supervised Manifold Alignment |
| **GlobalGeo** | `global_geo_align()` | Landmark geodesic alignment with optional coordinate-anchored column correspondence |
| **Multiscale** | `multiscale_manifold_align()` | Multiscale manifold alignment over a joint spectral hierarchy |
| **LSE** | `linear_sim_embed()` | Linear similarity embedding |
| **Multi-set UOT** | `multiset_uot_align()`, `uot_fit_pair()` | Unbalanced optimal transport of several datasets to a shared template |
| **Synchronisation** | `align_many()` | Lifts pairwise aligners to many domains by transform synchronisation |
| **SUMCOR** | `sumcor_align()`, `calibrate_domain()` | Event-keyed multi-domain SUMCOR with repeated observations and an optional event-feature metric; explicit per-domain encoders and frozen-reference calibration of new domains |
| **RAPID-MA** | `rapid_ma()`, `rapid_ma_feature_scores()` | Scalable structure-aware alignment of many unmatched domains with custom sparse relations; feature-side (dual) coordinates via barycentres |

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("bbuchsbaum/manifoldalign")
```

## Quick Start

```r
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

# Align with KEMA (as_hyperdesign() stores the labels in a column named `label`)
result <- kema(hd, y = label, ncomp = 5)

# Or with semi-supervised manifold alignment
result <- ssma_align(hd, y = label, ncomp = 5)
```

## Event-keyed SUMCOR alignment

`sumcor_align()` learns one linear encoder per domain from observations of
shared events: keyed items, conditions or time points that several domains
observed, possibly more than once. Correspondence is carried only by
`event_id` values, so domains may have different column counts, missing
events, or different row orders. Repeated observations of an event are
averaged in the agreement term but enter each domain's normalisation
individually. Optional event features add a metric over events; `beta = 0`
is the exact correspondence-only model, so the effect of the features can be
tested directly.

```r
fit <- sumcor_align(
  X_list, design = design_list,
  event_id = "event_id", repeat_id = "repeat_id",
  event_features = list(f1 = F1, f2 = F2),
  ncomp = 10,
  control = sumcor_control(beta = 0.5, covariance_shrinkage = 0.1)
)

# Encode new observations of one domain
Z_new <- oos_predict(fit, X_new, side = "domain3")

# Add a held-out domain against the frozen reference
fit2 <- calibrate_domain(fit, X_heldout, design = design_heldout,
                         domain_name = "heldout")
```

## Dual (feature-side) coordinates

Embedding methods place samples in the shared latent space. For RAPID-MA, which
is not a linear projector and therefore has no loadings, `rapid_ma_feature_scores()`
places each feature column at the weighted barycentre of the samples it loads
on, so features from different domains can be compared in the same space:

```r
fit  <- rapid_ma(X, relations = relations, ncomp = 16)
dual <- rapid_ma_feature_scores(fit, X)
dual$scores$domain1   # feature x component coordinates for one domain
dual$pooled           # equal-weight mean across domains
dual$dispersion       # how consistently each feature is placed across domains
```

## Unified Aligner Interface

Methods share an adapter interface. Aligner constructors take no arguments;
method settings are passed to `fit_pair()` or `fit_many()`, and
`aligner_capabilities()` reports what each aligner supports.

```r
# Native multi-domain aligners take a hyperdesign
fit <- fit_many(kema_aligner(), hd, y = label, ncomp = 5)

# Pairwise aligners take two matrices; PARROT needs the same columns in both
Y1 <- matrix(rnorm(40 * 10), 40, 10)
Y2 <- Y1[sample(40), ] + matrix(rnorm(400, sd = 0.1), 40, 10)
pair <- fit_pair(parrot_aligner(), Y1, Y2, ncomp = 3)

aligner_capabilities(parrot_aligner())
```

## Articles

Detailed tutorials are available as package vignettes:

- [Quickstart: Building a hyperdesign and running KEMA/GPA](https://bbuchsbaum.github.io/manifoldalign/articles/quickstart-hyperdesign.html)
- [Features, Correspondences, and Predictive Performance](https://bbuchsbaum.github.io/manifoldalign/articles/features-and-predictive-performance.html)
- [Kernel Manifold Alignment with KEMA](https://bbuchsbaum.github.io/manifoldalign/articles/kema-overview.html)
- [Generalized PCA Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/gpca-align.html)
- [Pairwise Graph Alignment with CONE](https://bbuchsbaum.github.io/manifoldalign/articles/cone-align.html)
- [Multi-Domain Graph Alignment with CONE](https://bbuchsbaum.github.io/manifoldalign/articles/cone-align-multiple.html)
- [Spectral Graph Alignment with GRASP](https://bbuchsbaum.github.io/manifoldalign/articles/grasp-align.html)
- [Position-Aware OT Alignment with PARROT](https://bbuchsbaum.github.io/manifoldalign/articles/parrot-align.html)
- [Fused-Partial Gromov-Wasserstein](https://bbuchsbaum.github.io/manifoldalign/articles/fpgw_tutorial.html)
- [Low-Rank Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/lowrank-align.html)
- [Multiset Manifold Alignment](https://bbuchsbaum.github.io/manifoldalign/articles/mma-align.html)
- [Coupled Diagonalization](https://bbuchsbaum.github.io/manifoldalign/articles/coupled-diagonalization.html)
- [Performance Overview](https://bbuchsbaum.github.io/manifoldalign/articles/performance-overview.html)

## Documentation

Full API reference and articles at <https://bbuchsbaum.github.io/manifoldalign/>.

## License

MIT

<!-- albersdown:theme-note:start -->
## Albers theme
This package uses the albersdown theme. Existing vignette theme hooks are replaced so `albers.css` and local `albers.js` render consistently on CRAN and GitHub Pages. The palette family is provided via `params$family` (default 'red'). The pkgdown site uses `template: { package: albersdown }`.
<!-- albersdown:theme-note:end -->
