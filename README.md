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

# Align with KEMA
result <- kema(hd, ncomp = 5)

# Or align with Generalized Procrustes
result <- generalized_procrustes(hd, ncomp = 5)
```

## Unified Aligner Interface

All pairwise alignment methods share a consistent adapter interface:

```r
# Create an aligner
aligner <- parrot_aligner(ncomp = 5)

# Fit a pair of domains
fit <- fit_pair(aligner, X1, X2, labels1, labels2)

# Align many domains at once
fit_many(aligner, list(X1, X2, X3), list(l1, l2, l3))
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
