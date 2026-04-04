# manifoldalign

**manifoldalign** is an R package for aligning data from multiple
domains into a shared latent space. It implements a broad family of
manifold alignment and domain adaptation methods — from kernel-based and
spectral approaches to optimal transport — preserving local manifold
structure while matching samples across domains.

## Alignment Methods

| Method                 | Function                                                                                                                                                                                     | Description                                                       |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| **KEMA**               | [`kema()`](https://bbuchsbaum.github.io/manifoldalign/reference/kema.md)                                                                                                                     | Kernel Manifold Alignment with exact and regression-based solvers |
| **GPCA**               | [`gpca_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/gpca_align.md)                                                                                                         | Generalized PCA alignment via shared subspace estimation          |
| **GPA**                | [`generalized_procrustes()`](https://bbuchsbaum.github.io/manifoldalign/reference/generalized_procrustes.md)                                                                                 | Generalized Procrustes Analysis for partial observations          |
| **CONE**               | [`cone_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/cone_align.md), [`cone_align_multiple()`](https://bbuchsbaum.github.io/manifoldalign/reference/cone_align_multiple.md) | Spectral graph alignment with anchor-based matching               |
| **GRASP**              | [`grasp()`](https://bbuchsbaum.github.io/manifoldalign/reference/grasp.md), [`grasp_multiset()`](https://bbuchsbaum.github.io/manifoldalign/reference/grasp_multiset.md)                     | Graph Spectral alignment using structural descriptors             |
| **PARROT**             | [`parrot()`](https://bbuchsbaum.github.io/manifoldalign/reference/parrot.md)                                                                                                                 | Position-Aware Regularized OT alignment                           |
| **FPGW**               | [`fpgw()`](https://bbuchsbaum.github.io/manifoldalign/reference/fpgw.md)                                                                                                                     | Fused Partial Gromov-Wasserstein transport                        |
| **Gromov-Wasserstein** | [`gromov_wasserstein()`](https://bbuchsbaum.github.io/manifoldalign/reference/gromov_wasserstein.md)                                                                                         | Structure-preserving optimal transport                            |
| **Low-Rank**           | [`lowrank_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/lowrank_align.md)                                                                                                   | Memory-efficient alignment via low-rank factorization             |
| **MMA**                | [`mma_align_multiple()`](https://bbuchsbaum.github.io/manifoldalign/reference/mma_align_multiple.md)                                                                                         | Multiset Manifold Alignment for \>2 domains                       |
| **Coupled Diag.**      | [`coupled_diagonalization()`](https://bbuchsbaum.github.io/manifoldalign/reference/coupled_diagonalization.md)                                                                               | Joint diagonalization across modalities                           |
| **Spectral MNN**       | [`spectral_mnn_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/spectral_mnn_align.md)                                                                                         | Mutual nearest neighbor spectral alignment                        |
| **Token-OT Graph**     | [`token_ot_graph_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/token_ot_graph_align.md)                                                                                     | Token-level OT with graph structure                               |
| **SSMA**               | [`ssma_align()`](https://bbuchsbaum.github.io/manifoldalign/reference/ssma_align.md)                                                                                                         | Semi-Supervised Manifold Alignment                                |

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

# Align with KEMA
result <- kema(hd, ncomp = 5)

# Or align with Generalized Procrustes
result <- generalized_procrustes(hd, ncomp = 5)
```

## Unified Aligner Interface

All pairwise alignment methods share a consistent adapter interface:

``` r
# Create an aligner
aligner <- parrot_aligner(ncomp = 5)

# Fit a pair of domains
fit <- fit_pair(aligner, X1, X2, labels1, labels2)

# Align many domains at once
fit_many(aligner, list(X1, X2, X3), list(l1, l2, l3))
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
