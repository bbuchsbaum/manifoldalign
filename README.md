# manifoldalign

<!-- badges: start -->
[![R-CMD-check](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/bbuchsbaum/manifoldalign/actions/workflows/pkgdown.yaml/badge.svg)](https://bbuchsbaum.github.io/manifoldalign/)
<!-- badges: end -->

**manifoldalign** is an R package for aligning data from multiple domains into a shared latent space. It implements a broad family of manifold alignment and domain adaptation methods — from kernel-based and spectral approaches to optimal transport — preserving local manifold structure while matching samples across domains.

## Alignment methods

No method is currently classified as stable or recommended. The table below
states the input model instead of treating every method as an observed-graph
aligner. See the [complete capability inventory](https://bbuchsbaum.github.io/manifoldalign/METHOD_COMPARISON_TABLE.html)
and [lifecycle decision record](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/ALIGNMENT_METHOD_SUNSET_REVIEW.md).

| Method | Function | Input contract | Lifecycle |
|---|---|---|---|
| **RAPID-MA** | `rapid_ma()`, `rapid_ma_feature_scores()` | Supplied sparse weighted relations, positions, and mixed attributes | Experimental |
| **KEMA** | `kema()` | Supervised feature-derived kernel manifolds | Experimental |
| **GPCA** | `gpca_align()` | Supervised feature manifolds | Experimental |
| **GPA** | `generalized_procrustes()` | Orthogonal feature-task alignment with partial observations | Experimental |
| **CONE** | `cone_align()`, `cone_align_multiple()` | Spectral graphs derived from features | Legacy/reference |
| **GlobalGeo** | `global_geo_align()` | Correspondence-linked geometry derived from features | Legacy/reference |
| **GRASP** | `grasp()`, `grasp_multiset()` | Functional maps on feature-derived graphs | Experimental |
| **PARROT** | `parrot()` | Anchor-guided dense OT on feature-derived networks | Experimental |
| **FPGW** | `fpgw()` | Partial metric-measure OT | Experimental |
| **Gromov-Wasserstein** | `gromov_wasserstein()` | Balanced metric-measure OT | Experimental |
| **Low-Rank** | `lowrank_align()` | Supervised low-rank feature manifolds | Experimental |
| **MMA** | `mma_align_multiple()` | Multiset feature-derived manifolds | Experimental |
| **Coupled Diag.** | `coupled_diagonalization()` | Joint feature-derived spectral bases | Experimental |
| **Spectral MNN** | `spectral_mnn_align()` | Feature-derived spectral graphs | Experimental |
| **Token-OT Graph** | `token_ot_graph_align()` | Sparse OT on feature-derived graphs | Experimental |
| **SSMA** | `ssma_align()` | Semi-supervised feature-derived manifolds | Experimental |
| **Multiset UOT** | `multiset_uot_align()` | Unbalanced metric-measure/template OT | Experimental |
| **SUMCOR** | `sumcor_align()`, `calibrate_domain()` | Repeated observations of shared keyed events; linear per-domain encoders with frozen-reference calibration | Experimental |
| **LSE** | `linear_sim_embed()` | Single-domain target-similarity embedding | Experimental utility |
| **Multiscale** | `multiscale_manifold_align()` | Incoherent experimental hierarchy retained for reproduction | Deprecated; removal 1.0.0 |

## Current evidence

The frozen alignment tournament contains 588 typed rows: 42 attributed-graph
cases by 14 methods, with zero missing result rows and zero validation
problems. It records failures, unsupported tasks, runtime, peak RSS, controls,
paired effects, and checksums for every published receipt file. Resumable binary
checkpoints and the debug log remain local operational state. Read the
[checked interpretation](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/INTERPRETATION.md)
or inspect the [raw results](https://github.com/bbuchsbaum/manifoldalign/blob/main/docs/evidence/alignment-tournament-v2/raw_results.csv).

RAPID-MA completed 42/42 and is the only tournament method that consumed the
supplied weighted relations together with positions and mixed attributes. It
is not a universal accuracy winner: a position-only control beat its node
matching in every non-oracle regime. Its frozen downstream classification probe
did beat source-only transfer in every non-oracle regime. The current
full-scoring scale receipt reaches 22,600 total nodes; an older runner reached
220,000, while the one-million-node and blocked-inductive promotion gates
failed. Those limitations are part of the public evidence contract.

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

# Align with KEMA, an experimental supervised method (as_hyperdesign() stores the
# labels in a column named `label`)
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
set.seed(1)
n_events <- 40
latent <- matrix(rnorm(n_events * 3), n_events, 3)
make_domain <- function(p) {
  rows <- rep(seq_len(n_events), times = 2)          # two repeats per event
  X <- latent[rows, ] %*% matrix(rnorm(3 * p), 3, p) +
    matrix(rnorm(length(rows) * p, sd = 0.5), length(rows), p)
  list(X = X, design = data.frame(event_id = paste0("e", rows),
                                  repeat_id = rep(1:2, each = n_events)))
}
domains <- list(a = make_domain(20), b = make_domain(25), c = make_domain(18))

fit <- sumcor_align(
  lapply(domains, `[[`, "X"),
  design = lapply(domains, `[[`, "design"),
  event_id = "event_id", repeat_id = "repeat_id",
  ncomp = 3,
  control = sumcor_control(beta = 0, basis_rank = 10)
)

# Encode new observations of one domain
Z_new <- oos_predict(fit, domains$c$X[1:5, ], side = "c")

# Add a held-out domain against the frozen reference
heldout <- make_domain(22)
fit2 <- calibrate_domain(fit, heldout$X, design = heldout$design,
                         domain_name = "heldout")
```

## Evaluating transfer to a new domain

`cv_alignment_transfer()` holds out one domain at a time, fits the method on
the source domains only, calibrates the held-out domain against that frozen
fit from calibration rows that share no event or group with the evaluation
rows, and scores a source-trained readout of a target variable plus linear
reconstructions in each domain's own columns (`recon_*_r2`). Methods that
cannot calibrate a new domain are reported as `unsupported`, never refitted
with the held-out domain.

```r
X_list <- lapply(domains, `[[`, "X")
design_list <- lapply(domains, function(d) {
  transform(d$design, label = latent[as.integer(sub("e", "", d$design$event_id)), 1] > 0)
})
cv <- cv_alignment_transfer(
  X_list,
  methods = list(sumcor = function(data, design) {
    sumcor_align(data, design = design, event_id = "event_id",
                 repeat_id = "repeat_id", ncomp = 3,
                 control = sumcor_control(beta = 0, basis_rank = 10))
  }),
  design = design_list,
  event_id = "event_id",
  target = "label"
)
cv
```

## Dual (feature-side) coordinates

Embedding methods place samples in the shared latent space. For RAPID-MA, which
is not a linear projector and therefore has no loadings, `rapid_ma_feature_scores()`
places each feature column at the weighted barycentre of the samples it loads
on, so features from different domains can be compared in the same space:

```r
set.seed(1)
X <- list(domain1 = matrix(abs(rnorm(40 * 6)), 40, 6),
          domain2 = matrix(abs(rnorm(35 * 6)), 35, 6))
fit  <- rapid_ma(X, ncomp = 3L)
dual <- rapid_ma_feature_scores(fit, X)
dual$scores$domain1   # feature x component coordinates for one domain
dual$pooled           # equal-weight mean across domains
dual$dispersion       # how consistently each feature is placed across domains
```

## Unified aligner interface

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
