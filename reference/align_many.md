# Align multiple datasets by synchronizing pairwise transforms

This meta-aligner lifts any pairwise method to N domains by: 1) choosing
a sparse dataset graph, 2) fitting pairwise models on edges, 3)
extracting relative transforms, and 4) running group-appropriate
synchronization to get per-domain transforms into a shared space.

## Usage

``` r
align_many(
  domains,
  algo,
  graph = c("auto", "mst", "complete", "star"),
  consensus = c("auto", "sync", "star", "tree"),
  k = NULL,
  parallel = TRUE,
  sync = NULL,
  ...
)
```

## Arguments

- domains:

  list of domain objects (often matrices), one per dataset

- algo:

  an aligner object created by new_aligner() with adapter methods

- graph:

  "auto" \| "mst" \| "complete" \| "star"

- consensus:

  "auto" \| "sync" \| "star" \| "tree" (currently: "sync"/"star")

- k:

  target latent dimension (optional)

- parallel:

  logical; reserved for future parallel pair fits

- sync:

  list of synchronization options (graph, refine, gl, perm settings)

- ...:

  forwarded to adapter methods

## Value

list with transforms, embeddings (if apply_transform succeeds), dataset
graph, pairwise fits, and diagnostics

## Details

If the adapter reports native multi-view support, dispatches to
fit_many().

## Examples

``` r
# \donttest{
set.seed(1)
X1 <- matrix(rnorm(50), 10, 5)
X2 <- matrix(rnorm(50), 10, 5)
X3 <- matrix(rnorm(50), 10, 5)
domains <- list(X1, X2, X3)
# align_many requires an aligner object with adapter methods
# See package vignettes for complete examples
# }
```
