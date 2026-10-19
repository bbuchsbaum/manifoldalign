# Bounded multi-view matching for a RAPID-MA fit

Bounded multi-view matching for a RAPID-MA fit

## Usage

``` r
rapid_ma_match(
  fit,
  from = 1L,
  to = 2L,
  anchors = NULL,
  prototype_per_node = 3L,
  candidate_cap = 32L,
  prototype_bucket_cap = 128L,
  candidate_views = c("position", "attribute", "latent"),
  view_candidate_k = 12L,
  assignment = c("one_to_one", "independent"),
  latent_weight = 1,
  structure_weight = 0.25,
  position_weight = 1,
  attribute_weight = 0.15,
  prototype_weight = 0.25
)
```

## Arguments

- fit:

  A fitted \`rapid_ma\` object.

- from, to:

  Source and target domains.

- anchors:

  Optional one-to-one anchor pairs or paired anchor-ID vectors.

- prototype_per_node:

  Number of strongest prototype buckets queried.

- candidate_cap:

  Maximum retained target candidates per source node.

- prototype_bucket_cap:

  Maximum target nodes retained per prototype.

- candidate_views:

  Retained comparable views that contribute bounded nearest-neighbour
  candidates in addition to prototype buckets.

- view_candidate_k:

  Number of candidates requested from each comparable view.

- assignment:

  Either one-to-one greedy assignment or independent retrieval. Pair
  transforms use the one-to-one default.

- latent_weight, structure_weight, position_weight, attribute_weight,
  prototype_weight:

  Nonnegative reranking weights.

## Value

A \`rapid_ma_matching\` object with a sparse assignment operator, match
table, coverage, confidence, and unmatched-node diagnostics.
