# Assess SSMA Alignment Quality on Synthetic Paired Domains

Runs \`ssma_align()\` on synthetic paired-domain problems and reports
diagnostics that separate embedding quality from decoding quality.

## Usage

``` r
assess_ssma(
  sizes = c(50L, 100L),
  d = 3L,
  noise_sd = 0.05,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 1,
  n_anchors = 8L,
  holdout_fraction = 0.25,
  n_reps = 3L,
  solvers = c("reduced", "operator"),
  serial = c(FALSE, TRUE),
  ssma_knn = 10L,
  ssma_rank_per_domain = 64L,
  ssma_lag_exclusion = 1L,
  ncomp = NULL,
  decode_k = 5L,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- sizes:

  Integer vector of node counts to assess.

- d:

  Integer feature dimension for synthetic node features.

- noise_sd:

  Numeric noise standard deviation added to target features.

- structure:

  Synthetic geometry for base coordinates: \`"ring"\`, \`"grid"\`,
  \`"random"\`, or \`"community"\`.

- permute_fraction:

  Fraction of nodes to permute (in \`(0,1\]\`).

- n_anchors:

  Number of synthetic anchor pairs to sample.

- holdout_fraction:

  Fraction of anchors reserved for holdout evaluation.

- n_reps:

  Integer number of replications per size.

- solvers:

  Character vector of SSMA solvers to assess: \`"reduced"\` and/or
  \`"operator"\`.

- serial:

  Logical vector indicating whether to run serial decontamination
  variants (\`FALSE\`, \`TRUE\`, or both).

- ssma_knn:

  Integer k used for SSMA within-domain graph construction.

- ssma_rank_per_domain:

  Integer rank cap for SSMA per-domain reduction.

- ssma_lag_exclusion:

  Integer lag exclusion for hard serial masking when \`serial = TRUE\`.

- ncomp:

  Optional integer embedding dimension. If \`NULL\`, uses \`min(10, n -
  2)\`.

- decode_k:

  Integer \`k\` used for top-k retrieval metrics.

- seed:

  Integer base random seed.

- verbose:

  Logical; print per-run progress.

## Value

A data.frame with per-run SSMA diagnostics including: runtime, embedding
finite-ness, raw/procrustes decoding metrics, holdout anchor recovery,
eigen diagnostics, serial diagnostics, and cross-solver subspace
distance when both solvers are run.
