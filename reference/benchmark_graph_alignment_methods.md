# Benchmark graph alignment methods on synthetic permuted graphs

Generates paired domains where the target is a permuted/noisy copy of
the source and compares multiple alignment methods available in this
package.

## Usage

``` r
benchmark_graph_alignment_methods(
  sizes = c(50L, 100L),
  d = 3L,
  noise_sd = 0.05,
  structure = c("ring", "grid", "random", "community"),
  permute_fraction = 1,
  n_anchors = 0L,
  methods = c("token_ot_graph", "fpgw", "cone_align", "grasp", "parrot", "ssma", "lra",
    "gpca"),
  decode_mode = c("native", "common_nn", "common_procrustes_lap"),
  ssma_procrustes = FALSE,
  candidate_k = 80L,
  coarsen_method = c("kmeans", "louvain"),
  token_mode = c("view_only", "view_plus_neighbors"),
  views = "raw",
  n_levels = 1L,
  prior_mode = c("none", "hard", "soft"),
  fpgw_omega1 = 0.5,
  fpgw_epsilon = 0.01,
  fpgw_max_iter = 50L,
  fpgw_inner_max_iter = 20L,
  fpgw_tol = 1e-06,
  ssma_solver = c("reduced", "operator"),
  ssma_knn = 10L,
  ssma_rank_per_domain = 64L,
  ssma_use_serial = FALSE,
  ssma_lag_exclusion = 1L,
  lra_mu = 0.5,
  lra_lambda = 0.01,
  lra_sv_thresh = 1,
  lra_solver = c("operator", "explicit"),
  gpca_u = 0.5,
  gpca_lambda = 0.1,
  n_reps = 3L,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- sizes:

  Integer vector of node counts to benchmark.

- d:

  Integer feature dimension for synthetic node features.

- noise_sd:

  Numeric noise standard deviation added to target features.

- structure:

  Character vector of synthetic geometries for base node coordinates:
  \`"ring"\`, \`"grid"\`, \`"random"\`, or \`"community"\`.

- permute_fraction:

  Fraction of nodes to permute (in \`(0,1\]\`).

- n_anchors:

  Integer number of anchor correspondences to add (0 allowed).

- methods:

  Character vector of methods to run. Options include:
  \`"token_ot_graph"\`, \`"cone_align"\`, \`"grasp"\`, \`"parrot"\`,
  \`"fpgw"\`, \`"ssma"\`, \`"lra"\`, \`"gpca"\`.

- decode_mode:

  Decoding mode for mapping source nodes to targets: \`"native"\` uses
  each method's native assignment/transport decoding; \`"common_nn"\`
  uses nearest-neighbor decoding on aligned scores for all methods;
  \`"common_procrustes_lap"\` uses anchor-Procrustes (when anchors are
  available) followed by LAP/Hungarian decoding on aligned scores.

- ssma_procrustes:

  Logical; when \`FALSE\` (default), SSMA ignores
  \`"common_procrustes_lap"\` and decodes in its own common space via
  nearest-neighbor (\`"common_nn"\`). Set \`TRUE\` to explicitly
  benchmark an \`SSMA+Procrustes\` decode variant.

- candidate_k:

  Integer number of candidates per source node for
  \`token_ot_graph_align()\`.

- coarsen_method:

  Coarsening method for multilevel Token-OT alignment: \`"kmeans"\` or
  \`"louvain"\`.

- token_mode:

  Tokenization mode passed to \[token_ot_graph_align_control()\].

- views:

  Views passed to \[token_ot_graph_align_control()\].

- n_levels:

  Integer number of multilevel hierarchy levels for
  \`token_ot_graph_align()\`.

- prior_mode:

  Prior lifting mode for multilevel (\`"none"\`, \`"hard"\`,
  \`"soft"\`).

- fpgw_omega1:

  Weight of the feature term for \`fpgw()\` (0–1).

- fpgw_epsilon:

  Entropic warm-start regularization for \`fpgw()\`.

- fpgw_max_iter:

  Maximum Frank-Wolfe iterations for \`fpgw()\`.

- fpgw_inner_max_iter:

  Inner iterations for \`fpgw()\` updates.

- fpgw_tol:

  Convergence tolerance for \`fpgw()\`.

- ssma_solver:

  Solver passed to \[ssma_align_control()\] for \`ssma_align()\`.

- ssma_knn:

  Integer k for within-domain graph construction in SSMA.

- ssma_rank_per_domain:

  Integer per-domain rank cap in SSMA reduction.

- ssma_use_serial:

  Logical; if \`TRUE\`, enable serial decontamination in SSMA.

- ssma_lag_exclusion:

  Non-negative lag exclusion radius for SSMA serial hard masking when
  \`ssma_use_serial = TRUE\`.

- lra_mu:

  Mixing parameter for \`lowrank_align()\` (\`mu\` argument).

- lra_lambda:

  Ridge parameter passed to \`lowrank_align()\`.

- lra_sv_thresh:

  Singular-value threshold passed to \`lowrank_align()\`.

- lra_solver:

  Solver backend passed to \`lowrank_align()\`.

- gpca_u:

  Within/between trade-off parameter passed to \`gpca_align()\`.

- gpca_lambda:

  Ridge regularization passed to \`gpca_align()\`.

- n_reps:

  Integer number of replications per size.

- seed:

  Integer seed for reproducibility.

- verbose:

  Logical; print per-run progress.

## Value

A data.frame with per-run runtime and accuracy metrics.

## Details

Benchmarks are intended for interactive profiling and method comparison
rather than automated CI.
