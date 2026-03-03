# Benchmark token_ot_graph_align scalability on synthetic permuted graphs

Generates paired domains where the target is a permuted/noisy copy of
the source and reports runtime and basic alignment accuracy as problem
size grows.

## Usage

``` r
benchmark_token_ot_graph_align_scalability(
  sizes = c(50L, 100L, 200L),
  d = 4L,
  noise_sd = 0.01,
  candidate_k = 80L,
  n_levels = 1L,
  prior_mode = c("none", "hard", "soft"),
  token_mode = c("view_only", "view_plus_neighbors"),
  views = "raw",
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

- candidate_k:

  Integer number of candidates per source node.

- n_levels:

  Integer number of multilevel hierarchy levels.

- prior_mode:

  Prior lifting mode for multilevel (\`"none"\`, \`"hard"\`,
  \`"soft"\`).

- token_mode:

  Tokenization mode passed to \[token_ot_graph_align_control()\].

- views:

  Views passed to \[token_ot_graph_align_control()\].

- n_reps:

  Integer number of replications per size.

- seed:

  Integer seed for reproducibility.

- verbose:

  Logical; print per-run progress.

## Value

A data.frame with per-run runtime and summary metrics.

## Details

Benchmarks are intended for interactive profiling rather than automated
CI.
