# Generate a reproducible structure-aware alignment benchmark fixture

Generate a reproducible structure-aware alignment benchmark fixture

## Usage

``` r
rapid_ma_benchmark_fixture(
  n_common = 120L,
  extra_target = 24L,
  source_features = 24L,
  target_features = 10L,
  labels_per_class = 4L,
  seed = 1L
)
```

## Arguments

- n_common:

  Number of corresponding nodes shared by the two domains.

- extra_target:

  Additional target-only nodes.

- source_features, target_features:

  Modality-specific feature widths.

- labels_per_class:

  Paired observed labels retained per class.

- seed:

  Deterministic seed.

## Value

A named fixture list accepted by \[benchmark_rapid_ma()\].
