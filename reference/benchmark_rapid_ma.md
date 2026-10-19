# Compare RAPID-MA with related package methods on identical splits

Compare RAPID-MA with related package methods on identical splits

## Usage

``` r
benchmark_rapid_ma(
  fixture = rapid_ma_benchmark_fixture(),
  methods = c("target_centroid", "rapid_ma", "rapid_ma_native", "rapid_ma_semisup",
    "lema", "ssma", "kema", "multiscale", "spectral_mnn", "token_ot"),
  ncomp = 8L,
  seed = 1L,
  controls = list(),
  keep_fits = FALSE
)
```

## Arguments

- fixture:

  A list from \[rapid_ma_benchmark_fixture()\], or a compatible list
  containing \`data\`, \`truth\`, \`labels\`, optional \`positions\` and
  \`attributes\`, and optional source-to-target \`correspondence\`.

- methods:

  Methods to run. Supported values are \`target_centroid\`,
  \`rapid_ma\`, \`rapid_ma_native\`, \`rapid_ma_semisup\`, \`lema\`,
  \`ssma\`, \`kema\`, \`multiscale\`, \`spectral_mnn\`, and
  \`token_ot\`. Embedding methods use one shared nearest-centroid head;
  explicitly named RAPID native and semi-supervised rows retain their
  own prediction heads.

- ncomp:

  Shared requested embedding width.

- seed:

  Deterministic method seed.

- controls:

  Optional named per-method control lists.

- keep_fits:

  Retain fitted objects in the result.

## Value

A \`rapid_ma_benchmark\` list containing a same-split results table,
failures, split/provenance metadata, and optionally fits.
