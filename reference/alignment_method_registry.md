# Canonical alignment method registry

Returns the machine-readable source of truth for public alignment
methods. Registry records distinguish supplied graphs from
feature-derived topology and metric-measure geometry; describe topology,
node-data, supervision, output, complexity, dependency, and scale
contracts; and embed lifecycle metadata validated by
\[validate_alignment_lifecycle()\].

## Usage

``` r
alignment_method_registry()

alignment_method_info(method, registry = alignment_method_registry())

alignment_method_table(registry = alignment_method_registry())

validate_alignment_method_registry(registry = alignment_method_registry())
```

## Arguments

- method:

  One canonical ID, exported entry point, adapter name, or benchmark
  alias.

- registry:

  A registry returned by \`alignment_method_registry()\`.

## Value

\`alignment_method_registry()\` returns a named list with class
\`alignment_method_registry\`. \`alignment_method_info()\` returns one
record. \`alignment_method_table()\` returns a data frame.
\`validate_alignment_method_registry()\` returns an
\`alignment_registry_problems\` data frame; zero rows means valid.

## Details

\`alignment_method_info()\` resolves a canonical method ID, exported
entry point, adapter name, or benchmark alias to exactly one record.
\`alignment_method_table()\` returns a flattened data frame suitable for
documentation and audits. \`validate_alignment_method_registry()\`
returns typed schema and contradiction problems.

## Examples

``` r
registry <- alignment_method_registry()
rapid <- alignment_method_info("rapid_ma", registry)
rapid$topology$supplied
#> [1] "native"
alignment_method_table(registry)[, c("id", "lifecycle", "primary_input")]
#>                                                  id        lifecycle
#> rapid_ma                                   rapid_ma     experimental
#> token_ot_graph_align           token_ot_graph_align     experimental
#> spectral_mnn_align               spectral_mnn_align     experimental
#> ssma_align                               ssma_align     experimental
#> global_geo_align                   global_geo_align legacy_reference
#> kema                                           kema     experimental
#> gpca_align                               gpca_align     experimental
#> lowrank_align                         lowrank_align     experimental
#> linear_sim_embed                   linear_sim_embed     experimental
#> generalized_procrustes       generalized_procrustes     experimental
#> grasp                                         grasp     experimental
#> grasp_multiset                       grasp_multiset     experimental
#> cone_align                               cone_align legacy_reference
#> cone_align_multiple             cone_align_multiple legacy_reference
#> parrot                                       parrot     experimental
#> gromov_wasserstein               gromov_wasserstein     experimental
#> fpgw                                           fpgw     experimental
#> coupled_diagonalization     coupled_diagonalization     experimental
#> mma_align_multiple               mma_align_multiple     experimental
#> multiscale_manifold_align multiscale_manifold_align       deprecated
#> multiset_uot_align               multiset_uot_align     experimental
#> ot_procrustes                         ot_procrustes     experimental
#>                                   primary_input
#> rapid_ma                         supplied_graph
#> token_ot_graph_align      feature_derived_graph
#> spectral_mnn_align        feature_derived_graph
#> ssma_align                feature_derived_graph
#> global_geo_align          feature_derived_graph
#> kema                      feature_derived_graph
#> gpca_align                     feature_manifold
#> lowrank_align                  feature_manifold
#> linear_sim_embed               feature_manifold
#> generalized_procrustes         feature_manifold
#> grasp                     feature_derived_graph
#> grasp_multiset            feature_derived_graph
#> cone_align                feature_derived_graph
#> cone_align_multiple       feature_derived_graph
#> parrot                    feature_derived_graph
#> gromov_wasserstein               metric_measure
#> fpgw                             metric_measure
#> coupled_diagonalization   feature_derived_graph
#> mma_align_multiple        feature_derived_graph
#> multiscale_manifold_align feature_derived_graph
#> multiset_uot_align               metric_measure
#> ot_procrustes                  feature_manifold
```
