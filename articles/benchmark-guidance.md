# Benchmark Guidance by Problem Regime

This package now tracks benchmark results in a regime-aware format
rather than a single pooled leaderboard. That matters because the right
method depends on the problem structure:

- are rows matched exactly or only weakly supervised?
- is overlap full or partial?
- are labels clean, noisy, or absent?
- do all domains share the same feature dimension?

This article reads the tracked synthetic benchmark outputs and turns
them into practical method guidance. The bundled tables are frozen
historical outputs; current
[`alignment_method_registry()`](https://bbuchsbaum.github.io/manifoldalign/reference/alignment_method_registry.md)
lifecycle states take precedence over historical rankings.

## Load benchmark outputs

``` r

bench_dir <- system.file("extdata", "benchmark-guidance", package = "manifoldalign",
                         mustWork = TRUE)

rank_path <- file.path(bench_dir, "toy_feature_aligners_regime_ranks.csv")
seed_path <- file.path(bench_dir, "toy_feature_seed_audit_summary.csv")

stopifnot(file.exists(rank_path), file.exists(seed_path))

regime_ranks <- utils::read.csv(rank_path, stringsAsFactors = FALSE)
seed_audit <- utils::read.csv(seed_path, stringsAsFactors = FALSE)
full_registry <- manifoldalign:::synthetic_alignment_scenarios(profile = "full")

required_rank_cols <- c(
  "scenario", "supervision_regime", "dimensionality_constraint",
  "method", "method_family", "overall_rank",
  "runtime_sec_mean", "max_abs_latent_cor_mean", "centroid_acc_mean"
)
required_seed_cols <- c(
  "scenario", "method", "method_family", "supervision_regime",
  "centroid_acc_mean", "centroid_acc_cv", "max_abs_latent_cor_mean"
)

stopifnot(all(required_rank_cols %in% names(regime_ranks)))
stopifnot(all(required_seed_cols %in% names(seed_audit)))
stopifnot(nrow(regime_ranks) > 0, nrow(seed_audit) > 0)
```

## Coverage

The checked-in fast benchmark and the checked-in seed-audit summary do
not have to cover exactly the same scenario set. The fast benchmark is
where we broaden regime coverage; the seed audit is where we focus on
stability diagnostics.

``` r

rank_scenarios <- sort(unique(regime_ranks$scenario))
seed_scenarios <- sort(unique(seed_audit$scenario))
full_scenarios <- sort(unique(full_registry$scenario))
full_only_scenarios <- setdiff(full_scenarios, rank_scenarios)

coverage <- data.frame(
  source = c(
    rep("fast_benchmark", length(rank_scenarios)),
    rep("seed_audit", length(seed_scenarios)),
    rep("full_registry_only", length(full_only_scenarios))
  ),
  scenario = c(rank_scenarios, seed_scenarios, full_only_scenarios),
  stringsAsFactors = FALSE
)

stopifnot(all(seed_scenarios %in% rank_scenarios | seed_scenarios %in% full_scenarios))
coverage
#>                source              scenario
#> 1      fast_benchmark     hard_nonisometric
#> 2      fast_benchmark       isometric_curve
#> 3      fast_benchmark         linear_affine
#> 4      fast_benchmark     noisy_label_shift
#> 5      fast_benchmark      partial_open_set
#> 6          seed_audit     hard_nonisometric
#> 7          seed_audit       isometric_curve
#> 8          seed_audit         linear_affine
#> 9  full_registry_only        highdim_stress
#> 10 full_registry_only many_domain_consensus
```

At the moment the checked-in fast benchmark covers the five lightweight
feature regimes:

- `linear_affine`
- `isometric_curve`
- `hard_nonisometric`
- `partial_open_set`
- `noisy_label_shift`

The full registry also defines `many_domain_consensus` and
`highdim_stress`, but those are kept in the `full` profile because they
behave more like stress tests than routine dev-loop benchmarks.

## Regime leaders

For each scenario and supervision regime we keep the top-ranked method
after combining runtime, latent correlation, and centroid recovery.

``` r

leader_rows <- !duplicated(regime_ranks[c(
  "scenario", "supervision_regime", "dimensionality_constraint"
)])

leaders <- regime_ranks[leader_rows, c(
  "scenario", "supervision_regime", "dimensionality_constraint",
  "method", "method_family", "overall_rank",
  "runtime_sec_mean", "max_abs_latent_cor_mean", "centroid_acc_mean"
)]

leaders
#>             scenario      supervision_regime dimensionality_constraint
#> 1      linear_affine    exact_correspondence              unrestricted
#> 10     linear_affine       label_supervision              unrestricted
#> 15     linear_affine   unsupervised_multiset              unrestricted
#> 16     linear_affine transport_side_geometry       shared_geometry_dim
#> 17     linear_affine    exact_correspondence         equal_feature_dim
#> 18     linear_affine   unsupervised_pairwise              unrestricted
#> 20     linear_affine          partial_anchor              unrestricted
#> 21   isometric_curve    exact_correspondence              unrestricted
#> 30   isometric_curve       label_supervision              unrestricted
#> 35   isometric_curve   unsupervised_multiset              unrestricted
#> 36   isometric_curve   unsupervised_pairwise              unrestricted
#> 38   isometric_curve          partial_anchor              unrestricted
#> 39 hard_nonisometric    exact_correspondence              unrestricted
#> 48 hard_nonisometric       label_supervision              unrestricted
#> 53 hard_nonisometric   unsupervised_multiset              unrestricted
#> 54 hard_nonisometric transport_side_geometry       shared_geometry_dim
#> 55 hard_nonisometric   unsupervised_pairwise              unrestricted
#> 56  partial_open_set    exact_correspondence              unrestricted
#> 63  partial_open_set   unsupervised_multiset              unrestricted
#> 64  partial_open_set transport_side_geometry       shared_geometry_dim
#> 65  partial_open_set   unsupervised_pairwise              unrestricted
#> 67  partial_open_set          partial_anchor              unrestricted
#> 68 noisy_label_shift    exact_correspondence              unrestricted
#> 77 noisy_label_shift       label_supervision              unrestricted
#> 82 noisy_label_shift   unsupervised_multiset              unrestricted
#> 83 noisy_label_shift transport_side_geometry       shared_geometry_dim
#> 84 noisy_label_shift    exact_correspondence         equal_feature_dim
#> 85 noisy_label_shift   unsupervised_pairwise              unrestricted
#> 87 noisy_label_shift          partial_anchor              unrestricted
#>                        method          method_family overall_rank
#> 1      global_geo_rowid_exact             global_geo     2.666667
#> 10          kema_linear_label                   kema     1.666667
#> 15 mma_consensus_unsupervised                    mma     1.000000
#> 16 multiset_uot_side_geometry           multiset_uot     1.000000
#> 17        genproc_rowid_exact generalized_procrustes     1.000000
#> 18    cone_align_unsupervised             cone_align     1.000000
#> 20      parrot_partial_anchor                 parrot     1.000000
#> 21     global_geo_rowid_exact             global_geo     1.500000
#> 30          kema_linear_label                   kema     1.666667
#> 35 mma_consensus_unsupervised                    mma     1.000000
#> 36    cone_align_unsupervised             cone_align     1.000000
#> 38      parrot_partial_anchor                 parrot     1.000000
#> 39     global_geo_rowid_exact             global_geo     1.500000
#> 48              lowrank_label                lowrank     1.833333
#> 53 mma_consensus_unsupervised                    mma     1.000000
#> 54 multiset_uot_side_geometry           multiset_uot     1.000000
#> 55         grasp_unsupervised                  grasp     1.000000
#> 56         spectral_mnn_exact           spectral_mnn     1.333333
#> 63 mma_consensus_unsupervised                    mma     1.000000
#> 64 multiset_uot_side_geometry           multiset_uot     1.000000
#> 65    cone_align_unsupervised             cone_align     1.166667
#> 67      parrot_partial_anchor                 parrot     1.000000
#> 68     global_geo_rowid_exact             global_geo     1.333333
#> 77          kema_linear_label                   kema     1.666667
#> 82 mma_consensus_unsupervised                    mma     1.000000
#> 83 multiset_uot_side_geometry           multiset_uot     1.000000
#> 84        genproc_rowid_exact generalized_procrustes     1.000000
#> 85    cone_align_unsupervised             cone_align     1.000000
#> 87      parrot_partial_anchor                 parrot     1.000000
#>    runtime_sec_mean max_abs_latent_cor_mean centroid_acc_mean
#> 1       0.070666667               0.9852159         0.9833333
#> 10      0.023666667               0.9662072         0.9833333
#> 15      0.087333333               0.9222068         0.9833333
#> 16      0.054666667               0.9679378         0.9833333
#> 17      0.045000000               0.9628587         0.9833333
#> 18      0.029000000               0.8989416         1.0000000
#> 20      0.164666667               0.7957893         0.9666667
#> 21      0.007333333               0.9910047         1.0000000
#> 30      0.020333333               0.9921191         1.0000000
#> 35      0.015333333               0.9242365         0.9866667
#> 36      0.010000000               0.7658405         1.0000000
#> 38      0.114333333               0.5081632         0.8266667
#> 39      0.007666667               0.9595887         1.0000000
#> 48      0.017000000               0.8769228         1.0000000
#> 53      0.014000000               0.8191251         0.9600000
#> 54      0.003000000               0.9961049         1.0000000
#> 55      0.015000000               0.4058587         1.0000000
#> 56      0.014666667               0.8977594         0.9305556
#> 63      0.016333333               0.8587978         0.8888889
#> 64      0.004000000               0.9922013         0.9305556
#> 65      0.010000000               0.8173933         1.0000000
#> 67      0.117000000               0.8914209         0.9861111
#> 68      0.007333333               0.9867950         0.9722222
#> 77      0.020666667               0.9894533         0.8888889
#> 82      0.015333333               0.9161527         0.9722222
#> 83      0.004000000               0.9770365         0.9722222
#> 84      0.006666667               0.9639595         0.9722222
#> 85      0.010666667               0.8941699         0.9861111
#> 87      0.119666667               0.8246086         0.9861111
```

Two patterns are immediate:

- `global_geo_rowid_exact` is now the strongest unrestricted
  exact-correspondence method on the full-overlap scenarios in the
  checked-in fast benchmark: `linear_affine`, `isometric_curve`,
  `hard_nonisometric`, and `noisy_label_shift`.
- `spectral_mnn_exact` still matters because it is the current leader on
  the partial/open-set exact-correspondence regime.
- In label-supervised regimes there is still no single winner:
  `kema_linear_label` leads on `linear_affine` and `isometric_curve`,
  `lowrank_label` leads on `hard_nonisometric`, and `gpca_label` leads
  on `noisy_label_shift`.
- The new multiset regimes are intentionally separate from the exact and
  label tables: `mma_consensus_unsupervised` leads the
  unsupervised-multiset blocks, and `multiset_uot_side_geometry` leads
  the transport-side-geometry blocks.

## Stability matters

Mean score alone is not enough. The seed audit tracks how much a method
moves across repeated random draws.

``` r

seed_audit <- seed_audit[order(
  -seed_audit$centroid_acc_cv,
  -seed_audit$max_abs_latent_cor_mean,
  seed_audit$method
), ]

most_variable <- seed_audit[seq_len(min(10L, nrow(seed_audit))), c(
  "scenario", "method", "method_family", "supervision_regime",
  "centroid_acc_mean", "centroid_acc_cv",
  "max_abs_latent_cor_mean"
)]

stopifnot(all(is.finite(most_variable$centroid_acc_cv)))
most_variable
#>             scenario                         method method_family
#> 7      linear_affine multiscale_spectral_supervised    multiscale
#> 13   isometric_curve              kema_linear_label          kema
#> 4      linear_affine            kema_rbf_full_label          kema
#> 5      linear_affine         kema_rbf_reduced_label          kema
#> 26 hard_nonisometric     multiscale_spectral_energy    multiscale
#> 14   isometric_curve            kema_rbf_full_label          kema
#> 27 hard_nonisometric multiscale_spectral_supervised    multiscale
#> 21 hard_nonisometric                     ssma_exact          ssma
#> 8      linear_affine     multiscale_spectral_hybrid    multiscale
#> 28 hard_nonisometric     multiscale_spectral_hybrid    multiscale
#>      supervision_regime centroid_acc_mean centroid_acc_cv
#> 7  exact_correspondence         0.7166667       0.2827454
#> 13    label_supervision         0.6366667       0.2764117
#> 4     label_supervision         0.7500000       0.2605356
#> 5     label_supervision         0.7833333       0.2572461
#> 26 exact_correspondence         0.7433333       0.2345128
#> 14    label_supervision         0.7033333       0.1976934
#> 27 exact_correspondence         0.7266667       0.1904436
#> 21 exact_correspondence         0.7000000       0.1748575
#> 8  exact_correspondence         0.8041667       0.1707059
#> 28 exact_correspondence         0.7366667       0.1660869
#>    max_abs_latent_cor_mean
#> 7                0.5865069
#> 13               0.2951984
#> 4                0.4916877
#> 5                0.5863469
#> 26               0.6253236
#> 14               0.2973345
#> 27               0.6015314
#> 21               0.9987103
#> 8                0.6292505
#> 28               0.6304174
```

These tracked runs show:

- `gpca_label` is the safest label-supervised baseline in the audited
  core trio.
- `kema_linear_label` is fast and often strong on latent correlation,
  but its centroid recovery is materially more seed-sensitive on the
  nonlinear curve case.
- The multiscale family is still clearly more variable than the best
  exact and label-supervised baselines.
- The checked-in seed audit is intentionally narrower than the fast
  benchmark, so it should be read as a stability diagnostic for the core
  trio rather than as a full coverage table for every benchmarked
  regime.

## Practical guidance

The current toy benchmark supports these rules of thumb:

- Start with `global_geo_rowid_exact` when exact row correspondences are
  available and the problem has full overlap with shared global
  geometry.
- Reach for `spectral_mnn_exact` earlier when overlap is
  partial/open-set, because that is where it overtakes the other
  exact-correspondence methods in the checked-in fast benchmark.
- Use `gpca_label` as the safest label-supervised baseline and
  `lowrank_label` as the strongest current label baseline on the harder
  nonlinear regime.
- Use `kema_linear_label` when you care about latent geometry and speed,
  but validate it with repeated seeds or resamples rather than trusting
  a single run.
- Treat `mma_consensus_unsupervised` as the current multiset baseline
  when you want a 3+ domain method without explicit side information.
- Treat `multiset_uot_side_geometry` as a different regime entirely: it
  is the right comparison point when a shared geometry/template is part
  of the problem specification.
- Treat generalized Procrustes as a useful exact-ID baseline only when
  feature dimensions match across domains.
- Treat coupled diagonalization and the multiscale methods as active
  research families rather than the current default recommendation on
  these synthetic tasks.
- Use `cone_align` as the unsupervised pairwise baseline when no
  correspondences are available between two domains. It leads the
  `unsupervised_pairwise` regime on most scenarios via iterative
  Procrustes on spectral embeddings.
- Use `grasp` when the unsupervised pairwise problem involves hard
  nonlinear distortion, where its diffusion descriptors provide an
  advantage.
- Use `parrot` when you have partial anchor correspondences (a subset of
  matched rows). It is the only method currently benchmarked in the
  `partial_anchor` regime.
- See `vignette(“method-selection-guide”)` for the full decision
  flowchart across all method families.

## What to benchmark next

The benchmark registry now separates lightweight fast scenarios from
heavier stress scenarios. The next high-value comparisons are:

- check in full-profile summaries for `many_domain_consensus` and
  `highdim_stress`,
- add truly adversarial side-information settings rather than only noisy
  labels,
- integrate `fpgw` and `token_ot_graph_align` into the benchmark (both
  return `multiblock_biprojector` objects but need input format
  adapters),
- add `gromov_wasserstein` with a transport-plan quality metric or
  barycentric projection step,
- track runtime and failure rates explicitly as first-class outputs.

Those settings are exactly where pooled “best method” tables tend to be
wrong. Use regime-specific summaries instead.
