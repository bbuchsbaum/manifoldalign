# Global Geometry Baseline Benchmark Report

## Scope
Benchmark run comparing:
- `global_geo`
- `kema`
- `lra`
- `gpca`
- `uot`
- `spectral_mnn`
- `cone_multi`
- `grasp_multi`

Using the same synthetic multi-domain benchmark input and `ncomp = 8`, with supervision modes:
- `id`
- `condition`
- `none` (for methods that do not use label/correspondence supervision)

Source files:
- `dev-scripts/benchmarks/benchmark_global_geo_vs_kema_lra.R`
- `dev-scripts/benchmarks/global_geo_vs_baselines_results.csv`
- `dev-scripts/benchmarks/global_geo_vs_baselines_summary.csv`

## Mean Results By Method/Supervision

| method | supervision | runtime_mean (s) | identity_mse_mean | geometry_spearman_mean | condition_knn_acc_mean |
|---|---|---:|---:|---:|---:|
| spectral_mnn | condition | 0.0493 | 0.1286 | 0.5213 | 0.9500 |
| global_geo | condition | 0.0610 | 0.0506 | 0.7424 | 0.8667 |
| gpca | condition | 0.0760 | 0.2861 | 0.8175 | 0.9750 |
| kema | condition | 0.1300 | 0.5936 | 0.1522 | 0.8875 |
| lra | condition | 0.1763 | 0.0500 | 0.6227 | 0.9958 |
| gpca | id | 0.0227 | 0.1341 | 0.9672 | 0.9500 |
| global_geo | id | 0.0307 | 0.0506 | 0.7424 | 0.8667 |
| spectral_mnn | id | 0.0427 | 0.0289 | 0.5213 | 0.9417 |
| kema | id | 0.1360 | 0.0622 | 0.1363 | 0.6500 |
| lra | id | 0.1663 | 0.0001 | 0.5445 | 0.9375 |
| spectral_mnn | none | 0.0280 | 0.1891 | 0.5213 | 0.8417 |
| cone_multi | none | 0.0330 | 0.0359 | 0.4345 | 0.9417 |
| uot | none | 0.0350 | 0.3293 | 0.9719 | 0.9625 |
| grasp_multi | none | 2.8420 | 0.4125 | 0.0314 | 0.6042 |

## Winners By Metric

### `id` supervision
- Fastest: `gpca` (0.0227 s)
- Best identity alignment (lowest MSE): `lra` (0.0001)
- Best geometry preservation (Spearman): `gpca` (0.9672)
- Best condition kNN accuracy: `gpca` (0.9500)

### `condition` supervision
- Fastest: `spectral_mnn` (0.0493 s)
- Best identity alignment (lowest MSE): `lra` (0.0500)
- Best geometry preservation (Spearman): `gpca` (0.8175)
- Best condition kNN accuracy: `lra` (0.9958)

### `none` supervision
- Fastest: `spectral_mnn` (0.0280 s)
- Best identity alignment (lowest MSE): `cone_multi` (0.0359)
- Best geometry preservation (Spearman): `uot` (0.9719)
- Best condition kNN accuracy: `uot` (0.9625)

## Interpretation
- `gpca` remains strongest on geometry rank-correlation in supervised settings on this benchmark.
- `lra` remains strongest on identity MSE in `id` mode and strongest condition separability in `condition` mode.
- `global_geo` stays competitive and balanced with strong geometry and low runtime.
- `spectral_mnn` remains very fast, improves identity alignment in `id` mode vs `global_geo` on this benchmark, and is the strongest geometry-preserving method in `none` mode.
- `kema` remains stable but trails on geometry for this setup.
- `grasp_multi` improved versus earlier runs but is still slower and less accurate than other baselines here.
- `uot` (unsupervised, anchored to domain 1) is very strong on geometry rank-correlation and condition kNN here, but does not optimize identity correspondences in this benchmark setup.

## Notes
- Metric values are identical across repetitions for each method/config because benchmark seeds are fixed; runtime is the primary varying quantity.
- Direct cross-method comparisons should consider objective mismatch across algorithms.
