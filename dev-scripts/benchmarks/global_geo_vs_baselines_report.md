# Global Geometry Baseline Benchmark Report

## Scope
Benchmark run comparing:
- `global_geo`
- `kema`
- `lra`
- `gpca`
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
| global_geo | condition | 0.0590 | 0.0506 | 0.7424 | 0.8667 |
| gpca | condition | 0.0750 | 0.2861 | 0.8175 | 0.9750 |
| kema | condition | 0.1337 | 0.5936 | 0.1522 | 0.8875 |
| lra | condition | 0.1760 | 0.0500 | 0.6227 | 0.9958 |
| gpca | id | 0.0213 | 0.1341 | 0.9672 | 0.9500 |
| global_geo | id | 0.0283 | 0.0506 | 0.7424 | 0.8667 |
| kema | id | 0.1357 | 0.0622 | 0.1363 | 0.6500 |
| lra | id | 0.1757 | 0.0001 | 0.5445 | 0.9375 |
| cone_multi | none | 0.0297 | 0.0359 | 0.4345 | 0.9417 |
| grasp_multi | none | 2.8693 | 0.4125 | 0.0314 | 0.6042 |

## Winners By Metric

### `id` supervision
- Fastest: `gpca` (0.021 s)
- Best identity alignment (lowest MSE): `lra` (0.0001)
- Best geometry preservation (Spearman): `gpca` (0.9672)
- Best condition kNN accuracy: `gpca` (0.9500)

### `condition` supervision
- Fastest: `global_geo` (0.059 s)
- Best identity alignment (lowest MSE): `lra` (0.0500)
- Best geometry preservation (Spearman): `gpca` (0.8175)
- Best condition kNN accuracy: `lra` (0.9958)

### `none` supervision
- Best unsupervised geometry: `cone_multi` (0.4345) vs `grasp_multi` (0.0314).

## Interpretation
- `gpca` remains strongest on geometry rank-correlation in supervised settings on this benchmark.
- `lra` remains strongest on identity MSE in `id` mode and strongest condition separability in `condition` mode.
- `global_geo` stays competitive and balanced with strong geometry and low runtime.
- `kema` remains stable but trails on geometry for this setup.
- `grasp_multi` improved versus earlier runs but is still slower and less accurate than other baselines here.

## Notes
- Metric values are identical across repetitions for each method/config because benchmark seeds are fixed; runtime is the primary varying quantity.
- Direct cross-method comparisons should consider objective mismatch across algorithms.
