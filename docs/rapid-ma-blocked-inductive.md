# RAPID-MA real-scene blocked-inductive validation

Status: **accuracy gate failed** on the frozen held-out geographic blocks.

## Bottom line

On nine previously untouched Chikusei blocks, feature-only RAPID-MA reached
57.88% macro-F1 and 49.94% mIoU. LeMA reached 56.18% macro-F1 and 51.97%
mIoU. RAPID-MA therefore led LeMA by 1.69 macro-F1 points but trailed it by
2.04 mIoU points. The preregistered target of at least three absolute points
on macro-F1 or mIoU was not met.

RAPID-MA is useful but not yet the claimed accuracy replacement:

- it improved over the direct target-only classifier by 3.61 macro-F1 and
  3.06 mIoU points;
- it was 1.79x faster than LeMA and retained 4.24x less fitted state;
- it did not dominate LeMA in OA, mIoU, rare-class recall, or calibration;
- KEMA completed no held-out fold under the common 30-second budget;
- fixed spatial and attribute smoothing substantially hurt untouched-region
  generalization, despite helping the development blocks.

These are real-scene results, but not independent-sensor results. The source is
real Chikusei hyperspectral reflectance and the target is a frozen ten-band
Gaussian approximation to Sentinel-2 VNIR responses generated from the same
scene.

## Frozen held-out results

Metrics are unweighted means across the nine held-out geographic blocks. The
blocks vary greatly in size and class composition, so standard deviations are
large and are retained in the summary artifact.

| Method | Completed | OA | Macro-F1 | mIoU | Rare recall | ECE | Runtime (s) | Retained fit (MB) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| RAPID-MA, feature OOS | 9/9 | 52.63 | **57.88** | 49.94 | 41.36 | 31.13 | 1.146 | 2.48 |
| LeMA | 9/9 | **58.05** | 56.18 | **51.97** | **59.64** | **26.87** | 2.049 | 10.50 |
| Target-only centroid | 9/9 | 51.61 | 54.27 | 46.88 | 36.23 | 39.27 | **0.096** | 0.00 |
| RAPID-MA, structured OOS | 9/9 | 41.77 | 45.96 | 39.07 | 34.72 | 27.02 | 1.392 | 2.48 |
| Target-only, structured OOS | 9/9 | 40.60 | 41.97 | 33.61 | 21.41 | 33.38 | 0.542 | 0.01 |
| SSMA | 9/9 | 29.11 | 34.61 | 32.01 | 28.06 | 42.18 | 0.459 | 0.16 |
| KEMA | 0/9 | NA | NA | NA | NA | NA | >30 each | NA |

RAPID-MA beat LeMA on macro-F1 in three of nine folds and on mIoU in two of
nine. The paired mean RAPID-minus-LeMA differences and 95% t intervals were:

| Metric | Mean difference | 95% interval |
|---|---:|---:|
| OA | -5.42 points | [-34.42, 23.57] |
| Macro-F1 | +1.69 points | [-19.74, 23.13] |
| mIoU | -2.04 points | [-22.84, 18.77] |

The intervals are descriptive with only nine heterogeneous blocks; they do
not support a superiority claim. Among test rows whose classes appeared in
training, RAPID-MA obtained 62.49% macro-F1 and 53.85% mIoU versus LeMA's
61.26% and 56.56%, respectively. The conclusion is unchanged.

## Leakage and budget controls

The benchmark uses a 4 by 4 equal-extent geographic grid. Development used
blocks 1, 2, 4, 5, and 6. The structured OOS weights were then frozen, the
implementation was committed as `d291f314fa9ea8a850404ee1aa963e10310796c8`,
and final evaluation used untouched blocks 7, 9, 10, 11, 12, 13, 14, 15, and
16.

For every held-out fold:

- the test block and a separate calibration block were absent from fitting,
  centering, scaling, relation construction, prototype construction, and OOS
  neighbor indexes;
- at most 750 fit rows were retained;
- every method received the same four labeled examples per available class,
  the same split, an eight-dimensional embedding budget where applicable, and
  a 30-second isolated-process limit;
- predictions for both calibration and test rows were generated before either
  label vector was accessed;
- calibration used only the separate calibration block; test labels were read
  only for final scoring;
- failures and timeouts remained in the raw table.

The target-only method is a nearest-centroid classifier on target features.
LeMA uses its learned low-modality linear projection. SSMA and KEMA use bounded
target-feature interpolation of their fitted embeddings. RAPID-MA uses bounded
target-feature interpolation of its learned relation-teacher probabilities.

## Alignment versus spatial smoothing

The diagnostic is a two-by-two comparison: target-only versus RAPID's
source-informed relation teacher, each with and without the same bounded OOS
structure smoother. The smoother combines target feature, position, and
attribute neighborhoods with fixed weights 0.50, 0.25, and 0.25.

On the five development blocks, the smoother raised RAPID-MA from 50.56% to
56.72% macro-F1 and from 44.64% to 48.25% mIoU. On the untouched blocks, it
reduced RAPID-MA by 11.91 macro-F1 and 10.87 mIoU points. The same smoother
reduced the target-only method by 12.30 macro-F1 and 13.27 mIoU points.

The similar damage with and without RAPID-MA identifies geographic smoothing
under region shift as the failure, rather than an OOS encoder bug. Position is
valuable within a scene, but a fixed nearest-neighbor mixture extrapolates
poorly across excluded tiles. Feature-only OOS remains the safe default. A
future structured predictor needs a training-only reliability gate that can
reject position and attribute views under geographic shift; development-block
gain is not sufficient evidence.

The feature-only RAPID-minus-target differences were +1.02 OA, +3.61
macro-F1, and +3.06 mIoU points. This is evidence that the source-informed
RAPID pipeline adds signal on average, but the paired 95% intervals include
zero and the comparison does not isolate manifold alignment from all other
RAPID teacher components.

## Calibration

Temperature was fit on the separate calibration block and never changes the
predicted class, so OA, macro-F1, and mIoU are identical before and after
calibration.

- RAPID-MA calibration-block NLL fell from 1.625 to 1.403. Held-out ECE fell
  only from 0.324 to 0.311 and improved in six of nine folds.
- LeMA calibration-block NLL fell from 2.525 to 1.181. Held-out ECE fell from
  0.511 to 0.269 and improved in seven of nine folds.
- Mean held-out ECE worsened for the target-only and SSMA models even though
  calibration-block NLL improved, demonstrating calibration-region shift.

The current RAPID probabilities are not well calibrated. A single global
temperature is too weak for spatially varying class priors and can transfer
poorly from one tile to another.

## Dataset and scope

The scene is the Chikusei airborne hyperspectral dataset distributed by Naoto
Yokoya under CC BY 3.0. The benchmark removes the bands marked noisy by the
source header, retains real reflectance and ground truth, and constructs the
target modality with fixed spectral response parameters. The downloaded
archive contained 1,066,530,806 bytes and had SHA-256:

```text
192825c63e764a16e9324470e005992f343c24317178e8a5546147aadde838dc
```

Official data page: <https://naotoyokoya.com/Download.html>

The result does not establish cross-scene transfer, a genuinely independent
sensor pair, open-set recognition, or broad remote-sensing superiority. Some
held-out blocks contain few classes; block 16 contains 166 labeled pixels from
one class. The raw artifact exposes each block rather than hiding this
heterogeneity in the mean.

## Reproduction

Download and acknowledge the official data terms:

```sh
RAPID_MA_ACCEPT_CHIKUSEI_LICENSE=yes \
  Rscript inst/benchmarks/download_rapid_ma_scene.R /path/to/data
```

Run the frozen held-out blocks:

```sh
RAPID_MA_PROFILE=exhaustive \
RAPID_MA_BLOCKS=7,9,10,11,12,13,14,15,16 \
RAPID_MA_MAX_FIT_ROWS=750 \
RAPID_MA_METHOD_TIMEOUT=30 \
Rscript inst/benchmarks/rapid_ma_blocked_inductive.R \
  /path/to/data/Chikusei_ENVI docs/rapid-ma-blocked-inductive
```

Artifacts:

- `rapid-ma-blocked-inductive-raw.csv`: every fold, method, metric, timeout,
  calibration result, runtime, and sampled RSS;
- `rapid-ma-blocked-inductive-summary.csv`: fold-level means and standard
  deviations;
- `rapid-ma-blocked-inductive-metadata.txt`: split, budget, dataset, platform,
  commit, and source hashes.

SHA-256:

```text
ebf34b378c0465d2077ae322432ee6d2f61d8b53158e2cae6cd0d42d4065c5cb  rapid-ma-blocked-inductive-raw.csv
4a9bcd3746e6de13e5816ab42da225824f9cf4ad3d0439024f556510ddc08885  rapid-ma-blocked-inductive-summary.csv
16bf06b7af41e70378b64751177d2b56abda8ac17164b94bb246387291008aa6  rapid-ma-blocked-inductive-metadata.txt
```

## Decision

Keep the new structure-aware OOS API because it correctly preserves spatial
and attribute information and is useful for explicit, gated experiments. Keep
feature-only interpolation as the default. Do not advertise the fixed
structured smoother or claim the accuracy target has been achieved.

The next accuracy improvement should operate inside the learned alignment,
not as unconditional post-fit geography: train a reliability-gated,
class-conditional OOS mixture using only fit-region cross-validation, include
an explicit target-only fallback, and evaluate it on a new scene or sensor pair
without reopening these held-out blocks.
