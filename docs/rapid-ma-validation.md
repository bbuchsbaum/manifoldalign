# RAPID-MA accuracy validation

> **Regenerated 2026-10-19** with the documented runners and frozen seeds
> (3001--3020 main, 4001--4020 ablation) at commit
> `78d33c9b34f0185b120920ae18ae5a150e17401a`, after the evaluation fixes
> listed in the previous status note: PCA-based target-only comparator and
> teacher seed, leave-one-label-out relation-safety folds without label
> leakage, macro-F1/mIoU over classes present in truth or predictions,
> undirected diffusion without reverse channels, and `timeout` status for
> isolated-runner timeouts. All 400 main and 200 ablation rows completed; no
> row failed, timed out, or was unsupported. Runtime and memory were measured
> on a different machine (4-core x86_64 Linux, R 4.5.3, OpenBLAS) from the
> superseded run and are not comparable with it. The runners were executed
> outside the package root, so their metadata files lack source hashes;
> `rapid-ma-regeneration-provenance.txt` records the commit, source md5s,
> package versions, CPU, and environment.
>
> Claims that no longer hold after regeneration: KEMA and Multiscale now
> complete every run; under misleading target structure the guarded result is
> *below* the target-only baseline (it was reported 1.33 points above); the
> attribute and propagation ablation effects are smaller; and the memory
> advantage over LeMA is not reproduced on this platform.

Status: synthetic transductive and structural-robustness evidence complete;
real-scene and blocked-inductive evidence not yet established.

## Bottom line

The complete RAPID-MA classification pipeline is highly accurate in this
controlled benchmark. Across 20 untouched seeds, its validated sparse relation
teacher reached 98.86% overall accuracy (OA) on the structured scenario and
96.89% after target-feature degradation. The independent dense LeMA reference
reached 92.39% and 69.43%, respectively.

KEMA is effectively tied on the easy structured case: over all 20 seeds
RAPID-MA was 0.15 OA, 0.36 macro-F1 and 0.64 mIoU points ahead (OA wins, ties
and losses 7/7/6). KEMA now completes every degraded-target run, and RAPID-MA
leads it there by 7.92 OA, 11.74 macro-F1 and 17.76 mIoU points (19 wins, one
tie).

The strongest conclusion is therefore robustness, not universal dominance.
RAPID-MA is substantially better than LeMA, SSMA, KEMA, and the target-only
baseline when target features degrade, and it remains competitive with KEMA
when the task is easy. This result applies to classification, not exact node
matching. The core embedding alone is much weaker, probabilities are poorly
calibrated, Hits@1 remains low, and the relation-safety gate does not protect
against misleading structure (see the ablation).

## Protocol

- Main evaluation seeds: 3001 through 3020. Ablation seeds: 4001 through
  4020. All development ranges are disjoint and recorded in metadata.
- `structured`: imbalanced forest/urban/water classes, heterogeneous feature
  widths, spatial positions, mixed attributes with missing values, known
  partial correspondence, and 24 target-only nodes.
- `degraded_target`: the same fixture plus deterministic Gaussian target
  feature noise with standard deviation 0.75; positions and attributes remain.
- Four paired observed labels per class; all other target labels are held out.
- Embedding methods use the same nearest-centroid classifier. Native-head rows
  are named explicitly.
- RAPID-MA uses a frozen 75% domain-local / 25% shared-embedding centroid seed,
  two sparse relation-propagation hops, and no pseudo-label rounds.
- Misleading-relation safety is selected only from out-of-fold probabilities at
  observed labels. The full relation set is retained unless a bounded none,
  single, or pair candidate improves balanced validation accuracy by at least
  1/6 (two balanced label decisions in this fixture).
- Main methods run in isolated processes with a five-second limit. Failures
  remain visible. Peak RSS is sampled every 20 ms.
- Held-out labels are not read during fitting, relation selection, or tuning.

## Main classification results

Values are successful-run means. Runs are out of 20. Lower ECE is better.

### Structured scenario

| Method | Runs | OA | Macro-F1 | mIoU | Rare recall | ECE | Runtime (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| RAPID-MA validated teacher | 20 | **0.9886** | **0.9851** | **0.9719** | 0.9969 | 0.3829 | 1.297 |
| KEMA | 20 | 0.9871 | 0.9815 | 0.9655 | 0.9828 | 0.5213 | 0.565 |
| Target-only centroid | 20 | 0.9795 | 0.9706 | 0.9454 | 0.9885 | 0.4020 | 0.007 |
| Dense LeMA reference | 20 | 0.9239 | 0.9012 | 0.8317 | 0.9594 | 0.3344 | 0.520 |
| SSMA | 20 | 0.9129 | 0.9057 | 0.8356 | 0.9353 | 0.4660 | 0.269 |
| Multiscale | 20 | 0.9064 | 0.8972 | 0.8233 | 0.9203 | 0.4581 | 0.412 |
| RAPID-MA embedding + centroid | 20 | 0.8269 | 0.8006 | 0.6973 | 0.8663 | 0.2287 | 1.104 |
| Spectral MNN | 20 | 0.7163 | 0.7131 | 0.5940 | 0.8410 | 0.2976 | 0.258 |
| Token-OT Graph | 20 | 0.6758 | 0.6905 | 0.5827 | 0.8168 | 0.2687 | 0.982 |
| RAPID-MA native UOT head | 20 | 0.5803 | 0.5114 | 0.4182 | **0.9780** | **0.1820** | 1.101 |

### Degraded-target scenario

| Method | Runs | OA | Macro-F1 | mIoU | Rare recall | ECE | Runtime (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| RAPID-MA validated teacher | 20 | **0.9689** | **0.9631** | **0.9316** | **0.9850** | 0.4281 | 1.333 |
| Target-only centroid | 20 | 0.8909 | 0.8608 | 0.7683 | 0.8818 | 0.3784 | 0.008 |
| KEMA | 20 | 0.8898 | 0.8457 | 0.7539 | 0.8167 | 0.4376 | 0.568 |
| Multiscale | 20 | 0.8220 | 0.7980 | 0.6736 | 0.8226 | 0.3768 | 0.425 |
| SSMA | 20 | 0.8091 | 0.7816 | 0.6513 | 0.7816 | 0.3759 | 0.251 |
| RAPID-MA embedding + centroid | 20 | 0.7477 | 0.7076 | 0.5836 | 0.7557 | 0.1943 | 1.097 |
| Dense LeMA reference | 20 | 0.6943 | 0.6573 | 0.5149 | 0.7052 | **0.1522** | 0.559 |
| Spectral MNN | 20 | 0.6511 | 0.6301 | 0.4794 | 0.7924 | 0.2443 | 0.257 |
| Token-OT Graph | 20 | 0.6394 | 0.6121 | 0.4834 | 0.6923 | 0.2414 | 1.058 |
| RAPID-MA native UOT head | 20 | 0.6144 | 0.5254 | 0.4270 | 0.9636 | 0.1790 | 1.114 |

Multiscale changed most from the superseded run (structured OA 0.573 to
0.906, rare recall 0 to 0.92): the old rows were degenerate single-class
predictions from a run in which half its fits failed. The native UOT head's
macro-F1 and mIoU fell (0.713 to 0.511 structured), consistent with the
corrected class set that now charges absent-class predictions as false
positives; it is not a headline method.

## Paired improvements

Against LeMA on the same 20 seeds, RAPID-MA reduced classification error by
85.1% in the structured scenario and 89.8% after target degradation. Its
macro-F1 gains were 8.39 and 30.58 percentage points; mIoU gains were 14.02 and
41.67 points.

Against the target-only centroid, RAPID-MA reduced error by 44.4% and 71.5%.
The structured gains were 0.91 OA, 1.46 macro-F1, and 2.65 mIoU points (OA
wins/ties/losses 12/4/4); the degraded-target gains were 7.80, 10.23, and 16.33
points (19/1/0). The target-only comparator is stronger than in the superseded
run because it now uses principal components of all standardized features.

Against KEMA, which now completes all 40 runs, RAPID-MA was +0.15 OA, +0.36
macro-F1, and +0.64 mIoU points in the structured scenario and +7.92, +11.74,
and +17.76 points after degradation.

## Structural and head ablation

The ablation uses 20 untouched degraded-target seeds (4001 through 4020).

| Variant | OA | Macro-F1 | mIoU | Interpretation |
|---|---:|---:|---:|---|
| Domain-local seed only | **0.9557** | **0.9491** | **0.9079** | Shared seed contributes nothing measurable here |
| Full validated RAPID-MA | 0.9553 | 0.9481 | 0.9070 | Frozen complete pipeline |
| No attributes | 0.9549 | 0.9436 | 0.9015 | Attributes add almost nothing (0.04 OA) |
| Aligned seed only | 0.9220 | 0.9165 | 0.8605 | Alignment evidence alone is insufficient |
| No relation propagation | 0.9038 | 0.8840 | 0.8026 | Propagation adds about 5 OA points |
| No position | 0.8936 | 0.8780 | 0.7889 | Position adds about 6 OA points |
| Target-only centroid | 0.8818 | 0.8514 | 0.7530 | No RAPID structure |
| Feature relations only | 0.8648 | 0.8313 | 0.7255 | Spatial/attribute channels matter |
| Misleading target structure | 0.8636 | 0.8189 | 0.7115 | Below target-only: the safety gate does not protect |
| Dense LeMA | 0.6655 | 0.6311 | 0.4837 | Dense reference |

For misleading target structure, the safety selector changed the relation set
on 17/20 seeds: feature-only on ten and no structural relations on seven. It
retained all relations on three. The guarded result is 1.82 OA and 3.25
macro-F1 points *below* the target-only baseline (below on 12 seeds, equal on
4, above on 4), although far above LeMA. The superseded run reported it 1.33
points above target-only; that run scored the gate's leave-one-label-out folds
while relation neighbours still carried the held-out label. With that leakage
removed, the gate avoids catastrophic use but does not make RAPID-MA safer than
ignoring the structure, so it should not be described as a safety guarantee.

The ablation also limits the modality-transfer claim: full RAPID-MA is 0.04 OA
and 0.10 macro-F1 points *below* the domain-local-seed variant. In this
fixture, the classification gain comes from learned structural propagation and
positions, not cross-domain embedding transfer. Improving the core alignment is
still necessary.

## Matching, memory, and reliability

RAPID-MA has 100% matching coverage, but exact retrieval is weak: structured
Hits@1 is 4.00% and bounded MRR@100 is 0.123. LeMA obtains 4.58% Hits@1 and
MRR@100 of 0.205 at only 10% coverage because its reference embeds paired
labeled source rows only. Under degradation, RAPID-MA Hits@1 is 3.17% versus
LeMA's 7.92% at 10% coverage. Classification success is not evidence of solved
node matching.

The follow-up bounded multi-view matcher fixes the largest candidate-recall
defect. On 20 untouched 300-source/360-target seeds it reaches 14.20% Hits@1
at 100% coverage, versus 3.00% for the prototype-only matcher and 1.35% for
RAPID embedding kNN. This is a substantial improvement, not a solved-matching
claim; see `rapid-ma-matching.md`.

At this small size, RAPID-MA is slower than dense LeMA (1.297 versus 0.520
seconds in the structured scenario). On this platform its mean endpoint RSS
increase is larger (105.6 versus 27.1 MB; 42.9 versus 24.3 MB after
degradation), and its sampled process peak is 339.7 versus 261.3 MB. The
superseded run on a different platform reported the opposite ordering; this
regeneration does not establish whether the platform, allocator, or code
accounts for the difference, so no small-n memory advantage is claimed.
Absolute RSS includes inherited R/package state. Large-n speed claims are
evaluated separately in `rapid-ma-scaling.md`.

Every method completed 20/20 runs in both scenarios. The slowest single run
took 1.6 seconds against the five-second limit.

## What is and is not established

Established:

- deterministic, leakage-checked synthetic classification accuracy;
- more than 80% relative error reduction over the dense LeMA oracle;
- resilience to target-feature degradation when spatial structure remains;
- bounded out-of-fold fallback under misleading structural channels that
  avoids catastrophic use but stays below the target-only baseline;
- strict blocked-inductive Chikusei evaluation with training-only calibration;
- a visible failed real-scene accuracy gate (+1.69 macro-F1 but -2.04 mIoU
  versus LeMA on untouched blocks); these two Chikusei results predate the
  evaluation fixes and have not been regenerated (see
  `rapid-ma-blocked-inductive.md`);
- bounded multi-view matcher improvement from 3.00% to 14.20% Hits@1;
- explicit calibration, matching, runtime, sampled peak RSS, and failures;
- frozen seeds, preprocessing, controls, source hashes, and raw records.

Not established:

- cross-scene or genuinely independent-sensor generalization;
- unseen-class/open-set accuracy;
- well-calibrated probabilities;
- high-accuracy exact node correspondence;
- million-node completion (the 1,000,002-node attempt timed out at 1,200
  seconds with 7.20 GB peak RSS) or a direct 10x comparison at 100,000 nodes
  against a surviving baseline.

The next accuracy work is a new cross-scene or independent-sensor benchmark
with a fit-only reliability gate for structural OOS views. The next performance
work is approximate/reusable relation indexing plus warm-started UOT. The
scaling benchmark establishes near-linear growth through 220,000 total nodes;
the measured million-node attempt is runtime-bound, not memory-bound.

## Artifacts

- `rapid-ma-validation-{raw,summary}.csv` and
  `rapid-ma-validation-metadata.txt`: main comparison.
- `rapid-ma-ablation-{raw,summary}.csv` and
  `rapid-ma-ablation-metadata.txt`: component and misleading-structure study.
- `rapid-ma-scaling.md`: complete-pipeline scaling comparison and limits.
- `rapid-ma-blocked-inductive.md`: real-scene held-out accuracy and calibration.
- `rapid-ma-matching.md`: frozen exact-node matching comparison.
- `rapid-ma-regeneration-provenance.txt`: commit, source md5s, package
  versions, CPU, and environment of the 2026-10-19 regeneration.

SHA-256:

```text
cd43302ff8b3622c2e3323ecb9f89c2542f2ede9995bd91e92c818df83573354  rapid-ma-validation-raw.csv
77531eb7950d049e552859300d332dcbabcc845ed8aeed21c1bc5447edb8f4d3  rapid-ma-validation-summary.csv
84e1a4de232fbbffa12d6d314ca1c3f1a85d700f7e48e14b076bbfe840f4809e  rapid-ma-validation-metadata.txt
f35c7d1cc8f4a319e044057169c925606933cb07854d39bce4a313c3d964c555  rapid-ma-ablation-raw.csv
24d0a972b655694561ade0feccde8252d5ed4789d81717ce0c8924f675617fcd  rapid-ma-ablation-summary.csv
5003ad04bbc3262e9aba2caba3bebc41cd355f53153b25b2167f66a8b991659d  rapid-ma-ablation-metadata.txt
8fcbca6d80391b28bd7e6bc9fa108bd65dbc25bb0fb4ad71716a1b14f1da73cf  rapid-ma-regeneration-provenance.txt
```
