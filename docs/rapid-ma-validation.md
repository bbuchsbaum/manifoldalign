# RAPID-MA accuracy validation

> **Regenerated 2026-10-19** with the documented runners and frozen seeds
> (3001--3020 main, 4001--4020 ablation) at commit
> `a6e1cd63b031b5c1f074813a48c5ca07d6479cb8`, after the evaluation fixes of
> the earlier status note (PCA-based target-only comparator and teacher seed,
> leave-one-label-out relation-safety folds without label leakage,
> macro-F1/mIoU over classes present in truth or predictions, undirected
> diffusion without reverse channels, `timeout` status), the relation-safety
> gate fix (candidates are compared with a no-transfer reference, the
> domain-local nearest centroid), deterministic `irlba` starts, and the
> exact RAPID-MA speedups. All 400 main and 200 ablation rows completed; no
> row failed, timed out, or was unsupported. A further 200 ablation rows on
> seeds 8001--8020, chosen after the gate was finalised and never inspected
> before this run, are reported as untouched evidence for the gate. Runtime
> and memory come from a 4-core x86_64 Linux machine (R 4.5.3, OpenBLAS) run
> with no other jobs, and are not comparable with the superseded macOS run.
> The runners were executed outside the package root, so their metadata
> files lack source hashes; `rapid-ma-regeneration-provenance.txt` records
> the commit, source md5s, package versions, CPU, and environment.
>
> Seed disclosure: while designing the gate fix, the frozen seeds were
> inspected, and its one-row non-inferiority tolerance was prompted by
> frozen seeds 3013 and 3014. The frozen-seed gate results are therefore not
> untouched evidence for that design; the 8001--8020 rows are.
>
> Claims that changed after regeneration: KEMA and Multiscale complete every
> run; the attribute and propagation ablation effects are smaller than first
> reported; the small-n memory advantage over LeMA is not reproduced on this
> platform; and the safety gate, after its fix, matches the target-only
> baseline under misleading structure instead of falling below it.

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
calibrated, and Hits@1 remains low. Under misleading structure the
relation-safety gate falls back to the no-transfer predictor and matches the
target-only baseline; it does not improve on it.

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
| RAPID-MA validated teacher | 20 | **0.9886** | **0.9851** | **0.9719** | 0.9969 | 0.3829 | 1.035 |
| KEMA | 20 | 0.9871 | 0.9815 | 0.9655 | 0.9828 | 0.5213 | 0.627 |
| Target-only centroid | 20 | 0.9795 | 0.9706 | 0.9454 | 0.9885 | 0.4020 | 0.011 |
| Dense LeMA reference | 20 | 0.9239 | 0.9012 | 0.8317 | 0.9594 | 0.3344 | 0.619 |
| SSMA | 20 | 0.9129 | 0.9057 | 0.8356 | 0.9353 | 0.4660 | 0.293 |
| Multiscale | 20 | 0.9049 | 0.8963 | 0.8217 | 0.9232 | 0.4566 | 0.497 |
| RAPID-MA embedding + centroid | 20 | 0.8269 | 0.8006 | 0.6973 | 0.8663 | 0.2287 | 0.756 |
| Spectral MNN | 20 | 0.7163 | 0.7131 | 0.5940 | 0.8410 | 0.2976 | 0.299 |
| Token-OT Graph | 20 | 0.6758 | 0.6905 | 0.5827 | 0.8168 | 0.2687 | 1.183 |
| RAPID-MA native UOT head | 20 | 0.5803 | 0.5114 | 0.4182 | **0.9780** | **0.1820** | 0.785 |

### Degraded-target scenario

| Method | Runs | OA | Macro-F1 | mIoU | Rare recall | ECE | Runtime (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| RAPID-MA validated teacher | 20 | **0.9689** | **0.9631** | **0.9316** | **0.9850** | 0.4281 | 1.034 |
| Target-only centroid | 20 | 0.8909 | 0.8608 | 0.7683 | 0.8818 | 0.3784 | 0.010 |
| KEMA | 20 | 0.8898 | 0.8457 | 0.7539 | 0.8167 | 0.4376 | 0.608 |
| Multiscale | 20 | 0.8220 | 0.7980 | 0.6736 | 0.8226 | 0.3768 | 0.476 |
| SSMA | 20 | 0.8091 | 0.7816 | 0.6513 | 0.7816 | 0.3759 | 0.287 |
| RAPID-MA embedding + centroid | 20 | 0.7477 | 0.7076 | 0.5836 | 0.7557 | 0.1943 | 0.754 |
| Dense LeMA reference | 20 | 0.6943 | 0.6573 | 0.5149 | 0.7052 | **0.1522** | 0.637 |
| Spectral MNN | 20 | 0.6511 | 0.6301 | 0.4794 | 0.7924 | 0.2443 | 0.293 |
| Token-OT Graph | 20 | 0.6394 | 0.6121 | 0.4834 | 0.6923 | 0.2414 | 1.175 |
| RAPID-MA native UOT head | 20 | 0.6144 | 0.5254 | 0.4270 | 0.9636 | 0.1790 | 0.722 |

Multiscale changed most from the run before the evaluation fixes (structured
OA 0.573 to 0.905, rare recall 0 to 0.92): the old rows were degenerate
single-class predictions from a run in which half its fits failed. Its rows
are now reproducible: its truncated SVDs previously drew their starting
vectors from the caller's random stream, which moved the structured OA by
about 0.002 between runs with the same seeds. The native UOT head's macro-F1
and mIoU fell (0.713 to 0.511 structured), consistent with the corrected
class set that charges absent-class predictions as false positives; it is
not a headline method.

## Paired improvements

Against LeMA on the same 20 seeds, RAPID-MA reduced classification error by
85.1% in the structured scenario and 89.8% after target degradation. Its
macro-F1 gains were 8.39 and 30.58 percentage points; mIoU gains were 14.02 and
41.67 points.

Against the target-only centroid, RAPID-MA reduced error by 44.4% and 71.5%.
The structured gains were 0.91 OA, 1.46 macro-F1, and 2.65 mIoU points (OA
wins/ties/losses 12/4/4); the degraded-target gains were 7.80, 10.23, and 16.33
points (19/1/0). The target-only comparator uses principal components of all
standardized features.

Against KEMA, which completes all 40 runs, RAPID-MA was +0.15 OA, +0.36
macro-F1, and +0.64 mIoU points in the structured scenario and +7.92, +11.74,
and +17.76 points after degradation.

## Structural and head ablation

The ablation uses degraded-target seeds. The frozen seeds 4001 through 4020
are shown with the untouched seeds 8001 through 8020 (see the seed
disclosure above).

| Variant | OA (4001--4020) | Macro-F1 | mIoU | OA (8001--8020) | Macro-F1 | mIoU | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---|
| Full validated RAPID-MA | 0.9553 | 0.9481 | 0.9070 | 0.9511 | 0.9434 | 0.9005 | Frozen complete pipeline |
| Domain-local seed only | 0.9557 | 0.9491 | 0.9079 | 0.9466 | 0.9380 | 0.8899 | Shared seed contributes little |
| No attributes | 0.9549 | 0.9436 | 0.9015 | 0.9553 | 0.9424 | 0.9004 | Attributes add almost nothing |
| Aligned seed only | 0.9220 | 0.9183 | 0.8619 | 0.9004 | 0.8948 | 0.8267 | Alignment evidence alone is insufficient |
| No relation propagation | 0.9038 | 0.8840 | 0.8026 | 0.8811 | 0.8603 | 0.7686 | Propagation adds 5 to 7 OA points |
| No position | 0.8936 | 0.8780 | 0.7889 | 0.8701 | 0.8591 | 0.7643 | Position adds 6 to 8 OA points |
| Target-only centroid | 0.8818 | 0.8514 | 0.7530 | 0.8405 | 0.8100 | 0.6993 | No RAPID structure |
| Misleading target structure | 0.8807 | 0.8486 | 0.7498 | 0.8424 | 0.8076 | 0.6966 | Gate falls back to no transfer |
| Feature relations only | 0.8689 | 0.8362 | 0.7330 | 0.8348 | 0.8085 | 0.6988 | Spatial/attribute channels matter |
| Dense LeMA | 0.6655 | 0.6311 | 0.4837 | 0.6652 | 0.6327 | 0.4859 | Dense reference |

The relation-safety gate compares each candidate relation set with a
no-transfer reference, the domain-local nearest centroid, using out-of-fold
balanced accuracy at the observed labels. The full relation set is kept when
it is within one held-out row of the reference; other candidates must beat
the comparator by 1/6 without a worse log-loss. Before this fix the gate
compared candidates only with the full relation set, never scored a
predictor equivalent to target-only, and fell 1.82 OA points below the
target-only baseline under misleading structure.

Under misleading target structure the gate now rejects transfer on 18 of the
frozen seeds and keeps the full set on two; the guarded result is equal to
target-only on 19 seeds and below it on one (mean -0.11 OA points). On the
untouched seeds it rejects transfer on 11, keeps the full set on 7 and picks
an alternative on 2; the guarded result is equal to target-only on 12 seeds,
above on 4 and below on 4 (mean +0.19 OA points). The gate therefore avoids
the cost of misleading structure but does not turn it into a gain, and
leave-one-label-out validation with four labels per class cannot resolve
one-row differences. Its fold scores are somewhat optimistic because the
shared-embedding part of each fold is not refitted without the held-out
label.

The ablation also limits the modality-transfer claim: full RAPID-MA is within
half an OA point of the domain-local-seed variant on both seed sets. In this
fixture, the classification gain comes from learned structural propagation
and positions, not cross-domain embedding transfer. Improving the core
alignment is still necessary.

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

At this small size, RAPID-MA is slower than dense LeMA (1.035 versus 0.619
seconds in the structured scenario). On this platform its mean endpoint RSS
increase is larger (66.1 versus 21.6 MB; 16.6 versus 6.7 MB after
degradation), and its sampled process peak is 317.0 versus 273.0 MB. The
superseded run on a different platform reported the opposite ordering; this
regeneration does not establish whether the platform, allocator, or code
accounts for the difference, so no small-n memory advantage is claimed.
Absolute RSS includes inherited R/package state. Large-n speed claims are
evaluated separately in `rapid-ma-scaling.md`.

Every method completed 20/20 runs in both scenarios. The slowest single run
took 1.9 seconds against the five-second limit.

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
regenerated scaling benchmark shows near-linear runtime (slope 0.98) and
retained state through 220,000 total nodes on a 4-core machine; the
million-node attempt on the previous machine was runtime-bound, not
memory-bound.

## Artifacts

- `rapid-ma-validation-{raw,summary}.csv` and
  `rapid-ma-validation-metadata.txt`: main comparison.
- `rapid-ma-ablation-{raw,summary}.csv` and
  `rapid-ma-ablation-metadata.txt`: component and misleading-structure study.
- `rapid-ma-ablation-holdout-{raw,summary}.csv` and
  `rapid-ma-ablation-holdout-metadata.txt`: the same ablation on untouched
  seeds 8001--8020.
- `rapid-ma-scaling.md`: complete-pipeline scaling comparison and limits.
- `rapid-ma-blocked-inductive.md`: real-scene held-out accuracy and calibration.
- `rapid-ma-matching.md`: frozen exact-node matching comparison.
- `rapid-ma-regeneration-provenance.txt`: commit, source md5s, package
  versions, CPU, and environment of the 2026-10-19 regeneration.

SHA-256:

```text
8b5224916968595ca3ee76c82099f8d5c9c126143f009192ba9e02b3f0a944cd  rapid-ma-validation-raw.csv
c064eb9490a7fb743bac31fbcb7eb6fa9eb4319d44396268905fb3b7fc1087ed  rapid-ma-validation-summary.csv
ebb1dbf4f219dc033868d726cc6c8ee0df1e30a4247167790ad0572c8ff6ed68  rapid-ma-validation-metadata.txt
22331c0640331ffe68a88e680592ebda2f0f3bc6a66aca6619df30bf9515aacf  rapid-ma-ablation-raw.csv
33ac1f193b16840d49b7685892c111962cd2d8aada112800ebfc3aea91573e36  rapid-ma-ablation-summary.csv
41b5825af5d69ccec94b2a4eea60b9520d90790195a4c031bc98765d9eb45033  rapid-ma-ablation-metadata.txt
46c721155fe44d46c597430d9f5b65b774858a0ec2ceedad11b8fd64109c3066  rapid-ma-ablation-holdout-raw.csv
27882224a9cd4c304bf05adc2376c05b121ff8390b5b26454a2106cca43e93e6  rapid-ma-ablation-holdout-summary.csv
8cbf07cb2687009fd5ab6578fd626b54138825391f90e39ab2ed3682b4bacfb0  rapid-ma-ablation-holdout-metadata.txt
5daae21f9402b7c594278cb1fb1bd14a513950816718075436f3426c917c5a96  rapid-ma-regeneration-provenance.txt
```
