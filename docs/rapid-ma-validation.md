# RAPID-MA accuracy validation

Status: synthetic transductive and structural-robustness evidence complete;
real-scene and blocked-inductive evidence not yet established.

## Bottom line

The complete RAPID-MA classification pipeline is highly accurate in this
controlled benchmark. Across 20 untouched seeds, its validated sparse relation
teacher reached 98.75% overall accuracy (OA) on the structured scenario and
96.21% after target-feature degradation. The independent dense LeMA reference
reached 92.39% and 69.20%, respectively.

KEMA is effectively tied on the easy structured case: on the 19 seeds where it
finished, its OA was 0.04 points higher, while RAPID-MA's macro-F1 and mIoU were
0.16 and 0.31 points higher. KEMA completed only 5/20 degraded-target runs;
RAPID-MA completed 20/20 and led by 6.21 OA points on those five paired runs.

The strongest conclusion is therefore robustness, not universal dominance.
RAPID-MA is substantially better than LeMA, SSMA, and the target-only baseline
when target features degrade, and it remains competitive with KEMA when the
task is easy. This result applies to classification, not exact node matching.
The core embedding alone is much weaker, probabilities are poorly calibrated,
and Hits@1 remains low.

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
| KEMA | 19 | **0.9904** | **0.9856** | **0.9724** | 0.9819 | 0.5194 | 0.182 |
| RAPID-MA validated teacher | 20 | 0.9875 | 0.9839 | 0.9696 | **0.9969** | 0.3838 | 0.272 |
| Target-only centroid | 20 | 0.9742 | 0.9628 | 0.9319 | 0.9708 | 0.3958 | 0.001 |
| Dense LeMA reference | 20 | 0.9239 | 0.9013 | 0.8319 | 0.9568 | 0.3339 | 0.139 |
| SSMA | 20 | 0.9163 | 0.9103 | 0.8430 | 0.9460 | 0.4690 | 0.057 |
| RAPID-MA embedding + centroid | 20 | 0.7955 | 0.7873 | 0.6805 | 0.9073 | **0.2046** | 0.254 |
| Spectral MNN | 20 | 0.6883 | 0.6968 | 0.5643 | 0.8488 | 0.2751 | 0.055 |
| Token-OT Graph | 20 | 0.6811 | 0.6661 | 0.5449 | 0.8087 | 0.2832 | 0.071 |
| RAPID-MA native UOT head | 20 | 0.6409 | 0.7131 | 0.5862 | **1.0000** | 0.2054 | 0.253 |
| Multiscale | 9 | 0.5732 | 0.2428 | 0.1911 | 0.0000 | 0.2399 | 0.099 |

### Degraded-target scenario

| Method | Runs | OA | Macro-F1 | mIoU | Rare recall | ECE | Runtime (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| RAPID-MA validated teacher | 20 | **0.9621** | **0.9566** | **0.9216** | **0.9911** | 0.4296 | 0.276 |
| KEMA | 5 | 0.8879 | 0.8566 | 0.7598 | 0.8655 | 0.4284 | 0.668 |
| Target-only centroid | 20 | 0.8553 | 0.8240 | 0.7147 | 0.8477 | 0.3483 | 0.001 |
| SSMA | 20 | 0.8072 | 0.7810 | 0.6506 | 0.7987 | 0.3740 | 0.058 |
| RAPID-MA embedding + centroid | 20 | 0.7614 | 0.7438 | 0.6133 | 0.8615 | 0.1727 | 0.256 |
| Dense LeMA reference | 20 | 0.6920 | 0.6543 | 0.5112 | 0.6995 | **0.1515** | 0.142 |
| Token-OT Graph | 20 | 0.6379 | 0.6273 | 0.4923 | 0.7852 | 0.2335 | 0.071 |
| Spectral MNN | 20 | 0.6129 | 0.6170 | 0.4858 | 0.7818 | 0.2399 | 0.055 |
| RAPID-MA native UOT head | 20 | 0.5917 | 0.6775 | 0.5529 | **0.9974** | 0.2402 | 0.259 |
| Multiscale | 10 | 0.5652 | 0.2407 | 0.1884 | 0.0000 | 0.2318 | 0.107 |

## Paired improvements

Against LeMA on the same 20 seeds, RAPID-MA reduced classification error by
83.6% in the structured scenario and 87.7% after target degradation. Its
macro-F1 gains were 8.26 and 30.24 percentage points; mIoU gains were 13.77 and
41.04 points.

Against the target-only centroid, RAPID-MA reduced error by 51.5% and 73.8%.
The structured gains were 1.33 OA, 2.10 macro-F1, and 3.77 mIoU points; the
degraded-target gains were 10.68, 13.26, and 20.69 points.

On KEMA's common successful seeds, RAPID-MA was -0.04 OA, +0.16 macro-F1, and
0.31 mIoU points ahead in the structured scenario. On the five common degraded
successes, RAPID-MA led by 6.21, 8.09, and 13.19 points. KEMA's 5/20 degraded
completion rate makes its successful-run mean selectively observed.

## Structural and head ablation

The ablation uses 20 untouched degraded-target seeds (4001 through 4020).

| Variant | OA | Macro-F1 | mIoU | Interpretation |
|---|---:|---:|---:|---|
| Full validated RAPID-MA | **0.9564** | **0.9498** | **0.9097** | Frozen complete pipeline |
| Domain-local seed only | 0.9561 | 0.9479 | 0.9065 | Shared seed contributes little here |
| No attributes | 0.9360 | 0.9152 | 0.8591 | Attributes add about 2 OA points |
| No position | 0.8860 | 0.8685 | 0.7764 | Position adds about 7 OA points |
| No relation propagation | 0.8667 | 0.8411 | 0.7385 | Propagation adds about 9 OA points |
| Misleading target structure | 0.8633 | 0.8283 | 0.7228 | Safety gate prevents catastrophic use |
| Aligned seed only | 0.8598 | 0.8556 | 0.7746 | Alignment evidence alone is insufficient |
| Feature relations only | 0.8595 | 0.8255 | 0.7199 | Spatial/attribute channels matter |
| Target-only centroid | 0.8500 | 0.8150 | 0.7010 | No RAPID structure |
| Dense LeMA | 0.6659 | 0.6314 | 0.4841 | Dense reference |

For misleading target structure, the safety selector changed the relation set
on 13/20 seeds: feature-only on 11, no propagation on one, and
feature-plus-spatial on one. It retained all relations on seven seeds. The
guarded result stays 1.33 OA and macro-F1 points above the target-only baseline,
but remains far below the full informative-structure result.

The ablation also limits the modality-transfer claim: full RAPID-MA exceeds the
domain-local-seed variant by only 0.04 OA and 0.19 macro-F1 points. In this
fixture, nearly all classification gain comes from learned structural
propagation, not cross-domain embedding transfer. Improving the core alignment
is still necessary.

## Matching, memory, and reliability

RAPID-MA has 100% matching coverage, but exact retrieval is weak: structured
Hits@1 is 4.50% and bounded MRR@100 is 0.131. LeMA obtains 4.17% Hits@1 and
MRR@100 of 0.203 at
only 10% coverage because its reference embeds paired labeled source rows only.
Under degradation, RAPID-MA Hits@1 is 3.42% versus LeMA's 7.92% at 10%
coverage. Classification success is not evidence of solved node matching.

The follow-up bounded multi-view matcher fixes the largest candidate-recall
defect. On 20 untouched 300-source/360-target seeds it reaches 14.13% Hits@1
at 100% coverage, versus 3.50% for the prototype-only matcher and 1.57% for
RAPID embedding kNN. This is a substantial improvement, not a solved-matching
claim; see `rapid-ma-matching.md`.

At this small size, RAPID-MA is slower than dense LeMA (0.272 versus 0.139
seconds). Its mean endpoint RSS increase is smaller (51.1 versus 124.1 MB), and
its absolute sampled process peak is 406.0 versus 474.1 MB in the structured
scenario. Absolute RSS includes inherited R/package state. Large-n speed claims
are evaluated separately in `rapid-ma-scaling.md`.

KEMA timed out on 1/20 structured and 15/20 degraded runs. Multiscale completed
9/20 and 10/20; failures were timeouts or `no right-hand side in 'b'`. All
other methods completed 20/20.

## What is and is not established

Established:

- deterministic, leakage-checked synthetic classification accuracy;
- more than 80% relative error reduction over the dense LeMA oracle;
- resilience to target-feature degradation when spatial structure remains;
- bounded out-of-fold fallback under misleading structural channels;
- strict blocked-inductive Chikusei evaluation with training-only calibration;
- a visible failed real-scene accuracy gate (+1.69 macro-F1 but -2.04 mIoU
  versus LeMA on untouched blocks);
- bounded multi-view matcher improvement from 3.50% to 14.13% Hits@1;
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

SHA-256:

```text
eec627ddcc1b1f3551e56aaa3fd4ab3c54590a0357a72852f952d6f59235510e  rapid-ma-validation-raw.csv
01990d4e1bd307f368d3a698ce6252658faf5252880327cb9ade1d79054c8887  rapid-ma-validation-summary.csv
eb3509475d65d7c71ab60c3452413fa580d2539ee4a20c3bddb14d5da497a5c6  rapid-ma-validation-metadata.txt
29cd9a12a8a62202ab68951dc521a0e0bd951607417d63264c723bd6481c581b  rapid-ma-ablation-raw.csv
b77954ea2be08ed9e7f3c5e1fed123dceee6568adf79a149d6a59d91ebae3b1f  rapid-ma-ablation-summary.csv
b7eb2c0061212842653373f5658c1649ba203e08a641726ea492db7bd9417021  rapid-ma-ablation-metadata.txt
```
