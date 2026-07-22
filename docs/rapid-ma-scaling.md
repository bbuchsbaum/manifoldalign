# RAPID-MA complete-pipeline scaling

Status: 220,000-node synthetic run complete; 1,000,002-node completion gate
failed on the frozen 1,200-second runtime ceiling.

## Bottom line

RAPID-MA's complete validated-teacher pipeline scales approximately linearly
over the measured range. It fit 100,000 source plus 120,000 target nodes
(220,000 total) in 108.3 seconds, retained a 482.9 MB result, and reached a
5.07 GB sampled process peak. The log-log time slope from 2,500 through 100,000
common nodes is 1.083; the retained-object slope is 0.993.

The million-node smoke did not complete. A frozen run with 454,546 source and
545,456 target nodes (1,000,002 total) hit its 1,200-second wall-clock limit.
It reached only 7.20 GB sampled peak RSS under a 30 GB kill ceiling, so runtime,
not memory, is the demonstrated blocker. The correspondence-ranking evaluator
was disabled for this large run and cannot explain the timeout.

The speed claim needs qualification:

- RAPID-MA is slower than SSMA, Spectral-MNN, and Token-OT at small sizes.
- It crosses those sparse baselines around 2,500 to 5,000 common nodes and is
  the only tested method that finishes the 10,000-common-node task under the
  shared ten-second budget.
- It is 2.48x faster than dense LeMA at 500 common nodes and 4.74x faster at
  1,000. LeMA exceeds ten seconds from 2,500 onward.
- It is 8.33x faster than KEMA at 500; KEMA exceeds ten seconds at 1,000, where
  RAPID takes 0.816 seconds (a greater-than-12.3x timeout lower bound).
- No direct 10x runtime comparison was completed at 100,000 common nodes,
  because dense LeMA and KEMA are already outside the time/memory envelope.

## Protocol

- Same deterministic synthetic fixture and fixed 12-label-per-domain budget at
  every size.
- Common-node sizes: 250, 500, 1,000, 2,500, 5,000, and 10,000. The target has
  20% additional unmatched nodes, so total nodes are 2.2 times the common size.
- Three repetitions through 1,000; one repetition at larger sizes.
- RAPID-only extensions at 50,000 and 100,000 common nodes.
- Every task creates its fixture and runs in an isolated process. RSS is sampled
  every 20 ms. Common tasks have a ten-second cap; RAPID large tasks have a
  four-minute cap. Timeout rows remain in the data.
- R 4.5.1 on arm64 macOS (Darwin 23.3), 14 logical cores, 36 GB system memory.
- Hits@1 is exact. Matching MRR is the declared bounded MRR@100 lower bound;
  the evaluator does not materialize an all-pairs distance matrix.

## RAPID-MA scaling

| Common nodes | Total nodes | Runtime (s) | Retained result (MB) | Sampled peak RSS (MB) | OA | Macro-F1 |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 550 | 0.316 | 1.56 | 436.8 | 0.9861 | 0.9804 |
| 500 | 1,100 | 0.475 | 2.71 | 464.9 | 0.9722 | 0.9712 |
| 1,000 | 2,200 | 0.816 | 5.18 | 483.9 | 0.9921 | 0.9895 |
| 2,500 | 5,500 | 2.085 | 12.45 | 757.7 | 0.9552 | 0.9517 |
| 5,000 | 11,000 | 3.742 | 24.47 | 884.9 | 0.9644 | 0.9673 |
| 10,000 | 22,000 | 7.311 | 48.63 | 1,590.2 | 0.9949 | 0.9937 |
| 50,000 | 110,000 | 46.239 | 241.65 | 3,687.5 | 0.9892 | 0.9856 |
| 100,000 | 220,000 | 108.298 | 482.86 | 5,074.0 | 0.9363 | 0.9400 |
| 454,546 | 1,000,002 | >1,200 (timeout) | NA | 7,199.3 | NA | NA |

Absolute peak RSS includes inherited R and loaded-package state. Retained state
is the cleaner scaling signal and grows essentially linearly.

## Runtime comparison

| Method | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---:|---:|---:|---:|---:|---:|
| RAPID-MA teacher | 0.316 | 0.475 | 0.816 | 2.085 | **3.742** | **7.311** |
| Dense LeMA | **0.308** | 1.181 | 3.864 | timeout | timeout | timeout |
| KEMA | 0.665 | 3.957 | timeout | timeout | timeout | timeout |
| SSMA | **0.064** | **0.096** | **0.220** | **1.160** | 4.170 | timeout |
| Spectral MNN | 0.069 | 0.112 | 0.282 | 1.295 | 4.538 | timeout |
| Token-OT Graph | 0.093 | 0.156 | 0.435 | 2.186 | 5.503 | timeout |
| Multiscale | failed | failed | failed | failed | failed | timeout |

Times through 1,000 are three-run means; larger entries are single runs.
“Timeout” means the isolated end-to-end task exceeded ten seconds, not that the
method can never complete with a larger budget. Multiscale's smaller cases fail
with its existing `no right-hand side in 'b'` solver error.

## Dense-memory contrast

The LeMA reference's learned graph is dense over target nodes plus paired
source landmarks. With 100,000 common nodes and 20% target-only expansion, that
graph has about 120,012 rows. One double-precision matrix of that shape requires
about 115.2 GB; retaining both `W` and `L` alone requires about 230.4 GB before
distance workspaces and embeddings. That exceeds the 36 GB machine without
running LeMA at this size.

At 10,000 common nodes, LeMA had already reached 6.65 GB sampled RSS when the
ten-second timeout stopped it. RAPID completed in 7.31 seconds with a 1.59 GB
sampled peak. At 100,000, RAPID's sparse retained result is 482.9 MB.

## Scaling contracts and limitations

Automated tests double fixture size and check that retained state grows by less
than 3x, transport nonzeros remain bounded by `q * n + K`, relation nonzeros
remain bounded by `degree_cap * n`, and label budget does not grow with size.

What this establishes:

- complete-pipeline operation at 220,000 total nodes;
- near-linear measured runtime and retained state;
- a clear crossover over dense LeMA and KEMA;
- a later crossover over the package's faster sparse baselines;
- no hidden quadratic matching evaluator in the benchmark.

What it does not establish:

- a completed million-node full fit (the frozen attempt timed out);
- distributed, disk-backed, or GPU execution;
- a direct 10x runtime win over LeMA at 100,000 common nodes;
- real-scene accuracy at large scale;
- stable accuracy from a fixed 12-label budget as sample size grows.

The measured million-node attempt changes the diagnosis. Memory did not
approach the machine ceiling: peak RSS was 7.20 GB. A same-control sampling
profile at 100,000 common nodes attributes 57.74% of sampled time to bounded
relation construction, including 54.41% in `RANN::nn2`, and 26.39% to sparse
UOT evaluation. These two paths account for 84.13% of sampled time.

The million task has 4.545x as many common nodes as the 100,000-common-node
profile but consumed more than 8.45x its profiled fit time. Because the larger
task is right-censored, the corresponding exponent of 1.410 is a lower bound.
The next implementation work is reusable approximate relation indexes (plus a
linear spatial-grid path) and warm-started UOT with cached candidate costs.
Both need frozen small-scale accuracy/oracle checks before the exact same
million-node gate is repeated.

## Artifacts

- `rapid-ma-scaling-raw.csv`: every task, timeout/failure, metric, and sampled
  peak.
- `rapid-ma-scaling-summary.csv`: per-method/per-size aggregation.
- `rapid-ma-scaling-metadata.txt`: hardware, budgets, seeds, and source hashes.
- `rapid-ma-million-{raw,summary}.csv` and
  `rapid-ma-million-metadata.txt`: the frozen 1,000,002-node timeout.
- `rapid-ma-million-profile.txt`: the compact 100,000-common-node function
  profile and filing-ready bottleneck targets.

SHA-256:

```text
08c1afe1235deb125c35c60f7ee0c42926d4644defbcfa10e0f44ddfb7cc47f7  rapid-ma-scaling-raw.csv
e71be1565d88563813976e28b39a98eb768647ebf793c691a3a82840140839d9  rapid-ma-scaling-summary.csv
20310590440ceab6c4ccf5160012e020f58ad919e37eae4380e0878a1a322bcf  rapid-ma-scaling-metadata.txt
db14e271af34a2ce1d110d2e729ca7466cac29a2a7b98c5d96e4de53e7c362de  rapid-ma-million-raw.csv
4753eca1622f6033e15204558f756c3a68d9e0511a8b6b7b2770b7406df86328  rapid-ma-million-summary.csv
3c6c6114d0f2732bd67a2700db69253f51c63a2ee19e260e5ae714c6d5370cc7  rapid-ma-million-metadata.txt
9c2074ba81eba487bfb51c16813ae7baa1e218c8d4861dcfe5012ddf698951d7  rapid-ma-million-profile.txt
```
