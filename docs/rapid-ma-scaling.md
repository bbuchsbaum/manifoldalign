# RAPID-MA complete-pipeline scaling

> **Regenerated 2026-10-19** with the documented runner and frozen budgets
> (`RAPID_MA_PROFILE=full`, ten-second common-task cap, four-minute cap for the
> RAPID-only sizes) at commit `91f23f03c5aaca9f4a55ab7b9c1fdfdc9765e3b7`, on a
> 4-core x86_64 Linux machine (R 4.5.3, OpenBLAS, 15.7 GB) with no other jobs
> running. All 84 tasks are recorded with status `ok`, `timeout`, or `failed`.
> The superseded run used a 14-logical-core arm64 macOS machine with 36 GB; its
> summary is kept as `rapid-ma-scaling-summary-previous-macos.csv`. Because the
> caps are wall-clock limits, which tasks finish depends on the hardware: the
> tables below are a within-machine comparison and are not compared number for
> number with the previous run.
>
> A same-machine control isolates the code changes from the hardware: the
> pre-fix commit `78fbc05` ran the RAPID-MA teacher in 2.30, 8.58, and 173.3
> seconds at 1,000, 5,000, and 50,000 common nodes, against 2.17, 8.21, and
> 180.5 seconds for the regenerated code
> (`rapid-ma-scaling-control-prefix-summary.csv`). The fixes therefore change
> runtime by at most about 5%; the slower absolute times here come from the
> machine.

Status: 110,000-node synthetic run complete on this machine; the
220,000-node fit exceeded the frozen four-minute cap here (it completed in
108 seconds on the previous, faster machine); the 1,000,002-node gate failed
on the previous machine and was not repeated.

## Bottom line

RAPID-MA's complete validated-teacher pipeline scales slightly
super-linearly in runtime and linearly in retained state over the measured
range. On this machine it fit 50,000 source plus 60,000 target nodes (110,000
total) in 180.5 seconds, retained a 231.6 MB result, and reached a 1.08 GB
sampled process peak. The log-log time slope from 2,500 through 50,000 common
nodes is 1.24; the retained-object slope is 0.99.

The 100,000-common-node (220,000 total) task did not finish within the frozen
240-second cap here (sampled peak 1.81 GB when stopped), and the
10,000-common-node task exceeded the ten-second common cap at 10.02 seconds.
Both are hardware-bound timeouts: the same-machine control shows the earlier
code is no faster.

The speed claims need qualification:

- RAPID-MA is slower than SSMA and Spectral MNN at every size where all three
  finish (up to 5,000 common nodes); none of the three finishes the
  10,000-common-node task in ten seconds on this machine.
- Token-OT Graph is slightly faster than RAPID-MA from 250 to 2,500 common
  nodes and exceeds ten seconds from 5,000, where RAPID takes 8.21 seconds.
- Dense LeMA is faster than RAPID-MA at 250 common nodes (1.18 versus 1.55
  seconds), slower at 500 (4.65 versus 1.83), and exceeds ten seconds from
  1,000 onward, where RAPID takes 2.17 seconds.
- KEMA completed no task here: at 250 common nodes its certified eigensolver
  rejected all three runs (backward error 6e-7 to 1e-5), from 500 to 2,500 it
  exceeded ten seconds, and from 5,000 it refused to allocate its dense
  kernel blocks (about 8.8 GB at 5,000 common nodes) under the default
  `max_dense_bytes`. (In the same-machine control the pre-fix benchmark arm
  failed for an unrelated reason: it passed the withdrawn
  `solver = "regression"`.)
- No direct 10x runtime comparison was completed at 100,000 common nodes.

## Protocol

- Same deterministic synthetic fixture and fixed 12-label-per-domain budget at
  every size.
- Common-node sizes: 250, 500, 1,000, 2,500, 5,000, and 10,000. The target has
  20% additional unmatched nodes, so total nodes are 2.2 times the common size.
- Three repetitions through 1,000; one repetition at larger sizes.
- RAPID-only extensions at 50,000 and 100,000 common nodes.
- Every task creates its fixture and runs in an isolated process. RSS is sampled
  every 20 ms. Common tasks have a ten-second cap; RAPID large tasks have a
  four-minute cap. Timeout and failure rows remain in the data.
- R 4.5.3 on x86_64 Linux, 4 logical cores, 15.7 GB system memory,
  `OMP_NUM_THREADS=1`.
- Hits@1 is exact. Matching MRR is the declared bounded MRR@100 lower bound;
  the evaluator does not materialize an all-pairs distance matrix.

## RAPID-MA scaling

| Common nodes | Total nodes | Runtime (s) | Retained result (MB) | Sampled peak RSS (MB) | OA | Macro-F1 |
|---:|---:|---:|---:|---:|---:|---:|
| 250 | 550 | 1.550 | 1.57 | 341.9 | 0.9826 | 0.9757 |
| 500 | 1,100 | 1.831 | 2.68 | 285.0 | 0.9683 | 0.9659 |
| 1,000 | 2,200 | 2.165 | 5.05 | 331.1 | 0.9888 | 0.9850 |
| 2,500 | 5,500 | 4.795 | 12.01 | 371.1 | 0.9719 | 0.9682 |
| 5,000 | 11,000 | 8.210 | 23.53 | 398.8 | 0.9531 | 0.9577 |
| 10,000 | 22,000 | >10 (timeout) | NA | 470.1 | NA | NA |
| 50,000 | 110,000 | 180.488 | 231.64 | 1,083.6 | 0.9892 | 0.9867 |
| 100,000 | 220,000 | >240 (timeout) | NA | 1,814.6 | NA | NA |

Absolute peak RSS includes inherited R and loaded-package state. Retained state
is the cleaner scaling signal and grows essentially linearly.

## Runtime comparison

| Method | 250 | 500 | 1,000 | 2,500 | 5,000 | 10,000 |
|---|---:|---:|---:|---:|---:|---:|
| RAPID-MA teacher | 1.550 | 1.831 | 2.165 | 4.795 | 8.210 | timeout |
| Dense LeMA | 1.181 | 4.647 | timeout | timeout | timeout | timeout |
| KEMA | failed | timeout | timeout | timeout | failed | failed |
| SSMA | **0.281** | **0.332** | **0.502** | **1.520** | **5.951** | timeout |
| Spectral MNN | 0.298 | 0.387 | 0.639 | 1.977 | 7.662 | timeout |
| Token-OT Graph | 1.208 | 1.290 | 1.968 | 4.560 | timeout | timeout |
| Multiscale | 0.497 | 0.556 | 0.752 | 4.390 | 7.975 | timeout |

Times through 1,000 are three-run means; larger entries are single runs.
"Timeout" means the isolated end-to-end task exceeded ten seconds on this
machine, not that the method can never complete with a larger budget.
"Failed" rows carry their error in the raw file: KEMA's certification
rejections at 250 and its dense-allocation refusals from 5,000. Multiscale,
which failed every size in the previous run with `no right-hand side in 'b'`,
now completes through 5,000 common nodes.

## Dense-memory contrast

The LeMA reference's learned graph is dense over target nodes plus paired
source landmarks. With 100,000 common nodes and 20% target-only expansion, that
graph has about 120,012 rows. One double-precision matrix of that shape requires
about 115.2 GB; retaining both `W` and `L` alone requires about 230.4 GB before
distance workspaces and embeddings, far beyond either machine.

At 10,000 common nodes, LeMA had reached 2.32 GB sampled RSS when the
ten-second timeout stopped it on this machine. RAPID's sparse retained result
is 231.6 MB at 50,000 common nodes.

## Scaling contracts and limitations

Automated tests double fixture size and check that retained state grows by less
than 3x, transport nonzeros remain bounded by `q * n + K`, relation nonzeros
remain bounded by `degree_cap * n`, and label budget does not grow with size.

What this establishes:

- complete-pipeline operation at 110,000 total nodes on a 4-core machine
  (and at 220,000 on the previous 14-core machine);
- slightly super-linear runtime (slope 1.24) and linear retained state;
- a crossover over dense LeMA between 250 and 500 common nodes;
- no crossover over SSMA or Spectral MNN within the ten-second budget on this
  machine (the previous machine showed one around 2,500 to 5,000 common
  nodes);
- no hidden quadratic matching evaluator in the benchmark.

What it does not establish:

- a completed million-node full fit (the frozen attempt timed out);
- distributed, disk-backed, or GPU execution;
- a direct 10x runtime win over LeMA at 100,000 common nodes;
- real-scene accuracy at large scale;
- stable accuracy from a fixed 12-label budget as sample size grows.

The following diagnosis comes from the previous machine and was not
repeated here. The measured million-node attempt changes the diagnosis. Memory did not
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
- `rapid-ma-scaling-metadata.txt`: hardware and budgets (the runner was
  executed outside the package root, so commit and source hashes are in
  `rapid-ma-regeneration-provenance.txt`).
- `rapid-ma-scaling-control-prefix-summary.csv`: same-machine control with
  the pre-fix commit.
- `rapid-ma-scaling-summary-previous-macos.csv`: the superseded summary from
  the previous machine.
- `rapid-ma-million-{raw,summary}.csv` and
  `rapid-ma-million-metadata.txt`: the frozen 1,000,002-node timeout
  (previous machine, not repeated).
- `rapid-ma-million-profile.txt`: the compact 100,000-common-node function
  profile and filing-ready bottleneck targets.

SHA-256:

```text
137b8dda967f6feb90e38ea6657790f6993523261797dab7c8d9005f558799dc  rapid-ma-scaling-raw.csv
e1f74abf6bebe996f1178775f38d23b218cf5aaae023b93812e263f7a56a69be  rapid-ma-scaling-summary.csv
abd001a90f7112ec1e812d8f58a3c5881da18c8ae6393763304f9001084ba28a  rapid-ma-scaling-metadata.txt
985e87c03c1b7412a24cde3e17fb6b2e1db5f1812cbd5f39358d4b66d3393428  rapid-ma-scaling-control-prefix-summary.csv
e71be1565d88563813976e28b39a98eb768647ebf793c691a3a82840140839d9  rapid-ma-scaling-summary-previous-macos.csv
fffd98245276f7c51872ef985d7ceaae23bdb693aed42bfd047218ffa25d50d7  rapid-ma-regeneration-provenance.txt
db14e271af34a2ce1d110d2e729ca7466cac29a2a7b98c5d96e4de53e7c362de  rapid-ma-million-raw.csv
4753eca1622f6033e15204558f756c3a68d9e0511a8b6b7b2770b7406df86328  rapid-ma-million-summary.csv
3c6c6114d0f2732bd67a2700db69253f51c63a2ee19e260e5ae714c6d5370cc7  rapid-ma-million-metadata.txt
9c2074ba81eba487bfb51c16813ae7baa1e218c8d4861dcfe5012ddf698951d7  rapid-ma-million-profile.txt
```
